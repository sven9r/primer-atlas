"""Source identity, taxonomy, sampling and immutable-output UNITE regressions.

All small test sequences are synthetic and do not establish marker coverage.
"""
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import random
import tarfile
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("unite", Path(__file__).resolve().parents[1] / "scripts/build_unite_reference_panel.py")
unite = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(unite)


def dna(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def record(number, phylum="PhylumA", class_name="ClassA", order="OrderA", family="FamilyA",
           genus="GenusA", representative="reps", length=600, sequence=None, accession=None, sh=None):
    accession = accession or "SYNTHETIC" + str(number)
    sh = sh or "SH_SYNTHETIC" + str(number)
    header = (f"Genus_species|{accession}|{sh}|{representative}|k__Fungi;p__{phylum};c__{class_name};"
              f"o__{order};f__{family};g__{genus};s__Genus_species")
    return header, sequence if sequence is not None else dna(length, number)


def write_fasta(path, records):
    path.write_text("".join(">" + header + "\n" + sequence + "\n" for header, sequence in records))


def bundle(directory, standard, developer=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    std, dev, archive = directory / "standard.fasta", directory / "developer.fasta", directory / "source.tgz"
    write_fasta(std, standard)
    write_fasta(dev, developer if developer is not None else standard)
    with tarfile.open(archive, "w:gz") as package:
        package.add(std, arcname=std.name)
        package.add(dev, arcname=dev.name)
    expected = {"archive": unite.sha_file(archive), "standard": unite.sha_file(std), "developer": unite.sha_file(dev)}
    return archive, std, dev, expected


class ImportTests(unittest.TestCase):
    def test_developer_identity_is_joined_by_sh_and_accession_not_file_order(self):
        first, second = record(1), record(2)
        developer = [(second[0], "ACGT" + second[1] + "TGCA"), (first[0], "GG" + first[1] + "CC")]
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, [first, second], developer)
            rows, sequences, issues, counts = unite.pair_sources(std, dev)
        by_accession = {r["accession"]: r for r in rows}
        self.assertEqual(counts, dict(standard_records=2, developer_records=2))
        self.assertEqual(by_accession["SYNTHETIC1"]["raw_length_bp"], 604)
        self.assertEqual(by_accession["SYNTHETIC2"]["raw_length_bp"], 608)
        self.assertEqual(sequences[by_accession["SYNTHETIC1"]["record_id"]], "GG" + first[1] + "CC")
        self.assertTrue(all(r["eligibility_status"] == "eligible" for r in rows))
        self.assertEqual(issues, [])
        self.assertEqual(by_accession["SYNTHETIC1"]["accession_version_status"], "not_reported")
        self.assertFalse(by_accession["SYNTHETIC1"]["accession"].endswith(".1"))

    def test_taxonomy_placeholders_and_incertae_sedis_are_not_resolved_orders(self):
        for value in ("", "unidentified", "unclassified_Fungi", "unknown", "uncultured_fungus", "NA"):
            self.assertEqual(unite.rank_status(value), "unresolved")
        for value in ("Fungi_ord_Incertae_sedis", "GS01_phy_Incertae_sedis", "Incertae sedis"):
            self.assertEqual(unite.rank_status(value), "incertae_sedis")
        self.assertEqual(unite.rank_status("Hypocreales"), "resolved")
        row = unite.parse_header(record(1, order="Fungi_ord_Incertae_sedis")[0])
        self.assertEqual(row["order_name_status"], "incertae_sedis")

    def test_hierarchy_ids_keep_homonymous_orders_and_families_separate(self):
        left = unite.parse_header(record(1, phylum="A")[0])
        right = unite.parse_header(record(2, phylum="B")[0])
        self.assertNotEqual(unite.group_identity(left, "order"), unite.group_identity(right, "order"))
        self.assertNotEqual(unite.group_identity(left, "family"), unite.group_identity(right, "family"))
        self.assertEqual(unite.group_identity(left, "class"), unite.group_identity(dict(left, accession="OTHER"), "class"))

    def test_uncertain_sequence_relationship_is_retained_but_not_selected(self):
        standard = record(1)
        developer = (standard[0], dna(590, 33))
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, [standard], [developer])
            rows, _, issues, _ = unite.pair_sources(std, dev)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["eligibility_status"], "withheld")
        self.assertIn("unverified_standard_developer_sequence_relationship", rows[0]["eligibility_reason"])
        self.assertEqual(rows[0]["source_header"], standard[0])
        self.assertEqual(rows[0]["source_developer_header"], developer[0])
        self.assertEqual(issues[0]["issue"], "withheld_reference")

    def test_reverse_complement_relationship_is_disclosed_without_reorienting_source(self):
        standard = record(1)
        reverse = standard[1].translate(unite.COMPLEMENT)[::-1]
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, [standard], [(standard[0], "AAAA" + reverse + "CCCC")])
            rows, sequences, _, _ = unite.pair_sources(std, dev)
        self.assertEqual(rows[0]["sequence_relationship"], "standard_reverse_complement_embedded_in_developer")
        self.assertEqual(sequences[rows[0]["record_id"]], "AAAA" + reverse + "CCCC")

    def test_missing_counterpart_is_not_paired_to_an_unrelated_record(self):
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, [record(1)], [record(2)])
            rows, _, _, _ = unite.pair_sources(std, dev)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r["eligibility_status"] == "withheld" for r in rows))
        self.assertTrue(any("missing_developer_record" in r["eligibility_reason"] for r in rows))
        self.assertTrue(any("missing_standard_record" in r["eligibility_reason"] for r in rows))

    def test_duplicate_occurrences_are_preserved_and_conflicts_withheld(self):
        original = record(1)
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, [original, original], [original])
            rows, _, issues, _ = unite.pair_sources(std, dev)
            self.assertEqual(rows[0]["standard_record_count"], 2)
            self.assertEqual(rows[0]["standard_record_ordinals"], "1|2")
            self.assertEqual(rows[0]["eligibility_status"], "eligible")
            self.assertEqual(issues[0]["issue"], "identical_duplicate")
        with tempfile.TemporaryDirectory() as directory:
            changed = (original[0], dna(600, 88))
            _, std, dev, _ = bundle(directory, [original, changed], [original])
            rows, _, issues, _ = unite.pair_sources(std, dev)
            self.assertEqual(rows[0]["eligibility_status"], "withheld")
            self.assertIn("standard_identity_conflict", rows[0]["eligibility_reason"])
            self.assertTrue(any(r["issue"] == "conflicting_duplicate" for r in issues))

    def test_one_sh_cannot_silently_change_accession_or_taxonomy(self):
        with tempfile.TemporaryDirectory() as directory:
            entries = [record(1, sh="SH_CONFLICT"), record(2, sh="SH_CONFLICT")]
            _, std, dev, _ = bundle(directory, entries)
            rows, _, _, _ = unite.pair_sources(std, dev)
            self.assertTrue(all("species_hypothesis_accession_conflict" in r["eligibility_reason"] for r in rows))
        with tempfile.TemporaryDirectory() as directory:
            original = record(1)
            dev = record(1, order="DifferentOrder", sequence=original[1])
            _, std, developer, _ = bundle(directory, [original], [dev])
            rows, _, _, _ = unite.pair_sources(std, developer)
            self.assertIn("standard_developer_header_conflict", rows[0]["eligibility_reason"])

    def test_invalid_alphabet_and_cross_order_identical_sequences_are_withheld(self):
        with tempfile.TemporaryDirectory() as directory:
            entries = [record(1, sequence="ACGT-X"), record(2, order="A", sequence=dna(600, 7)),
                       record(3, order="B", sequence=dna(600, 7))]
            _, std, dev, _ = bundle(directory, entries)
            rows, _, _, _ = unite.pair_sources(std, dev)
        self.assertTrue(all(r["eligibility_status"] == "withheld" for r in rows))
        self.assertIn("-", rows[0]["sequence_alphabet"])
        self.assertEqual(sum(r["eligibility_reason"] == "cross_order_identical_sequence_conflict" for r in rows), 2)


class SamplingTests(unittest.TestCase):
    def rows(self, entries):
        with tempfile.TemporaryDirectory() as directory:
            _, std, dev, _ = bundle(directory, entries)
            rows, _, _, _ = unite.pair_sources(std, dev)
        return rows

    def test_all_rare_orders_and_uncertain_order_groups_remain_visible(self):
        rows = self.rows([record(1, order="OrderRare", length=550), record(2, order="OrderRare", length=620),
                         record(3, order="unidentified", length=1000),
                         record(4, order="Fungi_ord_Incertae_sedis", length=250)])
        selected, _ = unite.choose_panel(rows, 100)
        coverage = unite.coverage_rows(rows, selected, 100)
        orders = [r for r in coverage if r["rank"] == "order"]
        self.assertEqual(len(selected), 4)
        self.assertEqual(len(orders), 3)
        self.assertEqual(sum(r["taxon_resolved"] for r in orders), 1)
        self.assertTrue(all(r["n_selected"] == r["n_eligible_unique"] for r in orders))
        self.assertTrue(all(r["shortfall_to_100"] == 100 - r["n_selected"] for r in orders))
        for row in selected:
            for rank in unite.COVERAGE_RANKS:
                match = [c for c in coverage if c["rank"] == rank and c["group_id"] == row[rank + "_group_id"]]
                self.assertEqual(len(match), 1)

    def test_selection_is_deterministic_and_retains_broad_length_strata(self):
        rows = self.rows([record(1, length=280), record(2, length=550), record(3, length=1100),
                         record(4, length=580), record(5, length=570)])
        a, _ = unite.choose_panel(rows, 3)
        b, _ = unite.choose_panel(list(reversed(rows)), 3)
        self.assertEqual([r["record_id"] for r in a], [r["record_id"] for r in b])
        self.assertEqual({r["length_stratum"] for r in a}, {"lt300", "500_649", "1000_1999"})

    def test_refs_lower_ambiguity_and_unique_sequence_policy(self):
        common = dna(600, 51)
        rows = self.rows([record(1, representative="reps", sequence=common),
                         record(2, representative="refs", sequence=common),
                         record(3, representative="refs", sequence=dna(600, 52)[:580] + "N" * 20)])
        selected, unique = unite.select_group(rows, 1)
        self.assertEqual(unique, 2)
        self.assertEqual(selected[0]["accession"], "SYNTHETIC2")
        self.assertEqual(selected[0]["ambiguity_count"], 0)


class BuildTests(unittest.TestCase):
    def build(self, directory, output_name="panel"):
        archive, std, dev, expected = bundle(Path(directory) / "source", [record(1), record(2, order="unidentified")])
        output = Path(directory) / output_name
        result = unite.build_panel(archive, std, dev, output, expected=expected, expected_records=2,
                                  source_doi="synthetic:test", source_version="synthetic-v1", source_url="https://example.invalid/source")
        return result, output, (archive, std, dev, expected)

    def test_atomic_output_idempotency_hashes_gzip_and_no_invented_accession_versions(self):
        with tempfile.TemporaryDirectory() as directory:
            result, output, inputs = self.build(directory)
            before = {p.name: (unite.sha_file(p), p.stat().st_mtime_ns) for p in output.iterdir()}
            archive, std, dev, expected = inputs
            replay = unite.build_panel(archive, std, dev, output, expected=expected, expected_records=2,
                                      source_doi="synthetic:test", source_version="synthetic-v1", source_url="https://example.invalid/source")
            self.assertEqual(replay, result)
            self.assertEqual(before, {p.name: (unite.sha_file(p), p.stat().st_mtime_ns) for p in output.iterdir()})
            for name, descriptor in result["artifacts"].items():
                self.assertEqual(unite.sha_file(output / name), descriptor["sha256"])
            with gzip.open(output / "selected_metadata.csv.gz", "rt") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual({r["accession"] for r in rows}, {"SYNTHETIC1", "SYNTHETIC2"})
            self.assertTrue(all(r["panel_version_id"] == result["version_id"] for r in rows))
            with gzip.open(output / "selected_raw.fasta.gz", "rt") as stream:
                headers = [line.strip()[1:] for line in stream if line.startswith(">")]
            self.assertEqual(headers, [r["record_id"] for r in rows])
            self.assertEqual(result["counts"]["resolved_order_groups"], 1)
            self.assertEqual(result["counts"]["unresolved_order_groups"], 1)
            with self.assertRaisesRegex(ValueError, "Different immutable panel"):
                unite.build_panel(archive, std, dev, output, per_group=1, expected=expected, expected_records=2,
                                  source_doi="synthetic:test", source_version="synthetic-v1")

    def test_existing_artifact_corruption_and_archive_member_mismatch_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            _, output, inputs = self.build(directory)
            archive, std, dev, expected = inputs
            (output / "source_coverage.csv").write_text("corrupt")
            with self.assertRaisesRegex(ValueError, "artifact checksum differs"):
                unite.build_panel(archive, std, dev, output, expected=expected, expected_records=2,
                                  source_doi="synthetic:test", source_version="synthetic-v1")
        with tempfile.TemporaryDirectory() as directory:
            archive, std, dev, expected = bundle(directory, [record(1)])
            std.write_text(std.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "source checksum differs"):
                unite.validate_bundle(archive, std, dev, expected, 1)
            expected["standard"] = unite.sha_file(std)
            with self.assertRaisesRegex(ValueError, "not the pinned archive member"):
                unite.validate_bundle(archive, std, dev, expected, 1)

    def test_incomplete_write_never_promotes_or_leaves_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, std, dev, expected = bundle(Path(directory) / "source", [record(1)])
            output = Path(directory) / "panel"
            with patch.object(unite, "write_csv", side_effect=OSError("synthetic disk error")):
                with self.assertRaisesRegex(OSError, "disk error"):
                    unite.build_panel(archive, std, dev, output, expected=expected, expected_records=1)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(directory).glob(".unite-panel-*")), [])

    def test_gzip_artifact_bytes_reproduce_in_another_output_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            result, output, inputs = self.build(directory)
            archive, std, dev, expected = inputs
            other = Path(directory) / "second"
            replay = unite.build_panel(archive, std, dev, other, expected=expected, expected_records=2,
                                      source_doi="synthetic:test", source_version="synthetic-v1", source_url="https://example.invalid/source")
            self.assertEqual(result["version_id"], replay["version_id"])
            self.assertEqual(result["artifacts"], replay["artifacts"])


if __name__ == "__main__":
    unittest.main()
