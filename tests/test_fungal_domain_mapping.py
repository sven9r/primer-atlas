"""Regression checks for sourced, feature-specific fungal LSU domain projection.

Fixtures are synthetic. They establish coordinate/support invariants rather than
claiming any biological domain boundary in a real fungal reference.
"""
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("domain_mapper", ROOT / "scripts/map_fungal_domains.py")
mapper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mapper)


def synthetic_sequence(length=200, seed=117):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def feature(location, product="28S ribosomal RNA", key="rRNA"):
    return ("<GBFeature><GBFeature_key>" + key + "</GBFeature_key>"
            "<GBFeature_location>" + escape(location) + "</GBFeature_location>"
            "<GBFeature_quals><GBQualifier><GBQualifier_name>product</GBQualifier_name>"
            "<GBQualifier_value>" + escape(product) + "</GBQualifier_value>"
            "</GBQualifier></GBFeature_quals></GBFeature>")


def records_xml(accession, sequence, features):
    return ("<GBSet><GBSeq><GBSeq_accession-version>" + accession + "</GBSeq_accession-version>"
            "<GBSeq_sequence>" + sequence.lower() + "</GBSeq_sequence>"
            "<GBSeq_feature-table>" + "".join(features) + "</GBSeq_feature-table>"
            "</GBSeq></GBSet>")


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.reference = synthetic_sequence()

    def project(self, reference=None, target=None, start=81, end=100, **kwargs):
        reference = self.reference if reference is None else reference
        target = reference if target is None else target
        return mapper.project_domain(reference, target, target.replace("-", ""), start, end, **kwargs)

    def test_inclusive_interval_is_raw_bases_not_alignment_width(self):
        plain = self.project()
        self.assertEqual(plain["status"], "supported")
        self.assertEqual((plain["start"], plain["end"]), (81, 100))
        reference = self.reference[:90] + "-----" + self.reference[90:]
        target = self.reference[:90] + "GATTA" + self.reference[90:]
        inserted = self.project(reference, target)
        self.assertEqual(inserted["status"], "supported")
        self.assertEqual((inserted["start"], inserted["end"]), (81, 105))
        self.assertEqual(inserted["end"] - inserted["start"] + 1, 25)

    def test_deleted_boundary_bases_project_to_retained_edges(self):
        target = list(self.reference)
        target[80] = target[99] = "-"
        proof = self.project(target="".join(target))
        self.assertEqual(proof["status"], "supported")
        self.assertEqual((proof["start"], proof["end"]), (81, 98))
        self.assertEqual(proof["end"] - proof["start"] + 1, 18)
        self.assertEqual(proof["boundary_method"], "alignment_column_edges_including_target_indels")

    def test_deleted_whole_domain_is_withheld(self):
        target = self.reference[:80] + "-" * 20 + self.reference[100:]
        proof = self.project(target=target)
        self.assertEqual(proof["status"], "no_target_bases_in_domain")
        self.assertLess(proof["end"], proof["start"])

    def test_actual_flank_lengths_control_coverage_and_scores(self):
        target = list(self.reference)
        for position in (34, 45, 60):
            target[position] = "-"
        proof = self.project(target="".join(target), left_bases=48, right_bases=33)
        self.assertEqual(proof["status"], "supported")
        self.assertEqual(proof["left_flank"]["source_flank_bases"], 48)
        self.assertEqual(proof["right_flank"]["source_flank_bases"], 33)
        self.assertEqual(proof["left_flank"]["paired_coverage"], 45 / 48)
        self.assertEqual(proof["left_flank"]["paired_identity"], 1)

    def test_observed_anchor_insertions_pay_the_same_penalty_as_competitors(self):
        reference = self.reference[:48] + "-----" + self.reference[48:]
        target = self.reference[:48] + "GATTA" + self.reference[48:]
        proof = self.project(reference, target)
        self.assertEqual(proof["status"], "supported")
        left = proof["left_flank"]
        self.assertEqual(left["paired_coverage"], 1)
        self.assertEqual(left["paired_identity"], 1)
        self.assertEqual(left["observed_alignment_score"], 2 * 64 - 2 * 5)
        self.assertEqual(left["observed_normalized_alignment_score"], 118 / 128)
        self.assertEqual((proof["start"], proof["end"]), (86, 105))

    def test_partial_source_or_too_short_conserved_context_is_withheld(self):
        self.assertEqual(self.project(start=32, end=45)["status"], "source_flank_unavailable")
        self.assertEqual(self.project(left_bases=31)["status"], "source_flank_unavailable")
        target = "-" * 90 + self.reference[90:]
        self.assertEqual(self.project(target=target)["status"], "boundary_flanks_unresolved")

    def test_ambiguous_target_domain_never_gets_known_boundaries(self):
        target = self.reference[:85] + "N" + self.reference[86:]
        self.assertEqual(self.project(target=target)["status"], "ambiguous_target_domain")

    def test_ambiguous_flank_is_withheld_even_with_high_identity(self):
        target = self.reference[:50] + "N" + self.reference[51:]
        proof = self.project(target=target)
        self.assertEqual(proof["status"], "boundary_flanks_unresolved")
        self.assertGreater(proof["left_flank"]["paired_identity"], 0.99)
        self.assertFalse(proof["left_flank"]["supported"])

    def test_exact_alternative_rDNA_copy_rejects_unique_anchor_claim(self):
        spacer = synthetic_sequence(96, seed=19)
        target = self.reference + spacer + self.reference
        reference = self.reference + "-" * (len(spacer) + len(self.reference))
        proof = self.project(reference, target)
        self.assertEqual(proof["status"], "boundary_flanks_unresolved")
        self.assertEqual(proof["left_flank"]["best_alternative_identity"], 1)
        self.assertEqual(proof["right_flank"]["best_alternative_identity"], 1)
        self.assertFalse(proof["left_flank"]["supported"])

    def test_indel_bearing_alternative_copy_also_rejects_unique_anchor_claim(self):
        spacer = synthetic_sequence(96, seed=19)
        duplicate = self.reference[:48] + "A" + self.reference[48:132] + "C" + self.reference[132:]
        target = self.reference + spacer + duplicate
        reference = self.reference + "-" * (len(spacer) + len(duplicate))
        proof = self.project(reference, target)
        self.assertEqual(proof["status"], "boundary_flanks_unresolved")
        for side in ("left_flank", "right_flank"):
            flank = proof[side]
            self.assertEqual(flank["anchor_margin"], 0.015625)
            self.assertLess(flank["anchor_margin"], mapper.MIN_ANCHOR_MARGIN)
            competition = flank["competitive_anchor_proof"]
            self.assertEqual(competition["method"], "full_query_semiglobal_DP")
            alternative = competition["best_alternative"]
            self.assertEqual(alternative["paired_identity"], 1)
            self.assertEqual(alternative["paired_coverage"], 1)
            self.assertEqual(alternative["matches"], 64)
            self.assertEqual(alternative["score"], 126)
            self.assertGreater(alternative["target_start"], len(self.reference) + len(spacer))
            self.assertEqual(alternative["target_end"] - alternative["target_start"] + 1, 65)

    def test_high_match_fraction_but_heavily_gapped_decoy_cannot_equal_a_better_anchor(self):
        own = list(self.reference)
        # The own left flank has51/64 matches and an ungapped score89.
        for position in range(16, 80, 5):
            own[position] = "ACGT"[("ACGT".index(own[position]) + 1) % 4]
        own = "".join(own)
        anchor = own[16:80]
        decoy = "".join(anchor[i:i + 4] + ("A" if i + 4 < 64 else "") for i in range(0, 64, 4))
        spacer = synthetic_sequence(96, seed=19)
        target = own + spacer + decoy
        reference = self.reference + "-" * (len(spacer) + len(decoy))
        proof = mapper.flank_proof(reference, target, range(17, 81), target)
        alternative = proof["competitive_anchor_proof"]["best_alternative"]
        self.assertGreater(alternative["matches"] / 64, proof["paired_identity"])
        self.assertEqual(proof["observed_alignment_score"], 89)
        self.assertLess(alternative["score"], proof["observed_alignment_score"])
        self.assertGreater(proof["anchor_margin"], mapper.MIN_ANCHOR_MARGIN)
        self.assertTrue(proof["supported"])
        self.assertEqual(proof["coordinate_system"], "one_based_inclusive_target_feature_bases")

    def test_partial_coverage_low_identity_alternative_is_not_competitive(self):
        proof = mapper.competitive_anchor(self.reference[:64], "T" * 200, 0, 64)
        self.assertIsNone(proof["best_alternative"])

    def test_alignment_width_disagreement_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "width differs"):
            self.project(target=self.reference[:-1])


class SourceFeatureTests(unittest.TestCase):
    def test_repeated_copies_remain_individual_features_and_partial_is_explicit(self):
        sequence = synthetic_sequence(440)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "records.xml"
            path.write_text(records_xml("COPY.1", sequence, [feature("<1..200"), feature("241..>440")]))
            rows = mapper.source_features(path, {"COPY.1": sequence})
        self.assertEqual([(r["start"], r["end"]) for r in rows], [(1, 200), (241, 440)])
        self.assertEqual([r["source_feature_number"] for r in rows], [0, 1])
        self.assertTrue(all(r["partial"] for r in rows))
        self.assertEqual([r["source_location"] for r in rows], ["<1..200", "241..>440"])

    def test_complement_compound_or_operon_labels_are_not_fabricated_as_forward_lsu(self):
        sequence = synthetic_sequence(200)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "records.xml"
            path.write_text(records_xml("COPY.1", sequence, [
                feature("complement(1..200)"), feature("join(1..80,121..200)"),
                feature("1..200", "contains small subunit, 5.8S, and large subunit", key="misc_feature")]))
            self.assertEqual(mapper.source_features(path, {"COPY.1": sequence}), [])

    def test_annotation_record_must_match_frozen_accession_sequence(self):
        sequence = synthetic_sequence(200)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "records.xml"
            path.write_text(records_xml("COPY.1", sequence, [feature("1..200")]))
            with self.assertRaisesRegex(ValueError, "differs from frozen sequence"):
                mapper.source_features(path, {"COPY.1": "T" + sequence[1:]})


class GeneratorTests(unittest.TestCase):
    def fixture(self, directory, partial=False):
        directory = Path(directory)
        snapshot, evidence, output = directory / "snapshot", directory / "evidence", directory / "output"
        snapshot.mkdir(); evidence.mkdir()
        reference = synthetic_sequence()
        target = reference + "N" * 40 + reference
        (snapshot / "sequences.fasta").write_text(">COPY.1 synthetic two-copy record\n" + target + "\n")
        (snapshot / "landmarks.csv").write_text("accession,landmark_id,label,start,end,source_note\n")
        locations = ["<1..200", "241..>440"] if partial else ["1..200", "241..440"]
        records = evidence / "features.xml"
        records.write_text(records_xml("COPY.1", target, [feature(location) for location in locations]))
        (snapshot / "manifest.json").write_text(json.dumps({"id": "synthetic-v2", "artifacts": {
            "sequences": {"sha256": mapper.sha(snapshot / "sequences.fasta")},
            "landmarks": {"sha256": mapper.sha(snapshot / "landmarks.csv")}}}))
        (snapshot / "run_manifest.json").write_text(json.dumps({"raw_responses": {"features.xml": mapper.sha(records)}}))
        (evidence / "STRUCTURAL.1.fasta").write_text(">STRUCTURAL.1 synthetic source\n" + reference + "\n")
        row = {"domain_id": "D1", "source_accession": "STRUCTURAL.1", "source_start": 81, "source_end": 100,
               "source_sequence_sha256": hashlib.sha256(reference.encode()).hexdigest(),
               "source_scheme": "synthetic classic-D fixture", "source_status": "source_annotated",
               "eligible_for_projection": True, "source_url": "https://example.invalid/synthetic-reference",
               "source_citation": "Synthetic source fixture; not a biological domain annotation"}
        withheld = dict(row, domain_id="D10", source_start=120, source_end=130,
                        eligible_for_projection=False, source_status="withheld_overlapping_source_annotation")
        domains = evidence / "domains.csv"
        mapper.write_csv(domains, [row, withheld])
        return dict(snapshot=snapshot, records=records, domains=domains, output=output, evidence=evidence)

    def run_fixture(self, fixture, projector=None, version_stdout="", version_stderr="synthetic-MAFFT-version"):
        argv = ["map_fungal_domains.py", "--snapshot", str(fixture["snapshot"]), "--records", str(fixture["records"]),
                "--domains", str(fixture["domains"]), "--output", str(fixture["output"])]
        def fake_mafft(command, **kwargs):
            if command == ["mafft", "--version"]:
                return subprocess.CompletedProcess(command, 0, stdout=version_stdout, stderr=version_stderr)
            pair = mapper.fasta(Path(command[-1]))
            self.assertEqual(pair["STRUCTURAL_REFERENCE"], pair["TARGET_FEATURE"])
            return subprocess.CompletedProcess(command, 0, stdout=Path(command[-1]).read_text(), stderr="")
        with patch.object(sys, "argv", argv), patch.object(mapper.subprocess, "run", side_effect=fake_mafft), contextlib.redirect_stdout(io.StringIO()):
            if projector is None:
                mapper.main()
            else:
                with patch.object(mapper, "project_domain", side_effect=projector):
                    mapper.main()

    def test_two_complete_copies_get_distinct_parent_genes_and_absolute_domain_offsets(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            inputs = [p for d in [fixture["snapshot"], fixture["evidence"]] for p in d.iterdir()]
            original = {str(p): mapper.sha(p) for p in inputs}
            self.run_fixture(fixture)
            rows = mapper.csv_rows(fixture["output"] / "landmarks.csv")
            domains = [r for r in rows if r["landmark_type"] == "divergent_domain"]
            parents = {(r["accession"], r["landmark_id"]): r for r in rows if r["landmark_type"] == "gene"}
            self.assertEqual([(int(r["start"]), int(r["end"])) for r in domains], [(81, 100), (321, 340)])
            self.assertEqual(len(parents), 2)
            self.assertNotEqual(domains[0]["parent_landmark"], domains[1]["parent_landmark"])
            for row in domains:
                parent = parents[row["accession"], row["parent_landmark"]]
                self.assertLessEqual(int(parent["start"]), int(row["start"]))
                self.assertLessEqual(int(row["end"]), int(parent["end"]))
                self.assertEqual(row["source_reference"], "STRUCTURAL.1")
                self.assertEqual(row["source_citation"], "Synthetic source fixture; not a biological domain annotation")
                self.assertEqual(row["evidence_status"], "mapped_reference_domain")
                self.assertEqual(row["boundary_status"], "projected")
            audit = json.loads((fixture["output"] / "domain_mapping_audit.json").read_text())
            self.assertEqual(sum(r["status"] == "withheld_overlapping_source_annotation" for r in audit), 2)
            for row in audit:
                if row["status"] != "mapped_reference_domain":
                    continue
                for projection in row["methods"].values():
                    self.assertEqual(projection["coordinate_system"], "one_based_inclusive_target_feature_bases")
                    self.assertEqual(row["start"], row["feature_start"] + projection["start"] - 1)
                    self.assertEqual(row["end"], row["feature_start"] + projection["end"] - 1)
                    self.assertEqual(projection["left_flank"]["coordinate_system"], projection["coordinate_system"])
            self.assertFalse(any(r["label"] == "D10" for r in domains))
            run = json.loads((fixture["output"] / "run_manifest.json").read_text())
            self.assertFalse(run["source_is_direct_fungal_annotation"])
            self.assertEqual(run["base_snapshot"], "synthetic-v2")
            self.assertEqual(run["mapped_domain_rows"], 2)
            self.assertEqual(run["mafft_version"], "synthetic-MAFFT-version")
            for relative, expected in run["artifacts_sha256"].items():
                self.assertEqual(mapper.sha(fixture["output"] / relative), expected)
            self.assertEqual({str(p): mapper.sha(p) for p in inputs}, original)

    def test_partial_parent_stays_partial_when_an_interior_domain_maps(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory, partial=True)
            self.run_fixture(fixture)
            rows = mapper.csv_rows(fixture["output"] / "landmarks.csv")
            parents = [r for r in rows if r["landmark_type"] == "gene"]
            domains = [r for r in rows if r["landmark_type"] == "divergent_domain"]
            self.assertEqual(len(parents), 2)
            self.assertEqual(len(domains), 2)
            self.assertTrue(all(r["boundary_status"] == "partial" for r in parents))
            self.assertTrue(all(r["evidence_status"] == "partial_source_annotation" for r in parents))
            self.assertEqual([r["source_location"] for r in parents], ["<1..200", "241..>440"])

    def test_unselected_source_domain_masks_flanks_without_creating_a_projection_or_audit_row(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            rows = mapper.csv_rows(fixture["domains"])
            # This source D2 ends only 24 bases before D1. Mapping D1 must be
            # withheld for unavailable conserved context, although D2 itself
            # is merely retained as boundary evidence and is not selected.
            context = dict(rows[0], domain_id="D2", source_start=30, source_end=56,
                           projection_selected=False)
            rows[0]["projection_selected"] = True
            rows[1]["projection_selected"] = True
            mapper.write_csv(fixture["domains"], rows + [context])
            self.run_fixture(fixture)
            audit = json.loads((fixture["output"] / "domain_mapping_audit.json").read_text())
            self.assertEqual(len(audit), 4)
            self.assertFalse(any(r["domain"] == "D2" for r in audit))
            selected = [r for r in audit if r["domain"] == "D1"]
            self.assertEqual(len(selected), 2)
            self.assertTrue(all(r["status"] == "boundary_mapping_unresolved" for r in selected))
            self.assertTrue(all(r["methods"]["localpair"]["status"] == "source_flank_unavailable" for r in selected))
            landmarks = mapper.csv_rows(fixture["output"] / "landmarks.csv")
            self.assertFalse(any(r["landmark_type"] == "divergent_domain" for r in landmarks))

    def test_one_base_method_disagreement_does_not_create_domain_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            calls = []
            def disagreement(*args, **kwargs):
                calls.append(1)
                return {"status": "supported", "start": 81 + (len(calls) % 2), "end": 100}
            self.run_fixture(fixture, projector=disagreement)
            rows = mapper.csv_rows(fixture["output"] / "landmarks.csv")
            self.assertFalse(any(r["landmark_type"] == "divergent_domain" for r in rows))
            audit = json.loads((fixture["output"] / "domain_mapping_audit.json").read_text())
            self.assertEqual(sum(r["status"] == "alignment_boundary_disagreement" for r in audit), 2)

    def test_unfrozen_source_response_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            fixture["records"].write_text(fixture["records"].read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "not the frozen acquisition source"):
                self.run_fixture(fixture)
            self.assertFalse((fixture["output"] / "landmarks.csv").exists())

    def test_version_printed_to_stdout_is_frozen_and_empty_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            self.run_fixture(fixture, version_stdout="synthetic-version-on-stdout", version_stderr="")
            run = json.loads((fixture["output"] / "run_manifest.json").read_text())
            self.assertEqual(run["mafft_version"], "synthetic-version-on-stdout")
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            with self.assertRaisesRegex(ValueError, "MAFFT version could not be recorded"):
                self.run_fixture(fixture, version_stdout="", version_stderr="")
            self.assertFalse((fixture["output"] / "landmarks.csv").exists())

    def test_structural_sequence_must_match_the_domain_evidence_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            source = fixture["evidence"] / "STRUCTURAL.1.fasta"
            source.write_text(">STRUCTURAL.1\n" + "A" * 200 + "\n")
            with self.assertRaisesRegex(ValueError, "differs from its domain evidence"):
                self.run_fixture(fixture)
            self.assertFalse((fixture["output"] / "landmarks.csv").exists())

    def test_named_structural_source_cannot_be_substituted_by_another_record(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            source = fixture["evidence"] / "STRUCTURAL.1.fasta"
            source.write_text(">OTHER.1\n" + synthetic_sequence() + "\n")
            with self.assertRaisesRegex(ValueError, "Missing named structural reference"):
                self.run_fixture(fixture)
            self.assertFalse((fixture["output"] / "landmarks.csv").exists())


if __name__ == "__main__":
    unittest.main()
