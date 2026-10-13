"""Immutable UNITE freezer regressions using explicitly synthetic evidence.

No synthetic feature interval, HMM call or MSA establishes biological coverage.
The fixtures exercise stage identity, exact raw slices and historical replay.
"""
from collections import Counter, defaultdict
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import random
import shutil
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_unite_reference_panel as source_builder
import freeze_unite_snapshot as freezer
import scan_unite_primer_sites as scanner

F = "TCCGTAGGTGAACCTGCGG"
R = "TCCTCCGCTTATTGATATGC"


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def gzip_text(path, text):
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as stream:
            stream.write(text.encode())


def describe(path, relative=None):
    return dict(path=relative or path.name, sha256=freezer.sha(path), bytes=path.stat().st_size)


def update_artifact_manifest(directory, filename):
    path = directory / filename
    state = freezer.load_json(path)
    for key, item in state["artifacts"].items():
        file = directory / item.get("path", key)
        item.update(sha256=freezer.sha(file), bytes=file.stat().st_size)
    json_write(path, state)
    return state


def mutate_csv(path, change):
    rows = scanner.read_csv(path)
    change(rows)
    scanner.write_gzip_csv(path, list(rows[0]), rows)


class Fixture:
    def __init__(self, directory, seed=1):
        self.root = Path(directory)
        self.panel, self.annotations, self.evidence, self.alignments = [self.root / name for name in ("panel", "annotations", "evidence", "alignments")]
        self.destination = self.root / "pinned"
        raw_dir = self.root / "raw"
        raw_dir.mkdir(parents=True)
        entries = []
        for number in range(3):
            rng = random.Random(seed * 100 + number)
            analysis = "".join(rng.choice("ACGT") for _ in range(260))
            if number < 2:
                analysis = analysis[:29] + F + analysis[29+len(F):228] + scanner.reverse_complement(R) + analysis[248:]
                # S/W ambiguity must be correctly reverse-complemented, not
                # silently changed by the upstream ITSx1.1.3 complement quirk.
                analysis = analysis[:100] + "SW" + analysis[102:]
            sequence = scanner.reverse_complement(analysis) if number == 1 else analysis
            order = ("OrderA", "OrderB", "Fungi_ord_Incertae_sedis")[number]
            class_name = "ClassA" if number < 2 else "unidentified"
            header = (f"Genus{number}_species|SYNTHETIC_{seed}_{number}|SH_SYNTHETIC_{seed}_{number}|refs|"
                      f"k__Fungi;p__PhylumA;c__{class_name};o__{order};f__Family{number};g__Genus{number};s__Genus{number}_species")
            entries.append((header, sequence))
        standard, developer, archive = [raw_dir / name for name in ("standard.fasta", "developer.fasta", "archive.tgz")]
        text = "".join(f">{name}\n{seq}\n" for name, seq in entries)
        standard.write_text(text); developer.write_text(text)
        with tarfile.open(archive, "w:gz") as package:
            package.add(standard, arcname=standard.name); package.add(developer, arcname=developer.name)
        expected = dict(archive=freezer.sha(archive), standard=freezer.sha(standard), developer=freezer.sha(developer))
        source_builder.build_panel(archive, standard, developer, self.panel, expected=expected, expected_records=3,
                                   source_doi="synthetic:fixture", source_version="Synthetic fixture; no biological coverage", source_url="https://example.invalid/synthetic")
        self.metadata = scanner.read_csv(self.panel / "selected_metadata.csv.gz")
        self.by_id = {row["record_id"]: row for row in self.metadata}
        self.source = dict(scanner.read_fasta(self.panel / "selected_raw.fasta.gz"))
        self.annotations.mkdir()
        orientations, regions, analyses = [], [], {}
        for row in self.metadata:
            key = row["record_id"]; number = int(row["accession"].split("_")[-1])
            strand = -1 if number == 1 else 1
            analysis = scanner.reverse_complement(self.source[key]) if strand == -1 else self.source[key]
            analyses[key] = analysis
            detected = number < 2
            state = "profile_detected" if detected else "not_detected"
            orientations.append(dict(record_id=key, source_strand=strand, orientation="profile_forward" if detected else "source_unverified",
                annotation_status=state, annotation_note="Explicitly synthetic fixture", source_sequence_sha256=freezer.bases_sha(self.source[key]),
                analysis_sequence_sha256=freezer.bases_sha(analysis)))
            for region, interval in zip(freezer.REGIONS, ((1, 50), (51, 90), (91, 130), (131, 210), (211, 260))):
                start, end = interval if detected else (None, None)
                original = (261-end, 261-start) if detected and strand == -1 else (start, end)
                regions.append(dict(record_id=key, region=region, start=start, end=end,
                    boundary_status="flank_detected" if detected and region in ("SSU", "LSU") else "complete" if detected else "not_detected",
                    source_strand=strand, annotation_source="Synthetic fixture only", coordinate_basis="1-based inclusive positions in oriented source sequence",
                    reason="Synthetic", original_start=original[0], original_end=original[1], itsx_reported_status="Synthetic",
                    annotation_status=state, annotation_note="Synthetic"))
        self.analyses, self.regions = analyses, regions
        gzip_text(self.annotations / "oriented_raw.fasta.gz", "".join(f">{key}\n{seq}\n" for key, seq in analyses.items()))
        scanner.write_gzip_csv(self.annotations / "orientations.csv.gz", list(orientations[0]), orientations)
        scanner.write_gzip_csv(self.annotations / "regions.csv.gz", list(regions[0]), regions)
        for region in freezer.REGIONS:
            records = {row["record_id"]: analyses[row["record_id"]][row["start"]-1:row["end"]]
                       for row in regions if row["region"] == region and row["start"] is not None}
            fasta = "".join(f">{key}\n{seq}\n" for key, seq in records.items())
            gzip_text(self.annotations / ("regions." + region.replace(".", "_") + ".fasta.gz"), fasta)
            if region in ("ITS1", "5.8S", "ITS2"):
                gzip_text(self.annotations / (region.replace(".", "_") + ".fasta.gz"), fasta)
        for filename in ("code_region_parser.py", "code_parser_scanner_dependency.py", "code_profile_run_annotator.py", "code_profile_run_scanner.py"):
            (self.annotations / filename).write_text("# Synthetic captured-code fixture; no profile execution.\n")
        tool = dict(itsx_version="1.1.3", hmmer_version="3.4", files={name: dict(path="/synthetic/"+name, sha256="a"*64)
                   for name in ("itsx_program", "itsx_wrapper", "hmmsearch", "fungal_profile")})
        json_write(self.annotations / "tool_identity_manifest.json", tool)
        semantics = dict(profile="F", preserve=True, complement=True, cpu=1)
        run = dict(schema_version="1.0.0", capture_mode="before_profile_execution", work_directory="/synthetic/work",
                   input_sha256=freezer.sha(self.panel / "selected_raw.fasta.gz"), n_input=3,
                   staged_input_sha256=freezer.bases_sha("".join(f">{key}\n{seq}\n" for key, seq in self.source.items())),
                   tool_manifest_sha256=freezer.sha(self.annotations / "tool_identity_manifest.json"), semantic_parameters=semantics,
                   executed_source={name+"_sha256": freezer.sha(self.annotations / ("code_profile_run_"+name+".py")) for name in ("annotator", "scanner")})
        json_write(self.annotations / "profile_run_identity.json", run)
        identity = dict(input_sha256=run["input_sha256"], metadata_sha256=freezer.sha(self.panel / "selected_metadata.csv.gz"),
                        script_sha256=freezer.sha(self.annotations / "code_region_parser.py"),
                        helper_sha256=freezer.sha(self.annotations / "code_parser_scanner_dependency.py"),
                        tool_manifest_sha256=run["tool_manifest_sha256"], semantic_parameters=semantics)
        manifest = dict(identity=identity, input_sha256=identity["input_sha256"], script_sha256=identity["script_sha256"],
                        helper_sha256=identity["helper_sha256"], tool_manifest_sha256=identity["tool_manifest_sha256"],
                        profile_run_manifest_sha256=freezer.sha(self.annotations / "profile_run_identity.json"),
                        parser_rule_version="itsx_regions_verified_oriented_slices_v2", n_input=3, n_detected=2, n_reverse=1, n_chimeric=0,
                        summary_proof=dict(real_detected_all_batches=2), region_counts=dict(Counter(r["region"]+":"+r["boundary_status"] for r in regions)),
                        artifacts={path.name: describe(path) for path in self.annotations.iterdir()})
        json_write(self.annotations / "annotation_manifest.json", manifest)
        catalog = [dict(catalog_id="FORWARD", primer_name="SYNTHETIC_F", sequence=F, direction="forward"),
                   dict(catalog_id="REVERSE", primer_name="SYNTHETIC_R", sequence=R, direction="reverse")]
        pairs = [dict(pair_id="PAIR", pair_label="Synthetic", forward_primer="SYNTHETIC_F", forward_sequence=F,
                      reverse_primer="SYNTHETIC_R", reverse_sequence=R, evidence_status="synthetic_fixture")]
        scanner.write_gzip_csv(self.root / "catalog.csv.gz", list(catalog[0]), catalog)
        scanner.write_gzip_csv(self.root / "pairs.csv.gz", list(pairs[0]), pairs)
        with contextlib.redirect_stdout(io.StringIO()):
            scanner.build(self.panel, self.root / "catalog.csv.gz", self.root / "pairs.csv.gz", self.evidence,
                          self.annotations / "oriented_raw.fasta.gz")
        self.alignments.mkdir(); (self.alignments / "alignments").mkdir()
        complete = defaultdict(list)
        for row in regions:
            if row["region"] in ("ITS1", "5.8S", "ITS2") and row["boundary_status"] == "complete":
                meta = self.by_id[row["record_id"]]
                for rank in ("class", "order"):
                    if meta[freezer.RANK_COLUMNS[rank]+"_status"] == "resolved":
                        complete[(rank, meta[rank+"_group_id"], row["region"])].append(row)
        index = []
        for number, ((rank, group_id, region), rows) in enumerate(sorted(complete.items())):
            path = self.alignments / "alignments" / f"synthetic_{number}.fasta.gz"
            records = {row["record_id"]: analyses[row["record_id"]][row["start"]-1:row["end"]] for row in rows}
            if len(records) > 1:
                for i, key in enumerate(records):
                    at = 5 if i == 0 else len(records[key])
                    records[key] = records[key][:at] + "-" + records[key][at:]
            gzip_text(path, "".join(f">{key}\n{seq}\n" for key, seq in records.items()))
            index.append(dict(rank=rank, group_id=group_id, region=region, file=str(path.relative_to(self.alignments)),
                n_source_regions=len(rows), n_input=len(rows), n_columns=len(next(iter(records.values()))),
                n_gap_cells=sum(value.count("-") for value in records.values()),
                n_ambiguous_cells=sum(sum(c not in "ACGT-" for c in value) for value in records.values()), alignment_scope="Synthetic fixture"))
        scanner.write_gzip_csv(self.alignments / "alignment_index.csv.gz", list(index[0]), index)
        state = dict(identity=dict(input_sha256=freezer.sha(self.panel / "selected_metadata.csv.gz"),
            regions_sha256=freezer.sha(self.annotations / "regions.csv.gz"), oriented_sha256=freezer.sha(self.annotations / "oriented_raw.fasta.gz"),
            script_sha256=freezer.sha(ROOT / "scripts/build_unite_regional_alignments.py"), limit=32, mafft_version="Synthetic fixture; no MAFFT execution"),
            n_alignments=len(index), n_alignment_records=sum(r["n_input"] for r in index),
            artifacts={str(path.relative_to(self.alignments)): describe(path, str(path.relative_to(self.alignments)))
                       for path in self.alignments.rglob("*") if path.is_file()})
        json_write(self.alignments / "alignment_manifest.json", state)

    def freeze(self, **kwargs):
        return freezer.freeze_snapshot(self.panel, self.annotations, self.evidence, self.alignments, self.destination, **kwargs)

    def rehash_annotation(self):
        return update_artifact_manifest(self.annotations, "annotation_manifest.json")

    def rehash_alignment(self):
        return update_artifact_manifest(self.alignments, "alignment_manifest.json")

    def rehash_search(self):
        path = self.evidence / "site_search_manifest.json"; state = freezer.load_json(path)
        for file, field in (("catalog.csv.gz", "catalog_artifact_sha256"), ("pairs.csv.gz", "pair_artifact_sha256"), ("primer_sites.csv.gz", "site_artifact_sha256")):
            state[field] = freezer.sha(self.evidence / file)
        json_write(path, state)
        return state


class SnapshotTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = Fixture(self.temporary.name)

    def tearDown(self): self.temporary.cleanup()

    def test_good_freeze_preserves_exact_sequence_bytes_source_hashes_and_unknown_evidence(self):
        f = self.fixture; result = f.freeze()
        root = f.destination / "versions" / result["id"]
        self.assertEqual(result["n_source"], 3)
        self.assertEqual(result["n_selected"], 3)
        self.assertTrue(result["id"].startswith("u1-"))
        self.assertEqual((root / "sequences.fasta.gz").read_bytes(), (f.annotations / "oriented_raw.fasta.gz").read_bytes())
        frozen = scanner.read_csv(root / "metadata.csv.gz")
        for row in frozen:
            self.assertEqual(row["sequence_sha256"], f.by_id[row["record_id"]]["sequence_sha256"])
            self.assertEqual(row["analysis_sequence_sha256"], freezer.bases_sha(f.analyses[row["record_id"]]))
            self.assertEqual(row["organism"], row["species"].replace("_", " "))
        self.assertEqual(sum(r["source_strand"] == "-1" for r in frozen), 1)
        self.assertEqual(sum(r["orientation"] == "source_unverified" for r in frozen), 1)
        self.assertIn("source_source_census_csv_gz", result["artifacts"])
        freezer.validate_frozen(root, result)
        self.assertEqual(freezer.load_json(f.destination / "latest.json")["sha256"], freezer.sha(root / "manifest.json"))

    def test_replay_never_rewrites_any_snapshot_or_pointer(self):
        f = self.fixture; first = f.freeze()
        before = {str(path.relative_to(f.destination)): (path.read_bytes(), path.stat().st_mtime_ns)
                  for path in f.destination.rglob("*") if path.is_file()}
        self.assertEqual(f.freeze()["id"], first["id"])
        after = {str(path.relative_to(f.destination)): (path.read_bytes(), path.stat().st_mtime_ns)
                 for path in f.destination.rglob("*") if path.is_file()}
        self.assertEqual(before, after)

    def test_ephemeral_run_paths_are_frozen_provenance_but_not_version_identity(self):
        f = self.fixture; first = f.freeze()
        run_path = f.annotations / "profile_run_identity.json"; run = freezer.load_json(run_path)
        run["work_directory"] = "/different/ephemeral/work"
        run["captured_at_utc"] = "different timestamp"
        json_write(run_path, run)
        state = f.rehash_annotation(); state["profile_run_manifest_sha256"] = freezer.sha(run_path)
        json_write(f.annotations / "annotation_manifest.json", state)
        self.assertEqual(f.freeze()["id"], first["id"])

    def test_new_evidence_appends_and_historical_bytes_survive_without_rewinding_latest(self):
        f = self.fixture; first = f.freeze()
        historical = f.destination / "versions" / first["id"]
        old = {str(path.relative_to(historical)): path.read_bytes() for path in historical.rglob("*") if path.is_file()}
        second = Fixture(Path(self.temporary.name) / "second", seed=2)
        second.destination = f.destination
        newer = second.freeze()
        self.assertTrue(newer["id"].startswith("u2-"))
        self.assertEqual(old, {str(path.relative_to(historical)): path.read_bytes() for path in historical.rglob("*") if path.is_file()})
        self.assertEqual(f.freeze()["id"], first["id"])
        self.assertEqual(freezer.load_json(f.destination / "latest.json")["id"], newer["id"])
        self.assertEqual(len(freezer.load_json(f.destination / "index.json")["versions"]), 2)

    def test_incomplete_or_foreign_search_proof_cannot_impute_zero_hits(self):
        f = self.fixture; path = f.evidence / "site_search_manifest.json"; state = freezer.load_json(path)
        state["complete"] = False; json_write(path, state)
        with self.assertRaisesRegex(ValueError, "Exhaustive site-search"):
            f.freeze()
        state["complete"] = True; state["sequence_sha256"] = "0"*64; json_write(path, state)
        with self.assertRaisesRegex(ValueError, "Exhaustive site-search"):
            f.freeze()
        self.assertFalse(f.destination.exists())

    def test_sparse_coordinates_duplicate_keys_and_count_loss_are_rejected(self):
        f = self.fixture; path = f.evidence / "primer_sites.csv.gz"
        original = path.read_bytes()
        mutate_csv(path, lambda rows: rows[0].update(start=str(int(rows[0]["start"])+1)))
        f.rehash_search()
        with self.assertRaisesRegex(ValueError, "source bases"):
            f.freeze()
        path.write_bytes(original)
        mutate_csv(path, lambda rows: rows.append(dict(rows[0])))
        state = f.rehash_search(); state["n_sparse_rows"] += 1; state["n_no_exact_hits"] -= 1
        json_write(f.evidence / "site_search_manifest.json", state)
        with self.assertRaisesRegex(ValueError, "Sparse site key"):
            f.freeze()

    def test_pair_counts_and_lengths_are_reconciled_from_exact_sites(self):
        f = self.fixture
        mutate_csv(f.evidence / "pairs.csv.gz", lambda rows: rows[0].update(amplicon_min_bp="999"))
        f.rehash_search()
        with self.assertRaisesRegex(ValueError, "length summary"):
            f.freeze()

    def test_taxonomy_group_alias_cannot_merge_unrelated_lineages(self):
        f = self.fixture
        mutate_csv(f.panel / "selected_metadata.csv.gz", lambda rows: rows[0].update(order_group_id=rows[1]["order_group_id"]))
        update_artifact_manifest(f.panel, "source_manifest.json")
        with self.assertRaisesRegex(ValueError, "Hierarchy identity"):
            f.freeze()

    def test_incorrect_reverse_complement_or_original_coordinate_projection_is_rejected(self):
        f = self.fixture
        mutate_csv(f.annotations / "regions.csv.gz", lambda rows: next(r for r in rows if r["source_strand"] == "-1").update(original_start="2"))
        f.rehash_annotation()
        with self.assertRaisesRegex(ValueError, "projection differs"):
            f.freeze()

    def test_retained_lsu_flank_cannot_be_called_a_complete_gene(self):
        f = self.fixture
        mutate_csv(f.annotations / "regions.csv.gz", lambda rows: next(r for r in rows if r["region"] == "LSU" and r["start"]).update(boundary_status="complete"))
        state = f.rehash_annotation()
        state["region_counts"] = dict(Counter(r["region"]+":"+r["boundary_status"] for r in scanner.read_csv(f.annotations / "regions.csv.gz")))
        json_write(f.annotations / "annotation_manifest.json", state)
        with self.assertRaisesRegex(ValueError, "complete gene"):
            f.freeze()

    def test_extracted_complete_region_must_equal_its_oriented_raw_slice(self):
        f = self.fixture; path = f.annotations / "ITS1.fasta.gz"
        rows = list(scanner.read_fasta(path)); key, seq = rows[0]
        rows[0] = (key, ("A" if seq[0] != "A" else "C") + seq[1:])
        gzip_text(path, "".join(f">{name}\n{value}\n" for name, value in rows)); f.rehash_annotation()
        with self.assertRaisesRegex(ValueError, "source slice"):
            f.freeze()

    def test_alignment_changed_base_wrong_owner_or_width_is_rejected_even_with_new_hash(self):
        f = self.fixture; index = scanner.read_csv(f.alignments / "alignment_index.csv.gz")
        row = next(r for r in index if int(r["n_input"]) == 2)
        path = f.alignments / row["file"]; original = path.read_bytes(); records = list(scanner.read_fasta(path))
        key, seq = records[0]
        records[0] = (key, ("A" if seq[0] != "A" else "C") + seq[1:])
        gzip_text(path, "".join(f">{name}\n{value}\n" for name, value in records)); f.rehash_alignment()
        with self.assertRaisesRegex(ValueError, "source base retention"):
            f.freeze()
        path.write_bytes(original); records = list(scanner.read_fasta(path)); records[0] = ("FOREIGN", records[0][1])
        gzip_text(path, "".join(f">{name}\n{value}\n" for name, value in records)); f.rehash_alignment()
        with self.assertRaisesRegex(ValueError, "membership/width"):
            f.freeze()

    def test_alignment_source_denominator_is_not_the_display_exemplar_count(self):
        f = self.fixture
        mutate_csv(f.alignments / "alignment_index.csv.gz", lambda rows: rows[0].update(n_source_regions="99"))
        f.rehash_alignment()
        with self.assertRaisesRegex(ValueError, "index counts"):
            f.freeze()

    def test_budget_and_interrupted_copy_leave_no_version_or_registry(self):
        f = self.fixture
        with self.assertRaisesRegex(ValueError, "byte budget"):
            f.freeze(max_bytes=1)
        self.assertFalse((f.destination / "index.json").exists())
        self.assertEqual(list((f.destination / "versions").iterdir()), [])
        with patch.object(freezer.shutil, "copyfile", side_effect=OSError("synthetic interrupted copy")):
            with self.assertRaises(OSError): f.freeze()
        self.assertFalse((f.destination / ".freeze.lock").exists())
        self.assertEqual(list((f.destination / "versions").iterdir()), [])

    def test_corrupt_history_or_path_escape_fails_before_promotion(self):
        f = self.fixture; result = f.freeze()
        pointer_before = (f.destination / "latest.json").read_bytes()
        historical = f.destination / "versions" / result["id"]
        (historical / "sequences.fasta.gz").write_bytes(b"changed historical bytes")
        with self.assertRaisesRegex(ValueError, "checksum differs"):
            f.freeze()
        self.assertEqual((f.destination / "latest.json").read_bytes(), pointer_before)
        for relative in ("../escape", "/absolute", "alignments/../../escape", "bad\\path"):
            with self.assertRaises(ValueError): freezer.safe_path(f.destination, relative)

    def test_completed_directory_can_recover_interrupted_registration(self):
        f = self.fixture; first = f.freeze()
        (f.destination / "index.json").unlink(); (f.destination / "latest.json").unlink()
        self.assertEqual(f.freeze()["id"], first["id"])
        self.assertEqual(len(freezer.load_json(f.destination / "index.json")["versions"]), 1)
        self.assertEqual(freezer.load_json(f.destination / "latest.json")["id"], first["id"])

    def test_index_written_before_latest_is_repaired_without_rewinding_history(self):
        f = self.fixture; first = f.freeze()
        second = Fixture(Path(self.temporary.name) / "second", seed=2); second.destination = f.destination
        newer = second.freeze()
        entries = freezer.load_json(f.destination / "index.json")["versions"]
        json_write(f.destination / "latest.json", entries[0])
        self.assertEqual(f.freeze()["id"], first["id"])
        self.assertEqual(freezer.load_json(f.destination / "latest.json")["id"], newer["id"])
        (f.destination / "latest.json").unlink()
        self.assertEqual(second.freeze()["id"], newer["id"])
        self.assertEqual(freezer.load_json(f.destination / "latest.json")["id"], newer["id"])

    def test_corrupt_latest_pointer_is_not_silently_repaired_or_promoted(self):
        f = self.fixture; f.freeze()
        path = f.destination / "latest.json"; pointer = freezer.load_json(path)
        pointer["sha256"] = "0"*64; json_write(path, pointer)
        before = path.read_bytes()
        with self.assertRaisesRegex(ValueError, "latest pointer/proof"):
            f.freeze()
        self.assertEqual(path.read_bytes(), before)


if __name__ == "__main__": unittest.main()
