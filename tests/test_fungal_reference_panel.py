"""Scientific boundary and retention checks for the bounded panel builder."""
import importlib.util
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

SPEC = importlib.util.spec_from_file_location("panel", Path(__file__).resolve().parents[1] / "scripts/build_fungal_reference_panel.py")
panel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(panel)


class BoundaryTests(unittest.TestCase):
    def feature(self, start, end, note, key="misc_RNA", qualifier="note", location=None):
        return (f"<GBFeature><GBFeature_key>{key}</GBFeature_key>"
                f"<GBFeature_location>{escape(location or str(start) + '..' + str(end))}</GBFeature_location>"
                f"<GBFeature_intervals><GBInterval><GBInterval_from>{start}</GBInterval_from>"
                f"<GBInterval_to>{end}</GBInterval_to></GBInterval></GBFeature_intervals>"
                f"<GBFeature_quals><GBQualifier><GBQualifier_name>{qualifier}</GBQualifier_name>"
                f"<GBQualifier_value>{note}</GBQualifier_value></GBQualifier></GBFeature_quals></GBFeature>")

    def record(self, *features):
        return ET.fromstring("<GBSeq><GBSeq_feature-table>" + "".join(features) + "</GBSeq_feature-table></GBSeq>")

    def test_entire_operon_is_not_a_58s_boundary(self):
        record = self.record(self.feature(1, 5000, "contains small subunit ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and large subunit ribosomal RNA"))
        self.assertEqual(panel.annotate(record, "X.1"), [])

    def test_explicit_genes_support_intervening_spacers(self):
        record = self.record(self.feature(1, 1800, "18S ribosomal RNA", "rRNA", "product"),
                             self.feature(2101, 2260, "5.8S ribosomal RNA", "rRNA", "product"),
                             self.feature(2501, 5000, "28S ribosomal RNA", "rRNA", "product"))
        rows = {x["landmark_id"]: x for x in panel.annotate(record, "X.1")}
        self.assertEqual((rows["ITS1"]["start"], rows["ITS1"]["end"]), (1801, 2100))
        self.assertEqual((rows["ITS2"]["start"], rows["ITS2"]["end"]), (2261, 2500))
        self.assertEqual(rows["5_8S"]["end"] - rows["5_8S"]["start"] + 1, 160)

    def test_complement_and_compound_features_are_not_force_placed(self):
        record = self.record(self.feature(1, 100, "18S ribosomal RNA", "rRNA", "product"))
        interval = record.find(".//GBInterval")
        ET.SubElement(interval, "GBInterval_iscomp", {"value": "true"})
        self.assertEqual(panel.annotate(record, "X.1"), [])
        compound = self.record(self.feature(1, 100, "18S ribosomal RNA", "rRNA", "product"))
        intervals = compound.find(".//GBFeature_intervals")
        second = ET.SubElement(intervals, "GBInterval")
        ET.SubElement(second, "GBInterval_from").text = "200"
        ET.SubElement(second, "GBInterval_to").text = "300"
        self.assertEqual(panel.annotate(compound, "X.1"), [])

    def test_fuzzy_gene_boundaries_do_not_derive_exact_spacers(self):
        record = self.record(self.feature(1, 1800, "18S ribosomal RNA", "rRNA", "product", "<1..>1800"),
                             self.feature(2101, 2260, "5.8S ribosomal RNA", "rRNA", "product"),
                             self.feature(2501, 5000, "28S ribosomal RNA", "rRNA", "product"))
        rows, note = panel.inspect_annotations(record, "X.1")
        ids = {x["landmark_id"] for x in rows}
        self.assertNotIn("18S", ids)
        self.assertNotIn("ITS1", ids)
        self.assertIn("ITS2", ids)
        self.assertIn("fuzzy", note)

    def test_fuzzy_explicit_spacer_is_not_filled_from_adjacent_genes(self):
        record = self.record(self.feature(1, 1800, "18S ribosomal RNA", "rRNA", "product"),
                             self.feature(1801, 2100, "internal transcribed spacer 1", location="<1801..2100"),
                             self.feature(2101, 2260, "5.8S ribosomal RNA", "rRNA", "product"))
        self.assertNotIn("ITS1", {x["landmark_id"] for x in panel.annotate(record, "X.1")})

    def test_repeated_gene_labels_are_withheld_instead_of_last_wins(self):
        record = self.record(self.feature(1, 1800, "18S ribosomal RNA", "rRNA", "product"),
                             self.feature(2101, 2260, "5.8S ribosomal RNA", "rRNA", "product"),
                             self.feature(2501, 3000, "28S ribosomal RNA", "rRNA", "product"),
                             self.feature(3501, 4000, "28S ribosomal RNA", "rRNA", "product"))
        rows, note = panel.inspect_annotations(record, "X.1")
        ids = {x["landmark_id"] for x in rows}
        self.assertNotIn("28S", ids)
        self.assertNotIn("ITS2", ids)
        self.assertIn("ITS1", ids)
        self.assertIn("repeated gene labels", note)


class AcquisitionTests(unittest.TestCase):
    def test_legacy_search_migration_requires_query_limit_and_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            legacy = cache / "search_Fungi.json"
            legacy.write_text(json.dumps({"esearchresult": {"retmax": "4", "count": "100", "idlist": ["A.1", "B.1", "C.1", "D.1"]}}))
            os.utime(legacy, (946684800, 946684800))
            provenance = {"searches": [{"clade": "Fungi", "query": "query", "sha256": panel.sha(legacy),
                                          "selected": ["A.1", "B.1", "C.1", "D.1"]}]}
            migrated = panel.search_cache_path(cache, "Fungi", "query", 4, provenance)
            self.assertTrue(migrated.exists())
            self.assertEqual(panel.sha(migrated), panel.sha(legacy))
            self.assertEqual(panel.acquisition_time(migrated), "2000-01-01T00:00:00Z")
            larger = panel.search_cache_path(cache, "Fungi", "query", 20, provenance)
            other_query = panel.search_cache_path(cache, "Fungi", "different query", 4, provenance)
            self.assertFalse(larger.exists())
            self.assertFalse(other_query.exists())
            self.assertNotEqual(migrated, larger)

    def test_same_count_wrong_accession_versions_are_rejected(self):
        def record(accession):
            return ET.fromstring(f"<GBSeq><GBSeq_accession-version>{accession}</GBSeq_accession-version></GBSeq>")
        panel.validate_accessions([record("A.1"), record("B.2")], ["B.2", "A.1"])
        with self.assertRaises(RuntimeError):
            panel.validate_accessions([record("A.2"), record("B.2")], ["A.1", "B.2"])
        with self.assertRaises(RuntimeError):
            panel.validate_accessions([record("A.1"), record("A.1")], ["A.1", "B.2"])

    def test_transport_exception_is_retried(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "raw.json"
            calls = []
            def run(*args, **kwargs):
                calls.append(1)
                if len(calls) == 1:
                    raise subprocess.CalledProcessError(1, "client")
                destination.write_text("{}")
                return subprocess.CompletedProcess("client", 0, stdout='{"ok":true}')
            with patch.object(panel.subprocess, "run", side_effect=run), patch.object(panel.time, "sleep"):
                self.assertEqual(panel.request("esearch", {}, destination), destination)
            self.assertEqual(len(calls), 2)


class RetentionTests(unittest.TestCase):
    def write_previous(self, directory, sequence):
        (directory / "panel.fasta").write_text(f">X.1\n{sequence}\n")
        row = {"accession": "X.1", "class": "OriginalClass", "phylum": "OriginalPhylum",
               "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
               "taxonomy_source_sha256": "a" * 64}
        panel.write_csv(directory / "panel_manifest.csv", [row], list(row))
        landmark = {"accession": "X.1", "landmark_id": "18S", "label": "18S",
                    "start": 1, "end": 1800, "source_note": "preceding annotation"}
        panel.write_csv(directory / "landmarks.csv", [landmark], list(landmark))
        (directory / "run_manifest.json").write_text(json.dumps({"taxonomy_snapshot_sha256": "a" * 64}))
        return row

    def test_restored_panel_appends_new_template_preserving_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            previous_sequence = "ACGT" * 750
            added_sequence = "ACGA" * 800
            original = self.write_previous(directory, previous_sequence)
            old_alignment_input_hash = panel.sha(directory / "panel.fasta")
            sequences = {"X.1": previous_sequence, "Y.1": added_sequence}
            metadata = [{"accession": "X.1", "class": "ReassignedClass", "taxonomy_source_sha256": "b" * 64},
                        {"accession": "Y.1", "class": "AddedClass", "taxonomy_source_sha256": "b" * 64}]
            landmarks = [{"accession": "X.1", "source_note": "new conflicting annotation"}]
            retained, lineage_hash = panel.retain_previous_panel(directory, sequences, metadata, landmarks)
            self.assertEqual(retained, ["X.1"])
            self.assertEqual(lineage_hash, panel.sha(directory / "run_manifest.json"))
            self.assertEqual(next(row for row in metadata if row["accession"] == "X.1"), original)
            self.assertEqual(hashlib.sha256(sequences["X.1"].encode()).hexdigest(), original["sequence_sha256"])
            self.assertEqual(landmarks[0]["source_note"], "preceding annotation")
            new_input = directory / "grown_alignment_input.fasta"
            new_input.write_text("".join(f">{accession}\n{sequence}\n" for accession, sequence in sequences.items()))
            self.assertNotEqual(panel.sha(new_input), old_alignment_input_hash)
            self.assertEqual(set(panel.read_fasta(new_input)), {"X.1", "Y.1"})

    def test_changed_sequence_for_same_accession_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_previous(directory, "ACGT" * 750)
            with self.assertRaises(RuntimeError):
                panel.retain_previous_panel(directory, {"X.1": "TGCA" * 750}, [], [])


if __name__ == "__main__":
    unittest.main()
