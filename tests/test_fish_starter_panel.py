#!/usr/bin/env python3
"""Source-to-artifact checks for the fish assay catalog and bounded panel."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "data/pinned/FISH_STARTER/v1"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class FishSourceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.primers = read_csv(ROOT / "data/catalog/fish_assay_primers.csv")
        cls.templates = read_csv(PANEL / "templates.csv")
        cls.manifest = json.loads((PANEL / "manifest.json").read_text())
        cls.records = {}
        for source in cls.manifest["sources"]:
            raw = gzip.decompress((PANEL / source["file"]).read_bytes())
            for record in ET.fromstring(raw).findall("./GBSeq"):
                cls.records[record.findtext("GBSeq_accession-version")] = record

    def test_catalog_identity_and_marker_scope(self):
        self.assertEqual(len(self.primers), 22)
        self.assertEqual(len({p["component_id"] for p in self.primers}), 22)
        self.assertEqual({p["marker_id"] for p in self.primers}, {"COI", "12S_MT", "18S", "28S"})
        for primer in self.primers:
            self.assertEqual(primer["target_group"], "FISH")
            self.assertIn(primer["direction"], {"forward", "reverse"})
            self.assertTrue(set(primer["sequence"]) <= set("ACGTRYSWKMBDHVN"))
            self.assertTrue(15 <= len(primer["sequence"]) <= 35)
            self.assertTrue(primer["source_url"].startswith("https://"))
            self.assertTrue(primer["verification_url"].startswith("https://"))

    def test_assay_sources_resolve_and_both_roles_are_retained(self):
        keys = {s["key"] for s in read_csv(ROOT / "data/catalog/sources.csv")}
        memberships = {}
        for primer in self.primers:
            self.assertIn(primer["source_key"], keys)
            for pair in primer["documented_pair_ids"].split("|"):
                memberships.setdefault(pair, []).append(primer)
        self.assertEqual(len(memberships), 12)
        for pair, primers in memberships.items():
            self.assertEqual(len(primers), 2, pair)
            self.assertEqual({p["direction"] for p in primers}, {"forward", "reverse"}, pair)
            self.assertEqual(len({p["marker_id"] for p in primers}), 1, pair)

    def test_18s_oligos_match_unchanged_pr2_snapshot(self):
        pr2 = read_csv(ROOT / "data/catalog/pr2_18s_primers.csv")
        original = {(p["primer_name"], p["sequence"], p["direction"]) for p in pr2}
        for primer in self.primers:
            if primer["marker_id"] == "18S":
                self.assertIn((primer["primer_name"], primer["sequence"], primer["direction"]), original)
                self.assertIn("Broad eukaryotes", primer["target_scope"])
                self.assertIn("not fish-specific", primer["evidence_note"])

    def test_archive_and_artifact_hashes(self):
        for name, metadata in self.manifest["artifacts"].items():
            raw = (PANEL / name).read_bytes()
            self.assertEqual(digest(raw), metadata["sha256"])
            self.assertEqual(len(raw), metadata["bytes"])
        for source in self.manifest["sources"]:
            raw = (PANEL / source["file"]).read_bytes()
            self.assertEqual(digest(raw), source["sha256"])
            self.assertEqual(digest(gzip.decompress(raw)), source["uncompressed_sha256"])
        builder = self.manifest["builder"]
        self.assertEqual(digest((ROOT / builder["path"]).read_bytes()), builder["sha256"])

    def test_all_templates_are_exact_original_source_feature_slices(self):
        self.assertEqual(len(self.templates), 13)
        self.assertEqual(len(self.records), 7)
        self.assertEqual(len({t["template_id"] for t in self.templates}), 13)
        for template in self.templates:
            record = self.records[template["accession_version"]]
            original_sequence = record.findtext("GBSeq_sequence").upper()
            start, end = int(template["source_start"]), int(template["source_end"])
            self.assertEqual(template["sequence"], original_sequence[start - 1:end])
            self.assertEqual(int(template["sequence_length"]), end - start + 1)
            self.assertEqual(digest(template["sequence"].encode("ascii")), template["sequence_sha256"])
            self.assertEqual(digest(original_sequence.encode("ascii")), template["source_sequence_sha256"])
            matching_features = [f for f in record.findall("./GBSeq_feature-table/GBFeature")
                                 if f.findtext("GBFeature_location") == template["source_feature_location"]]
            self.assertTrue(matching_features, template["template_id"])
            expected = {"COI": "cytochrome c oxidase subunit I", "18S": "18S", "28S": "28S"}
            if template["marker_id"] in expected:
                self.assertIn(expected[template["marker_id"]], template["source_feature"])
            elif template["marker_id"] == "16S_MT":
                self.assertTrue("16S" in template["source_feature"] or "l-rRNA" in template["source_feature"])
            else:
                self.assertTrue("12S" in template["source_feature"] or "s-rRNA" in template["source_feature"])

    def test_partial_and_predicted_sources_are_not_silently_upgraded(self):
        for template in self.templates:
            record = self.records[template["accession_version"]]
            definition = record.findtext("GBSeq_definition", "")
            if definition.startswith("PREDICTED:"):
                self.assertIn("MODEL", template["source_status"])
                self.assertIn("Predicted genomic", template["scope_note"])
            if "partial" in definition:
                self.assertEqual(template["partial"], "true")
                self.assertIn("missing primer sites", template["scope_note"])
        self.assertEqual(self.manifest["status"], "starter_mapping_panel")
        self.assertIn("not monthly maintained", self.manifest["scope"])


if __name__ == "__main__":
    unittest.main()
