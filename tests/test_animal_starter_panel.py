#!/usr/bin/env python3
"""Validate animal-group identity, archived source slices and bird oligos."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "data/pinned/ANIMAL_STARTER/v1"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class AnimalSourceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.templates = read_csv(PANEL / "templates.csv")
        cls.manifest = json.loads((PANEL / "manifest.json").read_text())
        cls.records = {}
        for source in cls.manifest["sources"]:
            for record in ET.fromstring(gzip.decompress((PANEL / source["file"]).read_bytes())).findall("./GBSeq"):
                cls.records[record.findtext("GBSeq_accession-version")] = record

    def test_each_requested_group_marker_has_its_own_source(self):
        expected = {("ARTHROPODS", "18S"), ("ARTHROPODS", "28S"), ("ARTHROPODS", "16S_MT"),
                    ("BIRDS", "12S_MT"), ("BIRDS", "COI"), ("BIRDS", "18S"), ("BIRDS", "28S"),
                    ("ZOOPLANKTON", "COI"), ("ZOOPLANKTON", "18S"), ("ZOOPLANKTON", "28S"), ("ZOOPLANKTON", "16S_MT")}
        self.assertEqual({(t["group_id"], t["marker_id"]) for t in self.templates}, expected)
        self.assertEqual(len(self.templates), 11)
        self.assertEqual(len({t["template_id"] for t in self.templates}), 11)
        for template in self.templates:
            self.assertEqual(template["target_group"], template["group_id"])
            required = "Arthropoda" if template["group_id"] in ("ARTHROPODS", "ZOOPLANKTON") else "Aves"
            self.assertIn(required, template["source_taxonomy"].split("; "))
            self.assertNotIn("Danio", template["organism"])
            self.assertNotIn("FISH", template["template_id"])

    def test_source_hashes_and_frozen_builder(self):
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

    def test_exact_feature_and_sequence_contract(self):
        for template in self.templates:
            record = self.records[template["accession_version"]]
            self.assertEqual(record.findtext("GBSeq_organism"), template["organism"])
            self.assertEqual(record.findtext("GBSeq_taxonomy"), template["source_taxonomy"])
            source = record.findtext("GBSeq_sequence").upper()
            start, end = int(template["source_start"]), int(template["source_end"])
            expected = source[start - 1:end]
            if template["source_strand"] == "-":
                expected = expected.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
            self.assertEqual(template["sequence"], expected)
            self.assertEqual(len(template["sequence"]), end - start + 1)
            self.assertEqual(len(template["sequence"]), int(template["sequence_length"]))
            self.assertEqual(digest(template["sequence"].encode("ascii")), template["sequence_sha256"])
            self.assertEqual(digest(source.encode("ascii")), template["source_sequence_sha256"])
            matches = [f for f in record.findall("./GBSeq_feature-table/GBFeature")
                       if f.findtext("GBFeature_location") == template["source_feature_location"]
                       and f.findtext("GBFeature_key") in ("CDS", "rRNA")]
            self.assertTrue(matches)

    def test_discontinuous_and_partial_drosophila_repeats_are_withheld(self):
        drosophila_18s = [t for t in self.templates if t["accession_version"] == "M21017.1"]
        self.assertEqual(len(drosophila_18s), 1)
        self.assertEqual(drosophila_18s[0]["source_feature_location"], "1..1995")
        self.assertEqual(drosophila_18s[0]["marker_id"], "18S")
        self.assertEqual(len(self.manifest["withheld_features"]), 2)
        for template in self.templates:
            self.assertNotIn("join", template["source_feature_location"])
            is_partial = "<" in template["source_feature_location"] or ">" in template["source_feature_location"]
            self.assertEqual(template["partial"], str(is_partial).lower())
            if is_partial:
                self.assertIn("missing primer site", template["scope_note"])

    def test_predicted_source_remains_explicit(self):
        models = [t for t in self.templates if t["source_status"] == "PREDICTED RefSeq MODEL"]
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]["accession_version"], "XR_005840241.1")
        self.assertIn("not an experimentally validated transcript", models[0]["scope_note"])
        self.assertEqual(self.manifest["status"], "starter_mapping_panel")
        self.assertIn("not monthly maintained", self.manifest["scope"])

    def test_bird_oligos_and_publication_identity(self):
        primers = [p for p in read_csv(ROOT / "data/catalog/animal_assay_primers.csv")
                   if p["marker_id"] == "COI" and p["target_group"] == "BIRDS"]
        expected = {"BirdF1": "TTCTCCAACCACAAAGACATTGGCAC", "BirdR1": "ACGTGGGAGATAATTCCAAATCCTG",
                    "BirdR2": "ACTACATGTGAGATGATTCCGAATCCAG", "BirdR3": "AGGAGTTTGCTAGTACGATGCC"}
        self.assertEqual({p["primer_name"]: p["sequence"] for p in primers}, expected)
        sources = {s["key"] for s in read_csv(ROOT / "data/catalog/sources.csv")}
        for primer in primers:
            self.assertEqual(primer["marker_id"], "COI")
            self.assertEqual(primer["target_group"], "BIRDS")
            self.assertIn(primer["source_key"], sources)
            self.assertEqual(primer["doi"], "10.1371/journal.pbio.0020312")
            self.assertIn("not evidence of water-eDNA specificity", primer["evidence_note"])
        memberships = {}
        for primer in primers:
            for pair in primer["documented_pair_ids"].split("|"):
                memberships.setdefault(pair, []).append(primer["direction"])
        self.assertEqual(len(memberships), 3)
        for directions in memberships.values():
            self.assertEqual(sorted(directions), ["forward", "reverse"])

    def test_bird_12s_and_arthropod_28s_keep_actual_primary_source_scope(self):
        primers = read_csv(ROOT / "data/catalog/animal_assay_primers.csv")
        self.assertEqual(len(primers), 10)
        self.assertEqual(len({p["component_id"] for p in primers}), 10)
        by_name = {p["primer_name"]: p for p in primers}
        expected = {"MiBird-U-F": "GGGTTGGTAAATCTTGTGCCAGC", "MiBird-U-R": "CATAGTGGGGTATCTAATCCCAGTTTG",
                    "D2F": "AGAGAGAGTTCAAGAGTACGTG", "D3R": "TAGTTCACCATCTTTCGGGTCCC"}
        for name, sequence in expected.items():
            self.assertEqual(by_name[name]["sequence"], sequence)
        for name in ("MiBird-U-F", "MiBird-U-R"):
            self.assertEqual(by_name[name]["target_group"], "BIRDS")
            self.assertEqual(by_name[name]["marker_id"], "12S_MT")
            self.assertEqual(by_name[name]["doi"], "10.1038/s41598-018-22817-5")
            self.assertIn("six random bases are excluded", by_name[name]["evidence_note"])
            self.assertIn("other vertebrates", by_name[name]["target_scope"])
        for name in ("D2F", "D3R"):
            self.assertEqual(by_name[name]["target_group"], "ARTHROPODS")
            self.assertEqual(by_name[name]["marker_id"], "28S")
            self.assertEqual(by_name[name]["doi"], "10.1371/journal.pone.0074837")
            self.assertIn("not universal arthropod recovery", by_name[name]["evidence_note"])
        sources = {s["key"] for s in read_csv(ROOT / "data/catalog/sources.csv")}
        self.assertTrue(all(p["source_key"] in sources for p in primers))
        mitochondrial_16s = [p for p in primers if p["marker_id"] == "16S_MT"]
        self.assertEqual({p["primer_name"]: p["sequence"] for p in mitochondrial_16s},
                         {"16Sar": "CGCCTGTTTATCAAAAACAT", "16Sbr": "CCGGTCTGAACTCAGATCACGT"})
        self.assertEqual({p["direction"] for p in mitochondrial_16s}, {"forward", "reverse"})
        for primer in mitochondrial_16s:
            self.assertEqual(primer["doi"], "10.3897/zookeys.876.37244")
            self.assertIn("not transferred to Daphnia", primer["evidence_note"])

    def test_drosophila_rna_matches_its_contiguous_genomic_copy(self):
        proof = self.manifest["genomic_correspondences"][0]
        template = next(t for t in self.templates if t["accession_version"] == proof["template_accession_version"])
        genomic = self.records[proof["genomic_accession_version"]].findtext("GBSeq_sequence").upper()
        sequence = genomic[proof["genomic_start"] - 1:proof["genomic_end"]]
        self.assertEqual(template["sequence"], sequence)
        self.assertEqual(digest(sequence.encode("ascii")), proof["sequence_sha256"])
        self.assertIn("contiguous annotated rDNA interval", template["scope_note"])

    def test_zooplankton_are_actual_source_supported_daphnia(self):
        rows = [t for t in self.templates if t["group_id"] == "ZOOPLANKTON"]
        self.assertEqual(len(rows), 4)
        self.assertEqual({r["organism"] for r in rows}, {"Daphnia pulex", "Daphnia magna"})
        self.assertEqual({r["marker_id"] for r in rows}, {"COI", "18S", "28S", "16S_MT"})
        sources = {r["key"] for r in read_csv(ROOT / "data/catalog/sources.csv")}
        for row in rows:
            self.assertIn("Daphnia", row["source_taxonomy"].split("; "))
            self.assertIn(row["ecology_source_key"], sources)
            self.assertIn("Freshwater zooplankton", row["ecology_scope"])
            self.assertTrue(row["ecology_source_url"].startswith("https://doi.org/"))
        mitochondrial_16s = next(r for r in rows if r["marker_id"] == "16S_MT")
        self.assertEqual(mitochondrial_16s["source_strand"], "-")
        self.assertEqual(mitochondrial_16s["source_feature_location"], "complement(12366..13738)")
        nuclear_28s = next(r for r in rows if r["marker_id"] == "28S")
        self.assertEqual(nuclear_28s["partial"], "true")


if __name__ == "__main__":
    unittest.main()
