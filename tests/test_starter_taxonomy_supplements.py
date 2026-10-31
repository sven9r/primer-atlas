#!/usr/bin/env python3
"""Check authoritative accession-to-taxid joins and rank denominators."""

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def record_taxid(record):
    values = []
    for feature in record.findall("./GBSeq_feature-table/GBFeature"):
        if feature.findtext("GBFeature_key") != "source":
            continue
        for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
            value = qualifier.findtext("GBQualifier_value", "")
            if qualifier.findtext("GBQualifier_name") == "db_xref" and value.startswith("taxon:"):
                values.append(value.split(":", 1)[1])
    return values


class StarterTaxonomyChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panels = {}
        for name in ("FISH_STARTER", "ANIMAL_STARTER"):
            directory = ROOT / "data/pinned" / name / "v1"
            manifest = json.loads((directory / "manifest.json").read_text())
            table = read_csv(directory / "taxonomy.csv")
            raw = gzip.decompress((directory / table[0]["raw_taxonomy_file"]).read_bytes())
            taxa = {taxon.findtext("TaxId"): taxon for taxon in ET.fromstring(raw).findall("./Taxon")}
            cls.panels[name] = (directory, manifest, table, taxa)

    def test_every_original_template_has_exactly_one_taxonomy_row(self):
        for name, (directory, manifest, table, _) in self.panels.items():
            originals = read_csv(directory / "templates.csv")
            ids = [row["template_id"] for row in table]
            self.assertEqual(len(ids), len(set(ids)))
            extras = manifest["taxonomy_supplement"]["external_template_ids"]
            self.assertEqual(set(ids), {row["template_id"] for row in originals} | set(extras))
            self.assertEqual(manifest["taxonomy_supplement"]["n_panel_template_rows"], len(originals))
            self.assertEqual(len(table), 13 if name == "FISH_STARTER" else 12)
            for original in originals:
                taxonomy = next(row for row in table if row["template_id"] == original["template_id"])
                for field in ("accession_version", "marker_id", "group_id"):
                    self.assertEqual(taxonomy[field], original[field])
                self.assertEqual(taxonomy["species"], original["organism"])

    def test_taxid_is_read_from_original_accession_not_inferred_from_name(self):
        for _, (directory, _, table, _) in self.panels.items():
            for row in table:
                compressed = (directory / row["source_record_file"]).read_bytes()
                self.assertEqual(digest(compressed), row["source_record_sha256"])
                records = ET.fromstring(gzip.decompress(compressed)).findall("./GBSeq")
                records = [r for r in records if r.findtext("GBSeq_accession-version") == row["accession_version"]]
                self.assertEqual(len(records), 1)
                self.assertEqual(record_taxid(records[0]), [row["taxid"]])
                self.assertEqual(row["source_record_taxid"], row["taxid"])
                self.assertEqual(row["source_record_organism"], records[0].findtext("GBSeq_organism"))

    def test_rank_names_are_exact_authoritative_rank_values(self):
        for _, (_, _, table, taxa) in self.panels.items():
            for row in table:
                taxon = taxa[row["taxid"]]
                nodes = [*taxon.findall("./LineageEx/Taxon"), taxon]
                for rank in RANKS:
                    names = {n.findtext("ScientificName") for n in nodes if n.findtext("Rank") == rank}
                    expected = next(iter(names)) if len(names) == 1 else "Unresolved"
                    self.assertEqual(row[rank], expected)
                self.assertEqual(row["taxonomy_status"], "source_taxid_verified")
                self.assertEqual(row["unresolved_ranks"], "")
                if row["group_id"] == "ARTHROPODS":
                    self.assertEqual(row["phylum"], "Arthropoda")
                elif row["group_id"] == "BIRDS":
                    self.assertEqual(row["class"], "Aves")
                elif row["group_id"] == "ZOOPLANKTON":
                    self.assertEqual(row["genus"], "Daphnia")
                    self.assertEqual(row["family"], "Daphniidae")

    def test_hierarchy_counts_reconcile_to_template_denominator(self):
        expected_family_counts = {
            "FISH_STARTER": {"Danionidae": 6, "Gadidae": 3, "Scyliorhinidae": 3, "Salmonidae": 1},
            "ANIMAL_STARTER": {"Drosophilidae": 4, "Phasianidae": 4, "Daphniidae": 4},
        }
        for name, (_, _, table, _) in self.panels.items():
            for rank in RANKS:
                counts = {label: sum(row[rank] == label for row in table) for label in {row[rank] for row in table}}
                self.assertEqual(sum(counts.values()), len(table))
                if rank == "family":
                    self.assertEqual(counts, expected_family_counts[name])

    def test_taxonomy_artifacts_and_raw_xml_are_checksummed(self):
        for _, (directory, manifest, table, _) in self.panels.items():
            artifact = manifest["artifacts"]["taxonomy.csv"]
            self.assertEqual(digest((directory / "taxonomy.csv").read_bytes()), artifact["sha256"])
            supplement = manifest["taxonomy_supplement"]
            self.assertEqual(supplement["n_rows"], len(table))
            for source in supplement["sources"]:
                raw = (directory / source["file"]).read_bytes()
                self.assertEqual(digest(raw), source["sha256"])
                self.assertEqual(digest(gzip.decompress(raw)), source["uncompressed_sha256"])
            for row in table:
                compressed = (directory / row["raw_taxonomy_file"]).read_bytes()
                self.assertEqual(digest(compressed), row["raw_taxonomy_sha256"])
                self.assertEqual(digest(gzip.decompress(compressed)), row["raw_taxonomy_xml_sha256"])
                self.assertEqual(row["taxonomy_source_sha256"], row["raw_taxonomy_xml_sha256"])
            builder = supplement["builder"]
            self.assertEqual(digest((ROOT / builder["path"]).read_bytes()), builder["sha256"])

    def test_existing_dyakuba_fasta_matches_archived_accession_source(self):
        directory, _, table, _ = self.panels["ANIMAL_STARTER"]
        row = next(r for r in table if r["template_id"] == "NC_001322.1_COX1")
        self.assertEqual(row["species"], "Drosophila yakuba")
        self.assertEqual(row["taxid"], "7245")
        raw = gzip.decompress((directory / row["source_record_file"]).read_bytes())
        record = ET.fromstring(raw).find("./GBSeq")
        source = record.findtext("GBSeq_sequence").upper()
        fasta = "".join(line.strip() for line in (ROOT / "data/reference/coi_reference_NC_001322.1.fasta").read_text().splitlines()
                        if not line.startswith(">"))
        self.assertEqual(fasta, source[1473:3009])

    def test_missing_or_conflicting_source_ranks_remain_unresolved(self):
        spec = importlib.util.spec_from_file_location("starter_taxonomy", ROOT / "scripts/build_starter_taxonomy_supplements.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        taxon = ET.fromstring("<Taxon><TaxId>1</TaxId><ScientificName>Example species</ScientificName><Rank>species</Rank>"
            "<LineageEx><Taxon><Rank>genus</Rank><ScientificName>First</ScientificName></Taxon>"
            "<Taxon><Rank>genus</Rank><ScientificName>Conflicting</ScientificName></Taxon></LineageEx></Taxon>")
        values = module.rank_values(taxon)
        self.assertEqual(values["species"], "Example species")
        self.assertEqual(values["genus"], "Unresolved")
        self.assertEqual(values["family"], "Unresolved")


if __name__ == "__main__":
    unittest.main()
