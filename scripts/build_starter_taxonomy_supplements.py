#!/usr/bin/env python3
"""Attach archived NCBI taxon hierarchies without changing template sequences."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
FIELDS = [
    "template_id", "accession_version", "marker_id", "group_id", "target_group", "taxid",
    *RANKS, "taxonomy_status", "unresolved_ranks", "taxonomy_source_url", "taxonomy_version",
    "taxonomy_source_sha256", "raw_taxonomy_file", "raw_taxonomy_sha256", "raw_taxonomy_xml_sha256",
    "source_record_file", "source_record_sha256", "source_record_taxid", "source_record_organism",
    "retrieved_date",
]


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def records_from_xml(raw):
    return ET.fromstring(raw).findall("./GBSeq")


def record_taxid(record):
    values = set()
    for feature in record.findall("./GBSeq_feature-table/GBFeature"):
        if feature.findtext("GBFeature_key") != "source":
            continue
        for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
            value = qualifier.findtext("GBQualifier_value", "")
            if qualifier.findtext("GBQualifier_name") == "db_xref" and value.startswith("taxon:"):
                values.add(value.split(":", 1)[1])
    if len(values) != 1:
        raise ValueError("A source record must identify exactly one NCBI taxon")
    return values.pop()


def rank_values(taxon):
    possibilities = {rank: set() for rank in RANKS}
    for node in [*taxon.findall("./LineageEx/Taxon"), taxon]:
        rank = node.findtext("Rank", "")
        if rank in possibilities:
            possibilities[rank].add(node.findtext("ScientificName", ""))
    # Missing or conflicting ranks stay unresolved; never infer a rank from a name.
    return {rank: next(iter(names)) if len(names) == 1 and all(names) else "Unresolved"
            for rank, names in possibilities.items()}


def build(root, taxonomy_xml, dyakuba_xml, date):
    taxonomy_raw = taxonomy_xml.read_bytes()
    taxa = {t.findtext("TaxId"): t for t in ET.fromstring(taxonomy_raw).findall("./Taxon")}
    if set(taxa) != {"7955", "8049", "7830", "8030", "7227", "9031", "7245", "6669", "35525"}:
        raise ValueError("Expected the nine explicitly requested source taxa")
    taxonomy_gz = gzip.compress(taxonomy_raw, mtime=0)
    dyakuba_raw = dyakuba_xml.read_bytes()
    dyakuba_records = records_from_xml(dyakuba_raw)
    if len(dyakuba_records) != 1 or dyakuba_records[0].findtext("GBSeq_accession-version") != "NC_001322.1":
        raise ValueError("Dyakuba source must be accession-version NC_001322.1")
    builder_hash = sha256(Path(__file__).read_bytes())
    for panel in ("FISH_STARTER", "ANIMAL_STARTER"):
        directory = root / "data/pinned" / panel / "v1"
        manifest_path = directory / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        # Immutable scientific template artifacts are read, never rebuilt here.
        before = {name: sha256((directory / name).read_bytes()) for name in ("templates.csv", "templates.fasta")}
        with (directory / "templates.csv").open(encoding="utf-8", newline="") as handle:
            templates = list(csv.DictReader(handle))
        records = {}
        for source in manifest["sources"]:
            compressed = (directory / source["file"]).read_bytes()
            if sha256(compressed) != source["sha256"]:
                raise ValueError(f"Corrupt original source archive: {source['file']}")
            for record in records_from_xml(gzip.decompress(compressed)):
                records[record.findtext("GBSeq_accession-version")] = (record, source["file"], sha256(compressed))
        tax_dir = directory / "taxonomy_source"
        tax_dir.mkdir(exist_ok=True)
        tax_filename = f"taxonomy_source/ncbi-taxonomy-{date.replace('-', '')}.xml.gz"
        (directory / tax_filename).write_bytes(taxonomy_gz)
        supplemental_sources = [{
            "file": tax_filename, "sha256": sha256(taxonomy_gz), "uncompressed_sha256": sha256(taxonomy_raw),
            "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            "params": {"db": "taxonomy", "id": "7955,8049,7830,8030,7227,9031,7245,6669,35525", "retmode": "xml"},
        }]
        external_template_ids = []
        if panel == "ANIMAL_STARTER":
            dyakuba_gz = gzip.compress(dyakuba_raw, mtime=0)
            dyakuba_filename = "taxonomy_source/ncbi-NC_001322.1.xml.gz"
            (directory / dyakuba_filename).write_bytes(dyakuba_gz)
            records["NC_001322.1"] = (dyakuba_records[0], dyakuba_filename, sha256(dyakuba_gz))
            templates.append({"template_id": "NC_001322.1_COX1", "accession_version": "NC_001322.1",
                              "marker_id": "COI", "group_id": "ARTHROPODS", "organism": "Drosophila yakuba"})
            external_template_ids.append("NC_001322.1_COX1")
            supplemental_sources.append({"file": dyakuba_filename, "sha256": sha256(dyakuba_gz),
                "uncompressed_sha256": sha256(dyakuba_raw),
                "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                "params": {"db": "nuccore", "id": "NC_001322.1", "rettype": "gb", "retmode": "xml"}})
        rows = []
        for template in templates:
            record, source_file, source_hash = records[template["accession_version"]]
            source_taxid = record_taxid(record)
            taxon = taxa.get(source_taxid)
            if taxon is None:
                raise ValueError(f"No taxonomy response for exact source taxid {source_taxid}")
            scientific_name = taxon.findtext("ScientificName", "")
            if record.findtext("GBSeq_organism") != scientific_name or template["organism"] != scientific_name:
                raise ValueError(f"Template, accession source and taxonomy names disagree: {template['template_id']}")
            ranks = rank_values(taxon)
            group = template["group_id"]
            if group == "ARTHROPODS" and ranks["phylum"] != "Arthropoda":
                raise ValueError("Arthropod template lacks authoritative Arthropoda membership")
            if group == "BIRDS" and ranks["class"] != "Aves":
                raise ValueError("Bird template lacks authoritative Aves membership")
            if group == "FISH" and ranks["class"] not in ("Actinopteri", "Chondrichthyes"):
                raise ValueError("Unexpected starter fish class; review source rather than guessing")
            if group == "ZOOPLANKTON" and (ranks["genus"] != "Daphnia" or not template.get("ecology_source_url")):
                raise ValueError("Zooplankton exemplar must retain Daphnia taxonomy and reviewed ecology source")
            missing = [rank for rank in RANKS if ranks[rank] == "Unresolved"]
            rows.append({"template_id": template["template_id"], "accession_version": template["accession_version"],
                "marker_id": template["marker_id"], "group_id": group, "target_group": group,
                "taxid": source_taxid, **ranks,
                "taxonomy_status": "source_taxid_verified" if not missing else "source_taxid_verified_unresolved_ranks",
                "unresolved_ranks": "|".join(missing),
                "taxonomy_source_url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={source_taxid}",
                "taxonomy_version": f"NCBI-Taxonomy-{date}-v1", "taxonomy_source_sha256": sha256(taxonomy_raw),
                "raw_taxonomy_file": tax_filename, "raw_taxonomy_sha256": sha256(taxonomy_gz),
                "raw_taxonomy_xml_sha256": sha256(taxonomy_raw),
                "source_record_file": source_file, "source_record_sha256": source_hash,
                "source_record_taxid": source_taxid, "source_record_organism": record.findtext("GBSeq_organism"),
                "retrieved_date": date})
        if len(rows) != len({r["template_id"] for r in rows}):
            raise ValueError("Taxonomy identities must remain one row per template")
        with (directory / "taxonomy.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, FIELDS, quoting=csv.QUOTE_ALL, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        artifact = directory / "taxonomy.csv"
        manifest["artifacts"]["taxonomy.csv"] = {"sha256": sha256(artifact.read_bytes()), "bytes": artifact.stat().st_size}
        manifest["taxonomy_supplement"] = {
            "version": f"NCBI-Taxonomy-{date}-v1", "retrieved_date": date, "n_rows": len(rows),
            "n_panel_template_rows": len(rows) - len(external_template_ids), "external_template_ids": external_template_ids,
            "scope": "Authoritative NCBI rank values joined through exact accession-source taxon IDs. No names inferred from FASTA headers, no unresolved rank filled from another rank, and no new sequence membership or scoring denominator.",
            "missing_rank_label": "Unresolved", "rank_fields": list(RANKS),
            "taxonomy_source_sha256_definition": "SHA-256 of the uncompressed original TaxaSet XML response; raw_taxonomy_sha256 hashes its stored gzip archive.",
            "builder": {"path": "scripts/build_starter_taxonomy_supplements.py", "sha256": builder_hash},
            "sources": supplemental_sources,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        after = {name: sha256((directory / name).read_bytes()) for name in before}
        if before != after:
            raise ValueError("Taxonomy supplementation changed immutable template artifacts")
        print(f"{panel}: {len(rows)} source-verified taxonomy rows; template bytes unchanged")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("taxonomy_xml", type=Path)
    parser.add_argument("dyakuba_xml", type=Path)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--date", default="2026-10-31")
    args = parser.parse_args()
    build(args.root, args.taxonomy_xml, args.dyakuba_xml, args.date)
