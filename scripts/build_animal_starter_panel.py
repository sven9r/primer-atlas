#!/usr/bin/env python3
"""Freeze source-backed arthropod, bird and zooplankton templates from NCBI XML.

Offline extraction only; no fish templates or inter-group aliases are used.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


SPECS = {
    ("M21017.1", "18S"): ("ARTHROPODS", "Drosophila melanogaster", 1, 1995),
    ("NR_133562.1", "28S"): ("ARTHROPODS", "Drosophila melanogaster", 1, 3970),
    ("NC_001323.1", "12S_MT"): ("BIRDS", "Gallus gallus", 1297, 2272),
    ("NC_001323.1", "COI"): ("BIRDS", "Gallus gallus", 6645, 8192),
    ("AF173612.1", "18S"): ("BIRDS", "Gallus gallus", 1, 1737),
    ("XR_005840241.1", "28S"): ("BIRDS", "Gallus gallus", 1, 4444),
    ("NC_026914.1", "COI"): ("ZOOPLANKTON", "Daphnia magna", 1313, 2849),
    ("NC_026914.1", "16S_MT"): ("ZOOPLANKTON", "Daphnia magna", 12366, 13738),
    ("AF014011.1", "18S"): ("ZOOPLANKTON", "Daphnia pulex", 1, 2293),
    ("EU370436.1", "28S"): ("ZOOPLANKTON", "Daphnia magna", 1, 3823),
    ("NC_001322.1", "16S_MT"): ("ARTHROPODS", "Drosophila yakuba", 12756, 14081),
}
LOCATION_OVERRIDES = {("NC_026914.1", "COI"): "<1313..2849",
                      ("NC_026914.1", "16S_MT"): "complement(12366..13738)",
                      ("EU370436.1", "28S"): "<1..>3823",
                      ("NC_001322.1", "16S_MT"): "complement(12756..14081)"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
FIELDS = [
    "template_id", "sequence", "marker_id", "organism", "source_url", "version",
    "group_id", "target_group", "accession_version", "source_feature", "source_feature_location",
    "source_start", "source_end", "source_strand", "sequence_length",
    "sequence_sha256", "source_sequence_sha256", "raw_source_file", "raw_source_sha256",
    "source_status", "partial", "scope_note", "retrieved_date", "source_taxonomy",
    "ecology_source_key", "ecology_source_url", "ecology_scope",
]


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def freeze(inputs, output, retrieved_date):
    output.mkdir(parents=True, exist_ok=True)
    (output / "source").mkdir(exist_ok=True)
    rows, sources, seen, archived_records = [], [], set(), {}
    for index, path in enumerate(inputs, 1):
        original = path.read_bytes()
        records = ET.fromstring(original).findall("./GBSeq")
        if not records:
            raise ValueError(f"No NCBI GBSeq records: {path}")
        filename = f"source/ncbi-records-{index:02d}.xml.gz"
        compressed = gzip.compress(original, mtime=0)
        (output / filename).write_bytes(compressed)
        accessions = [r.findtext("GBSeq_accession-version") for r in records]
        sources.append({
            "file": filename, "sha256": sha256(compressed),
            "uncompressed_sha256": sha256(original), "accessions": accessions,
            "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            "params": {"db": "nuccore", "id": ",".join(accessions),
                       "rettype": "gbwithparts" if accessions == ["NW_007931121.1"] else "gb", "retmode": "xml"},
        })
        for record in records:
            accession = record.findtext("GBSeq_accession-version")
            archived_records[accession] = record
            organism = record.findtext("GBSeq_organism", "")
            taxonomy = record.findtext("GBSeq_taxonomy", "")
            definition = record.findtext("GBSeq_definition", "")
            parent_sequence = record.findtext("GBSeq_sequence", "").upper()
            if len(parent_sequence) != int(record.findtext("GBSeq_length", "0")):
                raise ValueError(f"Source length mismatch: {accession}")
            if set(parent_sequence) - set("ACGTRYSWKMBDHVN"):
                raise ValueError(f"Invalid nucleotide alphabet: {accession}")
            for feature in record.findall("./GBSeq_feature-table/GBFeature"):
                kind = feature.findtext("GBFeature_key", "")
                qualifiers = {q.findtext("GBQualifier_name"): q.findtext("GBQualifier_value", "")
                              for q in feature.findall("./GBFeature_quals/GBQualifier")}
                product = qualifiers.get("product", "")
                marker = None
                if kind == "CDS" and qualifiers.get("gene") == "COX1":
                    marker = "COI"
                elif kind == "rRNA":
                    if product == "l-rRNA":
                        marker = "16S_MT"
                    for label, value in (("12S", "12S_MT"), ("16S", "16S_MT"), ("18S", "18S"), ("28S", "28S")):
                        if label in product:
                            marker = value
                            break
                spec = SPECS.get((accession, marker))
                if spec is None:
                    continue
                group, expected_organism, start, end = spec
                location = feature.findtext("GBFeature_location", "")
                # M21017 also contains a terminal partial 18S repeat; do not use it.
                if location != LOCATION_OVERRIDES.get((accession, marker), f"{start}..{end}"):
                    continue
                intervals = feature.findall("./GBFeature_intervals/GBInterval")
                if len(intervals) != 1 or organism != expected_organism:
                    raise ValueError(f"Unexpected species/feature geometry: {accession} {marker}")
                raw_bounds = sorted((int(intervals[0].findtext("GBInterval_from", "0")), int(intervals[0].findtext("GBInterval_to", "0"))))
                if raw_bounds != [start, end]:
                    raise ValueError(f"Source interval changed: {accession} {marker}")
                required_lineage = "Arthropoda" if group in ("ARTHROPODS", "ZOOPLANKTON") else "Aves"
                if required_lineage not in taxonomy.split("; "):
                    raise ValueError(f"Target group not supported by source taxonomy: {accession}")
                if group == "ZOOPLANKTON" and "Daphnia" not in taxonomy.split("; "):
                    raise ValueError("Zooplankton membership is restricted to the source-backed Daphnia exemplars")
                identity = (accession, marker)
                if identity in seen:
                    raise ValueError(f"Duplicate source template: {identity}")
                seen.add(identity)
                sequence = parent_sequence[start - 1:end]
                strand = "-" if location.startswith("complement(") else "+"
                if strand == "-":
                    sequence = sequence.translate(COMPLEMENT)[::-1]
                predicted = definition.startswith("PREDICTED:")
                partial = "<" in location or ">" in location
                ecology_key = ("colbourne2011_daphnia" if organism == "Daphnia pulex" else "orsini2016_daphnia") if group == "ZOOPLANKTON" else ""
                ecology_url = ("https://doi.org/10.1126/science.1197761" if organism == "Daphnia pulex" else "https://doi.org/10.1038/sdata.2016.30") if group == "ZOOPLANKTON" else ""
                rows.append({
                    "template_id": f"{group}_{accession}_{marker}", "sequence": sequence,
                    "marker_id": marker, "organism": organism,
                    "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                    "version": "animal-starter-v1", "group_id": group, "target_group": group,
                    "accession_version": accession, "source_feature": f"{kind}: {product}",
                    "source_feature_location": location, "source_start": start, "source_end": end,
                    "source_strand": strand, "sequence_length": len(sequence),
                    "sequence_sha256": sha256(sequence.encode("ascii")),
                    "source_sequence_sha256": sha256(parent_sequence.encode("ascii")),
                    "raw_source_file": filename, "raw_source_sha256": sha256(compressed),
                    "source_status": "PREDICTED RefSeq MODEL" if predicted else "NCBI source annotation",
                    "partial": str(partial).lower(),
                    "scope_note": "Selected source-backed mapping exemplar for this animal group; not representative taxonomic coverage or PCR validation." +
                        (" Predicted genomic rRNA model; not an experimentally validated transcript." if predicted else "") +
                        (" Source feature is explicitly partial; a missing primer site does not establish a primer mismatch or locus absence." if partial else "") +
                        (" Freshwater Daphnia zooplankton exemplar; ecology is source-supported for the named species, not inferred from Arthropoda membership." if group == "ZOOPLANKTON" else ""),
                    "retrieved_date": retrieved_date, "source_taxonomy": taxonomy,
                    "ecology_source_key": ecology_key, "ecology_source_url": ecology_url,
                    "ecology_scope": "Freshwater zooplankton; selected Daphnia exemplar" if group == "ZOOPLANKTON" else "",
                })
    if seen != set(SPECS):
        raise ValueError(f"Missing expected exact source features: {set(SPECS) - seen}")
    # A reviewed RNA reference is not by itself proof of a contiguous DNA span.
    # Verify this rDNA copy against its named genomic source before claiming it.
    genomic_record = archived_records.get("NW_007931121.1")
    if genomic_record is None:
        raise ValueError("The Drosophila genomic correspondence source is required")
    genomic_sequence = genomic_record.findtext("GBSeq_sequence", "").upper()
    drosophila_28s = next(r for r in rows if r["accession_version"] == "NR_133562.1")
    if genomic_sequence[70954:74924] != drosophila_28s["sequence"]:
        raise ValueError("Reviewed Drosophila RNA does not match its contiguous genomic source")
    correspondence_features = [f for f in genomic_record.findall("./GBSeq_feature-table/GBFeature")
        if f.findtext("GBFeature_key") == "rRNA" and f.findtext("GBFeature_location") == "70955..74924"
        and any(q.findtext("GBQualifier_name") == "transcript_id" and q.findtext("GBQualifier_value") == "NR_133562.1"
                for q in f.findall("./GBFeature_quals/GBQualifier"))]
    if len(correspondence_features) != 1:
        raise ValueError("Named genomic feature does not identify NR_133562.1 uniquely")
    drosophila_28s["scope_note"] += " Its 3,970 encoded bases exactly match the contiguous annotated rDNA interval NW_007931121.1:70955-74924; this verifies this copy's reference span, not PCR success or other repeat copies."
    rows.sort(key=lambda row: (row["group_id"], row["marker_id"], row["template_id"]))
    with (output / "templates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "templates.fasta").write_text("".join(
        f">{row['template_id']} {row['organism']}\n{row['sequence']}\n" for row in rows
    ), encoding="ascii")
    manifest = {
        "schema_version": 1, "version": "animal-starter-v1", "created_date": retrieved_date,
        "status": "starter_mapping_panel", "n_templates": len(rows),
        "n_source_accessions": len(archived_records),
        "group_marker_counts": {group: {marker: sum(r["group_id"] == group and r["marker_id"] == marker for r in rows)
                                        for marker in ("12S_MT", "16S_MT", "COI", "18S", "28S")}
                                for group in ("ARTHROPODS", "BIRDS", "ZOOPLANKTON")},
        "scope": "Selected arthropod, bird and source-supported freshwater zooplankton exemplars. Not a representative animal reference release, not monthly maintained, and not PCR performance evidence. Fish templates remain in their separate panel.",
        "coordinate_contract": "Exact contiguous source-feature slices oriented by the annotated strand, with local coordinates starting at one. source_start/source_end preserve original raw genomic bounds; source_strand records reverse complements. Partial source boundaries remain explicit. No inter-species alignment or SSU/LSU domain boundaries are inferred.",
        "ecology_sources": [{"source_key": "colbourne2011_daphnia", "species": "Daphnia pulex", "doi": "10.1126/science.1197761", "url": "https://doi.org/10.1126/science.1197761", "scope": "Primary genome paper identifies D. pulex as a freshwater lake/pond grazer and forage species."},
                            {"source_key": "orsini2016_daphnia", "species": "Daphnia magna", "doi": "10.1038/sdata.2016.30", "url": "https://doi.org/10.1038/sdata.2016.30", "scope": "Primary transcriptome study identifies D. magna as a pelagic freshwater lake/pond grazer and prey species."}],
        "withheld_features": [{"accession_version": "M21017.1", "location": "join(3288..5075,5121..7232)", "reason": "Discontinuous source 28S feature; not concatenated into genomic PCR geometry."},
                              {"accession_version": "M21017.1", "location": "11727..>12026", "reason": "Partial terminal repeat; complete 18S feature is available separately."}],
        "genomic_correspondences": [{"template_accession_version": "NR_133562.1",
            "genomic_accession_version": "NW_007931121.1", "genomic_start": 70955, "genomic_end": 74924,
            "strand": "+", "sequence_sha256": drosophila_28s["sequence_sha256"],
            "proof": "Exact sequence equality to the contiguous genomic rRNA feature with transcript_id NR_133562.1; repeat-specific reference span only."}],
        "acquisition": "NCBI Entrez skill ncbi_entrez.py efetch, nuccore, gb XML. Source organism taxonomy, accession versions, sequence lengths and feature coordinates validated before freezing.",
        "builder": {"path": "scripts/build_animal_starter_panel.py", "sha256": sha256(Path(__file__).read_bytes())},
        "artifacts": {name: {"sha256": sha256((output / name).read_bytes()), "bytes": (output / name).stat().st_size}
                      for name in ("templates.csv", "templates.fasta")},
        "sources": sources,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("data/pinned/ANIMAL_STARTER/v1"))
    parser.add_argument("--retrieved-date", default="2026-10-31")
    args = parser.parse_args()
    rows = freeze(args.inputs, args.output, args.retrieved_date)
    print(f"Frozen {len(rows)} group-specific animal starter templates in {args.output}")
