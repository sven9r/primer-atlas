#!/usr/bin/env python3
"""Freeze a small, source-traceable fish mapping panel from archived NCBI XML.

This is a starter geometry panel, not a representative fish reference release.
The network acquisition is recorded in manifest.json; this extractor is offline.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


MITO_ACCESSIONS = {"NC_002333.2", "NC_002081.1", "NC_001950.1"}
RNA_MARKERS = {
    "FJ915075.1": "18S",
    "FJ710886.1": "18S",
    "XR_012407110.1": "18S",
    "XR_012403044.1": "28S",
}
FIELDS = [
    "template_id", "sequence", "marker_id", "organism", "source_url", "version",
    "group_id", "accession_version", "source_feature", "source_feature_location",
    "source_start", "source_end", "source_strand", "sequence_length",
    "sequence_sha256", "source_sequence_sha256", "raw_source_file", "raw_source_sha256",
    "source_status", "partial", "scope_note", "retrieved_date",
]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def freeze(inputs: list[Path], output: Path, retrieved_date: str) -> list[dict]:
    output.mkdir(parents=True, exist_ok=True)
    raw_dir = output / "source"
    raw_dir.mkdir(exist_ok=True)
    rows: list[dict] = []
    sources: list[dict] = []
    seen: set[str] = set()
    for index, path in enumerate(inputs, start=1):
        original = path.read_bytes()
        records = ET.fromstring(original).findall("./GBSeq")
        if not records:
            raise ValueError(f"No complete GBSeq records in {path}")
        frozen_name = f"ncbi-records-{index:02d}.xml.gz"
        compressed = gzip.compress(original, mtime=0)
        (raw_dir / frozen_name).write_bytes(compressed)
        accessions = [r.findtext("GBSeq_accession-version") for r in records]
        sources.append({
            "file": f"source/{frozen_name}", "sha256": sha256(compressed),
            "uncompressed_sha256": sha256(original), "accessions": accessions,
            "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            "params": {"db": "nuccore", "id": ",".join(accessions), "rettype": "gb", "retmode": "xml"},
        })
        for record in records:
            accession = record.findtext("GBSeq_accession-version")
            if accession in seen:
                raise ValueError(f"Repeated source accession {accession}")
            seen.add(accession)
            source_sequence = record.findtext("GBSeq_sequence", "").upper()
            if len(source_sequence) != int(record.findtext("GBSeq_length", "0")):
                raise ValueError(f"Sequence-length mismatch in {accession}")
            if not source_sequence or set(source_sequence) - set("ACGTRYSWKMBDHVN"):
                raise ValueError(f"Unsupported sequence alphabet in {accession}")
            organism = record.findtext("GBSeq_organism")
            definition = record.findtext("GBSeq_definition", "")
            for feature in record.findall("./GBSeq_feature-table/GBFeature"):
                kind = feature.findtext("GBFeature_key")
                qualifiers = {
                    q.findtext("GBQualifier_name"): q.findtext("GBQualifier_value", "")
                    for q in feature.findall("./GBFeature_quals/GBQualifier")
                }
                marker = None
                if accession in MITO_ACCESSIONS:
                    if kind == "CDS" and qualifiers.get("gene") == "COX1":
                        marker = "COI"
                    elif kind == "rRNA" and qualifiers.get("product") in ("s-rRNA", "12S ribosomal RNA"):
                        marker = "12S_MT"
                    elif kind == "rRNA" and qualifiers.get("product") in ("l-rRNA", "16S ribosomal RNA"):
                        marker = "16S_MT"
                elif accession in RNA_MARKERS and kind == "rRNA":
                    marker = RNA_MARKERS[accession]
                if marker is None:
                    continue
                intervals = feature.findall("./GBFeature_intervals/GBInterval")
                if len(intervals) != 1:
                    raise ValueError(f"Joined/unknown feature intervals withheld: {accession} {marker}")
                interval = intervals[0]
                start = int(interval.findtext("GBInterval_from", "0"))
                end = int(interval.findtext("GBInterval_to", "0"))
                location = feature.findtext("GBFeature_location", "")
                if not (1 <= start <= end <= len(source_sequence)) or "complement" in location:
                    raise ValueError(f"Unsupported feature orientation/bounds: {accession} {location}")
                sequence = source_sequence[start - 1:end]
                partial = "<" in location or ">" in location or "partial" in definition.lower()
                predicted = definition.upper().startswith("PREDICTED:")
                rows.append({
                    "template_id": f"FISH_{accession}_{marker}", "sequence": sequence,
                    "marker_id": marker, "organism": organism,
                    "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                    "version": "fish-starter-v1", "group_id": "FISH",
                    "accession_version": accession,
                    "source_feature": f"{kind}: {qualifiers.get('product', '')}",
                    "source_feature_location": location, "source_start": start, "source_end": end,
                    "source_strand": "+", "sequence_length": len(sequence),
                    "sequence_sha256": sha256(sequence.encode("ascii")),
                    "source_sequence_sha256": sha256(source_sequence.encode("ascii")),
                    "raw_source_file": f"source/{frozen_name}", "raw_source_sha256": sha256(compressed),
                    "source_status": "PREDICTED RefSeq MODEL" if predicted else "NCBI source annotation",
                    "partial": str(partial).lower(),
                    "scope_note": "Small selected fish mapping panel; not population coverage or PCR validation. " +
                        ("Predicted genomic rRNA model; not an experimentally validated transcript. " if predicted else "") +
                        ("Source sequence is partial; missing primer sites do not establish locus absence." if partial else ""),
                    "retrieved_date": retrieved_date,
                })
    if len(rows) != 13 or len({row["template_id"] for row in rows}) != 13:
        raise ValueError("Expected exactly thirteen unique frozen starter templates")
    with (output / "templates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    fasta = "".join(f">{row['template_id']} {row['organism']}\n{row['sequence']}\n" for row in rows)
    (output / "templates.fasta").write_text(fasta, encoding="ascii")
    artifacts = {
        name: {"sha256": sha256((output / name).read_bytes()), "bytes": (output / name).stat().st_size}
        for name in ("templates.csv", "templates.fasta")
    }
    script = Path(__file__).resolve()
    manifest = {
        "schema_version": 1, "version": "fish-starter-v1", "created_date": retrieved_date,
        "status": "starter_mapping_panel", "n_templates": len(rows),
        "n_source_accessions": len(seen),
        "marker_counts": {m: sum(r["marker_id"] == m for r in rows) for m in ("12S_MT", "16S_MT", "COI", "18S", "28S")},
        "scope": "Selected NCBI records for interactive primer geometry. Not a representative fish reference release, not monthly maintained, and not PCR performance evidence.",
        "coordinate_contract": "Each template is an exact plus-strand slice of one archived accession-version feature. Local template coordinates start at one. source_start/source_end map back to the original record; no inter-template alignment or domain-coordinate projection is claimed.",
        "acquisition": "NCBI Entrez skill ncbi_entrez.py efetch, nuccore, gb XML. Returned accession versions and lengths were checked before freezing.",
        "builder": {"path": "scripts/build_fish_starter_panel.py", "sha256": sha256(script.read_bytes())},
        "artifacts": artifacts, "sources": sources,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("data/pinned/FISH_STARTER/v1"))
    parser.add_argument("--retrieved-date", default="2026-10-31")
    args = parser.parse_args()
    result = freeze(args.inputs, args.output, args.retrieved_date)
    print(f"Frozen {len(result)} starter templates in {args.output}")
