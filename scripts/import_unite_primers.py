#!/usr/bin/env python3
"""Create a versioned, runtime-independent snapshot of the UNITE primer table."""

from __future__ import annotations

import argparse
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote_plus

import pandas as pd


DEFAULT_URL = "https://unite.ut.ee/primers.php"
IUPAC = set("ACGTRYSWKMBDHVNI")


def reference_key(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return f"its_primary_{slug or 'not_reported'}"


def reference_url(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if not text or "unpublished" in text.lower():
        return DEFAULT_URL
    return "https://search.crossref.org/?q=" + quote_plus(text)


def clean_sequence(value: object) -> str:
    sequence = "" if pd.isna(value) else str(value).upper().replace(" ", "")
    return sequence if 15 <= len(sequence) <= 80 and set(sequence) <= IUPAC else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default=DEFAULT_URL)
    parser.add_argument("--output", default="data/catalog/unite_primers.csv")
    args = parser.parse_args()

    tables = pd.read_html(args.source)
    if not tables:
        raise RuntimeError("UNITE primer page contained no readable table")
    table = tables[0].copy()
    expected = {
        "Primer name", "Position in gene", "Sequence", "Gene/locus",
        "Direction", "Target", "Remarks", "Reference",
    }
    if not expected.issubset(table.columns):
        raise RuntimeError(f"UNITE table schema changed: {list(table.columns)}")

    table["sequence"] = table["Sequence"].map(clean_sequence)
    table["direction"] = table["Direction"].astype(str).str.lower().map(
        {"fwd": "forward", "rev": "reverse"}
    )
    table = table[table["sequence"].ne("") & table["direction"].notna()].copy()
    table["primer_name"] = table["Primer name"].astype(str).str.strip()
    table["reported_position"] = pd.to_numeric(
        table["Position in gene"], errors="coerce"
    ).astype("Int64")
    table["gene_locus"] = table["Gene/locus"].fillna("").astype(str).str.strip()
    table["target"] = table["Target"].fillna("").astype(str).str.strip()
    table["remarks"] = table["Remarks"].fillna("").astype(str).str.strip()
    table["reference"] = table["Reference"].fillna("").astype(str).str.strip()
    table["primary_reference_key"] = table["reference"].map(reference_key)
    table["primary_reference_url"] = table["reference"].map(reference_url)
    table["primary_reference_status"] = table["reference"].map(
        lambda value: "not_reported" if not value else (
            "unpublished" if "unpublished" in value.lower() else "publication_cited"
        )
    )
    table["marker_id"] = "ITS_FUNGAL"
    table["source_url"] = DEFAULT_URL
    table["source_snapshot_utc"] = datetime.now(timezone.utc).replace(
        microsecond=0
    ).isoformat()
    table["source_license_note"] = (
        "Factual primer metadata; retain UNITE attribution and cited primary references"
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "marker_id", "primer_name", "sequence", "direction",
        "reported_position", "gene_locus", "target", "remarks", "reference",
        "primary_reference_key", "primary_reference_url", "primary_reference_status",
        "source_url", "source_snapshot_utc", "source_license_note",
    ]
    table[columns].drop_duplicates(
        subset=["primer_name", "sequence", "direction"]
    ).to_csv(output, index=False)
    print(f"Wrote {len(table)} valid UNITE primer records to {output}")


if __name__ == "__main__":
    main()
