#!/usr/bin/env python3
"""Acquire and select a fixed-size, reproducible UNITE ITS pilot panel.

Raw release files and per-sequence results stay under ignored data/external and
data/derived directories. This script never invokes the public release builder.
"""

import csv
import hashlib
import json
import pathlib
import tarfile
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "data/external/its_fungal_pilot"
DERIVED = ROOT / "data/derived/its_fungal_pilot"
DOI = "10.15156/BIO/3301229"
VERSION = "UNITE 10.0; 2025-02-19; Fungi general FASTA dynamic"
URL = "https://s3.hpc.ut.ee/plutof-public/original/9489f7bc-7cc1-4e0a-84dc-c732476b9acd.tgz"
ARCHIVE_SHA256 = "a340864f947c517b671f1054a261e613a6619ef88dd8000de1699f5e1f21413b"
FASTA_SHA256 = "1bb35a664b20ef03484fe74e857929ccd537d9d51c9eef0305b99a1486e7386e"
ARCHIVE = EXTERNAL / "sh_general_release_19.02.2025.tgz"
FASTA = EXTERNAL / "sh_general_release_dynamic_19.02.2025.fasta"
PHYLA = (
    "Ascomycota", "Basidiomycota", "Chytridiomycota", "Glomeromycota",
    "Mortierellomycota", "Mucoromycota", "Rozellomycota", "Zoopagomycota",
)
PER_PHYLUM = 6


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_source():
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists():
        with urllib.request.urlopen(URL, timeout=120) as response, ARCHIVE.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError(f"UNITE archive checksum changed: {digest}")
    if not FASTA.exists():
        with tarfile.open(ARCHIVE, "r:gz") as archive:
            member = archive.getmember(FASTA.name)
            with archive.extractfile(member) as source, FASTA.open("wb") as output:
                while chunk := source.read(1024 * 1024):
                    output.write(chunk)
    fasta_digest = hashlib.sha256(FASTA.read_bytes()).hexdigest()
    if fasta_digest != FASTA_SHA256:
        raise ValueError(f"UNITE FASTA checksum changed: {fasta_digest}")
    return digest


def parse_fasta():
    with FASTA.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            if not header.startswith(">") or not sequence or sequence.startswith(">"):
                raise ValueError("Expected the two-line UNITE general FASTA format")
            yield header[1:].strip(), sequence.strip().upper()


def taxonomy(header):
    parts = header.split("|")
    if len(parts) != 5:
        return None
    ranks = dict(item.split("__", 1) for item in parts[4].split(";") if "__" in item)
    if ranks.get("k") != "Fungi":
        return None
    return dict(accession=parts[1], species_hypothesis=parts[2],
                representative_type=parts[3], phylum=ranks.get("p", ""),
                class_name=ranks.get("c", ""), order_name=ranks.get("o", ""),
                family=ranks.get("f", ""), genus=ranks.get("g", ""),
                species=ranks.get("s", ""))


def select(records, n):
    # Equal long/short strata expose missing terminal sites as well as retained
    # sites. Prefer different classes, then genera, then INSDC accessions.
    chosen = []
    n_long = min(4, sum(650 <= r["length_bp"] <= 1200 for r in records))
    for label, lo, hi, count in (("long", 650, 1200, n_long), ("short", 300, 549, n - n_long)):
        candidates = [r for r in records if lo <= r["length_bp"] <= hi]
        candidates.sort(key=lambda r: (
            r["accession"].startswith("UDB"),
            r["representative_type"] != "refs",
            hashlib.sha256(r["accession"].encode()).hexdigest(),
        ))
        for diversity_key in ("class_name", "genus", None):
            for row in candidates:
                if row["accession"] in {x["accession"] for x in chosen}:
                    continue
                if diversity_key and row[diversity_key] in {x[diversity_key] for x in chosen}:
                    continue
                picked = dict(row)
                picked["length_stratum"] = label
                chosen.append(picked)
                if sum(x["length_stratum"] == label for x in chosen) == count:
                    break
            if sum(x["length_stratum"] == label for x in chosen) == count:
                break
        if sum(x["length_stratum"] == label for x in chosen) != count:
            raise ValueError(f"Insufficient {label} reference sequences")
    if len(chosen) != n:
        raise ValueError("Panel selection count changed")
    return chosen


def main():
    digest = get_source()
    DERIVED.mkdir(parents=True, exist_ok=True)
    by_phylum = {p: [] for p in PHYLA}
    total = 0
    for header, sequence in parse_fasta():
        total += 1
        taxon = taxonomy(header)
        if taxon and taxon["phylum"] in by_phylum and set(sequence) <= set("ACGTRYSWKMBDHVN"):
            by_phylum[taxon["phylum"]].append(dict(
                taxon, source_header=header, sequence=sequence, length_bp=len(sequence)
            ))
    panel = []
    for phylum in PHYLA:
        panel.extend(select(by_phylum[phylum], PER_PHYLUM))
    panel.sort(key=lambda r: (PHYLA.index(r["phylum"]), r["length_stratum"], r["accession"]))
    for index, row in enumerate(panel, 1):
        row["panel_id"] = f"ITS{index:03d}"
        row["marker_region"] = "ITS; exact ITS1/5.8S/ITS2 boundaries not annotated in source FASTA"
        row["source_doi"] = DOI
        row["source_version"] = VERSION
        row["sequence_sha256"] = hashlib.sha256(row["sequence"].encode()).hexdigest()
    fields = ["panel_id", "accession", "species_hypothesis", "representative_type",
              "phylum", "class_name", "order_name", "family", "genus", "species",
              "marker_region", "length_stratum", "length_bp", "source_doi",
              "source_version", "source_header", "sequence_sha256"]
    write_csv(DERIVED / "panel_manifest.csv", panel, fields)
    with (DERIVED / "panel.fasta").open("w") as handle:
        for row in panel:
            handle.write(f'>{row["panel_id"]}\n{row["sequence"]}\n')
    (DERIVED / "source.json").write_text(json.dumps(dict(
        doi=DOI, version=VERSION, url=URL, archive_sha256=digest,
        archive_member=FASTA.name, fasta_sha256=FASTA_SHA256, source_sequences=total,
        panel_sequences=len(panel), phyla=list(PHYLA),
        selection="6 per phylum: up to 4 sequences 650-1200 bp, remainder 300-549 bp; diverse classes and genera; INSDC and RefS preferred",
    ), indent=2) + "\n")
    print(f"Selected {len(panel)} of {total} UNITE sequences across {len(PHYLA)} phyla")


if __name__ == "__main__":
    main()
