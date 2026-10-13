#!/usr/bin/env python3
"""Indexed exact-site evidence for a versioned UNITE panel.

The search mirrors fungal_explorer.R's exact IUPAC/site-specificity rules.
It stores sparse site evidence instead of a huge forward/reverse Cartesian
table. Missing hits are search results, never PCR failures or gene absence.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import tempfile

IUPAC = dict(A="A", C="C", G="G", T="T", U="T", R="AG", Y="CT",
             S="GC", W="AT", K="GT", M="AC", B="CGT", D="AGT",
             H="ACT", V="ACG", N="ACGT", I="ACGT")
COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNI", "TGCAAYRSWMKVHDBNI")

def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def read_fasta(path):
    with open_text(path) as handle:
        header, bases = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(bases).upper()
                header, bases = line[1:].strip().split()[0], []
            else:
                bases.append(line.strip())
        if header is not None:
            yield header, "".join(bases).upper()

def read_csv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle))

def write_gzip_csv(path, fields, rows):
    # Deterministic compressed artifacts make replays byte-identical.
    import io
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)

def prepare_primer(sequence, direction):
    sequence = sequence.upper()
    if not sequence or any(base not in IUPAC for base in sequence):
        return None
    query = reverse_complement(sequence) if direction == "reverse" else sequence
    pattern = "".join(base if len(IUPAC[base]) == 1 else "[" + IUPAC[base] + "]" for base in query)
    return dict(regex=re.compile("(?=(" + pattern + "))"), width=len(query),
                exact_probability=math.prod(len(IUPAC[base]) / 4 for base in query))

def exact_site(sequence, prepared):
    empty = dict(start=None, end=None, n_matches=0, status="no_exact_binding_site")
    if prepared is None:
        return dict(empty, status="invalid_oligo")
    width = prepared["width"]
    if len(sequence) < width:
        return dict(empty, status="template_shorter_than_oligo")
    hits = [match.start() + 1 for match in prepared["regex"].finditer(sequence)]
    if not hits:
        return empty
    if len(hits) > 1:
        return dict(empty, n_matches=len(hits), status="ambiguous_binding_site")
    expectation = max(1, len(sequence) - width + 1) * prepared["exact_probability"]
    if width < 7 or expectation > .05:
        return dict(empty, n_matches=1, status="low_specificity_site")
    return dict(start=hits[0], end=hits[0]+width-1, n_matches=1, status="exact_compatible_site")

def pair_result(forward, reverse, minimum=100):
    f = forward["status"] == "exact_compatible_site"
    r = reverse["status"] == "exact_compatible_site"
    status = "both_sites_missing"
    if f and not r: status = "reverse_site_missing"
    if r and not f: status = "forward_site_missing"
    if "low_specificity_site" in (forward["status"], reverse["status"]): status = "low_specificity_site"
    if "ambiguous_binding_site" in (forward["status"], reverse["status"]): status = "ambiguous_binding_site"
    length = None
    if f and r:
        status = "outward_or_overlapping_sites"
        if forward["end"] < reverse["start"]:
            length = reverse["end"] - forward["start"] + 1
            status = "scorable" if length >= minimum else "below_minimum_length"
    return status, length if status == "scorable" else None

def build(panel, catalog_path, pair_path, output, fasta_path=None):
    output.mkdir(parents=True, exist_ok=True)
    catalog = read_csv(catalog_path)
    pairs = read_csv(pair_path)
    fasta_path = fasta_path or panel / "selected_raw.fasta.gz"
    records = list(read_fasta(fasta_path))
    if len({name for name, _ in records}) != len(records):
        raise ValueError("Repeated UNITE record identity")
    if len({row["catalog_id"] for row in catalog}) != len(catalog):
        raise ValueError("Repeated primer identity")
    prepared = [prepare_primer(row["sequence"], row["direction"]) for row in catalog]
    lookup = {(row["primer_name"],row["sequence"],row["direction"]): i for i,row in enumerate(catalog)}
    indices = [(lookup[(p["forward_primer"],p["forward_sequence"],"forward")],
                lookup[(p["reverse_primer"],p["reverse_sequence"],"reverse")]) for p in pairs]
    lengths = [[] for _ in pairs]
    counts = [{} for _ in pairs]
    sparse = []
    n_missing = 0
    for n, (record_id, sequence) in enumerate(records, 1):
        sites = [exact_site(sequence, primer) for primer in prepared]
        for row, site in zip(catalog, sites):
            if site["status"] == "no_exact_binding_site":
                n_missing += 1
                continue
            sparse.append(dict(record_id=record_id,catalog_id=row["catalog_id"],
                               primer_name=row["primer_name"],direction=row["direction"],**site))
        for j,(fi,ri) in enumerate(indices):
            state,bp = pair_result(sites[fi],sites[ri])
            counts[j][state] = counts[j].get(state,0)+1
            if bp is not None: lengths[j].append(bp)
        if n % 1000 == 0: print(f"Searched {n}/{len(records)} UNITE records",flush=True)
    fields = ["record_id","catalog_id","primer_name","direction","start","end","n_matches","status"]
    write_gzip_csv(output / "primer_sites.csv.gz",fields,sparse)
    summaries = []
    for pair, state, values in zip(pairs, counts, lengths):
        row = dict(pair, n_panel=len(records), n_scorable=len(values),
                   amplicon_min_bp=min(values) if values else None,
                   amplicon_median_bp=statistics.median(values) if values else None,
                   amplicon_max_bp=max(values) if values else None)
        row.update({"n_"+key:state.get(key,0) for key in ("scorable","both_sites_missing","reverse_site_missing",
            "forward_site_missing","ambiguous_binding_site","low_specificity_site","outward_or_overlapping_sites","below_minimum_length")})
        summaries.append(row)
    extra = sorted(set(summaries[0]) - set(pairs[0])) if summaries else []
    write_gzip_csv(output / "pairs.csv.gz",list(pairs[0])+extra if pairs else [],summaries)
    write_gzip_csv(output / "catalog.csv.gz",list(catalog[0]),catalog)
    meta = dict(n_records=len(records),n_catalog=len(catalog),n_pairs=len(pairs),
                n_site_searches=len(records)*len(catalog),n_sparse_rows=len(sparse),n_no_exact_hits=n_missing,
                sequence_file=str(fasta_path),sequence_sha256=hashlib.sha256(fasta_path.read_bytes()).hexdigest(),
                catalog_sha256=hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
                pair_catalog_sha256=hashlib.sha256(pair_path.read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                search_rules="Unique exact IUPAC-compatible A/C/G/T template sites; overlapping hits counted; template ambiguity never establishes a site; exact chance expectation <=0.05; inward products >=100 encoded bp including primers",
                limitation="Missing exact hits can reflect truncated records or primer mismatches; no PCR probability or whole-taxon coverage inferred.")
    (output / "site_search_manifest.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps({k:v for k,v in meta.items() if k.startswith("n_")},indent=2))
    return meta

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--panel",type=Path,default=Path("data/derived/unite_reference_panel"))
    parser.add_argument("--catalog",type=Path,default=Path("data/pinned/ITS_FUNGAL/versions/v3-c6fa87f256e1/catalog.csv"))
    parser.add_argument("--pairs",type=Path,default=Path("data/pinned/ITS_FUNGAL/versions/v3-c6fa87f256e1/pairs.csv"))
    parser.add_argument("--output",type=Path,default=Path("data/derived/unite_reference_panel/evidence"))
    parser.add_argument("--fasta",type=Path)
    args=parser.parse_args()
    build(args.panel,args.catalog,args.pairs,args.output,args.fasta)
