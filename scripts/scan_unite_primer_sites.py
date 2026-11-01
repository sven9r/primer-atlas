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
RAW_ALPHABET = set("ACGTRYSWKMBDHVN")
SCAN_RULE_VERSION = "exact_iupac_unique_sites_v1"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def record_ids_sha256(ids):
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()

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
                words = line[1:].strip().split()
                if not words:
                    raise ValueError("Empty FASTA record identity")
                header, bases = words[0], []
            else:
                if line.strip() and header is None:
                    raise ValueError("Sequence data precedes a FASTA record")
                bases.append(line.strip())
        if header is not None:
            yield header, "".join(bases).upper()

def validate_raw_records(records):
    if not records:
        raise ValueError("The raw UNITE panel is empty")
    if len({name for name, _ in records}) != len(records):
        raise ValueError("Repeated UNITE record identity")
    for name, sequence in records:
        if not sequence or set(sequence) - RAW_ALPHABET:
            raise ValueError("UNITE raw template must be nonempty and ungapped IUPAC DNA: " + name)

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
    # Capture executed source and immutable inputs before the long scan.
    script_bytes = Path(__file__).read_bytes()
    script_sha256 = hashlib.sha256(script_bytes).hexdigest()
    output.mkdir(parents=True, exist_ok=True)
    catalog = read_csv(catalog_path)
    pairs = read_csv(pair_path)
    fasta_path = fasta_path or panel / "selected_raw.fasta.gz"
    sequence_sha256 = sha(fasta_path)
    catalog_sha256, pair_catalog_sha256 = sha(catalog_path), sha(pair_path)
    records = list(read_fasta(fasta_path))
    validate_raw_records(records)
    metadata_path = panel / "selected_metadata.csv.gz"
    metadata_sha256 = sha(metadata_path)
    metadata = read_csv(metadata_path)
    metadata_ids = [row["record_id"] for row in metadata]
    if len(set(metadata_ids)) != len(metadata_ids) or set(metadata_ids) != {name for name, _ in records}:
        raise ValueError("UNITE site scan must cover exactly the retained panel membership")
    if not catalog or not pairs:
        raise ValueError("The frozen primer/pair catalogs are empty")
    if len({row["catalog_id"] for row in catalog}) != len(catalog):
        raise ValueError("Repeated primer identity")
    if len({row["pair_id"] for row in pairs}) != len(pairs):
        raise ValueError("Repeated pair identity")
    if any(row["direction"] not in ("forward", "reverse") for row in catalog):
        raise ValueError("Primer direction must be forward or reverse")
    if len({(row["primer_name"],row["sequence"],row["direction"]) for row in catalog}) != len(catalog):
        raise ValueError("Repeated primer component identity")
    previous = output / "site_search_manifest.json"
    if previous.exists():
        saved = json.loads(previous.read_text())
        identity = dict(sequence_sha256=sequence_sha256, catalog_sha256=catalog_sha256,
                        pair_catalog_sha256=pair_catalog_sha256, script_sha256=script_sha256,
                        metadata_sha256=metadata_sha256)
        if not saved.get("complete") or any(saved.get(key) != value for key,value in identity.items()):
            raise ValueError("Different immutable site scan exists; choose a new output directory")
        for name,key in (("primer_sites.csv.gz","site_artifact_sha256"),("catalog.csv.gz","catalog_artifact_sha256"),("pairs.csv.gz","pair_artifact_sha256")):
            if sha(output/name) != saved[key]:
                raise ValueError("Site scan artifact checksum differs: " + name)
        return saved
    if any(output.iterdir()):
        raise ValueError("Incomplete site scan exists; choose a new output directory")
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
    # A reader may impute omitted rows only after verifying this complete proof.
    meta = dict(complete=True,scan_rule_version=SCAN_RULE_VERSION,
                n_records=len(records),n_catalog=len(catalog),n_pairs=len(pairs),
                n_site_searches=len(records)*len(catalog),n_sparse_rows=len(sparse),n_no_exact_hits=n_missing,
                sequence_file=str(fasta_path),sequence_sha256=sequence_sha256,
                record_ids_sha256=record_ids_sha256(name for name,_ in records),
                metadata_sha256=metadata_sha256,catalog_sha256=catalog_sha256,
                pair_catalog_sha256=pair_catalog_sha256,script_sha256=script_sha256,
                catalog_artifact_sha256=sha(output/"catalog.csv.gz"),
                pair_artifact_sha256=sha(output/"pairs.csv.gz"),
                site_artifact_sha256=sha(output/"primer_sites.csv.gz"),
                search_rules="Unique exact IUPAC-compatible A/C/G/T template sites; overlapping hits counted; template ambiguity never establishes a site; exact chance expectation <=0.05; inward products >=100 encoded bp including primers",
                limitation="Missing exact hits can reflect truncated records or primer mismatches; no PCR probability or whole-taxon coverage inferred.")
    if meta["n_sparse_rows"] + meta["n_no_exact_hits"] != meta["n_site_searches"]:
        raise ValueError("The sparse site proof does not cover every record/primer search")
    if sha(fasta_path)!=sequence_sha256 or sha(catalog_path)!=catalog_sha256 or sha(pair_path)!=pair_catalog_sha256 or sha(metadata_path)!=metadata_sha256:
        raise ValueError("An input changed while scanning; no completed proof can be emitted")
    (output/"code_scanner.py").write_bytes(script_bytes)
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
