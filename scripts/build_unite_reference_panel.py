#!/usr/bin/env python3
"""Build a provenance-preserving, order-balanced UNITE ITS reference panel.

Uses the pinned cached general release; never scrapes sequence search pages or
edits earlier fungal snapshots. Developer sequences are retained as published,
not padded into a complete ribosomal operon. No primer performance is inferred.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DOI = "10.15156/BIO/3301229"
SOURCE_VERSION = "UNITE 10.0; 2025-02-19; Fungi general FASTA dynamic"
SOURCE_URL = "https://s3.hpc.ut.ee/plutof-public/original/9489f7bc-7cc1-4e0a-84dc-c732476b9acd.tgz"
STANDARD_MEMBER = "sh_general_release_dynamic_19.02.2025.fasta"
DEVELOPER_MEMBER = "sh_general_release_dynamic_19.02.2025_dev.fasta"
PINNED_HASHES = {
    "archive": "a340864f947c517b671f1054a261e613a6619ef88dd8000de1699f5e1f21413b",
    "standard": "1bb35a664b20ef03484fe74e857929ccd537d9d51c9eef0305b99a1486e7386e",
    "developer": "e75ef3e66b3af93d0f545ee15cf4690c612802a679faab51ba915d46b6ff9c72",
}
RANKS = ("kingdom", "phylum", "class_name", "order_name", "family", "genus", "species")
PREFIXES = dict(zip("kpcofgs", RANKS))
COVERAGE_RANKS = ("phylum", "class", "order", "family")
RANK_COLUMN = {"phylum": "phylum", "class": "class_name", "order": "order_name", "family": "family"}
ALPHABET = set("ACGTRYSWKMBDHVN")
STRATA = ("lt300", "300_499", "500_649", "650_999", "1000_1999", "ge2000")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def sha_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def rank_status(value):
    normalized = re.sub(r"[_\s]+", " ", (value or "").strip()).lower()
    if "incertae sedis" in normalized:
        return "incertae_sedis"
    if not normalized or normalized in {"na", "n/a", "none", "not assigned", "not available"} or re.match(
            r"^(?:unidentified|unclassified|unknown|unresolved|uncultured)(?:\b|$)", normalized):
        return "unresolved"
    return "resolved"


def parse_header(header):
    pieces = header.split("|")
    if len(pieces) != 5 or not pieces[1] or not pieces[2]:
        raise ValueError("Expected organism|accession|SH|representative|taxonomy header")
    taxonomy = {}
    for token in pieces[4].split(";"):
        if "__" not in token:
            continue
        prefix, value = token.split("__", 1)
        if prefix in PREFIXES:
            column = PREFIXES[prefix]
            if column in taxonomy:
                raise ValueError("Repeated taxonomy rank in source header")
            taxonomy[column] = value
    result = dict(organism=pieces[0], accession=pieces[1], accession_raw=pieces[1],
                  species_hypothesis=pieces[2], representative_type=pieces[3])
    result["accession_version_status"] = "reported" if re.search(r"\.\d+$", pieces[1]) else "not_reported"
    for column in RANKS:
        result[column] = taxonomy.get(column, "")
        result[column + "_status"] = rank_status(result[column])
    result["record_id"] = "UNITE_" + hashlib.sha256((pieces[2] + "\0" + pieces[1]).encode()).hexdigest()[:20]
    return result


def group_identity(row, rank):
    columns = ("kingdom", "phylum", "class_name", "order_name", "family")
    end = columns.index(RANK_COLUMN[rank]) + 1
    identity = []
    for column in columns[:end]:
        value = row.get(column, "")
        identity.append("__UNRESOLVED__" if rank_status(value) == "unresolved" else value)
    payload = json.dumps([rank, identity], ensure_ascii=False, separators=(",", ":"))
    return "UNITE_" + rank.upper() + "_" + hashlib.sha256(payload.encode()).hexdigest()[:20]


def length_stratum(length):
    return STRATA[0 if length < 300 else 1 if length < 500 else 2 if length < 650 else
                  3 if length < 1000 else 4 if length < 2000 else 5]


def fasta_records(path):
    header, pieces, ordinal = None, [], 0
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield ordinal, header, "".join(pieces).upper()
                ordinal += 1
                header, pieces = line[1:], []
            elif line.strip():
                if header is None:
                    raise ValueError("Sequence data precedes the first FASTA header")
                pieces.append(line.strip())
        if header is not None:
            yield ordinal, header, "".join(pieces).upper()


def read_source(path, kind):
    records, issues = collections.defaultdict(list), []
    total = 0
    for ordinal, header, sequence in fasta_records(path):
        total += 1
        try:
            metadata = parse_header(header)
            key = (metadata["species_hypothesis"], metadata["accession"])
        except ValueError as error:
            metadata = {column: "" for column in RANKS}
            key = ("__MALFORMED__" + kind, str(ordinal))
            issues.append(dict(source=kind, ordinal=ordinal, issue="malformed_header", reason=str(error),
                               source_header=header, sequence_sha256=sequence_sha(sequence)))
        records[key].append(dict(metadata=metadata, header=header, sequence=sequence, ordinal=ordinal))
    return records, issues, total


def pair_sources(standard_path, developer_path):
    standard, issues, n_standard = read_source(standard_path, "standard")
    developer, developer_issues, n_developer = read_source(developer_path, "developer")
    issues.extend(developer_issues)
    sh_accessions = collections.defaultdict(set)
    for sh, accession in set(standard) | set(developer):
        if not sh.startswith("__MALFORMED__"):
            sh_accessions[sh].add(accession)
    rows, sequences = [], {}
    for key in sorted(set(standard) | set(developer)):
        std, dev = standard.get(key, []), developer.get(key, [])
        first = std[0] if std else dev[0]
        row = dict(first["metadata"])
        if "record_id" not in row:
            row.update(record_id="UNITE_INVALID_" + hashlib.sha256(repr(key).encode()).hexdigest()[:20],
                       accession="", accession_raw="", species_hypothesis="", representative_type="", organism="")
        sequence = dev[0]["sequence"] if dev else ""
        standard_sequence = std[0]["sequence"] if std else ""
        reasons = []
        if not std or not dev:
            reasons.append("missing_standard_record" if not std else "missing_developer_record")
        if key[0].startswith("__MALFORMED__"):
            reasons.append("malformed_header")
        elif len(sh_accessions[key[0]]) != 1:
            reasons.append("species_hypothesis_accession_conflict")
        for kind, entries in (("standard", std), ("developer", dev)):
            if len(entries) > 1:
                conflict = len({(r["header"], r["sequence"]) for r in entries}) > 1
                if conflict:
                    reasons.append(kind + "_identity_conflict")
                for entry in entries[1:]:
                    issues.append(dict(source=kind, ordinal=entry["ordinal"],
                        issue="conflicting_duplicate" if conflict else "identical_duplicate",
                        record_id=row["record_id"], source_header=entry["header"],
                        sequence_sha256=sequence_sha(entry["sequence"])))
        if std and dev and std[0]["metadata"] != dev[0]["metadata"]:
            reasons.append("standard_developer_header_conflict")
        if row.get("kingdom") != "Fungi":
            reasons.append("nonfungal_or_unresolved_kingdom")
        if not sequence or set(sequence) - ALPHABET:
            reasons.append("invalid_developer_sequence_alphabet")
        if not standard_sequence or set(standard_sequence) - ALPHABET:
            reasons.append("invalid_standard_sequence_alphabet")
        relationship = "unverified"
        if sequence and standard_sequence:
            if standard_sequence in sequence:
                relationship = "standard_embedded_in_developer"
            elif standard_sequence.translate(COMPLEMENT)[::-1] in sequence:
                relationship = "standard_reverse_complement_embedded_in_developer"
            if relationship == "unverified":
                reasons.append("unverified_standard_developer_sequence_relationship")
        row.update(raw_length_bp=len(sequence), standard_length_bp=len(standard_sequence),
            sequence_sha256=sequence_sha(sequence), standard_sequence_sha256=sequence_sha(standard_sequence),
            sequence_alphabet="".join(sorted(set(sequence))),
            standard_sequence_alphabet="".join(sorted(set(standard_sequence))),
            source_header=std[0]["header"] if std else "", source_developer_header=dev[0]["header"] if dev else "",
            standard_record_count=len(std), developer_record_count=len(dev),
            standard_record_ordinals="|".join(str(r["ordinal"]) for r in std),
            developer_record_ordinals="|".join(str(r["ordinal"]) for r in dev),
            sequence_relationship=relationship, ambiguity_count=sum(base not in "ACGT" for base in sequence),
            ambiguity_fraction=(sum(base not in "ACGT" for base in sequence) / len(sequence)) if sequence else 1,
            length_stratum=length_stratum(len(sequence)), eligibility_status="eligible" if not reasons else "withheld",
            eligibility_reason=" | ".join(reasons))
        for column in RANKS:
            row[column + "_status"] = rank_status(row.get(column, ""))
        for rank in COVERAGE_RANKS:
            row[rank + "_group_id"] = group_identity(row, rank)
        if reasons:
            issues.append(dict(source="paired", record_id=row["record_id"], issue="withheld_reference",
                reason=row["eligibility_reason"], source_header=row["source_header"],
                source_developer_header=row["source_developer_header"], sequence_sha256=row["sequence_sha256"]))
        sequences[row["record_id"]] = sequence
        rows.append(row)
    # The same raw sequence assigned to different orders is not silently
    # treated as independent support for both taxonomic paths.
    sequence_groups = collections.defaultdict(list)
    for row in rows:
        if row["eligibility_status"] == "eligible":
            sequence_groups[row["sequence_sha256"]].append(row)
    for duplicates in sequence_groups.values():
        if len({r["order_group_id"] for r in duplicates}) > 1:
            for row in duplicates:
                row["eligibility_status"] = "withheld"
                row["eligibility_reason"] = "cross_order_identical_sequence_conflict"
                issues.append(dict(source="paired", record_id=row["record_id"], issue="withheld_reference",
                    reason=row["eligibility_reason"], source_header=row["source_header"],
                    sequence_sha256=row["sequence_sha256"]))
    return rows, sequences, issues, dict(standard_records=n_standard, developer_records=n_developer)


def priority(row):
    return (not row["representative_type"].lower().startswith("refs"),
            row["ambiguity_fraction"], hashlib.sha256(row["record_id"].encode()).hexdigest())


def select_group(rows, cap):
    unique = {}
    for row in sorted(rows, key=priority):
        unique.setdefault(row["sequence_sha256"], row)
    buckets = {stratum: [] for stratum in STRATA}
    for row in unique.values():
        buckets[row["length_stratum"]].append(row)
    for bucket in buckets.values():
        bucket.sort(key=priority)
    chosen, seen_genera = [], set()
    while len(chosen) < min(cap, len(unique)):
        for stratum in STRATA:
            bucket = buckets[stratum]
            if not bucket or len(chosen) >= cap:
                continue
            best_reference_flag = priority(bucket[0])[0]
            index = next((i for i, r in enumerate(bucket)
                          if priority(r)[0] == best_reference_flag and r["genus_status"] == "resolved"
                          and r["genus"] not in seen_genera), 0)
            picked = bucket.pop(index)
            chosen.append(picked)
            if picked["genus_status"] == "resolved":
                seen_genera.add(picked["genus"])
    return chosen, len(unique)


def choose_panel(rows, per_group):
    groups = collections.defaultdict(list)
    for row in rows:
        if row["eligibility_status"] == "eligible":
            groups[row["order_group_id"]].append(row)
    selected, eligible_counts = [], {}
    for group_id in sorted(groups):
        chosen, eligible_counts[group_id] = select_group(groups[group_id], per_group)
        selected.extend(chosen)
    return sorted(selected, key=lambda row: (row["order_group_id"], row["record_id"])), eligible_counts


def coverage_rows(census, selected, per_group):
    result = []
    for rank in COVERAGE_RANKS:
        column = RANK_COLUMN[rank]
        groups = collections.defaultdict(list)
        selected_groups = collections.Counter(r[rank + "_group_id"] for r in selected)
        for row in census:
            groups[row[rank + "_group_id"]].append(row)
        ancestors = list(RANK_COLUMN.values())[:list(RANK_COLUMN).index(rank)]
        for group_id in sorted(groups):
            rows = groups[group_id]
            row = rows[0]
            eligible = [r for r in rows if r["eligibility_status"] == "eligible"]
            taxon_status = rank_status(row.get(column, ""))
            result.append(dict(rank=rank, taxon=row.get(column, "") or "unidentified",
                parent_path=" › ".join(row.get(c, "") or "unidentified" for c in ancestors),
                group_id=group_id, n_source=len(rows),
                n_eligible_unique=len({r["sequence_sha256"] for r in eligible}),
                n_selected=selected_groups[group_id], n_withheld=sum(r["eligibility_status"] != "eligible" for r in rows),
                target_per_order_group=per_group if rank == "order" else "",
                shortfall_to_100=max(0, per_group - selected_groups[group_id]) if rank == "order" else "",
                taxon_resolved=taxon_status == "resolved", taxon_status=taxon_status))
    return result


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    path = Path(path)
    if path.suffix == ".gz":
        with path.open("wb") as binary, gzip.GzipFile(filename="", fileobj=binary, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, extrasaction="ignore")
                writer.writeheader(); writer.writerows(rows)
    else:
        with path.open("w", encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(text, fieldnames=fields, extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)


def validate_bundle(archive, standard, developer, expected, expected_records):
    paths = dict(archive=Path(archive), standard=Path(standard), developer=Path(developer))
    hashes = {key: sha_file(path) for key, path in paths.items()}
    for key, wanted in expected.items():
        if hashes[key] != wanted:
            raise ValueError("Pinned " + key + " source checksum differs")
    with tarfile.open(archive, "r:gz") as package:
        for key in ("standard", "developer"):
            member = package.getmember(paths[key].name)
            if not member.isfile():
                raise ValueError("Source archive member is not a regular file")
            digest = hashlib.sha256()
            with package.extractfile(member) as stream:
                while chunk := stream.read(1024 * 1024):
                    digest.update(chunk)
            if digest.hexdigest() != hashes[key]:
                raise ValueError("Cached " + key + " FASTA is not the pinned archive member")
    return hashes


def build_panel(archive, standard, developer, output, per_group=100, expected=None, expected_records=102137,
                source_doi=SOURCE_DOI, source_version=SOURCE_VERSION, source_url=SOURCE_URL):
    if not isinstance(per_group, int) or per_group < 1:
        raise ValueError("Per-order target must be a positive integer")
    hashes = validate_bundle(archive, standard, developer, expected or PINNED_HASHES, expected_records)
    identity = dict(schema_version="1.0.0", source_doi=source_doi, source_version=source_version,
        input_sha256=hashes, generator_sha256=sha_file(__file__), per_order_group=per_group,
        expected_source_records=expected_records, selection_policy="v1_reference_genus_ambiguity_length_round_robin")
    version_id = "unite-" + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
    output = Path(output)
    if output.exists():
        manifest_path = output / "source_manifest.json"
        if not manifest_path.exists():
            raise ValueError("Existing output has no immutable panel manifest")
        existing = json.loads(manifest_path.read_text())
        if existing.get("version_id") != version_id:
            raise ValueError("Different immutable panel already exists; choose a new output directory")
        for relative, descriptor in existing["artifacts"].items():
            if sha_file(output / relative) != descriptor["sha256"]:
                raise ValueError("Existing panel artifact checksum differs: " + relative)
        return existing
    census, sequences, issues, counts = pair_sources(standard, developer)
    if any(counts[key] != expected_records for key in ("standard_records", "developer_records")):
        raise ValueError("Source record count differs from the pinned release")
    selected, eligible_counts = choose_panel(census, per_group)
    if not selected:
        raise ValueError("No source-linked developer references pass import QC")
    if len({r["record_id"] for r in selected}) != len(selected) or len({r["species_hypothesis"] for r in selected}) != len(selected):
        raise ValueError("Selected SH/reference identities are not unique")
    if len({r["sequence_sha256"] for r in selected}) != len(selected):
        raise ValueError("Selected raw sequence identities are not unique")
    for row in census:
        row.update(source_version=source_version, source_doi=source_doi, panel_version_id=version_id)
    coverage = coverage_rows(census, selected, per_group)
    order_rows = [r for r in coverage if r["rank"] == "order"]
    counts.update(source_reference_identities=len(census), selected_references=len(selected),
        eligible_reference_identities=sum(r["eligibility_status"] == "eligible" for r in census),
        withheld_reference_identities=sum(r["eligibility_status"] != "eligible" for r in census),
        resolved_order_groups=sum(r["taxon_resolved"] for r in order_rows),
        unresolved_order_groups=sum(not r["taxon_resolved"] for r in order_rows),
        order_groups_with_selected=sum(r["n_selected"] > 0 for r in order_rows),
        order_groups_below_target=sum(r["n_selected"] < per_group for r in order_rows),
        standard_embedded_in_developer=sum(r["sequence_relationship"] == "standard_embedded_in_developer" for r in census),
        developer_longer_than_standard=sum(r["raw_length_bp"] > r["standard_length_bp"] for r in census),
        developer_same_length_as_standard=sum(r["raw_length_bp"] == r["standard_length_bp"] for r in census),
        developer_shorter_than_standard=sum(r["raw_length_bp"] < r["standard_length_bp"] for r in census))
    manifest = dict(identity, version_id=version_id, source_url=source_url,
        archive_member_names=dict(standard=Path(standard).name, developer=Path(developer).name),
        built_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), counts=counts,
        taxonomy_source="UNITE release header; retained without translating to NCBI taxonomy",
        identifier_basis="Species hypothesis plus accession token exactly as published; no accession version inferred",
        sequence_basis="Developer FASTA bases exactly as published after uppercase/FASTA-line normalization",
        source_census_basis="One row per paired SH/accession identity; duplicate input occurrences retained as counts/ordinals/issues",
        selection_policy="Up to target per resolved or unresolved hierarchy-scoped order group; retain every rare eligible reference; unique SH and raw sequence; round-robin all length strata, source RefS flag then novel resolved genus then lower ambiguity then stable hash",
        limitations=["An ITS record's length does not establish complete 18S, ITS1, 5.8S, ITS2 or 28S coverage.",
          "Selected references are deliberately balanced, not a representative population sample or amplification validation.",
          "Source taxonomy placeholders and Incertae sedis remain distinct from resolved taxa.",
          "Unverified source pair relationships and conflicting identities are withheld, never silently reconciled."],
        issue_counts=dict(collections.Counter(r["issue"] for r in issues)))
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".unite-panel-", dir=output.parent))
    try:
        write_csv(staging / "source_census.csv.gz", census)
        write_csv(staging / "selected_metadata.csv.gz", selected)
        write_csv(staging / "source_coverage.csv", coverage)
        write_csv(staging / "import_issues.csv.gz", issues,
                  list(dict.fromkeys(["source", "ordinal", "record_id", "issue", "reason", "source_header",
                                     "source_developer_header", "sequence_sha256"])))
        with (staging / "selected_raw.fasta.gz").open("wb") as raw, gzip.GzipFile(
                filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            for row in selected:
                compressed.write((">" + row["record_id"] + "\n" + sequences[row["record_id"]] + "\n").encode())
        shutil.copyfile(__file__, staging / "code_unite_panel_builder.py")
        manifest["artifacts"] = {path.name: dict(sha256=sha_file(path), bytes=path.stat().st_size)
                                 for path in sorted(staging.iterdir())}
        for name, count in {"source_census.csv.gz": len(census), "selected_metadata.csv.gz": len(selected),
                            "source_coverage.csv": len(coverage), "import_issues.csv.gz": len(issues),
                            "selected_raw.fasta.gz": len(selected)}.items():
            manifest["artifacts"][name]["rows"] = count
        (staging / "source_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        if output.exists():
            raise ValueError("Output appeared during build; refusing to replace it")
        staging.rename(output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    cache = ROOT / "data/external/its_fungal_pilot"
    parser.add_argument("--archive", type=Path, default=cache / "sh_general_release_19.02.2025.tgz")
    parser.add_argument("--standard-fasta", type=Path, default=cache / STANDARD_MEMBER)
    parser.add_argument("--developer-fasta", type=Path, default=cache / DEVELOPER_MEMBER)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/derived/unite_reference_panel")
    parser.add_argument("--per-order", type=int, default=100)
    args = parser.parse_args()
    result = build_panel(args.archive, args.standard_fasta, args.developer_fasta, args.output_dir, args.per_order)
    print(json.dumps(dict(version_id=result["version_id"], counts=result["counts"]), indent=2))


if __name__ == "__main__":
    main()
