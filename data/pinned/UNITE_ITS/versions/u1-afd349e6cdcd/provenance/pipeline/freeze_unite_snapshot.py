#!/usr/bin/env python3
"""Validate and append an immutable broad-ITS UNITE snapshot.

No acquisition, annotation, site search, or alignment is performed here. A
sparse missing-site row can be interpreted only after its exhaustive search
proof is tied to these exact sequence/catalog bytes. SSU/LSU annotations are
profile-detected retained flanks, never evidence of complete rRNA genes.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import datetime
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import statistics
import tempfile

from build_unite_reference_panel import group_identity, parse_header, rank_status
from scan_unite_primer_sites import (SCAN_RULE_VERSION, exact_site, pair_result,
                                     prepare_primer, read_csv, read_fasta,
                                     record_ids_sha256, reverse_complement,
                                     validate_raw_records, write_gzip_csv)

SCHEMA_VERSION = "1.0.0"
FREEZE_RULE_VERSION = "unite-its-validated-immutable-v1"
REGIONS = ("SSU", "ITS1", "5.8S", "ITS2", "LSU")
RANK_COLUMNS = dict(phylum="phylum", **{"class": "class_name"},
                    order="order_name", family="family")
HASH = re.compile(r"[0-9a-f]{64}\Z")
SOURCE_LICENSE = dict(name="Creative Commons Attribution 4.0 International",
    spdx="CC-BY-4.0", url="https://creativecommons.org/licenses/by/4.0/",
    declared_at="https://unite.ut.ee/cite.php", verified_on="2026-10-13")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bases_sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def canonical_sha(value):
    return bases_sha(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False))


def load_json(path):
    with Path(path).open() as handle:
        return json.load(handle)


def safe_path(root, relative):
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise ValueError("Unsafe artifact path")
    value = Path(relative)
    if value.is_absolute() or any(part in (".", "..") for part in value.parts):
        raise ValueError("Unsafe artifact path: " + relative)
    result = Path(root) / value
    if result.is_symlink() or not result.resolve().is_relative_to(Path(root).resolve()):
        raise ValueError("Artifact escapes its declared owner")
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, name, optional=False):
    if optional and value in (None, ""):
        return None
    require(not isinstance(value, bool) and re.fullmatch(r"\d+", str(value)) is not None,
            "Invalid integer: " + name)
    return int(value)


def verify_artifacts(root, descriptors):
    require(isinstance(descriptors, dict) and bool(descriptors), "No declared artifacts")
    result = {}
    paths = set()
    for name, item in descriptors.items():
        relative = item.get("path", name)
        path = safe_path(root, relative)
        require(relative not in paths, "Repeated declared artifact path")
        paths.add(relative)
        require(path.is_file() and HASH.fullmatch(item.get("sha256", "")) is not None,
                "Missing artifact or invalid hash: " + relative)
        require(sha(path) == item["sha256"], "Artifact checksum differs: " + relative)
        if "bytes" in item:
            require(path.stat().st_size == integer(item["bytes"], relative), "Artifact size differs: " + relative)
        result[name] = path
    return result


def unique_rows(rows, key, name):
    values = [row[key] for row in rows]
    require(all(values) and len(set(values)) == len(values), name + " identities repeat or are absent")
    return dict(zip(values, rows))


def raw_records(path):
    records = list(read_fasta(path))
    validate_raw_records(records)
    return dict(records)


def validate_source(panel):
    manifest = load_json(panel / "source_manifest.json")
    verify_artifacts(panel, manifest["artifacts"])
    metadata = read_csv(panel / "selected_metadata.csv.gz")
    by_id = unique_rows(metadata, "record_id", "Selected metadata")
    sequences = raw_records(panel / "selected_raw.fasta.gz")
    require(set(sequences) == set(by_id), "Source metadata/sequence membership differs")
    require(len(by_id) == manifest["counts"]["selected_references"], "Selected source count differs")
    for key, row in by_id.items():
        require(bases_sha(sequences[key]) == row["sequence_sha256"] and
                len(sequences[key]) == integer(row["raw_length_bp"], "raw_length_bp"), "Original source bases differ")
        require(row["accession"] == row["accession_raw"], "Published accession token was changed")
        header = parse_header(row["source_developer_header"])
        require(row["eligibility_status"] == "eligible" and all(row[field] == value for field, value in header.items()),
                "Selected source header identity/taxonomy differs")
        for rank, column in RANK_COLUMNS.items():
            require(row[rank + "_group_id"] == group_identity(row, rank), "Hierarchy identity differs: " + rank)
            require(row[column + "_status"] == rank_status(row[column]), "Taxonomy status differs: " + rank)
    coverage = read_csv(panel / "source_coverage.csv")
    unique_rows(coverage, "group_id", "Coverage")
    census = read_csv(panel / "source_census.csv.gz")
    census_ids = unique_rows(census, "record_id", "Source census")
    require(len(census) == manifest["counts"]["source_reference_identities"], "Source census count differs")
    for key, row in by_id.items():
        require(key in census_ids and all(census_ids[key].get(field) == value for field, value in row.items()),
                "Selected metadata differs from source census")
    for rank, column in RANK_COLUMNS.items():
        groups = {row["group_id"]: row for row in coverage if row["rank"] == rank}
        source_groups, chosen_groups = defaultdict(list), Counter(row[rank + "_group_id"] for row in metadata)
        for row in census:
            require(row[rank + "_group_id"] == group_identity(row, rank), "Census hierarchy identity differs")
            source_groups[row[rank + "_group_id"]].append(row)
        require(set(groups) == set(source_groups), "Coverage omits or invents a hierarchy group")
        for key, rows in source_groups.items():
            row = groups[key]
            n_eligible = len({item["sequence_sha256"] for item in rows if item["eligibility_status"] == "eligible"})
            checks = dict(n_source=len(rows), n_eligible_unique=n_eligible,
                          n_selected=chosen_groups[key], n_withheld=sum(item["eligibility_status"] != "eligible" for item in rows))
            require(all(integer(row[field], field) == value for field, value in checks.items()), "Coverage denominator differs")
            ancestors = list(RANK_COLUMNS.values())[:list(RANK_COLUMNS).index(rank)]
            require(row["taxon"] == (rows[0][column] or "unidentified") and row["taxon_status"] == rank_status(rows[0][column]) and
                    row["parent_path"] == " › ".join(rows[0][field] or "unidentified" for field in ancestors),
                    "Coverage taxonomy label/status differs")
            require(row["taxon_resolved"].lower() in ("true", "false") and
                    (row["taxon_resolved"].lower() == "true") == (row["taxon_status"] == "resolved"), "Resolved flag differs")
            if rank == "order":
                target = integer(manifest["per_order_group"], "per_order_group")
                require(integer(row["target_per_order_group"], "target") == target and
                        integer(row["shortfall_to_100"], "shortfall") == max(0, target - chosen_groups[key]) and
                        chosen_groups[key] == min(target, n_eligible), "Order sampling/shortage differs")
    return manifest, metadata, sequences, coverage


def validate_annotations(panel, annotations, metadata, source_sequences):
    manifest = load_json(annotations / "annotation_manifest.json")
    declared = verify_artifacts(annotations, manifest["artifacts"])
    require(manifest.get("parser_rule_version") == "itsx_regions_verified_oriented_slices_v2",
            "Unsupported or unverified ITSx parser rules")
    require(manifest["identity"]["input_sha256"] == sha(panel / "selected_raw.fasta.gz") and
            manifest["identity"]["metadata_sha256"] == sha(panel / "selected_metadata.csv.gz") and
            manifest["input_sha256"] == sha(panel / "selected_raw.fasta.gz"), "Annotation source identity differs")
    require(manifest["n_input"] == len(metadata), "Annotation input count differs")
    require(sha(annotations / "code_region_parser.py") == manifest["script_sha256"] == manifest["identity"]["script_sha256"] and
            sha(annotations / "code_parser_scanner_dependency.py") == manifest["helper_sha256"] == manifest["identity"]["helper_sha256"],
            "Annotation executed parser/dependency differs")
    tool_path, run_path = annotations / "tool_identity_manifest.json", annotations / "profile_run_identity.json"
    require(sha(tool_path) == manifest["tool_manifest_sha256"] == manifest["identity"]["tool_manifest_sha256"] and
            sha(run_path) == manifest["profile_run_manifest_sha256"], "Profile proof hash differs")
    tool, run = load_json(tool_path), load_json(run_path)
    require(tool.get("itsx_version") == "1.1.3" and tool.get("hmmer_version") == "3.4", "Unsupported ITSx/HMMER tool identity")
    require({"itsx_program", "itsx_wrapper", "hmmsearch", "fungal_profile"} <= set(tool.get("files", {})), "Incomplete tool identity")
    require(all(HASH.fullmatch(item.get("sha256", "")) for item in tool["files"].values()), "Invalid tool file checksum")
    require(run.get("schema_version") == "1.0.0" and run.get("capture_mode") in
            ("before_profile_execution", "active_process_verified_before_edit"), "Unverified profile run identity")
    staged = "".join(f">{key}\n{value}\n" for key, value in source_sequences.items())
    require(run.get("input_sha256") == manifest["input_sha256"] and run.get("staged_input_sha256") == bases_sha(staged) and
            run.get("n_input") == len(metadata) and run.get("tool_manifest_sha256") == sha(tool_path) and
            run.get("semantic_parameters") == manifest["identity"]["semantic_parameters"], "Profile run input/tool parameters differ")
    for name in ("annotator", "scanner"):
        require(sha(annotations / ("code_profile_run_" + name + ".py")) == run["executed_source"].get(name + "_sha256"),
                "Captured profile-run code differs")
    orientations = unique_rows(read_csv(annotations / "orientations.csv.gz"), "record_id", "Orientation")
    sequences = raw_records(annotations / "oriented_raw.fasta.gz")
    require(set(orientations) == set(source_sequences) == set(sequences), "Annotation oriented membership differs")
    merged = []
    for row in metadata:
        key = row["record_id"]; orientation = orientations[key]
        require(orientation["source_strand"] in ("1", "-1"), "Invalid source strand")
        expected = reverse_complement(source_sequences[key]) if orientation["source_strand"] == "-1" else source_sequences[key]
        require(sequences[key] == expected and orientation["source_sequence_sha256"] == row["sequence_sha256"] and
                orientation["analysis_sequence_sha256"] == bases_sha(expected), "Oriented source bases differ")
        require(orientation["orientation"] in ("profile_forward", "source_unverified"), "Unsupported orientation evidence")
        if orientation["orientation"] == "source_unverified":
            require(orientation["source_strand"] == "1" and orientation["annotation_status"] == "not_detected",
                    "An unverified source orientation was inferred")
        merged.append(dict(row, **{field: value for field, value in orientation.items() if field != "record_id"},
                           organism=row["species"].replace("_", " ")))
    require(sum(row["source_strand"] == "-1" for row in orientations.values()) == manifest["n_reverse"] and
            sum(row["annotation_status"] == "chimeric_withheld" for row in orientations.values()) == manifest["n_chimeric"],
            "Annotation orientation/status counts differ")
    detected = sum(row["orientation"] == "profile_forward" for row in orientations.values())
    require(detected == manifest["n_detected"] == manifest["summary_proof"]["real_detected_all_batches"],
            "Owned profile-call completion count differs")
    regions = read_csv(annotations / "regions.csv.gz")
    keys = [(row["record_id"], row["region"]) for row in regions]
    require(len(set(keys)) == len(keys) and set(keys) == {(key, region) for key in sequences for region in REGIONS},
            "Region evidence must retain every selected identity and all five region states")
    counts = Counter(row["region"] + ":" + row["boundary_status"] for row in regions)
    require(dict(counts) == manifest["region_counts"], "Region state counts differ")
    by_key = dict(zip(keys, regions))
    for row in regions:
        key, region, status = row["record_id"], row["region"], row["boundary_status"]
        start, end = integer(row["start"], "region start", True), integer(row["end"], "region end", True)
        require((start is None) == (end is None), "Incomplete region coordinate pair")
        valid = start is not None and 1 <= start <= end <= len(sequences[key])
        require(status in ("complete", "partial", "flank_detected", "not_detected", "invalid_coordinates",
                           "chimeric_withheld", "extraction_mismatch"), "Unsupported boundary status")
        require(row["source_strand"] == orientations[key]["source_strand"] and
                row["annotation_status"] == orientations[key]["annotation_status"] and
                row["coordinate_basis"] == "1-based inclusive positions in oriented source sequence", "Region coordinate/strand evidence differs")
        if status in ("complete", "partial", "flank_detected"):
            require(valid, "Accepted region bounds leave the selected source")
        if status == "not_detected":
            require(start is None, "Undetected region was assigned coordinates")
        original = (integer(row["original_start"], "original start", True), integer(row["original_end"], "original end", True))
        expected_original = ((len(sequences[key]) - end + 1, len(sequences[key]) - start + 1)
                             if valid and row["source_strand"] == "-1" else (start, end) if valid else (None, None))
        require(original == expected_original, "Original/oriented region projection differs")
        if region in ("SSU", "LSU"):
            require(status != "complete", "A retained rRNA flank was labelled a complete gene")
        if status == "complete":
            require(orientations[key]["annotation_status"] == "profile_detected", "Withheld profile evidence became complete")
            left, right = {"ITS1": ("SSU", "5.8S"), "5.8S": ("ITS1", "ITS2"), "ITS2": ("5.8S", "LSU")}[region]
            lrow, rrow = by_key[(key, left)], by_key[(key, right)]
            require(integer(lrow["end"], "left anchor", True) == start - 1 and
                    integer(rrow["start"], "right anchor", True) == end + 1 and start > 1 and end < len(sequences[key]),
                    "Complete region lacks both contiguous anchors")
    for region in ("ITS1", "5.8S", "ITS2"):
        extracted_rows = list(read_fasta(annotations / (region.replace(".", "_") + ".fasta.gz")))
        require(len({key for key, _ in extracted_rows}) == len(extracted_rows), "Complete region extraction repeats an identity")
        extracted = dict(extracted_rows)
        complete = {row["record_id"]: sequences[row["record_id"]][int(row["start"]) - 1:int(row["end"])]
                    for row in regions if row["region"] == region and row["boundary_status"] == "complete"}
        require(extracted == complete, "Complete region extraction/source slice differs")
    # The original ITSx output additionally proves accepted partial and flank
    # intervals; withheld extraction mismatches remain present as audit evidence.
    for region in REGIONS:
        path = annotations / ("regions." + region.replace(".", "_") + ".fasta.gz")
        extracted_rows = list(read_fasta(path))
        require(len({key for key, _ in extracted_rows}) == len(extracted_rows) and
                {key for key, _ in extracted_rows} <= set(sequences), "Original ITSx extraction membership differs")
        extracted = dict(extracted_rows)
        for row in regions:
            if row["region"] == region and row["boundary_status"] in ("complete", "partial", "flank_detected"):
                expected = sequences[row["record_id"]][int(row["start"]) - 1:int(row["end"])]
                require(extracted.get(row["record_id"]) == expected, "Accepted original ITSx extraction/source slice differs")
    return manifest, merged, sequences, regions, declared


def validate_search(panel, evidence, sequence_path, metadata, sequences):
    proof = load_json(evidence / "site_search_manifest.json")
    catalog, pairs, sites = (read_csv(evidence / name) for name in ("catalog.csv.gz", "pairs.csv.gz", "primer_sites.csv.gz"))
    by_catalog = unique_rows(catalog, "catalog_id", "Primer catalog")
    unique_rows(pairs, "pair_id", "Pair catalog")
    require(bool(catalog) and bool(pairs), "An empty catalog does not establish primer evidence")
    checks = dict(sequence_sha256=sha(sequence_path), metadata_sha256=sha(panel / "selected_metadata.csv.gz"),
                  record_ids_sha256=record_ids_sha256(sequences), n_records=len(sequences), n_catalog=len(catalog), n_pairs=len(pairs),
                  n_site_searches=len(sequences) * len(catalog), n_sparse_rows=len(sites),
                  catalog_artifact_sha256=sha(evidence / "catalog.csv.gz"), pair_artifact_sha256=sha(evidence / "pairs.csv.gz"),
                  site_artifact_sha256=sha(evidence / "primer_sites.csv.gz"), script_sha256=sha(evidence / "code_scanner.py"))
    require(proof.get("complete") is True and proof.get("scan_rule_version") == SCAN_RULE_VERSION and
            all(proof.get(field) == expected for field, expected in checks.items()), "Exhaustive site-search proof differs")
    require(integer(proof.get("n_no_exact_hits"), "n_no_exact_hits") + len(sites) == proof["n_site_searches"],
            "Sparse search proof does not account for every primer/template search")
    components = {(row["primer_name"], row["sequence"], row["direction"]): row["catalog_id"] for row in catalog}
    require(len(components) == len(catalog) and all(row["direction"] in ("forward", "reverse") for row in catalog),
            "Catalog component identity/direction is ambiguous")
    prepared = {key: prepare_primer(row["sequence"], row["direction"]) for key, row in by_catalog.items()}
    sparse = {}
    for row in sites:
        key = (row["record_id"], row["catalog_id"])
        require(key not in sparse and key[0] in sequences and key[1] in by_catalog, "Sparse site key/ownership differs")
        require(row["primer_name"] == by_catalog[key[1]]["primer_name"] and row["direction"] == by_catalog[key[1]]["direction"],
                "Sparse site primer identity differs")
        expected = exact_site(sequences[key[0]], prepared[key[1]])
        actual = dict(start=integer(row["start"], "site start", True), end=integer(row["end"], "site end", True),
                      n_matches=integer(row["n_matches"], "site hits"), status=row["status"])
        require(expected == actual and actual["status"] != "no_exact_binding_site", "Sparse exact-site evidence does not match source bases")
        sparse[key] = actual
    # Reconcile every pair's denominator, statuses and length summary against
    # these sparse site rows. Omitted hits are permitted only by the proof above.
    absent = dict(start=None, end=None, n_matches=0, status="no_exact_binding_site")
    states = ("scorable", "both_sites_missing", "reverse_site_missing", "forward_site_missing", "ambiguous_binding_site",
              "low_specificity_site", "outward_or_overlapping_sites", "below_minimum_length")
    for row in pairs:
        fkey = (row["forward_primer"], row["forward_sequence"], "forward")
        rkey = (row["reverse_primer"], row["reverse_sequence"], "reverse")
        require(fkey in components and rkey in components, "Pair component is outside the catalog")
        counts, lengths = Counter(), []
        for record_id in sequences:
            state, bp = pair_result(sparse.get((record_id, components[fkey]), absent), sparse.get((record_id, components[rkey]), absent))
            counts[state] += 1
            if bp is not None: lengths.append(bp)
        require(integer(row["n_panel"], "pair n_panel") == len(sequences) and
                integer(row["n_scorable"], "pair n_scorable") == len(lengths) and
                all(integer(row["n_" + state], state) == counts[state] for state in states), "Pair support counts differ")
        for field, expected in (("amplicon_min_bp", min(lengths) if lengths else None),
                                ("amplicon_median_bp", statistics.median(lengths) if lengths else None),
                                ("amplicon_max_bp", max(lengths) if lengths else None)):
            actual = None if row[field] == "" else float(row[field])
            require(actual == expected, "Pair encoded length summary differs")
    return proof, catalog, pairs, sites


def expected_exemplars(rows, limit):
    rows = sorted(rows, key=lambda row: (float(row.get("ambiguity_fraction", 0)), row["representative_type"] != "refs", row["record_id"]))
    selected, used = [], set()
    for column in ("order_group_id", "genus", None):
        seen = set()
        for row in rows:
            if row["record_id"] in used or (column and row[column] in seen): continue
            selected.append(row["record_id"]); used.add(row["record_id"])
            if column: seen.add(row[column])
            if len(selected) >= limit: return set(selected)
    return set(selected)


def validate_alignments(panel, annotations, alignments, metadata, sequences, regions):
    manifest = load_json(alignments / "alignment_manifest.json")
    declared = verify_artifacts(alignments, manifest["artifacts"])
    identity = manifest["identity"]
    require(identity["input_sha256"] == sha(panel / "selected_metadata.csv.gz") and
            identity["regions_sha256"] == sha(annotations / "regions.csv.gz") and
            identity["oriented_sha256"] == sha(annotations / "oriented_raw.fasta.gz"), "Regional alignment input identity differs")
    require(isinstance(identity.get("mafft_version"), str) and bool(identity["mafft_version"].strip()), "MAFFT version is unrecorded")
    require(identity["script_sha256"] == sha(Path(__file__).with_name("build_unite_regional_alignments.py")),
            "Regional alignment generator differs from executed source")
    limit = integer(identity["limit"], "alignment limit")
    require(limit > 0, "Alignment limit must be positive")
    index = read_csv(alignments / "alignment_index.csv.gz")
    by_id = {row["record_id"]: row for row in metadata}
    by_region, expected_groups = {}, defaultdict(list)
    for row in regions:
        if row["region"] in ("ITS1", "5.8S", "ITS2") and row["boundary_status"] == "complete":
            by_region[(row["record_id"], row["region"])] = row
            meta = by_id[row["record_id"]]
            for rank in ("class", "order"):
                if meta[RANK_COLUMNS[rank] + "_status"] == "resolved":
                    expected_groups[(rank, meta[rank + "_group_id"], row["region"])].append(meta)
    keys = [(row["rank"], row["group_id"], row["region"]) for row in index]
    require(len(keys) == len(set(keys)) and set(keys) == set(expected_groups), "Regional alignment groups/regions differ")
    paths = set()
    for row in index:
        key = (row["rank"], row["group_id"], row["region"])
        relative = row["file"]
        require(relative.startswith("alignments/") and relative.endswith(".fasta.gz") and relative not in paths,
                "Alignment file path/identity differs")
        paths.add(relative)
        path = safe_path(alignments, relative)
        require(relative in declared, "Alignment file has no artifact proof")
        records = list(read_fasta(path)); ids = [name for name, _ in records]
        widths = {len(value) for _, value in records}
        expected = expected_exemplars(expected_groups[key], limit)
        require(len(ids) == len(set(ids)) and set(ids) == expected and len(widths) == 1 and all(widths),
                "Alignment membership/width differs from bounded regional exemplars")
        require(integer(row["n_source_regions"], "n_source_regions") == len(expected_groups[key]) and
                integer(row["n_input"], "n_input") == len(ids) and integer(row["n_columns"], "n_columns") == next(iter(widths)),
                "Alignment index counts differ")
        for record_id, value in records:
            region = by_region[(record_id, row["region"])]
            require(by_id[record_id][row["rank"] + "_group_id"] == row["group_id"] and
                    set(value) <= set("ACGTRYSWKMBDHVN-") and
                    value.replace("-", "") == sequences[record_id][int(region["start"]) - 1:int(region["end"])],
                    "Alignment source base retention/group ownership differs")
        require(integer(row["n_gap_cells"], "n_gap_cells") == sum(value.count("-") for _, value in records) and
                integer(row["n_ambiguous_cells"], "n_ambiguous_cells") == sum(sum(base not in "ACGT-" for base in value) for _, value in records),
                "Alignment encoded gap/ambiguity counts differ")
    require(manifest["n_alignments"] == len(index) and manifest["n_alignment_records"] == sum(int(row["n_input"]) for row in index),
            "Alignment manifest totals differ")
    # Undeclared output files could otherwise be accidentally shipped as proof.
    files = {str(path.relative_to(alignments)) for path in alignments.rglob("*") if path.is_file()}
    require(files == {str(path.relative_to(alignments)) for path in declared.values()} | {"alignment_manifest.json"},
            "Regional output contains undeclared artifacts")
    return manifest, index, declared


def atomic_json(path, value):
    path = Path(path)
    payload = (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode()
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=".pointer-", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name); handle.write(payload)
        temporary.replace(path)
    finally:
        if temporary and temporary.exists(): temporary.unlink()


@contextmanager
def exclusive_freeze(destination):
    lock = destination / ".freeze.lock"
    try: lock.mkdir()
    except FileExistsError: raise ValueError("Another snapshot freezer is active; investigate before retrying") from None
    try: yield
    finally: lock.rmdir()


def validate_frozen(directory, manifest):
    require(manifest.get("marker_id") == "UNITE_ITS" and manifest.get("schema_version") == SCHEMA_VERSION,
            "Historical UNITE manifest identity differs")
    verify_artifacts(directory, manifest["artifacts"])


def freeze_snapshot(panel, annotations, evidence, alignments, destination, label=None, max_bytes=70 * 1024 * 1024):
    panel, annotations, evidence, alignments, destination = map(Path, (panel, annotations, evidence, alignments, destination))
    # Validate every completed stage before creating a registry or staging output.
    source, source_meta, source_seq, coverage = validate_source(panel)
    annotation, metadata, sequences, regions, annotation_files = validate_annotations(panel, annotations, source_meta, source_seq)
    search, catalog, pairs, sites = validate_search(panel, evidence, annotations / "oriented_raw.fasta.gz", metadata, sequences)
    alignment, index, alignment_files = validate_alignments(panel, annotations, alignments, metadata, sequences, regions)
    inputs = dict(source_manifest=sha(panel / "source_manifest.json"), annotation_manifest=sha(annotations / "annotation_manifest.json"),
                  search_manifest=sha(evidence / "site_search_manifest.json"), alignment_manifest=sha(alignments / "alignment_manifest.json"))
    tool = load_json(annotations / "tool_identity_manifest.json")
    # Paths, timestamps, PIDs and MAFFT stderr are captured as provenance, but
    # cannot create another scientific identity for identical evidence bytes.
    identity_basis = dict(rule=FREEZE_RULE_VERSION, generator_sha256=sha(__file__),
        source=dict(source_version=source["source_version"], source_doi=source["source_doi"],
            input_sha256=source["input_sha256"], generator_sha256=source["generator_sha256"],
            per_order_group=source["per_order_group"],
            artifacts={name: item["sha256"] for name, item in source["artifacts"].items()}),
        annotation=dict(parser_rule_version=annotation["parser_rule_version"], semantic_parameters=annotation["identity"]["semantic_parameters"],
            script_sha256=annotation["script_sha256"], helper_sha256=annotation["helper_sha256"],
            tool_files={name: item["sha256"] for name, item in tool["files"].items()},
            data_sha256={name: sha(annotations / name) for name in ("oriented_raw.fasta.gz", "regions.csv.gz", "orientations.csv.gz")}),
        search=dict(scan_rule_version=search["scan_rule_version"], script_sha256=search["script_sha256"],
            artifacts={name: search[name] for name in ("sequence_sha256", "metadata_sha256", "catalog_artifact_sha256",
                                                       "pair_artifact_sha256", "site_artifact_sha256")}),
        alignment=dict(script_sha256=alignment["identity"]["script_sha256"], limit=alignment["identity"]["limit"],
            mafft_version=alignment["identity"]["mafft_version"], index_sha256=sha(alignments / "alignment_index.csv.gz"),
            alignments={name: item["sha256"] for name, item in alignment["artifacts"].items() if name.endswith(".fasta.gz")}))
    identity = canonical_sha(identity_basis)
    expected_hashes = {}
    for owner, descriptors in ((panel, source["artifacts"]), (annotations, annotation["artifacts"]), (alignments, alignment["artifacts"])):
        expected_hashes.update({safe_path(owner, item.get("path", name)).resolve(): item["sha256"] for name, item in descriptors.items()})
    expected_hashes.update({(evidence / name).resolve(): search[key] for name, key in (
        ("catalog.csv.gz", "catalog_artifact_sha256"), ("pairs.csv.gz", "pair_artifact_sha256"),
        ("primer_sites.csv.gz", "site_artifact_sha256"), ("code_scanner.py", "script_sha256"))})
    versions = destination / "versions"
    versions.mkdir(parents=True, exist_ok=True)
    with exclusive_freeze(destination):
        registry_path = destination / "index.json"
        registry = load_json(registry_path) if registry_path.exists() else dict(schema_version=SCHEMA_VERSION, marker_id="UNITE_ITS", versions=[])
        require(registry.get("marker_id") == "UNITE_ITS" and isinstance(registry.get("versions"), list), "Invalid historical registry")
        require(len({entry["id"] for entry in registry["versions"]}) == len(registry["versions"]), "Historical registry repeats identities")
        frozen = []
        for path in sorted(versions.glob("u*-*/manifest.json")):
            manifest = load_json(path); validate_frozen(path.parent, manifest)
            require(path.parent.name == manifest["id"] and re.fullmatch(r"u[1-9]\d*-[0-9a-f]{12}", manifest["id"]), "Invalid historical version path")
            frozen.append((path, manifest))
        by_id = {manifest["id"]: (path, manifest) for path, manifest in frozen}
        for entry in registry["versions"]:
            require(entry["id"] in by_id and safe_path(destination, entry["manifest"]) == by_id[entry["id"]][0] and
                    entry["sha256"] == sha(by_id[entry["id"]][0]), "Historical registry/proof differs")
        same = [(path, manifest) for path, manifest in frozen if manifest.get("identity_sha256") == identity]
        require(len(same) <= 1, "Immutable content identity repeats")
        if same:
            path, manifest = same[0]
            # Recover registration after a completed directory rename but an
            # interrupted pointer write. Existing registered history is untouched.
            if manifest["id"] not in {entry["id"] for entry in registry["versions"]}:
                registry["versions"].append(dict(id=manifest["id"], label=manifest["label"], manifest=str(path.relative_to(destination)),
                                                  sha256=sha(path), n_selected=manifest["n_selected"]))
                atomic_json(registry_path, registry)
                atomic_json(destination / "latest.json", registry["versions"][-1])
            return manifest
        ordinal = 1 + max((int(manifest["id"].split("-")[0][1:]) for _, manifest in frozen), default=0)
        version_id = f"u{ordinal}-{identity[:12]}"
        staging = Path(tempfile.mkdtemp(prefix=".snapshot-", dir=versions))
        try:
            artifacts = {}
            def copy(name, original, relative, rows=None):
                path = safe_path(staging, relative); path.parent.mkdir(parents=True, exist_ok=True)
                require(not path.exists(), "Snapshot artifact path repeats")
                shutil.copyfile(original, path)
                artifacts[name] = dict(path=relative, sha256=sha(path), bytes=path.stat().st_size)
                if original.resolve() in expected_hashes:
                    require(artifacts[name]["sha256"] == expected_hashes[original.resolve()], "Validated upstream bytes changed during copy")
                if rows is not None: artifacts[name]["rows"] = rows
            metadata_path = staging / "metadata.csv.gz"
            write_gzip_csv(metadata_path, list(metadata[0]), metadata)
            artifacts["metadata"] = dict(path=metadata_path.name, sha256=sha(metadata_path), bytes=metadata_path.stat().st_size, rows=len(metadata))
            main = dict(coverage=(panel / "source_coverage.csv", "source_coverage.csv", len(coverage)),
                sequences=(annotations / "oriented_raw.fasta.gz", "sequences.fasta.gz", len(sequences)),
                regions=(annotations / "regions.csv.gz", "regions.csv.gz", len(regions)),
                sites=(evidence / "primer_sites.csv.gz", "primer_sites.csv.gz", len(sites)),
                catalog=(evidence / "catalog.csv.gz", "catalog.csv.gz", len(catalog)),
                pairs=(evidence / "pairs.csv.gz", "pairs.csv.gz", len(pairs)),
                search_manifest=(evidence / "site_search_manifest.json", "site_search_manifest.json", None),
                alignment_index=(alignments / "alignment_index.csv.gz", "alignment_index.csv.gz", len(index)))
            original_destinations = {}
            for name, (original, relative, rows) in main.items():
                copy(name, original, relative, rows); original_destinations[original.resolve()] = relative
            for name in source["artifacts"]:
                copy("source_" + name.replace(".", "_"), panel / name, "provenance/source/" + name, source["artifacts"][name].get("rows"))
            copy("source_manifest", panel / "source_manifest.json", "provenance/source/source_manifest.json")
            for name, path in annotation_files.items():
                if path.resolve() not in original_destinations:
                    copy("annotation_" + name.replace(".", "_"), path, "provenance/annotation/" + str(path.relative_to(annotations)))
            copy("annotation_manifest", annotations / "annotation_manifest.json", "provenance/annotation/annotation_manifest.json")
            copy("scanner_code", evidence / "code_scanner.py", "provenance/evidence/code_scanner.py")
            for name, path in alignment_files.items():
                if path.resolve() not in original_destinations:
                    copy("alignment_" + name, path, str(path.relative_to(alignments)))
            copy("alignment_manifest", alignments / "alignment_manifest.json", "provenance/alignment/alignment_manifest.json")
            for name in ("freeze_unite_snapshot.py", "build_unite_regional_alignments.py", "build_unite_reference_panel.py"):
                copy("pipeline_" + name, Path(__file__).with_name(name), "provenance/pipeline/" + name)
            total = sum(item["bytes"] for item in artifacts.values())
            require(total <= max_bytes, "Snapshot exceeds the new-data byte budget")
            # Detect a completed upstream artifact changing between validation
            # and copying; all copied upstream proofs must remain those checked.
            require(all(sha(path) == inputs[name] for name, path in (
                ("source_manifest", panel / "source_manifest.json"), ("annotation_manifest", annotations / "annotation_manifest.json"),
                ("search_manifest", evidence / "site_search_manifest.json"), ("alignment_manifest", alignments / "alignment_manifest.json"))),
                "An upstream completion proof changed during freezing")
            require(artifacts["sequences"]["sha256"] == search["sequence_sha256"] and
                    artifacts["sites"]["sha256"] == search["site_artifact_sha256"] and
                    artifacts["catalog"]["sha256"] == search["catalog_artifact_sha256"] and
                    artifacts["pairs"]["sha256"] == search["pair_artifact_sha256"], "Copied site-search inputs changed")
            manifest = dict(schema_version=SCHEMA_VERSION, marker_id="UNITE_ITS", id=version_id, version_number=ordinal,
                label=label or f"UNITE ITS {ordinal} · {len(metadata):,} retained SH references", identity_sha256=identity,
                freeze_rule_version=FREEZE_RULE_VERSION, built_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                n_source=source["counts"]["source_reference_identities"], n_selected=len(metadata), n_catalog=len(catalog), n_pairs=len(pairs),
                n_sparse_sites=len(sites), n_alignments=len(index), n_alignment_records=alignment["n_alignment_records"],
                source_version=source["source_version"], source_doi=source["source_doi"], source_url=source["source_url"],
                source_input_sha256=source["input_sha256"], source_license=SOURCE_LICENSE, input_manifest_sha256=inputs,
                identity_basis=identity_basis,
                sequence_basis="Correctly oriented, otherwise unchanged developer-release bases; original source sequence hashes remain in metadata",
                taxonomy_source="Frozen UNITE release headers; hierarchy-specific unresolved and incertae groups remain explicit",
                limitations=["This pinned UNITE10.0 February2025 bulk release is distinct from a live UNITE service/version.",
                    "Selected references are balanced per hierarchy-scoped order, not a census sample, PCR performance estimate, or whole-taxon coverage.",
                    "Absent exact sites may reflect retained-fragment coverage or mismatches; only the exhaustive frozen search proof permits sparse zero-hit inference.",
                    "ITSx boundaries are computational; SSU/LSU intervals are detected retained flanks, not complete genes.",
                    "Regional MSAs compare bounded exemplars within each resolved class/order; full panel denominators are independent of display membership.",
                    "The original public source archive is identified by URL and checksum; normalized selected source bytes and the complete identity census are embedded."],
                artifacts=artifacts, artifact_bytes=total)
            atomic_json(staging / "manifest.json", manifest)
            validate_frozen(staging, manifest)
            target = versions / version_id
            require(not target.exists(), "Refusing to overwrite a historical snapshot")
            staging.rename(target)
            entry = dict(id=version_id, label=manifest["label"], manifest=f"versions/{version_id}/manifest.json",
                         sha256=sha(target / "manifest.json"), n_selected=len(metadata))
            registry["versions"].append(entry)
            atomic_json(registry_path, registry)
            atomic_json(destination / "latest.json", entry)
            return manifest
        finally:
            if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=Path("data/derived/unite_reference_panel"))
    parser.add_argument("--annotations", type=Path, default=Path("data/derived/unite_reference_panel/annotations_v2"))
    parser.add_argument("--evidence", type=Path, default=Path("data/derived/unite_reference_panel/evidence"))
    parser.add_argument("--alignments", type=Path, default=Path("data/derived/unite_reference_panel/regional_alignments"))
    parser.add_argument("--destination", type=Path, default=Path("data/pinned/UNITE_ITS"))
    parser.add_argument("--label")
    parser.add_argument("--max-bytes", type=int, default=70 * 1024 * 1024)
    args = parser.parse_args()
    result = freeze_snapshot(args.panel, args.annotations, args.evidence, args.alignments, args.destination, args.label, args.max_bytes)
    print(json.dumps({field: result[field] for field in ("id", "n_source", "n_selected", "n_catalog", "n_pairs", "n_alignments", "artifact_bytes")}, indent=2))
