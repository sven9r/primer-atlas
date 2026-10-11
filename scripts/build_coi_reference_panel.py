#!/usr/bin/env python3
"""Acquire an append-only, accession-versioned arthropod COI panel.

Original responses and legacy FASTA bytes are archived. New templates are
explicitly annotated COI features, not arbitrary complete mitochondrial records.
Query counts, retained records, marker QC and primer-site scoreability are
different quantities; this builder establishes only the first three.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
CLIENT = ROOT / "vendor/ncbi_entrez/ncbi_entrez.py"
SCRIPT_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
QUERY_SUFFIX = ("[Organism:exp] AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene])"
                " AND 500:20000[SLEN] NOT COII[Title]")
FIELDS = ("marker_id", "target_group", "accession", "entrez_uid", "taxid", "organism", "order", "family", "genus", "species",
          "taxonomy_status", "taxonomy_source_sha256", "taxonomy_source_file", "first_release_id", "added_at_utc", "header",
          "sequence_sha256", "header_sha256", "template_length_bp", "qc_status", "sequence_kind",
          "feature_start", "feature_end", "feature_strand", "feature_boundary_state", "marker_length_bp",
          "source_sha256", "source_file", "source_acquired_at_utc", "previous_run_manifest_sha256",
          "baseline_source", "original_raw_bytes_available", "published_release_id",
          "published_artifact_url", "published_artifact_sha256")


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def file_hash(path):
    return digest(Path(path).read_bytes())


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def read_csv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def accession_from_header(header):
    token = header.split()[0]
    token = re.sub(r"^_R_", "", token)
    if "|" in token:
        candidates = [part for part in token.split("|") if re.fullmatch(r"[A-Z][A-Z_]*[0-9]+\.[0-9]+", part)]
        if len(candidates) != 1:
            raise ValueError("Legacy header lacks one unique accession version")
        token = candidates[0]
    if not re.fullmatch(r"[A-Z][A-Z_]*[0-9]+\.[0-9]+", token):
        raise ValueError("Reference identity must contain an accession version: " + token)
    return token


def read_fasta(path):
    records, header, parts = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                accession = accession_from_header(header)
                if accession in records:
                    raise ValueError("Repeated FASTA accession version: " + accession)
                records[accession] = (header, "".join(parts))
            header, parts = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError("FASTA sequence lacks a header")
            parts.append(line.strip())
    if header is not None:
        accession = accession_from_header(header)
        if accession in records:
            raise ValueError("Repeated FASTA accession version: " + accession)
        records[accession] = (header, "".join(parts))
    return records


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for header, sequence in records.values():
            handle.write(">" + header + "\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def request(endpoint, params, output, offline=False, refresh=False):
    key = digest(json.dumps(dict(endpoint=endpoint, params=params), sort_keys=True))[:24]
    suffix = "json" if params.get("retmode") == "json" else "xml"
    destination = output / "raw" / "cache" / f"{endpoint}_{key}.{suffix}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists() or refresh:
        if offline:
            raise RuntimeError("Offline acquisition cache missing: " + destination.name)
        config = dict(endpoint=endpoint, params=params, max_items=4, max_depth=2,
                      timeout_sec=60, save_raw=True, raw_output_path=str(destination))
        last_error = "unknown"
        for attempt in range(3):
            time.sleep(0.4 if attempt == 0 else 2 ** attempt)
            try:
                response = subprocess.run([sys.executable, str(CLIENT)], input=json.dumps(config),
                                          capture_output=True, text=True, timeout=90, check=True)
                reply = json.loads(response.stdout)
                if reply.get("ok") and destination.exists():
                    break
                last_error = reply.get("error", {}).get("code", "missing_response")
            except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as error:
                last_error = type(error).__name__
        else:
            raise RuntimeError(f"NCBI {endpoint} failed after three attempts: {last_error}")
    # Immutable content addresses keep old source hashes valid across refreshes.
    source_hash = file_hash(destination)
    frozen = output / "raw" / "content" / f"{source_hash}.{suffix}"
    frozen.parent.mkdir(parents=True, exist_ok=True)
    if not frozen.exists():
        shutil.copy2(destination, frozen)
    return frozen


def qualifiers(feature):
    result = collections.defaultdict(list)
    for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
        result[qualifier.findtext("GBQualifier_name", "")].append(qualifier.findtext("GBQualifier_value", ""))
    return result


def is_coi_feature(feature):
    values = qualifiers(feature)
    genes = [re.sub(r"[^a-z0-9]", "", value.lower()) for value in values.get("gene", [])]
    products = values.get("product", [])
    return any(gene in {"coi", "co1", "cox1", "coxi", "mtco1", "mtcox1"} for gene in genes) or any(
        re.search(r"\bcytochrome\s+(?:c\s+)?oxidase\s+(?:subunit\s+)?(?:i|1)\b", product, re.I)
        for product in products)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN"))[::-1]


def extract_coi(record, minimum=500, maximum=2000, max_ambiguity=0.02):
    sequence = record.findtext("GBSeq_sequence", "").upper()
    if not sequence or re.search(r"[^ACGTURYSWKMBDHVN]", sequence):
        return None, "missing_or_invalid_record_sequence"
    candidates = [feature for feature in record.findall("./GBSeq_feature-table/GBFeature")
                  if feature.findtext("GBFeature_key") in {"CDS", "gene"} and is_coi_feature(feature)]
    cds = [feature for feature in candidates if feature.findtext("GBFeature_key") == "CDS"]
    candidates = cds or candidates
    if len(candidates) != 1:
        return None, "missing_coi_feature" if not candidates else "ambiguous_coi_features"
    feature = candidates[0]
    values = qualifiers(feature)
    if "pseudo" in values or "pseudogene" in values:
        return None, "annotated_pseudogene"
    intervals = feature.findall("./GBFeature_intervals/GBInterval")
    if len(intervals) != 1:
        return None, "compound_coi_feature_withheld"
    interval = intervals[0]
    left, right = interval.findtext("GBInterval_from"), interval.findtext("GBInterval_to")
    if not left or not right:
        return None, "missing_coi_feature_bounds"
    left, right = int(left), int(right)
    start, end = min(left, right), max(left, right)
    if start < 1 or end > len(sequence):
        return None, "coi_feature_outside_record"
    complement = interval.find("GBInterval_iscomp")
    strand = -1 if (complement is not None and complement.get("value", "true") != "false") or left > right else 1
    marker = sequence[start - 1:end]
    if strand == -1:
        marker = reverse_complement(marker)
    ambiguous = sum(base not in "ACGT" for base in marker)
    if len(marker) < minimum or len(marker) > maximum:
        return None, "coi_feature_length_outside_qc"
    if len(marker) - ambiguous < minimum or ambiguous / len(marker) > max_ambiguity:
        return None, "coi_feature_ambiguity_outside_qc"
    location = feature.findtext("GBFeature_location", "")
    # Fuzzy ends describe a partial marker, while its record coordinates still
    # identify real sequenced bases. Never interpret them as complete-gene ends.
    boundary = "partial" if "<" in location or ">" in location else "complete_annotated_interval"
    return dict(sequence=marker, feature_start=start, feature_end=end, feature_strand=strand,
                feature_boundary_state=boundary, marker_length_bp=len(marker)), "verified_marker_feature"


def record_identity(record):
    accession = record.findtext("GBSeq_accession-version", "")
    uid = next((value.text.split(":", 1)[1] for value in record.findall("./GBSeq_other-seqids/GBSeqid")
                if (value.text or "").startswith("gi:")), "")
    taxids = []
    for feature in record.findall("./GBSeq_feature-table/GBFeature"):
        if feature.findtext("GBFeature_key") == "source":
            taxids.extend(value.split(":", 1)[1] for value in qualifiers(feature).get("db_xref", []) if value.startswith("taxon:"))
    return accession, uid, taxids[0] if len(set(taxids)) == 1 else ""


def taxonomy_rows(document):
    result = {}
    for taxon in document.findall("./Taxon"):
        taxid = taxon.findtext("TaxId", "")
        nodes = list(taxon.findall("./LineageEx/Taxon")) + [taxon]
        row = dict(taxid=taxid, organism=taxon.findtext("ScientificName", ""), taxonomy_status="resolved")
        for rank in ("order", "family", "genus", "species"):
            values = [node.findtext("ScientificName", "") for node in nodes if node.findtext("Rank") == rank]
            row[rank] = values[0] if len(values) == 1 else ""
        row["lineage_names"] = [node.findtext("ScientificName", "") for node in nodes]
        result[taxid] = row
    return result


def plan_targets(groups, previous_counts, target_total, minimum, growth=0.01, monthly=False):
    if minimum < 1 or target_total < minimum * len(groups) or growth < 0:
        raise ValueError("Reference target cannot be below the per-group minimum")
    targets = {group: max(minimum, math.ceil(previous_counts.get(group, 0) * (1 + growth))
                         if monthly else previous_counts.get(group, 0)) for group in groups}
    total = max(target_total, sum(targets.values()))
    while sum(targets.values()) < total:
        group = min(groups, key=lambda value: (targets[value], groups.index(value)))
        targets[group] += 1
    return targets


def diverse_selection(rows, number):
    # Round-robin source organism labels, without claiming those labels are
    # accepted species identities. This reduces repeated-label dominance.
    bins = collections.defaultdict(list)
    for row in sorted(rows, key=lambda item: (-int(item.get("marker_length_bp") or 0), digest(item["accession"]))):
        bins[row.get("organism") or row["accession"]].append(row)
    queues = [collections.deque(bins[key]) for key in sorted(bins)]
    selected = []
    while queues and len(selected) < number:
        for queue in queues:
            if len(selected) == number:
                break
            selected.append(queue.popleft())
        queues = [queue for queue in queues if queue]
    return selected


def import_previous(directory, output):
    records = read_fasta(directory / "panel.fasta")
    metadata = read_csv(directory / "panel_manifest.csv")
    if len(metadata) != len(records) or {row["accession"] for row in metadata} != set(records):
        raise RuntimeError("Previous manifest and FASTA accession sets differ")
    previous_hash = file_hash(directory / "run_manifest.json")
    previous_manifest = output / "raw" / "lineage" / f"{previous_hash}.json"
    previous_manifest.parent.mkdir(parents=True, exist_ok=True)
    if not previous_manifest.exists():
        shutil.copy2(directory / "run_manifest.json", previous_manifest)
    for row in metadata:
        header, sequence = records[row["accession"]]
        if digest(sequence) != row["sequence_sha256"] or digest(header) != row["header_sha256"]:
            raise RuntimeError("Previous sequence or header hash differs: " + row["accession"])
        row["previous_run_manifest_sha256"] = previous_hash
    for folder in ("raw", "legacy_raw"):
        source = directory / folder
        if source.exists() and source.resolve() != (output / folder).resolve():
            shutil.copytree(source, output / folder, dirs_exist_ok=True)
    return records, {row["accession"]: row for row in metadata}, previous_hash


def verify_retained_marker_qc(records, metadata, output):
    """Recheck stored raw evidence without changing a retained sequence/header."""
    source_cache, taxon_cache = {}, {}
    for accession, row in metadata.items():
        source = row.get("source_file")
        if not source:
            row["qc_status"] = "retained_marker_source_unavailable"
            continue
        path = output / source
        expected_hash = row.get("source_sha256")
        if source not in source_cache:
            if not path.exists() or file_hash(path) != expected_hash:
                raise RuntimeError("Retained marker source checksum failed: " + accession)
            source_cache[source] = {node.findtext("GBSeq_accession-version"): node
                                    for node in ET.parse(path).getroot().findall("./GBSeq")}
        node = source_cache[source].get(accession)
        if node is None:
            raise RuntimeError("Retained marker source omits its accession version: " + accession)
        taxon_hash = row.get("taxonomy_source_sha256")
        if taxon_hash and taxon_hash not in taxon_cache:
            taxonomy_path = output / (row.get("taxonomy_source_file") or f"raw/content/{taxon_hash}.xml")
            if not taxonomy_path.exists() or file_hash(taxonomy_path) != taxon_hash:
                raise RuntimeError("Retained taxonomy source checksum failed: " + accession)
            taxon_cache[taxon_hash] = taxonomy_rows(ET.parse(taxonomy_path).getroot())
        taxon = taxon_cache.get(taxon_hash, {}).get(row.get("taxid"), {})
        feature, reason = extract_coi(node)
        if row["target_group"] not in taxon.get("lineage_names", []):
            row["qc_status"] = "retained_taxonomy_membership_unresolved"
        elif feature is None:
            row["qc_status"] = "retained_marker_unverified_" + reason
        else:
            expected = feature["sequence"] if row.get("sequence_kind") == "extracted_coi_feature" else node.findtext("GBSeq_sequence", "").upper()
            if records[accession][1].upper() != expected:
                row["qc_status"] = "retained_raw_differs_from_marker_source"
            else:
                row["qc_status"] = "verified_retained_marker_evidence"


def archive_legacy(directory, output, groups):
    result, archived = {}, []
    for path in sorted(directory.glob("*/*_GB.fasta")):
        group = path.parent.name
        destination = output / "legacy_raw" / group / path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and file_hash(destination) != file_hash(path):
            raise RuntimeError("Archived legacy FASTA changed: " + path.name)
        if not destination.exists():
            shutil.copy2(path, destination)
        archived.append(dict(target_group=group, path=str(destination.relative_to(output)), sha256=file_hash(path),
                             active_group=group in groups, records=len(read_fasta(path))))
        if group not in groups:
            continue
        for accession, (header, sequence) in read_fasta(path).items():
            if accession in result and result[accession][1:3] != (header, sequence):
                raise RuntimeError("Legacy accession has conflicting sequence/header/group evidence: " + accession)
            result[accession] = (group, header, sequence, file_hash(path))
    return result, archived


def run(args):
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    group_rows = [row for row in read_csv(args.groups_file) if row["order"] != "Bivalvia"]
    groups = [row["order"] for row in group_rows]
    built_at = utc_now()
    records, metadata, previous_hash = {}, {}, ""
    if args.previous_dir:
        records, metadata, previous_hash = import_previous(Path(args.previous_dir), output)
        verify_retained_marker_qc(records, metadata, output)
    legacy, legacy_files = archive_legacy(Path(args.legacy_raw_dir), output, groups) if args.legacy_raw_dir else ({}, [])
    if args.previous_dir:
        previous_manifest = json.loads((Path(args.previous_dir) / "run_manifest.json").read_text())
        inherited = {row["path"]: row for row in previous_manifest.get("legacy_archives", [])}
        inherited.update({row["path"]: row for row in legacy_files})
        legacy_files = list(inherited.values())
    published = read_csv(args.published_baseline) if args.published_baseline else []
    required = {row["accession"]: row for row in published if row["target_group"] in groups}
    for accession, (group, header, sequence, source_hash) in legacy.items():
        required.setdefault(accession, dict(target_group=group, accession=accession, template=header,
                                           baseline_source="local legacy raw FASTA"))
    counts = collections.Counter(row["target_group"] for row in metadata.values())
    for accession, row in required.items():
        if accession not in metadata:
            counts[row["target_group"]] += 1
    targets = plan_targets(groups, counts, args.target_total, args.minimum_per_group,
                           args.growth, monthly=bool(args.previous_dir or required))
    for accession, prior in required.items():
        if accession in metadata:
            if prior.get("published_release_id") and not metadata[accession].get("published_release_id"):
                metadata[accession]["first_release_id"] = prior["published_release_id"]
                if metadata[accession].get("sequence_kind") == "published_accession_record_recovery":
                    metadata[accession]["baseline_source"] = "published release accession recovery"
            for field in ("published_release_id", "published_artifact_url", "published_artifact_sha256"):
                if prior.get(field) and not metadata[accession].get(field):
                    metadata[accession][field] = prior[field]
    searches, rejections, taxonomy_cache = [], [], {}

    def fetch_taxonomy(taxids):
        missing = sorted(set(taxids) - set(taxonomy_cache) - {""})
        for start in range(0, len(missing), 100):
            ids = missing[start:start + 100]
            path = request("efetch", dict(db="taxonomy", id=",".join(ids), retmode="xml"),
                           output, args.offline)
            values = taxonomy_rows(ET.parse(path).getroot())
            for taxid, row in values.items():
                row["taxonomy_source_sha256"] = file_hash(path)
                row["taxonomy_source_file"] = str(path.relative_to(output))
                taxonomy_cache[taxid] = row
            for taxid in ids:
                taxonomy_cache.setdefault(taxid, dict(taxid=taxid, taxonomy_status="unresolved",
                                                     taxonomy_source_sha256=file_hash(path)))

    def acquire(ids, group, mandatory=False):
        candidates = []
        for start in range(0, len(ids), 100):
            selected = ids[start:start + 100]
            path = request("efetch", dict(db="nuccore", id=",".join(selected), rettype="gb", retmode="xml"),
                           output, args.offline)
            fetched = list(ET.parse(path).getroot().findall("./GBSeq"))
            response_hash = file_hash(path)
            if not fetched:
                raise RuntimeError("NCBI returned no sequence records for an acquisition batch")
            identities = [record_identity(record) for record in fetched]
            returned_accessions = [identity[0] for identity in identities]
            returned_uids = [identity[1] for identity in identities]
            requested = set(selected)
            returned = set(returned_accessions) if all("." in value for value in selected) else set(returned_uids)
            if len(fetched) != len(selected) or returned != requested or len(set(returned_accessions)) != len(fetched):
                raise RuntimeError("NCBI returned a different accession-version/UID set from the requested batch")
            fetch_taxonomy([identity[2] for identity in identities])
            for record, (accession, uid, taxid) in zip(fetched, identities):
                if accession in metadata:
                    continue
                feature, status = extract_coi(record)
                taxon = taxonomy_cache.get(taxid, {})
                membership = group in taxon.get("lineage_names", [])
                if not membership and not mandatory:
                    rejections.append(dict(target_group=group, accession=accession, reason="taxonomy_membership_unresolved", source_sha256=response_hash))
                    continue
                if feature is None and not mandatory:
                    rejections.append(dict(target_group=group, accession=accession, reason=status, source_sha256=response_hash))
                    continue
                legacy_value = legacy.get(accession)
                prior = required.get(accession, {})
                if legacy_value is not None:
                    header, sequence = legacy_value[1:3]
                    original_available = True
                    source_record = record.findtext("GBSeq_sequence", "").upper()
                    if sequence.upper() != source_record:
                        status = "legacy_raw_differs_from_accession_record"
                    elif feature is not None:
                        status = "verified_legacy_marker_record"
                    else:
                        status = "legacy_marker_unverified_" + status
                    kind = "legacy_raw_record"
                elif mandatory:
                    header = prior.get("template") or prior.get("header") or accession
                    sequence = record.findtext("GBSeq_sequence", "").upper()
                    original_available = False
                    status = "verified_published_marker_recovery" if feature is not None else "published_marker_unverified_" + status
                    kind = "published_accession_record_recovery"
                else:
                    sequence = feature["sequence"]
                    header = accession + " " + group + " annotated COI feature"
                    original_available, kind = True, "extracted_coi_feature"
                if not sequence:
                    raise RuntimeError("Retained accession has no recoverable sequence: " + accession)
                expected_hash = prior.get("sequence_sha256")
                if expected_hash and digest(sequence) != expected_hash:
                    raise RuntimeError("Recovered published sequence hash differs: " + accession)
                row = dict(marker_id="COI", target_group=group, accession=accession, entrez_uid=uid, taxid=taxid,
                           organism=taxon.get("organism") or record.findtext("GBSeq_organism", ""),
                           taxonomy_status=taxon.get("taxonomy_status", "unresolved"),
                           taxonomy_source_sha256=taxon.get("taxonomy_source_sha256", ""),
                           taxonomy_source_file=taxon.get("taxonomy_source_file", ""),
                           first_release_id=prior.get("first_release_id") or prior.get("published_release_id") or args.release_id,
                           added_at_utc=prior.get("added_at_utc") or built_at, header=header,
                           sequence_sha256=digest(sequence), header_sha256=digest(header), template_length_bp=len(sequence),
                           qc_status=status, sequence_kind=kind, source_sha256=response_hash,
                           source_file=str(path.relative_to(output)),
                           source_acquired_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(path.stat().st_mtime)),
                           previous_run_manifest_sha256=previous_hash,
                           baseline_source=prior.get("baseline_source") or ("published release accession recovery" if mandatory else "new NCBI annotated COI feature"),
                           original_raw_bytes_available=original_available)
                for rank in ("order", "family", "genus", "species"):
                    row[rank] = taxon.get(rank, "")
                if feature is not None:
                    row.update({key: value for key, value in feature.items() if key != "sequence"})
                if mandatory and not membership:
                    row["qc_status"] = "retained_taxonomy_membership_unresolved"
                for field in ("published_release_id", "published_artifact_url", "published_artifact_sha256"):
                    row[field] = prior.get(field, "")
                row["_sequence"] = sequence
                candidates.append(row)
            print(f"{group}: recovered/fetched {min(start + 100, len(ids))}/{len(ids)} records; accepted candidates={len(candidates)}", flush=True)
        return candidates

    # Recover every previously published identity before searching for additions.
    for group in groups:
        ids = sorted(accession for accession, row in required.items() if row["target_group"] == group and accession not in metadata)
        for row in acquire(ids, group, mandatory=True) if ids else []:
            sequence = row.pop("_sequence")
            records[row["accession"]] = (row["header"], sequence)
            metadata[row["accession"]] = row
    verified_counts = collections.Counter(row["target_group"] for row in metadata.values()
                                         if row["qc_status"].startswith("verified"))
    verified_targets = plan_targets(groups, verified_counts, args.target_total, args.minimum_per_group)
    for group in groups:
        retained_n = sum(row["target_group"] == group for row in metadata.values())
        targets[group] = max(targets[group], retained_n + max(0, verified_targets[group] - verified_counts[group]))
    qa = []
    for group in groups:
        previous_n = counts[group]
        retained = [row for row in metadata.values() if row["target_group"] == group]
        requested = max(0, targets[group] - len(retained))
        query = group + QUERY_SUFFIX
        candidate_limit = max(2000, targets[group] * 3)
        path = request("esearch", dict(db="nuccore", term=query, retmode="json", retmax=candidate_limit, idtype="acc"),
                       output, args.offline, args.refresh)
        search = json.loads(path.read_text())["esearchresult"]
        existing_accessions = set(metadata)
        uids = [uid for uid in search.get("idlist", []) if uid not in existing_accessions]
        search_row = dict(target_group=group, query=query, available_query_records=int(search["count"]),
                          fetched_candidate_limit=candidate_limit, source_sha256=file_hash(path),
                          source_file=str(path.relative_to(output)))
        searches.append(search_row)
        pool, consumed = [], 0
        while len(pool) < requested and consumed < len(uids):
            batch_size = max(100, math.ceil((requested - len(pool)) * 1.25))
            batch = uids[consumed:consumed + batch_size]
            pool.extend(acquire(batch, group))
            consumed += len(batch)
        chosen = diverse_selection(pool, requested)
        for row in chosen:
            sequence = row.pop("_sequence")
            records[row["accession"]] = (row["header"], sequence)
            metadata[row["accession"]] = row
        after = sum(row["target_group"] == group for row in metadata.values())
        verified = sum(row["target_group"] == group and row["qc_status"].startswith("verified") for row in metadata.values())
        qa.append(dict(target_group=group, previous_retained=previous_n, target_retained=targets[group],
                       added_new=after - previous_n, retained_after=after, verified_marker_records=verified,
                       minimum_per_group=args.minimum_per_group, target_met=after >= targets[group],
                       minimum_met=verified >= args.minimum_per_group, verified_target=verified_targets[group],
                       verified_target_met=verified >= verified_targets[group],
                       available_query_records=int(search["count"]),
                       candidate_records_examined=consumed,
                       target_shortfall_reason="" if after >= targets[group] else "bounded candidates exhausted before target; query counts are pre-QC",
                       release_id=args.release_id))
        print(f"{group}: retained={after}, target={targets[group]}, verified marker records={verified}", flush=True)
    rows = sorted(metadata.values(), key=lambda row: (groups.index(row["target_group"]), row["accession"]))
    write_fasta(output / "panel.fasta", {row["accession"]: records[row["accession"]] for row in rows})
    for name in ("panel_manifest.csv", "ledger.csv", "taxonomy.csv"):
        write_csv(output / name, rows, FIELDS)
    write_csv(output / "growth_qa.csv", qa, list(qa[0]))
    write_csv(output / "rejections.csv", rejections, ("target_group", "accession", "reason", "source_sha256"))
    artifacts = {name: dict(sha256=file_hash(output / name), path=name) for name in (
        "panel.fasta", "panel_manifest.csv", "ledger.csv", "taxonomy.csv", "growth_qa.csv", "rejections.csv")}
    manifest = dict(schema_version="1.0.0", marker_id="COI", release_id=args.release_id, built_at_utc=built_at,
                    target_total_minimum=args.target_total, minimum_per_group=args.minimum_per_group,
                    monthly_growth_fraction=args.growth, previous_run_manifest_sha256=previous_hash,
                    retained_count=len(rows), verified_marker_count=sum(row["qc_status"].startswith("verified") for row in rows),
                    published_baseline_count=len(published), archived_legacy_records=sum(row["records"] for row in legacy_files),
                    archived_legacy_active_records=sum(row["records"] for row in legacy_files if row.get("active_group")),
                    groups=groups, targets=targets, searches=searches, legacy_archives=legacy_files, artifacts=artifacts,
                    source_client_sha256=file_hash(CLIENT), acquisition_code_sha256=SCRIPT_SHA256,
                    limitations=["Query counts precede marker-feature and taxonomy QC.",
                                 "A retained COI fragment does not establish either primer binding site.",
                                 "Public accession recoveries cannot prove unavailable original raw-file byte identity.",
                                 "Sampling balances target groups and source organism labels within a bounded cached search; it is not random or exhaustive."])
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if any(not row["target_met"] or not row["verified_target_met"] for row in qa):
        raise RuntimeError("One or more groups did not meet the retained target; inspect growth_qa.csv")
    print(f"COI panel complete: {len(rows)} unique accession versions across {len(groups)} target groups", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "data/derived/coi_reference_panel"))
    parser.add_argument("--groups-file", default=str(ROOT / "data/orders_22.csv"))
    parser.add_argument("--previous-dir")
    parser.add_argument("--legacy-raw-dir")
    parser.add_argument("--published-baseline")
    parser.add_argument("--target-total", type=int, default=10000)
    parser.add_argument("--minimum-per-group", type=int, default=100)
    parser.add_argument("--growth", type=float, default=0.01)
    parser.add_argument("--release-id", default=time.strftime("%Y-%m-%d", time.gmtime()))
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
