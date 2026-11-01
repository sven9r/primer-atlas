#!/usr/bin/env python3
"""Project sourced LSU D domains onto retained, independently annotated LSU features.

This is a conservative homology projection, not direct fungal experimental
annotation. Unsupported boundaries stay in the audit, outside the landmark map.
Two LSU-only MAFFT strategies must agree on each interval, and both boundary
flanks must meet explicit identity, coverage and uniqueness checks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
FLANK = 64
MIN_FLANK = 32
MIN_PAIRED_COVERAGE = 0.75
MIN_PAIRED_IDENTITY = 0.70
MIN_ANCHOR_MARGIN = 0.08
METHODS = ("localpair", "genafpair")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta(path):
    records, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError("Duplicate FASTA identity")
            records[name] = ""
        elif name:
            records[name] += line.strip().upper()
    return records


def csv_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def source_features(records_path, sequences):
    features = []
    for record in ET.parse(records_path).getroot().findall("GBSeq"):
        accession = record.findtext("GBSeq_accession-version")
        if accession not in sequences:
            continue
        if record.findtext("GBSeq_sequence", "").upper() != sequences[accession]:
            raise ValueError("Source record differs from frozen sequence: " + accession)
        for number, feature in enumerate(record.findall("./GBSeq_feature-table/GBFeature")):
            qualifiers = {q.findtext("GBQualifier_name"): q.findtext("GBQualifier_value", "")
                          for q in feature.findall("./GBFeature_quals/GBQualifier")}
            product = qualifiers.get("product", "")
            location = feature.findtext("GBFeature_location", "")
            match = re.fullmatch(r"([<>]?)(\d+)\.\.([<>]?)(\d+)", location)
            if feature.findtext("GBFeature_key") != "rRNA" or not match:
                continue
            gene = "28S" if re.search(r"(?:25|26|28)S", product) else "18S" if "18S" in product else None
            if gene is None:
                continue
            start, end = int(match[2]), int(match[4])
            if not 1 <= start <= end <= len(sequences[accession]):
                raise ValueError("Source feature outside frozen sequence")
            features.append(dict(accession=accession, gene=gene, start=start, end=end,
                                 partial=bool(match[1] or match[3]), source_location=location,
                                 source_feature_number=number, product=product))
    return features


def alignment_projection(reference, target):
    if len(reference) != len(target):
        raise ValueError("Alignment width differs")
    reference_position = target_position = 0
    columns, before, after = {}, [], []
    for column, (a, b) in enumerate(zip(reference, target)):
        before.append(target_position)
        reference_position += a != "-"
        target_position += b != "-"
        after.append(target_position)
        if a != "-":
            columns[reference_position] = column
    return columns, before, after


def competitive_anchor(anchor, target_raw, own_start, own_end):
    """Find a competitive, gap-aware full-query match outside the own locus.

    Query bases must all be consumed; target prefixes/suffixes are free. Each
    end state retains its best alignment (match +2, mismatch -1, gap -2), with
    match count and paired-base count for the same coverage/identity gates used
    at the projected locus. This is a sequence-specificity check, not evidence
    that a competing sequence has the same ribosomal structural annotation.

    Coordinates are zero-based half-open internally and one-based inclusive
    in the returned proof. The own flank plus one query length on either side
    is excluded so shifted alignments of the same locus are not competitors.
    """
    length = len(anchor)
    if not 1 <= length <= FLANK:
        raise ValueError("Competitive anchor exceeds the bounded flank length")
    canonical = set("ACGT")
    # State: alignment score, matches, canonical pairs, target start (0-based).
    previous = [(0, 0, 0, position) for position in range(len(target_raw) + 1)]
    for query_position, query_base in enumerate(anchor, 1):
        current = [(-2 * query_position, 0, 0, 0)]
        for target_position, target_base in enumerate(target_raw, 1):
            paired = int(query_base in canonical and target_base in canonical)
            matched = int(bool(paired) and query_base == target_base)
            diagonal = previous[target_position - 1]
            deletion = previous[target_position]
            insertion = current[target_position - 1]
            candidates = (
                (diagonal[0] + (2 if matched else -1), diagonal[1] + matched,
                 diagonal[2] + paired, diagonal[3]),
                (deletion[0] - 2, deletion[1], deletion[2], deletion[3]),
                (insertion[0] - 2, insertion[1], insertion[2], insertion[3]),
            )
            # Ties prefer more evidence, then the shorter target interval.
            current.append(max(candidates))
        previous = current
    excluded_start = max(0, own_start - length)
    excluded_end = min(len(target_raw), own_end + length)
    best, best_key = None, None
    for end, (score, matches, paired, start) in enumerate(previous):
        if end <= start or not (end <= excluded_start or start >= excluded_end):
            continue
        coverage = paired / length
        identity = matches / paired if paired else 0
        if coverage < MIN_PAIRED_COVERAGE or identity < MIN_PAIRED_IDENTITY:
            continue
        key = (score / (2 * length), matches / length, identity, coverage)
        if best_key is None or key > best_key:
            best_key = key
            best = dict(score=score, matches=matches, paired_bases=paired,
                        query_bases=length, target_start=start + 1, target_end=end,
                        paired_coverage=round(coverage, 6),
                        paired_identity=round(identity, 6),
                        normalized_alignment_score=round(score / (2 * length), 6),
                        normalized_match_score=round(matches / length, 6))
    return dict(method="full_query_semiglobal_DP", scoring="match2_mismatch-1_gap-2",
                score_normalization="alignment_score_divided_by_2_times_query_bases",
                coordinate_system="one_based_inclusive_target_feature_bases",
                excluded_target_start=excluded_start + 1,
                excluded_target_end=excluded_end,
                best_alternative=best)


def flank_proof(reference, target, positions, target_raw):
    columns, before, after = alignment_projection(reference, target)
    indices = [columns[p] for p in positions]
    paired = [(reference[i], target[i]) for i in indices
              if reference[i] in "ACGT" and target[i] in "ACGT"]
    identity = sum(a == b for a, b in paired) / len(paired) if paired else 0
    coverage = len(paired) / len(indices)
    lo, hi = before[indices[0]], after[indices[-1]]
    anchor = "".join(reference[i] for i in indices)
    competition = competitive_anchor(anchor, target_raw, lo, hi)
    best_alternative = competition["best_alternative"]
    alternative = (best_alternative["matches"] / len(anchor)) if best_alternative else 0
    # Compare like-for-like alignment scores. Match fractions alone would
    # falsely elevate distant, heavily gapped decoys to the own flank's rank.
    observed_score = 0
    for column in range(indices[0], indices[-1] + 1):
        a, b = reference[column], target[column]
        if a == b == "-":
            continue
        if a == "-" or b == "-":
            observed_score -= 2
        elif a == b and a in "ACGT":
            observed_score += 2
        else:
            observed_score -= 1
    observed = observed_score / (2 * len(anchor))
    alternative_score = best_alternative["score"] / (2 * len(anchor)) if best_alternative else 0
    margin = observed - alternative_score
    return dict(source_flank_bases=len(indices), paired_coverage=round(coverage, 6), paired_identity=round(identity, 6),
                best_alternative_identity=round(alternative, 6),
                best_alternative_match_fraction=round(alternative, 6),
                observed_alignment_score=observed_score,
                observed_normalized_alignment_score=round(observed, 6),
                best_alternative_normalized_alignment_score=round(alternative_score, 6),
                anchor_margin=round(margin, 6),
                competitive_anchor_proof=competition,
                coordinate_system="one_based_inclusive_target_feature_bases",
                target_flank_start=lo + 1, target_flank_end=hi,
                supported=(coverage >= MIN_PAIRED_COVERAGE and identity >= MIN_PAIRED_IDENTITY
                           and margin >= MIN_ANCHOR_MARGIN
                           and hi > lo and set(target_raw[lo:hi]) <= set("ACGT")))


def project_domain(aligned_reference, aligned_target, target_raw, start, end,
                   left_bases=FLANK, right_bases=FLANK):
    columns, before, after = alignment_projection(aligned_reference, aligned_target)
    if min(left_bases, right_bases) < MIN_FLANK or start <= left_bases or end + right_bases > len(columns):
        return dict(status="source_flank_unavailable")
    left = flank_proof(aligned_reference, aligned_target, range(start - left_bases, start), target_raw)
    right = flank_proof(aligned_reference, aligned_target, range(end + 1, end + right_bases + 1), target_raw)
    lo, hi = before[columns[start]] + 1, after[columns[end]]
    status = "supported" if left["supported"] and right["supported"] else "boundary_flanks_unresolved"
    if hi < lo:
        status = "no_target_bases_in_domain"
    elif set(target_raw[lo - 1:hi]) - set("ACGT"):
        status = "ambiguous_target_domain"
    return dict(status=status, start=lo, end=hi, left_flank=left, right_flank=right,
                coordinate_system="one_based_inclusive_target_feature_bases",
                boundary_method="alignment_column_edges_including_target_indels")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--domains", type=Path, action="append", required=True)
    parser.add_argument("--source-proof", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    alignments = args.output / "alignments"
    alignments.mkdir(exist_ok=True)
    sequences = fasta(args.snapshot / "sequences.fasta")
    manifest = json.loads((args.snapshot / "manifest.json").read_text())
    if sha(args.snapshot / "sequences.fasta") != manifest["artifacts"]["sequences"]["sha256"]:
        raise ValueError("Frozen input sequence checksum differs")
    if sha(args.snapshot / "landmarks.csv") != manifest["artifacts"]["landmarks"]["sha256"]:
        raise ValueError("Frozen input landmark checksum differs")
    source_run = json.loads((args.snapshot / "run_manifest.json").read_text())
    if sha(args.records) not in source_run.get("raw_responses", {}).values():
        raise ValueError("NCBI feature response is not the frozen acquisition source")
    base_landmarks = csv_rows(args.snapshot / "landmarks.csv")
    features = source_features(args.records, sequences)
    domains = [row for path in args.domains for row in csv_rows(path)]
    version_result = subprocess.run(["mafft", "--version"], capture_output=True, text=True, check=True)
    mafft_version = version_result.stderr.strip() or version_result.stdout.strip()
    if not mafft_version:
        raise ValueError("MAFFT version could not be recorded")
    sources = {}
    for path in args.domains:
        for accession in {row["source_accession"] for row in csv_rows(path)}:
            source_path = path.parent / (accession + ".fasta")
            source = fasta(source_path)
            if accession not in source:
                raise ValueError("Missing named structural reference")
            expected = {row["source_sequence_sha256"] for row in csv_rows(path)
                        if row["source_accession"] == accession}
            if expected != {hashlib.sha256(source[accession].encode()).hexdigest()}:
                raise ValueError("Structural reference sequence differs from its domain evidence")
            sources[accession] = source[accession]
    extra_landmarks, audit = [], []
    for feature in features:
        existing = [r for r in base_landmarks if r.get("accession") == feature["accession"]
                    and r["label"] == feature["gene"] and int(r["start"]) == feature["start"]
                    and int(r["end"]) == feature["end"]]
        feature_id = existing[0]["landmark_id"] if existing else (
            feature["gene"] + "_FEATURE_" + str(feature["source_feature_number"]))
        if not existing:
            extra_landmarks.append(dict(accession=feature["accession"], landmark_id=feature_id,
                                        label=feature["gene"], start=feature["start"], end=feature["end"],
                                        source_note="NCBI " + ("partial" if feature["partial"] else "complete")
                                                    + " forward feature " + feature["source_location"],
                                        landmark_type="gene", evidence_status=("partial_source_annotation" if feature["partial"]
                                                                               else "annotated_reference"),
                                        boundary_status="partial" if feature["partial"] else "complete",
                                        source_reference=feature["accession"],
                                        source_location=feature["source_location"]))
        if feature["gene"] != "28S":
            continue
        target = sequences[feature["accession"]][feature["start"] - 1:feature["end"]]
        aligned = {}
        for accession, reference in sources.items():
            for method in METHODS:
                identity = hashlib.sha256((reference + "|" + target + "|" + method + "|" + mafft_version).encode()).hexdigest()[:16]
                input_path = alignments / (identity + ".input.fasta")
                output_path = alignments / (identity + ".aligned.fasta")
                input_path.write_text(">STRUCTURAL_REFERENCE\n" + reference + "\n>TARGET_FEATURE\n" + target + "\n")
                if not output_path.exists():
                    result = subprocess.run(["mafft", "--" + method, "--maxiterate", "1000",
                                             "--thread", "1", "--quiet", str(input_path)],
                                            check=True, capture_output=True, text=True, timeout=120)
                    output_path.write_text(result.stdout)
                pair = fasta(output_path)
                if pair["STRUCTURAL_REFERENCE"].replace("-", "") != reference or pair["TARGET_FEATURE"].replace("-", "") != target:
                    raise ValueError("Alignment changed an input sequence")
                aligned[accession, method] = pair
        for domain in domains:
            # Unselected source domains still bound conserved flank windows.
            # This permits a human D10-only projection while retaining the
            # independently published human D9/D12 boundaries as context.
            if domain.get("projection_selected", "true").lower() != "true":
                continue
            source = domain["source_accession"]
            proof = dict(accession=feature["accession"], feature_id=feature_id,
                         feature_start=feature["start"], feature_end=feature["end"],
                         domain=domain["domain_id"], source_reference=source,
                         source_start=int(domain["source_start"]), source_end=int(domain["source_end"]))
            if domain["eligible_for_projection"].lower() != "true":
                proof["status"] = domain["source_status"]
                audit.append(proof)
                continue
            projections = {}
            # Longer boundary windows add evidence without lowering identity
            # thresholds. Stop each window before an adjacent sourced domain.
            adjacent = [d for d in domains if d["source_accession"] == source
                        and d["eligible_for_projection"].lower() == "true"]
            prior = [int(d["source_end"]) for d in adjacent if int(d["source_end"]) < proof["source_start"]]
            following = [int(d["source_start"]) for d in adjacent if int(d["source_start"]) > proof["source_end"]]
            left_bases = min(FLANK, proof["source_start"] - max(prior, default=0) - 1)
            right_bases = min(FLANK, min(following, default=len(sources[source]) + 1) - proof["source_end"] - 1)
            for method in METHODS:
                pair = aligned[source, method]
                projections[method] = project_domain(pair["STRUCTURAL_REFERENCE"], pair["TARGET_FEATURE"],
                                                     target, proof["source_start"], proof["source_end"],
                                                     left_bases=left_bases, right_bases=right_bases)
            proof["methods"] = projections
            accepted = all(p["status"] == "supported" for p in projections.values())
            intervals = {(p.get("start"), p.get("end")) for p in projections.values()}
            if len(intervals) != 1:
                accepted = False
                proof["status"] = "alignment_boundary_disagreement"
            else:
                proof["status"] = "mapped_reference_domain" if accepted else "boundary_mapping_unresolved"
            if accepted:
                lo, hi = next(iter(intervals))
                start, end = feature["start"] + lo - 1, feature["start"] + hi - 1
                proof.update(start=start, end=end, encoded_bp=end - start + 1)
                extra_landmarks.append(dict(accession=feature["accession"],
                    landmark_id=feature_id + "_" + domain["domain_id"], label=domain["domain_id"],
                    start=start, end=end, source_note="Structural reference domain projected from " + source,
                    landmark_type="divergent_domain", parent_landmark=feature_id,
                    evidence_status="mapped_reference_domain", boundary_status="projected",
                    domain_scheme=domain["source_scheme"], source_reference=source,
                    source_start=proof["source_start"], source_end=proof["source_end"],
                    source_url=domain["source_url"], source_citation=domain["source_citation"],
                    source_feature_start=feature["start"], source_feature_end=feature["end"],
                    mapping_method="LSU-only MAFFT localpair/genafpair consensus; supported unique boundary flanks"))
            audit.append(proof)
    landmarks = base_landmarks + extra_landmarks
    if len({(r.get("accession"), r["landmark_id"]) for r in landmarks}) != len(landmarks):
        raise ValueError("Repeated landmark identity")
    write_csv(args.output / "landmarks.csv", landmarks)
    (args.output / "domain_mapping_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    inputs = args.output / "sources"
    inputs.mkdir(exist_ok=True)
    for path in [args.records, *args.domains, *args.source_proof, *(p.parent / (a + ".fasta") for p in args.domains
                 for a in {r["source_accession"] for r in csv_rows(p)})]:
        destination = inputs / path.name
        shutil.copy2(path, destination)
    shutil.copy2(Path(__file__), args.output / "code_domain_mapping.py")
    counts = {accession: sum(r.get("landmark_type") == "divergent_domain" and r["accession"] == accession
                            for r in extra_landmarks) for accession in sequences}
    run = dict(schema_version="1.0.0", base_snapshot=manifest["id"],
               source="Named structural LSU domain references projected onto frozen NCBI LSU features",
               coordinate_system="one_based_inclusive_raw_template_bases",
               evidence_status="mapped_reference_domain", source_is_direct_fungal_annotation=False,
               methods=list(METHODS), mafft_version=mafft_version, max_flank_bases=FLANK,
               min_flank_bases=MIN_FLANK, min_paired_coverage=MIN_PAIRED_COVERAGE,
               min_paired_identity=MIN_PAIRED_IDENTITY, min_anchor_margin=MIN_ANCHOR_MARGIN,
               anchor_specificity_method="full_query_semiglobal_DP",
               anchor_specificity_scoring="match2_mismatch-1_gap-2",
               anchor_specificity_normalization="alignment_score_divided_by_2_times_query_bases",
               mapped_domains_per_accession=counts, mapped_domain_rows=sum(counts.values()),
               partial_source_gene_rows=sum(r.get("evidence_status") == "partial_source_annotation"
                                            for r in extra_landmarks),
               unresolved_domain_rows=sum(r["status"] != "mapped_reference_domain" for r in audit),
               limitations=["Reference-projected boundaries are computational homology inferences.",
                            "Deleted boundary bases use alignment-edge projections; the source coordinates are not transplanted as raw fungal positions.",
                            "Only independently bounded NCBI LSU features are projected.",
                            "Unresolved/absent domains are withheld; absence of a label does not prove biological absence."])
    run["artifacts_sha256"] = {str(p.relative_to(args.output)): sha(p) for p in sorted(args.output.rglob("*"))
                                if p.is_file() and p.name != "run_manifest.json"}
    (args.output / "run_manifest.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps({key: run[key] for key in ("mapped_domain_rows", "mapped_domains_per_accession",
                                               "partial_source_gene_rows", "unresolved_domain_rows")}, indent=2))


if __name__ == "__main__":
    main()
