#!/usr/bin/env python3
"""Build a bounded, accession-versioned NCBI fungal rDNA evidence panel.

Cached original Entrez responses make an offline rebuild possible. Short or
unannotated records are retained with explicit missing-boundary metadata rather
than force-fitted to the 739 bp pilot. This is a selected panel, not a survey of
fungal diversity or a claim of assay validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
CLIENT = ROOT / "vendor/ncbi_entrez/ncbi_entrez.py"
BASELINE = ROOT / "data/reference/its_fungal_reference_FN812768.2.fasta"
STRATA = ("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota", "Glomeromycota")
QUERY_TEMPLATE = ('{clade}[Organism:exp] AND ("18S ribosomal RNA"[Title] OR '
                  '"small subunit ribosomal RNA"[Title]) AND '
                  '"internal transcribed spacer"[Title] AND 2500:12000[SLEN]')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fasta(path):
    records, header, bases = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(bases).upper().replace("-", "")
            header, bases = line[1:].split()[0], []
        else:
            bases.append(line.strip())
    if header is not None:
        records[header] = "".join(bases).upper().replace("-", "")
    return records


def request(endpoint, params, destination, offline=False, refresh=False):
    if destination.exists() and not refresh:
        return destination
    if offline:
        raise RuntimeError(f"Offline acquisition cache missing: {destination.name}")
    config = {"endpoint": endpoint, "params": params, "max_items": 4,
              "max_depth": 2, "timeout_sec": 60, "save_raw": True,
              "raw_output_path": str(destination)}
    last_error = "unknown"
    for attempt in range(3):
        time.sleep(0.4 if attempt == 0 else 2 ** attempt)
        try:
            result = subprocess.run([sys.executable, str(CLIENT)], input=json.dumps(config),
                                    capture_output=True, text=True, check=True, timeout=90)
            reply = json.loads(result.stdout)
            if reply.get("ok") and destination.exists():
                return destination
            last_error = reply.get("error", {}).get("code", "missing_raw_response")
        except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
            last_error = type(exc).__name__
    raise RuntimeError(f"NCBI {endpoint} failed after three attempts: {last_error}")


def acquisition_time(path):
    """Original cache file time is an acquisition-time proxy, not rebuild time."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(Path(path).stat().st_mtime))


def search_cache_path(cache, clade, query, per_clade, provenance=None):
    """Cache identities include the exact query and requested record limit."""
    identity = hashlib.sha256(json.dumps(dict(query=query, per_clade=per_clade),
                                        sort_keys=True).encode()).hexdigest()[:16]
    destination = cache / f"search_{clade}_{identity}.json"
    legacy = cache / f"search_{clade}.json"
    if destination.exists() or not legacy.exists() or provenance is None:
        return destination
    # Migration requires the preceding run's exact query, hash, selected IDs,
    # and the raw response's retmax. A four-record cache never satisfies a new
    # request for twenty, even when that request uses the same named clade.
    for search in provenance.get("searches", []):
        if search.get("clade") != clade or search.get("query") != query:
            continue
        raw = json.loads(legacy.read_text()).get("esearchresult", {})
        count = int(raw.get("count", -1))
        expected = min(per_clade, count)
        if (search.get("sha256") == sha(legacy)
                and int(raw.get("retmax", -1)) in {per_clade, expected}
                and len(raw.get("idlist", [])) == expected
                and raw.get("idlist") == search.get("selected")):
            shutil.copy2(legacy, destination)
            return destination
    return destination


def validate_accessions(records, requested):
    returned = [record.findtext("GBSeq_accession-version") for record in records]
    if len(returned) != len(requested) or len(set(returned)) != len(returned) or set(returned) != set(requested):
        raise RuntimeError("NCBI returned a different accession-version set from the selected records")


def qualifiers(feature):
    return {q.findtext("GBQualifier_name", ""): q.findtext("GBQualifier_value", "")
            for q in feature.findall("./GBFeature_quals/GBQualifier")}


def feature_bounds(feature):
    location = feature.findtext("GBFeature_location", "")
    if "<" in location or ">" in location:
        return None
    intervals = feature.findall("./GBFeature_intervals/GBInterval")
    if len(intervals) != 1:
        return None
    interval = intervals[0]
    start, end = interval.findtext("GBInterval_from"), interval.findtext("GBInterval_to")
    if not start or not end or interval.find("GBInterval_iscomp") is not None:
        return None
    return min(int(start), int(end)), max(int(start), int(end))


def inspect_annotations(record, accession):
    known, occurrences, notes = {}, {}, []
    for feature in record.findall("./GBSeq_feature-table/GBFeature"):
        q = qualifiers(feature)
        # A single 'contains SSU, ITS1, 5.8S, ITS2, LSU' misc_RNA feature
        # describes the entire operon and cannot annotate individual boundaries.
        key = feature.findtext("GBFeature_key", "")
        if key not in ("rRNA", "misc_RNA", "misc_feature"):
            continue
        text = (q.get("product") or q.get("gene") or q.get("note", "")).lower()
        groups = [bool("5.8s" in text),
                  bool("18s" in text or "small subunit ribosomal" in text),
                  bool("28s" in text or "large subunit ribosomal" in text),
                  bool("internal transcribed spacer 1" in text or "its1" in text),
                  bool("internal transcribed spacer 2" in text or "its2" in text)]
        if sum(groups) != 1:
            if sum(groups) > 1:
                notes.append("Multi-region feature lacks separate rDNA boundaries.")
            continue
        label = None
        if "5.8s" in text:
            label = "5_8S"
        elif "18s" in text or "small subunit ribosomal" in text:
            label = "18S"
        elif "28s" in text or "large subunit ribosomal" in text:
            label = "28S"
        elif "internal transcribed spacer 1" in text or "its1" == text.strip():
            label = "ITS1"
        elif "internal transcribed spacer 2" in text or "its2" == text.strip():
            label = "ITS2"
        if label:
            occurrences[label] = occurrences.get(label, 0) + 1
            bounds = feature_bounds(feature)
            if bounds is None:
                notes.append(f"{label} feature withheld: fuzzy, compound, or reverse-strand boundary.")
            else:
                known[label] = (*bounds, "NCBI annotated complete forward feature")
    for label, count in occurrences.items():
        if count > 1:
            known.pop(label, None)
            notes.append(f"{label} feature withheld: {count} repeated gene labels are ambiguous.")
    for name, left, right in (("ITS1", "18S", "5_8S"), ("ITS2", "5_8S", "28S")):
        if name not in known and occurrences.get(name, 0) == 0 and left in known and right in known:
            start, end = known[left][1] + 1, known[right][0] - 1
            if start <= end:
                known[name] = (start, end, "between NCBI annotated complete forward rRNA features")
    rows = [dict(accession=accession, landmark_id=name, label=name.replace("_", "."),
                 start=value[0], end=value[1], source_note=value[2])
            for name, value in known.items()]
    if not rows and not notes:
        notes.append("No separate complete forward rDNA feature boundaries are annotated.")
    return rows, " ".join(dict.fromkeys(notes))


def annotate(record, accession):
    return inspect_annotations(record, accession)[0]


def retain_previous_panel(previous_dir, sequences, metadata, landmarks):
    previous = read_fasta(previous_dir / "panel.fasta")
    with (previous_dir / "panel_manifest.csv").open() as handle:
        previous_meta = list(csv.DictReader(handle))
    with (previous_dir / "landmarks.csv").open() as handle:
        previous_landmarks = list(csv.DictReader(handle))
    manifest_path = previous_dir / "run_manifest.json"
    previous_manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    rows_by_accession = {row["accession"]: row for row in previous_meta}
    if set(rows_by_accession) != set(previous):
        raise RuntimeError("Previous panel sequence and metadata accession sets differ")
    for accession, sequence in previous.items():
        if accession in sequences and sequences[accession] != sequence:
            raise RuntimeError(f"Changed sequence for retained accession-version: {accession}")
        row = rows_by_accession[accession]
        if row.get("sequence_sha256") and row["sequence_sha256"] != hashlib.sha256(sequence.encode()).hexdigest():
            raise RuntimeError(f"Previous panel metadata sequence hash differs: {accession}")
        # Preserve the preceding taxonomy assignment even for accessions found
        # again today. Its source can differ from newly added templates.
        row.setdefault("taxonomy_source_sha256", previous_manifest.get("taxonomy_snapshot_sha256", ""))
        sequences[accession] = sequence
    retained = set(previous)
    metadata[:] = [row for row in metadata if row["accession"] not in retained] + previous_meta
    landmarks[:] = [row for row in landmarks if row["accession"] not in retained] + previous_landmarks
    return list(previous), sha(manifest_path) if manifest_path.exists() else None


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/derived/fungal_reference_panel")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "data/external/fungal_reference_panel")
    parser.add_argument("--per-clade", type=int, default=4)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--previous-dir", type=Path, help="Retain a preceding panel and its frozen metadata")
    args = parser.parse_args()
    if args.per_clade < 1:
        parser.error("--per-clade must be positive")
    out, cache = args.output_dir, args.cache_dir
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    preceding_manifest = out / "run_manifest.json"
    if not preceding_manifest.exists():
        preceding_manifest = ROOT / "data/derived/fungal_reference_panel/run_manifest.json"
    cache_provenance = json.loads(preceding_manifest.read_text()) if preceding_manifest.exists() else None
    accessions, searches = ["FN812768.2"], []
    for clade in STRATA:
        query = QUERY_TEMPLATE.format(clade=clade)
        search_path = search_cache_path(cache, clade, query, args.per_clade, cache_provenance)
        path = request("esearch", dict(db="nuccore", term=query, retmode="json",
                       retmax=args.per_clade, idtype="acc"), search_path,
                       args.offline, args.refresh)
        search = json.loads(path.read_text())["esearchresult"]
        expected_count = min(args.per_clade, int(search["count"]))
        if int(search.get("retmax", -1)) not in {args.per_clade, expected_count} or len(search["idlist"]) != expected_count:
            raise RuntimeError(f"Search cache does not match requested record limit: {clade}")
        accessions.extend(search["idlist"])
        searches.append(dict(clade=clade, query=query, count_available=int(search["count"]),
                             selected=search["idlist"], sha256=sha(path),
                             cache_file=path.name, acquired_at=acquisition_time(path)))
    accessions = list(dict.fromkeys(accessions))
    identity = hashlib.sha256("\n".join(accessions).encode()).hexdigest()[:16]
    raw = request("efetch", dict(db="nuccore", id=",".join(accessions), rettype="gb", retmode="xml"),
                  cache / f"records_{identity}.xml", args.offline, args.refresh)
    records = ET.parse(raw).getroot().findall("GBSeq")
    validate_accessions(records, accessions)
    taxids = []
    for record in records:
        for feature in record.findall("./GBSeq_feature-table/GBFeature"):
            if feature.findtext("GBFeature_key") == "source":
                for q in feature.findall("./GBFeature_quals/GBQualifier"):
                    value = q.findtext("GBQualifier_value", "")
                    if q.findtext("GBQualifier_name") == "db_xref" and value.startswith("taxon:"):
                        taxids.append(value.split(":", 1)[1])
    taxids = list(dict.fromkeys(taxids))
    tax_identity = hashlib.sha256("\n".join(taxids).encode()).hexdigest()[:16]
    raw_tax = request("efetch", dict(db="taxonomy", id=",".join(taxids), retmode="xml"),
                      cache / f"taxonomy_{tax_identity}.xml", args.offline, args.refresh)
    taxa = {}
    for taxon in ET.parse(raw_tax).getroot().findall("Taxon"):
        ranks = {x.findtext("Rank"): x.findtext("ScientificName")
                 for x in taxon.findall("./LineageEx/Taxon")}
        ranks[taxon.findtext("Rank")] = taxon.findtext("ScientificName")
        taxa[taxon.findtext("TaxId")] = ranks
    sequences, metadata, landmarks = {}, [], []
    taxonomy_source_sha = sha(raw_tax)
    for record in records:
        accession = record.findtext("GBSeq_accession-version")
        seq = record.findtext("GBSeq_sequence", "").upper()
        taxid = None
        for feature in record.findall("./GBSeq_feature-table/GBFeature"):
            q = qualifiers(feature)
            if feature.findtext("GBFeature_key") == "source":
                for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
                    value = qualifier.findtext("GBQualifier_value", "")
                    if value.startswith("taxon:"):
                        taxid = value.split(":", 1)[1]
        ranks = taxa.get(taxid, {})
        if "Fungi" not in record.findtext("GBSeq_taxonomy", ""):
            raise RuntimeError(f"Nonfungal record returned: {accession}")
        boundaries, annotation_note = inspect_annotations(record, accession)
        sequences[accession] = seq
        metadata.append(dict(accession=accession, organism=record.findtext("GBSeq_organism"),
            taxid=taxid, phylum=ranks.get("phylum", "unresolved"),
            **{"class": ranks.get("class", "unresolved")},
            order=ranks.get("order", "unresolved"), family=ranks.get("family", "unresolved"),
            genus=ranks.get("genus", "unresolved"), length_bp=len(seq),
            sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(), source="NCBI GenBank",
            source_url=f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
            definition=record.findtext("GBSeq_definition"),
            annotation_status="annotated_rDNA" if {x["landmark_id"] for x in boundaries} == {"18S", "ITS1", "5_8S", "ITS2", "28S"} and not annotation_note else "partial_or_unresolved_boundaries",
            annotation_note=annotation_note,
            taxonomy_source_sha256=taxonomy_source_sha,
            reference_role="original_coordinate_reference" if accession == "FN812768.2" else "new_NCBI_template"))
        landmarks.extend(boundaries)
    baseline_seq = next(iter(read_fasta(BASELINE).values()))
    if sequences.get("FN812768.2") != baseline_seq:
        raise RuntimeError("Pinned FN812768.2 sequence differs from NCBI accession-version")
    previous_accessions = []
    previous_run_manifest_sha = None
    if args.previous_dir:
        previous_accessions, previous_run_manifest_sha = retain_previous_panel(
            args.previous_dir, sequences, metadata, landmarks)
    with (out / "panel.fasta").open("w") as handle:
        for accession, sequence in sequences.items():
            handle.write(f">{accession}\n")
            handle.write("\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + "\n")
    fields = ["accession", "organism", "taxid", "phylum", "class", "order", "family", "genus",
              "length_bp", "sequence_sha256", "source", "source_url", "definition", "annotation_status",
              "annotation_note", "taxonomy_source_sha256", "reference_role"]
    write_csv(out / "panel_manifest.csv", metadata, fields)
    write_csv(out / "landmarks.csv", landmarks, ["accession", "landmark_id", "label", "start", "end", "source_note"])
    aligned = subprocess.run(["mafft", "--quiet", "--thread", "2", "--auto", str(out / "panel.fasta")],
                             capture_output=True, text=True, check=True)
    (out / "alignment.fasta").write_text(aligned.stdout)
    names = ("panel.fasta", "panel_manifest.csv", "landmarks.csv", "alignment.fasta")
    manifest = dict(schema_version="1.0.0", source="NCBI GenBank and NCBI Taxonomy",
        acquired_at=max(acquisition_time(raw), acquisition_time(raw_tax)),
        acquisition_time_basis="original raw-response cache file modification time; proxy when no first-download log is available",
        built_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        sampling_policy="bounded cached NCBI search per named phylum; not random, exhaustive, or guaranteed newest",
        panel_count=len(sequences), previous_count=len(previous_accessions),
        retained_previous_count=len(set(previous_accessions) & set(sequences)),
        added_count=len(set(sequences) - set(previous_accessions)), searches=searches,
        previous_run_manifest_sha256=previous_run_manifest_sha,
        raw_responses={raw.name:sha(raw), raw_tax.name:sha(raw_tax)},
        raw_response_acquired_at={raw.name:acquisition_time(raw), raw_tax.name:acquisition_time(raw_tax)},
        taxonomy_snapshot_sha256=taxonomy_source_sha,
        mafft_version=subprocess.run(["mafft", "--version"], capture_output=True, text=True).stderr.strip(),
        alignment_policy="full alignment without keeplength; product lengths use raw ungapped templates; cross-clade variable-region homology not assumed",
        artifacts_sha256={name:sha(out / name) for name in names})
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(panel_count=len(sequences), lengths=[min(map(len,sequences.values())),max(map(len,sequences.values()))],
        classes=sorted(set(x["class"] for x in metadata)), output=str(out)), indent=2))


if __name__ == "__main__":
    main()
