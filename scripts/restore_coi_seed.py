#!/usr/bin/env python3
"""Install a checksummed manual COI seed without acquiring or publishing data."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile

PANEL_PREFIX = PurePosixPath('data/derived/coi_reference_panel')
BASELINE_PATH = PurePosixPath('data/provenance/coi_published_baseline.csv')
MAX_ARCHIVE_BYTES = 4 * 1024**3
MAX_ARCHIVE_ENTRIES = 50000


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def safe_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or '\\' in value:
        raise ValueError('Unsafe COI archive path')
    return path


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def fasta(path):
    result, header, parts = {}, None, []
    def save():
        if header is None:
            return
        token = re.sub(r'^_R_', '', header.split()[0])
        if not re.fullmatch(r'[A-Z][A-Z_]*[0-9]+\.[0-9]+', token) or token in result:
            raise ValueError('Seed FASTA needs unique accession versions')
        sequence = ''.join(parts)
        if not sequence or not re.fullmatch(r'[ACGTURYSWKMBDHVNI]+', sequence):
            raise ValueError('Invalid seed sequence')
        result[token] = (header, sequence)
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            save(); header, parts = line[1:], []
            if not header.strip():
                raise ValueError('Missing seed FASTA header')
        elif line.strip():
            if header is None:
                raise ValueError('Seed sequence lacks a header')
            parts.append(line.strip())
    save()
    return result


def validate_panel(panel, groups, minimum_total=10000, minimum_per_group=100):
    hashes = {}
    def checked_hash(path):
        path = Path(path)
        if path not in hashes:
            hashes[path] = sha(path)
        return hashes[path]
    content_files = {}
    for path in (panel / 'raw/content').glob('*'):
        if path.is_file():
            content_files.setdefault(path.stem, []).append(path)
    required = {'panel.fasta', 'panel_manifest.csv', 'ledger.csv', 'taxonomy.csv',
                'growth_qa.csv', 'rejections.csv', 'run_manifest.json'}
    if not all((panel / name).is_file() for name in required):
        raise ValueError('Incomplete COI seed panel')
    run = json.loads((panel / 'run_manifest.json').read_text())
    if run.get('marker_id') != 'COI' or run.get('schema_version') != '1.0.0':
        raise ValueError('Invalid COI seed manifest')
    if set(run.get('groups', [])) != set(groups):
        raise ValueError('Seed target groups differ from the active arthropod groups')
    artifacts = run.get('artifacts', {})
    if not required.difference({'run_manifest.json'}) <= set(artifacts):
        raise ValueError('Seed manifest is missing artifact checksums')
    for name, artifact in artifacts.items():
        relative = safe_path(artifact.get('path', name))
        path = panel / relative
        if not path.is_file() or checked_hash(path) != artifact.get('sha256'):
            raise ValueError('COI seed artifact checksum differs: ' + name)
    sequences = fasta(panel / 'panel.fasta')
    records = rows(panel / 'panel_manifest.csv')
    if len(records) != len(sequences) or {row['accession'] for row in records} != set(sequences):
        raise ValueError('Seed sequence and metadata identities differ')
    counts = dict.fromkeys(groups, 0)
    verified_counts = dict.fromkeys(groups, 0)
    for row in records:
        group = row['target_group']
        if group not in counts:
            raise ValueError('Unknown seed target group')
        header, sequence = sequences[row['accession']]
        if hashlib.sha256(sequence.encode()).hexdigest() != row.get('sequence_sha256') or hashlib.sha256(header.encode()).hexdigest() != row.get('header_sha256'):
            raise ValueError('Seed sequence/header byte identity differs')
        counts[group] += 1
        verified_counts[group] += row.get('qc_status', '').startswith('verified')
        source_hash = row.get('source_sha256', '')
        if not re.fullmatch(r'[a-f0-9]{64}', source_hash):
            raise ValueError('Seed source response hash is missing')
        source_file = panel / safe_path(row.get('source_file', ''))
        if not source_file.is_file() or checked_hash(source_file) != source_hash:
            raise ValueError('Seed original source response differs')
        taxonomy_hash = row.get('taxonomy_source_sha256', '')
        if taxonomy_hash:
            candidates = content_files.get(taxonomy_hash, [])
            if not candidates or not any(checked_hash(path) == taxonomy_hash for path in candidates):
                raise ValueError('Seed taxonomy source response is missing or corrupt')
    if len(records) != run.get('retained_count') or len(records) < minimum_total or min(counts.values()) < minimum_per_group:
        raise ValueError('Seed retained count or group minimum failed')
    if sum(verified_counts.values()) != run.get('verified_marker_count') or sum(verified_counts.values()) < minimum_total or min(verified_counts.values()) < minimum_per_group:
        raise ValueError('Seed verified-marker policy failed')
    for name in ('ledger.csv', 'taxonomy.csv'):
        metadata = rows(panel / name)
        if len(metadata) != len(records) or {row['accession'] for row in metadata} != set(sequences):
            raise ValueError('Seed ledger/taxonomy accession coverage differs')
        by_id = {row['accession']: row for row in metadata}
        for record in records:
            for column in ('target_group', 'sequence_sha256', 'header_sha256', 'first_release_id', 'added_at_utc'):
                if record.get(column) != by_id[record['accession']].get(column):
                    raise ValueError('Seed frozen metadata differs: ' + column)
    for legacy in run.get('legacy_archives', []):
        path = panel / safe_path(legacy['path'])
        if not path.is_file() or checked_hash(path) != legacy['sha256']:
            raise ValueError('Seed legacy archive bytes differ')
    return dict(retained_count=len(records), group_counts=counts,
                verified_marker_count=sum(verified_counts.values()), run_manifest_sha256=checked_hash(panel / 'run_manifest.json'))


def install_seed(archive, expected_sha, root, groups):
    if not re.fullmatch(r'[a-fA-F0-9]{64}', expected_sha) or sha(archive) != expected_sha.lower():
        raise ValueError('COI seed archive SHA-256 differs')
    root = Path(root).resolve()
    parent = root / 'data/derived'
    parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix='coi_seed_', dir=parent))
    try:
        with tarfile.open(archive, 'r:gz') as bundle:
            members = bundle.getmembers()
            if len(members) > MAX_ARCHIVE_ENTRIES or any(member.size < 0 for member in members) or sum(member.size for member in members) > MAX_ARCHIVE_BYTES:
                raise ValueError('COI seed archive exceeds its bounded size')
            seen = set()
            for member in members:
                relative = safe_path(member.name)
                if relative in seen:
                    raise ValueError('Repeated COI seed archive path')
                seen.add(relative)
                # Parent directory entries are harmless; every file must live
                # in the one designated panel subtree, and links are forbidden.
                if member.isdir():
                    if relative != PANEL_PREFIX and relative not in PANEL_PREFIX.parents and PANEL_PREFIX not in relative.parents and relative not in BASELINE_PATH.parents:
                        raise ValueError('Unexpected COI seed archive directory')
                    continue
                if not member.isfile() or (PANEL_PREFIX not in relative.parents and relative != BASELINE_PATH):
                    raise ValueError('COI seed archive contains a link, special file or unexpected path')
                target = scratch / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = bundle.extractfile(member)
                if stream is None:
                    raise ValueError('Missing COI seed file payload')
                with stream, target.open('wb') as destination:
                    shutil.copyfileobj(stream, destination)
                os.utime(target, (member.mtime, member.mtime))
        panel = scratch / PANEL_PREFIX
        proof = validate_panel(panel, groups)
        baseline = scratch / BASELINE_PATH
        if baseline.is_file():
            panel_records = rows(panel / 'panel_manifest.csv')
            panel_keys = {(row['target_group'], row['accession']) for row in panel_records}
            published = rows(baseline)
            if not published or not all((row['target_group'], row['accession']) in panel_keys for row in published):
                raise ValueError('Seed dropped a previously published accession identity')
            proof['published_baseline_count'] = len(published)
            proof['published_baseline_sha256'] = sha(baseline)
            existing_baseline = root / BASELINE_PATH
            if existing_baseline.exists() and sha(existing_baseline) != sha(baseline):
                raise ValueError('Refusing to replace a different published COI baseline')
        destination = root / PANEL_PREFIX
        if destination.exists():
            raise ValueError('Refusing to overwrite an existing COI panel with a manual seed')
        os.replace(panel, destination)
        if baseline.is_file() and not (root / BASELINE_PATH).exists():
            (root / BASELINE_PATH).parent.mkdir(parents=True, exist_ok=True)
            os.replace(baseline, root / BASELINE_PATH)
        proof['archive_sha256'] = expected_sha.lower()
        evidence = root / 'data/provenance/coi_seed_restore.json'
        evidence.parent.mkdir(parents=True, exist_ok=True)
        evidence.write_text(json.dumps(proof, indent=2) + '\n')
        return proof
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--root', default='.')
    args = parser.parse_args()
    root = Path(args.root).resolve()
    groups = [row['order'] for row in rows(root / 'data/orders_22.csv') if row['order'] != 'Bivalvia']
    proof = install_seed(args.archive, args.sha256, root, groups)
    print('Verified manual COI seed:', proof['retained_count'], 'retained records;', proof['verified_marker_count'], 'verified markers')


if __name__ == '__main__':
    main()
