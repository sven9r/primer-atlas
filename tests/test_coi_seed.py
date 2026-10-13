import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('coi_seed', Path(__file__).resolve().parents[1] / 'scripts/restore_coi_seed.py')
seed = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(seed)


class SeedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / 'source'
        self.panel = self.source / seed.PANEL_PREFIX
        self.panel.mkdir(parents=True)
        self.groups = ['OrderA', 'OrderB']
        source = self.panel / 'raw/content/source.xml'
        source.parent.mkdir(parents=True)
        source.write_text('<fixture>source response bytes</fixture>')
        source_hash = seed.sha(source)
        taxonomy = source.with_name(source_hash + '.xml')
        taxonomy.write_bytes(source.read_bytes())
        records = []
        fasta = []
        for i, group in enumerate(self.groups):
            accession = 'AA00000' + str(i + 1) + '.1'
            header = accession + ' frozen original header'
            sequence = 'ACGT' * 150
            records.append(dict(target_group=group, accession=accession, header=header,
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                first_release_id='fixture-1', added_at_utc='2026-10-11T00:00:00Z',
                qc_status='verified_fixture_only', source_sha256=source_hash,
                source_file='raw/content/source.xml', taxonomy_source_sha256=source_hash))
            fasta += ['>' + header, sequence]
        (self.panel / 'panel.fasta').write_text('\n'.join(fasta) + '\n')
        for name in ('panel_manifest.csv', 'ledger.csv', 'taxonomy.csv'):
            self.csv(self.panel / name, records)
        self.csv(self.panel / 'growth_qa.csv', [dict(target_group=group, target_met=True) for group in self.groups])
        self.csv(self.panel / 'rejections.csv', [dict(reason='fixture only')])
        artifacts = {name: dict(path=name, sha256=seed.sha(self.panel / name)) for name in (
            'panel.fasta','panel_manifest.csv','ledger.csv','taxonomy.csv','growth_qa.csv','rejections.csv')}
        (self.panel / 'run_manifest.json').write_text(json.dumps(dict(marker_id='COI',schema_version='1.0.0',
            groups=self.groups,retained_count=2,verified_marker_count=2,artifacts=artifacts,legacy_archives=[])))
        baseline = self.source / seed.BASELINE_PATH
        baseline.parent.mkdir(parents=True)
        self.csv(baseline,[dict(target_group=records[0]['target_group'],accession=records[0]['accession'])])

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def csv(path, records):
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader(); writer.writerows(records)

    def archive(self, extra=None):
        path = self.root / 'seed.tar.gz'
        with tarfile.open(path,'w:gz') as bundle:
            for source in sorted(self.source.rglob('*')):
                if source.is_file(): bundle.add(source,arcname=str(source.relative_to(self.source)))
            if extra:
                member, payload = extra
                bundle.addfile(member,io.BytesIO(payload) if payload else None)
        return path

    def install(self, path, expected=None):
        target = self.root / 'worker'
        target.mkdir(exist_ok=True)
        validator = seed.validate_panel
        # Smaller fixtures exercise the same checks; the public CLI retains
        # the fixed100/group and10000verified total production thresholds.
        with patch.object(seed,'validate_panel',lambda panel,groups: validator(panel,groups,2,1)):
            return seed.install_seed(path,expected or seed.sha(path),target,self.groups)

    def test_verified_install_preserves_all_bytes_and_baseline(self):
        path = self.archive()
        proof = self.install(path)
        self.assertEqual(proof['retained_count'],2)
        self.assertEqual(proof['published_baseline_count'],1)
        installed = self.root / 'worker' / seed.PANEL_PREFIX
        for source in self.panel.rglob('*'):
            if source.is_file(): self.assertEqual(seed.sha(source),seed.sha(installed / source.relative_to(self.panel)))
        with self.assertRaises(ValueError): self.install(path)

    def test_archive_checksum_rejected_before_install(self):
        with self.assertRaises(ValueError): self.install(self.archive(),'0' * 64)
        self.assertFalse((self.root / 'worker' / seed.PANEL_PREFIX).exists())

    def test_traversal_and_outside_path_rejected(self):
        for name in ('../escape.txt','/absolute.txt','data/derived/coi_reference_panel/../escape.txt','data/provenance/other.csv'):
            member = tarfile.TarInfo(name);member.size=1
            with self.assertRaises(ValueError): self.install(self.archive((member,b'x')))
        self.assertFalse((self.root / 'worker' / seed.PANEL_PREFIX).exists())

    def test_symlink_rejected(self):
        member = tarfile.TarInfo('data/derived/coi_reference_panel/raw/link')
        member.type=tarfile.SYMTYPE;member.linkname='/tmp/outside'
        with self.assertRaises(ValueError): self.install(self.archive((member,b'')))

    def test_corrupt_artifact_and_original_response_rejected(self):
        original = (self.panel / 'panel.fasta').read_bytes()
        (self.panel / 'panel.fasta').write_bytes(original + b'corrupt')
        with self.assertRaises(ValueError): self.install(self.archive())
        (self.panel / 'panel.fasta').write_bytes(original)
        (self.panel / 'raw/content/source.xml').write_text('changed source response')
        with self.assertRaises(ValueError): self.install(self.archive())

    def test_production_policy_floor_is_enforced(self):
        with self.assertRaises(ValueError): seed.validate_panel(self.panel,self.groups)

    def test_published_membership_loss_rejected(self):
        baseline = self.source / seed.BASELINE_PATH
        self.csv(baseline,[dict(target_group='OrderA',accession='MISSING999.1')])
        with self.assertRaises(ValueError): self.install(self.archive())
        self.assertFalse((self.root / 'worker' / seed.PANEL_PREFIX).exists())


if __name__ == '__main__':
    unittest.main()
