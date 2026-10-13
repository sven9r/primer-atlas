import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import annotate_unite_regions as annotator
import scan_unite_primer_sites as scanner


def write_gzip(path,text):
    with path.open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as stream:
        stream.write(text.encode())

def write_csv(path,rows):
    text=io.StringIO();writer=csv.DictWriter(text,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if path.suffix=='.gz':write_gzip(path,text.getvalue())
    else:path.write_text(text.getvalue())

def position_line(key,length,bounds,note=''):
    return '\t'.join([key,str(length)+' bp.']+[region+': '+('-'.join(map(str,bounds[region])) if region in bounds else 'Not found') for region in annotator.REGIONS]+[note])

def detail_line(key,length,strand=1,chimera=''):
    return '\t'.join([key,str(length),'F','1' if strand==-1 else '0','4','1e-10','50','50','1',str(length),'1_SSU','4_LSU',chimera,'1_SSU,2_5.8,3_End,4_LSU','F: 4 1e-10 50'])

class AnnotationFixtures(unittest.TestCase):
    def make_fixture(self,root):
        panel=root/'panel';panel.mkdir();work=root/'work';work.mkdir()
        oriented='ACGTA'+'GATTACA'+'ACGTRYSWKMBD'+'TTTAGATTACATGCAA'+'GCGCGCTA'
        n=len(oriented)
        seqs={'F':oriented,'R':scanner.reverse_complement(oriented),'P':oriented,'C':oriented,'B':oriented,'U':oriented}
        text=''.join('>'+key+'\n'+seq+'\n' for key,seq in seqs.items())
        write_gzip(panel/'selected_raw.fasta.gz',text)
        write_csv(panel/'selected_metadata.csv.gz',[dict(record_id=key,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()) for key,seq in seqs.items()])
        (work/'panel.fasta').write_text(text)
        files={}
        for name in ('itsx_program','itsx_wrapper','hmmsearch','fungal_profile'):
            file=root/name;file.write_text('fixture '+name);files[name]=dict(path=str(file),sha256=annotator.sha(file))
        tool=root/'tool.json';tool.write_text(json.dumps(dict(itsx_version='1.1.3',hmmer_version='3.4',files=files)))
        for name,path in [('annotator',ROOT/'scripts/annotate_unite_regions.py'),('scanner',ROOT/'scripts/scan_unite_primer_sites.py')]:
            target=work/('old_'+name+'.py');target.write_bytes(path.read_bytes())
        run=dict(schema_version='1.0.0',capture_mode='active_process_verified_before_edit',work_directory=str(work),
            input_sha256=annotator.sha(panel/'selected_raw.fasta.gz'),staged_input_sha256=annotator.sha(work/'panel.fasta'),n_input=len(seqs),
            tool_manifest_sha256=annotator.sha(tool),semantic_parameters=dict(annotator.SEMANTICS,cpu=6),
            effective_command=[str(root/'itsx_wrapper'),'-i','panel.fasta','-o','regions','-t','F','--preserve','T','--save_regions','all','--detailed_results','T','--save_raw','T','--cpu','6','-p',str(root),'--temp',str(work/'temp')],
            executed_source={name+suffix:value for name in ('annotator','scanner') for suffix,value in (('_path',str(work/('old_'+name+'.py'))),('_sha256',annotator.sha(work/('old_'+name+'.py'))))})
        runfile=root/'run.json';runfile.write_text(json.dumps(run))
        full={'SSU':(1,5),'ITS1':(6,12),'5.8S':(13,24),'ITS2':(25,n-8),'LSU':(n-7,n)}
        partial={'ITS1':(1,12),'5.8S':(13,24),'ITS2':(25,n)}
        specs=[('F',full,'',1),('R',full,'',-1),('P',partial,'',1),('C',full,'Chimeric! ',1),('B',full,'Broken or partial sequence, only partial 5.8S! ',1)]
        pos=[];details=[];regional={region:[] for region in annotator.REGIONS}
        for i,(key,bounds,note,strand) in enumerate(specs):
            if i==2:
                pos.append('--END--\t bp.\tSSU: Not found\tITS1: Not found\t5.8S: Not found\tITS2: Not found\tLSU: Not found\tBroken or partial sequence, no 5.8S! ')
                details.append('--END--\t\t\t\t\t\t\t\t1\t\t***\t***\t\t\t:   ')
            pos.append(position_line(key,n,bounds,note))
            details.append(detail_line(key,n,strand,'Chimeric' if key=='C' else ''))
            if key=='C':continue
            for region,(start,end) in bounds.items():
                seq=oriented[start-1:end]
                if key=='R' and region=='5.8S':seq=seq.translate(str.maketrans('SW','WS'))
                regional[region].append('>'+key+'\n'+seq+'\n')
        (work/'regions.positions.txt').write_text('\n'.join(pos)+'\n')
        (work/'regions.extraction.results').write_text('\n'.join(details)+'\n')
        (work/'regions.summary.txt').write_text('ITSx run started at fixture.\nNumber of sequences in input file:\t6\nSequences detected as ITS by ITSx:\t4\nOn main strand:\t4\nOn complementary strand:\t0\nSequences detected as chimeric by ITSx:\t1\nFungi:\t3\nITSx run finished at fixture.\n')
        raw=work/'regions_ITSx_raw_output';raw.mkdir()
        (raw/'main.1.fasta').write_text(''.join('>'+key+' main\n'+seqs[key]+'\n' for key in ('F','R')))
        (raw/'main.2.fasta').write_text(''.join('>'+key+' main\n'+seqs[key]+'\n' for key in ('P','C','B','U')))
        for region,rows in regional.items():(work/('regions.'+region.replace('.','_')+'.fasta')).write_text(''.join(rows))
        return panel,work,tool,runfile,seqs

    def test_verified_reparse_preserves_orientation_partial_notes_and_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,work,tool,run,seqs=self.make_fixture(root)
            state=annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)
            self.assertEqual(state['n_input'],6);self.assertEqual(state['n_detected'],5)
            self.assertEqual(state['n_reverse'],1);self.assertEqual(state['n_chimeric'],1)
            self.assertEqual(state['summary_proof']['ignored_exact_batch_sentinels'],1)
            output=root/'output';regions=scanner.read_csv(output/'regions.csv.gz')
            rows={(row['record_id'],row['region']):row for row in regions}
            self.assertEqual(rows[('F','ITS1')]['boundary_status'],'complete')
            self.assertEqual(rows[('P','ITS1')]['boundary_status'],'partial')
            self.assertEqual(rows[('C','ITS1')]['boundary_status'],'chimeric_withheld')
            self.assertEqual(rows[('B','ITS1')]['boundary_status'],'partial')
            self.assertIn('Broken',rows[('B','ITS1')]['annotation_note'])
            self.assertEqual(rows[('R','5.8S')]['boundary_status'],'extraction_mismatch')
            self.assertEqual(rows[('R','ITS1')]['original_start'],str(len(seqs['R'])-12+1))
            orientations={row['record_id']:row for row in scanner.read_csv(output/'orientations.csv.gz')}
            self.assertEqual(orientations['U']['orientation'],'source_unverified')
            oriented=dict(scanner.read_fasta(output/'oriented_raw.fasta.gz'))
            self.assertEqual(oriented['R'],seqs['F']);self.assertEqual(set(oriented),set(seqs))
            for name in ('code_region_parser.py','code_parser_scanner_dependency.py','code_profile_run_annotator.py','code_profile_run_scanner.py','tool_identity_manifest.json','profile_run_identity.json'):
                self.assertEqual(annotator.sha(output/name),state['artifacts'][name]['sha256'])
            self.assertEqual(annotator.build(panel,root/'itsx_wrapper',output,6,work,tool,run),state)
            (output/'ITS1.fasta.gz').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'checksum'):annotator.build(panel,root/'itsx_wrapper',output,6,work,tool,run)

    def test_unverified_stale_work_and_semantic_or_tool_changes_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,work,tool,run,_=self.make_fixture(root)
            before=(work/'panel.fasta').read_bytes()
            with self.assertRaisesRegex(ValueError,'verified run identity'):
                annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool)
            self.assertEqual((work/'panel.fasta').read_bytes(),before)
            with self.assertRaisesRegex(ValueError,'semantic identity'):
                annotator.build(panel,root/'itsx_wrapper',root/'output',2,work,tool,run)
            (root/'fungal_profile').write_text('mutated model')
            with self.assertRaisesRegex(ValueError,'Tool identity checksum'):
                annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)

    def test_unknown_blank_strand_foreign_owner_and_bad_sentinels_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,work,tool,run,_=self.make_fixture(root)
            path=work/'regions.extraction.results';original=path.read_text()
            path.write_text(original.replace('F\t48\tF\t0','F\t48\tF\t'))
            with self.assertRaisesRegex(ValueError,'strand flag'):annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)
            path.write_text(original.replace('--END--','UNKNOWN'))
            with self.assertRaisesRegex(ValueError,'strand flag'):annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)
            path.write_text(original.replace('--END--\t\t\t\t','--END--\t\t\t0\t'))
            with self.assertRaisesRegex(ValueError,'sentinel'):annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)

    def test_invalid_anchor_cannot_establish_complete_region(self):
        self.assertEqual(annotator.boundary_status('ITS1',{'ITS1':(6,12),'SSU':(1,999),'5.8S':(13,24)},48),'partial')
        self.assertEqual(annotator.boundary_status('ITS2',{'ITS2':(25,40),'5.8S':(13,24),'LSU':(45,48)},48),'partial')
        self.assertEqual(annotator.boundary_status('ITS1',{'ITS1':(0,12)},48),'invalid_coordinates')

    def test_actual_command_override_and_same_ids_changed_bases_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,work,tool,run,seqs=self.make_fixture(root)
            original=json.loads(run.read_text())
            forged=copy.deepcopy(original);forged['effective_command']+=['--complement','F']
            run.write_text(json.dumps(forged))
            with self.assertRaisesRegex(ValueError,'biological command semantics'):
                annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)
            run.write_text(json.dumps(original))
            seqs['F']='T'+seqs['F'][1:]
            write_gzip(panel/'selected_raw.fasta.gz',''.join('>'+key+'\n'+seq+'\n' for key,seq in seqs.items()))
            write_csv(panel/'selected_metadata.csv.gz',[dict(record_id=key,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()) for key,seq in seqs.items()])
            before=(work/'panel.fasta').read_bytes()
            with self.assertRaisesRegex(ValueError,'input/tool/semantic identity'):
                annotator.build(panel,root/'itsx_wrapper',root/'output',6,work,tool,run)
            self.assertEqual((work/'panel.fasta').read_bytes(),before)

class ScanProofFixtures(unittest.TestCase):
    def fixture(self,root):
        panel=root/'panel';panel.mkdir()
        f='TCCGTAGGTGAACCTGCGG';r='TCCTCCGCTTATTGATATGC'
        seqs={'A':f+'N'*105+scanner.reverse_complement(r),'B':'A'*150}
        write_gzip(panel/'selected_raw.fasta.gz',''.join('>'+key+'\n'+seq+'\n' for key,seq in seqs.items()))
        write_csv(panel/'selected_metadata.csv.gz',[dict(record_id=key) for key in seqs])
        catalog=root/'catalog.csv';pairs=root/'pairs.csv'
        write_csv(catalog,[dict(catalog_id='F',primer_name='F',sequence=f,direction='forward'),dict(catalog_id='R',primer_name='R',sequence=r,direction='reverse')])
        write_csv(pairs,[dict(pair_id='F_R',pair_label='F + R',forward_primer='F',forward_sequence=f,reverse_primer='R',reverse_sequence=r)])
        return panel,catalog,pairs

    def test_completed_sparse_scan_proves_every_search_and_is_immutable(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,catalog,pairs=self.fixture(root)
            output=root/'output';state=scanner.build(panel,catalog,pairs,output)
            self.assertTrue(state['complete']);self.assertEqual(state['n_site_searches'],4)
            self.assertEqual(state['n_sparse_rows']+state['n_no_exact_hits'],4)
            self.assertEqual(state['record_ids_sha256'],hashlib.sha256(b'A\nB').hexdigest())
            for name,key in [('primer_sites.csv.gz','site_artifact_sha256'),('catalog.csv.gz','catalog_artifact_sha256'),('pairs.csv.gz','pair_artifact_sha256')]:
                self.assertEqual(scanner.sha(output/name),state[key])
            result=scanner.read_csv(output/'pairs.csv.gz')[0]
            self.assertEqual(result['n_scorable'],'1');self.assertEqual(result['amplicon_median_bp'],'144')
            self.assertEqual(scanner.build(panel,catalog,pairs,output),state)
            (output/'primer_sites.csv.gz').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'checksum'):scanner.build(panel,catalog,pairs,output)

    def test_gapped_templates_and_missing_panel_members_cannot_emit_proof(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);panel,catalog,pairs=self.fixture(root)
            raw=dict(scanner.read_fasta(panel/'selected_raw.fasta.gz'))
            raw['A']=raw['A'][:20]+'-'*100+raw['A'][20:]
            write_gzip(panel/'selected_raw.fasta.gz',''.join('>'+key+'\n'+seq+'\n' for key,seq in raw.items()))
            with self.assertRaisesRegex(ValueError,'ungapped'):scanner.build(panel,catalog,pairs,root/'output')
            self.assertFalse((root/'output/site_search_manifest.json').exists())
            write_gzip(panel/'selected_raw.fasta.gz','>A\nACGT\n')
            with self.assertRaisesRegex(ValueError,'exactly the retained'):scanner.build(panel,catalog,pairs,root/'output')

if __name__=='__main__':unittest.main()
