#!/usr/bin/env python3
"""Extract profile-supported ITS regions from verified immutable ITSx output."""
from __future__ import annotations
import argparse
from collections import Counter
import datetime
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from scan_unite_primer_sites import (read_fasta, read_csv, reverse_complement,
                                    validate_raw_records, write_gzip_csv)

REGIONS = ("SSU", "ITS1", "5.8S", "ITS2", "LSU")
SEMANTICS = dict(profile="F", preserve=True, save_regions="all", detailed_results=True,
                 complement=True, output_prefix="regions", min_domains=2,
                 allow_single_domain="1e-9,0", anchor=0, only_full=False, algorithm="hmmsearch")
PARSER_RULE_VERSION = "itsx_regions_verified_oriented_slices_v2"

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def bytes_sha(value): return hashlib.sha256(value).hexdigest()
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def canonical_sha(value): return bytes_sha(json.dumps(value,sort_keys=True,separators=(",",":")).encode())

def write_gzip_text(path, lines):
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as compressed:
            with io.TextIOWrapper(compressed,encoding="utf-8") as handle:
                for line in lines: handle.write(line)

def verified_output_fields(path,kind,expected_ids=None,batch_by_id=None):
    fields=[line.split("\t") for line in path.read_text().splitlines() if line.strip()]
    sentinel_details=["--END--","","","","","","","","1","","***","***","","",":   "]
    sentinel_positions=["--END--"," bp.","SSU: Not found","ITS1: Not found","5.8S: Not found","ITS2: Not found","LSU: Not found","Broken or partial sequence, no 5.8S! "]
    sentinel=sentinel_positions if kind=="positions" else sentinel_details
    lookup={key:i for i,key in enumerate(expected_ids or [])}
    seen_boundaries=set(); output=[]
    for i,row in enumerate(fields):
        if row[0]!="--END--": output.append(row);continue
        if row!=sentinel or not lookup or i==0 or i==len(fields)-1:
            raise ValueError("Unverified ITSx sentinel row")
        before,after=fields[i-1][0],fields[i+1][0]
        if before not in lookup or after not in lookup:
            raise ValueError("ITSx sentinel lacks owned batch neighbors")
        if batch_by_id is None or before not in batch_by_id or after not in batch_by_id:
            raise ValueError("ITSx sentinel lacks verified source batch ownership")
        batch=batch_by_id[before]
        if batch_by_id[after]!=batch+1 or batch in seen_boundaries:
            raise ValueError("ITSx sentinel is not at a unique internal batch boundary")
        seen_boundaries.add(batch)
    return output,seen_boundaries

def parse_position_details(path,expected_ids=None,batch_by_id=None):
    rows={}
    fields_list,_=verified_output_fields(path,"positions",expected_ids,batch_by_id)
    for fields in fields_list:
        if len(fields)<7: raise ValueError("ITSx position schema changed")
        record_id=fields[0].split()[0]
        if record_id in rows: raise ValueError("ITSx position identity repeats")
        length=re.fullmatch(r"(\d+) bp[.]",fields[1].strip())
        if not length: raise ValueError("ITSx input length is unavailable")
        intervals={}; statuses={}; notes=[]
        for field in fields[2:]:
            field=field.strip()
            match=re.fullmatch(r"(SSU|ITS1|5\.8S|ITS2|LSU): (.+)",field)
            if match:
                region,status=match.groups()
                if region in statuses: raise ValueError("ITSx position region repeats")
                statuses[region]=status
                bound=re.fullmatch(r"(\d+)-(\d+)",status)
                if bound: intervals[region]=tuple(map(int,bound.groups()))
            elif field: notes.append(field)
        if set(statuses)!=set(REGIONS): raise ValueError("ITSx position regions are missing")
        rows[record_id]=dict(intervals=intervals,statuses=statuses,input_bp=int(length[1]),notes=" ".join(notes))
    return rows

def parse_positions(path,expected_ids=None):
    return {key:row["intervals"] for key,row in parse_position_details(path,expected_ids).items()}

def parse_detailed_results(path,expected_ids=None,batch_by_id=None):
    rows={}
    fields_list,_=verified_output_fields(path,"details",expected_ids,batch_by_id)
    for fields in fields_list:
        if len(fields)<15: raise ValueError("ITSx detailed results schema changed")
        record_id=fields[0].split()[0]
        if fields[3] not in ("0","1"): raise ValueError("ITSx strand flag is not 0/1")
        if record_id in rows: raise ValueError("ITSx strand identity repeats")
        if fields[2]!="F": raise ValueError("ITSx output is not from the fungal profile")
        rows[record_id]=dict(strand=-1 if fields[3]=="1" else 1,input_bp=int(fields[1]),
                             chimera_flag=fields[12],domain_order=fields[13])
    return rows

def parse_strands(path,expected_ids=None):
    return {key:row["strand"] for key,row in parse_detailed_results(path,expected_ids).items()}

def valid_interval(interval,length):
    return interval is not None and all(value is not None for value in interval) and 1<=interval[0]<=interval[1]<=length

def boundary_status(region,intervals,length):
    if region not in intervals: return "not_detected"
    start,end=intervals[region]
    if not valid_interval((start,end),length): return "invalid_coordinates"
    if region in ("SSU","LSU"): return "flank_detected"
    def valid(other): return valid_interval(intervals.get(other),length)
    if region=="ITS1":
        complete=(valid("SSU") and valid("5.8S") and intervals["SSU"][1]+1==start and end+1==intervals["5.8S"][0])
    elif region=="ITS2":
        complete=(valid("5.8S") and valid("LSU") and intervals["5.8S"][1]+1==start and end+1==intervals["LSU"][0])
    else:
        complete=(start>1 and end<length and valid("ITS1") and valid("ITS2") and
                  intervals["ITS1"][1]+1==start and end+1==intervals["ITS2"][0])
    return "complete" if complete else "partial"

def load_tool_identity(path,itsx):
    path=Path(path)
    payload=path.read_bytes(); state=json.loads(payload)
    if state.get("itsx_version")!="1.1.3" or state.get("hmmer_version")!="3.4":
        raise ValueError("The verified ITSx1.1.3/HMMER3.4 tool identity is required")
    required={"itsx_program","itsx_wrapper","hmmsearch","fungal_profile"}
    if not required<=set(state.get("files",{})): raise ValueError("Tool identity lacks required files")
    for name,descriptor in state["files"].items():
        if sha(descriptor["path"])!=descriptor["sha256"]: raise ValueError("Tool identity checksum differs: "+name)
    profile=Path(state["files"]["fungal_profile"]["path"])
    for suffix in ("h3f","h3i","h3m","h3p"):
        pressed=Path(str(profile)+"."+suffix)
        if pressed.exists() and "fungal_profile_"+suffix not in state["files"]:
            raise ValueError("Pressed fungal HMM profile is not included in the tool identity")
        if pressed.exists() and Path(state["files"]["fungal_profile_"+suffix]["path"]).resolve()!=pressed.resolve():
            raise ValueError("Pressed fungal HMM profile ownership differs")
    if "profile_directory" in state and (Path(state["profile_directory"])/"F.hmm").resolve()!=profile.resolve():
        raise ValueError("Effective profile directory differs from verified fungal profile")
    launcher=Path(itsx).resolve()
    permitted={Path(state["files"][name]["path"]).resolve() for name in ("itsx_program","itsx_wrapper")}
    if launcher not in permitted: raise ValueError("ITSx launcher differs from verified tool identity")
    return state,payload

def validate_run_identity(state,work,source_hash,staged_hash,n_input,tool_hash,semantics,tool):
    if state.get("schema_version")!="1.0.0" or state.get("capture_mode") not in (
            "before_profile_execution","active_process_verified_before_edit"):
        raise ValueError("ITSx work has no recognized verified run provenance")
    if Path(state.get("work_directory","")).resolve()!=work.resolve(): raise ValueError("ITSx work ownership differs")
    checks=dict(input_sha256=source_hash,staged_input_sha256=staged_hash,n_input=n_input,
                tool_manifest_sha256=tool_hash,semantic_parameters=semantics)
    if any(state.get(key)!=value for key,value in checks.items()):
        raise ValueError("ITSx work input/tool/semantic identity differs")
    if sha(work/"panel.fasta")!=staged_hash: raise ValueError("ITSx staged input checksum differs")
    sources=state.get("executed_source",{})
    for name in ("annotator","scanner"):
        if name+"_path" not in sources or sha(sources[name+"_path"])!=sources.get(name+"_sha256"):
            raise ValueError("Original executed ITSx pipeline source is not preserved: "+name)
    command=state.get("effective_command",[])
    permitted={Path(tool["files"][name]["path"]).resolve() for name in ("itsx_program","itsx_wrapper")}
    prefix=command[:next((i for i,value in enumerate(command) if value.startswith("-")),len(command))]
    if not prefix or Path(prefix[-1]).resolve() not in permitted or (len(prefix)>1 and prefix[:-1]!=["perl"]):
        raise ValueError("Captured ITSx launcher differs from verified tools")
    options=command[len(prefix):]
    if len(options)%2 or len(set(options[::2]))!=len(options[::2]):
        raise ValueError("Captured ITSx options are malformed or repeated")
    values=dict(zip(options[::2],options[1::2]))
    biological_defaults={"-N":"2","-E":"1e-5","-S":"0","--allow_single_domain":"1e-9,0",
        "--complement":"T","--allow_reorder":"F","--selection_priority":"score",
        "--heuristics":"F","--nhmmer":"F","--anchor":"0","--require_anchor":"0",
        "--only_full":"F","--partial":"0","--truncate":"T","--concat":"F","--minlen":"0"}
    required={"-i":"panel.fasta","-o":"regions","-t":"F","--preserve":"T","--save_regions":"all",
        "--detailed_results":"T","--save_raw":"T","--cpu":str(semantics["cpu"])}
    harmless={"--graphical","--silent","--table","--fasta","--summary","--not_found","--positions"}
    if set(values)-set(biological_defaults)-set(required)-harmless-{"-p","--temp"}:
        raise ValueError("Captured ITSx command has unverified options")
    if any(values.get(flag,wanted)!=wanted for flag,wanted in biological_defaults.items()):
        raise ValueError("Captured ITSx biological command semantics differ")
    if "-p" not in values or Path(values["-p"]).resolve()!=Path(tool["files"]["fungal_profile"]["path"]).parent.resolve():
        raise ValueError("Captured ITSx profile directory differs")
    if "--temp" not in values or Path(values["--temp"]).resolve()!=(work/"temp").resolve():
        raise ValueError("Captured ITSx temporary work differs")
    for flag,wanted in (("-i","panel.fasta"),("-o","regions"),("-t","F"),("--preserve","T"),
                        ("--save_regions","all"),("--detailed_results","T"),("--save_raw","T"),("--cpu",str(semantics["cpu"]))):
        if flag not in command or command[command.index(flag)+1]!=wanted:
            raise ValueError("Captured ITSx command differs: "+flag)
    return state

def verified_batches(work,sequences):
    raw=work/"regions_ITSx_raw_output"
    files=sorted((path for path in raw.glob("main.*.fasta") if re.fullmatch(r"main[.]\d+[.]fasta",path.name)),
                 key=lambda path:int(path.name.split(".")[1]))
    if not files: raise ValueError("ITSx saved source batches are required for verified reuse")
    batches={};order=[];proof=[]
    for batch,path in enumerate(files):
        if int(path.name.split(".")[1])!=batch+1: raise ValueError("ITSx source batch numbering differs")
        records=list(read_fasta(path))
        if not records: raise ValueError("ITSx source batch is empty")
        for key,seq in records:
            if key in batches or key not in sequences or seq!=sequences[key]:
                raise ValueError("ITSx saved source batch ownership or bases differ")
            batches[key]=batch;order.append(key)
        proof.append(dict(file=path.name,sha256=sha(path),rows=len(records),first_id=records[0][0],last_id=records[-1][0]))
    if order!=list(sequences): raise ValueError("ITSx source batches do not retain the complete input order")
    return batches,proof

def check_finished_summary(path,expected_ids,details,sentinel_boundaries,batch_by_id):
    text=path.read_text()
    def number(label):
        match=re.search(re.escape(label)+r"\s*(\d+)",text)
        if not match: raise ValueError("ITSx completion summary lacks "+label)
        return int(match[1])
    if "ITSx run finished at" not in text: raise ValueError("ITSx run is not complete")
    if number("Number of sequences in input file:")!=len(expected_ids):
        raise ValueError("ITSx summary/input counts differ")
    # Original1.1.3 resets summary counters in each internal source batch. A
    # leaked previous-batch sentinel adds one blank-type main-strand call.
    last_batch=max(batch_by_id.values())
    last_ids={key for key,batch in batch_by_id.items() if batch==last_batch}
    last={key:row for key,row in details.items() if key in last_ids}
    spill=int(last_batch>0 and last_batch-1 in sentinel_boundaries)
    if number("Sequences detected as ITS by ITSx:")!=len(last)+spill or number("Fungi:")!=len(last):
        raise ValueError("ITSx summary/verified final-batch counts differ")
    if number("On main strand:")!=sum(row["strand"]==1 for row in last.values())+spill or number("On complementary strand:")!=sum(row["strand"]==-1 for row in last.values()):
        raise ValueError("ITSx summary strand counts differ")
    return dict(reported_detected=number("Sequences detected as ITS by ITSx:"),
                real_detected_all_batches=len(details),real_detected_final_batch=len(last),
                ignored_exact_batch_sentinels=len(sentinel_boundaries),
                scope="last_internal_source_batch" if last_batch else "all_records")

def build(panel,itsx,output,cpus=6,existing_work=None,tool_manifest=None,run_manifest=None):
    panel,itsx,output=Path(panel),Path(itsx),Path(output)
    tool_manifest=Path(tool_manifest or "tmp/unite-tools/tool_identity_manifest.json")
    # Freeze loaded algorithm/dependency/tool/input bytes before any long work.
    parser_bytes=Path(__file__).read_bytes()
    helper_path=Path(__file__).with_name("scan_unite_primer_sites.py")
    helper_bytes=helper_path.read_bytes()
    tool,tool_bytes=load_tool_identity(tool_manifest,itsx)
    source=panel/"selected_raw.fasta.gz"; metadata_path=panel/"selected_metadata.csv.gz"
    source_hash,metadata_hash=sha(source),sha(metadata_path)
    records=list(read_fasta(source)); validate_raw_records(records)
    sequences=dict(records); metadata=read_csv(metadata_path)
    ids=[row["record_id"] for row in metadata]
    if len(set(ids))!=len(ids) or set(sequences)!=set(ids): raise ValueError("UNITE metadata/sequence membership differs")
    for row in metadata:
        expected=row.get("sequence_sha256")
        if expected and expected!=bytes_sha(sequences[row["record_id"]].encode()): raise ValueError("UNITE original sequence hash differs")
    semantics=dict(SEMANTICS,cpu=int(cpus))
    if cpus<1: raise ValueError("ITSx CPU count must be positive")
    identity=dict(input_sha256=source_hash,metadata_sha256=metadata_hash,
        script_sha256=bytes_sha(parser_bytes),helper_sha256=bytes_sha(helper_bytes),
        tool_manifest_sha256=bytes_sha(tool_bytes),semantic_parameters=semantics,parser_rule_version=PARSER_RULE_VERSION)
    previous=output/"annotation_manifest.json"
    if previous.exists():
        state=json.loads(previous.read_text())
        if state.get("identity")!=identity: raise ValueError("Annotation identity differs; choose a new output directory")
        for descriptor in state["artifacts"].values():
            if sha(output/descriptor["path"])!=descriptor["sha256"]: raise ValueError("Annotation artifact checksum differs")
        return state
    if output.exists() and any(output.iterdir()): raise ValueError("Incomplete annotation output exists; choose a new directory")
    work=Path(existing_work) if existing_work else Path(tempfile.mkdtemp(prefix="primer-atlas-unite-annotation-"))
    if any(c.isspace() for c in str(work)): raise ValueError("ITSx work directory must not contain spaces")
    work.mkdir(parents=True,exist_ok=True)
    staged_bytes="".join(f">{record_id}\n{seq}\n" for record_id,seq in records).encode()
    staged_hash=bytes_sha(staged_bytes)
    provenance_path=Path(run_manifest) if run_manifest else work/"ITSx_run_identity.json"
    if provenance_path.exists():
        run_bytes=provenance_path.read_bytes(); run=json.loads(run_bytes)
        validate_run_identity(run,work,source_hash,staged_hash,len(records),bytes_sha(tool_bytes),semantics,tool)
    else:
        if run_manifest or any(work.iterdir()): raise ValueError("Existing ITSx work requires an explicitly verified run identity")
        (work/"panel.fasta").write_bytes(staged_bytes)
        (work/"code_annotator_executed.py").write_bytes(parser_bytes)
        (work/"code_scanner_executed.py").write_bytes(helper_bytes)
        profile_parent=tool.get("profile_directory",str(Path(tool["files"]["fungal_profile"]["path"]).parent))
        if any(c.isspace() for c in profile_parent): raise ValueError("ITSx profile path must not contain spaces")
        command=[str(itsx.resolve()),"-i","panel.fasta","-o","regions","-t","F","--preserve","T",
            "--save_regions","all","--detailed_results","T","--save_raw","T","--cpu",str(cpus),
            "--temp",str(work/"temp"),"-p",profile_parent]
        run=dict(schema_version="1.0.0",capture_mode="before_profile_execution",captured_at_utc=utc(),
            work_directory=str(work.resolve()),input_sha256=source_hash,staged_input_sha256=staged_hash,
            n_input=len(records),tool_manifest_sha256=bytes_sha(tool_bytes),semantic_parameters=semantics,
            effective_command=command,executed_source=dict(annotator_sha256=bytes_sha(parser_bytes),
                scanner_sha256=bytes_sha(helper_bytes),annotator_path=str((work/"code_annotator_executed.py").resolve()),
                scanner_path=str((work/"code_scanner_executed.py").resolve())))
        run_bytes=(json.dumps(run,indent=2)+"\n").encode();provenance_path.write_bytes(run_bytes)
    (work/"temp").mkdir(exist_ok=True)
    if not (work/"regions.positions.txt").exists():
        print(f"Annotating {len(sequences)} UNITE records with ITSx in {work}",flush=True)
        with (work/"ITSx.log").open("w") as log:
            subprocess.run(run["effective_command"],cwd=work,stdout=log,stderr=subprocess.STDOUT,check=True)
    required=[work/("regions"+suffix) for suffix in (".positions.txt",".extraction.results",".summary.txt",
                ".ITS1.fasta",".5_8S.fasta",".ITS2.fasta",".SSU.fasta",".LSU.fasta")]
    if any(not path.exists() for path in required): raise ValueError("ITSx did not produce required artifacts; inspect "+str(work/"ITSx.log"))
    expected_ids=list(sequences)
    batch_by_id,batch_proof=verified_batches(work,sequences)
    positions=parse_position_details(work/"regions.positions.txt",expected_ids,batch_by_id)
    details=parse_detailed_results(work/"regions.extraction.results",expected_ids,batch_by_id)
    if set(positions)!=set(details) or not set(positions)<=set(sequences): raise ValueError("ITSx output ownership differs")
    _,position_sentinels=verified_output_fields(work/"regions.positions.txt","positions",expected_ids,batch_by_id)
    _,detail_sentinels=verified_output_fields(work/"regions.extraction.results","details",expected_ids,batch_by_id)
    if position_sentinels!=detail_sentinels: raise ValueError("ITSx output sentinel identities differ")
    summary_proof=check_finished_summary(work/"regions.summary.txt",expected_ids,details,detail_sentinels,batch_by_id)
    extracted={}
    for region,suffix in (("SSU","SSU"),("ITS1","ITS1"),("5.8S","5_8S"),("ITS2","ITS2"),("LSU","LSU")):
        region_records=list(read_fasta(work/("regions."+suffix+".fasta")))
        if len({key for key,_ in region_records})!=len(region_records) or not {key for key,_ in region_records}<=set(positions):
            raise ValueError("ITSx extracted region ownership differs")
        extracted[region]=dict(region_records)
    oriented={};region_rows=[];orientation_rows=[];counts=Counter();issues=[]
    for record_id,seq in records:
        detected=record_id in details; detail=details.get(record_id,{})
        position=positions.get(record_id,dict(intervals={},statuses={},notes="",input_bp=len(seq)))
        if position["input_bp"]!=len(seq) or detail.get("input_bp",len(seq))!=len(seq): raise ValueError("ITSx record length differs from source")
        strand=detail.get("strand",1);analysis=reverse_complement(seq) if strand==-1 else seq
        note=position["notes"]
        chimeric=bool(detail.get("chimera_flag")) or "chimeric" in note.lower()
        broken=any(term in note.lower() for term in ("broken","overlap","too long"))
        annotation_status="chimeric_withheld" if chimeric else "broken_or_partial" if broken else "profile_detected" if detected else "not_detected"
        oriented[record_id]=analysis
        orientation_rows.append(dict(record_id=record_id,source_strand=strand,
            orientation="profile_forward" if detected else "source_unverified",annotation_status=annotation_status,
            annotation_note=note,source_sequence_sha256=bytes_sha(seq.encode()),analysis_sequence_sha256=bytes_sha(analysis.encode())))
        bounds=position["intervals"]
        for region in REGIONS:
            status=boundary_status(region,bounds,len(seq));interval=bounds.get(region,(None,None))
            if chimeric and region in bounds: status="chimeric_withheld"
            elif broken and status=="complete": status="partial"
            observed=extracted[region].get(record_id)
            if status in ("complete","partial","flank_detected"):
                expected=analysis[interval[0]-1:interval[1]]
                if observed!=expected:
                    status="extraction_mismatch"
                    issues.append(dict(record_id=record_id,region=region,reason="ITSx extraction differs from the correctly oriented raw slice; upstream S/W complement ambiguity is a known possible cause"))
            elif observed is not None and not chimeric:
                raise ValueError("ITSx extracted region has no valid reported source interval")
            counts[region+":"+status]+=1
            raw_start=raw_end=None
            if valid_interval(interval,len(seq)):
                raw_start,raw_end=(len(seq)-interval[1]+1,len(seq)-interval[0]+1) if strand==-1 else interval
            reason="Both bounding rRNA profiles detected" if status=="complete" and region in ("ITS1","ITS2") else status
            region_rows.append(dict(record_id=record_id,region=region,start=interval[0],end=interval[1],
                boundary_status=status,source_strand=strand,annotation_source="ITSx1.1.3 fungal HMM profiles",
                coordinate_basis="1-based inclusive positions in oriented source sequence",reason=reason,
                original_start=raw_start,original_end=raw_end,itsx_reported_status=position["statuses"].get(region,"No detection"),
                annotation_status=annotation_status,annotation_note=note))
    if sha(source)!=source_hash or sha(metadata_path)!=metadata_hash: raise ValueError("UNITE input changed during annotation")
    output.parent.mkdir(parents=True,exist_ok=True)
    staging=Path(tempfile.mkdtemp(prefix=".unite-annotation-",dir=output.parent))
    try:
        write_gzip_text(staging/"oriented_raw.fasta.gz",(f">{key}\n{value}\n" for key,value in oriented.items()))
        write_gzip_csv(staging/"regions.csv.gz",list(region_rows[0]),region_rows)
        write_gzip_csv(staging/"orientations.csv.gz",list(orientation_rows[0]),orientation_rows)
        write_gzip_csv(staging/"annotation_issues.csv.gz",["record_id","region","reason"],issues)
        for region in ("ITS1","5.8S","ITS2"):
            selected=[row for row in region_rows if row["region"]==region and row["boundary_status"]=="complete"]
            write_gzip_text(staging/(region.replace(".","_")+".fasta.gz"),
                (f'>{row["record_id"]}\n{oriented[row["record_id"]][row["start"]-1:row["end"]]}\n' for row in selected))
        for path in required: write_gzip_text(staging/(path.name+".gz"),iter(path.read_text().splitlines(keepends=True)))
        (staging/"code_region_parser.py").write_bytes(parser_bytes)
        (staging/"code_parser_scanner_dependency.py").write_bytes(helper_bytes)
        (staging/"tool_identity_manifest.json").write_bytes(tool_bytes)
        (staging/"profile_run_identity.json").write_bytes(run_bytes)
        for name in ("annotator","scanner"):
            shutil.copyfile(run["executed_source"][name+"_path"],staging/("code_profile_run_"+name+".py"))
        state=dict(identity=identity,input_sha256=source_hash,script_sha256=identity["script_sha256"],
            helper_sha256=identity["helper_sha256"],tool_manifest_sha256=identity["tool_manifest_sha256"],
            profile_run_manifest_sha256=bytes_sha(run_bytes),parser_rule_version=PARSER_RULE_VERSION,
            n_input=len(sequences),n_detected=len(positions),n_reverse=sum(row["strand"]==-1 for row in details.values()),
            n_chimeric=sum(row["annotation_status"]=="chimeric_withheld" for row in orientation_rows),
            region_counts=dict(sorted(counts.items())),summary_proof=summary_proof,source_batch_proof=batch_proof,
            command=run["effective_command"],work_directory=str(work),
            limitation="ITSx computational boundaries; SSU/LSU intervals are retained flanks, not complete genes. Complete regional alignments require valid adjacent anchors, unproblematic profile calls, and exact correctly oriented source slices. Original ITSx1.1.3 reverses S/W ambiguity incorrectly; mismatching outputs are withheld.",
            artifacts={path.name:dict(path=path.name,sha256=sha(path),bytes=path.stat().st_size) for path in sorted(staging.iterdir())})
        (staging/"annotation_manifest.json").write_text(json.dumps(state,indent=2)+"\n")
        if output.exists(): output.rmdir()
        staging.rename(output)
    finally:
        if staging.exists(): shutil.rmtree(staging)
    print(json.dumps({key:state[key] for key in ("n_input","n_detected","n_reverse","n_chimeric","region_counts")},indent=2),flush=True)
    return state

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--panel",type=Path,default=Path("data/derived/unite_reference_panel"))
    parser.add_argument("--itsx",type=Path,default=Path("tmp/unite-tools/bin/ITSx"))
    parser.add_argument("--output",type=Path,default=Path("data/derived/unite_reference_panel/annotations"))
    parser.add_argument("--cpu",type=int,default=6)
    parser.add_argument("--work",type=Path)
    parser.add_argument("--tool-manifest",type=Path,default=Path("tmp/unite-tools/tool_identity_manifest.json"))
    parser.add_argument("--run-manifest",type=Path,help="Explicit verified original profile-run identity for reusing frozen ITSx work")
    args=parser.parse_args();build(args.panel,args.itsx,args.output,args.cpu,args.work,args.tool_manifest,args.run_manifest)
