#!/usr/bin/env python3
"""Extract profile-supported ITS regions without forcing a universal MSA."""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from scan_unite_primer_sites import read_fasta, read_csv, reverse_complement, write_gzip_csv

REGIONS = ("SSU", "ITS1", "5.8S", "ITS2", "LSU")

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_gzip_text(path, lines):
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as compressed:
            with io.TextIOWrapper(compressed,encoding="utf-8") as handle:
                for line in lines: handle.write(line)

def parse_positions(path):
    rows = {}
    for line in path.read_text().splitlines():
        fields=line.split("\t")
        if not fields: continue
        record_id=fields[0].split()[0]
        if record_id in rows: raise ValueError("ITSx position identity repeats")
        intervals={}
        for field in fields[2:]:
            match=re.fullmatch(r"(SSU|ITS1|5\.8S|ITS2|LSU): (\d+)-(\d+)",field.strip())
            if match: intervals[match[1]]=(int(match[2]),int(match[3]))
        rows[record_id]=intervals
    return rows

def parse_strands(path):
    rows={}
    for line in path.read_text().splitlines():
        fields=line.split("\t")
        if len(fields)<4: raise ValueError("ITSx detailed results schema changed")
        record_id=fields[0].split()[0]
        if fields[3] not in ("0","1"): raise ValueError("ITSx strand flag is not 0/1")
        if record_id in rows: raise ValueError("ITSx strand identity repeats")
        rows[record_id] = -1 if fields[3] == "1" else 1
    return rows

def boundary_status(region,intervals,length):
    if region not in intervals: return "not_detected"
    start,end=intervals[region]
    if not (1<=start<=end<=length): return "invalid_coordinates"
    if region in ("SSU","LSU"): return "flank_detected"
    if region == "ITS1": complete="SSU" in intervals and "5.8S" in intervals
    elif region == "ITS2": complete="5.8S" in intervals and "LSU" in intervals
    else: complete=start>1 and end<length and "ITS1" in intervals and "ITS2" in intervals
    return "complete" if complete else "partial"

def build(panel,itsx,output,cpus=6,existing_work=None):
    source=panel/"selected_raw.fasta.gz"
    metadata=read_csv(panel/"selected_metadata.csv.gz")
    sequences=dict(read_fasta(source))
    if set(sequences)!={r["record_id"] for r in metadata}: raise ValueError("UNITE metadata/sequence membership differs")
    output.mkdir(parents=True,exist_ok=True)
    previous=output/"annotation_manifest.json"
    if previous.exists():
        state=json.loads(previous.read_text())
        if state["input_sha256"] != sha(source) or state["script_sha256"]!=sha(__file__):
            raise ValueError("Annotation identity differs; choose a new output directory")
        for artifact in state["artifacts"].values():
            if sha(output/artifact["path"])!=artifact["sha256"]: raise ValueError("Annotation artifact checksum differs")
        return state
    work=Path(existing_work) if existing_work else Path(tempfile.mkdtemp(prefix="primer-atlas-unite-annotation-"))
    if any(c.isspace() for c in str(work)): raise ValueError("ITSx work directory must not contain spaces")
    work.mkdir(parents=True,exist_ok=True)
    with (work/"panel.fasta").open("w") as handle:
        for record_id,seq in sequences.items(): handle.write(f">{record_id}\n{seq}\n")
    command=[str(itsx.resolve()),"-i","panel.fasta","-o","regions","-t","F","--preserve","T",
             "--save_regions","all","--detailed_results","T","--save_raw","T","--cpu",str(cpus),
             "--temp",str(work/"temp")]
    (work/"temp").mkdir(exist_ok=True)
    if not (work/"regions.positions.txt").exists():
        print(f"Annotating {len(sequences)} UNITE records with ITSx in {work}",flush=True)
        with (work/"ITSx.log").open("w") as log:
            subprocess.run(command,cwd=work,stdout=log,stderr=subprocess.STDOUT,check=True)
    required=[work/("regions"+suffix) for suffix in (".positions.txt",".extraction.results",".summary.txt",".ITS1.fasta",".5_8S.fasta",".ITS2.fasta")]
    if any(not p.exists() for p in required): raise ValueError("ITSx did not produce required artifacts; inspect "+str(work/"ITSx.log"))
    positions=parse_positions(work/"regions.positions.txt")
    strands=parse_strands(work/"regions.extraction.results")
    if set(positions)!=set(strands) or not set(positions)<=set(sequences): raise ValueError("ITSx output ownership differs")
    extracted={region:dict(read_fasta(work/("regions."+suffix+".fasta"))) for region,suffix in (("ITS1","ITS1"),("5.8S","5_8S"),("ITS2","ITS2"))}
    oriented={}; region_rows=[]; orientation_rows=[]; counts=Counter(); mismatch=[]
    for record_id,seq in sequences.items():
        strand=strands.get(record_id,1)
        analysis=reverse_complement(seq) if strand==-1 else seq
        oriented[record_id]=analysis
        orientation_rows.append(dict(record_id=record_id,source_strand=strand,
            orientation="profile_forward" if record_id in strands else "source_unverified",
            source_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),analysis_sequence_sha256=hashlib.sha256(analysis.encode()).hexdigest()))
        bounds=positions.get(record_id,{})
        for region in REGIONS:
            status=boundary_status(region,bounds,len(seq)); interval=bounds.get(region,(None,None))
            if region in extracted and status in ("complete","partial"):
                observed=extracted[region].get(record_id)
                expected=analysis[interval[0]-1:interval[1]]
                if observed != expected:
                    status="extraction_mismatch"
                    mismatch.append(dict(record_id=record_id,region=region,reason="ITSx extraction differs from the oriented raw slice"))
            counts[region+":"+status]+=1
            region_rows.append(dict(record_id=record_id,region=region,start=interval[0],end=interval[1],
                boundary_status=status,source_strand=strand,annotation_source="ITSx1.1.3 fungal HMM profiles",
                coordinate_basis="1-based inclusive positions in oriented source sequence",
                reason="Both bounding rRNA profiles detected" if status=="complete" and region in ("ITS1","ITS2") else status))
    write_gzip_text(output/"oriented_raw.fasta.gz",(f">{key}\n{value}\n" for key,value in oriented.items()))
    write_gzip_csv(output/"regions.csv.gz",list(region_rows[0]),region_rows)
    write_gzip_csv(output/"orientations.csv.gz",list(orientation_rows[0]),orientation_rows)
    write_gzip_csv(output/"annotation_issues.csv.gz",["record_id","region","reason"],mismatch)
    for region in ("ITS1","5.8S","ITS2"):
        selected=[r for r in region_rows if r["region"]==region and r["boundary_status"]=="complete"]
        write_gzip_text(output/(region.replace(".","_")+".fasta.gz"),
            (f'>{r["record_id"]}\n{oriented[r["record_id"]][r["start"]-1:r["end"]]}\n' for r in selected))
    for p in required:
        write_gzip_text(output/(p.name+".gz"),iter(p.read_text().splitlines(keepends=True)))
    files=[p for p in output.iterdir() if p.is_file()]
    state=dict(input_sha256=sha(source),script_sha256=sha(__file__),n_input=len(sequences),n_detected=len(positions),
        n_reverse=sum(v==-1 for v in strands.values()),region_counts=dict(sorted(counts.items())),command=command,work_directory=str(work),
        limitation="ITSx computational boundaries; SSU/LSU intervals are retained flanks, not claims of complete genes. Only both-flank-supported ITS1/ITS2 regions enter regional alignments.",
        artifacts={p.name:dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in files})
    previous.write_text(json.dumps(state,indent=2)+"\n")
    print(json.dumps({k:state[k] for k in ("n_input","n_detected","n_reverse","region_counts")},indent=2),flush=True)
    return state

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--panel",type=Path,default=Path("data/derived/unite_reference_panel"))
    parser.add_argument("--itsx",type=Path,default=Path("tmp/unite-tools/bin/ITSx"))
    parser.add_argument("--output",type=Path,default=Path("data/derived/unite_reference_panel/annotations"))
    parser.add_argument("--cpu",type=int,default=6)
    parser.add_argument("--work",type=Path)
    args=parser.parse_args();build(args.panel,args.itsx,args.output,args.cpu,args.work)
