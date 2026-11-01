#!/usr/bin/env python3
"""Build bounded regional comparison MSAs within frozen fungal clades."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from scan_unite_primer_sites import read_csv, read_fasta, write_gzip_csv
from annotate_unite_regions import write_gzip_text

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def exemplars(rows,limit):
    # Already sampled per order. Preserve several orders/genera in a class,
    # not just the first abundant lineage; every chosen identity stays visible.
    rows=sorted(rows,key=lambda r:(float(r.get("ambiguity_fraction",0)),r["representative_type"]!="refs",r["record_id"]))
    selected=[]; used=set()
    for column in ("order_group_id","genus",None):
        seen=set()
        for row in rows:
            if row["record_id"] in used: continue
            if column and row[column] in seen: continue
            selected.append(row);used.add(row["record_id"])
            if column: seen.add(row[column])
            if len(selected)>=limit:return selected
    return selected

def build(panel,annotations,output,limit=32,workers=4,mafft="mafft"):
    output.mkdir(parents=True,exist_ok=True)
    metadata=read_csv(panel/"selected_metadata.csv.gz")
    by_id={r["record_id"]:r for r in metadata}
    regions=read_csv(annotations/"regions.csv.gz")
    sequences=dict(read_fasta(annotations/"oriented_raw.fasta.gz"))
    if set(sequences)!=set(by_id):raise ValueError("Regional input membership differs")
    completed=[r for r in regions if r["region"] in ("ITS1","5.8S","ITS2") and r["boundary_status"]=="complete"]
    groups={}
    for row in completed:
        meta=by_id[row["record_id"]]
        for rank in ("class","order"):
            status_column={"class":"class_name_status","order":"order_name_status"}[rank]
            if meta[status_column]!="resolved":continue
            key=(rank,meta[rank+"_group_id"],row["region"])
            groups.setdefault(key,[]).append(dict(meta,region_start=int(row["start"]),region_end=int(row["end"])))
    command_identity=dict(input_sha256=sha(panel/"selected_metadata.csv.gz"),regions_sha256=sha(annotations/"regions.csv.gz"),
        oriented_sha256=sha(annotations/"oriented_raw.fasta.gz"),script_sha256=sha(__file__),limit=limit,mafft_version=subprocess.run([mafft,"--version"],capture_output=True,text=True,check=True).stderr.strip())
    manifest_path=output/"alignment_manifest.json"
    if manifest_path.exists():
        old=json.loads(manifest_path.read_text())
        if old["identity"]!=command_identity:raise ValueError("Regional alignment identity differs; use a new output")
        for item in old["artifacts"].values():
            if sha(output/item["path"])!=item["sha256"]:raise ValueError("Regional alignment checksum differs")
        return old
    (output/"alignments").mkdir(exist_ok=True)
    def align(item):
        (rank,group_id,region),rows=item
        chosen=exemplars(rows,limit)
        input_records={r["record_id"]:sequences[r["record_id"]][r["region_start"]-1:r["region_end"]] for r in chosen}
        token=hashlib.sha256((rank+"|"+group_id+"|"+region).encode()).hexdigest()[:20]
        relative=Path("alignments")/(rank+"_"+token+"_"+region.replace(".","_")+".fasta.gz")
        command=[mafft,"--auto","--thread","1"]
        with tempfile.TemporaryDirectory(prefix="unite-mafft-") as temporary:
            source=Path(temporary)/"input.fasta"
            source.write_text("".join(f">{key}\n{value}\n" for key,value in input_records.items()))
            if len(chosen)==1:
                alignment=input_records;log="Single retained regional record; no multi-row alignment inferred."
            else:
                run=subprocess.run(command+[str(source)],capture_output=True,text=True,check=True)
                result=Path(temporary)/"aligned.fasta";result.write_text(run.stdout)
                alignment=dict(read_fasta(result));log=run.stderr
        if set(alignment)!=set(input_records) or len({len(v) for v in alignment.values()})!=1:
            raise ValueError("MAFFT regional membership or widths differ")
        for key,value in alignment.items():
            if value.replace("-","")!=input_records[key]:raise ValueError("MAFFT changed a regional source base")
        write_gzip_text(output/relative,(f">{r['record_id']}\n{alignment[r['record_id']]}\n" for r in chosen))
        write_gzip_text(output/relative.with_suffix(".log.gz"),[log])
        width=len(next(iter(alignment.values())))
        return dict(rank=rank,group_id=group_id,region=region,file=str(relative),n_source_regions=len(rows),n_input=len(chosen),
            n_columns=width,n_gap_cells=sum(v.count("-") for v in alignment.values()),
            n_ambiguous_cells=sum(sum(c not in "ACGT-" for c in v) for v in alignment.values()),
            alignment_scope="Complete profile-bounded region in one resolved "+rank+"; bounded diverse exemplars",
            command=" ".join(command),sequence_retention="all extracted source bases retained")
    items=sorted(groups.items())
    print(f"Building {len(items)} regional class/order alignments; up to {limit} exemplars each",flush=True)
    index=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for number,row in enumerate(pool.map(align,items),1):
            index.append(row)
            if number%50==0:print(f"Verified {number}/{len(items)} regional alignments",flush=True)
    if not index:raise ValueError("No complete profile-bounded regions available")
    write_gzip_csv(output/"alignment_index.csv.gz",list(index[0]),index)
    files=[p for p in output.rglob("*") if p.is_file()]
    manifest=dict(identity=command_identity,n_alignments=len(index),n_alignment_records=sum(r["n_input"] for r in index),
        limitation="MSAs describe computational regional comparisons within clades, not universal positional homology across fungi. The full retained panel supplies primer-site denominators; each display alignment uses at most32 exemplars.",
        artifacts={str(p.relative_to(output)):dict(path=str(p.relative_to(output)),sha256=sha(p),bytes=p.stat().st_size) for p in files})
    manifest_path.write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--panel",type=Path,default=Path("data/derived/unite_reference_panel"))
    parser.add_argument("--annotations",type=Path,default=Path("data/derived/unite_reference_panel/annotations_v2"))
    parser.add_argument("--output",type=Path,default=Path("data/derived/unite_reference_panel/regional_alignments"))
    parser.add_argument("--limit",type=int,default=32)
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--mafft",default="mafft")
    args=parser.parse_args();build(args.panel,args.annotations,args.output,args.limit,args.workers,args.mafft)
