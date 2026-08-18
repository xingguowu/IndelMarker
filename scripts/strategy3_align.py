#!/usr/bin/env python3
"""Region-driven CDS localisation and sequence alignment for strategy 3."""
import argparse, os, re, subprocess
from pathlib import Path

def run(cmd, stdout=None): subprocess.run(cmd, check=True, stdout=stdout)
def attrs(s): return dict(x.split("=", 1) for x in s.split(";") if "=" in x)
def read_fasta(path):
    d, k, b = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if k: d[k] = "".join(b)
            k, b = line[1:].split()[0], []
        else: b.append(line.strip())
    if k: d[k] = "".join(b)
    return d

p=argparse.ArgumentParser(); p.add_argument("--gff3",required=True); p.add_argument("--cds",required=True); p.add_argument("--genomes",required=True); p.add_argument("--chromosome",required=True); p.add_argument("--start",type=int,required=True); p.add_argument("--end",type=int,required=True); p.add_argument("--upstream",type=int,required=True); p.add_argument("--downstream",type=int,required=True); p.add_argument("--outdir",required=True); p.add_argument("--exonerate",required=True); p.add_argument("--bedtools",required=True); p.add_argument("--mafft",required=True); p.add_argument("--threads",type=int,default=1); a=p.parse_args()
out=Path(a.outdir); (out/"genes").mkdir(parents=True,exist_ok=True); (out/"alignment").mkdir(exist_ok=True)
genes=[]
for line in open(a.gff3):
    if line.startswith("#"): continue
    f=line.rstrip().split("\t")
    # Input can be conventional 9-column GFF3 or the supplied 6-column BED.
    if len(f) == 9:
        if f[0]!=a.chromosome or f[2] not in {"gene", "mRNA", "transcript"}: continue
        start, end, strand = int(f[3]), int(f[4]), f[6]
        d=attrs(f[8]); gid=d.get("ID", d.get("Parent", "")).split(":")[-1]
    elif len(f) >= 4:
        if f[0] != a.chromosome: continue
        # BED start is 0-based and end-exclusive; report standard 1-based coordinates.
        start, end, strand, gid = int(f[1]) + 1, int(f[2]), f[5] if len(f) > 5 else ".", f[3]
    else:
        continue
    if start <= a.end and end >= a.start and gid:
        genes.append((gid, f[0], str(start), str(end), strand))
seen=set(); genes=[x for x in genes if not (x[0] in seen or seen.add(x[0]))]
with open(out/"alignment/pos.tsv","w") as fh:
    fh.write("gene\tchromosome\tstart\tend\tstrand\n"); [fh.write("\t".join(x)+"\n") for x in genes]
cds=read_fasta(a.cds); genomes=dict(x.split("=",1) for x in a.genomes.split(","))
with open(out/"alignment/manifest.tsv","w") as manifest:
    manifest.write("gene\talignment\n")
    for gid,*_ in genes:
        query=out/"genes"/f"{gid}.cds.fa"
        if gid not in cds: continue
        query.write_text(f">{gid}\n{cds[gid]}\n")
        fasta_files=[query]
        for label, genome in genomes.items():
            exo=out/"genes"/f"{gid}.{label}.exonerate"
            with open(exo,"w") as log: run([a.exonerate,str(query),genome,"--model","est2genome","--showtargetgff","yes","--showvulgar","no","--showalignment","no","--bestn","1","--cores",str(a.threads)],log)
            target_name=None; target_range=None
            for line in open(exo):
                m=re.search(r"Target:\s+(\S+)", line)
                if m: target_name=m.group(1)
                m=re.search(r"Target range:\s+(\d+)\s+->\s+(\d+)",line)
                if m: target_range=(int(m.group(1)),int(m.group(2)))
                if target_name and target_range: break
            if not (target_name and target_range): continue
            chrom=target_name; s,e=target_range; lo=max(0,min(s,e)-1-a.upstream); hi=max(s,e)+a.downstream
            bed=out/"genes"/f"{gid}.{label}.bed"; bed.write_text(f"{chrom}\t{lo}\t{hi}\t{gid}_{label}\t0\t{'+' if s<=e else '-'}\n")
            seq=out/"genes"/f"{gid}.{label}.fa"
            run([a.bedtools,"getfasta","-s","-name","-fi",genome,"-bed",str(bed),"-fo",str(seq)]); fasta_files.append(seq)
        aln=out/"alignment"/f"{gid}.fa"
        with open(aln,"w") as fh: run([a.mafft,"--auto",*[str(x) for x in fasta_files]],fh)
        manifest.write(f"{gid}\t{aln}\n")
