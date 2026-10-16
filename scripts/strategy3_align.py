#!/usr/bin/env python3
"""Region-driven CDS localisation and sequence alignment for strategy 3."""
import argparse, os, re, subprocess
from concurrent.futures import ThreadPoolExecutor
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

p=argparse.ArgumentParser(); p.add_argument("--gff3",required=True); p.add_argument("--cds",required=True); p.add_argument("--genomes",required=True); p.add_argument("--chromosome"); p.add_argument("--start",type=int); p.add_argument("--end",type=int); p.add_argument("--regions-file"); p.add_argument("--max-genes",type=int); p.add_argument("--upstream",type=int,required=True); p.add_argument("--downstream",type=int,required=True); p.add_argument("--outdir",required=True); p.add_argument("--exonerate",required=True); p.add_argument("--minimap2"); p.add_argument("--bedtools",required=True); p.add_argument("--mafft",required=True); p.add_argument("--threads",type=int,default=1); a=p.parse_args()
if not a.regions_file and not (a.chromosome and a.start is not None and a.end is not None): p.error("provide --regions-file or --chromosome/--start/--end")
out=Path(a.outdir); (out/"genes").mkdir(parents=True,exist_ok=True); (out/"alignment").mkdir(exist_ok=True)
regions=[]
if a.regions_file:
    for line in open(a.regions_file):
        if line.strip() and not line.startswith("#"):
            chrom,start,end,*_ = line.split()
            regions.append((chrom,int(start),int(end)))
else: regions.append((a.chromosome,a.start,a.end))
genes=[]
for line in open(a.gff3):
    if line.startswith("#"): continue
    f=line.rstrip().split("\t")
    # Input can be conventional 9-column GFF3 or the supplied 6-column BED.
    if len(f) == 9:
        if f[2] not in {"gene", "mRNA", "transcript"}: continue
        start, end, strand = int(f[3]), int(f[4]), f[6]
        d=attrs(f[8]); gid=d.get("ID", d.get("Parent", "")).split(":")[-1]
    elif len(f) == 5:
        # B73v5_gff_db.txt: chromosome, 1-based start, end, gene ID, strand.
        start, end, gid, strand = int(f[1]), int(f[2]), f[3], f[4]
    elif len(f) >= 4:
        # BED start is 0-based and end-exclusive; report standard 1-based coordinates.
        start, end, strand, gid = int(f[1]) + 1, int(f[2]), f[5] if len(f) > 5 else ".", f[3]
    else:
        continue
    if any(f[0] == chrom and start <= region_end and end >= region_start for chrom,region_start,region_end in regions) and gid:
        genes.append((gid, f[0], str(start), str(end), strand))
seen=set(); genes=[x for x in genes if not (x[0] in seen or seen.add(x[0]))]
if a.max_genes is not None: genes=genes[:a.max_genes]
with open(out/"alignment/pos.tsv","w") as fh:
    fh.write("gene\tchromosome\tstart\tend\tstrand\n"); [fh.write("\t".join(x)+"\n") for x in genes]
cds=read_fasta(a.cds); genomes=dict(x.split("=",1) for x in a.genomes.split(","))
selected={}
for gid,*_ in genes:
    cds_id = gid if gid in cds else next((key for key in cds if key.startswith(gid + "_")), None)
    if cds_id: selected[gid]=cds[cds_id]

# Mapping all selected CDSs in one minimap2 invocation per genome is orders of
# magnitude faster than running exonerate against a multi-gigabase genome for
# every gene. The legacy exonerate path remains available when --minimap2 is
# not supplied.
if a.minimap2:
    query_all=out/"genes"/"selected.cds.fa"
    query_all.write_text("".join(f">{gid}\n{seq}\n" for gid,seq in selected.items()))
    hits={label:{} for label in genomes}
    for label,genome in genomes.items():
        sam=out/"genes"/f"{label}.cds.sam"
        with open(sam,"w") as handle: run([a.minimap2,"-ax","splice","--secondary=no","-t",str(a.threads),genome,str(query_all)],handle)
        for line in open(sam):
            if line.startswith("@"):
                continue
            f=line.rstrip().split("\t")
            if len(f)<11 or int(f[1]) & 4:
                continue
            gid,flag,chrom,pos,cigar,mapq=f[0],int(f[1]),f[2],int(f[3]),f[5],int(f[4])
            span=sum(int(n) for n,op in re.findall(r"(\d+)([MDN=X])",cigar))
            candidate=(mapq,chrom,pos,pos+span-1,"-" if flag & 16 else "+")
            if gid not in hits[label] or candidate[0]>hits[label][gid][0]: hits[label][gid]=candidate
    alignment_jobs=[]
    for gid in selected:
        fasta_files=[out/"genes"/f"{gid}.cds.fa"]
        fasta_files[0].write_text(f">{gid}\n{selected[gid]}\n")
        for label,genome in genomes.items():
            if gid not in hits[label]:
                continue
            _,chrom,start,end,strand=hits[label][gid]
            bed=out/"genes"/f"{gid}.{label}.bed"; bed.write_text(f"{chrom}\t{max(0,start-1-a.upstream)}\t{end+a.downstream}\t{gid}_{label}\t0\t{strand}\n")
            seq=out/"genes"/f"{gid}.{label}.fa"; run([a.bedtools,"getfasta","-s","-name","-fi",genome,"-bed",str(bed),"-fo",str(seq)]); fasta_files.append(seq)
        aln=out/"alignment"/f"{gid}.fa"
        merged=out/"genes"/f"{gid}.merged.fa"
        merged.write_text("".join(Path(item).read_text() for item in fasta_files))
        alignment_jobs.append((gid, merged, aln))

    # One long interval may contain dozens of genes. Running the independent
    # MAFFT jobs concurrently is substantially faster than serial alignment;
    # each job receives one thread so the total never exceeds --threads.
    def align(job):
        _, merged, aln = job
        with open(aln,"w") as fh:
            run([a.mafft,"--thread","1","--auto",str(merged)],fh)
    with ThreadPoolExecutor(max_workers=max(1,a.threads)) as pool:
        list(pool.map(align,alignment_jobs))
    with open(out/"alignment/manifest.tsv","w") as manifest:
        manifest.write("gene\talignment\n")
        for gid,_,aln in alignment_jobs:
            manifest.write(f"{gid}\t{aln}\n")
    raise SystemExit(0)
with open(out/"alignment/manifest.tsv","w") as manifest:
    manifest.write("gene\talignment\n")
    for gid,*_ in genes:
        query=out/"genes"/f"{gid}.cds.fa"
        if gid not in selected: continue
        query.write_text(f">{gid}\n{selected[gid]}\n")
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
        merged=out/"genes"/f"{gid}.merged.fa"
        merged.write_text("".join(Path(item).read_text() for item in fasta_files))
        with open(aln,"w") as fh: run([a.mafft,"--thread",str(a.threads),"--auto",str(merged)],fh)
        manifest.write(f"{gid}\t{aln}\n")
