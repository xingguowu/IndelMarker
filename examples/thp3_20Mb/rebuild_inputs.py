#!/usr/bin/env python3
"""Re-extract THP3 +/-20 Mb inputs from the named public assemblies and B73 GFF3.
Both FASTA inputs may contain a complete genome or only chr3; gzip is supported.
Coordinates are 1-based inclusive. The packaged FASTA ID is chr3_THP3.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
REGIONS={'B73':(180402538,220408213),'Ames21814':(187896121,227901230)}

def op(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def extract(path,start,end):
 chunks=[];position=0;active=False
 with op(path) as f:
  for line in f:
   if line.startswith('>'):
    if active:break
    active=line[1:].split()[0]=='chr3';continue
   if not active:continue
   x=line.strip();n=position+len(x)
   if n>=start and position<end:chunks.append(x[max(0,start-1-position):min(len(x),end-position)])
   position=n
   if position>=end:break
 seq=''.join(chunks).upper()
 if len(seq)!=end-start+1:raise ValueError(f'Cannot extract chr3:{start}-{end} from {path}')
 return seq

def attrs(s):return dict(x.split('=',1) for x in s.split(';') if '=' in x)
def fasta(path,name,seq):path.write_text('>'+name+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--b73-genome',type=Path,required=True);p.add_argument('--ames21814-genome',type=Path,required=True);p.add_argument('--b73-gff3',type=Path,required=True);p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'data');a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=True);seqs={}
 for label,path in [('B73',a.b73_genome),('Ames21814',a.ames21814_genome)]:
  seqs[label]=extract(path,*REGIONS[label]);fasta(a.output/f'{label}.THP3_flanks20Mb.fa','chr3_THP3',seqs[label])
 # Prealigned BAMs retain original chromosome coordinates and lengths.
 # Reproduce samtools faidx's 80-column wrapping for the complete chr3 record.
 fasta(a.output/'B73.chr3.fa','chr3',extract(a.b73_genome,1,238017767))
 s,e=REGIONS['B73'];genes={};transcripts={};coding={};thp3=None
 with op(a.b73_gff3) as f:
  for line in f:
   if line.startswith('#'):continue
   f=line.rstrip().split('\t')
   if len(f)!=9 or f[0]!='chr3':continue
   meta=attrs(f[8])
   if f[2]=='gene' and meta.get('ID')=='Zm00001eb152450':thp3=(int(f[3]),int(f[4]))
   if int(f[3])<s+200 or int(f[4])>e-200:continue
   if f[2]=='gene':genes[meta['ID']]=(int(f[3]),int(f[4]),f[6])
   elif f[2] in {'mRNA','transcript'}:transcripts[meta['ID']]=meta.get('Parent','')
   elif f[2]=='CDS':
    for parent in meta.get('Parent','').split(','):coding.setdefault(parent,[]).append((int(f[3]),int(f[4]),f[6]))
 if thp3!=(200402538,200408213):raise ValueError(f'THP3 coordinates differ from the verified annotation: {thp3}')
 longest={};chosen={};complement=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
 for transcript,pieces in sorted(coding.items()):
  gene=transcripts.get(transcript)
  if gene not in genes:continue
  bases=''.join(seqs['B73'][lo-s:hi-s+1] for lo,hi,_ in sorted(pieces))
  if pieces[0][2]=='-':bases=bases.translate(complement)[::-1]
  if len(bases)>len(longest.get(gene,'')):longest[gene]=bases;chosen[gene]=transcript
 order=sorted(longest,key=lambda g:genes[g][0])
 (a.output/'B73.genes.bed').write_text(''.join(f'chr3_THP3\t{genes[g][0]-s}\t{genes[g][1]-s+1}\t{g}\t0\t{genes[g][2]}\n' for g in order))
 (a.output/'B73.longest_cds.fa').write_text(''.join(f'>{g}\n{longest[g]}\n' for g in order))
 (a.output/'B73.original_gene_coordinates.tsv').write_text('gene\tchromosome\tstart\tend\tstrand\tselected_transcript\tcds_bases\n'+''.join(f'{g}\tchr3\t{genes[g][0]}\t{genes[g][1]}\t{genes[g][2]}\t{chosen[g]}\t{len(longest[g])}\n' for g in order))
 print(json.dumps(dict(genes=len(order),thp3=thp3,sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(a.output.iterdir()) if p.suffix in {'.fa','.bed','.tsv'}}),indent=2))
if __name__=='__main__':main()
