#!/usr/bin/env python3
"""Content checks for the expanded THP3 case; accepts genuine missing query hits."""
import argparse,csv,hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile
ROOT=Path(__file__).resolve().parent

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()

def fasta(path):
 result={};key=None;chunks=[]
 with Path(path).open() as h:
  for line in h:
   if line.startswith('>'):
    if key is not None:result[key]=''.join(chunks)
    key=line[1:].strip();chunks=[]
   else:chunks.append(line.strip())
 if key is not None:result[key]=''.join(chunks)
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--strategy',choices=['1','2','3'],required=True);p.add_argument('--outdir',type=Path);a=p.parse_args()
 out=a.outdir or ROOT/'results'/f'strategy{a.strategy}'
 info=json.loads((ROOT/'provenance.json').read_text())
 assemblies={'data/B73.THP3_flanks20Mb.fa','data/Ames21814.THP3_flanks20Mb.fa'}
 required=assemblies if a.strategy=='1' else assemblies|{'data/B73.genes.bed','data/B73.longest_cds.fa','data/B73.original_gene_coordinates.tsv'} if a.strategy=='3' else {'data/B73.chr3.fa'}|{f'data/{sample}.THP3_flanks20Mb.bam{suffix}' for sample in ['B73','Ames21814'] for suffix in ['', '.bai']}
 for name in sorted(required):
  expected=info['sha256'][name]
  path=ROOT/name
  assert path.is_file() and sha(path)==expected,f'Input changed: {path}'
 report=dict(strategy=a.strategy,input_sha256='PASS',checked_input_files=sorted(required),scope='Computational outputs only; no PCR or full-genome primer-specificity validation')
 if a.strategy in {'1','2'}:
  offset=0;minimum,maximum=100,500
  if a.strategy=='1':
   ref=next(iter(fasta(ROOT/'data/B73.THP3_flanks20Mb.fa').values()))
   vcf=out/'strategy1/filtered_indels.vcf'
  else:
   ref=fasta(ROOT/'data/B73.chr3.fa')['chr3']
   minimum,maximum=20,100
   assert len(ref)==238017767
   import pysam
   read_info=json.loads((ROOT/'sources/regional_bam_provenance.json').read_text())
   assert read_info['status']=='PASS' and read_info['region']=='chr3:180402538-220408213'
   for label in ['B73','Ames21814']:
    metadata=read_info['samples'][label]
    assert metadata['source_region_alignment_records']==metadata['exported_alignment_records']
    assert metadata['source_depth']==metadata['exported_depth']
    with pysam.AlignmentFile(str(ROOT/'data'/f'{label}.THP3_flanks20Mb.bam'),'rb') as bam:
     assert bam.check_index() and bam.get_reference_length('chr3')==len(ref)
     count=bam.count('chr3',180402537,220408213,read_callback='nofilter')
     assert count==metadata['exported_alignment_records']['records']
     assert {rg['SM'] for rg in bam.header.to_dict()['RG']}=={'B73' if label=='B73' else 'Teosinte'}
   vcf=out/'strategy2/contrasting_indels.vcf'
  variants=[l.rstrip().split('\t') for l in vcf.open() if not l.startswith('#') and l.strip()]
  assert variants,'Header-only VCF is not a positive demonstration'
  for v in variants:
   assert minimum<=abs(len(v[3])-len(v[4]))<=maximum
   position=int(v[1])-1-offset
   assert 0<=position<len(ref) and ref[position:position+len(v[3])]==v[3]
   if a.strategy=='2':
    assert v[0]=='chr3' and 180402538<=int(v[1])<=220408213
    gt_index=v[8].split(':').index('GT')
    assert {sample.split(':')[gt_index].replace('|','/') for sample in v[9:11]}=={'0/0','1/1'}
  ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
  with ZipFile(out/'primer3/primer_markers.xlsx') as z:sheet=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
  rows=[[ ''.join(c.itertext()) for c in row.findall('s:c',ns)] for row in sheet.findall('.//s:row',ns)][1:]
  assert len(rows)==len(variants)
  pairs=0;rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
  for row in rows:
   template=row[5].replace('[','').replace(']','')
   for k in (6,15,24):
    if not row[k] or not row[k+3]:continue
    left,right=int(row[k+6]),int(row[k+7]);f,r=row[k],row[k+3]
    assert template[left:left+len(f)]==f
    assert template[right-len(r)+1:right+1]==rc(r)
    assert right-left+1==int(row[k+8]);pairs+=1
  assert pairs>0
  report.update(filtered_indel_loci=len(variants),primer_pairs_matching_reference_template=pairs)
 else:
  genes={l.split('\t')[3] for l in (ROOT/'data/B73.genes.bed').read_text().splitlines() if l}
  rows=list(csv.DictReader((out/'strategy3/alignment/manifest.tsv').open(),delimiter='\t'));assert {x['gene'] for x in rows}==genes
  both=0;missing=[];statuses=[]
  for row in rows:
   alignment=Path(row['alignment'])
   if not alignment.exists():alignment=out/'strategy3/alignment'/alignment.name
   seqs=fasta(alignment);lengths={len(s) for s in seqs.values()}
   assert len(lengths)==1 and next(iter(lengths))>0 and row['gene'] in seqs
   present=[label for label in ['B73','Ames21814'] if any(f'{row["gene"]}_{label}::' in name for name in seqs)]
   assert 'B73' in present,f'B73 sequence missing: {row["gene"]}'
   if len(present)==2:both+=1
   else:missing.append(row['gene'])
   if row['gene']=='Zm00001eb152450':assert len(present)==2,'THP3 must include both genomic sequences'
   statuses.append([row['gene'],','.join(present),','.join(x for x in ['B73','Ames21814'] if x not in present),next(iter(lengths))])
  with (out/'strategy3/alignment/genome_coverage.tsv').open('w') as h:
   w=csv.writer(h,delimiter='\t');w.writerow(['gene','genomes_in_alignment','genomes_without_extracted_sequence','alignment_columns']);w.writerows(statuses)
  report.update(aligned_genes=len(rows),both_genomes=both,without_Ames21814_sequence=len(missing),missing_query_genes=missing,missing_hit_interpretation='No extracted hit within this Ames21814 interval; not evidence of biological gene absence')
 report['status']='PASS';(out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='missing_query_genes'},indent=2))
if __name__=='__main__':main()
