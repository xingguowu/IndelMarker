#!/usr/bin/env python3
"""Fetch the complete B73 v5 reference for WGS mapping (not the 40-Mb extract)."""
import gzip,hashlib,json,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent
URL='https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz?download=thp3_whole_reference'

def main():
 data=ROOT/'data';data.mkdir(exist_ok=True);out=data/'B73.whole.fa'
 if out.exists():raise SystemExit(f'{out} exists; preserve it and verify its provenance before reusing.')
 archive=data/'B73.whole.fa.gz.download';partial=archive.with_suffix('.part')
 offset=partial.stat().st_size if partial.exists() else 0
 req=urllib.request.Request(URL,headers={'Range':f'bytes={offset}-','User-Agent':'CornMarker-public-example/1.0'})
 with urllib.request.urlopen(req,timeout=120) as response:
  if offset and (response.status!=206 or not response.headers.get('Content-Range','').startswith(f'bytes {offset}-')):raise RuntimeError('Resume range was not honored; partial file retained.')
  if not offset and response.status==206 and not response.headers.get('Content-Range','').startswith('bytes 0-'):raise RuntimeError('Server returned an unexpected cached byte range.')
  with partial.open('ab' if offset else 'wb') as f:
   while True:
    block=response.read(4*1024*1024)
    if not block:break
    f.write(block)
 lengths={};key=None;h=hashlib.sha256();temporary=out.with_suffix('.fa.part')
 with gzip.open(partial,'rb') as src,temporary.open('wb') as dst:
  for line in src:
   dst.write(line);h.update(line)
   if line.startswith(b'>'):key=line[1:].split()[0].decode();lengths[key]=0
   else:lengths[key]+=len(line.strip())
 if lengths.get('chr3')!=238017767:raise RuntimeError('Reference chr3 length differs from B73 NAM-5.0; temporary output retained.')
 temporary.rename(out);partial.rename(archive)
 info=dict(assembly='Zm-B73-REFERENCE-NAM-5.0',source=URL,chromosome_lengths=lengths,uncompressed_file_sha256=h.hexdigest(),gzip_validation='Full decompression including gzip CRC checks completed')
 out.with_suffix('.fa.provenance.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
if __name__=='__main__':main()
