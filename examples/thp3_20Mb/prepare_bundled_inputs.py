#!/usr/bin/env python3
"""Restore the bundled large inputs without network access; check every SHA-256."""
import gzip,hashlib,json,os,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()

def main():
 manifest=json.loads((ROOT/'bundled_inputs.json').read_text())
 restored=0
 for entry in manifest['files']:
  dest=ROOT/entry['destination']
  if dest.exists():
   if dest.stat().st_size!=entry['bytes'] or sha(dest)!=entry['sha256']:raise ValueError(f'Existing input differs from the bundled version; preserve your file and restore from a fresh download: {dest.name}')
   continue
  for item in entry['payloads']:
   src=ROOT/item['path']
   if src.stat().st_size!=item['bytes'] or sha(src)!=item['sha256']:raise ValueError(f'Corrupted bundled payload: {src.name}')
  dest.parent.mkdir(parents=True,exist_ok=True)
  temporary=None
  try:
   with tempfile.NamedTemporaryFile(prefix=dest.name+'.',suffix='.tmp',dir=dest.parent,delete=False) as out:
    temporary=Path(out.name)
    for item in entry['payloads']:
     src=ROOT/item['path'];op=gzip.open if entry['encoding']=='gzip' else open
     with op(src,'rb') as inp:
      for block in iter(lambda:inp.read(4*1024*1024),b''):out.write(block)
   if temporary.stat().st_size!=entry['bytes'] or sha(temporary)!=entry['sha256']:raise ValueError(f'Reconstructed input checksum failed: {dest.name}')
   os.replace(temporary,dest);restored+=1
  finally:
   if temporary is not None and temporary.exists():temporary.unlink()
 print(f'Bundled input preparation: PASS ({restored} restored; all {len(manifest["files"])} checked).')
if __name__=='__main__':main()
