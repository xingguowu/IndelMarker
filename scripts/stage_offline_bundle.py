#!/usr/bin/env python3
"""Preserve application files, exact package records, and upstream notices."""
import json
import platform
import shutil
import sys
from pathlib import Path

source, runtime, stage = map(lambda x: Path(x).resolve(), sys.argv[1:])
stage.mkdir()
manifest = source / 'envs/offline-source-files.txt'
for name in manifest.read_text().splitlines():
    relative = Path(name)
    path = (source / relative).resolve()
    if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(source) or not path.is_file():
        raise ValueError(f'Invalid or missing application source: {name}')
    destination = stage / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
(stage / 'offline/linux-64').mkdir(parents=True)
(stage / 'software').mkdir()
records = []
for metadata in sorted((runtime / 'conda-meta').glob('*.json')):
    package = json.loads(metadata.read_text())
    records.append({key: package.get(key) for key in ('name', 'version', 'build', 'subdir', 'url', 'sha256', 'md5', 'license', 'depends')})
    cached = Path((package.get('link') or {}).get('source', ''))
    package_dir = stage / 'software/upstream-notices' / metadata.stem
    for name in ('licenses', 'recipe'):
        origin = cached / 'info' / name
        if origin.is_dir():
            shutil.copytree(origin, package_dir / name)
(stage / 'software/packages.json').write_text(json.dumps(records, indent=2) + '\n')
(stage / 'software/platform.json').write_text(json.dumps({
    'os': 'Linux', 'architecture': 'x86_64', 'minimum_glibc': '2.28',
    'builder': platform.platform(), 'package_count': len(records),
    'runtime': 'offline/linux-64/indelmarker-runtime.tar.gz',
    'source': 'Official conda-forge and Bioconda package channels',
    'offline_validation': 'Must pass a fresh-container --network none test before release.'
}, indent=2) + '\n')
shutil.copy2(stage / 'docs/OFFLINE_INSTALL.md', stage / 'OFFLINE_README.md')
print(f'Staged application and {len(records)} exact package records.')
