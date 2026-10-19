#!/usr/bin/env python3
"""Download the fixed Zenodo THP3 dataset using only Python's standard library."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'examples/thp3_20Mb/zenodo_downloads.json'
TOP = 'InDelMarker-offline-linux-x86_64-glibc2.28/'
PREFIXES = tuple('examples/thp3_20Mb/' + s + '/' for s in
                 ('data', 'data_payload', 'reference_results'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def matches(path, size, sha):
    return path.is_file() and path.stat().st_size == size and digest(path) == sha


def fetch(entry, directory):
    name = entry['name']
    if Path(name).name != name:
        raise ValueError('Invalid download name')
    target = directory / name
    if target.is_symlink():
        raise ValueError('Refusing symlink: ' + str(target))
    if matches(target, entry['bytes'], entry['sha256']):
        return name
    if target.exists():
        raise ValueError('Existing file differs; move it aside before retrying: ' + str(target))
    if not entry['url'].startswith('https://zenodo.org/api/records/23474945/files/'):
        raise ValueError('Unexpected source URL')
    for attempt in range(4):
        temporary = None
        try:
            request = urllib.request.Request(entry['url'], headers={
                'User-Agent': 'InDelMarker-THP3-downloader/1.0'})
            with tempfile.NamedTemporaryFile(dir=directory, prefix='.' + name,
                                             suffix='.tmp', delete=False) as out:
                temporary = Path(out.name)
                with urllib.request.urlopen(request, timeout=180) as response:
                    shutil.copyfileobj(response, out, 1024 * 1024)
            if not matches(temporary, entry['bytes'], entry['sha256']):
                raise ValueError('Downloaded size or SHA-256 mismatch: ' + name)
            os.replace(temporary, target)
            return name
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()


def assemble(manifest, directory):
    archive = directory / manifest['archive']
    if archive.is_symlink():
        raise ValueError('Refusing symlink: ' + str(archive))
    if matches(archive, manifest['archive_bytes'], manifest['archive_sha256']):
        return archive
    if archive.exists():
        raise ValueError('Existing archive differs; move it aside: ' + str(archive))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=directory, prefix='.assembled-',
                                         suffix='.tmp', delete=False) as out:
            temporary = Path(out.name)
            for index in range(1, 89):
                with (directory / (manifest['archive'] + '.part%03d' % index)).open('rb') as source:
                    shutil.copyfileobj(source, out, 4 * 1024 * 1024)
        if not matches(temporary, manifest['archive_bytes'], manifest['archive_sha256']):
            raise ValueError('Assembled archive size or SHA-256 mismatch')
        os.replace(temporary, archive)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return archive


def install_data(archive, directory, repository):
    """Copy only manifest-listed dataset files; preserve all code/configuration."""
    inventory = json.loads((directory / 'public_case_file_manifest.json').read_text())
    selected = {e['path']: e for e in inventory if e['path'].startswith(PREFIXES)}
    destinations = {}
    for name, entry in selected.items():
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Invalid archive path')
        dest = repository.joinpath(*relative.parts)
        if dest.is_symlink() or not dest.resolve().is_relative_to(repository.resolve()):
            raise ValueError('Unsafe destination: ' + str(dest))
        if dest.exists() and not matches(dest, entry['bytes'], entry['sha256']):
            raise ValueError('Preserve your modified file; move it aside first: ' + str(dest))
        destinations[name] = dest
    seen = set()
    with tarfile.open(archive, 'r|gz') as bundle:
        for member in bundle:
            name = member.name.removeprefix(TOP)
            if name not in selected:
                continue
            if not member.name.startswith(TOP) or not member.isfile() or name in seen:
                raise ValueError('Invalid or repeated dataset member: ' + member.name)
            seen.add(name)
            dest, entry = destinations[name], selected[name]
            if dest.exists():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=dest.parent, prefix='.' + dest.name,
                                                 suffix='.tmp', delete=False) as out:
                    temporary = Path(out.name)
                    with bundle.extractfile(member) as source:
                        shutil.copyfileobj(source, out, 4 * 1024 * 1024)
                if not matches(temporary, entry['bytes'], entry['sha256']):
                    raise ValueError('Dataset member checksum failed: ' + name)
                os.replace(temporary, dest)
            finally:
                if temporary is not None and temporary.exists():
                    temporary.unlink()
    if seen != set(selected):
        raise ValueError('Archive is missing listed dataset members')
    print('Installed/verified %d dataset files; code and configurations preserved.' % len(seen), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download-dir', type=Path, default=ROOT / '.downloads/thp3_20261019')
    parser.add_argument('--jobs', type=int, default=4, help='Concurrent public downloads (1–8)')
    parser.add_argument('--install-data', action='store_true', help='Install data and expected outputs into this checkout')
    args = parser.parse_args()
    if not 1 <= args.jobs <= 8:
        parser.error('--jobs must be between 1 and 8')
    manifest = json.loads(MANIFEST.read_text())
    directory = args.download_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    print('Dataset: ' + manifest['record_url'], flush=True)
    print('Allow about 5 GB free disk space for parts, archive and installed inputs.', flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(fetch, entry, directory) for entry in manifest['files']]
        for count, future in enumerate(as_completed(futures), 1):
            print('[%d/%d] verified %s' % (count, len(futures), future.result()), flush=True)
    archive = assemble(manifest, directory)
    print('Verified archive: ' + str(archive), flush=True)
    if args.install_data:
        install_data(archive, directory, ROOT)
    print('PASS. Completed files are reused on retry; no account or token is needed.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        sys.exit('ERROR: ' + str(error))
