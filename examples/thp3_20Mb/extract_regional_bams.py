#!/usr/bin/env python3
"""Export all alignments overlapping a region, without filtering or subsampling.

Run where the original, indexed BAMs are available. Only command-line/path
metadata in headers are removed; all alignment records and coordinates remain
unchanged. Original files are opened read-only. No credentials are recorded.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def run(args, **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def digest_command(args, depth=False):
    h = hashlib.sha256()
    count = total = covered = 0
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    for line in p.stdout:
        h.update(line)
        count += 1
        if depth:
            value = int(line.rstrip().split(b'\t')[-1])
            total += value
            covered += value > 0
    if p.wait():
        raise RuntimeError('Verification command failed')
    result = {'sha256': h.hexdigest(), 'records': count}
    if depth:
        result.update(depth_sum=total, covered_positions=covered,
                      mean_depth=total / count if count else 0)
    return result


def file_sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--b73', type=Path, required=True)
    p.add_argument('--ames21814', type=Path, required=True)
    p.add_argument('--region', default='chr3:180402538-220408213')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--samtools', default='samtools')
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=False)
    version = subprocess.check_output([a.samtools, '--version'], text=True).splitlines()[0]
    report = {'region': a.region, 'coordinates': 'B73 NAM-5.0, 1-based inclusive',
              'samtools': version, 'extraction': 'samtools view -b --no-PG SOURCE REGION',
              'subsampling': False, 'mapping_quality_filter': None,
              'alignment_flag_filter': None, 'duplicate_removal': False,
              'header_changes': 'Remove @PG CL, @SQ UR and @CO only; retain sequence dictionary, RG and alignment records',
              'depth_command': 'samtools depth -a -d 0 -q 0 -Q 0 -g 0x700 -r REGION BAM',
              'depth_definition': 'Per-position aligned-base depth, including duplicate, secondary and QC-failed records; no MAPQ/base-quality threshold; deletion positions not counted',
              'samples': {}}
    for label, source in [('B73', a.b73), ('Ames21814', a.ames21814)]:
        print('Extracting ' + label, flush=True)
        dest = a.outdir / (label + '.THP3_flanks20Mb.bam')
        raw = a.outdir / (label + '.unsanitized.tmp.bam')
        header = subprocess.check_output([a.samtools, 'view', '--no-PG', '-H', str(source)], text=True)
        clean = []
        for line in header.splitlines():
            fields = line.split('\t')
            if fields[0] == '@CO':
                continue
            remove = 'CL:' if fields[0] == '@PG' else 'UR:' if fields[0] == '@SQ' else None
            clean.append('\t'.join(x for x in fields if not (remove and x.startswith(remove))))
        header_file = a.outdir / (label + '.header.sam')
        header_file.write_text('\n'.join(clean) + '\n')
        run([a.samtools, 'view', '--no-PG', '-@', '2', '-b', '-o', str(raw), str(source), a.region])
        with dest.open('wb') as out:
            run([a.samtools, 'reheader', '-P', str(header_file), str(raw)], stdout=out)
        raw.unlink()
        header_file.unlink()
        run([a.samtools, 'index', '-@', '2', str(dest)])
        run([a.samtools, 'quickcheck', '-v', str(dest)])
        sam_args = [a.samtools, 'view', '--no-PG']
        before = digest_command(sam_args + [str(source), a.region])
        after = digest_command(sam_args + [str(dest), a.region])
        assert before == after, 'Alignment records changed for ' + label
        depth_args = [a.samtools, 'depth', '-a', '-d', '0', '-q', '0', '-Q', '0', '-g', '0x700', '-r', a.region]
        depth_before = digest_command(depth_args + [str(source)], depth=True)
        depth_after = digest_command(depth_args + [str(dest)], depth=True)
        assert depth_before == depth_after, 'Depth changed for ' + label
        entry = {'material': label, 'source_bam_basename': source.name,
                 'source_bytes': source.stat().st_size,
                 'source_region_alignment_records': before, 'exported_alignment_records': after,
                 'source_depth': depth_before, 'exported_depth': depth_after,
                 'files': [{'name': f.name, 'bytes': f.stat().st_size, 'sha256': file_sha(f)}
                           for f in [dest, Path(str(dest) + '.bai')]],
                 'status': 'PASS: all regional alignment records and per-position depth identical'}
        report['samples'][label] = entry
        (a.outdir / 'regional_bam_provenance.json').write_text(json.dumps(report, indent=2) + '\n')
        print(label + ': ' + entry['status'], flush=True)
    report['status'] = 'PASS'
    (a.outdir / 'regional_bam_provenance.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
