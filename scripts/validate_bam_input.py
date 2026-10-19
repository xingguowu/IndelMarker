#!/usr/bin/env python3
"""Validate prealigned BAM metadata without altering or filtering reads."""
import argparse
import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference', required=True)
p.add_argument('--bam', required=True)
p.add_argument('--index', required=True)
a = p.parse_args()
with pysam.AlignmentFile(a.bam, 'rb', index_filename=a.index) as bam, pysam.FastaFile(a.reference) as ref:
    if not bam.check_index():
        raise ValueError('BAM index is missing or unreadable')
    header = bam.header.to_dict()
    if header.get('HD', {}).get('SO') != 'coordinate':
        raise ValueError('BAM must be coordinate sorted')
    samples = {rg.get('SM') for rg in header.get('RG', [])}
    if len(samples) != 1 or not next(iter(samples)):
        raise ValueError('Each parental BAM must have exactly one nonempty RG sample name')
    sequences = dict(zip(bam.references, bam.lengths))
    for name, length in zip(ref.references, ref.lengths):
        if sequences.get(name) != length:
            raise ValueError(f'Reference/BAM contig mismatch: {name}, FASTA={length}, BAM={sequences.get(name)}')
    print(f'Prealigned BAM metadata: PASS; sample={next(iter(samples))}; no realignment, filtering or duplicate marking applied')
