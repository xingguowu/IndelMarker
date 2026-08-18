#!/usr/bin/env python3
"""Convert retained VCF indels to Primer3 Boulder-IO input on reference sequence."""
import argparse

def fasta(path):
    data, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: data[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else: chunks.append(line.strip())
    if name: data[name] = "".join(chunks).upper()
    return data

p = argparse.ArgumentParser(); p.add_argument("--vcf", required=True); p.add_argument("--reference", required=True); p.add_argument("--output", required=True); p.add_argument("--flank", type=int, default=400); p.add_argument("--product-size-range", default="250-700"); a = p.parse_args()
seqs = fasta(a.reference)
with open(a.vcf) as src, open(a.output, "w") as dst:
    for line in src:
        if line.startswith("#"): continue
        chrom, pos, _, ref, alt, *_ = line.rstrip().split("\t")
        if chrom not in seqs or "," in alt: continue
        pos = int(pos); left = max(0, pos - 1 - a.flank); right = min(len(seqs[chrom]), pos - 1 + len(ref) + a.flank)
        template = seqs[chrom][left:right]; target_start = pos - 1 - left
        dst.write(f"SEQUENCE_ID={chrom}_{pos}\nSEQUENCE_TEMPLATE={template}\nSEQUENCE_TARGET={target_start},{len(ref)}\nPRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\nPRIMER_NUM_RETURN=3\nPRIMER_PRODUCT_SIZE_RANGE={a.product_size_range}\n=\n")
