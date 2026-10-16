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

p = argparse.ArgumentParser()
p.add_argument("--vcf", required=True); p.add_argument("--reference", required=True); p.add_argument("--output", required=True)
p.add_argument("--flank", type=int, default=400); p.add_argument("--product-size-range", default="250-700")
p.add_argument("--num-return", type=int, default=3)
p.add_argument("--optimal-primer-size", type=int, default=24); p.add_argument("--min-primer-size", type=int, default=18); p.add_argument("--max-primer-size", type=int, default=26)
p.add_argument("--optimal-tm", type=float, default=60); p.add_argument("--min-tm", type=float, default=57); p.add_argument("--max-tm", type=float, default=60)
p.add_argument("--min-gc-percent", type=float, default=40); p.add_argument("--max-gc-percent", type=float, default=60)
p.add_argument("--max-ns-accepted", type=int, default=0); p.add_argument("--max-poly-x", type=int, default=4)
a = p.parse_args()
if not (0 < a.min_primer_size <= a.optimal_primer_size <= a.max_primer_size):
    p.error("primer sizes must satisfy 0 < min <= optimal <= max")
if not (0 < a.min_tm <= a.optimal_tm <= a.max_tm):
    p.error("Tm values must satisfy 0 < min <= optimal <= max")
if not (0 <= a.min_gc_percent <= a.max_gc_percent <= 100):
    p.error("GC percentages must satisfy 0 <= min <= max <= 100")
seqs = fasta(a.reference)
with open(a.vcf) as src, open(a.output, "w") as dst:
    for line in src:
        if line.startswith("#"): continue
        chrom, pos, _, ref, alt, *_ = line.rstrip().split("\t")
        if chrom not in seqs or "," in alt: continue
        pos = int(pos); left = max(0, pos - 1 - a.flank); right = min(len(seqs[chrom]), pos - 1 + len(ref) + a.flank)
        template = seqs[chrom][left:right]; target_start = pos - 1 - left
        dst.write(
            f"SEQUENCE_ID={chrom}_{pos}\nSEQUENCE_TEMPLATE={template}\nSEQUENCE_TARGET={target_start},{len(ref)}\n"
            f"PRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\n"
            f"PRIMER_NUM_RETURN={a.num_return}\nPRIMER_PRODUCT_SIZE_RANGE={a.product_size_range}\n"
            f"PRIMER_OPT_SIZE={a.optimal_primer_size}\nPRIMER_MIN_SIZE={a.min_primer_size}\nPRIMER_MAX_SIZE={a.max_primer_size}\n"
            f"PRIMER_OPT_TM={a.optimal_tm}\nPRIMER_MIN_TM={a.min_tm}\nPRIMER_MAX_TM={a.max_tm}\n"
            f"PRIMER_MIN_GC={a.min_gc_percent}\nPRIMER_MAX_GC={a.max_gc_percent}\n"
            f"PRIMER_MAX_NS_ACCEPTED={a.max_ns_accepted}\nPRIMER_MAX_POLY_X={a.max_poly_x}\n=\n"
        )
