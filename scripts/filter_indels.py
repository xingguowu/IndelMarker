#!/usr/bin/env python3
"""Keep biallelic insertion/deletion records by allele length and genotype."""
import argparse, gzip

def op(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def gt_ok(fields, mode):
    if not mode:
        return True
    # The workflow fixes sample order as p1,p2. A: 0/0 vs 1/1; B also allows het vs homozygous alt.
    gts = [x.split(":", 1)[0].replace("|", "/") for x in fields[9:11]]
    if len(gts) != 2: return False
    if mode == "A": return set(gts) == {"0/0", "1/1"}
    if mode == "B": return ("1/1" in gts and any(x in {"0/0", "0/1", "1/0"} for x in gts))
    raise ValueError("parental genotypes must be A or B")

def indel_length(fields):
    """Return absolute InDel length for both ordinary and symbolic SV VCFs."""
    ref, alt, info = fields[3], fields[4], fields[7] if len(fields) > 7 else ""
    meta = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
    svtype = meta.get("SVTYPE")
    if svtype and svtype not in {"INS", "DEL"}:
        return None
    if "SVLEN" in meta:
        try:
            return abs(int(meta["SVLEN"].split(",", 1)[0]))
        except ValueError:
            return None
    return abs(len(ref) - len(alt)) if len(ref) != len(alt) else None

p = argparse.ArgumentParser()
p.add_argument("--input", required=True); p.add_argument("--output", required=True)
p.add_argument("--min-length", type=int, required=True); p.add_argument("--max-length", type=int, required=True)
p.add_argument("--parental-genotypes", choices=["A", "B"])
a = p.parse_args()
with op(a.input) as src, open(a.output, "w") as dst:
    for line in src:
        if line.startswith("#"):
            dst.write(line); continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or "," in f[4] or f[4] in {".", "*"}: continue
        diff = indel_length(f)
        if diff is not None and a.min_length <= diff <= a.max_length and gt_ok(f, a.parental_genotypes): dst.write(line)
