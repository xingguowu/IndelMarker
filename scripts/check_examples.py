#!/usr/bin/env python3
"""Validate actual example content, not just the existence of output files."""
import argparse
import csv
import hashlib
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from itertools import zip_longest
from pathlib import Path
from zipfile import ZipFile

import yaml

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fasta(path):
    result, name = {}, None
    for line in Path(path).open():
        if line.startswith(">"):
            name = line[1:].strip()
            require(name not in result, f"Duplicate FASTA identifier in {path}: {name}")
            result[name] = ""
        else:
            require(name is not None, f"Missing FASTA header: {path}")
            result[name] += line.strip()
    return result


def fastq(path):
    with Path(path).open() as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            sequence, plus, quality = [handle.readline().rstrip() for _ in range(3)]
            require(header.startswith("@") and plus.startswith("+") and sequence,
                    f"Malformed FASTQ record in {path}: {header}")
            require(len(sequence) == len(quality), f"FASTQ sequence/quality length differs: {path}")
            yield header[1:].split()[0], sequence


def check_pairs(paths):
    count = 0
    for first, second in zip_longest(fastq(paths["r1"]), fastq(paths["r2"])):
        require(first is not None and second is not None, "R1/R2 have different record counts")
        require(first[0].removesuffix("/1") == second[0].removesuffix("/2"),
                f"R1/R2 identifiers are not paired: {first[0]}, {second[0]}")
        count += 1
    require(count > 100, "The example requires sufficient read pairs to exercise variant calling")
    return count


def variants(path):
    samples, records = [], []
    for line in Path(path).open():
        if line.startswith("#CHROM"):
            samples = line.rstrip().split("\t")[9:]
        elif not line.startswith("#") and line.strip():
            fields = line.rstrip().split("\t")
            require(len(fields) >= 8, f"Malformed VCF record: {path}")
            records.append(fields)
    require(records, f"No variant records in {path}; a header-only VCF is not a passing example")
    return samples, records


def workbook(path):
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as archive:
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sheet.findall(".//s:row", ns):
        values = ["".join(cell.itertext()) for cell in row.findall("s:c", ns)]
        rows.append(values)
    require(len(rows) > 1, f"Primer workbook has no data rows: {path}")
    successful = [row for row in rows[1:] if len(row) >= 15 and row[6] and row[9] and row[14]]
    require(successful, f"Primer3 returned no usable primer pairs: {path}")
    return len(rows) - 1, len(successful)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["demo", "maize"], required=True)
    parser.add_argument("--strategy", choices=["1", "2", "3"], required=True)
    parser.add_argument("--inputs-only", action="store_true")
    args = parser.parse_args()
    directory = ROOT / ("demo" if args.suite == "demo" else "examples/maize")
    config = yaml.safe_load((directory / f"strategy{args.strategy}/config.yaml").read_text())
    # Repository-relative paths are the workflow's documented contract.
    import os
    os.chdir(ROOT)
    s = config[f"strategy{args.strategy}"]
    required = {
        "1": ["samtools", "minimap2", "syri", "primer3_core"],
        "2": ["samtools", "bwa", "gatk", "java", "fastp", "primer3_core"],
        "3": ["minimap2", "bedtools", "mafft"],
    }[args.strategy]
    for tool in required:
        require(shutil.which(tool), f"Missing executable: {tool}. Activate the complete Conda environment.")
    if args.suite == "maize":
        provenance = json.loads((directory / "provenance.json").read_text())
        for file, checksum in provenance["sha256"].items():
            path = directory / file
            require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == checksum,
                    f"Bundled input is missing or changed: {path}")
        public_sources = json.loads((directory / "source_verification.json").read_text())
        for file, details in public_sources.items():
            sequence = "".join(fasta(directory / "data" / file).values()).upper()
            require(hashlib.sha256(sequence.encode()).hexdigest() == details["sequence_sha256"],
                    f"Sequence differs from the verified public-source extract: {file}")
    if args.strategy == "1":
        for field in ("reference", "query"):
            records = fasta(s[field])
            require(s["chromosome"] in records, f"Requested chromosome missing in {s[field]}")
    elif args.strategy == "2":
        for parent, paths in s["parents"].items():
            check_pairs(paths)
    else:
        require(Path(s["annotation"]).is_file() and Path(s["longest_cds_fasta"]).is_file(),
                "Regional annotation or CDS file missing")
    if args.inputs_only:
        print(f"PASS: {args.suite} Strategy {args.strategy} input data and tools")
        return
    out = Path(config["outdir"])
    result = {"suite": args.suite, "strategy": args.strategy, "status": "PASS"}
    expected = json.loads((directory / "expected.json").read_text())[f"strategy{args.strategy}"]
    if args.strategy in {"1", "2"}:
        vcf = out / ("strategy1/filtered_indels.vcf" if args.strategy == "1" else "strategy2/contrasting_indels.vcf")
        samples, records = variants(vcf)
        if "position" in expected:
            matching = [f for f in records if f[0] == expected["chromosome"]
                        and abs(int(f[1]) - expected["position"]) <= 2
                        and abs(len(f[3]) - len(f[4])) == expected["indel_length"]]
            require(matching, f"Known example InDel was not recovered: {expected}")
            if args.strategy == "2":
                keys = matching[0][8].split(":")
                genotypes = {sample: fields.split(":")[keys.index("GT")].replace("|", "/")
                             for sample, fields in zip(samples, matching[0][9:])}
                require(genotypes == expected["genotypes"], f"Unexpected parental calls: {genotypes}")
        result["variant_records"] = len(records)
        result["workbook_rows"], result["loci_with_primers"] = workbook(out / "primer3/primer_markers.xlsx")
        if "variant_records" in expected:
            require(result["variant_records"] == expected["variant_records"], "Unexpected variant count")
    else:
        with (out / "strategy3/alignment/manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        require(rows, "No genes in the regional alignment manifest")
        require({row["gene"] for row in rows} == set(expected["genes"]), "Unexpected genes in manifest")
        for row in rows:
            aligned = fasta(row["alignment"])
            require(len({len(seq) for seq in aligned.values()}) == 1, "Unequal alignment lengths")
            for label in expected["genomes"]:
                require(any(f"_{label}::" in name for name in aligned),
                        f"Missing genome {label} in {row['alignment']}")
            if args.suite == "demo":
                sequences = [seq.replace("-", "") for name, seq in aligned.items() if "::" in name]
                require(abs(len(sequences[0]) - len(sequences[1])) == expected["indel_length"],
                        "The known regional InDel is absent from the genome alignments")
        result["aligned_genes"] = len(rows)
    (out / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
