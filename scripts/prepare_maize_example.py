#!/usr/bin/env python3
"""Rebuild the bundled B73/Ames21814 extracts from the published assemblies."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "demo"))
from generate_demo_inputs import fasta, paired_reads, write


def region(path, chromosome, start, end):
    chunks, position, active = [], 0, False
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if active:
                    break
                active = line[1:].split()[0] == chromosome
                continue
            if not active:
                continue
            sequence = line.strip()
            next_position = position + len(sequence)
            if next_position >= start and position < end:
                chunks.append(sequence[max(0, start - 1 - position):min(len(sequence), end - position)])
            position = next_position
            if position >= end:
                break
    sequence = "".join(chunks).upper()
    if len(sequence) != end - start + 1:
        raise ValueError(f"Cannot extract {chromosome}:{start}-{end} from {path}")
    return sequence


def attributes(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--b73", required=True)
    p.add_argument("--ames21814", required=True)
    p.add_argument("--gff3", required=True)
    p.add_argument("--output", default="examples/maize")
    a = p.parse_args()
    out = Path(a.output)
    start, end = 194750001, 195050000
    query_start, query_end = 199705286, 200083834
    reference = region(a.b73, "chr3", start, end)
    query = region(a.ames21814, "chr3", query_start, query_end)
    data = out / "data"
    fasta(data / "B73.region.fa", "chr3_example", reference)
    fasta(data / "Ames21814.region.fa", "chr3_example", query)
    genes, transcripts, coding = {}, {}, {}
    for line in Path(a.gff3).open():
        if line.startswith("#"):
            continue
        f = line.rstrip().split("\t")
        if len(f) != 9 or f[0] != "chr3" or int(f[3]) < start + 200 or int(f[4]) > end - 200:
            continue
        meta = attributes(f[8])
        if f[2] == "gene":
            genes[meta["ID"]] = (int(f[3]), int(f[4]), f[6])
        elif f[2] in {"mRNA", "transcript"}:
            transcripts[meta["ID"]] = meta.get("Parent", "")
        elif f[2] == "CDS":
            for parent in meta.get("Parent", "").split(","):
                coding.setdefault(parent, []).append((int(f[3]), int(f[4]), f[6]))
    longest = {}
    complement = str.maketrans("ACGTN", "TGCAN")
    for transcript, pieces in coding.items():
        gene = transcripts.get(transcript)
        if gene not in genes:
            continue
        bases = "".join(reference[s - start:e - start + 1] for s, e, _ in sorted(pieces))
        if pieces[0][2] == "-":
            bases = bases.translate(complement)[::-1]
        if len(bases) > len(longest.get(gene, "")):
            longest[gene] = bases
    if not longest:
        raise ValueError("No complete annotated genes were extracted")
    bed, cds = [], []
    original = ["gene\tchromosome\tstart\tend\tstrand\n"]
    for gene in sorted(longest, key=lambda g: genes[g][0]):
        s, e, strand = genes[gene]
        bed.append(f"chr3_example\t{s - start}\t{e - start + 1}\t{gene}\t0\t{strand}\n")
        cds.append(f">{gene}\n{longest[gene]}\n")
        original.append(f"{gene}\tchr3\t{s}\t{e}\t{strand}\n")
    write(data / "B73.genes.bed", "".join(bed))
    write(data / "B73.longest_cds.fa", "".join(cds))
    write(data / "B73.original_gene_coordinates.tsv", "".join(original))
    # Small read-based demonstration in the same real genomic background.
    # These are simulated reads, never described as sequenced WGS data.
    read_start, read_end = 1, 12000
    query_read_start, query_read_end = 1, 12000
    fasta(data / "B73.read_region.fa", "chr3_reads", reference[read_start - 1:read_end])
    read_sequences = {
        "B73_simulated": reference[read_start - 1:read_end],
        "Ames21814_simulated": query[query_read_start - 1:query_read_end],
    }
    read_pairs = {}
    for number, (label, bases) in enumerate(read_sequences.items()):
        fasta(data / f"{label}.haplotype.fa", "chr3_reads", bases)
        read_pairs[label] = paired_reads(bases, data / label, 20261015 + number)
    checksums = {f"data/{file.name}": hashlib.sha256(file.read_bytes()).hexdigest()
                 for file in sorted(data.iterdir())
                 if file.is_file() and file.suffix in {".fa", ".bed", ".tsv", ".fastq"}}
    provenance = {
        "reference": {"assembly": "Zm-B73-REFERENCE-NAM-5.0", "chromosome": "chr3",
                      "start": start, "end": end,
                      "source": "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/"},
        "query": {"assembly": "Zv-Ames21814-REFERENCE-CEMPS-1.0", "chromosome": "chr3",
                  "start": query_start, "end": query_end,
                  "source": "https://download.maizegdb.org/Zv-Ames21814-REFERENCE-CEMPS-1.0/"},
        "coordinate_system": "Original intervals: 1-based inclusive; bundled BED: 0-based half-open. Both extracted FASTAs renamed chr3_example.",
        "annotation": "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3; longest CDS per complete gene, extracted from B73 sequence.",
        "genes": list(longest), "sha256": checksums,
        "read_example": {
            "data_type": "simulated paired reads from unmodified public assembly extracts, not experimental WGS",
            "read_length": 150, "fragment_length_range": [350, 550], "step": 7,
            "read_pairs": read_pairs,
            "reference_original_interval": [start + read_start - 1, start + read_end - 1],
            "query_original_interval": [query_start + query_read_start - 1, query_start + query_read_end - 1],
            "seeds": [20261015, 20261016],
        },
        "interpretation": "Published genomic sequence extracts for computational demonstration; not PCR validation or a reproduction of the THP3 mapping experiment.",
    }
    write(out / "provenance.json", json.dumps(provenance, indent=2) + "\n")
    expected = {
        "strategy1": {"variant_records": 2},
        "strategy2": {"chromosome": "chr3_reads", "position": 3792, "indel_length": 42,
                      "genotypes": {"p1": "0/0", "p2": "1/1"}},
        "strategy3": {"genes": list(longest), "genomes": ["B73", "Ames21814"]},
    }
    write(out / "expected.json", json.dumps(expected, indent=2) + "\n")
    print(f"Prepared two real genome regions and {len(longest)} annotated genes.")


if __name__ == "__main__":
    main()
