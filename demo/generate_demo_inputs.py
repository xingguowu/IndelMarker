#!/usr/bin/env python3
"""Create deterministic, non-empty examples; no biological data are simulated as real."""
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sequence(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def fasta(path, name, bases):
    write(path, f">{name}\n" + "\n".join(bases[i:i + 80] for i in range(0, len(bases), 80)) + "\n")


def paired_reads(bases, prefix, seed, step=7, read_length=150):
    """Tile unique, variable-length fragments; mates have matching identifiers."""
    rng = random.Random(seed)
    complement = str.maketrans("ACGT", "TGCA")
    mates = [[], []]
    for number, start in enumerate(range(0, len(bases) - 550, step), 1):
        insert = rng.randint(350, 550)
        first = bases[start:start + read_length]
        second = bases[start + insert - read_length:start + insert].translate(complement)[::-1]
        for mate, read in enumerate((first, second), 1):
            mates[mate - 1].append(f"@{prefix.name}_{number}/{mate}\n{read}\n+\n{'I' * len(read)}\n")
    for mate in (1, 2):
        write(Path(f"{prefix}.R{mate}.fastq"), "".join(mates[mate - 1]))
    return len(mates[0])


def main():
    reference = sequence(50_000, 20260818)
    inserted = sequence(360, 20260819)
    # Prevent a left-normalisation ambiguity at the insertion boundary.
    inserted = inserted[:-1] + next(b for b in "ACGT" if b != reference[24999])
    fasta(ROOT / "strategy1/reference.fa", "chr01", reference)
    query = reference[:25000] + inserted + reference[25000:]
    # Include SNP background as well as the target InDel, as in real assemblies.
    # SyRI 1.8.2 drops a terminal isolated InDel when no later variant flushes it.
    changed = next(b for b in "ACGT" if b != query[40360])
    query = query[:40360] + changed + query[40361:]
    fasta(ROOT / "strategy1/query.fa", "chr01", query)

    reference = sequence(12_000, 20261015)
    position = 6000
    while reference[position - 1] == reference[position + 39]:
        position += 1
    alternate = reference[:position] + reference[position + 40:]
    fasta(ROOT / "strategy2/reference.fa", "chr01", reference)
    pairs = {}
    for parent, bases, seed in (("p1", reference, 17), ("p2", alternate, 29)):
        pairs[parent] = paired_reads(bases, ROOT / "strategy2" / parent, seed)

    # A single unambiguous, 1.2-kb CDS in a 12-kb region, with a real alignment gap.
    genome = sequence(12_000, 314159)
    cds = "ATG" + sequence(1194, 271828) + "TAA"
    genome = genome[:4000] + cds + genome[5200:]
    alt = genome[:4500] + sequence(30, 12345) + genome[4500:]
    fasta(ROOT / "strategy3/genome_p1.fa", "chr01", genome)
    fasta(ROOT / "strategy3/genome_p2.fa", "chr01", alt)
    fasta(ROOT / "strategy3/longest_cds.fa", "gene_demo", cds)
    write(ROOT / "strategy3/annotation.bed", "chr01\t4000\t5200\tgene_demo\t0\t+\n")
    expected = {
        "data_type": "deterministic synthetic installation examples; not experimental validation",
        "strategy1": {"chromosome": "chr01", "position": 25000, "indel_length": 360},
        "strategy2": {"chromosome": "chr01", "position": position, "indel_length": 40,
                      "genotypes": {"p1": "0/0", "p2": "1/1"}, "read_pairs": pairs},
        "strategy3": {"genes": ["gene_demo"], "genomes": ["P1", "P2"], "indel_length": 30},
    }
    write(ROOT / "expected.json", json.dumps(expected, indent=2) + "\n")
    print("Demo inputs are ready (three deterministic synthetic strategies).")


if __name__ == "__main__":
    main()
