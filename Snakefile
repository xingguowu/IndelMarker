"""Unified CornMarker-style InDel workflow (strategies 1–3)."""
from pathlib import Path

configfile: "config.yaml"

STRATEGY = str(config["strategy"])
OUT = config["outdir"].rstrip("/")
TOOLS = config.get("tools", {})
SAMTOOLS = TOOLS.get("samtools", "samtools")
MINIMAP2 = TOOLS.get("minimap2", "minimap2")
SYRI = TOOLS.get("syri", "syri")
BCFTOOLS = TOOLS.get("bcftools", "bcftools")
BWA = TOOLS.get("bwa", "bwa")
GATK = TOOLS.get("gatk", "gatk")
FASTP = TOOLS.get("fastp", "fastp")
PRIMER3 = TOOLS.get("primer3_core", "primer3_core")
EXONERATE = TOOLS.get("exonerate", "exonerate")
BEDTOOLS = TOOLS.get("bedtools", "bedtools")
MAFFT = TOOLS.get("mafft", "mafft")
PYTHON = TOOLS.get("python", "python3")
THREADS = int(config.get("threads", 8))

if STRATEGY not in {"1", "2", "3"}:
    raise ValueError("config.strategy must be 1, 2, or 3")

def primer_targets(vcf):
    return [f"{OUT}/primer3/primer3_input.txt", f"{OUT}/primer3/primer3_output.txt"]

if STRATEGY == "1":
    S1 = config["strategy1"]
    REF, QUERY, CHR = S1["reference"], S1["query"], S1["chromosome"]
    FINAL = f"{OUT}/strategy1/filtered_indels.vcf"
elif STRATEGY == "2":
    S2 = config["strategy2"]
    REF = S2["reference"]
    PARENTS = S2["parents"]
    if set(PARENTS) != {"p1", "p2"}:
        raise ValueError("strategy2.parents must contain exactly p1 and p2")
    FINAL = f"{OUT}/strategy2/contrasting_indels.vcf"
else:
    S3 = config["strategy3"]
    FINAL = f"{OUT}/strategy3/alignment/pos.tsv"

rule all:
    input:
        FINAL if STRATEGY == "3" else [FINAL] + primer_targets(FINAL)

# ---------------------------------------------------------------------------
# Strategy 1: chromosome-scale assemblies -> SyRI -> InDel primers
# ---------------------------------------------------------------------------
if STRATEGY == "1":
    rule s1_extract_chromosome:
        input: ref=REF, query=QUERY
        output: ref=f"{OUT}/strategy1/fasta/ref.{CHR}.fa", query=f"{OUT}/strategy1/fasta/query.{CHR}.fa"
        shell:
            "mkdir -p {OUT}/strategy1/fasta; {SAMTOOLS} faidx {input.ref} {CHR} > {output.ref}; {SAMTOOLS} faidx {input.query} {CHR} > {output.query}"

    rule s1_align:
        input: rules.s1_extract_chromosome.output
        output: f"{OUT}/strategy1/alignment/{CHR}.sam"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy1/alignment; {MINIMAP2} -ax {S1.get('preset', 'asm5')} --eqx -t {threads} {input.ref} {input.query} > {output}"

    rule s1_syri:
        input: sam=rules.s1_align.output, ref=rules.s1_extract_chromosome.output.ref, query=rules.s1_extract_chromosome.output.query
        output: vcf=f"{OUT}/strategy1/syri/{CHR}_syri.vcf", report=f"{OUT}/strategy1/syri/{CHR}_syri.out"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy1/syri; {SYRI} -c {input.sam} -r {input.ref} -q {input.query} -F S --dir {OUT}/strategy1/syri --prefix {CHR}_ --nc {threads}"

    rule s1_filter_indels:
        input: rules.s1_syri.output.vcf
        output: FINAL
        params: min_len=S1.get("min_indel_length", 300), max_len=S1.get("max_indel_length", 500)
        shell:
            "mkdir -p {OUT}/strategy1; {PYTHON} scripts/filter_indels.py --input {input} --output {output} --min-length {params.min_len} --max-length {params.max_len}"

# ---------------------------------------------------------------------------
# Strategy 2: paired-end parental WGS -> GATK -> contrasting InDel primers
# ---------------------------------------------------------------------------
if STRATEGY == "2":
    rule s2_reference_indexes:
        input: REF
        output: fai=REF + ".fai", bwa=REF + ".bwt", dictionary=str(Path(REF).with_suffix(".dict"))
        shell:
            "{SAMTOOLS} faidx {input}; {BWA} index {input}; {GATK} CreateSequenceDictionary -R {input} -O {output.dictionary}"

    rule s2_fastp:
        input: r1=lambda wc: PARENTS[wc.sample]["r1"], r2=lambda wc: PARENTS[wc.sample]["r2"]
        output: r1=f"{OUT}/strategy2/clean/{{sample}}.R1.fq.gz", r2=f"{OUT}/strategy2/clean/{{sample}}.R2.fq.gz", json=f"{OUT}/strategy2/qc/{{sample}}.json", html=f"{OUT}/strategy2/qc/{{sample}}.html"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy2/clean {OUT}/strategy2/qc; {FASTP} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} -j {output.json} -h {output.html}"

    rule s2_map:
        input: r1=rules.s2_fastp.output.r1, r2=rules.s2_fastp.output.r2, index=rules.s2_reference_indexes.output
        output: bam=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bam", bai=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bam.bai"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy2/bam {OUT}/strategy2/tmp; {BWA} mem -t {threads} {REF} {input.r1} {input.r2} | {SAMTOOLS} sort -@ {threads} -o {OUT}/strategy2/tmp/{{wildcards.sample}}.sorted.bam; {GATK} MarkDuplicates -I {OUT}/strategy2/tmp/{{wildcards.sample}}.sorted.bam -O {output.bam} -M {OUT}/strategy2/bam/{{wildcards.sample}}.metrics.txt --CREATE_INDEX true"

    rule s2_haplotypecaller:
        input: bam=rules.s2_map.output.bam, bai=rules.s2_map.output.bai, index=rules.s2_reference_indexes.output
        output: f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy2/gvcf; {GATK} HaplotypeCaller -R {REF} -I {input.bam} -ERC GVCF -O {output} --native-pair-hmm-threads {threads}"

    rule s2_joint_genotype:
        input: gvcfs=expand(f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz", sample=["p1", "p2"]), index=rules.s2_reference_indexes.output
        output: f"{OUT}/strategy2/joint.raw.vcf.gz"
        shell:
            "mkdir -p {OUT}/strategy2; {GATK} CombineGVCFs -R {REF} --variant {input.gvcfs[0]} --variant {input.gvcfs[1]} -O {OUT}/strategy2/joint.g.vcf.gz; {GATK} GenotypeGVCFs -R {REF} -V {OUT}/strategy2/joint.g.vcf.gz -O {output}"

    rule s2_filter_indels:
        input: rules.s2_joint_genotype.output
        output: FINAL
        params: min_len=S2.get("min_indel_length", 20), max_len=S2.get("max_indel_length", 100), genotypes=S2.get("genotype_mode", "A")
        shell:
            "{PYTHON} scripts/filter_indels.py --input {input} --output {output} --min-length {params.min_len} --max-length {params.max_len} --parental-genotypes {params.genotypes}"

# Shared primer3 steps for strategies 1 and 2.
if STRATEGY in {"1", "2"}:
    rule make_primer3_input:
        input: vcf=FINAL, reference=REF
        output: f"{OUT}/primer3/primer3_input.txt"
        params: flank=config.get("primer3", {}).get("flank", 400), product_range=config.get("primer3", {}).get("product_size_range", "250-700")
        shell:
            "mkdir -p {OUT}/primer3; {PYTHON} scripts/vcf_to_primer3.py --vcf {input.vcf} --reference {input.reference} --output {output} --flank {params.flank} --product-size-range {params.product_range}"

    rule run_primer3:
        input: rules.make_primer3_input.output
        output: f"{OUT}/primer3/primer3_output.txt"
        shell:
            "{PRIMER3} < {input} > {output}"

# ---------------------------------------------------------------------------
# Strategy 3: region genes -> CDS-guided genome extraction -> alignment only
# ---------------------------------------------------------------------------
if STRATEGY == "3":
    rule s3_region_alignments:
        input: gff=S3["annotation"] if "annotation" in S3 else S3["gff3"], cds=S3["longest_cds_fasta"], genomes=lambda wc: list(S3["genomes"].values())
        output: positions=f"{OUT}/strategy3/alignment/pos.tsv", manifest=f"{OUT}/strategy3/alignment/manifest.tsv"
        params: chromosome=S3["chromosome"], start=S3["start"], end=S3["end"], upstream=S3.get("upstream", 2000), downstream=S3.get("downstream", 400), genomes=lambda wc: ",".join(f"{k}={v}" for k,v in S3["genomes"].items())
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy3/alignment; {PYTHON} scripts/strategy3_align.py --gff3 {input.gff} --cds {input.cds} --genomes '{params.genomes}' --chromosome {params.chromosome} --start {params.start} --end {params.end} --upstream {params.upstream} --downstream {params.downstream} --outdir {OUT}/strategy3 --exonerate '{EXONERATE}' --bedtools '{BEDTOOLS}' --mafft '{MAFFT}' --threads {threads}"
