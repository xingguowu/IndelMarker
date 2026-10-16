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

# ``threads`` is the total core budget passed to Snakemake.  Do not allocate
# that whole budget to every parental sample: doing so serializes P1 and P2,
# even on a large server.  Users may override the derived values in
# ``strategy2`` when their scheduler has a different allocation model.
S2_FASTP_THREADS = max(1, THREADS // 2)
S2_MAP_THREADS = max(1, THREADS // 2)
S2_CALL_THREADS = max(1, min(8, THREADS // 2))

if STRATEGY not in {"1", "2", "3"}:
    raise ValueError("config.strategy must be 1, 2, or 3")

def primer_targets(vcf):
    return [f"{OUT}/primer3/primer_markers.xlsx"]

if STRATEGY == "1":
    S1 = config["strategy1"]
    REF, QUERY, CHR = S1["reference"], S1["query"], S1["chromosome"]
    S1_PRESET = S1.get("preset", "asm5")
    FINAL = f"{OUT}/strategy1/filtered_indels.vcf"
elif STRATEGY == "2":
    S2 = config["strategy2"]
    REF = S2["reference"]
    PARENTS = S2["parents"]
    if set(PARENTS) != {"p1", "p2"}:
        raise ValueError("strategy2.parents must contain exactly p1 and p2")
    S2_FASTP_THREADS = int(S2.get("fastp_threads", S2_FASTP_THREADS))
    S2_MAP_THREADS = int(S2.get("map_threads", S2_MAP_THREADS))
    S2_CALL_THREADS = int(S2.get("call_threads", S2_CALL_THREADS))
    S2_TMPDIR = S2.get("tmpdir", f"{OUT}/strategy2/tmp")
    FINAL = f"{OUT}/strategy2/contrasting_indels.vcf"
else:
    S3 = config["strategy3"]
    if "genomes" not in S3:
        legacy_genomes = {
            label: S3[key]
            for label, key in (("P1", "genome_qry1"), ("P2", "genome_qry2"))
            if key in S3
        }
        if len(legacy_genomes) != 2:
            raise ValueError("strategy3 requires a genomes mapping or both genome_qry1 and genome_qry2")
        S3["genomes"] = legacy_genomes
    S3_REGION_ARGS = f"--regions-file {S3['regions_file']}" if "regions_file" in S3 else f"--chromosome {S3['chromosome']} --start {S3['start']} --end {S3['end']}"
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
            "mkdir -p {OUT}/strategy1/fasta; {SAMTOOLS} faidx {input.ref}; {SAMTOOLS} faidx {input.query}; {SAMTOOLS} faidx {input.ref} {CHR} > {output.ref}; {SAMTOOLS} faidx {input.query} {CHR} > {output.query}"

    rule s1_align:
        input:
            ref=rules.s1_extract_chromosome.output.ref,
            query=rules.s1_extract_chromosome.output.query
        output: f"{OUT}/strategy1/alignment/{CHR}.sam"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy1/alignment; {MINIMAP2} -ax {S1_PRESET} --eqx -t {threads} {input.ref} {input.query} > {output}"

    rule s1_syri:
        input: sam=rules.s1_align.output, ref=rules.s1_extract_chromosome.output.ref, query=rules.s1_extract_chromosome.output.query
        output: vcf=f"{OUT}/strategy1/syri/{CHR}_syri.vcf", report=f"{OUT}/strategy1/syri/{CHR}_syri.out"
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy1/syri; {SYRI} -c {input.sam} -r {input.ref} -q {input.query} -F S --cigar --dir {OUT}/strategy1/syri --prefix {CHR}_ --nc {threads}"

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
        threads: S2_FASTP_THREADS
        shell:
            "mkdir -p {OUT}/strategy2/clean {OUT}/strategy2/qc; {FASTP} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} -j {output.json} -h {output.html}"

    rule s2_map:
        input: r1=rules.s2_fastp.output.r1, r2=rules.s2_fastp.output.r2, index=rules.s2_reference_indexes.output
        output: bam=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bam", bai=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bai"
        threads: S2_MAP_THREADS
        params: tmpdir=S2_TMPDIR
        shell:
            "mkdir -p {OUT}/strategy2/bam {params.tmpdir}; {BWA} mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' -t {threads} {REF} {input.r1} {input.r2} | {SAMTOOLS} sort -@ {threads} -T {params.tmpdir}/{wildcards.sample}.sort -o {params.tmpdir}/{wildcards.sample}.sorted.bam; {GATK} MarkDuplicates -I {params.tmpdir}/{wildcards.sample}.sorted.bam -O {output.bam} -M {OUT}/strategy2/bam/{wildcards.sample}.metrics.txt --CREATE_INDEX true"

    rule s2_all_bams_complete:
        input:
            expand(f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bam", sample=["p1", "p2"]),
            expand(f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bai", sample=["p1", "p2"])
        output: touch(f"{OUT}/strategy2/.all_bams_complete")

    rule s2_haplotypecaller:
        input: bam=rules.s2_map.output.bam, bai=rules.s2_map.output.bai, index=rules.s2_reference_indexes.output, all_bams=rules.s2_all_bams_complete.output
        output: f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz"
        threads: S2_CALL_THREADS
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
        params:
            flank=config.get("primer3", {}).get("flank", 400),
            product_range=config.get("primer3", {}).get("product_size_range", "250-700"),
            num_return=config.get("primer3", {}).get("num_return", 3),
            optimal_size=config.get("primer3", {}).get("optimal_primer_size", 24),
            min_size=config.get("primer3", {}).get("min_primer_size", 18),
            max_size=config.get("primer3", {}).get("max_primer_size", 26),
            optimal_tm=config.get("primer3", {}).get("optimal_tm", 60),
            min_tm=config.get("primer3", {}).get("min_tm", 57),
            max_tm=config.get("primer3", {}).get("max_tm", 60),
            min_gc=config.get("primer3", {}).get("min_gc_percent", 40),
            max_gc=config.get("primer3", {}).get("max_gc_percent", 60),
            max_ns=config.get("primer3", {}).get("max_ns_accepted", 0),
            max_poly_x=config.get("primer3", {}).get("max_poly_x", 4)
        shell:
            "mkdir -p {OUT}/primer3; {PYTHON} scripts/vcf_to_primer3.py --vcf {input.vcf} --reference {input.reference} --output {output} --flank {params.flank} --product-size-range {params.product_range} --num-return {params.num_return} --optimal-primer-size {params.optimal_size} --min-primer-size {params.min_size} --max-primer-size {params.max_size} --optimal-tm {params.optimal_tm} --min-tm {params.min_tm} --max-tm {params.max_tm} --min-gc-percent {params.min_gc} --max-gc-percent {params.max_gc} --max-ns-accepted {params.max_ns} --max-poly-x {params.max_poly_x}"

    rule run_primer3:
        input: rules.make_primer3_input.output
        output: f"{OUT}/primer3/primer3_output.txt"
        shell:
            "if [ -s {input} ]; then {PRIMER3} < {input} > {output}; else : > {output}; fi"

    rule export_primer_excel:
        input:
            vcf=FINAL,
            reference=REF,
            primer3=rules.run_primer3.output
        output:
            f"{OUT}/primer3/primer_markers.xlsx"
        params:
            flank=config.get("primer3", {}).get("flank", 400)
        shell:
            "{PYTHON} scripts/primer3_to_excel.py --vcf {input.vcf} --reference {input.reference} --primer3-output {input.primer3} --output {output} --flank {params.flank}"

# ---------------------------------------------------------------------------
# Strategy 3: region genes -> CDS-guided genome extraction -> alignment only
# ---------------------------------------------------------------------------
if STRATEGY == "3":
    rule s3_region_alignments:
        input: gff=S3["annotation"] if "annotation" in S3 else S3["gff3"], cds=S3["longest_cds_fasta"], genomes=lambda wc: list(S3["genomes"].values()), regions=lambda wc: [S3["regions_file"]] if "regions_file" in S3 else []
        output: positions=f"{OUT}/strategy3/alignment/pos.tsv", manifest=f"{OUT}/strategy3/alignment/manifest.tsv"
        params: region_args=S3_REGION_ARGS, upstream=S3.get("upstream", 2000), downstream=S3.get("downstream", 2000), genomes=lambda wc: ",".join(f"{k}={v}" for k,v in S3["genomes"].items())
        threads: THREADS
        shell:
            "mkdir -p {OUT}/strategy3/alignment; {PYTHON} scripts/strategy3_align.py --gff3 {input.gff} --cds {input.cds} --genomes '{params.genomes}' {params.region_args} --upstream {params.upstream} --downstream {params.downstream} --outdir {OUT}/strategy3 --exonerate '{EXONERATE}' --minimap2 '{MINIMAP2}' --bedtools '{BEDTOOLS}' --mafft '{MAFFT}' --threads {threads}"
