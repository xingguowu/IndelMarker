"""Unified CornMarker-style InDel workflow (strategies 1–3)."""
from pathlib import Path
import shlex

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
    S2_ALL_BAMS = all("bam" in parent for parent in PARENTS.values())
    S2_INPUT_TYPE = S2.get("input_type", "bam" if S2_ALL_BAMS else "fastq")
    if S2_INPUT_TYPE not in {"bam", "fastq"}:
        raise ValueError("strategy2.input_type must be bam or fastq")
    S2_BAM_MODE = S2_INPUT_TYPE == "bam"
    if (any("bam" in parent for parent in PARENTS.values()) and not S2_ALL_BAMS) or (S2_BAM_MODE and not S2_ALL_BAMS):
        raise ValueError("strategy2 parents must both use BAM or both use paired FASTQ")
    # Snakemake recursively merges config.yaml defaults with --configfile.
    # Explicit BAM mode must ignore inherited placeholder FASTQ paths.
    S2_DOWNSAMPLING_ARGS = ""
    if "max_reads_per_alignment_start" in S2:
        limit = int(S2["max_reads_per_alignment_start"])
        if limit < 0:
            raise ValueError("max_reads_per_alignment_start must be nonnegative; 0 disables downsampling")
        S2_DOWNSAMPLING_ARGS = f"--max-reads-per-alignment-start {limit}"
    S2_FASTP_THREADS = int(S2.get("fastp_threads", S2_FASTP_THREADS))
    S2_MAP_THREADS = int(S2.get("map_threads", S2_MAP_THREADS))
    S2_CALL_THREADS = int(S2.get("call_threads", S2_CALL_THREADS))
    S2_TMPDIR = S2.get("tmpdir", f"{OUT}/strategy2/tmp")
    # Align to the complete reference, but optionally restrict variant calling
    # to selected regions. This retains genome-wide mapping competition.
    S2_INTERVALS = S2.get("intervals", [])
    if isinstance(S2_INTERVALS, str):
        S2_INTERVALS = [S2_INTERVALS]
    S2_INTERVAL_ARGS = " ".join("-L " + shlex.quote(str(region)) for region in S2_INTERVALS)
    S2_CHUNKS = {}
    if "calling_chunk_bp" in S2:
        chunk_bp = int(S2["calling_chunk_bp"])
        if chunk_bp < 1 or not S2_INTERVALS:
            raise ValueError("calling_chunk_bp requires a positive size and explicit contig:start-end intervals")
        for region in S2_INTERVALS:
            try:
                contig, bounds = region.rsplit(":", 1)
                start, end = map(int, bounds.split("-"))
                if start < 1 or start > end:
                    raise ValueError()
            except ValueError:
                raise ValueError("Chunked calling requires contig:start-end intervals with 1-based coordinates")
            for lo in range(start, end + 1, chunk_bp):
                S2_CHUNKS[f"chunk{len(S2_CHUNKS):05}"] = f"{contig}:{lo}-{min(end, lo + chunk_bp - 1)}"
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
    if S2_BAM_MODE:
        rule s2_reference_indexes:
            input: REF
            output: fai=REF + ".fai", dictionary=str(Path(REF).with_suffix(".dict"))
            shell:
                "{SAMTOOLS} faidx {input:q}; {GATK} CreateSequenceDictionary -R {input:q} -O {output.dictionary:q}"

        rule s2_map:
            input:
                bam=lambda wc: PARENTS[wc.sample]["bam"],
                bai=lambda wc: PARENTS[wc.sample].get("bai", PARENTS[wc.sample]["bam"] + ".bai"),
                index=rules.s2_reference_indexes.output
            output: bam=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bam", bai=f"{OUT}/strategy2/bam/{{sample}}.sorted.dedup.bai"
            shell:
                "{PYTHON} scripts/validate_bam_input.py --reference {REF:q} --bam {input.bam:q} --index {input.bai:q}; mkdir -p {OUT}/strategy2/bam; cp {input.bam:q} {output.bam:q}; cp {input.bai:q} {output.bai:q}"
    else:
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

    if S2_CHUNKS:
        rule s2_haplotypecaller_chunk:
            input: bam=rules.s2_map.output.bam, bai=rules.s2_map.output.bai, index=rules.s2_reference_indexes.output, all_bams=rules.s2_all_bams_complete.output
            output: vcf=f"{OUT}/strategy2/gvcf/{{sample}}/{{chunk}}.g.vcf.gz", tbi=f"{OUT}/strategy2/gvcf/{{sample}}/{{chunk}}.g.vcf.gz.tbi"
            threads: S2_CALL_THREADS
            params: interval=lambda wc: S2_CHUNKS[wc.chunk], downsampling=S2_DOWNSAMPLING_ARGS
            shell:
                "mkdir -p {OUT}/strategy2/gvcf/{wildcards.sample}; {GATK} HaplotypeCaller -R {REF} -I {input.bam} -ERC GVCF -O {output.vcf} --native-pair-hmm-threads {threads} -L {params.interval:q} {params.downsampling}"

        rule s2_haplotypecaller:
            input: lambda wc: expand(f"{OUT}/strategy2/gvcf/{{sample}}/{{chunk}}.g.vcf.gz", sample=wc.sample, chunk=list(S2_CHUNKS))
            output: vcf=f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz", tbi=f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz.tbi"
            params: variants=lambda wc, input: " ".join("-I " + shlex.quote(str(f)) for f in input)
            shell:
                "{GATK} GatherVcfs {params.variants} -O {output.vcf}.gathered.vcf; {BCFTOOLS} view --no-version -Oz -o {output.vcf} {output.vcf}.gathered.vcf; {BCFTOOLS} index -t {output.vcf}; rm {output.vcf}.gathered.vcf"
    else:
        rule s2_haplotypecaller:
            input: bam=rules.s2_map.output.bam, bai=rules.s2_map.output.bai, index=rules.s2_reference_indexes.output, all_bams=rules.s2_all_bams_complete.output
            output: vcf=f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz", tbi=f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz.tbi"
            threads: S2_CALL_THREADS
            params: intervals=S2_INTERVAL_ARGS, downsampling=S2_DOWNSAMPLING_ARGS
            shell:
                "mkdir -p {OUT}/strategy2/gvcf; {GATK} HaplotypeCaller -R {REF} -I {input.bam} -ERC GVCF -O {output.vcf} --native-pair-hmm-threads {threads} {params.intervals} {params.downsampling}"

    rule s2_joint_genotype:
        input: gvcfs=expand(f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz", sample=["p1", "p2"]), tbi=expand(f"{OUT}/strategy2/gvcf/{{sample}}.g.vcf.gz.tbi", sample=["p1", "p2"]), index=rules.s2_reference_indexes.output
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
