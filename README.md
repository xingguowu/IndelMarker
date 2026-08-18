# InDelMarker: a unified Snakemake workflow for InDel marker development

InDelMarker implements three complementary workflows for identifying polymorphic InDels and preparing PCR marker designs or targeted sequence alignments.

## Installation

The supplied Conda environment provides all command-line programs used by the workflow. Python is pinned to 3.10 because SyRI 1.8.2 requires it; Snakemake is kept on the compatible 8.x release line. SyRI, Primer3, and Minimap2 provide `osx-arm64` builds through Bioconda ([SyRI](https://bioconda.github.io/recipes/syri/README.html), [Primer3](https://bioconda.github.io/recipes/primer3/README.html), [Minimap2](https://bioconda.github.io/recipes/minimap2/README.html)).

```bash
git clone https://github.com/xingguowu/IndelMarker.git
cd IndelMarker
conda config --set channel_priority strict
conda env create -f envs/indelmarker.yaml
conda activate indelmarker
```

## Quick start

Edit `config.yaml`, set `strategy` to `"1"`, `"2"`, or `"3"`, and provide the input paths required by that strategy. Then run:

```bash
snakemake --cores 16 --use-conda
```

Before starting a full analysis, verify that the environment and workflow can be resolved:

```bash
for program in snakemake samtools minimap2 syri bcftools bwa gatk fastp primer3_core exonerate bedtools mafft; do command -v "$program"; done
snakemake --snakefile Snakefile --configfile config.yaml --dry-run --cores 1
```

## Workflows

### Strategy 1: assembly-based InDel discovery

This workflow is designed for chromosome-scale genome assemblies. It extracts a matching chromosome from the reference and query assemblies, aligns them with `minimap2`, calls structural variation with SyRI, filters InDels by length, and prepares Primer3 input.

Outputs:

- `results/strategy1/filtered_indels.vcf`
- `results/primer3/primer3_input.txt`
- `results/primer3/primer3_output.txt`

The recommended InDel length range for gel-based markers is 300–500 bp.

### Strategy 2: short-read parental variant calling

This workflow processes paired-end reads from two parents with `fastp`, `bwa`, and GATK joint genotyping. It retains length-filtered InDels with contrasting parental genotypes and prepares Primer3 input.

Outputs:

- `results/strategy2/contrasting_indels.vcf`
- `results/primer3/primer3_input.txt`
- `results/primer3/primer3_output.txt`

Parents must be named `p1` and `p2` in `config.yaml`. Genotype mode `A` retains only `0/0` versus `1/1`; mode `B` also retains heterozygous versus `1/1` records. The recommended InDel length range is 20–100 bp.

### Strategy 3: targeted regional sequence alignment

This workflow supports fine-mapping intervals. It identifies B73 genes overlapping the selected interval, retrieves their longest CDS sequences, locates homologous regions in one or more selected genomes with Exonerate, extracts flanking sequences, and creates one MAFFT alignment per gene.

The `annotation` input accepts standard GFF3 or BED files, including the supplied `B73v5.gene.bed` format. Strategy 3 does not run Primer3; its alignments are intended for manual marker inspection and primer design.

Outputs:

- `results/strategy3/alignment/pos.tsv`
- `results/strategy3/alignment/manifest.tsv`
- `results/strategy3/alignment/<gene>.fa`

## Notes

- Strategy 2 automatically creates `samtools faidx`, BWA, and GATK sequence-dictionary indexes for the reference FASTA.
- Tool executables can be overridden in the `tools` section of `config.yaml` when a site-specific installation must be used.
- Large input and output files are excluded by `.gitignore`; use `config.local.yaml` for local paths if you maintain site-specific settings.
