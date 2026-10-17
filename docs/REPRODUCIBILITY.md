# Reproducibility

## Downloaded examples

The repository includes all inputs, configurations, deterministic synthetic-data generation code, and output-content checks. A ZIP download has the same inputs as a Git clone.

```bash
bash scripts/create_conda_env.sh
conda activate indelmarker
bash scripts/run_examples.sh --suite all --cores 2
```

Software is version-constrained in `envs/indelmarker.yaml` and `envs/indelmarker-linux-64.txt`. The prebuilt Linux runtime additionally supplies `envs/indelmarker-linux-64.lock`, an exact URL/build list with package hashes generated in a glibc-2.28 Linux environment. The online installer uses this list to skip solving; the offline release includes the complete relocatable runtime and `software/packages.json`. Record `conda list --explicit` with any manuscript analysis. PuLP and setuptools constraints retain compatibility with Snakemake 7.32.4; Java 17 is required by GATK. The offline runtime uses `snakemake-minimal` for local workflows without optional cloud integrations.

## Input provenance and coordinates

- `demo/expected.json`: truth for the synthetic examples.
- `examples/maize/provenance.json`: source assembly versions, original genomic intervals, annotation origin, coordinate conventions, and SHA-256 hashes of bundled files.
- `examples/maize/data/B73.original_gene_coordinates.tsv`: original B73 chromosome coordinates for the genes retained in the regional case.
- `examples/maize/expected.json`: verified outputs required by the maize example validator.

The real assembly extracts retain unmodified nucleotide sequences and are renamed `chr3_example` to give the two local regions a matching FASTA identifier. For a local position `p`, the original coordinate is `original_start + p - 1`. The extracted reference and query have different original starts; apply the correct offset to each. These local extracts are suitable for demonstrating the workflow, not for assessing genome-wide rearrangement context or primer uniqueness.

## Verification

The runner validates input presence, R1/R2 pairing, read count, tool availability, and maize input checksums before starting. Afterwards it requires:

1. Non-empty InDel VCF records and recovery of known synthetic variants.
2. The expected parental genotypes for the synthetic read case.
3. Actual primer sequences and predicted product sizes in an Excel workbook.
4. The expected gene set and both genome sequences in each regional alignment.

The checks write `validation.json` to each strategy output directory and exit nonzero on failure. These checks verify computational examples; they do not measure experimental PCR success.

## Regenerate public maize extracts

The complete source assemblies and B73 GFF3 are available from the links in the maize guide. After downloading them, the bundled extracts can be rebuilt with:

```bash
python scripts/prepare_maize_example.py \
  --b73 /path/to/Zm-B73-REFERENCE-NAM-5.0.fa \
  --ames21814 /path/to/Zv-Ames21814-REFERENCE-CEMPS-1.0.fa \
  --gff3 /path/to/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
```

Full assemblies are not required to run the already bundled examples. The maize read example uses simulated reads generated from the specified public haplotype extracts and must not be described as raw experimental WGS.
