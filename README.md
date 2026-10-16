# InDelMarker: three workflows with runnable example data

InDelMarker develops candidate InDel markers from genome assemblies (Strategy I), paired-end reads (Strategy II), or local gene alignments (Strategy III). Strategies I and II produce an Excel primer workbook; Strategy III produces alignments for manual marker selection. Candidate designs still require experimental PCR validation.

## Download, install, and run

Use **Linux x86_64** with Conda installed for the complete three-strategy environment. The small examples use two cores; allow approximately 4 GB RAM and space for the Conda packages. The example runner defaults to a 2-GB Java heap unless `JAVA_TOOL_OPTIONS` is already set. Production maize analyses require substantially more memory, storage, and time.

```bash
git clone https://github.com/xingguowu/IndelMarker.git
cd IndelMarker
bash scripts/create_conda_env.sh
conda activate indelmarker
bash scripts/run_examples.sh --suite all --cores 2
```

Alternatively, download the repository with **Code → Download ZIP**, unzip it, open a terminal in the extracted folder, and start at `bash scripts/create_conda_env.sh`. No Git checkout, external example-data download, private server, or manual configuration editing is needed. Installation downloads software dependencies; the example inputs themselves are included in the repository.

The final message must be:

```text
All requested example strategies completed successfully.
```

Success means the output-content checks passed, including actual variant records, usable primer pairs, and alignments containing both genomes. A header-only VCF or an empty primer workbook is a failure.

## Two example suites

| Suite | Input | Purpose |
| --- | --- | --- |
| `demo` | Small deterministic synthetic assemblies, paired reads, and a gene model | Installation and end-to-end tests with known InDels |
| `maize` | Public B73/Ames21814 genome extracts and B73 gene annotations; the read-based example is explicitly labelled simulated | A realistic maize sequence example that runs without full-genome downloads |

Run either suite, or only one strategy:

```bash
bash demo/run_all.sh --cores 2
bash scripts/run_examples.sh --suite maize --cores 2
bash scripts/run_examples.sh --suite maize --strategy 1 --cores 2
```

The scripts also work when invoked from outside the repository. Repeat a command to resume an interrupted analysis; completed jobs are retained. There is no automatic deletion of earlier results.

Read [the synthetic example guide](demo/README.md) and [the real maize sequence guide](examples/maize/README.md) for provenance, coordinates, expected results, and interpretation limits.

## Where to find results

| Strategy | Synthetic suite | Maize suite |
| --- | --- | --- |
| I: primer workbook | `demo/work/strategy1/primer3/primer_markers.xlsx` | `examples/maize/results/strategy1/primer3/primer_markers.xlsx` |
| II: primer workbook | `demo/work/strategy2/primer3/primer_markers.xlsx` | `examples/maize/results/strategy2/primer3/primer_markers.xlsx` |
| III: alignment manifest | `demo/work/strategy3/strategy3/alignment/manifest.tsv` | `examples/maize/results/strategy3/strategy3/alignment/manifest.tsv` |

Each strategy also writes `validation.json` under its output directory. The workbook contains the locus, alleles, InDel length, template sequence, and up to three primer pairs with predicted reference product sizes. It does not establish PCR success or genome-wide primer specificity.

## Run your own data

Copy and edit the configuration corresponding to your data, then run from the repository root:

```bash
cp config.yaml config.local.yaml
# Edit config.local.yaml with your own input paths, strategy, and output directory.
snakemake --snakefile Snakefile --configfile config.local.yaml \
  --cores 16 --printshellcmds --rerun-incomplete --latency-wait 60
```

- **Strategy I:** two sufficiently contiguous assemblies with the same requested chromosome identifier. Default InDel filter: 300–500 bp. Alignment quality and assembly gaps affect calls.
- **Strategy II:** paired reads for two diploid parents, named `p1` and `p2`, against a common reference. Default filter: 20–100 bp and `0/0` versus `1/1`. This configuration is not a pooled BSA calling model. Variant quality, depth, primer-binding conservation, and PCR performance require additional assessment for biological use.
- **Strategy III:** BED/GFF3 annotation, longest-CDS FASTA, two genome FASTAs, and a target interval. It writes per-gene MAFFT alignments and does not automatically design primers.

The default Primer3 settings are 400-bp flanks, 250–700-bp reference products, up to three pairs, 18–26-bp primers (optimum 24), Tm 57–60 °C, GC 40–60%, no Ns, and maximum poly-X of four. Settings use the snake_case keys shown in `config.yaml`.

## Troubleshooting

| Problem | Action |
| --- | --- |
| `conda` not found | Install Miniforge/Miniconda for your Linux architecture, reopen the shell, and repeat installation. |
| Invalid inherited Conda channel | Use `scripts/create_conda_env.sh`, which explicitly sets channels for this installation. |
| `snakemake`, `gatk`, or `java` missing | Run `conda activate indelmarker`; use the complete environment rather than mixing installations. |
| Environment already exists | Use `bash scripts/create_conda_env.sh indelmarker-examples`, then activate that name. |
| Input missing after ZIP download | Extract the complete repository, including `demo/` and `examples/maize/`. |
| Interrupted run | Repeat the same example command; `--rerun-incomplete` is already included. |
| Changed or corrupted maize inputs | Restore them from a fresh download; the validator checks SHA-256 hashes. |
| No variants or no usable primers | Validation fails explicitly. Keep the log and check the configuration and software versions. |

The complete environment targets Linux. On Apple Silicon, use a Linux host for the supported installation; native local runs require an independently provisioned compatible environment. Never force-install Linux binary packages into a macOS Conda environment.

More detail: [reproducibility guide](docs/REPRODUCIBILITY.md).
