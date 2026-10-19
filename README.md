# InDelMarker: three workflows with runnable example data

InDelMarker develops candidate InDel markers from genome assemblies (Strategy I), paired-end reads or prealigned BAMs (Strategy II), or local gene alignments (Strategy III). Strategies I and II produce an Excel primer workbook; Strategy III produces alignments for manual marker selection. Candidate designs still require experimental PCR validation.

## Download, install, and run

### Offline installation

Download the complete **Linux x86_64 / glibc >= 2.28 offline bundle** and its checksum from [Releases](https://github.com/xingguowu/IndelMarker/releases). Copy them to the offline computer, verify the checksum, and extract the bundle. Inside its directory, run:

```bash
bash scripts/install_offline_linux.sh
bash scripts/run_offline_examples.sh --suite all --cores 2
```

The bundle includes Python, Java, every workflow program, both example datasets, exact package/build records, and upstream notices. The offline computer needs no Conda installation, dependency solving, or Internet connection. See [offline installation and compatibility](docs/OFFLINE_INSTALL.md). The ordinary **Code → Download ZIP** source archive contains no binary runtime.

### Online installation

Use **Linux x86_64** with Conda installed for the complete three-strategy environment. The small examples use two cores; allow approximately 4 GB RAM and space for the Conda packages. The example runner defaults to a 2-GB Java heap unless `JAVA_TOOL_OPTIONS` is already set. Production maize analyses require substantially more memory, storage, and time.

```bash
git clone https://github.com/xingguowu/IndelMarker.git
cd IndelMarker
bash scripts/create_conda_env.sh
# If Conda was not initialized in this terminal, first run the source command
# printed by the installer, then activate the environment.
conda activate indelmarker
bash scripts/run_examples.sh --suite all --cores 2
```

Alternatively, download the repository with **Code → Download ZIP**, unzip it, open a terminal in the extracted folder, and start at `bash scripts/create_conda_env.sh`. For the original demo and maize installation tests, no Git checkout, external data download, private server, or manual configuration editing is needed. Their inputs are included in the repository. The formal THP3 case is downloaded separately from Zenodo; its large inputs are excluded from the default source archive. Installation downloads software dependencies unless the offline runtime is used.

On a server, Conda may be installed but absent from the login shell's `PATH`. The installer also checks `CONDA_EXE`, `CONDA_PREFIX`, common installations under your home directory, and `/opt/conda`. It prints the exact `source .../etc/profile.d/conda.sh` command needed before `conda activate`. For a site-specific installation, pass its executable explicitly:

```bash
INDELMARKER_CONDA_EXE=/actual/installation/bin/conda bash scripts/create_conda_env.sh
```

If Conda is still missing, check `module avail` and load the Conda/Anaconda/Miniforge module provided by your cluster. If none is available, install [Miniforge from its official instructions](https://github.com/conda-forge/miniforge#install) in a new directory you own, then retry. No administrator privileges are needed.

The online installer uses `envs/indelmarker-linux-64.lock` when supplied. This exact package/build list skips dependency solving and downloads the already selected Linux packages. The locked runtime requires glibc 2.28 or newer.

The final message must be:

```text
All requested example strategies completed successfully.
```

Success means the output-content checks passed, including actual variant records, usable primer pairs, and alignments containing both genomes. A header-only VCF or an empty primer workbook is a failure.

## Two small example suites

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

An additional [THP3-centred case](examples/thp3_20Mb/README.md) uses B73/Ames21814 sequence extracts with 20-Mb flanks, 1,097 B73 gene models and original-depth regional BAMs. All three strategies passed on native Linux with the published offline runtime. The formal data, including reference sequences, BAMs, indexes and selected expected outputs, are archived separately in [Zenodo (DOI: 10.5281/zenodo.23474945)](https://doi.org/10.5281/zenodo.23474945), version 2026.10.19.

**Complete single-file download:** [THP3 full-depth case archive, 1.37 GB](https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz) and [SHA-256 file](https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256). This is a direct public HTTPS download; no account, API token or manual part assembly is required. It is byte-for-byte identical to the archive preserved on Zenodo. For interrupted downloads, use:

```bash
curl -fL -C - -o InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz \
  https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz
curl -fL -o InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256 \
  https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256
sha256sum -c InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256
```

For offline use, extract the base runtime first and this complete case archive second, following the [case installation guide](examples/thp3_20Mb/README.md#download-and-offline-installation). For an existing source checkout, the complete case archive can be extracted into a separate folder and its `examples/thp3_20Mb/data/`, `data_payload/` and `reference_results/` directories copied into the checkout. The case runner restores and verifies the three reference FASTAs.

Alternatively, retrieve the archived Zenodo copy using Python 3.9 or newer:


```bash
python3 scripts/download_thp3_data.py --install-data
# After activating the workflow environment:
bash scripts/run_thp3_examples.sh all 32
```

The downloader retrieves all 88 lossless archive parts and seven companion files, checks SHA-256 values, assembles the original complete archive, and installs only dataset files. It preserves code and configurations, reuses verified files on retry, and needs no account or API token. Allow approximately 5 GB of free disk space. For an offline machine, download and assemble the case beforehand and follow the [offline case instructions](examples/thp3_20Mb/README.md#download-and-offline-installation).

**Repository layout:** GitHub's default branch contains code, small installation tests, case configurations, provenance and validation summaries. Large THP3 inputs and result tables are distributed separately through the complete case release and the Zenodo archive. The [offline software environment](https://github.com/xingguowu/IndelMarker/releases/tag/offline-2026.10.17) remains on GitHub Releases. Earlier Git commits and the fixed scientific source snapshot still contain the historical reference payloads; ordinary `git clone` may download that history. For a lightweight initial download, use **Code → Download ZIP** or `git clone --depth 1 https://github.com/xingguowu/IndelMarker.git`.

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

With the offline bundle, use its command runner instead of activating Conda:

```bash
bash scripts/run_offline.sh snakemake --snakefile Snakefile --configfile config.local.yaml \
  --cores 16 --printshellcmds --rerun-incomplete --latency-wait 60
```

- **Strategy I:** two sufficiently contiguous assemblies with the same requested chromosome identifier. Default InDel filter: 300–500 bp. Alignment quality and assembly gaps affect calls.
- **Strategy II:** paired reads or existing coordinate-sorted, indexed BAMs for two diploid parents, named `p1` and `p2`, against a common reference. For existing alignments, set `strategy2.input_type: bam` and provide `bam`/`bai` paths for both parents; this skips FASTQ cleaning, mapping and duplicate marking. Default filter: 20–100 bp and `0/0` versus `1/1`. Optional `strategy2.intervals` restricts variant calling while retaining mapping against the complete reference. This configuration is not a pooled BSA calling model. Variant quality, depth, primer-binding conservation, and PCR performance require additional assessment for biological use.
- **Strategy III:** BED/GFF3 annotation, longest-CDS FASTA, two genome FASTAs, and a target interval. It writes per-gene MAFFT alignments and does not automatically design primers.

The default Primer3 settings are 400-bp flanks, 250–700-bp reference products, up to three pairs, 18–26-bp primers (optimum 24), Tm 57–60 °C, GC 40–60%, no Ns, and maximum poly-X of four. Settings use the snake_case keys shown in `config.yaml`.

## Troubleshooting

| Problem | Action |
| --- | --- |
| `conda` not found | Load your cluster's Conda module, specify `INDELMARKER_CONDA_EXE`, or install Miniforge following the linked official instructions. |
| Invalid inherited Conda channel | Use `scripts/create_conda_env.sh`, which explicitly sets channels for this installation. |
| Stuck at `Solving environment` | Use the updated installer with the exact `.lock` file, or the complete offline release bundle. Installing from version specifications alone still requires solving. |
| `snakemake`, `gatk`, or `java` missing | Run `conda activate indelmarker`; use the complete environment rather than mixing installations. |
| Environment already exists | Use `bash scripts/create_conda_env.sh indelmarker-examples`, then activate that name. |
| Input missing after ZIP download | Extract the complete repository, including `demo/` and `examples/maize/`. |
| Interrupted run | Repeat the same example command; `--rerun-incomplete` is already included. |
| Changed or corrupted maize inputs | Restore them from a fresh download; the validator checks SHA-256 hashes. |
| No variants or no usable primers | Validation fails explicitly. Keep the log and check the configuration and software versions. |

The complete environment targets Linux. On Apple Silicon, use a Linux host for the supported installation; native local runs require an independently provisioned compatible environment. Never force-install Linux binary packages into a macOS Conda environment.

More detail: [reproducibility guide](docs/REPRODUCIBILITY.md).
