# THP3 regions and full-depth B73/Ames21814 BAMs

This additional case uses actual public genome sequences spanning THP3 and 20,000,000 bp beyond each gene boundary. It is separate from the small `demo` and `maize` installation tests. Those original tests remain unchanged.

| Material | Assembly | THP3 gene model | Gene coordinates, chr3 | Included interval, chr3 | Extract length |
| --- | --- | --- | --- | --- | --- |
| B73 | Zm-B73-REFERENCE-NAM-5.0 | Zm00001eb152450 | 200,402,538–200,408,213 | 180,402,538–220,408,213 | 40,005,676 bp |
| Ames21814 | Zv-Ames21814-REFERENCE-CEMPS-1.0 | Zv00003aa016390, corresponding gene supported by B73 gene alignment | 207,896,121–207,901,230 | 187,896,121–227,901,230 | 40,005,110 bp |

Assembly coordinates are **1-based inclusive**. Both extracted FASTA records are named `chr3_THP3`; their positions start at 1. Convert a local position to its original assembly position by adding 180,402,537 for B73 or 187,896,120 for Ames21814. BED coordinates are **0-based half-open**. The two equally wide flanking windows need not contain identical sets of genes because their sequence lengths and local genome structures differ.

Source genomes and annotations: [B73 MaizeGDB directory](https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/) and [Ames21814 MaizeGDB directory](https://download.maizegdb.org/Zv-Ames21814-REFERENCE-CEMPS-1.0/). Exact input hashes, source coordinates and the extraction method are recorded in `provenance.json` and `source_verification.json`. The THP3 gene alignment used to locate the Ames21814 region is retained in `sources/THP3_B73_to_Ames21814.paf`.

## Download and offline installation

The complete `InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz` archive contains the source code, original-depth regional BAMs, indexes and all reference inputs for this case. It requires no access to the authors' server. It uses the same top-level folder as the earlier offline release so that the two archives can be extracted in the same directory.

**Direct complete download:** [InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz (1.37 GB)](https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz), with its [SHA-256 checksum](https://github.com/xingguowu/IndelMarker/releases/download/thp3-full-depth-2026.10.19/InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256). No account or manual part assembly is needed. The complete archive is identical to the Zenodo-preserved archive and includes all inputs, indexes, selected outputs and the tested software snapshot. Copy it with the base runtime to an offline machine, then use the extraction commands below.

The formal dataset is [Zenodo, DOI: 10.5281/zenodo.23474945](https://doi.org/10.5281/zenodo.23474945), version 2026.10.19. Its large reference sequences, BAMs, indexes and result tables are excluded from GitHub's default source archive. On a connected machine, run from the checkout root (Python >=3.9, standard library only):

```bash
python3 scripts/download_thp3_data.py --install-data
```

This downloads and verifies all 95 Zenodo files, restores the complete archive and installs only data, compressed reference payloads and reference result files into this checkout. Code and configurations are preserved. The default download folder is `.downloads/thp3_20261019/`; verified files are reused on retry. Allow about 5 GB free space. Alternatively, download all 88 parts and the seven companion files from Zenodo and run the supplied `bash assemble_thp3_case.sh`.

For a machine without the runtime, also obtain the verified base offline archive and checksum from [offline-2026.10.17](https://github.com/xingguowu/IndelMarker/releases/tag/offline-2026.10.17). Copy the base archive, assembled case archive and both checksums to the offline Linux x86_64 / glibc >=2.28 machine. In the folder containing those files, run:

```bash
sha256sum -c InDelMarker-offline-linux-x86_64-glibc2.28.tar.gz.sha256
sha256sum -c InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz.sha256
tar -xzf InDelMarker-offline-linux-x86_64-glibc2.28.tar.gz
tar -xzf InDelMarker-THP3-full-depth-case-2026.10.19.tar.gz
cd InDelMarker-offline-linux-x86_64-glibc2.28
bash scripts/install_offline_linux.sh
bash scripts/run_offline.sh bash scripts/run_thp3_examples.sh all 32
```

If this base runtime is already installed, extract the THP3 archive beside the existing top-level folder and run the last command from that folder. With an independently installed compatible environment, the THP3 archive alone supplies the code and data; use `bash scripts/run_thp3_examples.sh all 32`.

## Run Strategies I and III after installing the data

First install the Zenodo dataset with the command above. Activate the complete InDelMarker environment, then run from the repository root:

```bash
bash scripts/run_thp3_examples.sh assemblies 8
```

To run only one strategy:

```bash
bash scripts/run_thp3_examples.sh 1 8
bash scripts/run_thp3_examples.sh 3 8
```

The offline runtime can run the same commands with `bash scripts/run_offline.sh bash scripts/run_thp3_examples.sh assemblies 8`, provided this source revision and its regional inputs have been copied into the offline bundle. The previously published offline archive predates this case and does not already contain it.

Strategy I compares the two 40-Mb extracts, retains 100–500-bp InDels and designs reference-based primers. Strategy III selects complete annotated B73 protein-coding genes with at least 200 bp of margin at either regional boundary, uses the longest CDS per gene, and extracts CDS-guided genomic segments with 200-bp flanks for alignment. **The 20-Mb value defines the input region; it is not a 20-Mb flank applied separately to every gene.**

`check_outputs.py` checks hashes, VCF reference alleles, reference-template primer coordinates, equal-length alignment records and THP3 sequence coverage. `genome_coverage.tsv` explicitly lists missing parental sequences; absence of a mapped sequence within this extract does not establish biological gene absence.

## Strategy II: full-depth regional BAMs

**Lite means a smaller genomic interval, with the original alignment depth retained.** The two BAMs contain every source alignment overlapping **B73 NAM-5.0 chr3:180,402,538–220,408,213**, the THP3 gene plus 20 Mb upstream and downstream (40,005,676 bp). There is no random subsampling, MAPQ threshold, alignment-flag filter or duplicate removal during extraction. Reads overlapping a boundary are retained as complete alignment records; mates outside the interval are not independently fetched. These are experimental sequencing data, not simulated reads.

| Parent | Material | Original BAM sample name | Delivered input |
| --- | --- | --- | --- |
| p1 | B73 | B73 | `data/B73.THP3_flanks20Mb.bam` and `.bam.bai` |
| p2 | Ames21814, identity confirmed by the authors | Teosinte | `data/Ames21814.THP3_flanks20Mb.bam` and `.bam.bai` |

Both BAMs use the **B73 reference coordinate system**, including the Ames21814 reads. Their coordinates must not be confused with positions in the Ames21814 assembly used by Strategies I and III. The original genome-wide alignments and duplicate flags are preserved. B73 contains 8,903,524 alignment records and Ames21814 contains 7,858,180; mean aligned-base depth across all 40,005,676 positions is 33.2769× and 26.1815×, respectively, including duplicate, secondary and QC-failed records under the documented depth command. These depth values are identical before and after extraction. Private command-line/path metadata are removed from headers; sequence dictionaries, read groups and alignment records remain unchanged. `sources/regional_bam_provenance.json` records the before/after SAM-record digests, counts, per-position depth digests, file sizes and SHA-256 checksums. All of those comparisons must pass before distributing the BAMs.

The supplied **complete B73 chr3 reference** retains its original `chr3` identifier and length of 238,017,767 bp; it is not the renamed 40-Mb assembly extract. Thus marker positions from Strategy II are already original B73 NAM-5.0 coordinates. Only the THP3 ±20-Mb interval is called. The example divides this interval into nonoverlapping 1-Mb chunks, calls them independently with identical settings, then gathers each parent's gVCFs in coordinate order before joint genotyping. Chunking partitions the computation; it does not sample reads or remove any part of the target interval. Up to 32 allocated CPU cores are suggested for this full-depth case; fewer cores run fewer chunks concurrently. Reference/BAM dictionary validation remains enabled. The BAM input path skips FASTQ cleaning, remapping and duplicate marking. GATK HaplotypeCaller is configured with `--max-reads-per-alignment-start 0` to disable its alignment-start downsampling, as specified in the [official GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/27007962724507-HaplotypeCaller). Standard caller quality/read filters still apply; preservation of all records in the input BAM does not mean that every read supports a reported genotype.

```bash
bash scripts/run_thp3_examples.sh 2 32
# Run all three strategies:
bash scripts/run_thp3_examples.sh all 32
```

The complete THP3 case archive supplies BAMs and indexes directly in `data/`. Large reference FASTAs are restored losslessly from `data_payload/` by the runner. The default source-code ZIP omits the formal input datasets and reference result tables; install them from Zenodo or use the complete case archive for execution. All input data and an installed runtime must be present before taking the case to an offline machine. No sequencing-data download is performed during the analysis.

The workflow retains 20–100-bp biallelic InDels with contrasting `0/0` and `1/1` calls, then designs reference-based primers. The checker verifies input hashes, indexed BAM record counts, reference alleles, contrasting genotypes and primer-template matches. This regional case does not estimate unbiased genome-wide calling performance; PCR and genome-wide primer specificity remain separate validation steps.

## Observed local results

All three expanded cases completed on native Linux x86_64 / glibc 2.28 using the published offline runtime, six allocated CPU cores and 24 GB requested memory. The combined job elapsed time was 1 h 29 min 24 s, including runtime installation and input preparation. See `native_linux_validation_20261019.json`. These allocation settings do not establish a minimum hardware requirement. The assembly-based cases also passed on macOS arm64 with four threads; the Strategy I primer sequences and coordinates were identical across platforms. The new reports are separate from the original six-case Linux offline validation.

- Strategy I: 155 filtered InDel loci and 465 reference-template-verified primer pairs.
- A further exact-binding screen within the two 40-Mb extracts identified 138 pairs at 58 loci with one predicted 50–2,000-bp amplicon in each parent and different product lengths. These are regional computational predictions, not full-genome specificity checks or PCR validation.
- Strategy III: 1,097 nonempty gene alignments. All contain B73 sequence; 916 also contain Ames21814 sequence, including THP3. The remaining 181 have no extracted Ames21814 sequence within the supplied interval. Do not interpret those missing hits as gene deletions.
- Strategy II: 2,487 contrasting 20–100-bp InDels and 7,443 reference-template-verified primer pairs. Primer3 returned three pairs at each of 2,481 loci; six loci yielded no pairs under the supplied settings. These are computational designs, not PCR success counts.

Selected reference outputs are under `reference_results/`. The author's separately reported 36.55% experimental primer success rate concerns the 10-Mb THP3 laboratory test and is not derived from these computational counts.

## Regenerate regional inputs

If the complete source assemblies and B73 annotation are available:

```bash
python examples/thp3_20Mb/rebuild_inputs.py \
  --b73-genome /path/to/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
  --ames21814-genome /path/to/Zv-Ames21814-REFERENCE-CEMPS-1.0.fa.gz \
  --b73-gff3 /path/to/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz \
  --output /path/to/regenerated_inputs
```

Compare regenerated file SHA-256 values with `provenance.json` before using them. No synthetic reads are created by this script.

To regenerate the read-based inputs where the original indexed BAMs are available:

```bash
python examples/thp3_20Mb/extract_regional_bams.py \
  --b73 /path/to/B73.sorted.mkdup.bam \
  --ames21814 /path/to/Teosinte.sorted.mkdup.bam \
  --region chr3:180402538-220408213 \
  --outdir /path/to/new_empty_output_directory
```

The export script opens the originals read-only, creates indexes for the new BAMs, and checks exact alignment-record and per-position depth preservation. It never reduces sequencing depth.
