# Offline installation: Linux x86_64, glibc 2.28 or newer

Use the complete `InDelMarker-offline-linux-x86_64-glibc2.28.tar.gz` release bundle and its SHA-256 file. The ordinary GitHub source ZIP contains application code and examples; it does not contain the binary runtime.

Download on a computer with Internet access, then copy both files to the offline Linux computer. CentOS/RHEL/Rocky Linux 8 is the build compatibility target. This bundle is not a native macOS, Windows, ARM64, or glibc-2.17 runtime.

```bash
sha256sum -c InDelMarker-offline-linux-x86_64-glibc2.28.tar.gz.sha256
tar -xzf InDelMarker-offline-linux-x86_64-glibc2.28.tar.gz
cd InDelMarker-offline-linux-x86_64-glibc2.28
bash scripts/install_offline_linux.sh
bash scripts/run_offline_examples.sh --suite all --cores 2
```

The offline computer needs Bash, GNU core utilities, tar, and gzip. It needs no preinstalled Conda, Python, dependency solver, administrator privileges, or network connection. Installation verifies the runtime checksum, extracts the complete environment, repairs prefixes with the bundled Python, and checks the workflow imports and executables. The runner automatically selects this runtime without changing shell startup files.

The default installation directory is `.runtime` inside the extracted bundle. To select another empty directory:

```bash
bash scripts/install_offline_linux.sh offline/linux-64/indelmarker-runtime.tar.gz /path/to/indelmarker-runtime
INDELMARKER_RUNTIME_PREFIX=/path/to/indelmarker-runtime bash scripts/run_offline_examples.sh --suite all --cores 2
```

Keep an initialized runtime at its installation path. If it must move, install the original archive again into a new empty directory. Repeating installation at the same path and checksum retains a completed installation. Existing unrelated files are never deleted.

For your own data, copy `config.yaml` to `config.local.yaml`, edit the input and output paths, and run the bundled workflow through the same runtime:

```bash
bash scripts/run_offline.sh snakemake --snakefile Snakefile --configfile config.local.yaml --cores 16 --rerun-incomplete
```

The command runner also supports individual bundled programs, for example `bash scripts/run_offline.sh samtools --version`.

Successful execution ends with `All requested example strategies completed successfully.` Each strategy writes its own `validation.json`; empty VCFs, empty primer workbooks, and incomplete gene alignments fail validation. The maize Strategy II reads are simulated and are labelled accordingly.

`envs/indelmarker-linux-64.lock` records exact package URLs/builds and MD5 hashes. `software/packages.json` records versions, builds, upstream URLs, available SHA-256 hashes, dependencies, and license identifiers. Upstream notices and recipes are retained under `software/upstream-notices` when included in the original packages. The runtime uses `snakemake-minimal`, which supplies the same workflow executable without optional cloud storage integrations.

## Maintainer build and verification

Build with `envs/offline-build.Dockerfile`, which pins a Rocky Linux 8.10 image and verifies the Miniforge installer checksum. Run `scripts/build_offline_bundle.sh` inside that Linux x86_64 builder with a new runtime prefix and an output directory. Before release, install the resulting bundle in a new Linux container with `--network none`, no Conda installation or package cache, and a different installation path; run all six example workflows there. Retain the output-content checks with the release.
