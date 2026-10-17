#!/usr/bin/env bash
# Install the prebuilt runtime without Conda, Python, a solver, or the network.
set -euo pipefail
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
archive=${1:-"$project_dir/offline/linux-64/indelmarker-runtime.tar.gz"}
prefix=${2:-"$project_dir/.runtime"}
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "ERROR: this offline runtime requires Linux x86_64." >&2
    exit 1
fi
glibc=$(getconf GNU_LIBC_VERSION 2>/dev/null) || {
    echo "ERROR: a glibc-based Linux system is required." >&2; exit 1;
}
glibc=${glibc##* }
IFS=. read -r glibc_major glibc_minor _ <<< "$glibc"
if (( glibc_major < 2 || (glibc_major == 2 && glibc_minor < 28) )); then
    echo "ERROR: glibc >= 2.28 is required; detected $glibc." >&2
    exit 1
fi
for tool in tar gzip sha256sum; do
    command -v "$tool" >/dev/null || { echo "ERROR: missing system utility: $tool" >&2; exit 1; }
done
if [[ ! -f "$archive" || ! -f "$archive.sha256" ]]; then
    echo "ERROR: offline runtime or its checksum is missing: $archive" >&2
    echo "Download the complete offline release bundle; the source ZIP alone contains no binary runtime." >&2
    exit 1
fi
archive_dir=$(CDPATH= cd -- "$(dirname -- "$archive")" && pwd)
archive_name=$(basename "$archive")
read -r expected _ < "$archive.sha256"
[[ "$expected" =~ ^[0-9a-f]{64}$ ]] || { echo "ERROR: invalid runtime checksum." >&2; exit 1; }
actual=$(sha256sum < "$archive")
actual=${actual%% *}
[[ "$actual" == "$expected" ]] || { echo "ERROR: runtime checksum mismatch; copy or download the archive again." >&2; exit 1; }
mkdir -p "$(dirname "$prefix")"
prefix_parent=$(CDPATH= cd -- "$(dirname -- "$prefix")" && pwd)
prefix="$prefix_parent/$(basename "$prefix")"
if [[ -d "$prefix" ]] && [[ -n "$(ls -A "$prefix")" ]]; then
    if [[ -f "$prefix/.indelmarker-ready" && -f "$prefix/.indelmarker-prefix" ]] &&
       [[ "$(cat "$prefix/.indelmarker-ready")" == "$expected" && "$(cat "$prefix/.indelmarker-prefix")" == "$prefix" ]]; then
        echo "Offline runtime is already installed: $prefix"
        exit 0
    fi
    echo "ERROR: installation directory is not empty: $prefix" >&2
    echo "Choose another empty directory as the second argument. Existing files are retained." >&2
    exit 1
fi
mkdir -p "$prefix"
tar -xzf "$archive_dir/$archive_name" -C "$prefix"
"$prefix/bin/python" -I "$prefix/bin/conda-unpack"
export PATH="$prefix/bin:$PATH"
export CONDA_PREFIX="$prefix"
export JAVA_HOME="$prefix/lib/jvm"
export JAVA_LD_LIBRARY_PATH="$JAVA_HOME/lib/server"
export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH
"$prefix/bin/python" -I -c 'import sys, yaml, snakemake, pulp, numpy, pysam; assert sys.version_info[:2] == (3, 10); print("Bundled Python and workflow imports: PASS")'
for tool in snakemake samtools bcftools minimap2 syri primer3_core bwa fastp gatk java bedtools mafft exonerate; do
    [[ -x "$prefix/bin/$tool" ]] || { echo "ERROR: runtime executable missing: $tool" >&2; exit 1; }
done
printf '%s\n' "$prefix" > "$prefix/.indelmarker-prefix"
printf '%s\n' "$expected" > "$prefix/.indelmarker-ready"
echo "Offline runtime installed successfully: $prefix"
printf 'INDELMARKER_RUNTIME_PREFIX=%q bash %q --suite all --cores 2\n' "$prefix" "$script_dir/run_offline_examples.sh"
