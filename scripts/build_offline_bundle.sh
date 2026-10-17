#!/usr/bin/env bash
# Maintainer build: run in the supplied Linux x86_64/glibc-2.28 container.
set -euo pipefail
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
runtime_prefix=${1:?Usage: build_offline_bundle.sh EMPTY_RUNTIME_PREFIX OUTPUT_DIRECTORY}
output_dir=${2:?Usage: build_offline_bundle.sh EMPTY_RUNTIME_PREFIX OUTPUT_DIRECTORY}
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || { echo "Build requires Linux x86_64." >&2; exit 1; }
[[ "$(getconf GNU_LIBC_VERSION)" == 'glibc 2.28' ]] || { echo "Use the supplied glibc-2.28 builder for the CentOS 8 bundle." >&2; exit 1; }
[[ ! -e "$runtime_prefix" ]] || { echo "Runtime build prefix already exists; choose a new empty path." >&2; exit 1; }
mkdir -p "$output_dir"
output_dir=$(CDPATH= cd -- "$output_dir" && pwd)
bundle_name=InDelMarker-offline-linux-x86_64-glibc2.28
stage="$output_dir/$bundle_name"
[[ ! -e "$stage" ]] || { echo "Bundle staging directory already exists: $stage" >&2; exit 1; }
specification="$project_dir/envs/indelmarker-offline-linux-64.txt"
[[ ! -f "$project_dir/envs/indelmarker-linux-64.lock" ]] || specification="$project_dir/envs/indelmarker-linux-64.lock"
conda create --yes --prefix "$runtime_prefix" --solver libmamba \
    --override-channels --strict-channel-priority \
    --channel https://conda.anaconda.org/conda-forge \
    --channel https://conda.anaconda.org/bioconda \
    --file "$specification"
python "$script_dir/stage_offline_bundle.py" "$project_dir" "$runtime_prefix" "$stage"
conda list --prefix "$runtime_prefix" --explicit --md5 > "$stage/envs/indelmarker-linux-64.lock"
conda-pack --prefix "$runtime_prefix" --output "$stage/offline/linux-64/indelmarker-runtime.tar.gz" --n-threads 2
(
    cd "$stage/offline/linux-64"
    sha256sum indelmarker-runtime.tar.gz > indelmarker-runtime.tar.gz.sha256
)
tar -czf "$output_dir/$bundle_name.tar.gz" -C "$output_dir" "$bundle_name"
(
    cd "$output_dir"
    sha256sum "$bundle_name.tar.gz" > "$bundle_name.tar.gz.sha256"
)
echo "Built offline bundle: $output_dir/$bundle_name.tar.gz"
echo "Validate in a fresh container with --network none before publishing."
