#!/usr/bin/env bash
# Create the workflow environment while ignoring broken global .condarc channels.
# Conda retains downloaded packages in its persistent local package cache.
set -euo pipefail

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "ERROR: the complete supplied environment targets Linux x86_64." >&2
    echo "Use a Linux x86_64 host; do not install Linux binaries into a native macOS environment." >&2
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found in PATH." >&2
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
environment_name=${1:-indelmarker}
# By default use Conda's persistent local cache.  This avoids duplicating large
# platform-specific archives into the Git checkout.  Set
# INDELMARKER_CONDA_PKGS_DIRS to use a project- or site-specific cache instead.
conda_base=$(conda info --base)
package_cache=${INDELMARKER_CONDA_PKGS_DIRS:-"$conda_base/pkgs"}
offline_mode=${INDELMARKER_OFFLINE:-0}

mkdir -p "$package_cache"

conda_args=(
    create --yes --name "$environment_name"
    --override-channels
    --strict-channel-priority
    --channel conda-forge
    --channel bioconda
    --file "$project_dir/envs/indelmarker-linux-64.txt"
)

if [[ "$offline_mode" == "1" ]]; then
    conda_args+=(--offline)
fi

CONDA_PKGS_DIRS="$package_cache" conda "${conda_args[@]}"

echo "Environment '$environment_name' created. Activate it with: conda activate $environment_name"
echo "Package cache: $package_cache"
