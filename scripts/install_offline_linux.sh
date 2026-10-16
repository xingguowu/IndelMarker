#!/usr/bin/env bash
# Install the bundled Linux environment without contacting Conda channels.
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
package_cache="$project_dir/offline/linux-64/conda-pkgs"

if [[ ! -d "$package_cache" ]] || ! find "$package_cache" -maxdepth 1 -type f \( -name '*.conda' -o -name '*.tar.bz2' \) -print -quit | grep -q .; then
    echo "ERROR: Offline package cache is missing: $package_cache" >&2
    exit 1
fi

INDELMARKER_CONDA_PKGS_DIRS="$package_cache" \
INDELMARKER_OFFLINE=1 \
    "$script_dir/create_conda_env.sh" "${1:-indelmarker}"
