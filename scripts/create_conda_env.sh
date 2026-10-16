#!/usr/bin/env bash
# Create the workflow environment while ignoring broken global .condarc channels.
# Conda retains downloaded packages in its persistent local package cache.
set -euo pipefail

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "ERROR: the complete supplied environment targets Linux x86_64." >&2
    echo "Use a Linux x86_64 host; do not install Linux binaries into a native macOS environment." >&2
    exit 1
fi

# Cluster login shells may not have initialized an existing Conda installation.
# Honour an explicit executable before searching PATH and common user locations.
conda_command=${INDELMARKER_CONDA_EXE:-}
if [[ -n "$conda_command" ]]; then
    if [[ ! -x "$conda_command" ]]; then
        echo "ERROR: INDELMARKER_CONDA_EXE is not executable: $conda_command" >&2
        exit 1
    fi
elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
    # Keep the active distribution when a different system Conda is on PATH.
    conda_command=$CONDA_EXE
elif command -v conda >/dev/null 2>&1; then
    conda_command=conda
else
    conda_candidates=(
        "${CONDA_PREFIX:-}/bin/conda"
        "$HOME/miniforge3/bin/conda"
        "$HOME/miniconda3/bin/conda"
        "$HOME/anaconda3/bin/conda"
        "$HOME/mambaforge/bin/conda"
        /opt/conda/bin/conda
    )
    for candidate in "${conda_candidates[@]}"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            conda_command=$candidate
            break
        fi
    done
fi

if [[ -z "$conda_command" ]]; then
    cat >&2 <<'HELP'
ERROR: Conda was not found in PATH, CONDA_EXE, or common installation locations.
On a cluster, first check the site's environment modules:
    module avail
Load the Conda/Anaconda/Miniforge module shown by your site, then retry.
For an existing installation in another location, use its actual path:
    INDELMARKER_CONDA_EXE=/path/to/bin/conda bash scripts/create_conda_env.sh
If Conda is not installed, follow the official Miniforge instructions:
    https://github.com/conda-forge/miniforge#install
The installer does not require sudo. Do not overwrite an existing installation.
HELP
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
environment_name=${1:-indelmarker}
# By default use Conda's persistent local cache.  This avoids duplicating large
# platform-specific archives into the Git checkout.  Set
# INDELMARKER_CONDA_PKGS_DIRS to use a project- or site-specific cache instead.
conda_base=$("$conda_command" info --base)
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

CONDA_PKGS_DIRS="$package_cache" "$conda_command" "${conda_args[@]}"

echo "Environment '$environment_name' created. Run these commands in your current terminal:"
printf 'source %q\n' "$conda_base/etc/profile.d/conda.sh"
printf 'conda activate %q\n' "$environment_name"
echo "Package cache: $package_cache"
