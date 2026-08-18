#!/usr/bin/env bash
# Create the workflow environment while ignoring broken global .condarc channels.
set -euo pipefail

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found in PATH." >&2
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
environment_name=${1:-indelmarker}

conda create --yes --name "$environment_name" \
    --override-channels \
    --strict-channel-priority \
    --channel conda-forge \
    --channel bioconda \
    --file "$project_dir/envs/indelmarker-linux-64.txt"

echo "Environment '$environment_name' created. Activate it with: conda activate $environment_name"
