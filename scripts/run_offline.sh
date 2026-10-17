#!/usr/bin/env bash
set -euo pipefail
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(dirname "$script_dir")
prefix=${INDELMARKER_RUNTIME_PREFIX:-"$project_dir/.runtime"}
if [[ ! -f "$prefix/.indelmarker-ready" || ! -f "$prefix/.indelmarker-prefix" ]]; then
    echo "ERROR: install the offline runtime first with bash scripts/install_offline_linux.sh." >&2
    exit 1
fi
prefix=$(CDPATH= cd -- "$prefix" && pwd)
if [[ "$(cat "$prefix/.indelmarker-prefix")" != "$prefix" ]]; then
    echo "ERROR: the initialized runtime was moved. Install the original archive in a new empty directory." >&2
    exit 1
fi
export PATH="$prefix/bin:$PATH"
export CONDA_PREFIX="$prefix"
export JAVA_HOME="$prefix/lib/jvm"
export JAVA_LD_LIBRARY_PATH="$JAVA_HOME/lib/server"
# MAFFT's packaged default prefix assignment is unquoted after relocation.
# Its supported environment override preserves paths containing spaces.
export MAFFT_BINARIES="$prefix/libexec/mafft"
export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH
if (( $# == 0 )); then
    echo "Usage: bash scripts/run_offline.sh COMMAND [ARGUMENT ...]" >&2
    exit 2
fi
exec "$@"
