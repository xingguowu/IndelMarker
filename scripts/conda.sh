#!/usr/bin/env bash
# Backwards-compatible entry point for the supported installer.
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec bash "$script_dir/create_conda_env.sh" "$@"
