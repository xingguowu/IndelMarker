#!/usr/bin/env bash
set -euo pipefail
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec bash "$script_dir/run_offline.sh" bash "$script_dir/run_examples.sh" "$@"
