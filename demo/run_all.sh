#!/usr/bin/env bash
set -euo pipefail
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
exec bash "$project_dir/scripts/run_examples.sh" --suite demo "$@"
