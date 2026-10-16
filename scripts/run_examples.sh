#!/usr/bin/env bash
# Run directly from any directory, including ZIP downloads without a .git folder.
set -euo pipefail
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$project_dir"
suite=demo
strategy=all
cores=2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --suite|--strategy|--cores)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --suite) suite=$2;;
                --strategy) strategy=$2;;
                --cores) cores=$2;;
            esac
            shift 2;;
        -h|--help)
            echo "Usage: bash scripts/run_examples.sh [--suite demo|maize|all] [--strategy 1|2|3|all] [--cores N]"
            exit 0;;
        *) echo "Unknown option: $1" >&2; exit 2;;
    esac
done
[[ "$suite" =~ ^(demo|maize|all)$ ]] || { echo "Invalid suite: $suite" >&2; exit 2; }
[[ "$strategy" =~ ^(1|2|3|all)$ ]] || { echo "Invalid strategy: $strategy" >&2; exit 2; }
[[ "$cores" =~ ^[1-9][0-9]*$ ]] || { echo "--cores must be a positive integer" >&2; exit 2; }
# Bound Java heap for the small examples; retain user-supplied JVM settings.
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:--Xmx2g}"
for program in python snakemake; do
    command -v "$program" >/dev/null || { echo "Missing $program. Activate the indelmarker Conda environment." >&2; exit 1; }
done
suites=("$suite")
[[ "$suite" != all ]] || suites=(demo maize)
strategies=("$strategy")
[[ "$strategy" != all ]] || strategies=(1 2 3)
for current in "${suites[@]}"; do
    if [[ "$current" == demo ]]; then
        directory=demo
        python demo/generate_demo_inputs.py
    else
        directory=examples/maize
    fi
    for number in "${strategies[@]}"; do
        python scripts/check_examples.py --suite "$current" --strategy "$number" --inputs-only
        echo "Running $current Strategy $number with $cores core(s)."
        snakemake --snakefile Snakefile --configfile "$directory/strategy${number}/config.yaml" \
            --cores "$cores" --config "threads=$cores" --printshellcmds --rerun-incomplete --latency-wait 10
        python scripts/check_examples.py --suite "$current" --strategy "$number"
    done
done
echo "All requested example strategies completed successfully."
