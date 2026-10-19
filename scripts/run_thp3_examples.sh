#!/usr/bin/env bash
set -euo pipefail
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$project_dir"
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:--Xmx2g}"
python examples/thp3_20Mb/prepare_bundled_inputs.py
strategy=${1:-all}
cores=${2:-8}
[[ "$cores" =~ ^[1-9][0-9]*$ ]] || { echo 'Core count must be a positive integer.' >&2; exit 2; }
case "$strategy" in
  all) strategies=(1 2 3);;
  assemblies) strategies=(1 3);;
  1|3) strategies=("$strategy");;
  2)
    strategies=(2);;
  *) echo 'Usage: bash scripts/run_thp3_examples.sh [all|assemblies|1|2|3] [cores]' >&2; exit 2;;
esac
for number in "${strategies[@]}"; do
  if [[ "$number" == 2 ]]; then
    for f in examples/thp3_20Mb/data/B73.chr3.fa examples/thp3_20Mb/data/{B73,Ames21814}.THP3_flanks20Mb.bam{,.bai}; do
      [[ -s "$f" ]] || { echo "Missing: $f. See examples/thp3_20Mb/README.md." >&2; exit 1; }
    done
  fi
  snakemake --snakefile Snakefile --configfile "examples/thp3_20Mb/strategy${number}/config.yaml" --cores "$cores" --config "threads=$cores" --printshellcmds --rerun-incomplete --latency-wait 10
  python examples/thp3_20Mb/check_outputs.py --strategy "$number"
done
