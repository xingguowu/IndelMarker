#!/usr/bin/env bash
# Maintainer regression test. The end user does not need Docker.
set -euo pipefail
archive=${1:?Usage: verify_offline_bundle.sh BUNDLE_TAR_GZ REPORT_DIRECTORY}
report_dir=${2:?Usage: verify_offline_bundle.sh BUNDLE_TAR_GZ REPORT_DIRECTORY}
archive_dir=$(CDPATH= cd -- "$(dirname -- "$archive")" && pwd)
archive_name=$(basename "$archive")
[[ -f "$archive_dir/$archive_name.sha256" ]] || { echo "Bundle checksum missing." >&2; exit 1; }
mkdir -p "$report_dir"
report_dir=$(CDPATH= cd -- "$report_dir" && pwd)
docker run --interactive --rm --network none --platform linux/amd64 --cpus 2 --memory 6g \
    --mount "type=bind,source=$archive_dir,target=/payload,readonly" \
    --mount "type=bind,source=$report_dir,target=/validation" \
    rockylinux/rockylinux:8.10@sha256:e8a49c5403b687db05d4d67333fa45808fbe74f36e683cec7abb1f7d0f2338c6 \
    bash -s -- "$archive_name" <<'VERIFY'
set -euo pipefail
archive_name=$1
test "$(getconf GNU_LIBC_VERSION)" = 'glibc 2.28'
test "$(uname -m)" = x86_64
# Some Docker kernels expose unused tunnel interfaces even with --network none.
# Require an empty IPv4 routing table in this isolated network namespace.
test -z "$(awk 'NR > 1 {print}' /proc/net/route)"
if command -v conda >/dev/null || command -v python >/dev/null || command -v python3 >/dev/null; then
    echo "Unexpected preinstalled Conda/Python on PATH." >&2
    exit 1
fi
cd /payload
sha256sum -c "$archive_name.sha256"
mkdir -p '/tmp/IndelMarker offline verification'
tar -xzf "$archive_name" -C '/tmp/IndelMarker offline verification'
cd '/tmp/IndelMarker offline verification/InDelMarker-offline-linux-x86_64-glibc2.28'
printf 'corrupt archive\n' > /tmp/indelmarker-corrupt.tar.gz
cp offline/linux-64/indelmarker-runtime.tar.gz.sha256 /tmp/indelmarker-corrupt.tar.gz.sha256
if bash scripts/install_offline_linux.sh /tmp/indelmarker-corrupt.tar.gz /tmp/invalid-runtime > /validation/checksum-rejection.log 2>&1; then
    echo "Corrupted runtime was not rejected." >&2
    exit 1
fi
test ! -e /tmp/invalid-runtime
bash scripts/install_offline_linux.sh | tee /validation/installation.log
bash scripts/install_offline_linux.sh | tee /validation/repeated-installation.log
bash scripts/run_offline_examples.sh --suite all --cores 2 2>&1 | tee /validation/examples.log
.runtime/bin/python -I - "$archive_name" <<'PY'
import json
import platform
import sys
from pathlib import Path
root = Path.cwd()
results = []
for suite, directory in [('demo', 'demo/work'), ('maize', 'examples/maize/results')]:
    for strategy in ['1', '2', '3']:
        result = json.loads((root / directory / ('strategy' + strategy) / 'validation.json').read_text())
        assert result['status'] == 'PASS', result
        results.append(result)
report = {
    'status': 'PASS', 'network_mode': 'none',
    'network_interfaces': [line.split(':', 1)[0].strip() for line in Path('/proc/net/dev').read_text().splitlines()[2:]],
    'ipv4_routes': Path('/proc/net/route').read_text().splitlines()[1:],
    'archive_sha256': (Path('/payload') / (sys.argv[1] + '.sha256')).read_text().split()[0],
    'preinstalled_conda_on_path': False, 'preinstalled_python_on_path': False,
    'glibc': platform.libc_ver(), 'machine': platform.machine(),
    'installation_path_contains_spaces': True, 'corrupted_runtime_rejected': True,
    'repeated_installation_passed': True, 'results': results,
}
Path('/validation/offline-validation.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
PY
VERIFY
