#!/usr/bin/env bash
set -euo pipefail
manuscript_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$manuscript_dir"
if [[ "${1:-}" == "--figures" ]]; then
    MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/lesto-manuscript-mpl}" python3 prepare_results.py
elif [[ -n "${1:-}" ]]; then
    printf 'Usage: bash manuscript/build.sh [--figures]\n' >&2
    exit 2
fi
mkdir -p build
if [[ -n "${TECTONIC_BIN:-}" ]]; then
    compiler="$TECTONIC_BIN"
elif command -v tectonic >/dev/null 2>&1; then
    compiler="$(command -v tectonic)"
elif [[ -x /home/jan/.local/opt/tectonic-0.17.0/tectonic ]]; then
    compiler=/home/jan/.local/opt/tectonic-0.17.0/tectonic
else
    printf 'Install Tectonic or compile main.tex with latexmk -pdf.\n' >&2
    exit 1
fi
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
    "$compiler" --keep-logs --keep-intermediates --outdir build main.tex
cp build/main.pdf lesto-iodine-methods.pdf
printf 'Built %s/lesto-iodine-methods.pdf\n' "$manuscript_dir"
