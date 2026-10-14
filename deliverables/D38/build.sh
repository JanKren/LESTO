#!/usr/bin/env bash
# Build LESTO deliverable D38 with Tectonic.  Figures, generated tables and
# the bibliography are read from ../../manuscript (not copied).
set -euo pipefail
d38_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$d38_dir"
compiler="${TECTONIC_BIN:-/home/jan/.local/opt/tectonic-0.17.0/tectonic}"
if [[ ! -x "$compiler" ]]; then
    if command -v tectonic >/dev/null 2>&1; then
        compiler="$(command -v tectonic)"
    else
        printf 'Tectonic not found; set TECTONIC_BIN.\n' >&2
        exit 1
    fi
fi
mkdir -p build
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
    "$compiler" --keep-logs --keep-intermediates --outdir build main.tex
cp build/main.pdf LESTO-D38-GEMS-CFD-coupling.pdf
printf 'Built %s/LESTO-D38-GEMS-CFD-coupling.pdf\n' "$d38_dir"
