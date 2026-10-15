#!/usr/bin/env bash
# Build the Physics of Fluids (REVTeX 4-2, aip/pof) version of the manuscript.
# Usage: bash manuscript/pof/build.sh            # two-column reprint
#        bash manuscript/pof/build.sh --preprint # single-column review copy
set -euo pipefail
pof_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$pof_dir"
mode=reprint
if [[ "${1:-}" == "--preprint" ]]; then
    mode=preprint
elif [[ -n "${1:-}" ]]; then
    printf 'Usage: bash manuscript/pof/build.sh [--preprint]\n' >&2
    exit 2
fi
if [[ -n "${TECTONIC_BIN:-}" ]]; then
    compiler="$TECTONIC_BIN"
elif command -v tectonic >/dev/null 2>&1; then
    compiler="$(command -v tectonic)"
elif [[ -x /home/jan/.local/opt/tectonic-0.17.0/tectonic ]]; then
    compiler=/home/jan/.local/opt/tectonic-0.17.0/tectonic
else
    printf 'Install Tectonic or compile lesto-pof.tex with latexmk -pdf.\n' >&2
    exit 1
fi
src=lesto-pof.tex
if [[ "$mode" == preprint ]]; then
    src=lesto-pof-preprint.tex
    sed 's/\[aip,pof,reprint,/[aip,pof,preprint,/' lesto-pof.tex > "$src"
fi
mkdir -p build
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
    "$compiler" --keep-logs --keep-intermediates --outdir build "$src"
out="${src%.tex}.pdf"
cp "build/$out" "$out"
[[ "$mode" == preprint ]] && rm -f "$src"
log="build/${src%.tex}.log"
printf 'Built %s/%s\n' "$pof_dir" "$out"
printf 'Overfull boxes: %s; undefined references/citations: %s\n' \
    "$(grep -c '^Overfull' "$log" || true)" \
    "$(grep -ciE 'undefined' "$log" || true)"
# Supplementary material: the shared ../supplementary.tex, built from the
# manuscript root so that its generated/ and references.bib paths resolve.
if [[ "$mode" == reprint ]]; then
    (cd .. && XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
        "$compiler" --keep-logs --keep-intermediates --outdir pof/build supplementary.tex)
    cp build/supplementary.pdf lesto-pof-supplementary.pdf
    slog=build/supplementary.log
    printf 'Built %s/lesto-pof-supplementary.pdf\n' "$pof_dir"
    printf 'Supplement overfull boxes: %s; undefined references/citations: %s\n' \
        "$(grep -c '^Overfull' "$slog" || true)" \
        "$(grep -ciE 'undefined' "$slog" || true)"
fi
