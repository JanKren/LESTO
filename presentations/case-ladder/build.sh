#!/usr/bin/env bash
# Build the case-ladder deck with Tectonic.  Figures are redrawn only when
# asked (they need runs/ and the tracked repository files).
#
#   bash build.sh              # main.tex -> lesto-case-ladder.pdf
#   bash build.sh --refresh    # make_figures.py first
#
# Environment: TECTONIC_BIN (default the workstation's Tectonic 0.17.0),
# XDG_CACHE_HOME (default /tmp/lesto-beamer-cache), LESTO_CPUS (taskset
# list for the build, default 4,5,11), LESTO_BEAMER_OFFLINE=1 (cache only).
set -euo pipefail
cd -- "$(dirname -- "$0")"
PAGES=32
mkdir -p build
cpus=${LESTO_CPUS:-4,5,11}
if [[ "${1:-}" == "--refresh" ]]; then
    MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/lesto-matplotlib}" taskset -c "$cpus" python3 make_figures.py
elif [[ -n "${1:-}" ]]; then
    printf 'Usage: bash build.sh [--refresh]\n' >&2
    exit 2
fi
engine=${TECTONIC_BIN:-/home/jan/.local/opt/tectonic-0.17.0/tectonic}
[[ -x "$engine" ]] || engine=$(command -v tectonic || true)
if [[ -z "$engine" ]]; then
    printf 'Install Tectonic or set TECTONIC_BIN.\n' >&2
    exit 1
fi
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}"
flags=(--keep-logs --keep-intermediates --outdir build)
if [[ "${LESTO_BEAMER_OFFLINE:-0}" == 1 ]]; then flags+=(--only-cached); fi
taskset -c "$cpus" nice -n 10 "$engine" "${flags[@]}" main.tex > build/compile.stdout 2>&1 || {
    tail -40 build/compile.stdout >&2
    exit 1
}
cp build/main.pdf lesto-case-ladder.pdf
overfull=$(grep -c 'Overfull' build/main.log || true)
pages=$(pdfinfo lesto-case-ladder.pdf | awk '/^Pages:/{print $2}')
if [[ "$pages" != "$PAGES" ]]; then
    printf 'Expected %s pages, found %s.\n' "$PAGES" "$pages" >&2
    exit 1
fi
printf 'Built %s/lesto-case-ladder.pdf (%s pages, %s overfull boxes).\n' "$PWD" "$pages" "$overfull"
