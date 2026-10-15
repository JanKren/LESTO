#!/usr/bin/env bash
# Build the case-ladder discussion report with Tectonic.  Figures and PSI
# assets are read from ../case-ladder (run its build.sh --refresh first to
# redraw them).  TECTONIC_BIN overrides the engine.
set -euo pipefail
cd -- "$(dirname -- "$0")"
engine=${TECTONIC_BIN:-$(command -v tectonic || echo /home/jan/.local/opt/tectonic-0.17.0/tectonic)}
[[ -x "$engine" ]] || { printf 'Install Tectonic or set TECTONIC_BIN.\n' >&2; exit 1; }
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}"
mkdir -p build
"$engine" --keep-logs --outdir build main.tex > build/compile.stdout 2>&1 || { tail -30 build/compile.stdout >&2; exit 1; }
cp build/main.pdf lesto-case-ladder-explained.pdf
printf 'Built %s/lesto-case-ladder-explained.pdf\n' "$PWD"
