#!/usr/bin/env bash
# Build the reviewed, frozen deck. Refresh plots/status only when requested.
set -euo pipefail
cd -- "$(dirname -- "$0")"
mkdir -p build
if [[ "${1:-}" == "--refresh" ]]; then
    MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/lesto-matplotlib}" python3 prepare_figures.py
elif [[ -n "${1:-}" ]]; then
    printf 'Usage: bash build.sh [--refresh]\n' >&2
    exit 2
fi

if command -v latexmk >/dev/null 2>&1; then
    latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex > build/compile.stdout 2>&1
elif command -v pdflatex >/dev/null 2>&1; then
    pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > build/compile.stdout 2>&1
    pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex >> build/compile.stdout 2>&1
else
    deck_engine="${TECTONIC_BIN:-}"
    if [[ -z "$deck_engine" ]] && command -v tectonic >/dev/null 2>&1; then
        deck_engine="$(command -v tectonic)"
    fi
    if [[ -z "$deck_engine" && -x /tmp/lesto-beamer-tools/tectonic ]]; then
        deck_engine=/tmp/lesto-beamer-tools/tectonic
        if [[ -d /tmp/lesto-beamer-cache ]]; then
            export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}"
        fi
    fi
    if [[ -z "$deck_engine" ]]; then
        printf 'Install a Beamer-capable LaTeX engine or set TECTONIC_BIN.\n' >&2
        exit 1
    fi
    deck_flags=(--keep-logs --keep-intermediates --outdir build)
    if [[ "${LESTO_BEAMER_OFFLINE:-0}" == 1 ]]; then deck_flags+=(--only-cached); fi
    nice -n 10 "$deck_engine" "${deck_flags[@]}" main.tex > build/compile.stdout 2>&1
fi
cp build/main.pdf lesto-m10-review.pdf
if command -v pdfinfo >/dev/null 2>&1; then
    deck_pages="$(pdfinfo lesto-m10-review.pdf | awk '/^Pages:/{print $2}')"
    if [[ "$deck_pages" != 20 ]]; then
        printf 'Expected 20 physical pages, found %s.\n' "$deck_pages" >&2
        exit 1
    fi
fi
printf 'Built %s/lesto-m10-review.pdf (20 slides).\n' "$PWD"
