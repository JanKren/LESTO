#!/usr/bin/env bash
# Build the International Journal of Heat and Mass Transfer (Elsevier,
# elsarticle) version of the manuscript.
# Usage: bash manuscript/ijhmt/build.sh              # working draft: arxiv.sty, line numbers
#        bash manuscript/ijhmt/build.sh --elsarticle # Elsevier review copy (12 pt, 1.5 spacing)
#        bash manuscript/ijhmt/build.sh --print      # Elsevier two-column 5p, length estimate
set -euo pipefail
ijhmt_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$ijhmt_dir"
mode=draft
case "${1:-}" in
    '') ;;
    --elsarticle) mode=elsarticle ;;
    --print) mode=print ;;
    *) printf 'Usage: bash manuscript/ijhmt/build.sh [--elsarticle|--print]\n' >&2
       exit 2 ;;
esac
if [[ -n "${TECTONIC_BIN:-}" ]]; then
    compiler="$TECTONIC_BIN"
elif command -v tectonic >/dev/null 2>&1; then
    compiler="$(command -v tectonic)"
elif [[ -x /home/jan/.local/opt/tectonic-0.17.0/tectonic ]]; then
    compiler=/home/jan/.local/opt/tectonic-0.17.0/tectonic
else
    printf 'Install Tectonic or compile lesto-ijhmt.tex with latexmk -pdf.\n' >&2
    exit 1
fi
src=lesto-ijhmt.tex
[[ "$mode" == elsarticle ]] && src=lesto-ijhmt-elsarticle.tex
if [[ "$mode" == print ]]; then
    # Elsevier's final two-column layout, single spacing, no line numbers.
    src=lesto-ijhmt-print.tex
    sed -e 's/\[preprint,12pt\]{elsarticle}/[final,5p,times,twocolumn]{elsarticle}/' \
        -e '/^\\usepackage\[margin=1in\]{geometry}$/d' -e 's/^\\onehalfspacing$//' -e 's/^\\linenumbers$//' lesto-ijhmt-elsarticle.tex > "$src"
fi
mkdir -p build
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
    "$compiler" --keep-logs --keep-intermediates --outdir build "$src"
out="${src%.tex}.pdf"
cp "build/$out" "$out"
[[ "$mode" == print ]] && rm -f "$src"
log="build/${src%.tex}.log"
printf 'Built %s/%s (%s pages)\n' "$ijhmt_dir" "$out" \
    "$(pdfinfo "$out" | awk '/^Pages:/{print $2}')"
printf 'Overfull boxes: %s; undefined references/citations: %s\n' \
    "$(grep -c '^Overfull' "$log" || true)" \
    "$(grep -ciE 'undefined' "$log" || true)"
# Journal limits that the source can check: abstract <= 250 words (inline
# math counts as one word, % comments excluded) and 3-5 highlights of at
# most 85 characters including spaces.
python3 - <<'EOF'
import re
t = open('../sections/abstract-journal.tex').read()
t = '\n'.join(l for l in t.splitlines() if not l.lstrip().startswith('%'))
t = re.sub(r'\$[^$]*\$', 'M', t)
words = len(t.split())
h = [l.rstrip('\n') for l in open('highlights.txt') if l.strip()]
long_h = [len(l) for l in h if len(l) > 85]
print(f'Abstract: {words} words (limit 250); highlights: {len(h)} '
      f'(3-5), longest {max(map(len, h))} characters (limit 85)')
if words > 250 or not 3 <= len(h) <= 5 or long_h:
    raise SystemExit('IJHMT limit exceeded')
EOF
# Supplementary material: the shared ../supplementary.tex, built from the
# manuscript root so that its generated/ and references.bib paths resolve.
if [[ "$mode" == draft ]]; then
    (cd .. && XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/lesto-beamer-cache}" \
        "$compiler" --keep-logs --keep-intermediates --outdir ijhmt/build supplementary.tex)
    cp build/supplementary.pdf lesto-ijhmt-supplementary.pdf
    slog=build/supplementary.log
    printf 'Built %s/lesto-ijhmt-supplementary.pdf\n' "$ijhmt_dir"
    printf 'Supplement overfull boxes: %s; undefined references/citations: %s\n' \
        "$(grep -c '^Overfull' "$slog" || true)" \
        "$(grep -ciE 'undefined' "$slog" || true)"
fi
