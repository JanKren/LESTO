#!/bin/bash
#------------------------------------------------------------------------------
#
# tests/regress.sh
#
# PURPOSE
#
# Compares a rhoFixedFlowFoam log with a reference, line by line, after
# removing what legitimately differs between two runs of the same solver.
#
# Usage
#
#   regress.sh <reference log> <log>
#       Complete logs, normally of two builds on the same OpenFOAM.  Every
#       filtered line must be identical, including the residual history and
#       the final "End".
#
#   regress.sh -excerpt <out_excerpt> <log>
#       Reference is the out_excerpt of a run/ case, written by another
#       OpenFOAM version on another host.  Its first block, from "Create
#       time" up to its first "..." line, is compared with the same number of
#       filtered lines of the log, from "Create time" on.
#
#   regress.sh [-excerpt] -noMockReport <reference> <log>
#       As above, after removing from both files the report that the
#       coupled mock of tag 0.1 prints in every step: the lines
#       "Integrated source <species> = ... kg/s" and "Net thermochemistry
#       source = ... kg/s", and the balance blocks from the blank line
#       before the first "Species balance <species>:" to the four lines of
#       "Tracked species total:".  For the out_excerpts of 008-010, which
#       were recorded before tag 0.1 added the report.
#
#   regress.sh [-excerpt] [-noMockReport] -ignore <regex> <reference> <log>
#       As above, after removing from both files every line that matches the
#       extended regular expression <regex> (e.g. the lines of the species
#       whose results tag 0.1 changed in the out_excerpt of 010).  The
#       number of removed reference lines is printed with the result.
#
# Filter (applied to both files)
#
#   - removed: the header lines that identify the run (Exec, Date, Time,
#     Host, PID, Case) and every ExecutionTime line;
#   - removed: the compilation chatter of dynamicCode, which appears only
#     when a case's coded boundary-condition library does not exist yet
#     (Could not load ..., cannot open shared object ..., Creating new
#     library ..., Invoking wmake libso ..., wmake libso ..., and the
#     indented ln:/dep:/Ctoo:/link: lines of wmake);
#   - replaced: the case directory, taken from the log's own "Case :" line,
#     by <case>.  The only remaining lines that contain it are the "Using
#     dynamicCode for patch ..." lines, which are thus compared too.
#
# Prints "IDENTICAL: <n> lines (<m> time steps)", or "DIFFERENT" followed by
# the diff (reference <, log >).  Exit status 0 when identical, 1 when
# different, 2 on a usage error.
#
#------------------------------------------------------------------------------

# the per-step report of the coupled mock of tag 0.1 (-noMockReport)
dropMockReport() {
  awk '
    /^Integrated source .* kg\/s$/ { next }
    /^Net thermochemistry source = .* kg\/s$/ { next }
    skip == 1 {
      if ($0 == "Tracked species total:") { skip = 2; n = 0 }
      next
    }
    skip == 2 { if (++n >= 4) skip = 0; next }
    blank {
      blank = 0
      if ($0 ~ /^Species balance .*:$/) { skip = 1; next }
      print ""
    }
    $0 == "" { blank = 1; next }
    { print }
    END { if (blank) print "" }
  '
}

filterLog() {
  local caseDir
  caseDir=$(sed -n 's/^Case   *: *//p' "$1" | head -1)

  sed -e '/^\(Exec\|Date\|Time\|Host\|PID\|Case\) *: /d' \
      -e '/^ExecutionTime = /d' \
      -e '/^Could not load "/d' \
      -e '/: cannot open shared object file: /d' \
      -e '/^Creating new library in "/d' \
      -e '/^Invoking wmake libso /d' \
      -e '/^wmake libso /d' \
      -e '/^    \(ln\|dep\|Ctoo\|link\): /d' \
      "$1" \
  | if [ -n "$caseDir" ]
    then
      sed -e "s|$caseDir|<case>|g"
    else
      cat
    fi \
  | if $noMockReport
    then
      dropMockReport
    else
      cat
    fi
}

# the lines that do not match the -ignore expression (all without one)
dropIgnored() {
  if [ -n "$ignore" ]
  then
    grep -v -E -- "$ignore"
  else
    cat
  fi
}

excerpt=false
noMockReport=false
ignore=
while [ $# -gt 0 ]
do
  case "$1" in
    -excerpt) excerpt=true ;;
    -noMockReport) noMockReport=true ;;
    -ignore) ignore=${2:-}; shift ;;
    *) break ;;
  esac
  shift
done

if [ $# -ne 2 ] || [ ! -f "$1" ] || [ ! -f "$2" ]
then
  echo "Usage: regress.sh [-excerpt] [-noMockReport] [-ignore <regex>]" \
       "<reference> <log>" >&2
  exit 2
fi

ref=$(mktemp)
new=$(mktemp)
all=$(mktemp)
trap 'rm -f "$ref" "$new" "$all"' EXIT

if $excerpt
then
  filterLog "$1" | sed -n '/^Create time/,$p' | sed '/^\.\.\.$/,$d' > "$all"
  dropIgnored < "$all" > "$ref"
  filterLog "$2" | sed -n '/^Create time/,$p' | dropIgnored \
    | head -n "$(wc -l < "$ref")" > "$new"
else
  filterLog "$1" > "$all"
  dropIgnored < "$all" > "$ref"
  filterLog "$2" | dropIgnored > "$new"
fi

ignored=
if [ -n "$ignore" ]
then
  ignored=", $(( $(wc -l < "$all") - $(wc -l < "$ref") )) reference lines\
 matching '$ignore' left out"
fi

if [ ! -s "$ref" ]
then
  echo "Empty reference after filtering: $1" >&2
  exit 2
fi

if cmp -s "$ref" "$new"
then
  echo "IDENTICAL: $(wc -l < "$ref") lines" \
       "($(grep -c '^Time = ' "$ref") time steps$ignored)"
  exit 0
fi

echo "DIFFERENT"
diff "$ref" "$new"
exit 1

#------------------------------------------------------------------------------
