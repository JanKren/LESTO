#!/bin/bash
#------------------------------------------------------------------------------
#
# tests/cmpFields.sh
#
# PURPOSE
#
# Byte-for-byte comparison of the files written by two runs.
#
# Usage
#
#   cmpFields.sh [-x <path>]... <case A> <case B> [<time> ...]
#
# For every listed time (default: every time directory of case A except 0)
# the script compares the files below <time>/ and below each
# processor*/<time>/ of case A with the same files of case B, using cmp.
# The two cases must hold exactly the same set of files there.
#
#   -x <path>   leave out files whose path relative to the time directory
#               is <path>, e.g. -x uniform/time (may be repeated)
#
# Prints "IDENTICAL: <n> files in <m> time directories", or the files that
# differ or exist in only one case.  Exit status 0 when identical, 1 when
# not, 2 on a usage error.
#
#------------------------------------------------------------------------------

export LC_ALL=C

excludes=()
while [ "${1:-}" = "-x" ]
do
  excludes+=("$2")
  shift 2
done

if [ $# -lt 2 ] || [ ! -d "$1" ] || [ ! -d "$2" ]
then
  echo "Usage: cmpFields.sh [-x <path>]... <case A> <case B> [<time> ...]" >&2
  exit 2
fi

caseA=$1
caseB=$2
shift 2

if [ $# -gt 0 ]
then
  times=("$@")
else
  mapfile -t times < <(cd "$caseA" && ls -d [0-9]* 2>/dev/null \
    | grep -v -x '0' | sort -g)
fi

if [ ${#times[@]} -eq 0 ]
then
  echo "No time directories to compare in $caseA" >&2
  exit 2
fi

# files below the listed times of one case, relative to the case, sorted
listFiles() {
  local t d
  for t in "${times[@]}"
  do
    for d in "$1/$t" "$1"/processor*/"$t"
    do
      [ -d "$d" ] && (cd "$1" && find "${d#$1/}" -type f)
    done
  done | sort | while read -r f
  do
    local keep=true x
    for x in "${excludes[@]}"
    do
      case "$f" in
        */"$x") keep=false ;;
      esac
    done
    $keep && echo "$f"
  done
}

filesA=$(mktemp)
filesB=$(mktemp)
trap 'rm -f "$filesA" "$filesB"' EXIT

listFiles "$caseA" > "$filesA"
listFiles "$caseB" > "$filesB"

status=0
nDirs=0
for t in "${times[@]}"
do
  for d in "$caseA/$t" "$caseA"/processor*/"$t"
  do
    [ -d "$d" ] && nDirs=$((nDirs + 1))
  done
done

if [ ! -s "$filesA" ]
then
  echo "No files below ${times[*]} in $caseA" >&2
  exit 2
fi

onlyA=$(comm -23 "$filesA" "$filesB")
onlyB=$(comm -13 "$filesA" "$filesB")
if [ -n "$onlyA" ]
then
  echo "Only in $caseA:"
  echo "$onlyA" | sed 's/^/  /'
  status=1
fi
if [ -n "$onlyB" ]
then
  echo "Only in $caseB:"
  echo "$onlyB" | sed 's/^/  /'
  status=1
fi

nSame=0
while read -r f
do
  if cmp -s "$caseA/$f" "$caseB/$f"
  then
    nSame=$((nSame + 1))
  else
    echo "Differs: $f"
    status=1
  fi
done < <(comm -12 "$filesA" "$filesB")

if [ $status -eq 0 ]
then
  echo "IDENTICAL: $nSame files in $nDirs time directories"
fi
exit $status

#------------------------------------------------------------------------------
