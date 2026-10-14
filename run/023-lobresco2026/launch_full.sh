#!/bin/bash
# Full runs of the five Lobresco cases.  Usage:
#   bash launch_full.sh <level: prod|coarse> <first core> [cases...]
# Each case is one serial process pinned to its own core (core, core+1, ...).
# Cases start from the prepared 0/ (startFrom latestTime: a stopped run resumes).
set -euo pipefail
G=$(cd "$(dirname "$0")" && pwd)
level=${1:?level}; core=${2:?first core}; shift 2
cases=${*:-T106_6 T106_60 T106_300 T5_60 T540_60}
set +u; source ~/OpenFOAM/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1; set -u
echo "solver sha256: $(sha256sum "$G/bin/rhoFixedFlowFoam" | cut -d' ' -f1)"
for c in $cases; do
  d="$G/cases/$level/$c"
  nohup taskset -c "$core" "$G/bin/rhoFixedFlowFoam" -case "$d" >> "$d/log.rhoFixedFlowFoam" 2>&1 < /dev/null &
  echo "$c on core $core pid $!"
  core=$((core + 1))
done
# After completion:  python3 "$G/analyse.py" "$G"/cases/$level/T*
