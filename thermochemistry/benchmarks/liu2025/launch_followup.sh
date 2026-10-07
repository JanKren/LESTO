#!/bin/bash
# Durable pipeline: audit, native smoke, full cases, then automatic reporting.
set -euo pipefail
study_root=${1:?Study directory required}
study_jobs=${2:-16}
report_out=${3:?Report directory required}
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
if test -z "${WM_PROJECT_DIR:-}"; then
    set +eu
    source /home/jan/OpenFOAM/OpenFOAM-v2412/etc/bashrc
    set -eu
fi
export MPLCONFIGDIR=/tmp/lesto-matplotlib
while ! test -f "$study_root/study.json"; do
    if test -f /tmp/lesto-m10e-followup-prepare.log && rg -q '^Traceback' /tmp/lesto-m10e-followup-prepare.log; then
        echo 'Preparation failed; inspect /tmp/lesto-m10e-followup-prepare.log' >&2
        exit 1
    fi
    sleep 10
done
if ! test -f "$study_root/smoke.json"; then
    python3 "$script_dir/followup_study.py" finalise --out "$study_root"
    python3 "$script_dir/followup_study.py" smoke --out "$study_root" --jobs 8
fi
python3 "$script_dir/followup_study.py" run --out "$study_root" --jobs "$study_jobs" --report-out "$report_out"
