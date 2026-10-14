#!/usr/bin/env bash
set -euo pipefail
package_dir=$(cd -- "$(dirname -- "$0")" && pwd)
if [[ "${WM_PROJECT_VERSION:-}" != "v2412" ]]; then
    foam_setup="${LESTO_FOAM_BASHRC:-$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc}"
    if [[ ! -f "$foam_setup" ]]; then
        echo 'Load OpenFOAM v2412 first, or set LESTO_FOAM_BASHRC to its etc/bashrc.' >&2
        exit 1
    fi
    set +eu
    source "$foam_setup"
    set -eu
fi
export MPLBACKEND=Agg
python3 "$package_dir/run_study.py" check "$@"
echo 'Checks passed. Starting the study in this terminal; use the nohup command in README to detach.'
exec python3 "$package_dir/run_study.py" run "$@"
