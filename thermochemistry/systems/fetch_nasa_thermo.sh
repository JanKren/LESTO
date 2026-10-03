#!/bin/bash
# NASA Glenn thermo.inp (Gurvich 1991 Pb-I data), Apache-2.0, pinned commit.
set -euo pipefail
curl -sSL -o nasa_thermo.inp \
  https://raw.githubusercontent.com/nasa/cea/3f4441d2/data/thermo.inp
ls -l nasa_thermo.inp
