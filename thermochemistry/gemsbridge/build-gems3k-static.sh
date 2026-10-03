#!/bin/bash
# Static GEMS3K for the C bridge: no ThermoFun, header-only spdlog/fmt, hidden
# symbols.  Needs the conda-forge toolchain (GEMS3K requires C++20).
set -euo pipefail
GEMS3K_SRC=${GEMS3K_SRC:-$HOME/src/GEMS3K}
GEMS_ENV=${GEMS_ENV:-$HOME/opt/gems-env}
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate "$GEMS_ENV"
mkdir -p "$GEMS3K_SRC/build-static-bridge"
cd "$GEMS3K_SRC/build-static-bridge"
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DUSE_THERMOFUN=OFF -DUSE_SPDLOG_PRECOMPILED=OFF \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON \
  -DBUILD_SOLMOD=OFF -DBUILD_SOLMOD_PYTHON=OFF -DBUILD_TOOLS=OFF -DBUILD_TESTING=OFF \
  -DCMAKE_CXX_FLAGS="-fvisibility=hidden -fvisibility-inlines-hidden -DFMT_HEADER_ONLY"
ninja
ls -l lib/libGEMS3K-static.a
