#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

cmake -S . -B build-codespaces \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_OPENMP=ON \
    -DENABLE_MPI=OFF \
    -DENABLE_SPGLIB=ON \
    -DENABLE_MARCH_NATIVE=OFF \
    -DBUILD_TESTING=OFF
cmake --build build-codespaces --parallel 2

# Re-running setup must preserve students' inputs and results.
mkdir -p runs
if [ ! -e runs/bccFe ]; then
    cp -r example/bulk/bccFe runs/bccFe
fi

printf '\nReady. Open runs/bccFe/input.nml, then run:\n'
printf '  cd runs/bccFe\n  rslmto.x input.nml > run.log 2>&1\n'
