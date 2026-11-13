#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# f90nml is not an Ubuntu 24.04 apt package. Keep pip out of system Python.
python3 -m venv --system-site-packages .venv-tutorial
.venv-tutorial/bin/python -m pip install 'f90nml==1.5.0'
export PATH="$repo_root/.venv-tutorial/bin:$PATH"

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
if [ ! -e runs/bccfe-scf ]; then
    cp -r tutorials/bccFe runs/bccfe-scf
fi

printf '\nReady. Follow docs/rs-tutorial.md. Your Fe working deck is runs/bccfe-scf/input.nml.\n'

# Install the release-prep wheel bundled with this tutorial branch.
bash .devcontainer/install-rsutils.sh
