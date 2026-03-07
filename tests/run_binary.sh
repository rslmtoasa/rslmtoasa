#!/bin/bash
# Thin wrapper for binary invocation.
# Usage: run_binary.sh <binary> [mpi_procs]
# Must be run from the test scratch directory (cwd). Output is written to testrun.log.

set -euo pipefail

binary="$1"
mpi_procs="${2:-1}"

if [ ! -x "$binary" ]; then
    echo "ERROR: binary not executable: $binary" >&2
    exit 1
fi

if [ "$mpi_procs" -gt 1 ]; then
    OMP_NUM_THREADS=1 mpirun -n "$mpi_procs" "$binary" > testrun.log 2>&1
else
    OMP_NUM_THREADS=1 "$binary" > testrun.log 2>&1
fi
