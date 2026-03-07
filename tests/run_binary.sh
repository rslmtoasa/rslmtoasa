#!/bin/bash
# Thin wrapper for binary invocation.
# Usage: run_binary.sh <binary>
# Must be run from the test scratch directory (cwd). Output is written to testrun.log.

set -euo pipefail

binary="$1"

if [ ! -x "$binary" ]; then
    echo "ERROR: binary not executable: $binary" >&2
    exit 1
fi

OMP_NUM_THREADS=1 "$binary" > testrun.log 2>&1
