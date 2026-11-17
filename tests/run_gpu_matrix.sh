#!/bin/bash
# Manual pre-release GPU check - not run in CI (GitHub-hosted runners have no GPU).
#
# Run this by hand on a workstation with an NVIDIA GPU before a release, after
# building with -DENABLE_CUDA_PLUGIN=ON. It runs:
#   1. The CPU regression backend matrix (ctest --label-regex backend), same
#      as CI, against the CUDA-enabled binary.
#   2. A CPU-vs-GPU consistency pass over every recur='chebyshev' regression
#      case (tests/regression/run_gpu_matrix.py), since gpu_plugin bypasses
#      cheb_backend and dispatches straight to CUDA.
#
# Usage:
#   tests/run_gpu_matrix.sh [build_dir]
#
# build_dir defaults to "build" and must already be configured+built with
# -DENABLE_CUDA_PLUGIN=ON (see docs/source/getting_started.rst).

set -euo pipefail
cd "$(dirname "$0")/.."

BUILD_DIR="${1:-build}"
BINARY="${BUILD_DIR}/bin/rslmto.x"

if [ ! -x "$BINARY" ]; then
    echo "ERROR: $BINARY not found or not executable."
    echo "Configure and build with -DENABLE_CUDA_PLUGIN=ON -DRUN_REG_TESTS=ON first, e.g.:"
    echo "  cmake -S . -B $BUILD_DIR -DENABLE_CUDA_PLUGIN=ON -DRUN_REG_TESTS=ON"
    echo "  cmake --build $BUILD_DIR --parallel"
    exit 1
fi

echo "=== 1/2: CPU regression backend matrix (${BINARY}) ==="
ctest --test-dir "$BUILD_DIR" --label-regex backend --output-on-failure --parallel 1
backend_status=$?

echo
echo "=== 2/2: CPU-vs-GPU consistency (chebyshev regression cases) ==="
python3 tests/regression/run_gpu_matrix.py --binary "$BINARY"
gpu_status=$?

echo
if [ "$backend_status" -eq 0 ] && [ "$gpu_status" -eq 0 ]; then
    echo "RESULT: all GPU matrix checks passed."
    exit 0
else
    echo "RESULT: GPU matrix FAILED (backend_status=$backend_status, gpu_status=$gpu_status)."
    exit 1
fi
