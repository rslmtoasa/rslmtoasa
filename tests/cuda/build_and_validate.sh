#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cuda_root="${CUDA_HOME:-/usr/local/cuda}"
out_dir="${1:-/tmp/rslmto-acc01-lowlevel}"

mkdir -p "$out_dir"
export LD_LIBRARY_PATH="$cuda_root/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

"$cuda_root/bin/nvcc" -O3 -arch=native -Xcompiler -fPIC -shared \
  "$repo_dir/source/cuda/rsrec_gpu.cu" \
  -o "$out_dir/librsrec_gpu.so" -lcublas -lcufft -llapack

cd "$repo_dir"
python3 tests/cuda/rsrec_validate.py --library "$out_dir/librsrec_gpu.so"
