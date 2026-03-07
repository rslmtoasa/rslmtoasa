#!/bin/bash

set -euo pipefail

if [ "$#" -ne 16 ]; then
    echo "Usage: $0 <python_exec> <runner> <binary> <case_dir> <scratch_root> <test_name> <nsp> <recur> <hoh> <lld> <nstep> <manifest> <references_dir> <keys> <abs_tol> <rel_tol>"
    exit 2
fi

python_exec="$1"
runner="$2"
binary="$3"
case_dir="$4"
scratch_root="$5"
test_name="$6"
nsp="$7"
recur="$8"
hoh="$9"
shift 9
lld="$1"
nstep="$2"
manifest="$3"
references_dir="$4"
keys="$5"
abs_tol="$6"
rel_tol="$7"

script_dir="$(cd "$(dirname "$0")" && pwd)"
compare_script="${script_dir}/compare_reference_data.py"

if [ ! -x "$runner" ]; then
    echo "ERROR: runner not executable: $runner"
    exit 1
fi

if [ ! -x "$python_exec" ]; then
    echo "ERROR: python executable not found: $python_exec"
    exit 1
fi

if [ ! -f "$compare_script" ]; then
    echo "ERROR: compare script not found: $compare_script"
    exit 1
fi

/bin/bash "$runner" "$binary" "$case_dir" "$scratch_root" "$test_name" "$nsp" "$recur" "$hoh" "$lld" "$nstep"

"$python_exec" "$compare_script" \
    --manifest "$manifest" \
    --references-dir "$references_dir" \
    --scratch-root "$scratch_root" \
    --pattern "$test_name" \
    --keys "$keys" \
    --abs-tol "$abs_tol" \
    --rel-tol "$rel_tol" \
    --include-cheb
