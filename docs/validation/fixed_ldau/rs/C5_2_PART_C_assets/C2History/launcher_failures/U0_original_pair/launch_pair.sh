#!/bin/bash
set -u
BIN=/private/tmp/c5-b1-frozen/a2112df3224a376c280ec5398eebeb68b6984f68dfef1138dc9b5d40432298d3/rslmto.x
ROOT=/private/tmp/c52-c2-stage/launch_u0
start_run() {
  run_id="$1"
  run_dir="$ROOT/$run_id"
  (
    cd "$run_dir" || exit 90
    echo "$BASHPID" > wrapper.pid
    env OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=4 \
      OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
      BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      RSLMTO_FIXED_U_TRACE="$run_dir/fixed_u_trace.csv" \
      RSLMTO_FIXED_U_RUN_ID="$run_id" \
      RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=1 RSLMTO_FIXED_U_TRACE_QL_ONLY=0 \
      RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0 \
      RSLMTO_FIXED_U_A2_SCAN=0 "$BIN" > testrun.log 2>&1 &
    solver_pid=$!
    echo "$solver_pid" > solver.pid
    echo "solver_pid=$solver_pid wrapper_pid=$BASHPID"
    wait "$solver_pid"
    rc=$?
    echo "$rc" > exit_code.txt
    echo "exit_code=$rc"
    exit 0
  ) &
  echo "$!" > "$run_dir/launcher_subshell.pid"
}
start_run c2_u0_block
start_run c2_u0_chebyshev
wait
