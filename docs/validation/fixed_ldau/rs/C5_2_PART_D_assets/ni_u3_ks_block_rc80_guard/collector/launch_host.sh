#!/bin/bash
set -eu
ROOT=/private/tmp/c52-rd-ni-followups/ni_u3_ks_block_rc80
RUN="$ROOT/run"
BIN=/private/tmp/c5-b1-frozen/a2112df3224a376c280ec5398eebeb68b6984f68dfef1138dc9b5d40432298d3/rslmto.x
EXPECTED_BIN=a2112df3224a376c280ec5398eebeb68b6984f68dfef1138dc9b5d40432298d3
EXPECTED_INPUT=4a9f5009f2d25becb20af6ae2fb7482c5372af96b282734fc06de2728eacc16d
EXPECTED_SEED=bc3340a8ad20d2467eaeccd43929a659d4cfff2a9f28b748fb6554edc4ee3ac9
EXPECTED_ELEMENT=4824dcfb4e01db54fa0ca1e4fe34b31cc398e05382463807504c21b3c5c18d57
GUARD="$ROOT/collector/disk_guard.sh"
EXPECTED_GUARD=3049362e51a3165866c73bb997667a0c74662fc673328085b2636c25dae74cdd
FLOOR=1073741824
MIN_PREEXEC=2684354560
sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
available_bytes() { df -Pk "$RUN" | awk 'NR==2 {print $4 * 1024}'; }
if [ "${C52_RD_NI_U3_KS_AUTHORIZED:-0}" != 1 ]; then echo 'refusing launch: explicit authorization env is unset' >&2; exit 77; fi
if [ -e "$RUN/solver.pid" ] || [ -e "$RUN/exit_code.txt" ] || [ -e "$RUN/lifecycle.txt" ]; then echo 'refusing launch: run directory has prior lifecycle metadata' >&2; exit 78; fi
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ] || [ "$(sha256 "$GUARD")" != "$EXPECTED_GUARD" ]; then echo 'frozen executable or guard hash mismatch' >&2; exit 79; fi
if [ "$(sha256 "$ROOT/immutable/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$ROOT/immutable/Ni-imp_out.nml")" != "$EXPECTED_SEED" ] || [ "$(sha256 "$ROOT/immutable/Ni-imp.nml")" != "$EXPECTED_ELEMENT" ]; then echo 'immutable input or seed hash mismatch' >&2; exit 80; fi
if [ "$(sha256 "$RUN/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$RUN/Ni-imp_out.nml")" != "$EXPECTED_SEED" ] || [ "$(sha256 "$RUN/Ni-imp.nml")" != "$EXPECTED_ELEMENT" ]; then echo 'working copy hash mismatch' >&2; exit 81; fi
for f in "$RUN/input.nml" "$RUN/Ni-imp_out.nml" "$RUN/Ni-imp.nml"; do
  if [ ! -w "$f" ]; then echo "working file is not writable: $f" >&2; exit 82; fi
  perm=$(stat -f '%Lp' "$f"); owner_digit=$(printf '%s' "$perm" | cut -c1)
  case "$owner_digit" in 2|3|6|7) ;; *) echo "working file lacks owner-write bit: $f mode=$perm" >&2; exit 83 ;; esac
done
available=$(available_bytes)
if [ "$available" -lt "$MIN_PREEXEC" ]; then echo "prelaunch disk check: ${available} bytes free, below $MIN_PREEXEC" >&2; exit 84; fi
{
  printf '%s\n' "run_id=rd_ni_u3_ks_bulk_rc80" "wrapper_pid=$$" "preexec_free_bytes=$available" "binary_sha256=$EXPECTED_BIN" "input_sha256=$EXPECTED_INPUT" "seed_sha256=$EXPECTED_SEED" "element_sha256=$EXPECTED_ELEMENT" "omp_num_threads=1" "omp_dynamic=FALSE" "rslmto_omp_threads_serial=1" "openblas_num_threads=1" "mkl_num_threads=1" "veclib_maximum_threads=1" "blis_num_threads=1" "numexpr_num_threads=1" "fixed_u_trace=$RUN/fixed_u_trace.csv" "started_local=$(date '+%Y-%m-%dT%H:%M:%S%z')"
} > "$RUN/lifecycle.txt"
cd "$RUN"
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=1 \
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  RSLMTO_FIXED_U_TRACE="$RUN/fixed_u_trace.csv" \
  RSLMTO_FIXED_U_RUN_ID=rd_ni_u3_ks_bulk_rc80 \
  RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=1 RSLMTO_FIXED_U_TRACE_QL_ONLY=0 \
  RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0 \
  RSLMTO_FIXED_U_A2_SCAN=0 "$BIN" > testrun.log 2>&1 &
solver_pid=$!
printf '%s\n' "$solver_pid" > "$RUN/solver.pid"
printf '%s\n' "$solver_pid" > "$RUN/solver_pid.actual"
"$GUARD" "$RUN" "$solver_pid" > "$RUN/disk_guard.log" 2>&1 &
guard_pid=$!
printf '%s\n' "$guard_pid" > "$RUN/disk_guard.pid"
rc=0
wait "$solver_pid" || rc=$?
printf '%s\n' "$rc" > "$RUN/exit_code.txt"
printf '%s\n' "solver_pid=$solver_pid" "disk_guard_pid=$guard_pid" "exit_code=$rc" "ended_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "$RUN/lifecycle.txt"
kill "$guard_pid" 2>/dev/null || true
wait "$guard_pid" 2>/dev/null || true
exit "$rc"
