#!/bin/bash
set -eu
ROOT=/private/tmp/c52-rd-ni-u0-imp-control
RUN="$ROOT/run"
WATCHER="$ROOT/collector/watch_trace_and_state.py"
GUARD="$ROOT/collector/disk_guard.sh"
BIN=/private/tmp/c5-rd1-frozen/4bb76cc2bbd251d95cdc434eec0348b9082069af2deeb4f4859f6d8056779a4f/rslmto.x
EXPECTED_BIN=4bb76cc2bbd251d95cdc434eec0348b9082069af2deeb4f4859f6d8056779a4f
EXPECTED_INPUT=008bce3b5160ef76b40201b2d03797db346da9eef9a07046b5150907b90538b7
EXPECTED_HOST=736a3c4e8fd2e1cdc7550502c4788e1b9f0a031c06a44be2ee77a8fd5b1ac1c4
EXPECTED_IMP=0213bf83b07bd82446f293bf19be8ba37b0a2bc816011d29f2c7e8d5c483738a
EXPECTED_BULK_ATOM=7385d1f467c336795f71bfa1a8dff292cd3ec500bc00a9182324c060ab2563ac
EXPECTED_IMP_ATOM=4824dcfb4e01db54fa0ca1e4fe34b31cc398e05382463807504c21b3c5c18d57
RUN_ID=rd_ni_i_u0_local_from_rs_bulk_u0_rc80
FLOOR=1073741824
PYTHON=/opt/homebrew/opt/python@3.11/bin/python3.11
sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
available_bytes() { df -Pk "$RUN" | awk 'NR==2 {print $4 * 1024}'; }
if [ "${C52_RD_NI_I_U0_CONTROL_AUTHORIZED:-0}" != 1 ]; then echo 'refusing launch: root authorization env unset' >&2; exit 77; fi
if [ -e "$RUN/solver.pid" ] || [ -e "$RUN/exit_code.txt" ] || [ -e "$RUN/lifecycle.txt" ]; then echo 'refusing launch: nonempty prior lifecycle' >&2; exit 78; fi
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ]; then echo 'frozen executable hash mismatch' >&2; exit 79; fi
if [ "$(sha256 "$ROOT/immutable/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$ROOT/immutable/Ni-bulk_out.nml")" != "$EXPECTED_HOST" ] || [ "$(sha256 "$ROOT/immutable/Ni-imp_out.nml")" != "$EXPECTED_IMP" ]; then echo 'immutable inputs/restart mismatch' >&2; exit 80; fi
for spec in "input.nml:$EXPECTED_INPUT" "Ni-bulk_out.nml:$EXPECTED_HOST" "Ni-imp_out.nml:$EXPECTED_IMP" "Ni-bulk.nml:$EXPECTED_BULK_ATOM" "Ni-imp.nml:$EXPECTED_IMP_ATOM"; do
  f=${spec%%:*}; want=${spec#*:}
  if [ "$(sha256 "$RUN/$f")" != "$want" ]; then echo "working copy hash mismatch: $f" >&2; exit 81; fi
  perm=$(stat -f '%Lp' "$RUN/$f"); owner_digit=$(printf '%s' "$perm" | cut -c1)
  case "$owner_digit" in 2|3|6|7) ;; *) echo "working input lacks owner-write bit: $f mode=$perm" >&2; exit 88 ;; esac
done
available=$(available_bytes); if [ "$available" -lt "$FLOOR" ]; then echo "disk below floor: $available" >&2; exit 85; fi
printf '%s\n' "launcher_pid=$$" "preexec_free_bytes=$available" "started_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$RUN/lifecycle.txt"
"$PYTHON" "$WATCHER" --run-id "$RUN_ID" --solver-pid-file "$RUN/solver.pid" --run-dir "$RUN" --state-file Ni-imp_scf.nml --output-dir "$ROOT/captures" --ready-file "$RUN/watcher.ready" > "$RUN/watcher.log" 2>&1 &
watcher_pid=$!; printf '%s\n' "$watcher_pid" > "$RUN/watcher.pid"
tries=0
while [ ! -s "$RUN/watcher.ready" ]; do
  if ! kill -0 "$watcher_pid" 2>/dev/null; then wait "$watcher_pid" 2>/dev/null || true; echo 'watcher exited before readiness' >&2; exit 82; fi
  tries=$((tries+1)); if [ "$tries" -gt 1000 ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo 'watcher readiness timeout' >&2; exit 83; fi
  sleep .1
done
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ] || [ "$(sha256 "$RUN/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$RUN/Ni-bulk_out.nml")" != "$EXPECTED_HOST" ] || [ "$(sha256 "$RUN/Ni-imp_out.nml")" != "$EXPECTED_IMP" ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo 'pre-exec hash mismatch' >&2; exit 84; fi
available=$(available_bytes); if [ "$available" -lt "$FLOOR" ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo "disk below floor before exec: $available" >&2; exit 85; fi
cd "$RUN"
env OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=4 \
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  RSLMTO_FIXED_U_TRACE="$RUN/fixed_u_trace.csv" \
  RSLMTO_FIXED_U_RUN_ID="$RUN_ID" \
  RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=1 RSLMTO_FIXED_U_TRACE_QL_ONLY=0 \
  RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0 \
  RSLMTO_FIXED_U_A2_SCAN=0 "$BIN" > testrun.log 2>&1 &
solver_pid=$!; printf '%s\n' "$solver_pid" > "$RUN/solver.pid"; printf '%s\n' "$solver_pid" > "$RUN/solver_pid.actual"
printf '%s\n' "solver_pid=$solver_pid" "watcher_pid=$watcher_pid" "env_OMP_NUM_THREADS=4" "env_RSLMTO_OMP_THREADS_SERIAL=4" "env_library_thread_caps=1" "exec_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "$RUN/lifecycle.txt"
"$GUARD" "$RUN" "$solver_pid" > "$RUN/disk_guard.log" 2>&1 &
guard_pid=$!; printf '%s\n' "$guard_pid" > "$RUN/disk_guard.pid"
rc=0; wait "$solver_pid" || rc=$?
printf '%s\n' "$rc" > "$RUN/exit_code.txt"
printf '%s\n' "exit_code=$rc" "ended_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "$RUN/lifecycle.txt"
kill "$guard_pid" 2>/dev/null || true; wait "$guard_pid" 2>/dev/null || true
wait "$watcher_pid" || true
exit "$rc"
