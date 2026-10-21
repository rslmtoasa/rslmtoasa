#!/bin/bash
set -eu
ROOT=/private/tmp/c52-rd-ni-followups/ni_u3_rs_block_rc80
RUN="$ROOT/run"
WATCHER="$ROOT/collector/watch_trace_and_state.py"
GUARD="$ROOT/collector/disk_guard.sh"
BIN=/private/tmp/c5-b1-frozen/a2112df3224a376c280ec5398eebeb68b6984f68dfef1138dc9b5d40432298d3/rslmto.x
EXPECTED_BIN=a2112df3224a376c280ec5398eebeb68b6984f68dfef1138dc9b5d40432298d3
EXPECTED_INPUT=$(shasum -a 256 "$ROOT/immutable/input.nml" | awk '{print $1}')
EXPECTED_SEED=0213bf83b07bd82446f293bf19be8ba37b0a2bc816011d29f2c7e8d5c483738a
EXPECTED_ELEMENT=4824dcfb4e01db54fa0ca1e4fe34b31cc398e05382463807504c21b3c5c18d57
FLOOR=1073741824
PYTHON=/opt/homebrew/opt/python@3.11/bin/python3.11
sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
available_bytes() { df -Pk "$RUN" | awk 'NR==2 {print $4 * 1024}'; }
if [ "${C52_RD_NI_U3_AUTHORIZED:-0}" != 1 ]; then echo 'refusing launch: explicit root authorization env is unset' >&2; exit 77; fi
if [ -e "$RUN/solver.pid" ] || [ -e "$RUN/exit_code.txt" ] || [ -e "$RUN/lifecycle.txt" ]; then echo 'refusing launch: run directory has prior lifecycle metadata' >&2; exit 78; fi
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ]; then echo 'frozen binary hash mismatch' >&2; exit 79; fi
if [ "$(sha256 "$ROOT/immutable/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$ROOT/immutable/Ni-imp_out.nml")" != "$EXPECTED_SEED" ]; then echo 'immutable input or seed hash mismatch' >&2; exit 80; fi
if [ "$(sha256 "$RUN/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$RUN/Ni-imp_out.nml")" != "$EXPECTED_SEED" ] || [ "$(sha256 "$RUN/Ni-imp.nml")" != "$EXPECTED_ELEMENT" ]; then echo 'working copy hash mismatch' >&2; exit 81; fi
if [ ! -w "$RUN" ]; then echo 'run directory is not writable' >&2; exit 86; fi
for f in "$RUN/input.nml" "$RUN/Ni-imp_out.nml" "$RUN/Ni-imp.nml"; do
  if [ ! -w "$f" ]; then echo "working input is not writable: $f" >&2; exit 87; fi
  perm=$(stat -f '%Lp' "$f"); owner_digit=$(printf '%s' "$perm" | cut -c1)
  case "$owner_digit" in 2|3|6|7) ;; *) echo "working file lacks owner-write bit: $f mode=$perm" >&2; exit 88 ;; esac
done
available=$(available_bytes)
if [ "$available" -lt "$FLOOR" ]; then echo "disk guard: ${available} bytes free, below $FLOOR" >&2; exit 85; fi
printf '%s\n' "launcher_pid=$$" "preexec_free_bytes=$available" "input_mode=$(stat -f '%Lp' "$RUN/input.nml")" "seed_mode=$(stat -f '%Lp' "$RUN/Ni-imp_out.nml")" "element_mode=$(stat -f '%Lp' "$RUN/Ni-imp.nml")" "all_working_files_writable=yes" "started_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$RUN/lifecycle.txt"
"$PYTHON" "$WATCHER" --run-id rd_ni_u3_bulk_rc80_block --solver-pid-file "$RUN/solver.pid" --run-dir "$RUN" --state-file Ni-imp_scf.nml --output-dir "$ROOT/captures" --ready-file "$RUN/watcher.ready" > "$RUN/watcher.log" 2>&1 &
watcher_pid=$!
printf '%s\n' "$watcher_pid" > "$RUN/watcher.pid"
tries=0
while [ ! -s "$RUN/watcher.ready" ]; do
  if ! kill -0 "$watcher_pid" 2>/dev/null; then wait "$watcher_pid" 2>/dev/null || true; echo 'collector exited before readiness' >&2; exit 82; fi
  tries=$((tries + 1)); if [ "$tries" -gt 1000 ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo 'collector readiness timeout' >&2; exit 83; fi
  sleep 0.1
done
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ] || [ "$(sha256 "$RUN/input.nml")" != "$EXPECTED_INPUT" ] || [ "$(sha256 "$RUN/Ni-imp_out.nml")" != "$EXPECTED_SEED" ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo 'pre-exec hash mismatch' >&2; exit 84; fi
available=$(available_bytes)
if [ "$available" -lt "$FLOOR" ]; then kill "$watcher_pid" 2>/dev/null || true; wait "$watcher_pid" 2>/dev/null || true; echo "disk guard: ${available} bytes free before exec" >&2; exit 85; fi
cd "$RUN"
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=1 \
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  RSLMTO_FIXED_U_TRACE="$RUN/fixed_u_trace.csv" \
  RSLMTO_FIXED_U_RUN_ID=rd_ni_u3_bulk_rc80_block \
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
printf '%s\n' "exit_code=$rc" "ended_local=$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "$RUN/lifecycle.txt"
kill "$guard_pid" 2>/dev/null || true
wait "$guard_pid" 2>/dev/null || true
wait "$watcher_pid" || true
exit "$rc"
