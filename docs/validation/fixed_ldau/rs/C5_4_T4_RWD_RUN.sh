#!/bin/bash
set -euo pipefail
run=$1 bin=$2 expected_bin=$3 ledger=$4 watcher=$5 disk_guard=$6 runid=$7 record=$8
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
test -d "$run" && test -w "$run" && test -x "$bin" && test -f "$watcher" && test -f "$disk_guard"
test "$(sha "$bin")" = "$expected_bin" || { echo binary_hash_mismatch >&2; exit 79; }
test "$(find "$run" -maxdepth 1 -type f | wc -l | tr -d ' ')" = 5 || { echo run_directory_not_pristine >&2; exit 78; }
for spec in 'input.nml:1f63c85f66c328811d9c0ef99f3e82454491ec40dbffcffdff455ec48e26e14c' 'Ni-bulk.nml:7385d1f467c336795f71bfa1a8dff292cd3ec500bc00a9182324c060ab2563ac' 'Ni-bulk_out.nml:736a3c4e8fd2e1cdc7550502c4788e1b9f0a031c06a44be2ee77a8fd5b1ac1c4' 'Ni-imp.nml:4824dcfb4e01db54fa0ca1e4fe34b31cc398e05382463807504c21b3c5c18d57' 'Ni-imp_out.nml:0213bf83b07bd82446f293bf19be8ba37b0a2bc816011d29f2c7e8d5c483738a'; do
  f=${spec%%:*}; want=${spec#*:}; test "$(sha "$run/$f")" = "$want" || { echo "input_hash_mismatch $f" >&2; exit 80; }
done
free=$(df -Pk "$run" | awk 'NR==2 {print $4 * 1024}'); floor=1073741824
test "$free" -ge "$floor" || { echo "disk_below_floor $free" >&2; exit 85; }
for key in OMP_NUM_THREADS OMP_DYNAMIC RSLMTO_OMP_THREADS_SERIAL OPENBLAS_NUM_THREADS MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS NUMEXPR_NUM_THREADS; do
  case "$key" in OMP_DYNAMIC) export "$key=FALSE" ;; *) export "$key=1" ;; esac
done
python3 - "$ledger" "$record" <<'PY'
import json,sys
with open(sys.argv[1],"a",encoding="utf-8") as f:
 f.write(json.dumps(json.loads(sys.argv[2]),sort_keys=True,separators=(",",":"))+"\n"); f.flush()
PY
mkdir -p "$run/captures"; : > "$run/fixed_u_trace.csv"
"$(command -v python3)" "$watcher" --run-id "$runid" --solver-pid-file "$run/solver.pid" --run-dir "$run" --state-file Ni-imp_scf.nml --output-dir "$run/captures" --ready-file "$run/watcher.ready" > "$run/watcher.log" 2>&1 &
wpid=$!; tries=0
while [ ! -s "$run/watcher.ready" ]; do
  if ! kill -0 "$wpid" 2>/dev/null; then wait "$wpid" || true; echo watcher_not_ready >&2; exit 82; fi
  tries=$((tries+1)); if [ "$tries" -gt 1000 ]; then kill "$wpid" 2>/dev/null || true; wait "$wpid" || true; exit 83; fi
  sleep .1
done
cd "$run"; ulimit -c 0
export RSLMTO_FIXED_U_TRACE="$run/fixed_u_trace.csv" RSLMTO_FIXED_U_RUN_ID="$runid"
export RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=0 RSLMTO_FIXED_U_TRACE_QL_ONLY=0
export RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0 RSLMTO_FIXED_U_A2_SCAN=0
"$bin" > testrun.log 2>&1 &
pid=$!; printf '%s\n' "$pid" > solver.pid
printf '%s\n' "solver_pid=$pid" "watcher_pid=$wpid" "OMP_NUM_THREADS=1" "RSLMTO_OMP_THREADS_SERIAL=1" "free_bytes=$free" > lifecycle.txt
bash "$disk_guard" "$run" "$pid" > disk_guard.log 2>&1 &
gpid=$!; printf '%s\n' "$gpid" > disk_guard.pid
rc=0; bounded_stop=0
while kill -0 "$pid" 2>/dev/null; do
  loops=$(rg -c 'Perform recursion at step' testrun.log 2>/dev/null || true); loops=${loops:-0}
  if [ "$loops" -ge 37 ]; then kill -TERM "$pid" 2>/dev/null || true; bounded_stop=1; break; fi
  sleep 1
done
wait "$pid" || rc=$?
printf '%s\n' "$rc" > exit_code.txt
printf '%s\n' "exit_code=$rc" "bounded_after_36_completed_loops=$bounded_stop" >> lifecycle.txt
kill "$gpid" 2>/dev/null || true; wait "$gpid" 2>/dev/null || true
wait "$wpid" || true
exit "$rc"
