#!/bin/bash
set -euo pipefail
run=$1 bin=$2 ledger=$3 watcher=$4 captures=$5 runid=$6 record=$7
test -d "$run" && test -w "$run" && test -x "$bin" && test -f "$watcher"
test ! -e "$run/solver.pid" && test ! -e "$run/exit_code.txt"
free=$(df -Pk "$run" | awk 'NR==2 {print $4 * 1024}')
floor=1073741824
if [ "$free" -lt "$floor" ]; then echo "disk below floor: $free" >&2; exit 85; fi
for key in OMP_NUM_THREADS OMP_DYNAMIC RSLMTO_OMP_THREADS_SERIAL OPENBLAS_NUM_THREADS MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS NUMEXPR_NUM_THREADS; do
  case "$key" in OMP_DYNAMIC) export "$key=FALSE" ;; *) export "$key=1" ;; esac
done
python3 - "$ledger" "$record" <<'PY'
import json,sys
with open(sys.argv[1],"a",encoding="utf-8") as f:
 f.write(json.dumps(json.loads(sys.argv[2]),sort_keys=True,separators=(",",":"))+"\n"); f.flush()
PY
mkdir -p "$captures"
"$(command -v python3)" "$watcher" --run-id "$runid" --solver-pid-file "$run/solver.pid" --run-dir "$run" --state-file Ni-imp_scf.nml --output-dir "$captures" --ready-file "$run/watcher.ready" > "$run/watcher.log" 2>&1 &
wpid=$!
tries=0
while [ ! -s "$run/watcher.ready" ]; do
  if ! kill -0 "$wpid" 2>/dev/null; then wait "$wpid" || true; echo 'watcher failed before readiness' >&2; exit 82; fi
  tries=$((tries+1)); if [ "$tries" -gt 1000 ]; then kill "$wpid" 2>/dev/null || true; wait "$wpid" || true; exit 83; fi
  sleep .1
done
cd "$run"
export RSLMTO_FIXED_U_TRACE="$run/fixed_u_trace.csv" RSLMTO_FIXED_U_RUN_ID="$runid"
export RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=1 RSLMTO_FIXED_U_TRACE_QL_ONLY=0
export RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0 RSLMTO_FIXED_U_A2_SCAN=0
"$bin" > testrun.log 2>&1 &
pid=$!; printf '%s\n' "$pid" > solver.pid
printf '%s\n' "solver_pid=$pid" "watcher_pid=$wpid" "env_OMP_NUM_THREADS=1" "free_bytes=$free" > lifecycle.txt
rc=0; wait "$pid" || rc=$?
printf '%s\n' "$rc" > exit_code.txt
kill "$wpid" 2>/dev/null || true; wait "$wpid" || true
exit "$rc"
