#!/bin/bash
set -euo pipefail
run=$1
binary=$2
ledger=$3
record=$4
runid=$5
trace=${6:-}
test -d "$run" && test -w "$run"
test -x "$binary"
test ! -e "$run/solver.pid" && test ! -e "$run/exit_code.txt"
free=$(df -Pk "$run" | awk 'NR==2 {print $4 * 1024}')
floor=1073741824
if [ "$free" -lt "$floor" ]; then echo "disk guard: $free bytes free, below $floor" >&2; exit 85; fi
for key in OMP_NUM_THREADS OMP_DYNAMIC RSLMTO_OMP_THREADS_SERIAL OPENBLAS_NUM_THREADS MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS NUMEXPR_NUM_THREADS; do
  case "$key" in OMP_DYNAMIC) export "$key=FALSE" ;; *) export "$key=1" ;; esac
done
python3 - "$ledger" "$record" <<'PY'
import json,sys
with open(sys.argv[1],"a",encoding="utf-8") as f:
    f.write(json.dumps(json.loads(sys.argv[2]),sort_keys=True,separators=(",",":"))+"\n")
    f.flush()
PY
cd "$run"
if [ -n "$trace" ]; then
  export RSLMTO_FIXED_U_TRACE="$run/fixed_u_trace.csv"
  export RSLMTO_FIXED_U_RUN_ID="$runid"
  export RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=1 RSLMTO_FIXED_U_TRACE_QL_ONLY=0
  export RSLMTO_FIXED_U_TRACE_ITERATION=0 RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=0
  export RSLMTO_FIXED_U_A2_SCAN=0
fi
"$binary" > testrun.log 2>&1 &
pid=$!
printf '%s\n' "$pid" > solver.pid
printf '%s\n' "run_id=$runid" "pid=$pid" "binary=$binary" "free_bytes=$free" > lifecycle.txt
rc=0
wait "$pid" || rc=$?
printf '%s\n' "$rc" > exit_code.txt
printf '%s\n' "exit_code=$rc" >> lifecycle.txt
exit "$rc"
