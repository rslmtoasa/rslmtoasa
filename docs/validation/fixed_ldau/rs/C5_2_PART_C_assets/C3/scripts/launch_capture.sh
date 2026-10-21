#!/bin/bash
set -eu
if [ "$#" -ne 2 ]; then echo 'usage: launch_capture.sh <accepted-C2-endpoint-lock.json> <fresh-capture-directory>' >&2; exit 64; fi
LOCK_PATH="$1"
CAPTURE_DIR="$2"
TEMPLATE=/private/tmp/c52-c2-stage/c3_debug_capture_template
BIN=/private/tmp/c52-c2-stage/crash_triage_preflight/frozen/debug_bounds/b8b00eff0f1a4c3b0536006f7593c53e229cac7f0c77cbfed29126125af26b69/rslmto.x
EXPECTED_BIN=b8b00eff0f1a4c3b0536006f7593c53e229cac7f0c77cbfed29126125af26b69
sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
if [ "${C52_EXPLICIT_RUN_AUTHORIZED:-0}" != 1 ]; then echo 'refusing launch: explicit root authorization env is unset' >&2; exit 77; fi
if [ ! -f "$LOCK_PATH" ]; then echo 'missing explicit accepted endpoint lock' >&2; exit 78; fi
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ]; then echo 'frozen binary hash mismatch' >&2; exit 79; fi
if [ -e "$CAPTURE_DIR" ] && [ "$(ls -A "$CAPTURE_DIR" 2>/dev/null)" != '' ]; then echo 'capture directory must be fresh and empty' >&2; exit 80; fi
mkdir -p "$CAPTURE_DIR"
/opt/homebrew/opt/python@3.11/bin/python3.11 "$TEMPLATE/bind_accepted_endpoint.py" --lock "$LOCK_PATH" --capture-dir "$CAPTURE_DIR"
cp "$TEMPLATE/input.nml" "$CAPTURE_DIR/input.nml"
if [ "$(sha256 "$BIN")" != "$EXPECTED_BIN" ]; then echo 'pre-exec frozen binary hash mismatch' >&2; exit 81; fi
if [ "$(sha256 "$CAPTURE_DIR/input.nml")" != "$(sha256 "$TEMPLATE/input.nml")" ]; then echo 'capture deck hash mismatch' >&2; exit 82; fi
if [ ! -w "$CAPTURE_DIR/input.nml" ] || [ ! -w "$CAPTURE_DIR/Fe_out.nml" ]; then echo 'working input/restart must be writable' >&2; exit 83; fi
printf '%s\n' "$$" > "$CAPTURE_DIR/wrapper.pid"
cd "$CAPTURE_DIR"
env OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=4 \
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  RSLMTO_FIXED_U_TRACE="$CAPTURE_DIR/records.csv" \
  RSLMTO_FIXED_U_RUN_ID=c3_fe_u2_cheb_fullcapture_i1 \
  RSLMTO_FIXED_U_TRACE_ITERATION=1 \
  RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY=1 \
  RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY=0 \
  RSLMTO_FIXED_U_TRACE_QL_ONLY=0 \
  "$BIN" > capture.stdout.log 2>&1 &
solver_pid=$!
printf '%s\n' "$solver_pid" > "$CAPTURE_DIR/solver.pid"
rc=0
wait "$solver_pid" || rc=$?
printf '%s\n' "$rc" > "$CAPTURE_DIR/exit_code.txt"
exit "$rc"
