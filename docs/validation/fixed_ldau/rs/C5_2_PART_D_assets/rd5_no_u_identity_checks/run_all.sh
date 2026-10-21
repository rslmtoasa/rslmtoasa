#!/bin/bash
set -u
set -o pipefail
ROOT=/private/tmp/c52-rd5-no-u
BASE=/private/tmp/c5-clean8b-release-build/bin/rslmto.x
REV=/private/tmp/c5-rd1-frozen/4bb76cc2bbd251d95cdc434eec0348b9082069af2deeb4f4859f6d8056779a4f/rslmto.x
RS=/private/tmp/c5-2-partA/a1_bccFe/source_snapshot
KS=/private/tmp/c5-3-b2-prep/identity_ks/clean8b_baseline
run_one() {
  name="$1"; snap="$2"; bin="$3"; expected="$4"
  d="$ROOT/$name"
  test ! -e "$d" || { echo "refusing existing run dir $d"; return 90; }
  mkdir "$d" || return 91
  cp "$snap/Fe.nml" "$snap/input.nml" "$d/" || return 92
  (cd "$d" && env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE RSLMTO_OMP_THREADS_SERIAL=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 "$bin" > testrun.log 2>&1) &
  p=$!
  printf '%s\n' "$p" > "$d/solver_pid.txt"
  printf 'started %s pid=%s binary=%s\n' "$name" "$p" "$bin" | tee -a "$ROOT/lifecycle.txt"
  wait "$p"; rc=$?
  printf '%s\n' "$rc" > "$d/exit_code.txt"
  printf 'finished %s pid=%s exit=%s\n' "$name" "$p" "$rc" | tee -a "$ROOT/lifecycle.txt"
  printf '%s\n' "$rc" > "$d/return_code.txt"
  if [ "$rc" -ne 0 ]; then return "$rc"; fi
  [ "$expected" = x ] || return 0
}
run_one rs_clean8b "$RS" "$BASE" x || exit $?
run_one rs_rd1 "$RS" "$REV" x || exit $?
run_one ks_clean8b "$KS" "$BASE" x || exit $?
run_one ks_rd1 "$KS" "$REV" x || exit $?
