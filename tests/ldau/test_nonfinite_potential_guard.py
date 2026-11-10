#!/usr/bin/env python3
"""LdaUNonfinitePotentialGuard: a k-space run whose potential parameters are NaN stops with a message naming the atom, shell and field, before the eigensolver.
Cause of the guard: a failed atomic radial solve (RSEQSR; log line 'RSEQ NIT GT 400 AND BAD NODES') leaves NaN potential parameters, and the run used to stop
in 'reciprocal backend: LAPACK eigensolver failed' (KNOWN_ISSUES.md, spectral radial moments under PBE).
  clean   bcc Fe k-space deck, one step, unchanged: exit 0 (the guard does not fire on finite parameters).
  nan_c   the same deck with c(0,1) = NaN in Fe.nml: non-zero exit and the log contains 'B1 non-finite guard: atom=Fe' and 'potential parameter C'.
          Bug caught: the guard removed or not reached (the run then stops elsewhere with another message, or continues with NaN).

usage: test_nonfinite_potential_guard.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    src = rr.CASES / "k_space_scf/bccFe"
    deck = lambda t: rr.use_kspace(rr.set_nstep(t, 1))
    failed = False

    clean = rr.prepare(src, base / "clean", deck)
    code = rr.run(binary, clean)
    ok = code == 0
    print(f"{'PASS' if ok else 'FAIL'} clean exit code {code} (need 0)")
    failed |= not ok

    bad = rr.prepare(src, base / "nan_c", deck)
    par = bad / "Fe.nml"
    text = par.read_text()
    new = re.sub(r"(?m)^( c\(:, 1\) = )[^,]*,", r"\1NaN,", text, count=1)
    if new == text:
        print("FAIL nan_c could not edit c(:, 1) in Fe.nml")
        return 1
    par.write_text(new)
    code = rr.run(binary, bad)
    log = " ".join(rr.log_text(bad).split())
    ok = code != 0 and "B1 non-finite guard: atom=Fe" in log and "potential parameter C" in log
    print(f"{'PASS' if ok else 'FAIL'} nan_c exit code {code} (need non-zero), guard message {'found' if ok else 'not found'}")
    failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
