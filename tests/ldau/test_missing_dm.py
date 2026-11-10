#!/usr/bin/env python3
"""LdaUMissingDm: a fixed-U restart whose potential file names a density-matrix file that is missing continues instead of stopping (bcc Fe,
nsp = 1, U = 2 eV, J = 0, block recursion, rc = 40). A short run writes Fe_out.nml and Fe_dm.nml; two restarts of one step follow:
  with_dm      Fe_out.nml + Fe_dm.nml: the restart reads the density matrix (control).
  missing_dm   Fe_out.nml alone (the file named by `dm_file` is removed): the run must finish (exit 0), log the warning 'Referenced
               density-matrix file is missing ... continuing without it', give finite numbers, and differ from the control by more than
               1e-6 muB in the first-step moment (the density matrix is used when present, so the two starts are different states).
Bug caught: the fatal on a missing sidecar file (the run stops), or a fallback that silently reads stale data (no difference to the control).
Not covered: the orbital state reached afterwards (a restart without the density matrix loses the orbital occupation pattern),
impurity hosts (their provenance checks are separate).

usage: test_missing_dm.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_hoh_u as hu  # noqa: E402


def run(binary: Path, work: Path) -> tuple[int, str]:
    status, log = gh.run_binary(binary, work)
    return status, re.sub(r"\x1b\[[0-9;]*m", "", log)


def deck(scratch: Path, name: str, nstep: int) -> Path:
    work = gh.fresh(hu.BULK_FIXTURE, scratch, name)
    gh.set_line(work / "input.nml", r"^nstep\s*=.*$", f"nstep = {nstep}")
    gh.set_line(work / "input.nml", r"^rc\s*=.*$", "rc = 040")
    gh.set_line(work / "input.nml", r"^hubbard_u_general\(1,3\)\s*=.*$", "hubbard_u_general(1,3) = 2.0")
    return work


def first_moment(log: str) -> float | None:
    m = re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", log)
    return float(m[0]) if m else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, scratch = args.binary.resolve(), args.scratch_root.resolve()

    base = deck(scratch, "base", 2)
    status, log = run(binary, base)
    if status != 0 or not (base / "Fe_dm.nml").exists():
        print(f"FAIL setup the base run did not finish or wrote no density-matrix file\n{log[-1200:]}")
        return 1
    results = {}
    for name, keep_dm in (("with_dm", True), ("missing_dm", False)):
        work = deck(scratch, name, 1)
        shutil.copy2(base / "Fe_out.nml", work / "Fe.nml")
        if keep_dm:
            shutil.copy2(base / "Fe_dm.nml", work / "Fe_dm.nml")
        results[name] = run(binary, work)

    failed = False
    st, log = results["with_dm"]
    ok = st == 0 and first_moment(log) is not None
    print(f"{'PASS' if ok else 'FAIL'} with_dm restart with the density matrix: exit {st}, first-step moment {first_moment(log)}")
    failed |= not ok
    st, log = results["missing_dm"]
    mom, ref = first_moment(log), first_moment(results["with_dm"][1])
    warned = "Referenced density-matrix file is missing" in log and "continuing without it" in log
    finite = mom is not None and math.isfinite(mom)
    differs = finite and ref is not None and abs(mom - ref) > 1.0e-6
    ok = st == 0 and warned and finite and differs
    print(f"{'PASS' if ok else 'FAIL'} missing_dm exit {st}, warning logged {warned}, first-step moment {mom} against {ref} with the file "
          f"(need exit 0, the warning, a finite value, and a difference > 1e-6)")
    failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
