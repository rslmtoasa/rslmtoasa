#!/usr/bin/env python3
"""LdaUKspaceGbt: k-space gbt_single_q at q = 0 against periodic_nc, with U (two-atom bcc Fe cell, nsp = 3, 8^3 k, 3 steps). At q = 0 and
theta = 0 the two magnetic representations build the same operator, so E_F and both site moments must agree at every step. The two
routes assemble H(k) in different routines (reciprocal_fourier: from the real-space blocks ee; reciprocal_gbt: from structure constants and
potential parameters). Cases and the bug each one catches:
  fm_second   ferromagnet, hoh (kspace_ham_order = 'second'): V_U missing from the onsite enu, or added to h where it would enter h o h
  af_second   antiferromagnet (the same potential file on both sites, mom = +z, -z), hoh: V_U swapped between the spin channels of the
              mom(3) < 0 site together with E_nu (V_U is built from the rotating-frame P and is already in the majority-first order);
              the two sublattices then stop being equivalent after step 1
  af_first    antiferromagnet, first order: the same swap, with V_U added to the onsite block of h
  af_ccor     antiferromagnet, hoh + ccor_2c (ccor_vmt_mode = vmad_scalar; the Fe seed has no vrmax): V_U with the CCOR term present
Not covered: q != 0, constraining fields with U, kspace_ham_order = 'orthogonal' (stops with a fatal when U is set), real-space
gbt_single_q with a mom(3) < 0 site (KNOWN_ISSUES.md, F1; see test_gbt_hoh_u.py for the ferromagnet).

usage: test_kspace_gbt_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL = 2.0e-6  # Ry and muB: the logs print six decimals
NSTEP = 3


def deck(scratch: Path, name: str, rep: str, *, af: bool, order: str, ccor: bool) -> Path:
    work = t3.af_deck(scratch, name, nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, -1.0 if af else 1.0), hoh=order == "second",
                      nstep=NSTEP, kspace=True)
    path = work / "input.nml"
    text = path.read_text()
    text = text.replace("&control\n", "&control\ndensity_policy = 'constrained_spiral'\n", 1)
    sdot = ".true.\nscreening = 'fitted'" if ccor else ".false."
    text = text.replace("strux_backend = 'strux_lib'\n", f"strux_backend = 'strux_lib'\nstrux_want_sdot = {sdot}\n", 1)
    ham = f"magnetic_representation = '{rep}'\n" + ("ccor_2c = .true.\nccor_strict = .true.\nccor_vmt_mode = 'vmad_scalar'\n" if ccor else "")
    text = text.replace("&hamiltonian\n", "&hamiltonian\n" + ham, 1)
    path.write_text(text.replace("&reciprocal\n", f"&reciprocal\nkspace_ham_order = '{order}'\n", 1))
    return work


def numbers(binary: Path, work: Path) -> tuple[list[float], list[float], list[float]] | str:
    status, log = gh.run_binary(binary, work)
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if status != 0:
        return f"run failed\n{log[-1200:]}"
    ef = [float(x) for x in re.findall(r"N\(E_F=\s*(-?\d+\.\d+)", log)]
    m = {a: [float(x) for x in re.findall(rf"Spin moment of atom\s+{a} is\s+(\S+)", log)] for a in (1, 2)}
    return ef, m[1], m[2]


def case(binary: Path, scratch: Path, name: str, *, af: bool, order: str, ccor: bool) -> tuple[bool, str]:
    ref = numbers(binary, deck(scratch, f"{name}_pnc", "periodic_nc", af=af, order=order, ccor=ccor))
    gbt = numbers(binary, deck(scratch, f"{name}_gbt", "gbt_single_q", af=af, order=order, ccor=ccor))
    for label, res in (("periodic_nc", ref), ("gbt_single_q", gbt)):
        if isinstance(res, str):
            return False, f"{label}: {res}"
    if any(len(x) != NSTEP for x in (*ref, *gbt)):
        return False, f"expected {NSTEP} steps, got {[len(x) for x in (*ref, *gbt)]} (E_F, m1, m2 of periodic_nc, then gbt_single_q)"
    d_ef = max(abs(a - b) for a, b in zip(ref[0], gbt[0]))
    d_m = max(abs(a - b) for r, g in zip(ref[1:], gbt[1:]) for a, b in zip(r, g))
    ok = d_ef < TOL and d_m < TOL
    return ok, (f"gbt_single_q against periodic_nc over {NSTEP} steps: max |dE_F| {d_ef:.1e} Ry, max |dm| {d_m:.1e} muB (need < {TOL}); "
                f"|m1|, |m2| at the last step {gbt[1][-1]:.6f}, {gbt[2][-1]:.6f}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    failed = False
    for name, af, order, ccor in (("fm_second", False, "second", False), ("af_second", True, "second", False),
                                  ("af_first", True, "first", False), ("af_ccor", True, "second", True)):
        ok, info = case(args.binary.resolve(), args.scratch_root.resolve(), name, af=af, order=order, ccor=ccor)
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
