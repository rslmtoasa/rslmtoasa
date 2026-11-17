#!/usr/bin/env python3
"""LdaUGbtHoh: real-space gbt_single_q at q = 0 against periodic_nc, with U and hoh (two-atom ferromagnetic bcc Fe cell, nsp = 3,
block recursion, 3 steps). At q = 0 and theta = 0 the two magnetic representations build the same operator, so E_F and the site
moment must agree at every step. The reference is the periodic_nc route, which builds its Hamiltonian in a different routine
(build_bulkham; gbt_single_q uses build_gbt_bulkham). Cases and the bug each one catches:
  hoh_u0       hoh, U = 0:        the onsite term of h taken from cx instead of cex (the collinear builder uses cex and adds enim = cx - cex)
  hoh_u2       hoh, U = 2 eV:     V_U added to ee, where it also enters h o h, instead of to enim (H = h - h o h + enim + V_U)
  hoh_ccor_u2  hoh + ccor_2c, U = 2: the same, with the CCOR term present (ccor_vmt_mode = vmad_scalar; the Fe seed has no vrmax)
  first_u2     no hoh, U = 2 eV:  V_U in ee at first order (path unchanged; checks the deck and the comparison itself)
Not covered: antiferromagnets (a site with mom(3) < 0 differs between the two routes even at U = 0 without hoh, KNOWN_ISSUES.md, F1),
k-space gbt_single_q, q != 0.

usage: test_gbt_hoh_u.py --binary rslmto.x --scratch-root DIR
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


def deck(scratch: Path, name: str, rep: str, *, hoh: bool, ccor: bool, u_ev: float, nstep: int) -> Path:
    work = t3.af_deck(scratch, name, nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, 1.0), hoh=hoh, nstep=nstep, recur="block")
    path = work / "input.nml"
    text = path.read_text()
    sdot = ".true.\nscreening = 'fitted'" if ccor else ".false."
    text = text.replace("rc = 040", "rc = 040\npbc = .true.\nb1 = .true.\nb2 = .true.\nb3 = .true.\nn1 = 20\nn2 = 20\nn3 = 20\n"
                        f"strux_backend = 'strux_lib'\nstrux_want_sdot = {sdot}", 1)
    text = text.replace("&control\n", "&control\ndensity_policy = 'constrained_spiral'\n", 1)
    text = re.sub(r"hubbard_u_general\((\d),3\) = [\d.]+", lambda m: f"hubbard_u_general({m.group(1)},3) = {u_ev}", text)
    ham = f"magnetic_representation = '{rep}'\n" + ("ccor_2c = .true.\nccor_strict = .true.\nccor_vmt_mode = 'vmad_scalar'\n" if ccor else "")
    path.write_text(text.replace("&hamiltonian\n", "&hamiltonian\n" + ham, 1))
    return work


def numbers(binary: Path, work: Path) -> tuple[list[float], list[float]] | str:
    status, log = gh.run_binary(binary, work)
    if status != 0:
        return f"run failed\n{log[-1200:]}"
    ef, mom = t3.af_numbers(re.sub(r"\x1b\[[0-9;]*m", "", log))
    return ef, [m[2] for m in mom[1]]


def case(binary: Path, scratch: Path, name: str, *, hoh: bool, ccor: bool, u_ev: float, nstep: int) -> tuple[bool, str]:
    ref = numbers(binary, deck(scratch, f"{name}_pnc", "periodic_nc", hoh=hoh, ccor=ccor, u_ev=u_ev, nstep=nstep))
    gbt = numbers(binary, deck(scratch, f"{name}_gbt", "gbt_single_q", hoh=hoh, ccor=ccor, u_ev=u_ev, nstep=nstep))
    for label, res in (("periodic_nc", ref), ("gbt_single_q", gbt)):
        if isinstance(res, str):
            return False, f"{label}: {res}"
    if len(ref[0]) != nstep or len(gbt[0]) != nstep or len(ref[1]) != nstep or len(gbt[1]) != nstep:
        return False, f"expected {nstep} steps, got E_F {len(ref[0])}/{len(gbt[0])}, moments {len(ref[1])}/{len(gbt[1])}"
    d_ef = max(abs(a - b) for a, b in zip(ref[0], gbt[0]))
    d_m = max(abs(a - b) for a, b in zip(ref[1], gbt[1]))
    ok = d_ef < TOL and d_m < TOL
    return ok, (f"gbt_single_q against periodic_nc over {nstep} steps: max |dE_F| {d_ef:.1e} Ry, max |dm| {d_m:.1e} muB (need < {TOL}); "
                f"E_F {', '.join(f'{x:.6f}' for x in gbt[0])}; m {', '.join(f'{x:.6f}' for x in gbt[1])}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    failed = False
    for name, hoh, ccor, u_ev, nstep in (("hoh_u0", True, False, 0.0, 1), ("hoh_u2", True, False, 2.0, 3),
                                          ("hoh_ccor_u2", True, True, 2.0, 3), ("first_u2", False, False, 2.0, 3)):
        ok, info = case(args.binary.resolve(), args.scratch_root.resolve(), name, hoh=hoh, ccor=ccor, u_ev=u_ev, nstep=nstep)
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
