#!/usr/bin/env python3
"""LdaUKspaceLocalNsp3: hubbard_radial_moments = 'local' in k-space at nsp = 3 against the nsp = 1 route (two-atom bcc Fe cell, U = 2 eV,
J = 0, 8^3 k, 3 steps). For collinear moments the two spin representations describe the same state, so E_F and the site moments must
agree at every step; nsp = 1 (the second site's file with swapped spin channels) uses the scalar correction, nsp = 3 the spinor correction
along the site axis. Cases and the bug each one catches:
  af_local    antiferromagnet (mom = +z, -z at nsp = 3): the spinor correction applied in the wrong spin frame for the mom(3) < 0 site,
              or V_U, P or M1 of the two spin channels exchanged
  fm_local    ferromagnet: the spinor correction itself (no spin-frame subtlety)
  af_hoh      antiferromagnet with hoh (kspace_ham_order = 'second'): V_U additive in enim, H0 = H - V_U in the moments
  af_engaged  nsp = 3 'local' against nsp = 3 'spectral': the two must differ by more than 1e-3 Ry in E_F at some step, so the
              correction is acting (a skipped correction would make the first two cases pass trivially)
Not covered: nsp = 4, non-collinear axes (stop with a fatal), calctype other than B.

usage: test_kspace_local_nsp3.py --binary rslmto.x --scratch-root DIR
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
ENGAGED = 1.0e-3  # Ry


def numbers(binary: Path, scratch: Path, name: str, *, nsp: int, mode: str, af: bool, hoh: bool = False) -> tuple[list[float], list[float], list[float]] | str:
    mom2 = (0.0, 0.0, -1.0) if af else (0.0, 0.0, 1.0)
    work = t3.af_deck(scratch, name, nsp=nsp, mom2=mom2, swap2=(nsp == 1 and af), nstep=NSTEP, kspace=True, radial=mode, hoh=hoh)
    status, log = gh.run_binary(binary, work)
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if status != 0:
        return f"run failed\n{log[-1200:]}"
    ef = [float(x) for x in re.findall(r"N\(E_F=\s*(-?\d+\.\d+)", log)]
    m = {a: [float(x) for x in re.findall(rf"Spin moment of atom\s+{a} is\s+(\S+)", log)] for a in (1, 2)}
    return ef, m[1], m[2]


def parity(binary: Path, scratch: Path, name: str, *, af: bool, hoh: bool = False) -> tuple[bool, str]:
    ref = numbers(binary, scratch, f"{name}_nsp1", nsp=1, mode="local", af=af, hoh=hoh)
    res = numbers(binary, scratch, f"{name}_nsp3", nsp=3, mode="local", af=af, hoh=hoh)
    for label, x in (("nsp = 1", ref), ("nsp = 3", res)):
        if isinstance(x, str):
            return False, f"{label}: {x}"
    if any(len(x) != NSTEP for x in (*ref, *res)):
        return False, f"expected {NSTEP} steps, got {[len(x) for x in (*ref, *res)]}"
    d_ef = max(abs(a - b) for a, b in zip(ref[0], res[0]))
    d_m = max(abs(a - b) for r, g in zip(ref[1:], res[1:]) for a, b in zip(r, g))
    return d_ef < TOL and d_m < TOL, (f"nsp = 3 against nsp = 1 over {NSTEP} steps: max |dE_F| {d_ef:.1e} Ry, max |dm| {d_m:.1e} muB "
                                      f"(need < {TOL}); |m1|, |m2| at the last step {res[1][-1]:.6f}, {res[2][-1]:.6f}")


def engaged(binary: Path, scratch: Path) -> tuple[bool, str]:
    loc = numbers(binary, scratch, "eng_local", nsp=3, mode="local", af=True)
    spec = numbers(binary, scratch, "eng_spectral", nsp=3, mode="spectral", af=True)
    for label, x in (("local", loc), ("spectral", spec)):
        if isinstance(x, str):
            return False, f"{label}: {x}"
    d = max(abs(a - b) for a, b in zip(loc[0], spec[0]))
    return d > ENGAGED, f"nsp = 3 'local' against 'spectral': max |dE_F| {d:.1e} Ry (need > {ENGAGED})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, scratch = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, ok_info in (("af_local", parity(binary, scratch, "af", af=True)), ("fm_local", parity(binary, scratch, "fm", af=False)),
                          ("af_hoh", parity(binary, scratch, "afh", af=True, hoh=True)),
                          ("af_engaged", engaged(binary, scratch))):
        ok, info = ok_info
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
