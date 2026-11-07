#!/usr/bin/env python3
"""Fixed U at nsp = 3 and 4 (spinor U in the noncollinear Hamiltonians), checked for collinear +z/-z setups.

Checks (each prints PASS/FAIL by name; the test fails if any fails):
  shift_nsp3     bcc Fe ferromagnet, one real-space Chebyshev step, no hoh: E_F(U=2) - E_F(U=0) at nsp = 3 equals the
                 nsp = 1 value within 1e-4 Ry (the shift is above 0.05 Ry). The nsp = 1 route is collinear and has no spinor V.
  shift_nsp4_hoh the same with SOC and hoh at nsp = 4 against nsp = 2 (spinor U with V in enim).
  af_fe          two-atom Fe antiferromagnet (G-type, simple cubic cell), U = 2 eV, 2 steps, no hoh: nsp = 3 with the
                 noncollinear setup (the same potential file on both sites, mom = +z and -z) against nsp = 1 with the
                 spin-swapped file on site 2. E_F and the site-1 moment agree within 1e-5, and the two site moments are
                 opposite at nsp = 3.

Bugs caught: V_U dropped at nsp >= 3 (the shifts are 0); V_U built from the wrong spin channel for the site with
mom = -z (the AF run leaves the nsp = 1 result). Not caught: the spinor against the collinear V at nsp = 4 (they differ by
7e-5 Ry in E_F at step 2 and not at step 1).

The code rewrites fermi in input.nml after a run, so every deck is built fresh.

usage: test_nsp34_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_hoh_u as hu  # noqa: E402

SHIFT_TOL = 1.0e-4  # Ry, nsp = 3 against 1 and 4 against 2
AF_TOL = 1.0e-5


def swap_spin_channels(text: str) -> str:
    """Swap the spin channels of a potential file: name(:, s), ql(k, :, s), the two-entry xi_p, xi_d, rac, vrmax = (avg, diff)."""
    out = []
    for line in text.splitlines():
        m_ql = re.match(r"(\s*)(\w+)\((:|\d+), :, ([12])\)(.*)", line)
        m_ch = re.match(r"(\s*)(\w+)\(:, ([12])\)(.*)", line)
        m_two = re.match(r"(\s*)(xi_p|xi_d|rac|vrmax)\s*=\s*(\S+),\s*(\S+)\s*$", line)
        if m_ql:
            line = f"{m_ql.group(1)}{m_ql.group(2)}({m_ql.group(3)}, :, {3 - int(m_ql.group(4))}){m_ql.group(5)}"
        elif m_ch:
            line = f"{m_ch.group(1)}{m_ch.group(2)}(:, {3 - int(m_ch.group(3))}){m_ch.group(4)}"
        elif m_two:
            a, b = m_two.group(3), m_two.group(4)
            if m_two.group(2) == "vrmax":
                b = b[1:] if b.startswith("-") else "-" + b
            else:
                a, b = b, a
            line = f"{m_two.group(1)}{m_two.group(2)} = {a}, {b}"
        out.append(line)
    return "\n".join(out) + "\n"


def af_deck(scratch: Path, name: str, *, nsp: int) -> Path:
    """Two-atom Fe AF cell from the bcc Fe seed: nsp = 1 swaps the spin channels of site 2, nsp = 3 sets mom = -z on it."""
    work = scratch / name
    work.mkdir(parents=True)
    seed = (hu.BULK_FIXTURE / "Fe.nml").read_text()
    site2 = swap_spin_channels(seed) if nsp == 1 else re.sub(r"(?m)^ mom = .*$", " mom = 0.0, 0.0, -1.0", seed)
    (work / "Fe1.nml").write_text(seed.replace("symbol = 'Fe'", "symbol = 'Fe1'"))
    (work / "Fe2.nml").write_text(site2.replace("symbol = 'Fe'", "symbol = 'Fe2'"))
    (work / "lattice.nml").write_text(
        "&lattice\nnbulk_bulk = 2\nntot = 2\nnbas = 2\nnrec = 2\n"
        "a(:, 1) = 1.0, 0.0, 0.0\na(:, 2) = 0.0, 1.0, 0.0\na(:, 3) = 0.0, 0.0, 1.0\n"
        "crd(:, 1) = 0.0, 0.0, 0.0\ncrd(:, 2) = 0.5, 0.5, 0.5\n"
        "izp(1) = 1\nizp(2) = 2\nno(1) = 1\nno(2) = 2\niu(1) = 1\niu(2) = 2\nib(1) = 1\nib(2) = 2\nirec(1) = 1\nirec(2) = 2\n/\n")
    (work / "input.nml").write_text(
        "&calculation\npre_processing = 'bravais'\nverbose = T\n/\n"
        "&lattice\nrc = 040\nalat = 2.86120\ncrystal_sym = 'file'\nwav = 1.40880\nntype = 2\nct(1) = 3.0d0\nct(2) = 3.0d0\nr2 = 9.00d0\n/\n"
        "&atoms\ndatabase = './'\nlabel(1) = 'Fe1'\nlabel(2) = 'Fe2'\n/\n"
        "&self\nnstep = 2\n/\n"
        "&energy\nfermi = -0.042265\nenergy_min = -3.0\nenergy_max = 3.0\nchannels_ldos = 2500\n/\n"
        f"&control\ncalctype = 'B'\nnsp = {nsp}\nlld = 200\nrecur = 'chebyshev'\n/\n"
        "&mix\nbeta = 0.1\nmixtype = 'broyden'\n/\n"
        "&hamiltonian\nhoh = .false.\nhubbard_u_general(1,3) = 2.0\nhubbard_u_general(2,3) = 2.0\n/\n")
    return work


def af_numbers(log: str) -> tuple[list[float], list[float], list[float]]:
    ef = [float(x) for x in re.findall(r"Free Fermi energy:\s*(-?\d+\.\d+)", log)]
    mz = {1: [], 2: []}
    for atom, _, _, z in re.findall(r"Spin moment projections of atom\s+(\d)\s+is\s+(\S+)\s+(\S+)\s+(\S+)", log):
        mz[int(atom)].append(float(z))
    return ef, mz[1], mz[2]


def shift_nsp3(binary: Path, scratch: Path) -> tuple[bool, str]:
    ref = hu.shift(binary, scratch, "n1", hoh=False)
    new = hu.shift(binary, scratch, "n3", hoh=False, nsp=3)
    for r in (ref, new):
        if isinstance(r, str):
            return False, r
    ok = new > hu.SHIFT_MIN and abs(new - ref) < SHIFT_TOL
    return ok, f"E_F shift, U=2 eV, no hoh: nsp=3 {new:+.6f} Ry, nsp=1 {ref:+.6f} Ry (need > {hu.SHIFT_MIN} and |difference| < {SHIFT_TOL})"


def shift_nsp4_hoh(binary: Path, scratch: Path) -> tuple[bool, str]:
    ref = hu.shift(binary, scratch, "n2", hoh=True, fixture=hu.SOC_FIXTURE, nsp=2)
    new = hu.shift(binary, scratch, "n4", hoh=True, fixture=hu.SOC_FIXTURE, nsp=4)
    for r in (ref, new):
        if isinstance(r, str):
            return False, r
    ok = new > hu.SHIFT_MIN and abs(new - ref) < SHIFT_TOL
    return ok, f"E_F shift, U=2 eV, hoh, SOC: nsp=4 {new:+.6f} Ry, nsp=2 {ref:+.6f} Ry (need > {hu.SHIFT_MIN} and |difference| < {SHIFT_TOL})"


def af_fe(binary: Path, scratch: Path) -> tuple[bool, str]:
    res = {}
    for nsp in (1, 3):
        status, log = gh.run_binary(binary, af_deck(scratch, f"af_nsp{nsp}", nsp=nsp))
        if status != 0:
            return False, f"run failed (nsp={nsp})\n{log[-1500:]}"
        res[nsp] = af_numbers(log)
    (ef1, m11, _), (ef3, m13, m23) = res[1], res[3]
    if min(len(ef1), len(ef3), len(m11), len(m13), len(m23)) < 2:
        return False, "fewer than two SCF steps in the log"
    d_ef = max(abs(a - b) for a, b in zip(ef1, ef3))
    d_m = max(abs(abs(a) - abs(b)) for a, b in zip(m11, m13))
    opposite = max(abs(a + b) for a, b in zip(m13, m23))
    ok = d_ef < AF_TOL and d_m < AF_TOL and opposite < AF_TOL
    return ok, (f"nsp=3 against nsp=1: max |dE_F| {d_ef:.1e} Ry, max ||m1| difference| {d_m:.1e} muB, "
                f"max |m1 + m2| at nsp=3 {opposite:.1e} muB (need < {AF_TOL})")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    failed = False
    for name, check in (("shift_nsp3", shift_nsp3), ("shift_nsp4_hoh", shift_nsp4_hoh), ("af_fe", af_fe)):
        ok, info = check(args.binary.resolve(), args.scratch_root.resolve())
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
