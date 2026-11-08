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
  af_fe_local    the same pair with hubbard_radial_moments = 'local': nsp = 3 (spinor correction along each site axis)
                 against nsp = 1 (the scalar route, validated against Quantum ESPRESSO for MnO).

Bugs caught: V_U dropped at nsp >= 3 (the shifts are 0); V_U built from the wrong spin channel for the site with
mom = -z (the AF run leaves the nsp = 1 result); the local radial-moment correction taken in the global instead of the
site spin frame for the site with mom = -z (af_fe_local leaves the nsp = 1 result). Not caught: the spinor against the collinear V at nsp = 4 (they differ by
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


def moment_line(text: str, mom: tuple[float, float, float]) -> str:
    return re.sub(r"(?m)^ mom = .*$", f" mom = {mom[0]!r}, {mom[1]!r}, {mom[2]!r}", text)


def af_deck(scratch: Path, name: str, *, nsp: int, mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, -1.0), swap2: bool | None = None,
            hoh: bool = False, nstep: int = 2, recur: str = "chebyshev", kspace: bool = False, self_extra: str = "",
            radial: str = "auto") -> Path:
    """Two-atom Fe cell (simple cubic, sites at 0 and (1/2,1/2,1/2)) from the bcc Fe seed, U = 2 eV on both sites.
    Site moments mom1, mom2 (the noncollinear convention: the same potential file on both sites); swap2 swaps the spin
    channels of the site-2 file (the nsp = 1 convention, default at nsp = 1, where mom is ignored)."""
    work = scratch / name
    work.mkdir(parents=True)
    seed = (hu.BULK_FIXTURE / "Fe.nml").read_text()
    if swap2 is None:
        swap2 = nsp == 1
    if nsp == 1:  # the spin lives in the channels; mom stays +z
        mom1 = mom2 = (0.0, 0.0, 1.0)
    site2 = swap_spin_channels(seed) if swap2 else seed
    (work / "Fe1.nml").write_text(moment_line(seed, mom1).replace("symbol = 'Fe'", "symbol = 'Fe1'"))
    (work / "Fe2.nml").write_text(moment_line(site2, mom2).replace("symbol = 'Fe'", "symbol = 'Fe2'"))
    periodic = "pbc = .true.\nb1 = .true.\nb2 = .true.\nb3 = .true.\nn1 = 20\nn2 = 20\nn3 = 20\nstrux_backend = 'strux_lib'\n" if kspace else ""
    reciprocal = ""
    if kspace:
        text = (gh.ROOT / "tests/scf/cases/k_space_scf/bccFe/input.nml").read_text()
        block = re.search(r"(?ms)^&reciprocal\b.*?^/\s*$", text).group(0)
        reciprocal = re.sub(r"(?m)^nk([123])\s*=.*$", r"nk\1 = 8", block) + "\n"
    (work / "lattice.nml").write_text(
        "&lattice\nnbulk_bulk = 2\nntot = 2\nnbas = 2\nnrec = 2\n"
        "a(:, 1) = 1.0, 0.0, 0.0\na(:, 2) = 0.0, 1.0, 0.0\na(:, 3) = 0.0, 0.0, 1.0\n"
        "crd(:, 1) = 0.0, 0.0, 0.0\ncrd(:, 2) = 0.5, 0.5, 0.5\n"
        "izp(1) = 1\nizp(2) = 2\nno(1) = 1\nno(2) = 2\niu(1) = 1\niu(2) = 2\nib(1) = 1\nib(2) = 2\nirec(1) = 1\nirec(2) = 2\n/\n")
    (work / "input.nml").write_text(
        "&calculation\npre_processing = 'bravais'\nverbose = T\n/\n"
        f"&lattice\nrc = 040\nalat = 2.86120\ncrystal_sym = 'file'\nwav = 1.40880\nntype = 2\nct(1) = 3.0d0\nct(2) = 3.0d0\nr2 = 9.00d0\n{periodic}/\n"
        "&atoms\ndatabase = './'\nlabel(1) = 'Fe1'\nlabel(2) = 'Fe2'\n/\n"
        f"&self\nnstep = {nstep}\n{'use_kspace = .true.' + chr(10) if kspace else ''}{self_extra}/\n"
        "&energy\nfermi = -0.042265\nenergy_min = -3.0\nenergy_max = 3.0\nchannels_ldos = 2500\n/\n"
        f"&control\ncalctype = 'B'\nnsp = {nsp}\nlld = {200 if recur == 'chebyshev' else 21}\nrecur = '{recur}'\n/\n"
        "&mix\nbeta = 0.1\nmixtype = 'broyden'\n/\n"
        f"{reciprocal}"
        f"&hamiltonian\nhoh = {'.true.' if hoh else '.false.'}\nhubbard_u_general(1,3) = 2.0\nhubbard_u_general(2,3) = 2.0\n"
        f"{'' if radial == 'auto' else 'hubbard_radial_moments = ' + repr(radial) + chr(10)}/\n")
    return work


def af_numbers(log: str) -> tuple[list[float], dict[int, list[tuple[float, float, float]]]]:
    """E_F per step and the moment vector of each site per step."""
    ef = [float(x) for x in re.findall(r"Free Fermi energy:\s*(-?\d+\.\d+)", log)]
    mom: dict[int, list[tuple[float, float, float]]] = {1: [], 2: []}
    for atom, x, y, z in re.findall(r"Spin moment projections of atom\s+(\d)\s+is\s+(\S+)\s+(\S+)\s+(\S+)", log):
        mom[int(atom)].append((float(x), float(y), float(z)))
    return ef, mom


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


def af_fe(binary: Path, scratch: Path, radial: str = "auto") -> tuple[bool, str]:
    res = {}
    for nsp in (1, 3):
        status, log = gh.run_binary(binary, af_deck(scratch, f"af_{radial}_nsp{nsp}", nsp=nsp, radial=radial))
        if status != 0:
            return False, f"run failed (nsp={nsp})\n{log[-1500:]}"
        res[nsp] = af_numbers(log)
    (ef1, m1), (ef3, m3) = res[1], res[3]
    if min(len(ef1), len(ef3), len(m1[1]), len(m3[1]), len(m3[2])) < 2:
        return False, "fewer than two SCF steps in the log"
    d_ef = max(abs(a - b) for a, b in zip(ef1, ef3))
    d_m = max(abs(abs(a[2]) - abs(b[2])) for a, b in zip(m1[1], m3[1]))
    opposite = max(abs(a[2] + b[2]) for a, b in zip(m3[1], m3[2]))
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
    for name, check in (("shift_nsp3", shift_nsp3), ("shift_nsp4_hoh", shift_nsp4_hoh), ("af_fe", af_fe),
                        ("af_fe_local", lambda b, s: af_fe(b, s, "local"))):
        ok, info = check(args.binary.resolve(), args.scratch_root.resolve())
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
