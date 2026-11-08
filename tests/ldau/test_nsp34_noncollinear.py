#!/usr/bin/env python3
"""Fixed U at nsp = 3 and 4 with noncollinear moment directions: symmetry and consistency checks.

There is no external oracle for these states. The expected values are exact symmetries of the problem or identities
between two routes; each check compares runs of the same code, so it is a consistency check, not an independent reference.
All runs are real-space Chebyshev (block recursion is not rotation invariant, see tests/KNOWN_ISSUES.md), U = 2 eV on Fe d, 2 steps
(the spinor density matrix has spin coherences from step 2).

Checks (each prints PASS/FAIL by name; the test fails if any fails):
  rot_fm          FM Fe, nsp = 3, moment along z and along (60 deg, azimuth 45 deg): E_F, |m| and the eigenvalues of the 10x10
                  spinor density matrix are equal (a global spin rotation is a symmetry without SOC); the tilted run has large
                  spin-off-diagonal blocks in the global frame, so the check is not trivial.
  af_angle        two-atom Fe, nsp = 3, site moments at 90 deg: rotating both moments together, and exchanging the two sites,
                  leave E_F and the site |m| unchanged; moments parallel (nsp = 3) equal the nsp = 1 ferromagnet.
  soc_cubic       FM Fe with SOC at nsp = 4, moment along z and along x (related by a cubic rotation): equal; with SOC held at
                  zero nsp = 4 equals nsp = 3.
  ks_af_mirror    two-atom Fe AF in k-space at nsp = 3 (the same potential file on both sites, mom = +z and -z): |m1| = |m2| at
                  every step (the sublattice mirror symmetry; the nsp = 1 spin-swapped route breaks it, see KNOWN_ISSUES).
  restart         tilted AF: a restart from <label>_out.nml and <label>_dm.nml after step 2 reproduces step 3 of the continuous run.

Bugs caught: a spinor V or density matrix built in the wrong spin frame (rot_fm, af_angle), a tilted-axis error that only the
site-to-site coupling shows (af_angle), V_U lost in the SOC path (soc_cubic), the k-space route dropping the mirror symmetry
(ks_af_mirror), moment axes or density matrix not persisted in the restart files (restart).

usage: test_nsp34_noncollinear.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_hoh_u as hu  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL = 2.0e-6  # Ry and muB: the logs print six decimals
EIG_TOL = 1.0e-8
SQ2 = 1 / math.sqrt(2)


def clean(log: str) -> str:
    return re.sub(r"\x1b\[[0-9;]*m", "", log)


def moments(log: str, atom: int = 1) -> list[float]:
    return [float(x) for x in re.findall(rf"Spin moment of atom\s+{atom} is\s+(\S+)", log)]


def run(binary: Path, work: Path) -> tuple[int, str]:
    status, log = gh.run_binary(binary, work)
    return status, clean(log)


def spinor_eigenvalues(dm_file: Path) -> tuple[np.ndarray, float]:
    """Eigenvalues of the 10x10 d spinor density matrix and its largest spin-off-diagonal element."""
    text = dm_file.read_text()
    p = np.zeros((10, 10), complex)
    for s in (1, 2):
        for u in (1, 2):
            for m in range(1, 6):
                re_ = re.search(rf"dm_re\(3,{s},{u},{m},:\)\s*=\s*(.*)", text)
                im_ = re.search(rf"dm_im\(3,{s},{u},{m},:\)\s*=\s*(.*)", text)
                if re_ and im_:
                    p[(s - 1) * 5 + m - 1, (u - 1) * 5:u * 5] = [float(a) + 1j * float(b) for a, b in zip(re_.group(1).split(), im_.group(1).split())]
    return np.sort(np.linalg.eigvalsh((p + p.conj().T) / 2)), float(np.abs(p[:5, 5:]).max())


def fm_deck(scratch: Path, name: str, mom: tuple[float, float, float], *, nsp: int = 3, soc_off: bool = False,
            fixture: Path = hu.BULK_FIXTURE) -> Path:
    work = hu.bulk_deck(scratch, name, hoh=False, u_ev=2.0, nsp=nsp, fixture=fixture)
    gh.set_line(work / "input.nml", r"^nstep\s*=.*$", "nstep = 2")
    if soc_off:
        path = work / "input.nml"
        path.write_text(path.read_text().replace("&self", "&self\nsoc_scale = 0.0\nfix_soc = .true.", 1))
    fe = work / "Fe.nml"
    fe.write_text(t3.moment_line(fe.read_text(), mom))
    return work


def rot_fm(binary: Path, scratch: Path) -> tuple[bool, str]:
    res = {}
    for name, mom in (("z", (0.0, 0.0, 1.0)), ("tilt", (0.6123724357, 0.6123724357, 0.5))):
        work = fm_deck(scratch, f"rot_{name}", mom)
        status, log = run(binary, work)
        if status != 0:
            return False, f"run failed ({name})\n{log[-1500:]}"
        eig, offdiag = spinor_eigenvalues(work / "Fe_dm.nml")
        res[name] = (t3.af_numbers(log)[0], moments(log), eig, offdiag)
    (ef0, m0, e0, o0), (ef1, m1, e1, o1) = res["z"], res["tilt"]
    d_ef = max(abs(a - b) for a, b in zip(ef0, ef1))
    d_m = max(abs(a - b) for a, b in zip(m0, m1))
    d_eig = float(np.abs(e0 - e1).max())
    ok = len(ef0) == len(ef1) == 2 and d_ef < TOL and d_m < TOL and d_eig < EIG_TOL and o1 > 0.1 and o0 < 1e-6
    return ok, (f"z against tilted: max |dE_F| {d_ef:.1e} Ry, max ||m| difference| {d_m:.1e} muB, max spinor-P eigenvalue difference "
                f"{d_eig:.1e}; max |P_up,dn| in the global frame {o0:.3f} (z) and {o1:.3f} (tilted, need > 0.1)")


def af_angle(binary: Path, scratch: Path) -> tuple[bool, str]:
    setups = {
        "A": dict(nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0)),    # 90 deg
        "B": dict(nsp=3, mom1=(1.0, 0.0, 0.0), mom2=(0.0, 0.0, -1.0)),   # A rotated about y by 90 deg
        "C": dict(nsp=3, mom1=(1.0, 0.0, 0.0), mom2=(0.0, 0.0, 1.0)),    # A with the two sites exchanged
        "D": dict(nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, 1.0)),    # parallel
        "E": dict(nsp=1, swap2=False),                                    # nsp = 1 ferromagnet
    }
    out = {}
    for key, kw in setups.items():
        status, log = run(binary, t3.af_deck(scratch, f"ang_{key}", **kw))
        if status != 0:
            return False, f"run failed ({key})\n{log[-1500:]}"
        out[key] = (t3.af_numbers(log)[0], moments(log, 1), moments(log, 2))
    def diff(x: str, y: str) -> float:
        (ef_x, a_x, b_x), (ef_y, a_y, b_y) = out[x], out[y]
        pairs = list(zip(ef_x, ef_y)) + list(zip(a_x, a_y)) + list(zip(b_x, b_y))
        return max(abs(p - q) for p, q in pairs)
    d_rot, d_swap, d_fm = diff("A", "B"), diff("A", "C"), diff("D", "E")
    ok = max(d_rot, d_swap, d_fm) < TOL
    return ok, (f"E_F and site |m| over 2 steps: rotated {d_rot:.1e}, sites exchanged {d_swap:.1e}, parallel nsp=3 against nsp=1 FM "
                f"{d_fm:.1e} (need < {TOL})")


def soc_cubic(binary: Path, scratch: Path) -> tuple[bool, str]:
    runs = {"z": dict(mom=(0.0, 0.0, 1.0), nsp=4), "x": dict(mom=(1.0, 0.0, 0.0), nsp=4),
            "x_soc0": dict(mom=(1.0, 0.0, 0.0), nsp=4, soc_off=True), "x_nsp3": dict(mom=(1.0, 0.0, 0.0), nsp=3)}
    out = {}
    for key, kw in runs.items():
        status, log = run(binary, fm_deck(scratch, f"soc_{key}", fixture=hu.SOC_FIXTURE, **kw))
        if status != 0:
            return False, f"run failed ({key})\n{log[-1500:]}"
        out[key] = (t3.af_numbers(log)[0], moments(log))
    def diff(x: str, y: str) -> float:
        return max(abs(p - q) for p, q in list(zip(out[x][0], out[y][0])) + list(zip(out[x][1], out[y][1])))
    d_cubic, d_soc0 = diff("z", "x"), diff("x_soc0", "x_nsp3")
    ok = max(d_cubic, d_soc0) < 5 * TOL
    return ok, (f"nsp=4 with SOC, moment along z against x: {d_cubic:.1e}; nsp=4 with SOC off against nsp=3: {d_soc0:.1e} "
                f"(E_F in Ry and |m| in muB over 2 steps, need < {5 * TOL:.0e})")


def ks_af_mirror(binary: Path, scratch: Path) -> tuple[bool, str]:
    status, log = run(binary, t3.af_deck(scratch, "ks_af", nsp=3, kspace=True, nstep=3))
    if status != 0:
        return False, f"run failed\n{log[-1500:]}"
    m1, m2 = moments(log, 1), moments(log, 2)
    if len(m1) < 3 or len(m2) < 3:
        return False, "fewer than three SCF steps in the log"
    worst = max(abs(a - b) for a, b in zip(m1, m2))
    return worst < TOL, f"|m1| - |m2| over 3 k-space steps: max {worst:.1e} muB (need < {TOL}); |m1| = {', '.join(f'{x:.6f}' for x in m1)}"


def restart(binary: Path, scratch: Path) -> tuple[bool, str]:
    tilted = dict(nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0))
    status, log3 = run(binary, t3.af_deck(scratch, "rs_cont", nstep=3, **tilted))
    if status != 0:
        return False, f"continuous run failed\n{log3[-1500:]}"
    first = t3.af_deck(scratch, "rs_first", nstep=2, **tilted)
    status, log = run(binary, first)
    if status != 0:
        return False, f"first run failed\n{log[-1500:]}"
    second = scratch / "rs_second"
    second.mkdir()
    for f in ("lattice.nml", "Fe1_out.nml", "Fe2_out.nml", "Fe1_dm.nml", "Fe2_dm.nml"):
        shutil.copy(first / f, second / f)
    (second / "input.nml").write_text(re.sub(r"(?m)^nstep\s*=.*$", "nstep = 1", (first / "input.nml").read_text()))
    status, log = run(binary, second)
    if status != 0:
        return False, f"restart run failed\n{log[-1500:]}"
    ef_c, ef_r = t3.af_numbers(log3)[0][2], t3.af_numbers(log)[0][0]
    m_c, m_r = moments(log3)[2], moments(log)[0]
    ok = abs(ef_c - ef_r) < TOL and abs(m_c - m_r) < TOL
    return ok, (f"step 3 of the continuous run against step 1 of the restart: E_F {ef_c:+.6f} / {ef_r:+.6f} Ry, "
                f"|m| {m_c:.6f} / {m_r:.6f} muB (need < {TOL})")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    failed = False
    for name, check in (("rot_fm", rot_fm), ("af_angle", af_angle), ("soc_cubic", soc_cubic),
                        ("ks_af_mirror", ks_af_mirror), ("restart", restart)):
        ok, info = check(args.binary.resolve(), args.scratch_root.resolve())
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
