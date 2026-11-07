#!/usr/bin/env python3
"""Fixed U with hoh: V_U is an additive on-site term (enim), checked against the no-hoh route (V_U in ee).

Checks (each prints PASS/FAIL by name; the test fails if any fails):
  fermi_shift       bcc Fe, one real-space Chebyshev step, U = 2 eV: E_F(U=2) - E_F(U=0) with hoh must exceed 0.05 Ry and
                    match the same difference without hoh within 0.01 Ry. The reference is the no-hoh route, where V_U is
                    added to ee: a different place in the code from enim. Measured on the deck: +0.0679 (hoh) and +0.0668 (no hoh).
  fermi_shift_ccor  the same with ccor_2c on (deck with vrmax, which the correction needs), both shifts above 0.05 Ry and
                    within 0.01 Ry of each other. This is the only check of ccor_2c with U without hoh (V_U in the folded
                    ee + eecc blocks). Measured: +0.0688 (hoh) and +0.0709 (no hoh).
  fermi_shift_soc   the same at nsp = 2 with spin-orbit coupling (spinor U): both shifts above 0.05 Ry and within 0.01 Ry of each
                    other. Only the spin-diagonal placement of V_U is tested: |P_up,dn| is 1.3e-2 on this deck. Measured: +0.0680 (hoh), +0.0668 (no hoh).
  guard_kfirst      hoh with U in k-space with kspace_ham_order = 'first' stops with "kspace_ham_order='second'"

Bugs caught: V_U not reaching the hoh Hamiltonian (the hoh shift is then 0, first check fails); V_U added inside h
(the hoh shift moves away from the no-hoh one); V_U lost in the ccor_2c operator blocks (the shift is 0 or differs
between hoh and no hoh); V_U dropped at nsp = 2 with hoh; a guard that stops firing (hoh with U would run on an Hamiltonian without V_U).

The code rewrites fermi in input.nml after a run, so every deck is built fresh from the committed fixture.

usage: test_hoh_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402

BULK_FIXTURE = gh.ROOT / "tests/scf/cases/bulk/bccFe_hubbard_nsp1"
CCOR_FIXTURE = gh.ROOT / "tests/scf/cases/bulk/bccFe_hubbard_ccor"
SOC_FIXTURE = gh.ROOT / "tests/scf/cases/bulk/bccFe_hubbard"  # nsp = 2 with SOC
SHIFT_MIN = 0.05  # Ry
SHIFT_TOL = 0.01  # Ry


def fermi_of(log: str) -> float:
    found = re.findall(r"Free Fermi energy:\s*(-?\d+\.\d+)", log)
    if not found:
        raise RuntimeError("no 'Free Fermi energy' line in the log")
    return float(found[-1])


def bulk_deck(scratch: Path, name: str, *, hoh: bool, u_ev: float, nsp: int = 1, fixture: Path = BULK_FIXTURE,
              window: float = 3.0, lld: int = 200) -> Path:
    work = gh.fresh(fixture, scratch, name)
    path = work / "input.nml"
    gh.set_line(path, r"^hoh\s*=.*$", f"hoh = {'.true.' if hoh else '.false.'}")
    gh.set_line(path, r"^hubbard_u_general\(1,3\)\s*=.*$", f"hubbard_u_general(1,3) = {u_ev}")
    gh.set_line(path, r"^recur\s*=.*$", "recur = 'chebyshev'")
    gh.set_line(path, r"^rc\s*=.*$", "rc = 040")
    gh.set_line(path, r"^lld\s*=.*$", f"lld = {lld}")
    gh.set_line(path, r"^energy_min\s*=.*$", f"energy_min = -{window}")
    gh.set_line(path, r"^energy_max\s*=.*$", f"energy_max = {window}")
    gh.set_line(path, r"^nsp\s*=.*$", f"nsp = {nsp}")
    return work


def shifts(binary: Path, scratch: Path, tag: str, **deck) -> tuple[float, float] | str:
    """E_F(U=2) - E_F(U=0) in Ry with hoh and without, or the failure text."""
    out = {}
    for hoh in (True, False):
        ef = {}
        for u_ev in (0.0, 2.0):
            work = bulk_deck(scratch, f"{tag}_hoh{int(hoh)}_u{int(u_ev)}", hoh=hoh, u_ev=u_ev, **deck)
            status, log = gh.run_binary(binary, work)
            if status != 0:
                return f"run failed (hoh={hoh}, U={u_ev})\n{log[-1500:]}"
            ef[u_ev] = fermi_of(log)
        out[hoh] = ef[2.0] - ef[0.0]
    return out[True], out[False]


def check_shifts(result: tuple[float, float] | str) -> tuple[bool, str]:
    if isinstance(result, str):
        return False, result
    hoh, plain = result
    ok = min(hoh, plain) > SHIFT_MIN and abs(hoh - plain) < SHIFT_TOL
    return ok, f"E_F shift, U=2 eV: hoh {hoh:+.4f} Ry, no hoh {plain:+.4f} Ry " \
               f"(need both > {SHIFT_MIN} and |difference| < {SHIFT_TOL})"


def fermi_shift(binary: Path, scratch: Path) -> tuple[bool, str]:
    return check_shifts(shifts(binary, scratch, "shift"))


def fermi_shift_ccor(binary: Path, scratch: Path) -> tuple[bool, str]:
    # vrmax widens the spectrum: a +-3 Ry window is not enough with ccor_2c (Chebyshev moments diverge)
    return check_shifts(shifts(binary, scratch, "shift_ccor", fixture=CCOR_FIXTURE, window=5.0, lld=300))


def expect_fatal(status: int, log: str, message: str) -> tuple[bool, str]:
    ok = status != 0 and message in log.lower()
    return ok, "" if ok else f"status={status}, message {message!r} not found\n{log[-1500:]}"


def fermi_shift_soc(binary: Path, scratch: Path) -> tuple[bool, str]:
    return check_shifts(shifts(binary, scratch, "shift_soc", fixture=SOC_FIXTURE, nsp=2))


def guard_kfirst(binary: Path, scratch: Path) -> tuple[bool, str]:
    work = gh.fresh(gh.HOST_FIXTURE, scratch, "guard_kfirst")
    gh.configure_host(work)
    path = work / "input.nml"
    text = path.read_text().replace("hubbard_radial_moments = 'local'\n", "hoh = .true.\n")
    text = re.sub(r"(?m)^&reciprocal\s*$", "&reciprocal\nkspace_ham_order = 'first'", text, count=1)
    path.write_text(text)
    return expect_fatal(*gh.run_binary(binary, work), "kspace_ham_order='second'")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    failed = False
    for name, check in (("fermi_shift", fermi_shift), ("fermi_shift_ccor", fermi_shift_ccor),
                        ("fermi_shift_soc", fermi_shift_soc), ("guard_kfirst", guard_kfirst)):
        ok, info = check(args.binary.resolve(), args.scratch_root.resolve())
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
