#!/usr/bin/env python3
"""LdaULocalEnergy: the Hubbard energy of hubbard_radial_moments = 'local' at nsp = 3 and 4 (real-space bulk, Chebyshev, U = 2 eV on Fe d, J = 0, 2 steps,
single thread). Without the routing the printed E_Hub_P_OUT is 0 and physical_total_energy is base_physical_total_energy; with it the energy of
the spinor P_out is added. For collinear setups the spin representations describe the same state, so the printed E_Hub_P_OUT (Ry, 16 digits) and
physical_total_energy (Ry) must agree at every step:
  af_nsp3_nsp1  antiferromagnet, nsp = 3 (same file, mom = +z, -z) against nsp = 1 (scalar energy; the second site's file with swapped channels)
  soc0_nsp4_nsp3  (z, x) pair, nsp = 4 with SOC held at zero against nsp = 3
  fm_nsp4_nsp2  ferromagnet along z with SOC, nsp = 4 against nsp = 2 (the spinor energy of the bulk)
  engaged       the nsp = 3 antiferromagnet has E_Hub_P_OUT above 1e-2 Ry (a skipped routing prints 0 and makes the parities hold trivially)
Bug caught: the energy not routed at nsp >= 3, or built from the wrong spin frame or site. Not covered: k-space, impurities, J != 0, an independent
reference for the energy itself (the nsp = 1 energy is the reference), a finite-difference dE/dU check.

usage: test_local_energy_u.py --binary rslmto.x --scratch-root DIR
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
import test_nsp34_noncollinear as nc  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL_HUB = 1.0e-8  # Ry
TOL_TOT = 1.0e-6  # Ry, on totals of about -5084
ENGAGED = 1.0e-2  # Ry
SOC_OFF = "soc_scale = 0.0\nfix_soc = .true.\n"


def energies(status: int, log: str):
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if status != 0:
        return f"run failed\n{log[-1200:]}"
    hub = [float(x) for x in re.findall(r"E_Hub_P_OUT=\s*(\S+)", log)]
    tot = [float(x) for x in re.findall(r"(?<!base_)physical_total_energy=\s*(\S+)", log)]
    return (hub, tot) if len(hub) == len(tot) == 2 else f"expected 2 energies, got {len(hub)}, {len(tot)}"


def deck_run(binary: Path, work: Path):
    return energies(*gh.run_binary(binary, work))


def compare(a, b, text: str):
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    dh = max(abs(x - y) for x, y in zip(a[0], b[0]))
    dt = max(abs(x - y) for x, y in zip(a[1], b[1]))
    return dh < TOL_HUB and dt < TOL_TOT, f"{text}: max |dE_Hub| {dh:.1e} Ry (need < {TOL_HUB:.0e}), max |dE_total| {dt:.1e} Ry (need < {TOL_TOT:.0e})"


def fm(binary: Path, scratch: Path, name: str, nsp: int):
    work = nc.fm_deck(scratch, name, (0.0, 0.0, 1.0), nsp=nsp, fixture=hu.SOC_FIXTURE)
    path = work / "input.nml"
    path.write_text(path.read_text().replace("&hamiltonian", "&hamiltonian\nhubbard_radial_moments = 'local'\n", 1))
    return energies(*nc.run(binary, work))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    af = dict(radial="local", mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, -1.0))
    pair = dict(radial="local", mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0))
    a1 = deck_run(binary, t3.af_deck(base, "af_nsp1", nsp=1, **{**af, "mom2": (0.0, 0.0, 1.0)}, swap2=True))
    a3 = deck_run(binary, t3.af_deck(base, "af_nsp3", nsp=3, **af))
    p3 = deck_run(binary, t3.af_deck(base, "pair_nsp3", nsp=3, **pair))
    p4 = deck_run(binary, t3.af_deck(base, "pair_nsp4_soc0", nsp=4, self_extra=SOC_OFF, **pair))
    f2, f4 = fm(binary, base, "fm_nsp2", 2), fm(binary, base, "fm_nsp4", 4)
    results = [("af_nsp3_nsp1", compare(a3, a1, "nsp = 3 against nsp = 1, AF")),
               ("soc0_nsp4_nsp3", compare(p4, p3, "nsp = 4 SOC off against nsp = 3, (z, x) pair")),
               ("fm_nsp4_nsp2", compare(f4, f2, "nsp = 4 against nsp = 2, FM with SOC")),
               ("engaged", (False, a3) if isinstance(a3, str) else (min(a3[0]) > ENGAGED, f"E_Hub_P_OUT at nsp = 3: {a3[0]} Ry (need > {ENGAGED})"))]
    failed = False
    for name, (ok, info) in results:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
