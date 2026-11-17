#!/usr/bin/env python3
"""LdaUHohEnergy: the Hubbard energy of hubbard_radial_moments = 'local' with hoh = .true. (U = 2 eV on Fe d, single thread; V_U in enim). The routing adds the
energy of P_out to physical_total_energy as without hoh; for collinear setups the spin representations describe the same state, so E_Hub_P_OUT (Ry, 16
digits) and physical_total_energy (Ry) must agree at every step. This is a consistency check between nsp = 3 and nsp = 1, not a check of the energy value
(the finite-difference dE/dU of bulk bcc Fe is in docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md):
  ks_af_nsp3_nsp1   k-space two-atom antiferromagnet, 3 steps (decks of test_kspace_energy_u.py)
  rs_af_nsp3_nsp1   the same antiferromagnet in real space (Chebyshev), 2 steps (decks of test_local_energy_u.py)
  imp_nsp3_nsp1     B2FeCo impurity (calctype I, decks of test_impurity_nsp3_u.py, J = 0), 4 steps
  surf_nsp3_nsp1    Ni(001) surface (calctype S, decks of test_surface_nsp3_u.py, U = 3 eV, J = 0.7), 4 steps
  engaged           |E_Hub_P_OUT| of the four nsp = 1 runs is above 1e-2 Ry at every step (a skipped routing prints 0 and makes the parities hold trivially)
Bug caught: the energy not routed with hoh, or built from the wrong spin frame or site there. Not covered: nsp = 4, SOC, an independent reference for the
energy itself, ccor_2c.

usage: test_hoh_energy_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_cluster_energy_u as ce  # noqa: E402
import test_impurity_nsp3_u as imp  # noqa: E402
import test_kspace_energy_u as ke  # noqa: E402
import test_local_energy_u as le  # noqa: E402
import test_local_nsp4_workflows_u as wf  # noqa: E402
import test_nsp34_u as t3  # noqa: E402
import test_surface_nsp3_u as surf  # noqa: E402


def af(binary: Path, base: Path, name: str, nsp: int, kspace: bool):
    steps = wf.NSTEP if kspace else 2
    work = t3.af_deck(base, name, nsp=nsp, mom2=(0.0, 0.0, -1.0), swap2=(nsp == 1), nstep=steps, kspace=kspace, radial="local", hoh=True)
    status, log = gh.run_binary(binary, work)
    return f"{name}: run failed\n{log[-1200:]}" if status else ke.energies(log, steps)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    r = {"k1": af(binary, base, "ks_af_nsp1", 1, True), "k3": af(binary, base, "ks_af_nsp3", 3, True),
         "r1": af(binary, base, "rs_af_nsp1", 1, False), "r3": af(binary, base, "rs_af_nsp3", 3, False),
         "i1": ce.energies(imp.host_and_impurity(binary, base / "imp", "imp_nsp1", 1, "local", True)),
         "i3": ce.energies(imp.host_and_impurity(binary, base / "imp", "imp_nsp3", 3, "local", True)),
         "s1": ce.energies(surf.surface(binary, base / "surf", "surf_nsp1", 1, "local", True)),
         "s3": ce.energies(surf.surface(binary, base / "surf", "surf_nsp3", 3, "local", True))}
    cases = [("ks_af_nsp3_nsp1", le.compare(r["k3"], r["k1"], "k-space AF nsp=3 against nsp=1")),
             ("rs_af_nsp3_nsp1", le.compare(r["r3"], r["r1"], "real-space AF nsp=3 against nsp=1")),
             ("imp_nsp3_nsp1", le.compare(r["i3"], r["i1"], "impurity nsp=3 against nsp=1")),
             ("surf_nsp3_nsp1", le.compare(r["s3"], r["s1"], "surface nsp=3 against nsp=1"))]
    bad = next((v for v in (r["k1"], r["r1"], r["i1"], r["s1"]) if isinstance(v, str)), None)
    if bad:
        cases.append(("engaged", (False, bad)))
    else:
        low = {k: min(abs(x) for x in r[k][0]) for k in ("k1", "r1", "i1", "s1")}
        cases.append(("engaged", (min(low.values()) > le.ENGAGED, "smallest |E_Hub_P_OUT| of the nsp=1 runs (Ry, need > " + f"{le.ENGAGED:.0e}): " +
                                  ", ".join(f"{k} {v:.2e}" for k, v in low.items()))))
    failed = False
    for name, (ok, info) in cases:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
