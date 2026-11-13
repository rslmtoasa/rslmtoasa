#!/usr/bin/env python3
"""LdaUClusterEnergy: the Hubbard energy of hubbard_radial_moments = 'local' on impurity and surface clusters (first order, 4 steps, single thread).
Without the routing the printed E_Hub_P_OUT is 0 and physical_total_energy is base_physical_total_energy; with it the energy of the cluster sites' P_out
(the host types are not included) is added. Both representations describe the same collinear state, so E_Hub_P_OUT (Ry, 16 digits) and
physical_total_energy (Ry) must agree at every step:
  imp_nsp3_nsp1   B2FeCo impurity (calctype I, tests/ldau/test_impurity_nsp3_u.py decks, U = 2 eV, J = 0): nsp = 3 against nsp = 1 (scalar energy)
  surf_nsp3_nsp1  Ni(001) surface (calctype S, tests/ldau/test_surface_nsp3_u.py decks, U = 3 eV, J = 0.7, constrained_spiral): nsp = 3 against nsp = 1
  engaged         |E_Hub_P_OUT| of both nsp = 1 runs is above 1e-2 Ry at every step (a skipped routing prints 0 and makes the parities hold trivially)
Bug caught: the energy not routed on a cluster, taken from the wrong type (the recursion site index is not the type), or built from the wrong spin frame.
Not covered: hoh (not routed anywhere), nsp = 4, k-space, an independent reference for the energy itself (the nsp = 1 energy is the reference), a
finite-difference dE/dU check.

usage: test_cluster_energy_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402
import test_impurity_nsp3_u as imp  # noqa: E402
import test_local_energy_u as le  # noqa: E402
import test_surface_nsp3_u as surf  # noqa: E402


def energies(d):
    if isinstance(d, str):
        return d
    log = rr.log_text(d)
    hub = [float(x) for x in re.findall(r"E_Hub_P_OUT=\s*(\S+)", log)]
    tot = [float(x) for x in re.findall(r"(?<!base_)physical_total_energy=\s*(\S+)", log)]
    return (hub, tot) if len(hub) == len(tot) == rr.NSTEP else f"expected {rr.NSTEP} energies, got {len(hub)}, {len(tot)}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    runs = {"i1": energies(imp.host_and_impurity(binary, base / "imp", "imp_nsp1", 1, "local")),
            "i3": energies(imp.host_and_impurity(binary, base / "imp", "imp_nsp3", 3, "local")),
            "s1": energies(surf.surface(binary, base / "surf", "surf_nsp1", 1, "local", False)),
            "s3": energies(surf.surface(binary, base / "surf", "surf_nsp3", 3, "local", False))}
    cases = [("imp_nsp3_nsp1", le.compare(runs["i3"], runs["i1"], "impurity nsp=3 against nsp=1")),
             ("surf_nsp3_nsp1", le.compare(runs["s3"], runs["s1"], "surface nsp=3 against nsp=1"))]
    bad = next((v for v in (runs["i1"], runs["s1"]) if isinstance(v, str)), None)
    if bad:
        cases.append(("engaged", (False, bad)))
    else:
        low = min(abs(x) for k in ("i1", "s1") for x in runs[k][0])
        cases.append(("engaged", (low > le.ENGAGED, f"smallest |E_Hub_P_OUT| of the nsp=1 runs over {rr.NSTEP} steps {low:.1e} Ry (need > {le.ENGAGED:.0e})")))
    failed = False
    for name, (ok, info) in cases:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
