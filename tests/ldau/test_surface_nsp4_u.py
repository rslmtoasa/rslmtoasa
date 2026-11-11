#!/usr/bin/env python3
"""LdaUSurfaceNsp4: hubbard_radial_moments = 'local' at nsp = 4 on the reduced Ni(001) surface (calctype S, deck tests/ldau/decks/ni001_soc_u3, U = 3,
J = 0.7, Chebyshev recursion, 4 steps, single thread, moments along +z, density_policy = 'constrained_spiral'). SOC is held at zero (soc_scale = 0,
fix_soc), so the nsp = 4 spinor correction (V by type) must reproduce the nsp = 3 one: per-iteration Diff, Fermi energy, spin moment and spin-diagonal
d density-matrix entries agree to SURF_SOC0_TOL. The decks are those of test_surface_nsp3_u.py (surface() there takes nsp = 3 or 4).
  nsp4_vs_nsp3      first order
  nsp4_vs_nsp3_hoh  hoh (V_U in enim)
  engaged           nsp = 4 local against spectral: spin moments differ by more than 1e-3 muB at some step
Bug caught: the correction missing or in the wrong spin frame on a surface cluster at nsp = 4 (a skipped correction passes the parity cases only if
'engaged' also fails). Not covered: SOC on, noncollinear axes, an independent reference for the corrected moments (nsp = 3 is the same code).

usage: test_surface_nsp4_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402
import test_surface_hoh_u as sh  # noqa: E402
import test_surface_nsp3_u as s3  # noqa: E402


def compare(a, b):
    for x in (a, b):
        if isinstance(x, str):
            return False, x
    return rr.compare_signatures(rr.signature(b, True), rr.signature(a, True), rr.SURF_SOC0_TOL)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    n3, n4 = s3.surface(binary, base, "nsp3", 3, "local", False), s3.surface(binary, base, "nsp4", 4, "local", False)
    h3, h4 = s3.surface(binary, base, "nsp3_hoh", 3, "local", True), s3.surface(binary, base, "nsp4_hoh", 4, "local", True)
    spec = s3.surface(binary, base, "nsp4_spectral", 4, "spectral", False)
    cases = [("nsp4_vs_nsp3", compare(n3, n4)), ("nsp4_vs_nsp3_hoh", compare(h3, h4))]
    failed_run = next((x for x in (n4, spec) if isinstance(x, str)), None)
    cases.append(("engaged", (False, failed_run) if failed_run else sh.engaged(n4, spec)))
    failed = False
    for name, (ok, info) in cases:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
