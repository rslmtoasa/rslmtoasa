#!/usr/bin/env python3
"""WP9 Battery C -- multi-sublattice Goldstone measurement (bcc FeCo).

WHAT IS BEING TESTED
=====================
`tests/KNOWN_ISSUES.md` records a historical finding (pre-WP1 architecture):
the multi-sublattice frozen-magnon branch (`&frozen_magnon branch_mode =
'auto'`, `source/calculation.f90::post_processing_frozen_magnon_auto`) failed
the Goldstone theorem on a two-sublattice bcc FeCo system -- the acoustic
branch had a finite omega(Gamma) (~0.28 Ry, k-space route) instead of going to
zero, as it must for a uniform rotation of every sublattice moment.

This script is a MEASUREMENT, not a gate. Its job is to report omega(Gamma)
on the current architecture exactly as measured -- it must NOT be "fixed" by
shifting, renormalizing, or tuning anything toward zero (explicit instruction
from the WP9 task). A large, nonzero result here is an ACCEPTABLE and
EXPECTED possible outcome; a passing assertion on this script would be
inappropriate. The only assertion is a loose sanity bound to catch a totally
broken run (NaN, a wildly unphysical energy scale), not Goldstone-satisfaction.

MEASURED RESULT (2026-08-07, this exact deck, real-space route, `recur=
'block'`, `lld=21`, `strux_backend='strux_lib'`, `nsp=3`, `theta_probe=20deg`):
acoustic branch omega(Gamma) = -2.15e-22 Ry -- zero to numerical noise, NOT
the ~0.28 Ry violation the KNOWN_ISSUES entry records. See
`tests/KNOWN_ISSUES.md`'s "multi-sublattice acoustic magnon not gapless at
Gamma" entry for the full comparison and caveats (in particular: this is the
real-space route; the original violation was measured on the k-space route,
which this script does not exercise or re-validate).

`base/Fe.nml` and `base/Co.nml` are seeded from an already-converged
reference potential (not the flat `mom=(0,0,1)` starting guess) -- starting
cold made the Broyden mixing oscillate for several steps before settling
(see git history / session notes if this ever needs re-deriving); seeding
from a converged state keeps this test's runtime and stability reasonable.

USAGE
=====
    run_wp9_goldstone.py --binary <rslmto.x> --deck-root <dir>
                          --scratch-root <dir>
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Loose sanity bound: catches a broken run (NaN, or an omega at the scale of
# a full band energy, ~Ry), not a physics judgement. The historical violation
# itself was ~0.28 Ry, so this must stay well above that.
SANITY_BOUND_RY = 5.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--deck-root", type=Path, required=True)
    ap.add_argument("--scratch-root", type=Path, required=True)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()

    base = args.deck_root / "base"
    workdir = args.scratch_root / "run"
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)

    result = subprocess.run(
        [str(args.binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=args.timeout,
    )
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(f"ERROR: rslmto.x exited {result.returncode}")

    branches = workdir / "frozen_magnon_branches.dat"
    if not branches.exists():
        raise SystemExit(f"ERROR: {branches} was not produced")
    text = branches.read_text()

    rows = [line.split() for line in text.splitlines() if line and not line.startswith("#")]
    if not rows:
        raise SystemExit("ERROR: no data rows in frozen_magnon_branches.dat")

    # Format: q1 q2 q3 branch omega nactive
    at_gamma = [r for r in rows if abs(float(r[0])) < 1e-9 and abs(float(r[1])) < 1e-9 and abs(float(r[2])) < 1e-9]
    if not at_gamma:
        raise SystemExit("ERROR: no Gamma (q=0) row found -- deck must include q=0 as the first q-point")

    branch_omegas = {int(r[3]): float(r[4]) for r in at_gamma}
    acoustic_branch = min(branch_omegas, key=lambda b: abs(branch_omegas[b]))
    omega_acoustic = branch_omegas[acoustic_branch]

    print("=== WP9 Battery C: multi-sublattice Goldstone measurement (bcc FeCo, real-space route) ===")
    print(f"  omega(Gamma) per branch: {branch_omegas}")
    print(f"  smallest-|omega| branch at Gamma (the acoustic candidate): "
          f"branch {acoustic_branch}, omega = {omega_acoustic:.6e} Ry")
    print("  Historical (pre-WP1, k-space route) violation: ~0.28 Ry.")
    print("  This measurement is reported, not gated, per the WP9 task instructions --")
    print("  see tests/KNOWN_ISSUES.md for interpretation and caveats.")

    if math.isnan(omega_acoustic) or math.isinf(omega_acoustic):
        print("FAIL: omega(Gamma) is NaN/Inf -- run is broken, not just Goldstone-violating.")
        return 1
    if abs(omega_acoustic) > SANITY_BOUND_RY:
        print(f"FAIL: |omega(Gamma)| = {abs(omega_acoustic):.3e} Ry exceeds the loose sanity "
              f"bound {SANITY_BOUND_RY} Ry -- this is not a Goldstone-quality judgement, "
              "it means the run produced an unphysical energy scale.")
        return 1

    print("PASS (sanity bound only -- see printed omega(Gamma) for the actual physics result)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
