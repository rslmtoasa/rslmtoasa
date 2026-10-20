#!/usr/bin/env python3
"""ValGbtHarmonic. Relation: Heisenberg form of the rigid-rotation energy, DeltaOmega(q, theta) = Omega(q, theta) - Omega(0, 0)
= sin^2(theta) * c(q) + O(theta^4) (Halilov et al., PRB 58, 293 (1998); Pajda et al., PRB 64, 174402 (2001)).
frozen_magnon with kspace_ham_order = 'orthogonal' on tests/spin_response/fe_bcc, zero SCF steps, 64^3 Gamma-centred mesh, 300 K,
q1 = (0.125,-0.125,0.125) and q2 = (0.25,-0.25,0.25) direct, theta = 1, 2, 5, 10, 20 deg. Prints DeltaOmega/sin^2(theta) per theta and the
intercept and slope of a linear fit in sin^2(theta). Numbers only, no thresholds; exit 0 unless a run fails.
usage: val_gbt_harmonic.py <rslmto.x> <scratch directory>"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gbt"))
from gbt_app_common import run_orthogonal

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
Q = [(0.125, -0.125, 0.125), (0.25, -0.25, 0.25)]
THETA = [1.0, 2.0, 5.0, 10.0, 20.0]
dom = np.array([run_orthogonal(binary, scratch / f"theta{t:g}", 64, 300.0, t, Q)[2] for t in THETA])     # (theta, q)
s2 = np.sin(np.radians(THETA)) ** 2
print("DeltaOmega/sin^2(theta) [mRy], 64^3, 300 K")
print(f"{'theta[deg]':>10s} {'q1':>12s} {'q2':>12s}")
for t, row in zip(THETA, dom):
    print(f"{t:10.1f} {row[0]:12.4f} {row[1]:12.4f}")
slope, intercept = np.polyfit(s2, dom, 1)
print(f"{'intercept':>10s} {intercept[0]:12.4f} {intercept[1]:12.4f}   [mRy]")
print(f"{'slope':>10s} {slope[0]:12.4f} {slope[1]:12.4f}   [mRy per unit sin^2(theta)]")
