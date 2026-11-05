#!/usr/bin/env python3
"""ValGbtConvergence. Relation R12: the adiabatic spectrum converges with k-sampling and smearing (Halilov et al., PRB 58, 293 (1998);
Pajda et al., PRB 64, 174402 (2001)). frozen_magnon with kspace_ham_order = 'orthogonal' on tests/spin_response/fe_bcc, zero SCF steps,
theta = 5 deg, q1 = (0.125,-0.125,0.125), q2 = (0.25,-0.25,0.25), q3 = (0.0625,-0.0625,0.0625) direct; 32^3, 64^3, 96^3 at 300 K and
100 K, 600 K at 64^3. Prints omega, DeltaOmega/sin^2(theta), mu and M per run. Numbers only, no thresholds; exit 0 unless a run fails.
usage: val_gbt_convergence.py <rslmto.x> <scratch directory>"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gbt"))
from gbt_app_common import run_orthogonal

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
Q = [(0.125, -0.125, 0.125), (0.25, -0.25, 0.25), (0.0625, -0.0625, 0.0625)]
RUNS = [(32, 300.0), (64, 300.0), (96, 300.0), (64, 100.0), (64, 600.0)]
print("theta = 5 deg; q1 = (0.125,-0.125,0.125), q2 = (0.25,-0.25,0.25), q3 = (0.0625,-0.0625,0.0625) direct")
print(f"{'mesh':>5s} {'T[K]':>6s} {'q':>6s} {'omega[meV]':>12s} {'DOm/sin2[mRy]':>14s} {'mu[Ry]':>12s} {'M[muB]':>10s}")
for nk, temp in RUNS:
    mu, moment, dom, omega = run_orthogonal(binary, scratch / f"nk{nk}_T{temp:g}", nk, temp, 5.0, Q)
    for i, (d, w) in enumerate(zip(dom, omega)):
        print(f"{nk:4d}^3 {temp:6.0f} {'q' + str(i + 1):>6s} {w:12.4f} {d:14.4f} {mu:12.6f} {moment:10.5f}")
