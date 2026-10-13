"""Cr: Omega_RPA, Omega_force = N U (I + chi0(q,0) U) and the LSWT of oracles/cr_lswt.dat from static run_path.py runs.

Usage: report_force.py <oracles/cr_lswt.dat> <dir[,dir...]>
The directories are the --part slices of one static run (new _adiabatic.dat with the force columns). Prints, for the q with
largest cubic-equivalent component >= 1/8 and for the rest: min / median / max of Omega_RPA / LSWT, Omega_force / LSWT and
Omega_RPA / Omega_force (positive branches), and the share of Omega_force / LSWT inside 0.75 to 1.3. Then the Gamma-X ladder
(Omega_RPA / Omega_force next to 1/sqrt(1 - lambda_u)), with lambda_u = 2 + U tr chi0(0, 0): the sum of the two eigenvalues of
I + chi0 U at q = 0 is 2 + U tr chi0 for U_1 = U_2, and one of them is the Goldstone zero.
"""
import math
import pathlib
import statistics
import sys

canon = lambda q: tuple(sorted((round(abs(float(x)), 6) for x in q), reverse=True))
lswt = {canon(l.split()[:3]): float(l.split()[5]) / 13605.693 for l in pathlib.Path(sys.argv[1]).read_text().splitlines() if not l.startswith("#")}
rows = {}   # canonical q -> (Omega_RPA+, Omega_force+)
gamma = None
for d in map(pathlib.Path, sys.argv[2].split(",")):
    qs = [tuple(float(x) for x in l.split()) for l in (d / "qpath.dat").read_text().splitlines()]
    for k, (q, l) in enumerate(zip(qs, (d / "spin_response_adiabatic.dat").read_text().splitlines()[1:]), 1):
        a = [float(x) for x in l.split()]
        n = (len(a) - 8) // 2
        rows[q] = (max(a[7:7 + n]), max(a[8 + n:]))
        if not any(q):
            u = float(next(s.split()[1] for s in (d / "spin_response_state.dat").read_text().splitlines()[-2:]))
            tr0 = float((d / f"spin_response_q{k:03d}.dat").read_text().splitlines()[1].split()[1])
            gamma = 2.0 + u * tr0
stats = lambda v: f"{min(v):.3f} {statistics.median(v):.3f} {max(v):.3f}"
for name, sel in (("max component >= 1/8", lambda q: canon(q)[0] >= 0.125), ("max component < 1/8", lambda q: 0 < canon(q)[0] < 0.125)):
    s = [q for q in rows if sel(q) and lswt[canon(q)] > 0]
    f = [rows[q][1] / lswt[canon(q)] for q in s]
    print(f"{name}: {len(s)} q | RPA/LSWT min/median/max {stats([rows[q][0] / lswt[canon(q)] for q in s])} | force/LSWT {stats(f)}, "
          f"{sum(0.75 <= x <= 1.3 for x in f)} of {len(f)} in [0.75, 1.3] | RPA/force {stats([rows[q][0] / rows[q][1] for q in s])}")
print(f"lambda_u = {gamma:.5f}; 1/sqrt(1 - lambda_u) = {1 / math.sqrt(1 - gamma):.4f}")
print("# Gamma-X ladder: q | RPA (Ry) | force (Ry) | RPA/force | LSWT (Ry) | force/LSWT")
for q in sorted(q for q in rows if q[1] == 0 and q[2] == 0 and 0 < q[0] <= 0.25):
    r, f = rows[q]
    print(f"{q[0]:.5f} | {r:.6e} | {f:.6e} | {r / f:.4f} | {lswt[canon(q)]:.6e} | {f / lswt[canon(q)]:.3f}")
