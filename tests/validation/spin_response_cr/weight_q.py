"""Weight of the positive-omega tr L peak of dynamic Cr runs (run_path.py --gx --dynamic) against the moments and the LSWT.

Usage: weight_q.py <oracles/cr_lswt.dat> <run directory> [<run directory> ...]
Per q and run: the omega > 0 maximum of tr L (omega_p, height h); the weight W_b = pi eta h (a Lorentzian of width eta) and W_c =
the integral over |omega - omega_p| <= 10 eta divided by (2/pi) arctan(10); omega_p over the adiabatic omega_ad of the same run;
W over the Juelich moment of the run (site 1); and the effective moment M_LSWT omega_LSWT / omega_ad that would make LSWT equal
the adiabatic frequency (omega_LSWT from the oracle at the cubic-equivalent q, M_LSWT from its header). In LSWT the weight of the
positive pole of the global-frame tr L is one sublattice moment.
"""
import math
import pathlib
import re
import sys

MEV_PER_RY = 13605.693
canon = lambda q: tuple(sorted((round(abs(float(x)), 6) for x in q), reverse=True))

oracle = pathlib.Path(sys.argv[1]).read_text().splitlines()
m_lswt = float(re.search(r"\|M\| = ([0-9.]+)", oracle[2]).group(1))
lswt = {canon(l.split()[:3]): float(l.split()[5]) / MEV_PER_RY for l in oracle if not l.startswith("#")}
print(f"# M_LSWT = {m_lswt} mu_B; columns: q_direct(3) | eta (Ry) | omega_p (Ry) | omega_p/omega_ad | W_b | W_c | W_b/M_J | W_c/M_J | M_eff (mu_B)")
for d in map(pathlib.Path, sys.argv[2:]):
    eta = float(re.search(r"^eta = (\S+)", (d / "input_driver.nml").read_text(), flags=re.M).group(1))
    state = [l.split() for l in (d / "spin_response_state.dat").read_text().splitlines()]
    m_j = abs(float(next(l for l in state if l and l[0] == "1")[2]))
    qs = [tuple(float(x) for x in l.split()) for l in (d / "qpath.dat").read_text().splitlines()]
    ad = [l.split() for l in (d / "spin_response_adiabatic.dat").read_text().splitlines() if not l.startswith("#")]
    for k, (q, a) in enumerate(zip(qs, ad), 1):
        rows = [[float(x) for x in l.split()] for l in (d / f"spin_response_q{k:03d}.dat").read_text().splitlines() if not l.startswith("#")]
        w, L = [r[0] for r in rows], [r[5] for r in rows]
        ip = max((i for i in range(len(w)) if w[i] > 0), key=lambda i: L[i])
        near = [i for i in range(len(w)) if abs(w[i] - w[ip]) <= 10 * eta]
        area = sum((w[i + 1] - w[i]) * (L[i] + L[i + 1]) / 2 for i in near[:-1])
        wb, wc = math.pi * eta * L[ip], area / (2 / math.pi * math.atan(10))
        ad_w = max(float(x) for x in a[7:-1])
        print(f"{q[0]:8.5f} {q[1]:8.5f} {q[2]:8.5f} | {eta:7.1e} | {w[ip]:.5f} | {w[ip] / ad_w:6.3f} | {wb:7.4f} | {wc:7.4f} | "
              f"{wb / m_j:6.3f} | {wc / m_j:6.3f} | {m_lswt * lswt[canon(q)] / ad_w:6.3f}   (M_J = {m_j:.5f}, edge: "
              f"{'yes' if ip in (0, len(w) - 1) else 'no'})")
