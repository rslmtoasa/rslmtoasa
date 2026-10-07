"""omega_s.py <runs dir> [--weight W] [--out FILE]: static delta_q = 1 + U Re tr chi0(q, 0), omega_s = U M delta_q and the stiffness D.
Reads every finished <runs dir>/st_N<N>_* directory (run3a.py --static). U = U_Juelich, M = Juelich moment, from that run's state.
--weight W (default 2.5447, the q = 0 pole weight) adds omega_s^W = U W delta_q. D = omega_s/|q|^2 (meV Angstrom^2): raw at the smallest
xi, and the intercept D0 of a quadratic fit of D(|q|) over xi <= 0.2 (a model, not a limit; needs >= 4 points)."""
import glob, os, re, sys
import numpy as np
from state import RY_MEV, read_state, read_q, read_dispersion

W = float(sys.argv[sys.argv.index("--weight") + 1]) if "--weight" in sys.argv else 2.5447
rows = {}
for d in sorted(glob.glob(sys.argv[1] + "/st_N*")):
    if not os.path.exists(f"{d}/spin_response_state.dat"): continue
    n = int(re.search(r"st_N(\d+)", d).group(1))
    st, disp = read_state(d), read_dispersion(d)
    for iq in range(1, len(disp) + 1):
        xi = np.linalg.norm(disp[iq - 1, 3:6]); dl = 1.0 + st["u_juelich"] * read_q(d, iq)[1]
        rows.setdefault(n, []).append((xi, disp[iq - 1, 6], dl, st["u_juelich"] * st["m_juelich"] * dl * RY_MEV, st))
lines = [f"# N xi |q|(1/A) delta_q omega_s(meV) D(meV A^2) omega_s^W(meV) D^W(meV A^2), W = {W}"]
for n in sorted(rows):
    r = sorted((x for x in rows[n] if x[0] > 0), key=lambda x: x[0])
    for xi, q, dl, w, st in r:
        k = W / st["m_juelich"]
        lines.append(f"{n} {xi:.4f} {q:.6f} {dl:.6e} {w:.5f} {w/q**2:.3f} {w*k:.5f} {w*k/q**2:.3f}")
    s = [x for x in r if x[0] <= 0.2 + 1e-9]
    txt = f"# N={n} U={r[0][4]['u_juelich']:.8f} M={r[0][4]['m_juelich']:.6f}: raw D(xi={r[0][0]:.3f}) = {r[0][3]/r[0][1]**2:.2f} (D^W {r[0][3]/r[0][1]**2*W/r[0][4]['m_juelich']:.2f})"
    if len(s) >= 4:
        c = np.polyfit([x[1] for x in s], [x[3] / x[1] ** 2 for x in s], 2)[::-1]
        txt += f"; D0 = {c[0]:.2f} (D0^W {c[0]*W/r[0][4]['m_juelich']:.2f}) from {len(s)} points xi <= 0.2"
    lines.append(txt)
print("\n".join(lines))
if "--out" in sys.argv: open(sys.argv[sys.argv.index("--out") + 1], "w").write("\n".join(lines) + "\n")
