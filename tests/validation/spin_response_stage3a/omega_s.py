"""omega_s.py <runs dir> [--out FILE]: static delta_q = 1 + U Re tr chi0(q, 0), omega_s = U M delta_q and the stiffness D.
Reads every <runs dir>/st_N<N>_* directory (run3a.py --static). U = U_Juelich, M = Juelich moment, both from that run's state.
D = omega_s/|q|^2 at the smallest q, and the least-squares a in omega_s = a |q|^2 + b |q|^4 over xi <= 0.2 (meV Angstrom^2)."""
import glob, re, sys
import numpy as np
from state import RY_MEV, read_state, read_q, read_dispersion

rows = {}
for d in sorted(glob.glob(sys.argv[1] + "/st_N*")):
    n = int(re.search(r"st_N(\d+)", d).group(1))
    st, disp = read_state(d), read_dispersion(d)
    for iq in range(1, len(disp) + 1):
        xi = np.linalg.norm(disp[iq - 1, 3:6])
        delta = 1.0 + st["u_juelich"] * read_q(d, iq)[1]
        rows.setdefault(n, []).append((xi, disp[iq - 1, 6], delta, st["u_juelich"] * st["m_juelich"] * delta * RY_MEV, st))
lines = ["# N xi |q|(1/A) delta_q omega_s(meV) omega_s/|q|^2(meV A^2)"]
for n in sorted(rows):
    r = sorted((x for x in rows[n] if x[0] > 0), key=lambda x: x[0])
    for xi, q, dl, w, st in r:
        lines.append(f"{n} {xi:.4f} {q:.6f} {dl:.6e} {w:.5f} {w/q**2:.3f}")
    s = [x for x in r if x[0] <= 0.2 + 1e-9]
    if len(s) >= 3:
        q2 = np.array([x[1] ** 2 for x in s]); w = np.array([x[3] for x in s])
        a, b = np.linalg.lstsq(np.stack([q2, q2 ** 2], 1), w, rcond=None)[0]
        lines.append(f"# N={n} U={r[0][4]['u_juelich']:.8f} M={r[0][4]['m_juelich']:.6f}: raw D(xi={r[0][0]:.3f}) = {r[0][3]/r[0][1]**2:.2f}, "
                     f"fit D = {a:.2f} meV A^2, b = {b:.1f} meV A^4 over {len(s)} points xi <= 0.2")
print("\n".join(lines))
if "--out" in sys.argv: open(sys.argv[sys.argv.index("--out") + 1], "w").write("\n".join(lines) + "\n")
