"""classify.py <runs dir> [--eta 1e-3] [--weight W]: classify every interior local maximum of tr L in the 60^3 scans.
Scan directories <runs dir>/sc60_e<eta>_*, run3a.py with a window; omega_s from <runs dir>/st_N60_* (static, eta = 0).
--weight W replaces M by the measured q = 0 pole weight W (omega_s = U W delta_q).
A maximum counts if its height is >= 10% of the dominant one (parabolic refinement); it is the RPA branch if it lies within 10% of
omega_s = U M delta_q, otherwise 'other' (not diagnosed). The nearest maximum and its deviation are printed for every xi."""
import glob, sys
import numpy as np
from state import RY_MEV, read_state, read_q, read_dispersion

eta = sys.argv[sys.argv.index("--eta") + 1] if "--eta" in sys.argv else "1e-3"
HEIGHT, BAND = 0.10, 0.10
weight = float(sys.argv[sys.argv.index("--weight") + 1]) if "--weight" in sys.argv else None
ws = {}
for d in glob.glob(sys.argv[1] + "/st_N60_*"):
    st, disp = read_state(d), read_dispersion(d)
    for iq in range(1, len(disp) + 1):
        ws[round(float(np.linalg.norm(disp[iq - 1, 3:6])), 4)] = st["u_juelich"] * (weight or st["m_juelich"]) * (1.0 + st["u_juelich"] * read_q(d, iq)[1]) * RY_MEV
res = []
for d in glob.glob(f"{sys.argv[1]}/sc60_e{eta}_*"):
    disp = read_dispersion(d)
    for iq in range(1, len(disp) + 1):
        xi = round(float(np.linalg.norm(disp[iq - 1, 3:6])), 4)
        q = read_q(d, iq); om, trl = q[:, 0], q[:, 5]
        idx = [i for i in range(1, len(om) - 1) if trl[i] > trl[i - 1] and trl[i] > trl[i + 1]]
        pk = []
        for i in idx:
            h = 0.5 * (om[i + 1] - om[i - 1]); den = 2 * (trl[i - 1] - 2 * trl[i] + trl[i + 1])
            pk.append((om[i] + h * (trl[i - 1] - trl[i + 1]) / den, trl[i]))
        top = max(p[1] for p in pk)
        pk = [(w * RY_MEV, h / top) for w, h in pk if h >= HEIGHT * top]
        res.append((xi, ws[xi], pk))
n_rpa = n_other = 0
print(f"# eta={eta}  xi omega_s(meV) | maxima: omega(meV) height/dominant class | nearest-to-omega_s deviation")
for xi, w, pk in sorted(res):
    cls = ["RPA" if abs(p[0] - w) <= BAND * w else "other" for p in pk]
    n_rpa += cls.count("RPA"); n_other += cls.count("other")
    near = min(pk, key=lambda p: abs(p[0] - w))
    print(f"{xi:.4f} {w:8.2f} | " + "  ".join(f"{p[0]:.1f}({p[1]:.2f}){c}" for p, c in zip(pk, cls)) + f" | {near[0]/w-1:+.3f}")
print(f"# maxima: {n_rpa} RPA branch, {n_other} other, {len(res)} q")
