"""wq.py <runs dir>: pole weight W(q) = pi eta tr L(omega_peak) per eta, and the fit 1/(pi h) = (eta + Gamma)/W over the eta ladder.
Peak and height from the parabola through the maximum and its two neighbours (interior maximum required, else n/a).
W(0) = pi eta tr L(0) at eta = 1e-3 from <runs dir>/q0_e1e-3 (weight_q0.py estimator (b))."""
import glob, math, sys
import numpy as np
from state import RY_MEV, read_q, read_state, omega_s_table

ws = omega_s_table(sys.argv[1], 60)
q0 = read_q(f"{sys.argv[1]}/q0_e1e-3", 1)
w0 = math.pi * 1e-3 * q0[int(np.argmin(abs(q0[:, 0]))), 5]
print(f"# W(0) = {w0:.5f} (eta 1e-3); columns: eta(Ry) peak(meV) peak/omega_s(M) W=pi*eta*h W/W(0); fit; and W_pos = M peak/omega_s(M) = peak/(U delta_q), from the peak position only")
for xi in (0.1, 0.2, 0.2333, 0.3333):
    pts = []
    for d in sorted(glob.glob(f"{sys.argv[1]}/wq_x{xi}_e*"), key=lambda d: -float(d.split("_e")[-1])):
        eta = float(d.split("_e")[-1]); q = read_q(d, 1); om, t = q[:, 0], q[:, 5]
        i = int(np.argmax(t))
        if i in (0, len(t) - 1):
            print(f"xi={xi} eta={eta:g}  peak on the window edge: n/a"); continue
        den = t[i - 1] - 2 * t[i] + t[i + 1]; sh = 0.5 * (t[i - 1] - t[i + 1]) / den
        h = t[i] - 0.25 * (t[i - 1] - t[i + 1]) * sh; wp = om[i] + sh * (om[i + 1] - om[i])
        pts.append((eta, 1 / (math.pi * h)))
        print(f"xi={xi} eta={eta:g}  {wp*RY_MEV:8.3f} {wp/ws[xi]:6.3f} {math.pi*eta*h:8.4f} {math.pi*eta*h/w0:6.3f}  W_pos/W(0) = {read_state(d)['m_juelich']*wp/ws[xi]/w0:.3f}")
    if len(pts) >= 3:
        e, y = np.array(pts).T; a, b = np.polyfit(e, y, 1)
        print(f"xi={xi} fit over {len(pts)} eta: W = {1/a:.4f} (W/W(0) = {1/a/w0:.3f}), Gamma = {b/a*RY_MEV:.3f} meV")
