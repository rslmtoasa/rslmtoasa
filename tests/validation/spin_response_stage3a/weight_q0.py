"""weight_q0.py <runs dir>: Goldstone-pole weight at q = 0 from tr L, per eta, against the Juelich moment.
(b) M = pi * eta * tr L(omega = 0): the pole is a Lorentzian of width eta and weight M.
(c) M = integral of tr L over |omega| <= w, divided by (2/pi) arctan(w/eta), the Lorentzian fraction inside the window."""
import glob, math, sys
import numpy as np
from state import read_state, read_q

for d in sorted(glob.glob(sys.argv[1] + "/q0_e*"), reverse=True):
    eta = float(d.split("_e")[-1])
    st, q = read_state(d), read_q(d, 1)
    om, trl = q[:, 0], q[:, 5]
    i0 = int(np.argmin(abs(om)))
    assert om[i0] == 0.0
    mb = math.pi * eta * trl[i0]
    cols = [f"eta={eta:g}", f"M_juelich={st['m_juelich']:.5f}", f"(b) {mb:.5f} ({mb/st['m_juelich']-1:+.2%})"]
    for w in (0.0025, 0.005, 0.01):
        s = abs(om) <= w + 1e-12
        mc = np.trapezoid(trl[s], om[s]) / (2 / math.pi * math.atan(w / eta))
        cols.append(f"(c) w={w:g}: {mc:.5f} ({mc/st['m_juelich']-1:+.2%})")
    print("  ".join(cols))
