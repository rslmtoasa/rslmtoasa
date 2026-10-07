"""moment_q0.py <runs dir>: zeroth moment of the Juelich chi0 at q = 0, integral of -Im tr chi0 / pi over the window (trapezoid), split
at omega = 0, from <runs dir>/m0_N* (run3a.py, window -w..w, eta, spacing eta/4). The Lorentzian tail outside the window is estimated
as 2 eta/(pi w) for transitions near omega = 0 (+-10% for transitions out to 1.5 Ry), applied as 1/(1 - 2 eta/(pi w))."""
import glob, math, sys
import numpy as np
from state import read_q, read_state

for d in sorted(glob.glob(sys.argv[1] + "/m0_N*")):
    q = read_q(d, 1); om = q[:, 0]; t = -q[:, 2] / math.pi
    eta, w = float(d.split("_e")[-1]), om[-1]
    tot = np.trapezoid(t, om); pos = np.trapezoid(t[om >= 0], om[om >= 0]); neg = np.trapezoid(t[om <= 0], om[om <= 0])
    c = 1 / (1 - 2 * eta / (math.pi * w))
    print(f"{d.split('/')[-1]}: raw {tot:.5f} (omega>0: {pos:.5f}, omega<0: {neg:.5f}), tail-corrected {tot*c:.5f} (factor {c:.5f}); M_juelich {read_state(d)['m_juelich']:.5f}")
