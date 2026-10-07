"""wq_runs.py <runs dir>: runs for the pole weight W(q) at 60^3, Gamma-H, xi = 0.1, 0.2, 7/30, 1/3, one run per (xi, eta).
Window [0.5, 1.6] omega_s (omega_s = U M delta_q from the static 60^3 runs), spacing eta/8. eta from 2e-3 down to 1.25e-4 (halving)
with omega_s/25 <= eta <= omega_s, so that the peak is resolved and no finer than needed. Four runs at a time."""
import pathlib, subprocess, sys
from state import RY_MEV, omega_s_table

here = pathlib.Path(__file__).parent
ws = omega_s_table(sys.argv[1], 60)
jobs = []
for xi in (0.1, 0.2, 0.2333, 0.3333):
    w = ws[xi]
    xq = {0.1: 3, 0.2: 6, 0.2333: 7, 0.3333: 10}[xi] / 30
    for eta in (2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4):
        if w / 25 <= eta <= w:
            lo, hi = 0.5 * w, 1.6 * w
            jobs.append(["python3", str(here / "run3a.py"), "60", repr(xq), "--tag", f"wq_x{xi}_e{eta:g}", "--eta", repr(eta),
                         f"--window={lo:.6g},{hi:.6g},{round((hi - lo) / (eta / 8)) + 1}"])
running = []
for j in jobs:
    while len(running) >= 4:
        running = [p for p in running if p.poll() is None]
        if len(running) >= 4: running[0].wait()
    running.append(subprocess.Popen(j))
for p in running: p.wait()
