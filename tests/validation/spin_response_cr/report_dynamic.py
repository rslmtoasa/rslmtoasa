"""Compare the tr L peaks of a dynamic run_path.py run with the static adiabatic frequencies.

Usage: report_dynamic.py <dynamic scratch directory> <static scratch directory[,directory...]>
Per q of the dynamic run: the omega > 0 of the largest tr L, whether it is the last grid point, the static w+ (from the adiabatic
files of the static run, same q, same mesh), their ratio, and max |L(w) + L(-w)| / max |L| over the grid.
"""
import pathlib
import sys

dyn = pathlib.Path(sys.argv[1])
static = {}
for d in (pathlib.Path(s) for s in sys.argv[2].split(",")):
    q = [tuple(float(x) for x in l.split()) for l in (d / "qpath.dat").read_text().splitlines()]
    a = [l.split() for l in (d / "spin_response_adiabatic.dat").read_text().splitlines() if not l.startswith("#")]
    static.update({qi: max(float(x) for x in r[7:-1]) for qi, r in zip(q, a)})
print("# q_direct(3) | w(max tr L) (Ry) | at grid edge | static w+ (Ry) | ratio | max|L(w)+L(-w)|/max|L|")
for k, l in enumerate((dyn / "qpath.dat").read_text().splitlines(), 1):
    q = tuple(float(x) for x in l.split())
    rows = [[float(x) for x in r.split()] for r in (dyn / f"spin_response_q{k:03d}.dat").read_text().splitlines() if not r.startswith("#")]
    w, L = [r[0] for r in rows], [r[5] for r in rows]
    n = len(w)
    ip = max((i for i in range(n) if w[i] > 0), key=lambda i: L[i])
    anti = max(abs(L[i] + L[n - 1 - i]) for i in range(n)) / max(abs(x) for x in L)
    ws = static[q]
    print(f"{q[0]:8.5f} {q[1]:8.5f} {q[2]:8.5f} | {w[ip]:8.4f} | {'yes' if ip == n - 1 else 'no':3s} | {ws:10.6f} | "
          f"{w[ip] / ws if ws > 1e-6 else float('nan'):8.3f} | {anti:.2e}")
