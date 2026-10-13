"""Tabulate static Cr runs made by run_path.py: adiabatic frequency, smallest eigenvalue of I + chi0(q, 0) U, the Gamma-X exponent.

Usage: report_path.py <dir[,dir...]> [<dir[,dir...]> ...]
One argument per mesh; the directories of one argument are the --part slices of that run. Prints one row per q (path order of
qpath_labels.dat of the first mesh) with the positive adiabatic branch and the smallest eigenvalue for every mesh, then per
mesh: the largest |w+ + w-| / w+ over q, the smallest eigenvalue and where it is, and the local exponents
ln(w2/w1) / ln(q2/q1) of successive points on Gamma-X with q <= 1/8.
"""
import math
import pathlib
import sys


def rows(path):
    return [[float(x) for x in l.replace("n/a", " n/a").replace("n/a", "nan").split()] for l in path.read_text().splitlines()
            if not l.startswith("#")]


def load(dirs):
    out = {}   # q (direct, rounded) -> (eigenvalues of omega, min eig of I + chi0 U)
    for d in dirs:
        d = pathlib.Path(d)
        q = [tuple(round(float(x), 9) for x in l.split()) for l in (d / "qpath.dat").read_text().splitlines()]
        adia, disp = rows(d / "spin_response_adiabatic.dat"), rows(d / "spin_response_dispersion.dat")
        assert len(q) == len(adia) == len(disp), f"{d}: row counts differ"
        for qi, a, s in zip(q, adia, disp):
            out[qi] = (a[7:-1], s[-1])
    return out


meshes = [load(arg.split(",")) for arg in sys.argv[1:]]
labels = [l.split() for l in pathlib.Path(sys.argv[1].split(",")[0], "qpath_labels.dat").read_text().splitlines()[1:]]
order = [(tuple(round(float(x), 9) for x in l[1:4]), l[4]) for l in labels]
print("# q_direct(3) | segment | w+ (Ry) per mesh | min eig(I+chi0 U) per mesh")
for q, lab in order:
    cells = [m.get(q) for m in meshes]
    w = " ".join(f"{max(c[0]):12.6e}" if c else f"{'n/a':>12}" for c in cells)
    e = " ".join(f"{c[1]:12.5e}" if c else f"{'n/a':>12}" for c in cells)
    print(f"{q[0]:8.5f} {q[1]:8.5f} {q[2]:8.5f} {lab:7s} | {w} | {e}")
for k, m in enumerate(meshes):
    pair = max(abs(min(e) + max(e)) / max(e) for q, (e, _) in m.items() if any(q))
    qmin = min(m, key=lambda q: m[q][1])
    print(f"# mesh {k + 1}: Gamma w = {m[(0.0, 0.0, 0.0)][0]}; max |w+ + w-|/w+ (q /= 0) = {pair:.2e}; min eig = {m[qmin][1]:.5e} at q = {qmin}")
    gx = sorted((q, m[q][0]) for q in m if q[1] == 0 and q[2] == 0 and 0 < q[0] <= 0.125)
    for (q1, e1), (q2, e2) in zip(gx, gx[1:]):
        print(f"#   p({q1[0]:.5f} -> {q2[0]:.5f}) = {math.log(max(e2) / max(e1)) / math.log(q2[0] / q1[0]):.3f}")
