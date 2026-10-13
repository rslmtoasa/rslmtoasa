"""bcc Cr antiferromagnet spin response along Gamma-X-M-Gamma-R-X | M-R (simple cubic cell, direct = Cartesian in 2 pi / a).

Usage: run_path.py <rslmto.x> <deck directory> <scratch directory> --mesh N [--part i/n] [--dynamic OMEGA_MAX N_OMEGA ETA]
The deck directory (tests/spin_response/cr_bcc) is only read. Static mode (default): n_omega = 1, omega = 0, eta = 0.
--dynamic: omega from -OMEGA_MAX to OMEGA_MAX (Ry), N_OMEGA points, eta in Ry. The q points are steps of 1/32 along the path,
plus 1/64 steps on Gamma-X up to 1/8 when N is a multiple of 64. --part i/n runs the i-th of n contiguous slices of the q list
(single-thread processes for the static mode). Writes <scratch>/qpath_labels.dat (index, direct q, segment label) for the whole list.
"""
import argparse
import pathlib
import re
import shutil
import subprocess
import sys
from fractions import Fraction as F

ap = argparse.ArgumentParser()
ap.add_argument("binary")
ap.add_argument("deck")
ap.add_argument("scratch")
ap.add_argument("--mesh", type=int, required=True)
ap.add_argument("--part", default="1/1")
ap.add_argument("--dynamic", nargs=3, type=float, metavar=("OMEGA_MAX", "N_OMEGA", "ETA"))
a = ap.parse_args()

G, X, M, R = (F(0),) * 3, (F(1, 2), F(0), F(0)), (F(1, 2), F(1, 2), F(0)), (F(1, 2),) * 3
SEGMENTS = [("G-X", G, X), ("X-M", X, M), ("M-G", M, G), ("G-R", G, R), ("R-X", R, X), ("M-R", M, R)]
points = []   # (q, label); a repeated q keeps its first label
for name, p0, p1 in SEGMENTS:
    for i in range(17):
        q = tuple(x0 + (x1 - x0) * i / 16 for x0, x1 in zip(p0, p1))
        if all(q != p[0] for p in points):
            points.append((q, name))
if a.mesh % 64 == 0:
    for j in range(1, 8):
        q = (F(j, 64), F(0), F(0))
        if all(q != p[0] for p in points):
            points.append((q, "ladder"))
assert all((x * a.mesh).denominator == 1 for q, _ in points for x in q), "q not commensurate with the mesh"

i, n = (int(s) for s in a.part.split("/"))
mine = points[(i - 1) * len(points) // n: i * len(points) // n]

shutil.rmtree(a.scratch, ignore_errors=True)
scratch = pathlib.Path(a.scratch)
scratch.mkdir(parents=True)
for f in pathlib.Path(a.deck).glob("*.nml"):
    shutil.copy(f, scratch)
text = pathlib.Path(a.deck, "input.nml").read_text()
old = "pre_processing = 'bravais'"
if old not in text:
    sys.exit("run_path.py: " + old + " not found in " + a.deck + "/input.nml")
for key in ("nk1", "nk2", "nk3"):
    text, k = re.subn(rf"^{key} = \d+$", f"{key} = {a.mesh}", text, flags=re.M)
    if k != 1:
        sys.exit("run_path.py: " + key + " not found in " + a.deck + "/input.nml")
scratch.joinpath("qpath.dat").write_text("".join(f"{float(x)!r} {float(y)!r} {float(z)!r}\n" for (x, y, z), _ in mine))
scratch.joinpath("qpath_labels.dat").write_text("# index qx qy qz (direct) segment\n" + "".join(
    f"{k + 1} {float(x)!r} {float(y)!r} {float(z)!r} {lab}\n" for k, ((x, y, z), lab) in enumerate(points)))
if a.dynamic:
    w, nw, eta = a.dynamic
    window = f"omega_min = {-w!r}\nomega_max = {w!r}\nn_omega = {int(nw)}\neta = {eta!r}\n"
else:
    window = "n_omega = 1\neta = 0.0\n"
nml = scratch / "input_driver.nml"
nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'")
               + f"\n&spin_response\nmethod = 'juelich'\nq_file = 'qpath.dat'\n{window}/\n")
subprocess.run([a.binary, nml.name], cwd=scratch, check=True, stdout=open(scratch / "driver.log", "w"), stderr=subprocess.STDOUT)
