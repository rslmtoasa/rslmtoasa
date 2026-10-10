"""fcc Ni spin-response dispersion against the LKAG + LSWT reference, along the reference q path.

Usage: run_dispersion.py <rslmto.x> <deck directory> <scratch directory> <ni_lswt.dat> [--path full|short] [--mesh N] [--eta E]
The deck directory (tests/spin_response/ni_fcc) is only read. The q points are the direct coordinates of the rows of
ni_lswt.dat: all 89 for --path full, every fourth row (23, with every labelled point) for --path short.
The window is 0 to 0.05 Ry in steps of eta/4. Prints, and writes to <scratch>/dispersion_compare.dat, the LSWT energy and
the dominant tr L peak per q. No tolerance is applied.
"""
import argparse
import math
import pathlib
import re
import shutil
import subprocess
import sys

MEV_PER_RY = 13605.693
OMEGA_MAX = 0.05

ap = argparse.ArgumentParser()
ap.add_argument("binary")
ap.add_argument("deck")
ap.add_argument("scratch")
ap.add_argument("oracle")
ap.add_argument("--path", choices=["full", "short"], default="short")
ap.add_argument("--mesh", type=int, default=24)
ap.add_argument("--eta", type=float, default=2.0e-3)
a = ap.parse_args()

rows = [l.split() for l in pathlib.Path(a.oracle).read_text().splitlines() if not l.startswith("#")]
if a.path == "short":
    rows = rows[::4]
q = [[float(x) for x in r[0:3]] for r in rows]
lswt_mev = [float(r[7]) for r in rows]
label = [r[9] if len(r) > 9 else "" for r in rows]

shutil.rmtree(a.scratch, ignore_errors=True)
shutil.copytree(a.deck, a.scratch)
scratch = pathlib.Path(a.scratch)
text = pathlib.Path(a.deck, "input.nml").read_text()
old = "pre_processing = 'bravais'"
if old not in text:
    sys.exit("run_dispersion.py: " + old + " not found in " + a.deck + "/input.nml")
for key in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{key} = \d+$", f"{key} = {a.mesh}", text, flags=re.M)
    if n != 1:
        sys.exit("run_dispersion.py: " + key + " not found in " + a.deck + "/input.nml")
n_omega = math.ceil(OMEGA_MAX / (a.eta / 4) - 1e-9) + 1
scratch.joinpath("qpath.dat").write_text("".join(f"{x!r} {y!r} {z!r}\n" for x, y, z in q))
block = (f"\n&spin_response\nmethod = 'juelich'\nq_file = 'qpath.dat'\nomega_min = 0.0\nomega_max = {OMEGA_MAX}\n"
         f"n_omega = {n_omega}\neta = {a.eta:g}\n/\n")
nml = scratch / "input_driver.nml"
nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + block)
subprocess.run([a.binary, nml.name], cwd=scratch, check=True, stdout=open(scratch / "driver.log", "w"), stderr=subprocess.STDOUT)

# the dispersion file has no space before an n/a in the last two columns
disp = [l.replace("n/a", " n/a").split() for l in (scratch / "spin_response_dispersion.dat").read_text().splitlines()
        if not l.startswith("#")]
if len(disp) != len(q):
    sys.exit(f"run_dispersion.py: {len(disp)} dispersion rows for {len(q)} q")
out = ["# q_direct(3) | LSWT [meV] | peak [meV] | (peak - LSWT)/LSWT | label ; mesh %d^3, eta %g Ry, %s path" % (a.mesh, a.eta, a.path)]
for qi, ref, d, lab in zip(q, lswt_mev, disp, label):
    peak = float(d[7]) * MEV_PER_RY if d[7] != "n/a" else None
    dev = f"{(peak - ref) / ref:+10.3f}" if peak is not None and ref > 0.0 else f"{'n/a':>10}"
    out.append(f"{qi[0]:9.5f} {qi[1]:9.5f} {qi[2]:9.5f} {ref:10.3f} {'n/a' if peak is None else format(peak, '10.3f'):>10} {dev}  {lab}")
(scratch / "dispersion_compare.dat").write_text("\n".join(out) + "\n")
print("\n".join(out))
