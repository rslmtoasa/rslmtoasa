"""Inversion null test: with SOC, centrosymmetric bcc Fe has omega(q) = omega(-q) on a Gamma-centred full-zone mesh.

Usage: inversion_odd.py <rslmto.x> <deck directory> <scratch directory> <tolerances.nml>
Runs nsp = 2 on the Fe deck at 16^3, static (omega = 0, eta = 0), for +-q pairs along [100] and [110], and checks that the odd
part of omega_adiabatic and of omega_force is below inversion_odd_abs (Ry). The q are in units of 2 pi/a; the primitive cell is
a1 = (-1, 1, 1)/2, a2 = (1, -1, 1)/2, a3 = (1, 1, -1)/2, so q_direct = (q . a1, q . a2, q . a3).
Catches an error that makes the response depend on the sign of k (a conjugation or phase-convention slip in H(k) or in the k+q
endpoint) and a mesh that is not symmetric under k -> -k.
"""
import pathlib
import re
import shutil
import subprocess
import sys

DIRECTIONS = {"100": (1, 0, 0), "110": (1, 1, 0)}
H = {"100": (0.125, 0.25, 0.1), "110": (0.0625, 0.125, 0.1)}
CELL = ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5))


def tolerance(path, key):
    m = re.search(rf"^\s*{key}\s*=\s*([0-9.eEdD+-]+)", pathlib.Path(path).read_text(), flags=re.M)
    if not m or not float(m.group(1).replace("d", "e").replace("D", "e")) > 0:
        sys.exit(f"inversion_odd.py: {path}: key {key} missing or <= 0")
    return float(m.group(1).replace("d", "e").replace("D", "e"))


binary, deck, scratch, tolfile = sys.argv[1:5]
tol = tolerance(tolfile, "inversion_odd_abs")
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
text = pathlib.Path(deck, "input.nml").read_text()
old = "pre_processing = 'bravais'"
if old not in text:
    sys.exit("inversion_odd.py: " + old + " not found in " + deck + "/input.nml")
for key in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{key} = \d+$", f"{key} = 16", text, flags=re.M)
    if n != 1:
        sys.exit("inversion_odd.py: " + key + " not found in " + deck + "/input.nml")
text, n = re.subn(r"^nsp = 1", "nsp = 2", text, flags=re.M)
if n != 1:
    sys.exit("inversion_odd.py: nsp = 1 not found in " + deck + "/input.nml")

pairs = [(name, h) for name, hs in H.items() for h in hs]
rows = []
for name, h in pairs:
    for sign in (1, -1):
        rows.append([sign * h * sum(c * d for c, d in zip(a, DIRECTIONS[name])) for a in CELL])
block = "\n&spin_response\nmethod = 'juelich'\nn_q = " + str(len(rows)) + "\n"
block += "".join(f"q_list(:, {i + 1}) = {q[0]!r}, {q[1]!r}, {q[2]!r}\n" for i, q in enumerate(rows))
block += "omega_min = 0.0\nn_omega = 1\neta = 0.0\n/\n"
nml = pathlib.Path(scratch, "input_driver.nml")
nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + block)
subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, "driver.log"), "w"),
               stderr=subprocess.STDOUT)

if "u_m0_nosoc" not in pathlib.Path(scratch, "spin_response_state.dat").read_text():
    sys.exit("inversion_odd.py: the run did not take the SOC branch (no u_m0_nosoc in spin_response_state.dat)")
data = [ln.split() for ln in pathlib.Path(scratch, "spin_response_adiabatic.dat").read_text().splitlines()[1:]]
if len(data) != len(rows):
    sys.exit(f"inversion_odd.py: {len(data)} rows in spin_response_adiabatic.dat, expected {len(rows)}")
print("direction  h[2pi/a]  omega_RPA(+q)[Ry]  odd_RPA[Ry]  odd_force[Ry]  (tolerance " + f"{tol:.1e} Ry)")
worst = 0.0
for i, (name, h) in enumerate(pairs):
    p, m = data[2 * i], data[2 * i + 1]
    wp, wm, fp, fm = float(p[7]), float(m[7]), float(p[9]), float(m[9])
    if not abs(wp) > 1.0e-6:
        sys.exit(f"inversion_odd.py: omega_RPA = {wp:.3e} Ry at {name} h = {h}: not a magnon frequency, the check would be vacuous")
    worst = max(worst, abs(wp - wm), abs(fp - fm))
    print(f"{name:>9} {h:9.4f} {wp:+.9e} {wp - wm:+.2e} {fp - fm:+.2e}")
print(f"max |odd| over omega_RPA and omega_force: {worst:.2e} Ry")
sys.exit(0 if worst <= tol else 1)
