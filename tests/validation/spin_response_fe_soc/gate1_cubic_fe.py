"""Cubic bcc Fe null test: the q = 0 static adiabatic frequency with SOC against the SOC strength.

Usage: gate1_cubic_fe.py <rslmto.x> <deck directory> <scratch directory> <nk> [scale ...]
Runs nsp = 2 on the Fe deck with xi_p and xi_d times each scale (default 1, 1/2, 1/4) on an nk^3 mesh, static, q = 0, and prints
omega_adiabatic and omega_force of _adiabatic.dat. Cubic anisotropy is fourth order in xi, so omega(xi) is expected to fall as xi^4.
"""
import math
import pathlib
import re
import shutil
import subprocess
import sys

BLOCK = "\n&spin_response\nn_q = 1\nq_list(:, 1) = 0.0, 0.0, 0.0\nomega_min = 0.0\nn_omega = 1\neta = 0.0\n/\n"

binary, deck, scratch, nk = sys.argv[1:5]
scales = [float(eval(s)) for s in sys.argv[5:]] or [1.0, 0.5, 0.25]
text = pathlib.Path(deck, "input.nml").read_text()
old = "pre_processing = 'bravais'"
if old not in text:
    sys.exit("gate1_cubic_fe.py: " + old + " not found in " + deck + "/input.nml")
for key in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{key} = \d+$", f"{key} = {nk}", text, flags=re.M)
    if n != 1:
        sys.exit("gate1_cubic_fe.py: " + key + " not found")
text, n = re.subn(r"^nsp = 1", "nsp = 2", text, flags=re.M)
if n != 1:
    sys.exit("gate1_cubic_fe.py: nsp = 1 not found")
print(f"mesh {nk}^3, q = 0, omega = 0, eta = 0   columns: scale  omega_adiabatic[Ry]  omega_force[Ry]  slope_adiabatic  slope_force")
prev = None
for scale in scales:
    run = pathlib.Path(f"{scratch}_{scale:g}")
    shutil.rmtree(run, ignore_errors=True)
    shutil.copytree(deck, run)
    atom = run / "Fe.nml"
    body = atom.read_text()
    for key in ("xi_p", "xi_d"):
        m = re.search(rf"^ {key} = (\S+), (\S+)", body, flags=re.M)
        body = body.replace(m.group(0), f" {key} = {float(m.group(1)) * scale!r}, {float(m.group(2)) * scale!r}")
    atom.write_text(body)
    nml = run / "input_driver.nml"
    nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + BLOCK)
    subprocess.run([binary, nml.name], cwd=run, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    row = (run / "spin_response_adiabatic.dat").read_text().splitlines()[1].split()
    w, wf = float(row[7]), float(row[9])
    slope = ["n/a", "n/a"] if prev is None else [f"{math.log(abs(a / b)) / math.log(scale / prev[0]):.3f}" for a, b in ((w, prev[1]), (wf, prev[2]))]
    print(f"{scale:8g} {w:+.6e} {wf:+.6e}  {slope[0]:>8} {slope[1]:>8}")
    prev = (scale, w, wf)
