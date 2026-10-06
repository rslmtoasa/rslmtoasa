"""Rerun the frozen 24^3 Fe response and compare it with tests/spin_response/baseline/.

Usage: baseline_repro.py [--freeze] <rslmto.x> <deck directory> <scratch directory> <baseline directory> <tolerances.nml>
The deck directory is only read. --freeze writes the run's q files and state file into the baseline directory.
Every column of every q file is compared normwise, every numeric scalar of the state file relatively.
"""
import math
import pathlib
import re
import shutil
import subprocess
import sys

# xi = 0.25, 1/3, 0.5 as direct (t, -t, t) with t = 3/24, 4/24, 6/24: mesh points of the 24^3 mesh. eta = 2e-3 Ry, 101 points.
Q_T = [3 / 24, 4 / 24, 6 / 24]
FILES = ["spin_response_q001.dat", "spin_response_q002.dat", "spin_response_q003.dat", "spin_response_state.dat"]


def run(binary, deck, scratch):
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.copytree(deck, scratch)
    text = pathlib.Path(deck, "input.nml").read_text()
    old = "pre_processing = 'bravais'"
    if old not in text:
        sys.exit("baseline_repro.py: " + old + " not found in " + deck + "/input.nml")
    for key in ("nk1", "nk2", "nk3"):
        text, n = re.subn(rf"^{key} = 12$", f"{key} = 24", text, flags=re.M)
        if n != 1:
            sys.exit("baseline_repro.py: " + key + " = 12 not found in " + deck + "/input.nml")
    q = "".join(f"q_list(:, {i + 1}) = {t!r}, {-t!r}, {t!r}\n" for i, t in enumerate(Q_T))
    block = f"\n&spin_response\nmethod = 'juelich'\nn_q = 3\n{q}omega_min = 0.0\nomega_max = 0.05\nn_omega = 101\neta = 2.0e-3\n/\n"
    nml = pathlib.Path(scratch, "input_driver.nml")
    nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + block)
    subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, "driver.log"), "w"),
                   stderr=subprocess.STDOUT)


def tolerance(path, key):
    m = re.search(rf"^\s*{key}\s*=\s*([0-9.eEdD+-]+)", pathlib.Path(path).read_text(), flags=re.M)
    if not m or not float(m.group(1).replace("d", "e").replace("D", "e")) > 0:
        sys.exit(f"baseline_repro.py: {path}: key {key} missing or <= 0")
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def rel(a, b):
    return abs(a - b) / abs(b) if b != 0.0 else abs(a - b)


def compare_q(new, ref):
    a = [[float(x) for x in l.split()] for l in new.read_text().splitlines() if not l.startswith("#")]
    b = [[float(x) for x in l.split()] for l in ref.read_text().splitlines() if not l.startswith("#")]
    if len(a) != len(b):
        return math.inf
    worst = 0.0
    for c in range(len(b[0])):
        num = math.sqrt(sum((x[c] - y[c]) ** 2 for x, y in zip(a, b)))
        den = math.sqrt(sum(y[c] ** 2 for y in b))
        worst = max(worst, num / den if den != 0.0 else num)
    return worst


def compare_state(new, ref):
    a, b = new.read_text().splitlines(), ref.read_text().splitlines()
    if len(a) != len(b):
        return math.inf
    worst = 0.0
    for la, lb in zip(a, b):
        ta, tb = la.split(), lb.split()
        if len(ta) != len(tb):
            return math.inf
        for x, y in zip(ta, tb):
            try:
                worst = max(worst, rel(float(x), float(y)))
            except ValueError:
                if x != y:
                    return math.inf
    return worst


def main():
    args = sys.argv[1:]
    freeze = "--freeze" in args
    binary, deck, scratch, baseline, tolfile = [a for a in args if a != "--freeze"]
    run(binary, deck, scratch)
    if freeze:
        pathlib.Path(baseline).mkdir(exist_ok=True)
        for f in FILES:
            shutil.copy(pathlib.Path(scratch, f), pathlib.Path(baseline, f))
        return
    tol = tolerance(tolfile, "baseline_repro_rel")
    failed = False
    print("file | measured | tolerance")
    for f in FILES:
        new, ref = pathlib.Path(scratch, f), pathlib.Path(baseline, f)
        d = compare_state(new, ref) if f.endswith("state.dat") else compare_q(new, ref)
        print(f"{f} | {d:.3e} | {tol:.1e} | {'ok' if d <= tol else 'FAIL'}")
        failed = failed or not d <= tol
    sys.exit(1 if failed else 0)


main()
