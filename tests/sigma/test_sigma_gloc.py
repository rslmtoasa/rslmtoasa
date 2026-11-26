#!/usr/bin/env python3
"""do_sigma deck test on bcc Fe (nb = 18, one atom type): G_ii(z) from sigma.dat against closed forms.

A: Sigma = sigma*I, sigma = 0.02 - 0.05i Ry, at offset 0.0 + 0.50i.   B: Sigma = 0 at offset A - sigma = -0.02 + 0.55i.
Closed form: (z - H - sigma)^-1 = (z - sigma - H)^-1, so block A must equal block B.
  Bug caught: a wrong on-site term (missing -i Im b, missing 1/a, Sigma on the wrong sites) or a wrong sign of the offset
  (z = E_F - offset) makes A differ from B, because A runs with a complex centre b and on-site blocks, B with real b and none.
C: Sigma = 0 at offset 0.0 + 20.0i: G -> I/z, so |z G - I|_max < 0.1 with z = E_F + 20i.
  Bug caught: a missing 1/a factor in the sum or a wrong E_F offset (E_F = -0.06 Ry here; only a gross error is caught at this tolerance) moves z G away from I.
"""
import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL_AB, TOL_C = 1e-8, 0.1


def read_gloc(path):
    lines = path.read_text().splitlines()
    nz, nrec, nb = map(int, lines[1].split())
    pos, blocks = 2, []
    for _ in range(nz):
        pos += 1
        rows = []
        for _ in range(nb):
            rows.append(lines[pos])
            pos += 1
        blocks.append(None if any("n/a" in r for r in rows) else
                      np.array([[float(v[0]) + 1j * float(v[1]) for v in zip(r.split()[0::2], r.split()[1::2])] for r in rows]))
    return blocks


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    work = Path(a.workdir)
    shutil.rmtree(work, ignore_errors=True)
    shutil.copytree(HERE, work, ignore=shutil.ignore_patterns("__pycache__", "*.py"))
    t0 = time.time()
    run = subprocess.run([a.binary], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300)
    elapsed = time.time() - t0
    (work / "run.log").write_text(run.stdout)
    if run.returncode != 0:
        print(run.stdout[-3000:])
        raise SystemExit(f"ERROR: run exited {run.returncode}")
    fermi = float(re.findall(r"Free Fermi energy:\s*([-+0-9.eE]+)", run.stdout)[-1])
    ga, gb, gc = read_gloc(work / "gloc_sigma.dat")
    if ga is None or gb is None or gc is None:
        raise SystemExit("ERROR: gloc_sigma.dat holds n/a for a block that must converge")
    ab = np.abs(ga - gb).max() / np.abs(gb).max()
    c = np.abs((fermi + 20j) * gc - np.eye(gc.shape[0])).max()
    print(f"A vs B: max|G_A - G_B|/max|G_B| = {ab:.3e} (tol {TOL_AB:g})")
    print(f"C: max|z G_C - I| = {c:.3e} (tol {TOL_C:g}); runtime {elapsed:.1f} s")
    if not (ab <= TOL_AB and c < TOL_C):
        raise SystemExit("FAIL")
    print("PASS")


if __name__ == "__main__":
    sys.exit(main())
