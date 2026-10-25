#!/usr/bin/env python3
"""AppGbtKspaceScSweep: frozen_magnon mode = 'scf' on the k-space path (unconstrained cone SCF at every q of a q file, bcc Fe,
tests/spin_response/fe_bcc potential, 8^3 mesh, 300 K, first order, theta = 30 deg, q = 0 and two q along Gamma-H with q/2 on the mesh).
Consistency, not validation: E_tot(q) of the sweep (frozen_magnon_kspace.dat) equals the energy of a separate single-q gbt_single_q SCF of the
same deck at each q (report.out), within scf_identity_ry; both sides converged; the sweep energies differ between q by more than that
tolerance (non-trivial). Bugs caught: q not applied inside the loop (every row at the first q), a state of the previous q leaking into the next SCF, an energy summed over the wrong atoms, a converged flag that is not set.
usage: test_kspace_scsweep.py <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import SR, lattice_vectors, tolerance

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
Q = [(0.0, 0.0, 0.0), (0.25, -0.25, 0.25), (0.5, -0.5, 0.5)]       # direct; q/2 on the 8^3 mesh
NK, THETA = 8, 30.0


def deck(dst, sweep, q=None):
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "fe_bcc", dst)
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["self"].update(nstep=150, conv_thr=5.0e-10)
    nml["control"]["nsp"] = 3
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK)
    nml["hamiltonian"] = {"theta_ss": THETA}
    if sweep:
        nml["calculation"]["post_processing"] = "frozen_magnon"
        nml["frozen_magnon"] = {"mode": "scf", "q_file": "q.dat", "q_coordinates": "direct"}
        (dst / "q.dat").write_text(f"{len(Q)}\n" + "".join(" ".join(map(str, p)) + "\n" for p in Q))
    else:
        nml["hamiltonian"].update(magnetic_representation="gbt_single_q", q_ss=list(np.linalg.inv(lattice_vectors("bccFe_prim")) @ np.array(q)))
    nml.write(str(dst / "input.nml"), force=True)
    r = subprocess.run([str(binary), "input.nml"], cwd=dst, capture_output=True, text=True)
    (dst / "run.log").write_text(re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr))
    assert r.returncode == 0, f"{dst}: exit code {r.returncode}"
    return dst


tol = tolerance("scf_identity_ry")
sweep = np.loadtxt(deck(scratch / "sweep", True) / "frozen_magnon_kspace.dat")      # q direct (3), q Cart (3), E, E - E(q1), M, omega (2), converged
single = [float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (deck(scratch / f"q{i}", False, q) / "report.out").read_text()).group(1))
          for i, q in enumerate(Q)]
dE = sweep[:, 6] - np.array(single)
for i in range(len(Q)):
    print(f"q = {Q[i]}: sweep E = {sweep[i, 6]:.10f} Ry, single-q E = {single[i]:.10f} Ry, difference {dE[i]:.3e}, converged {int(sweep[i, 11])}, "
          f"E - E(q1) = {sweep[i, 7]:.6e} Ry, M = {sweep[i, 8]:.6f} muB (tolerance {tol:.1e})")
ok = bool(np.all(np.abs(dE) < tol) and np.all(sweep[:, 11] == 1) and np.abs(sweep[1:, 7]).min() > 100 * tol)
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
