#!/usr/bin/env python3
"""AppGbtKspaceSymmetry. Relations R2, R4, R5, R6 (Sandratskii, Adv. Phys. 47, 91 (1998)) on the frozen_magnon_kspace
output for the one-site bcc Fe cell (tests/spin_response/fe_bcc potential), 8^3 Gamma-centred mesh, one cone angle.
R2: Omega(q,0) = Omega(0,0) and N(q,0) = N(0,0) for every q (q/2 on the mesh), within omega_identity_ry, n_identity.
R4: Omega(q+G,theta) = Omega(q,theta), G a reciprocal lattice vector, q/2 and (q+G)/2 on the mesh.
R5: Omega(-q,theta) = Omega(q,theta); weak in bcc (it also follows from inversion).
R6: Omega(0,theta) = Omega(0,0) holds by construction (global rotation); reported, not counted.
Bugs caught: a q/2 shift that is not a mesh translation (R2, R4), a mesh mapping that is off by one point (R2, R4),
a lost q -> -q or q -> q+G symmetry in the q handling (R4, R5). A swapped sign of the +-q/2 assignment is not seen
here (R5, inversion); UnitGbtKspaceOneBand and AppGbtKspaceGauge guard it.
usage: test_kspace_symmetry.py <rslmto.x> <scratch directory>"""
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import SR, tolerance

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
Q = [(0, 0, 0), (0, 0, 0.25), (0.25, 0, 0.5), (0, 0, -0.25), (0, 0, 1.25), (-0.25, 0, -0.5)]     # direct; G = (0, 0, 1)
PAIRS_R4, PAIRS_R5 = [(1, 4)], [(1, 3), (2, 5)]                                                  # row indices of (q, q+G), (q, -q)
THETA_DEG, NK = 10.0, 8

shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(SR / "fe_bcc", scratch)
nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
nml["calculation"].update(pre_processing="none", post_processing="frozen_magnon")
nml["lattice"].update(rc=80, ct=[3.0], r2=9.0)
for key in ("n1", "n2", "n3", "b1", "b2", "b3", "pbc"):
    nml["lattice"].pop(key, None)
nml["self"].update(nstep=3)
nml["control"].update(nsp=3)
nml["mix"].update(beta=0.01)
nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, use_symmetry_reduction=False, use_time_reversal=False, auto_find_fermi=False)
nml["hamiltonian"] = {"hoh": False, "theta_ss": THETA_DEG}
nml["frozen_magnon"] = {"mode": "mft", "q_file": "q.dat", "q_coordinates": "direct"}
nml.write(str(scratch / "input.nml"), force=True)
(scratch / "q.dat").write_text(f"{len(Q)}\n" + "\n".join(" ".join(map(str, q)) for q in Q) + "\n")
subprocess.run([str(binary), "input.nml"], cwd=scratch, check=True, stdout=open(scratch / "run.log", "w"), stderr=subprocess.STDOUT)

row = np.loadtxt(scratch / "frozen_magnon_kspace.dat")        # q direct (3), q Cartesian (3), theta, Omega(q,theta), Omega(q,0), N(q,theta), omega (2), N(q,0)
omega_t, omega_0, n_t, n_0 = row[:, 7], row[:, 8], row[:, 9], row[:, 12]
tol_omega, tol_n = tolerance("omega_identity_ry"), tolerance("n_identity")
r2_omega, r2_n = np.abs(omega_0 - omega_0[0]).max(), np.abs(n_0 - n_0[0]).max()
r4 = max(abs(omega_t[a] - omega_t[b]) for a, b in PAIRS_R4)
r5 = max(abs(omega_t[a] - omega_t[b]) for a, b in PAIRS_R5)
spread = np.abs(omega_t[1:] - omega_0[0]).min()                # the checks must have something to compare: Omega(q,theta) moves with q
print(f"R2 max |Omega(q,0) - Omega(0,0)| = {r2_omega:.3e} Ry (tol {tol_omega:.1e}); max |N(q,0) - N(0,0)| = {r2_n:.3e} (tol {tol_n:.1e})")
print(f"R4 max |Omega(q+G,theta) - Omega(q,theta)| = {r4:.3e} Ry (tol {tol_omega:.1e})")
print(f"R5 max |Omega(-q,theta) - Omega(q,theta)| = {r5:.3e} Ry (tol {tol_omega:.1e}); weak in bcc (also follows from inversion)")
print(f"R6 |Omega(0,theta) - Omega(0,0)| = {abs(omega_t[0] - omega_0[0]):.3e} Ry, by construction, not counted")
print(f"min over q != 0 of |Omega(q,theta) - Omega(0,0)| = {spread:.3e} Ry (the identities are not trivially zero)")
ok = r2_omega < tol_omega and r2_n < tol_n and r4 < tol_omega and r5 < tol_omega and spread > 100 * tol_omega
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
