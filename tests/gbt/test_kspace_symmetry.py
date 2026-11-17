#!/usr/bin/env python3
"""AppGbtKspaceSymmetry. Relations R2, R4, R5, R6 (Sandratskii, Adv. Phys. 47, 91 (1998)) on the frozen_magnon_kspace
output for the one-site bcc Fe cell (tests/spin_response/fe_bcc potential), 8^3 Gamma-centred mesh, one cone angle.
R2: Omega(q,0) = Omega(0,0) and N(q,0) = N(0,0) for every q (q/2 on the mesh), within omega_identity_ry, n_identity.
R4: Omega(q+G,theta) = Omega(q,theta), G a reciprocal lattice vector, q/2 and (q+G)/2 on the mesh.
R5: Omega(-q,theta) = Omega(q,theta); weak in bcc (it also follows from inversion).
R6: Omega(0,theta) = Omega(0,0) holds by construction (global rotation); reported, not counted.
Zero SCF steps (&self nstep = 0): every potential parameter the run writes (Fe_scf.nml), or none is written, equals the input Fe.nml bit for bit.
Bugs caught: a q/2 shift that is not a mesh translation (R2, R4), a mesh mapping that is off by one point (R2, R4),
a lost q -> -q or q -> q+G symmetry in the q handling (R4, R5), a hidden SCF step that changes the reference potential
(zero-step check), the same identities for kspace_ham_order = 'orthogonal'. A swapped sign of the +-q/2 assignment is not seen
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

def make_deck(deck, nstep, order="auto"):
    shutil.rmtree(deck, ignore_errors=True)
    shutil.copytree(SR / "fe_bcc", deck)
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["calculation"].update(pre_processing="none", post_processing="frozen_magnon")
    nml["lattice"].update(rc=80, ct=[3.0], r2=9.0)
    for key in ("n1", "n2", "n3", "b1", "b2", "b3", "pbc"):
        nml["lattice"].pop(key, None)
    nml["self"].update(nstep=nstep)
    nml["control"].update(nsp=3)
    nml["mix"].update(beta=0.01)
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, use_symmetry_reduction=False, use_time_reversal=False, auto_find_fermi=False,
                             kspace_ham_order=order)
    nml["hamiltonian"] = {"hoh": False, "theta_ss": THETA_DEG}
    nml["frozen_magnon"] = {"mode": "mft", "q_file": "q.dat", "q_coordinates": "direct"}
    nml.write(str(deck / "input.nml"), force=True)
    (deck / "q.dat").write_text(f"{len(Q)}\n" + "\n".join(" ".join(map(str, q)) for q in Q) + "\n")
    subprocess.run([str(binary), "input.nml"], cwd=deck, check=True, stdout=open(deck / "run.log", "w"), stderr=subprocess.STDOUT)


def identities(deck, label):
    row = np.loadtxt(deck / "frozen_magnon_kspace.dat")        # q direct (3), q Cartesian (3), theta, Omega(q,theta), Omega(q,0), N(q,theta), omega (2), N(q,0)
    omega_t, omega_0, n_0 = row[:, 7], row[:, 8], row[:, 12]
    tol_omega, tol_n = tolerance("omega_identity_ry"), tolerance("n_identity")
    r2_omega, r2_n = np.abs(omega_0 - omega_0[0]).max(), np.abs(n_0 - n_0[0]).max()
    r4 = max(abs(omega_t[a] - omega_t[b]) for a, b in PAIRS_R4)
    r5 = max(abs(omega_t[a] - omega_t[b]) for a, b in PAIRS_R5)
    spread = np.abs(omega_t[1:] - omega_0[0]).min()            # the checks must have something to compare: Omega(q,theta) moves with q
    print(f"[{label}] R2 max |Omega(q,0) - Omega(0,0)| = {r2_omega:.3e} Ry (tol {tol_omega:.1e}); max |N(q,0) - N(0,0)| = {r2_n:.3e} (tol {tol_n:.1e})")
    print(f"[{label}] R4 max |Omega(q+G,theta) - Omega(q,theta)| = {r4:.3e} Ry (tol {tol_omega:.1e})")
    print(f"[{label}] R5 max |Omega(-q,theta) - Omega(q,theta)| = {r5:.3e} Ry (tol {tol_omega:.1e}); weak in bcc (also follows from inversion)")
    print(f"[{label}] R6 |Omega(0,theta) - Omega(0,0)| = {abs(omega_t[0] - omega_0[0]):.3e} Ry, by construction, not counted")
    print(f"[{label}] min over q != 0 of |Omega(q,theta) - Omega(0,0)| = {spread:.3e} Ry (the identities are not trivially zero)")
    return r2_omega < tol_omega and r2_n < tol_n and r4 < tol_omega and r5 < tol_omega and spread > 100 * tol_omega


make_deck(scratch, 3)
ok = identities(scratch, "first order, nstep 3")
zero = scratch.parent / (scratch.name + "_nstep0")
make_deck(zero, 0)
out, given = zero / "Fe_scf.nml", f90nml.read(str(zero / "Fe.nml"))
written = f90nml.read(str(out)) if out.exists() else {}
changed = [f"{g}.{k}" for g in written for k in written[g] if np.any(np.asarray(written[g][k]) != np.asarray(given[g][k]))]
print(f"zero SCF steps: {sum(len(written[g]) for g in written)} potential parameters written, {len(changed)} differ from the input: {changed[:4]}")
ok = ok and not changed
ortho = scratch.parent / (scratch.name + "_orthogonal")        # exact-ASA H^gamma, fixed potential (kspace_ham_order = 'orthogonal')
make_deck(ortho, 0, "orthogonal")
ok = identities(ortho, "orthogonal, nstep 0") and ok
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
