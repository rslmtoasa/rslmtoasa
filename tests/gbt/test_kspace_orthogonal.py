#!/usr/bin/env python3
"""AppGbtKspaceOrthogonal. Relations R8 and R11 (Andersen & Jepsen, PRL 53, 2571 (1984); Mermin 1965) on the exact-ASA H^gamma:
frozen_magnon_kspace with kspace_ham_order = 'orthogonal' on tests/spin_response/fe_bcc, zero SCF steps, 32^3 Gamma-centred mesh, 300 K,
theta_ss = 5 deg, q = (0.125,-0.125,0.125) and (0.25,-0.25,0.25) direct. Omega(0,0), Omega(q,0), Omega(q,5 deg), N and mu are compared with
oracle O3 (tests/gbt/oracles/o3_hgamma, group o3_driver_mesh32, direct sums, no Richardson) at omega_identity_ry and n_identity.
M is reported only. Bugs caught: wrong wow exponent or alpha source in the glue, C without vmad, dele used as Delta, Gamma of the wrong sign
or side, a hidden SCF step, a k -+ q/2 swap, a wrong mu (occupation or smearing).
usage: test_kspace_orthogonal.py <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import ROOT, SR, tolerance

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
Q = [(0.125, -0.125, 0.125), (0.25, -0.25, 0.25)]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(SR / "fe_bcc", scratch)
nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
nml["calculation"].update(post_processing="frozen_magnon")
nml["self"].update(nstep=0)
nml["control"].update(nsp=3)
nml["reciprocal"].update(nk1=32, nk2=32, nk3=32, temperature=300.0, auto_find_fermi=True, kspace_ham_order="orthogonal")
nml["hamiltonian"] = {"hoh": False, "theta_ss": 5.0}
nml["frozen_magnon"] = {"mode": "mft", "q_file": "q.dat", "q_coordinates": "direct"}
nml.write(str(scratch / "input.nml"), force=True)
(scratch / "q.dat").write_text("3\n0 0 0\n" + "\n".join(" ".join(map(str, q)) for q in Q) + "\n")
subprocess.run([str(binary), "input.nml"], cwd=scratch, check=True, stdout=open(scratch / "run.log", "w"), stderr=subprocess.STDOUT)

group = re.search(r"&o3_driver_mesh32(.*?)\n/", (ROOT / "tests/gbt/oracles/o3_hgamma/o3_fe_bcc.nml").read_text(), re.S).group(1)
o3 = {k: float(v) for k, v in re.findall(r"^\s*(\w+)\s*=\s*([-+0-9.eE]+)\s*$", group, re.M)}
header = (scratch / "frozen_magnon_kspace.dat").read_text().splitlines()[0]
mu, moment = (float(re.search(rf"{k}\s*=\s*([-+0-9.eE]+)", header).group(1)) for k in (r"mu \[Ry\]", r"M \[muB\]"))
row = np.loadtxt(scratch / "frozen_magnon_kspace.dat")      # q direct (3), q Cartesian (3), theta, Omega(q,theta), Omega(q,0), N(q,theta), omega (2), N(q,0)
tol_omega, tol_n = tolerance("omega_identity_ry"), tolerance("n_identity")
dev = {"mu": abs(mu - o3["mu_ry"]), "Omega(0,0)": abs(row[0, 8] - o3["omega_00_ry"]), "N(0,0)": abs(row[0, 12] - o3["n_00"])}
for iq in (1, 2):
    dev[f"Omega(q{iq},5deg)"] = abs(row[iq, 7] - o3[f"omega_q{iq}_theta5_ry"])
    dev[f"Omega(q{iq},0)"] = abs(row[iq, 8] - o3[f"omega_q{iq}_theta0_ry"])
    dev[f"N(q{iq},5deg)"] = abs(row[iq, 9] - o3[f"n_q{iq}_theta5"])
for name, d in dev.items():
    print(f"{name:16s} |code - O3| = {d:.3e}  (tol {tol_n if name.startswith('N') else tol_omega:.1e})")
print(f"mu = {mu:.9f} Ry (O3 {o3['mu_ry']:.9f}); M = {moment:.6f} muB (reported only)")
ok = all(d < (tol_n if name.startswith("N") else tol_omega) for name, d in dev.items())
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
