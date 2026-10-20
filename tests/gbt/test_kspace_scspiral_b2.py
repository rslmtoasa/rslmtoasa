#!/usr/bin/env python3
"""AppGbtKspaceScSpiralB2: I1 (equivalent sites). bcc Fe as the primitive cell (16^3 mesh) and as the doubled primitive cell with two identical
sites (a1, a2, 2*a3, the second atom at a3; fixture bccFe_dbl, mesh 16x16x8): the constrained GBT SCF gives the same energy per atom, E_F, moment
per site and constraining field per site, at theta = 0 and 30 deg, q = (0, 0.25, 0) 2pi/alat (q/2 on both meshes; q.tau = 0.125 is a nontrivial
phase between the two sites). The two Gamma-centred meshes sample the same k set (the primitive k and k + a3-direction half step fold onto one
doubled-cell k), which is the precondition: the collinear SCF of the two cells must agree in E per atom and E_F within scf_identity_ry, else
the test stops. The CsCl cell (bccFe_sc2) has no equivalent mesh pair: a primitive n^3 mesh is the even-parity half of the cubic n^3 mesh.
Deck: tests/spin_response/fe_bcc potential, first order, 300 K. Tolerances: scf_identity_ry (E per atom, E_F, B), scf_moment_mub (M).
Bug caught: a missing or wrong q.tau phase between the sites, a frame mix-up between the sites, or the constraining field indexed by type instead
of by site (a site left without its field: its B and m_perp differ from the primitive cell).
usage: test_kspace_scspiral_b2.py <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import FIX, SR, tolerance

Q = (0.0, 0.25, 0.0)
CELLS = {"prim": ("bccFe_prim", (16, 16, 16), 1), "dbl": ("bccFe_dbl", (16, 16, 8), 2)}


def deck(dst, cell, nsp, theta=0.0, constrained=True, nstep=80):
    fixture, nk, n = CELLS[cell]
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    shutil.copy(FIX / fixture / "lattice.nml", dst / "lattice.nml")
    text = (SR / "fe_bcc" / "Fe.nml").read_text()
    for i in range(1, n + 1):
        (dst / f"Fe{i}.nml").write_text(text.replace("symbol = 'Fe'", f"symbol = 'Fe{i}'"))
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["lattice"].update(crystal_sym="file", ntype=n, ct=[4.0] * n)
    nml["atoms"] = {"database": "./", "label": [f"Fe{i}" for i in range(1, n + 1)]}
    nml["self"]["nstep"] = nstep
    nml["control"]["nsp"] = nsp
    nml["reciprocal"].update(nk1=nk[0], nk2=nk[1], nk3=nk[2], kspace_ham_order="first")
    nml["hamiltonian"] = {"hoh": False}
    if nsp == 3:
        nml["hamiltonian"].update(magnetic_representation="gbt_single_q", q_ss=list(Q), theta_ss=theta)
    nml.write(str(dst / "input.nml"), force=True)
    if nsp == 3 and constrained:
        with open(dst / "input.nml", "a") as f:
            f.write("\n&constraints\n constraints_enable = .true.\n constraints_i_cons = 3\n constraints_tolerance = 1e-6\n constraints_diagnostics = .true.\n")
            f.write("".join(f" constraints_mom_ref({i},{a}) = {v}\n" for a in range(1, n + 1) for i, v in enumerate(("0.0", "0.0", "1.0"), 1)))
            f.write("/\n")
    return dst


def scf(binary, d, n):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    assert r.returncode == 0, f"{d}: exit code {r.returncode}"
    out = dict(E=float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1)) / n,
               EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log)
    if (d / "constraint_diagnostics.dat").exists():
        rows = [[float(x) for x in l.split()] for l in (d / "constraint_diagnostics.dat").read_text().splitlines() if not l.startswith("#")]
        out["M"] = np.array([r[8] for r in rows[-n:]])                       # moment magnitude per site, last iteration
        out["B"] = np.array([float(x) for x in re.findall(r"B_next_xy_per_site \(Ry\) =(.*)", log)[-1].split()]).reshape(n, 2)
    return out


def main(binary, scratch):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    prim = scf(binary, deck(scratch / "pre_prim", "prim", 1), 1)
    dbl = scf(binary, deck(scratch / "pre_dbl", "dbl", 1), 2)
    pre = (abs(dbl["E"] - prim["E"]), abs(dbl["EF"] - prim["EF"]))
    print(f"precondition (collinear): |dE/atom| = {pre[0]:.3e} Ry, |dE_F| = {pre[1]:.3e} Ry, tolerance {tol:.1e}")
    if max(pre) >= tol:
        sys.exit("precondition failed: the two meshes do not sample the same k set; stopping")
    ok = True
    for theta in (0.0, 30.0):
        p = scf(binary, deck(scratch / f"prim_{theta:g}", "prim", 3, theta), 1)
        d = scf(binary, deck(scratch / f"dbl_{theta:g}", "dbl", 3, theta), 2)
        dev = dict(E=abs(d["E"] - p["E"]), EF=abs(d["EF"] - p["EF"]), M=np.abs(d["M"] - p["M"][0]).max(), B=np.abs(d["B"] - p["B"][0]).max())
        print(f"theta = {theta:g}: converged {p['conv']}, {d['conv']}; |dE/atom| {dev['E']:.3e} Ry, |dE_F| {dev['EF']:.3e} Ry, max |dM| {dev['M']:.3e} muB, "
              f"max |dB| {dev['B']:.3e} Ry; M = {p['M'][0]:.6f} muB, B_x = {p['B'][0][0]:.4e} Ry")
        ok = ok and p["conv"] and d["conv"] and dev["E"] < tol and dev["EF"] < tol and dev["B"] < tol and dev["M"] < tol_m
    return ok


if __name__ == "__main__":
    binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    scratch.mkdir(parents=True, exist_ok=True)
    ok = main(binary, scratch)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
