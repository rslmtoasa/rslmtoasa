#!/usr/bin/env python3
"""AppGbtKspaceScSpiralAfm: GBT SCF on two antiparallel sites (bcc Cr, CsCl cell, 12^3 mesh, first order). The collinear
reference is the nsp = 1 SCF of the deck as it stands; the GBT deck has the Cr2 potential (majority spin in channel 1) on both sites, with mom = (0,0,+1)
on site 1 and (0,0,-1) on site 2 (the frame carries the orientation, see reciprocal_gbt.f90): the mirror-image AFM state of the same energy.
Default deck tests/spin_response/cr_bcc; a developer deck directory can be given as 4th argument.
  identity  I2 (Sandratskii 1998): theta = 0, q = 0, without constraints: |dE|, |dE_F| < scf_identity_ry against the collinear AFM SCF;
            the site moments are reported. Bug: the up/down projection axis taken from potential%mom (-z on site 2) swaps the spin channels of site 2: the moments
            collapse from 0.887 to 3e-4 muB, dE = 7.6e-4 Ry.
usage: test_kspace_scspiral_afm.py <mode> <rslmto.x> <scratch directory> [deck directory]"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

from gbt_app_common import SR, tolerance

NK = 12


def deck(dst, src, gbt, theta=0.0, constrained=True, nstep=80):
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    nml = f90nml.read(str(src / "input.nml"))
    del nml["calculation"]["post_processing"]
    nml["calculation"]["pre_processing"] = "bravais"          # runs the SCF of the deck (crystal_sym = 'file')
    nml["self"]["nstep"] = nstep
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, kspace_ham_order="first")
    if gbt:
        text = (dst / "Cr2.nml").read_text()      # channel 1 = majority spin, as the local frame of every GBT site needs
        (dst / "Cr1.nml").write_text(text.replace("symbol = 'Cr2'", "symbol = 'Cr1'"))
        (dst / "Cr2.nml").write_text(re.sub(r"^ mom = .*$", " mom = 0.0, 0.0, -1.0", text, flags=re.M))
        nml["control"]["nsp"] = 3
        nml["hamiltonian"] = {"hoh": False, "magnetic_representation": "gbt_single_q", "q_ss": [0.0, 0.0, 0.0], "theta_ss": theta}
    nml.write(str(dst / "input.nml"), force=True)
    if gbt and constrained:
        with open(dst / "input.nml", "a") as f:
            f.write("\n&constraints\n constraints_enable = .true.\n constraints_i_cons = 3\n constraints_tolerance = 1e-6\n")
            f.write("".join(f" constraints_mom_ref({i},{a}) = {v}\n" for a in (1, 2) for i, v in enumerate(("0.0", "0.0", "1.0"), 1)))
            f.write("/\n")
    return dst


def scf(binary, d):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    assert r.returncode == 0, f"{d}: exit code {r.returncode}"
    return dict(E=float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1)),
                EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log,
                moments=[float(m) for _, m in re.findall(r"Spin moment of atom\s+(\d+) is\s+(\S+)", log)[-2:]])


def mode_identity(binary, scratch, src):
    tol = tolerance("scf_identity_ry")
    ref = scf(binary, deck(scratch / "collinear", src, False))
    ok = ref["conv"]
    for tag, constrained in (("free", False),):
        gbt = scf(binary, deck(scratch / tag, src, True, constrained=constrained))
        dE, dEF = gbt["E"] - ref["E"], gbt["EF"] - ref["EF"]
        print(f"{tag}: converged {gbt['conv']}; dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, tolerance {tol:.1e}; "
              f"site moments [muB]: collinear {ref['moments']}, GBT {gbt['moments']}")
        ok = ok and gbt["conv"] and abs(dE) < tol and abs(dEF) < tol
    return ok


if __name__ == "__main__":
    mode, binary, scratch = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3])
    src = Path(sys.argv[4]).resolve() if len(sys.argv) > 4 else SR / "cr_bcc"
    scratch.mkdir(parents=True, exist_ok=True)
    ok = globals()[f"mode_{mode}"](binary, scratch, src)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
