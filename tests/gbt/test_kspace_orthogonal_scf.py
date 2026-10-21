#!/usr/bin/env python3
"""AppKspaceOrthogonalScf: self-consistent bcc Fe (tests/spin_response/fe_bcc potential, 12^3 mesh, 300 K) with kspace_ham_order = 'orthogonal'
outside GBT (periodic_nc nsp = 3 and collinear nsp = 1). The exact-ASA H^gamma itself is anchored by oracle O3 (eigenvalues) through
AppGbtKspaceOrthogonal; the moments by O3b. The checks here are of the SCF wiring, and the first two modes run the same kernel on both sides:
  rotation   global spin rotation (no SOC): the SCF with every moment rotated to (theta, phi) gives the same E and E_F (scf_identity_ry), the
             same |m| and the moment direction of the start, as |m| * (final direction - start direction) (scf_moment_mub), the direction read from the written potential Fe_scf.nml at full precision. Eigenvalues are invariant under the frame
             rotation anyway, so what this can catch is the eigenvectors: a local-frame instead of lab spin basis, or a frame built from the
             wrong moment components. Bug caught there: dropping the local-to-lab rotation of H (the moment then comes out along the wrong axis).
  collinear  consistency, not validation: the nsp = 1 SCF equals the nsp = 3 SCF with mom = +z (E, E_F, |m|). Bug caught: a spin channel
             of the collinear path taken from the moment sign instead of the lab order.
  afm        consistency, not validation: collinear AFM bcc Cr (tests/spin_response/cr_bcc, CsCl cell, 12^3 mesh): the nsp = 1 SCF equals the nsp = 3 SCF (periodic_nc,
             mom = +z / -z per site, the mirror-image state with channel 1 = majority on both sites) in E, E_F (scf_identity_ry) and the two |m| (scf_moment_mub).
             Bug caught: the nsp = 3 frame ignoring the moment sign (site 2 not flipped: the nsp = 3 moments collapse to 0, dE = 1.2e-3 Ry). An nsp = 1 frame
             taken from the moment sign is not reachable here: the nsp = 1 SCF leaves mom = +z on both sites, and mom = -z in an nsp = 1 file collapses the moment.
  moments    identity against oracle O3b (tests/gbt/oracles/o3_hgamma, H^gamma moments of bcc Fe): the nsp = 1 SCF of the fixed fe_bcc potential at 300 K on the full
             Gamma-centred 8^3 and 16^3 meshes (no symmetry reduction, weights 1/N^3), read after one iteration with linear mixing beta = 1, which returns the raw
             moments: ql(1) = n0, gravity_center + vmad (input file) = centre, ql(3) = central q2, per l and spin (1 = majority). Bug caught: a moment taken from
             another eigensystem than the H^gamma one (e.g. the first-order H), a weight or Fermi-level error, a wrong l or spin assignment.
  guards     'orthogonal' stops with its message for ccor_2c, for hoh = .true. and for constraining fields. Bug caught: a guard removed, the run
             then solves a Hamiltonian that is not the one asked for (CCOR dropped silently, h - h o h on top of H^gamma).
usage: test_kspace_orthogonal_scf.py <mode> <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import ROOT, SR, tolerance

NK = 12


def deck(dst, nsp, direction=None, nstep=80):
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "fe_bcc", dst)
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["self"]["nstep"] = nstep
    nml["control"]["nsp"] = nsp
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, kspace_ham_order="orthogonal")
    nml["hamiltonian"] = {"hoh": False}
    nml.write(str(dst / "input.nml"), force=True)
    if direction is not None:
        f = dst / "Fe.nml"
        f.write_text(re.sub(r"^ mom = .*$", " mom = " + ", ".join(repr(float(x)) for x in direction), f.read_text(), flags=re.M))
    return dst


def scf(binary, d):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    assert r.returncode == 0, f"{d}: exit code {r.returncode}"
    return dict(E=float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1)),
                EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log,
                m=float(re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", log)[-1]),
                moments=np.array([float(m) for _, m in re.findall(r"Spin moment of atom\s+(\d+) is\s+(\S+)", log)[-2:]]))


def final_direction(d):
    """Unit moment direction of atom 1 from the potential written by the run (Fe_scf.nml, full precision; the log prints f10.6)."""
    return np.array([float(x) for x in re.search(r"^\s*mom\s*=\s*(.*)$", (d / "Fe_scf.nml").read_text(), re.M).group(1).split(",")])


def mode_rotation(binary, scratch):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    ref = scf(binary, deck(scratch / "rot_0_0", 3, (0.0, 0.0, 1.0)))
    print(f"unrotated: converged {ref['conv']}, E = {ref['E']:.10f} Ry, E_F = {ref['EF']:.10f} Ry, |m| = {ref['m']:.6f} muB")
    ok = ref["conv"]
    for theta, phi in ((30, 0), (45, 45), (60, 120), (180, 0)):
        t, p = np.radians(theta), np.radians(phi)
        start = np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
        d = deck(scratch / f"rot_{theta}_{phi}", 3, start)
        r = scf(binary, d)
        dE, dEF, dm = r["E"] - ref["E"], r["EF"] - ref["EF"], abs(r["m"] - ref["m"])
        dvec = r["m"] * np.abs(final_direction(d) - start).max()
        print(f"theta {theta}, phi {phi}: converged {r['conv']}; dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, d|m| = {dm:.3e} muB, "
              f"|m| max |dir - start| = {dvec:.3e} muB (tolerances {tol:.1e}, {tol_m:.1e})")
        ok = ok and r["conv"] and abs(dE) < tol and abs(dEF) < tol and dm < tol_m and dvec < tol_m
    return ok


def mode_collinear(binary, scratch):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    col = scf(binary, deck(scratch / "nsp1", 1))
    nc = scf(binary, deck(scratch / "nsp3", 3, (0.0, 0.0, 1.0)))
    dE, dEF, dm = nc["E"] - col["E"], nc["EF"] - col["EF"], abs(nc["m"] - col["m"])
    print(f"nsp = 1: E = {col['E']:.10f} Ry, E_F = {col['EF']:.10f} Ry, |m| = {col['m']:.6f} muB; nsp = 3 minus nsp = 1: "
          f"dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, d|m| = {dm:.3e} muB (tolerances {tol:.1e}, {tol_m:.1e})")
    return col["conv"] and nc["conv"] and abs(dE) < tol and abs(dEF) < tol and dm < tol_m


def afm_deck(dst, nsp):
    """cr_bcc as the deck of test_kspace_scspiral_afm. nsp = 1: the files as they stand (Cr1 holds the minority values in channel 1, Cr2 the majority; mom = +z in both,
    the SCF then sets mom from the sign of the site moment). nsp = 3: the Cr2 potential (channel 1 = majority) on both sites with mom = +z, -z."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "cr_bcc", dst)
    nml = f90nml.read(str(SR / "cr_bcc" / "input.nml"))
    del nml["calculation"]["post_processing"]
    nml["calculation"]["pre_processing"] = "bravais"
    nml["self"]["nstep"] = 80
    nml["control"]["nsp"] = nsp
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, kspace_ham_order="orthogonal")
    nml["hamiltonian"] = {"hoh": False}
    nml.write(str(dst / "input.nml"), force=True)
    cr2 = (dst / "Cr2.nml").read_text()
    setmom = lambda text, z: re.sub(r"^ mom = .*$", f" mom = 0.0, 0.0, {z}", text, flags=re.M)
    if nsp == 3:
        (dst / "Cr1.nml").write_text(setmom(cr2.replace("symbol = 'Cr2'", "symbol = 'Cr1'"), 1.0))
        (dst / "Cr2.nml").write_text(setmom(cr2, -1.0))
    return dst


def mode_afm(binary, scratch):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    col, nc = (scf(binary, afm_deck(scratch / f"nsp{n}", n)) for n in (1, 3))
    dE, dEF, dm = nc["E"] - col["E"], nc["EF"] - col["EF"], np.abs(nc["moments"] - col["moments"]).max()
    print(f"nsp = 1: converged {col['conv']}, E = {col['E']:.10f} Ry, E_F = {col['EF']:.10f} Ry, |m| = {col['moments']} muB; nsp = 3 minus nsp = 1: "
          f"dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, max d|m| = {dm:.3e} muB (nsp = 3 |m| = {nc['moments']}; tolerances {tol:.1e}, {tol_m:.1e})")
    if col["moments"].min() < 0.1:
        print("the collinear AFM moment collapsed (below 0.1 muB); the identity would be trivial")
        return False
    return col["conv"] and nc["conv"] and abs(dE) < tol and abs(dEF) < tol and dm < tol_m


def par_array(text, key):
    """(l, spin) values of a &par line such as 'ql(1, :, 1) = a, b, c' (key 'ql(1, :, {s})') for spin 1, 2."""
    return np.array([[float(x) for x in re.search(rf"^\s*{re.escape(key.format(s=s))}\s*=\s*(.*)$", text, re.M).group(1).split(",")] for s in (1, 2)]).T


def mode_moments(binary, scratch):
    tol = tolerance("n_identity")
    oracle = f90nml.read(str(ROOT / "tests/gbt/oracles/o3_hgamma/o3_fe_bcc.nml"))
    vmad = float(re.search(r"^\s*vmad\s*=\s*(\S+)", (SR / "fe_bcc" / "Fe.nml").read_text(), re.M).group(1).replace("D", "E"))
    ok = True
    for n in (8, 16):
        d = deck(scratch / f"mesh{n}", 1, nstep=1)
        nml = f90nml.read(str(d / "input.nml"))
        nml["reciprocal"].update(nk1=n, nk2=n, nk3=n, temperature=300.0)
        nml["mix"] = {"beta": 1.0, "mixtype": "linear"}
        nml.write(str(d / "input.nml"), force=True)
        r = scf(binary, d)
        text = (d / "Fe_scf.nml").read_text()
        ref = oracle[f"o3b_moments_mesh{n}"]
        pairs = {"n0": (par_array(text, "ql(1, :, {s})"), ref["n0"]),
                 "centre": (par_array(text, "gravity_center(:, {s})") + vmad, ref["centre_ry"]),
                 "central q2": (par_array(text, "ql(3, :, {s})"), ref["central_q2_ry2"])}
        dev = {k: np.abs(v[0] - np.array(v[1]).reshape((3, 2), order="F")).max() for k, v in pairs.items()}
        print(f"{n}^3 (one iteration, mixing beta = 1): E_F code {r['EF']:.10f} Ry, O3b mu {ref['mu_ry']:.10f} Ry; max |code - O3b| over l, spin: "
              + ", ".join(f"{k} {v:.3e}" for k, v in dev.items()) + f" (tolerance {tol:.1e})")
        ok = ok and all(v < tol for v in dev.values())
    return ok


def mode_guards(binary, scratch):
    cases = (("ccor", "does not support ccor_2c"), ("hoh", "hoh must be .false."), ("constraints", "does not support constraining fields"))
    ok = True
    for tag, message in cases:
        d = deck(scratch / tag, 3, (0.0, 0.0, 1.0), nstep=1)
        nml = f90nml.read(str(d / "input.nml"))
        if tag == "ccor":
            nml["lattice"].update(strux_want_sdot=True, screening="fitted")
            nml["hamiltonian"].update(ccor_2c=True, ccor_strict=True)
        if tag == "hoh":
            nml["hamiltonian"]["hoh"] = True
        nml.write(str(d / "input.nml"), force=True)
        if tag == "constraints":
            with open(d / "input.nml", "a") as f:
                f.write("\n&constraints\n constraints_enable = .true.\n constraints_i_cons = 3\n constraints_tolerance = 1e-6\n")
                f.write("".join(f" constraints_mom_ref({i},1) = {v}\n" for i, v in enumerate(("0.0", "0.0", "1.0"), 1)) + "/\n")
        r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
        found = message in re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
        print(f"{tag}: exit code {r.returncode}, message found: {found}")
        ok = ok and r.returncode != 0 and found
    return ok


if __name__ == "__main__":
    mode, binary, scratch = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3])
    scratch.mkdir(parents=True, exist_ok=True)
    ok = globals()[f"mode_{mode}"](binary, scratch)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
