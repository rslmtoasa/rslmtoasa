#!/usr/bin/env python3
"""AppGbtKspaceScSpiralI3: I3 (Sandratskii 1998), flat spirals (theta = 90 deg) without constraints. The k-space GBT SCF on the two-site CsCl cell
(fixture bccFe_sc2) equals the non-collinear periodic_nc SCF of the explicit supercell (fixtures bccFe_super_q025 / _q033: 2L sites, a3 = (0,0,L))
with every type the same potential and mom = the spiral direction, same order and temperature. q = (0,0,1/L) 2pi/alat along Gamma-H. Both sides keep
the axis of every site fixed and build the potential from the longitudinal moment (density_policy = constrained_spiral). Cr: site b antiparallel to
site a at every cell (theta_b = pi - theta, phi_b = pi). Gate: E per atom, E_F (scf_identity_ry), |m_par| and |m_perp| per site (scf_moment_mub), compared site by site.
Equivalent k-sampling: supercell Gamma-centred mesh (N1, N2, N3) and CsCl cell mesh (N1, N2, N3 L) with N3 even, q/2 on the CsCl mesh: each supercell k is
the fold of L CsCl points. A one-site primitive mesh has no equivalent (parity argument, see test_kspace_scspiral_b2.py). Precondition: the collinear SCF
of the two cells on these meshes (Fe ferromagnetic, Cr antiferromagnetic, nsp = 1) agrees in E per atom and E_F within scf_identity_ry, else the test stops.
The spiral state must keep |m_par| > 0.1 muB on every site, else the test stops (a collapsed moment makes the identity trivial).
Reported, not gated: the GBT |m_perp| per site (vanishes by symmetry without a field), the misalignment angle atan(|m_perp|/|m_par|) per site of both
sides (maximum, spread over the supercell sites), and the drift of a supercell run with density_policy = relaxed_reference (angle between each written mom
and its spiral direction).
Bug caught: a wrong q.tau phase or site frame of the GBT operator at finite q and theta = 90 deg, a wrong site-resolved density (channel assignment).
With the last argument ccor every deck carries ccor_2c (ccor_strict, Sdot, screening fitted). The supercell side builds Hcc from explicit moment
vectors, the GBT side by linking S and Sdot; the SCF densities are built from the eigenvectors of H + Hcc on both sides. Bug caught there:
the link or the onsite term of Hcc missing in the SCF assembly of the spiral at finite q and theta = 90 deg.
usage: test_kspace_scspiral_i3.py <fe|cr> <3|4> <first|second> <rslmto.x> <scratch directory> [N1 N2 N3] [ccor]"""
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import FIX, SR, tolerance

SUPER = {3: "bccFe_super_q033", 4: "bccFe_super_q025"}
CCOR = False         # set by the optional last argument ccor: every deck then carries ccor_2c (strict, Sdot, fitted screening)
MIN_MOMENT = 0.1      # muB: non-triviality guard (not an identity tolerance); a collapsed moment would make every comparison hold trivially


def potential(system, kind, i, site):
    """Text of the potential file of type i (1-based) at CsCl sublattice site (0 or 1)."""
    if system == "fe":
        text, symbol = (SR / "fe_bcc" / "Fe.nml").read_text(), "Fe"
    else:     # Cr: the collinear deck stores channels by lab spin (Cr1, Cr2); GBT and non-collinear need channel 1 = majority (Cr2)
        text, symbol = (SR / "cr_bcc" / ("Cr2.nml" if kind != "col" or site == 1 else "Cr1.nml")).read_text(), "Cr"
        text = text.replace("symbol = 'Cr2'", f"symbol = '{symbol}'").replace("symbol = 'Cr1'", f"symbol = '{symbol}'")
        return text.replace(f"symbol = '{symbol}'", f"symbol = '{symbol}{i}'")
    return text.replace(f"symbol = '{symbol}'", f"symbol = '{symbol}{i}'")


def deck(dst, system, L, cell, kind, order, nk, relaxed=False, nstep=300):
    """kind: col (nsp = 1), gbt (sc2 cell only) or nc (supercell only, periodic_nc with the spiral moments)."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    n = 2 if cell == "sc2" else 2 * L
    shutil.copy(FIX / ("bccFe_sc2" if cell == "sc2" else SUPER[L]) / "lattice.nml", dst / "lattice.nml")
    q_z = 1.0 / L
    for i in range(1, n + 1):
        site = (i - 1) % 2
        text = potential(system, kind, i, site)
        if kind == "gbt":
            mom = (0.0, 0.0, -1.0 if system == "cr" and site == 1 else 1.0)
        elif kind == "nc":
            phi = 2.0 * math.pi * q_z * 0.5 * (i - 1)
            sign = -1.0 if system == "cr" and site == 1 else 1.0
            mom = (sign * math.cos(phi), sign * math.sin(phi), 0.0)
        else:
            mom = None
        if mom is not None:
            text = re.sub(r"^ mom = .*$", " mom = " + ", ".join(repr(float(x)) for x in mom), text, flags=re.M)
        (dst / f"{'Fe' if system == 'fe' else 'Cr'}{i}.nml").write_text(text)
    nml = f90nml.read(str(SR / ("fe_bcc" if system == "fe" else "cr_bcc") / "input.nml"))
    nml["calculation"] = {"pre_processing": "bravais", "verbose": True}
    n3 = -(-20 // L) if cell == "super" else 20
    nml["lattice"].update(crystal_sym="file", ntype=n, ct=[4.0] * n, n3=n3, ndim=20 * 20 * n3 * n + 1000)
    nml["atoms"] = {"database": "./", "label": [f"{'Fe' if system == 'fe' else 'Cr'}{i}" for i in range(1, n + 1)]}
    nml["self"].update(nstep=nstep)
    nml["mix"].update(beta=0.02, mixtype="broyden")      # the default beta = 0.1 (Fe) / 0.15 (Cr) diverges on the 2L-site supercells after reaching 1e-5
    nml["control"].update(nsp=1 if kind == "col" else 3)
    if kind != "col":
        nml["control"].update(gbt_scf_diagnostics=True, density_policy="relaxed_reference" if relaxed else "constrained_spiral")
    nml["reciprocal"].update(nk1=nk[0], nk2=nk[1], nk3=nk[2], kspace_ham_order=order)
    nml["hamiltonian"] = {"hoh": order == "second"}
    if CCOR:
        nml["lattice"].update(strux_want_sdot=True, screening="fitted")
        nml["hamiltonian"].update(ccor_2c=True, ccor_strict=True)
    if kind == "gbt":
        nml["hamiltonian"].update(magnetic_representation="gbt_single_q", q_ss=[0.0, 0.0, q_z], theta_ss=90.0)
    nml.write(str(dst / "input.nml"), force=True)
    return dst


def scf(binary, d, n):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    assert r.returncode == 0, f"{d}: exit code {r.returncode}"
    out = dict(E=float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1)) / n,
               EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log)
    if (d / "gbt_scf_diagnostics.dat").exists():
        rows = [l.split() for l in (d / "gbt_scf_diagnostics.dat").read_text().splitlines() if not l.startswith("#")][-n:]
        out["par"] = np.array([abs(float(r[4 + 23])) for r in rows])                 # |longitudinal moment| per site, last iteration
        out["perp"] = np.array([float(r[4 + 24]) for r in rows])                     # transverse moment per site
    return out


def angles(r):
    return np.degrees(np.arctan2(r["perp"], r["par"]))


def main(binary, scratch, system, L, order, nk_super):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    nk_cell = (nk_super[0], nk_super[1], nk_super[2] * L)
    n = 2 * L
    pre = [scf(binary, deck(scratch / "pre_cell", system, L, "sc2", "col", order, nk_cell), 2),
           scf(binary, deck(scratch / "pre_super", system, L, "super", "col", order, nk_super), n)]
    d_pre = (abs(pre[0]["E"] - pre[1]["E"]), abs(pre[0]["EF"] - pre[1]["EF"]))
    print(f"precondition (collinear, {order} order): |dE/atom| = {d_pre[0]:.3e} Ry, |dE_F| = {d_pre[1]:.3e} Ry, tolerance {tol:.1e}")
    if max(d_pre) >= tol:
        sys.exit("precondition failed: the two meshes do not sample the same k set; stopping")
    gbt = scf(binary, deck(scratch / "gbt", system, L, "sc2", "gbt", order, nk_cell), 2)
    nc = scf(binary, deck(scratch / "nc", system, L, "super", "nc", order, nk_super), n)
    if min(gbt["par"].min(), nc["par"].min()) < MIN_MOMENT:
        sys.exit(f"the spiral state is (nearly) non-magnetic (|m_par| = {gbt['par']} muB); the identity would be trivial; stopping")
    dE, dEF = abs(gbt["E"] - nc["E"]), abs(gbt["EF"] - nc["EF"])
    site = np.arange(n) % 2           # CsCl sublattice of each supercell site
    dpar = np.abs(nc["par"] - gbt["par"][site]).max()
    dperp = np.abs(nc["perp"] - gbt["perp"][site]).max()
    print(f"{system} q = 1/{L}, {order} order, meshes CsCl {nk_cell} / supercell {tuple(nk_super)}: converged {gbt['conv']}, {nc['conv']}; "
          f"|dE/atom| = {dE:.3e} Ry, |dE_F| = {dEF:.3e} Ry, max |d m_par| = {dpar:.3e} muB, max |d m_perp| = {dperp:.3e} muB (tolerances {tol:.1e}, {tol_m:.1e})")
    print(f"GBT |m_par| per site {np.round(gbt['par'], 6)} muB, |m_perp| per site {gbt['perp']} muB")
    a_g, a_n = angles(gbt), angles(nc)
    print(f"misalignment angle [deg]: GBT {a_g}, supercell max {a_n.max():.3e}, spread over sites {a_n.max() - a_n.min():.3e}")
    d = scratch / "nc_relaxed"
    rel = deck(d, system, L, "super", "nc", order, nk_super, relaxed=True, nstep=60)
    scf(binary, rel, n)
    drift = []
    for i in range(1, n + 1):
        written = np.array([float(x) for x in re.search(r"^ mom = (.*)$", (rel / f"{'Fe' if system == 'fe' else 'Cr'}{i}_scf.nml").read_text(), flags=re.M).group(1).split(",")])
        start = np.array([float(x) for x in re.search(r"^ mom = (.*)$", (rel / f"{'Fe' if system == 'fe' else 'Cr'}{i}.nml").read_text(), flags=re.M).group(1).split(",")])
        drift.append(math.degrees(math.acos(max(-1.0, min(1.0, float(written @ start) / np.linalg.norm(written) / np.linalg.norm(start))))))
    print(f"relaxed_reference supercell (report only): drift of the site directions [deg]: max {max(drift):.3e}, spread {max(drift) - min(drift):.3e}")
    return gbt["conv"] and nc["conv"] and dE < tol and dEF < tol and dpar < tol_m and dperp < tol_m


if __name__ == "__main__":
    CCOR = sys.argv[-1] == "ccor"
    if CCOR:
        sys.argv.pop()
    system, L, order = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    binary, scratch = Path(sys.argv[4]).resolve(), Path(sys.argv[5])
    nk_super = tuple(int(x) for x in sys.argv[6:9]) if len(sys.argv) > 8 else (4, 4, 2)
    scratch.mkdir(parents=True, exist_ok=True)
    ok = main(binary, scratch, system, L, order, nk_super)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
