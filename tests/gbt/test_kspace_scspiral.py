#!/usr/bin/env python3
"""AppGbtKspaceScSpiral: one-site constrained-spiral SCF on the k-space GBT path (bcc Fe, tests/spin_response/fe_bcc potential,
16^3 mesh, 300 K, first order). Modes, each with its relation and the bug it catches:
  identity  I2 (spiral SCF at theta = 0 equals the collinear SCF): |dE|, |dE_F| < scf_identity_ry. Bug: the field leaking into H at
            B = 0, a wrong k -+ q/2 or a swapped spin channel in the SCF density.
  identity_orthogonal  the same at kspace_ham_order = 'orthogonal'. Consistency, not validation: both sides run the same H^gamma kernel family
            (the physics anchors are O3 and O3b). Bug caught: a k -+ q/2 shift that is not a mesh translation (q
            scaled by 0.3 in the assembly: dE = 1.1e-4 Ry). Not caught: errors shared with the nsp = 1 path (a swapped spin channel of C moves both sides; that one is caught by
            AppKspaceOrthogonalScf_moments (O3b), not by I2), or
            a shift by a multiple of q/2 (a mesh translation).
  identity_ccor, identity_ccor_second  the same with ccor_2c on both sides (first and second order). Bug: Hcc missing or wrong in the spiral SCF
            assembly (also in the SCF density: the moments then differ), Hcc inside h - h o h at second order.
  converge  S3, S8: theta = 30 deg, q != 0: the SCF converges, max |m_perp| falls below scf_moment_mub and sum_a |B_a||m_a_perp| below scf_identity_ry, N stays
            neutral every iteration (|dN| < n_identity). Bug: a field that does not enter H (m_perp stays 3e-4 muB).
  converge_orthogonal  the same at kspace_ham_order = 'orthogonal'. Bug: the field not entering H there (m_perp stays at its unconstrained value).
  chi       S3 (Pauli response): with H += +B.sigma, m = n_up - n_dn is pushed against B, so chi = dm_perp/dB at fixed potential is
            negative on both local axes (finite differences, b = 1e-4 Ry). Bug: a sign error of the field, or a conjugation error in the
            Bx -+ i By term (chi_xx and chi_yy of opposite sign).
  symmetry  S7: E(q) = E(-q) (weak in centrosymmetric bcc) and E(q + G) = E(q), G = (1,1,0) 2pi/alat (even mesh), constrained SCF at theta = 30 deg.
  chi_orthogonal  the same at kspace_ham_order = 'orthogonal' (see its docstring): negative definite and symmetric.
  symmetry_orthogonal  the same at kspace_ham_order = 'orthogonal'. The E(-q) = E(q) half is weak in bcc (inversion); the E(q + G) half is the
            informative one. Bug caught: q not handled lattice-periodically (injected: q reduced by 10% of
            its nearest-integer direct part, dE(q + G) = -5.3e-4 Ry). Not caught: the lattice matrix transposed (bcc is cubic: E unchanged).
  guards    the k-space GBT constraint guard stops the run with its message (a target other than (0,0,1)).
usage: test_kspace_scspiral.py <mode> <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

from gbt_app_common import SR, tolerance

Q = (0.0, 0.25, 0.0)           # Cartesian, 2 pi/alat: direct (0.125, -0.125, 0.125) of the primitive cell, q/2 on the 16^3 mesh
NK = 16


def deck(dst, nsp, q=None, theta=0.0, constrained=False, nstep=60, bfield=None, order="first", target=(0.0, 0.0, 1.0), ccor=False):
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "fe_bcc", dst)
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["self"]["nstep"] = nstep
    nml["control"]["nsp"] = nsp
    nml["reciprocal"].update(nk1=NK, nk2=NK, nk3=NK, kspace_ham_order=order)
    nml["hamiltonian"] = {"hoh": order == "second"}
    if ccor:       # two-centre CCOR on both sides: Sdot, fitted screening, strict
        nml["lattice"].update(strux_want_sdot=True, screening="fitted")
        nml["hamiltonian"].update(ccor_2c=True, ccor_strict=True)
    if q is not None:
        nml["hamiltonian"].update(magnetic_representation="gbt_single_q", q_ss=list(q), theta_ss=theta)
    nml.write(str(dst / "input.nml"), force=True)
    if constrained:
        with open(dst / "input.nml", "a") as f:
            f.write("\n&constraints\n constraints_enable = .true.\n constraints_i_cons = 3\n constraints_tolerance = 1e-6\n")
            f.write(" constraints_diagnostics = .true.\n" if bfield is not None else "")
            for i, v in enumerate(target, 1):
                f.write(f" constraints_mom_ref({i},1) = {v!r}\n")
            for i, v in enumerate(bfield or (), 1):
                f.write(f" constraints_bfield({i},1) = {v!r}\n")
            f.write("/\n")
    return dst


def run(binary, d):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    return r.returncode, log


def results(d, log):
    energy = float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1))
    return dict(E=energy, EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log,
                dN=max(abs(float(x)) for x in re.findall(r"dN=\s*(\S+),", log)),
                m_perp=float(re.findall(r"max \|m_perp\| \(muB\) =\s*(\S+)", log)[-1]) if "GBT_CONSTRAINT" in log else None,
                residual=float(re.findall(r"residual sum_a \|B_a\|\|m_a_perp\| \(Ry\) =\s*(\S+)", log)[-1]) if "GBT_CONSTRAINT" in log else None)


def scf(binary, d, **kw):
    d = deck(d, **kw)
    code, log = run(binary, d)
    assert code == 0, f"{d}: exit code {code}"
    return results(d, log)


def mode_identity(binary, scratch, order="first", ccor=False):
    tol = tolerance("scf_identity_ry")
    ref = scf(binary, scratch / "collinear", nsp=1, order=order, ccor=ccor)
    spiral = scf(binary, scratch / "spiral", nsp=3, q=Q, theta=0.0, constrained=True, order=order, ccor=ccor)
    dE, dEF = spiral["E"] - ref["E"], spiral["EF"] - ref["EF"]
    print(f"collinear E = {ref['E']:.10f} Ry, E_F = {ref['EF']:.10f} Ry; spiral dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, tolerance {tol:.1e}")
    return spiral["conv"] and abs(dE) < tol and abs(dEF) < tol


def mode_identity_orthogonal(binary, scratch):
    return mode_identity(binary, scratch, order="orthogonal")


def mode_identity_ccor(binary, scratch):
    return mode_identity(binary, scratch, ccor=True)


def mode_identity_ccor_second(binary, scratch):
    return mode_identity(binary, scratch, order="second", ccor=True)


def mode_converge(binary, scratch, order="first"):
    tol, tol_m, tol_n = tolerance("scf_identity_ry"), tolerance("scf_moment_mub"), tolerance("n_identity")
    r = scf(binary, scratch / "cone", nsp=3, q=Q, theta=30.0, constrained=True, order=order)
    print(f"converged {r['conv']}, max |m_perp| = {r['m_perp']:.3e} muB, residual = {r['residual']:.3e} Ry, max |dN| = {r['dN']:.2e}")
    return r["conv"] and r["m_perp"] < tol_m and r["residual"] < tol and r["dN"] < tol_n


def chi_tensor(binary, scratch, q, order, b=1.0e-4):
    """chi[a, c] = dm_perp_a/dB_c [muB/Ry] at fixed potential (one iteration, central differences in the local x, y components)."""
    m_perp = {}
    for tag, field in (("px", (b, 0, 0)), ("mx", (-b, 0, 0)), ("py", (0, b, 0)), ("my", (0, -b, 0))):
        d = deck(scratch / tag, nsp=3, q=q, theta=30.0, constrained=True, nstep=1, bfield=field, order=order)
        assert run(binary, d)[0] == 0
        row = [float(x) for x in [l for l in (d / "constraint_diagnostics.dat").read_text().splitlines() if not l.startswith("#")][0].split()]
        m_perp[tag] = row[8] * np.array(row[5:7])            # moment magnitude * local direction (x, y)
    return np.array([[(m_perp["px"][a] - m_perp["mx"][a]) / (2 * b), (m_perp["py"][a] - m_perp["my"][a]) / (2 * b)] for a in range(2)])


def mode_converge_orthogonal(binary, scratch):
    return mode_converge(binary, scratch, order="orthogonal")


def mode_chi(binary, scratch):
    chi = chi_tensor(binary, scratch, Q, "first")
    print(f"chi [muB/Ry]: xx {chi[0, 0]:.4f}, xy {chi[0, 1]:.4f}, yx {chi[1, 0]:.4f}, yy {chi[1, 1]:.4f}")
    return chi[0, 0] < 0 and chi[1, 1] < 0


def mode_chi_orthogonal(binary, scratch):
    """Same at kspace_ham_order = 'orthogonal', q = (0, 0, 0.125): chi negative definite (both eigenvalues of its symmetric part < 0) and
    symmetric, |chi_xy - chi_yx| b < scf_moment_mub (H is linear in B, so the response is a second derivative of the grand potential).
    The symmetry check holds by construction here: chi_xy = chi_yx = 0 by the mirror symmetry of the cone (measured 0 to 1e-17 muB); only a
    bug that makes chi_xx and chi_yy differ or flip sign is caught, through the eigenvalues.
    Bug caught: a sign error or conjugation error of the field in H (chi_xx, chi_yy of opposite sign, or chi_xy != chi_yx)."""
    b, tol_m = 1.0e-4, tolerance("scf_moment_mub")
    chi = chi_tensor(binary, scratch, (0.0, 0.0, 0.125), "orthogonal", b)
    eig = np.linalg.eigvalsh((chi + chi.T) / 2)
    asym = abs(chi[0, 1] - chi[1, 0]) * b
    print(f"chi [muB/Ry]: xx {chi[0, 0]:.4f}, xy {chi[0, 1]:.4f}, yx {chi[1, 0]:.4f}, yy {chi[1, 1]:.4f}; eigenvalues of the symmetric part "
          f"{eig[0]:.4f}, {eig[1]:.4f}; |chi_xy - chi_yx| b = {asym:.3e} muB, tolerance {tol_m:.1e}")
    return eig[1] < 0 and asym < tol_m


def mode_symmetry(binary, scratch, order="first"):
    tol = tolerance("scf_identity_ry")
    E = {}
    for tag, q in (("q", Q), ("minus_q", (0.0, -Q[1], 0.0)), ("q_plus_G", (1.0, 1.0 + Q[1], 0.0))):
        E[tag] = scf(binary, scratch / tag, nsp=3, q=q, theta=30.0, constrained=True, order=order)["E"]
    d1, d2 = E["minus_q"] - E["q"], E["q_plus_G"] - E["q"]
    print(f"E(q) = {E['q']:.10f} Ry; E(-q) - E(q) = {d1:.3e} Ry (weak: inversion); E(q+G) - E(q) = {d2:.3e} Ry; tolerance {tol:.1e}")
    return abs(d1) < tol and abs(d2) < tol


def mode_symmetry_orthogonal(binary, scratch):
    return mode_symmetry(binary, scratch, order="orthogonal")


def mode_guards(binary, scratch):
    cases = (("target", dict(target=(0.0, 1.0, 0.0)), "needs constraints_mom_ref = 0, 0, 1"),)
    ok = True
    for tag, kw, message in cases:
        d = deck(scratch / tag, nsp=3, q=Q, theta=30.0, constrained=True, nstep=1, **kw)
        code, log = run(binary, d)
        print(f"{tag}: exit code {code}, message found: {message in log}")
        ok = ok and code != 0 and message in log
    return ok


if __name__ == "__main__":
    mode, binary, scratch = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3])
    scratch.mkdir(parents=True, exist_ok=True)
    ok = globals()[f"mode_{mode}"](binary, scratch)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
