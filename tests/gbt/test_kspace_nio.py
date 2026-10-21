#!/usr/bin/env python3
"""AppGbtKspaceNio: antiferromagnetic NiO (fixture nio_af, four sites, hoh + ccor_2c strict, screening fitted) on the k-space GBT path.
The deck encodes the down-spin Ni as 'mirrored' (Nid: mom = +z, spin 1 holds the minority values). The B1 convention of the
k-space glue (a site with mom(3) < 0, spin 1 holds the majority values, the frame carries the orientation) is the file Nid_B1.nml
made by make_nid_b1.py. Both use ccor_vmt_mode = 'surface_scalar'; 'pair_surface' is refused by the k-space GBT glue. Modes, each with its
relation and the bug it catches:
  r2        R2 (Sandratskii, Adv. Phys. 47, 91 (1998)) at theta = 0, fixed potentials, both orders. (a) q != 0: the GBT eigenvalues at
            k equal the collinear periodic_nc (nsp = 1, mirrored Nid) sector eigenvalues H_up(k - q/2) and H_dn(k + q/2), with the
            mirrored Nid and with the B1 Nid. (b) q = 0: GBT with the B1 Nid equals periodic_nc nsp = 3 with the same B1 file.
            Bug caught: a sign +1 convention of W, C, ccd, obar or E_nu that a mom(3) < 0 site breaks with CCOR, a mom, center_band,
            width_band, shifted_band or obar field of the B1 file left unswapped, Hcc inside h - h o h.
  i2        I2 (Sandratskii 1998): self-consistent, unconstrained, theta = 0, q != 0, second order. The GBT SCF equals the collinear
            nsp = 1 SCF of the deck (mirrored Nid) in E and E_F (scf_identity_ry) and in the four site moments (scf_moment_mub); GBT with
            the mirrored Nid and with the B1 Nid. Bug caught: Hcc missing or misplaced in the spiral SCF, a spin channel of a mom(3) < 0
            site taken the wrong way round in the density or in Hcc.
  guard     ccor_vmt_mode = 'pair_surface' with gbt_single_q stops with the message of the glue. Bug: the guard removed, the run then
            gives eigenvalues from an undefined channel-wise pair average.
usage: test_kspace_nio.py <mode> <rslmto.x> <scratch directory>"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from gbt_app_common import make_nio_deck, make_nio_scf_deck, mesh_index, read_export, run, tolerance

NK = (8, 8, 8)
NK0 = (4, 4, 4)
Q = (0.25, 0.5, 0.25)          # direct; q/2 times NK is an integer vector


def sector_deviation(e_gbt, k_gbt, spin, index, shift):
    m = np.rint(k_gbt * np.array(NK)).astype(int)
    dev = 0.0
    for ik in range(len(k_gbt)):
        up = spin[0, :, index[tuple((m[ik] - shift) % NK)]]
        dn = spin[1, :, index[tuple((m[ik] + shift) % NK)]]
        dev = max(dev, np.abs(np.sort(e_gbt[:, ik]) - np.sort(np.concatenate([up, dn]))).max())
    return dev


def mode_r2(binary, scratch):
    tol = tolerance("eig_identity_ry")
    shift = 0.5 * np.array(Q) * np.array(NK)
    assert np.allclose(shift, np.rint(shift)), "q/2 is not a mesh translation"
    shift = np.rint(shift).astype(int)
    worst = 0.0
    for order in ("first", "second"):
        ref = make_nio_deck(scratch / f"ref_{order}", "mirrored", "periodic_nc", NK, order=order, spin_resolved=True)
        run(binary, ref)
        k_ref, spin = read_export(ref, spin=True)
        index = mesh_index(k_ref, NK)
        shutil.rmtree(ref)
        for nid in ("mirrored", "b1"):
            gbt = make_nio_deck(scratch / f"gbt_{nid}_{order}", nid, "gbt_single_q", NK, q_dir=Q, order=order)
            run(binary, gbt)
            k_gbt, e_gbt = read_export(gbt)
            shutil.rmtree(gbt)
            dev = sector_deviation(e_gbt, k_gbt, spin, index, shift)
            print(f"(a) q != 0, {order}, GBT Nid {nid}: {len(k_gbt)} k-points, max |E_gbt - E_nc(k -+ q/2)| = {dev:.3e} Ry")
            worst = max(worst, dev)
        ref = make_nio_deck(scratch / f"ref0_{order}", "b1", "periodic_nc", NK0, order=order)
        gbt = make_nio_deck(scratch / f"gbt0_{order}", "b1", "gbt_single_q", NK0, order=order)
        run(binary, ref)
        run(binary, gbt)
        dev = np.abs(np.sort(read_export(ref)[1], axis=0) - np.sort(read_export(gbt)[1], axis=0)).max()
        shutil.rmtree(ref)
        shutil.rmtree(gbt)
        print(f"(b) q = 0, {order}, Nid B1 on both sides: max |E_gbt - E_nc| = {dev:.3e} Ry")
        worst = max(worst, dev)
    print(f"max over all cases {worst:.3e} Ry, tolerance {tol:.1e}")
    return worst < tol


def scf(binary, d):
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    log = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    (d / "run.log").write_text(log)
    assert r.returncode == 0, f"{d}: exit code {r.returncode}"
    return dict(E=float(re.search(r"Physical DFT total energy of system:\s+(\S+)", (d / "report.out").read_text()).group(1)),
                EF=float(re.findall(r"EF=\s*(\S+) Ry", log)[-1]), conv="Converged!" in log,
                moments=np.array([float(m) for _, m in re.findall(r"Spin moment of atom\s+(\d+) is\s+(\S+)", log)[-4:]]))


def mode_i2(binary, scratch):
    tol, tol_m = tolerance("scf_identity_ry"), tolerance("scf_moment_mub")
    q = (1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0)               # q/2 times (6, 6, 6) is an integer vector
    nk = (6, 6, 6)
    ref = scf(binary, make_nio_scf_deck(scratch / "collinear", "mirrored", False, nk))
    print(f"collinear: converged {ref['conv']}, E = {ref['E']:.10f} Ry, E_F = {ref['EF']:.10f} Ry, moments {ref['moments']} muB")
    ok = ref["conv"]
    if ref["moments"][:2].min() < 0.1:
        sys.exit("the collinear NiO state is (nearly) non-magnetic; the identity would be trivial; stopping")
    for nid in ("mirrored", "b1"):
        gbt = scf(binary, make_nio_scf_deck(scratch / f"gbt_{nid}", nid, True, nk, q_dir=q))
        dE, dEF, dm = gbt["E"] - ref["E"], gbt["EF"] - ref["EF"], np.abs(gbt["moments"] - ref["moments"]).max()
        print(f"GBT Nid {nid}: converged {gbt['conv']}; dE = {dE:.3e} Ry, dE_F = {dEF:.3e} Ry, max |dm| = {dm:.3e} muB "
              f"(tolerances {tol:.1e}, {tol_m:.1e}); moments {gbt['moments']} muB")
        ok = ok and gbt["conv"] and abs(dE) < tol and abs(dEF) < tol and dm < tol_m
    return ok


def mode_guard(binary, scratch):
    d = make_nio_deck(scratch / "pair_surface", "b1", "gbt_single_q", NK0, vmt_mode="pair_surface")
    r = subprocess.run([str(binary), "input.nml"], cwd=d, capture_output=True, text=True)
    found = "channel-wise pair average" in r.stdout + r.stderr
    print(f"pair_surface: exit code {r.returncode}, message found: {found}")
    return r.returncode != 0 and found


MODES = {"r2": mode_r2, "i2": mode_i2, "guard": mode_guard}

if __name__ == "__main__":
    mode, binary, scratch = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    sys.exit(0 if MODES[mode](binary, scratch) else 1)
