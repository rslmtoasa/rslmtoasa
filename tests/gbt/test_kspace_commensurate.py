#!/usr/bin/env python3
"""AppGbtKspaceCommensurate. Relation R3 (Sandratskii, Adv. Phys. 47, 91 (1998)): if the spiral is periodic in an n-site
supercell, the GBT eigenvalues at the primitive k that fold onto a supercell point K equal the eigenvalues of the ordinary
non-collinear supercell Hamiltonian at K, same potential, moments m(R + tau). K_z = k_z + M/(2L) (M = q_z L, L the
supercell length in alat), because the translation by L carries D(2 pi M) = (-1)^M. Reference: periodic_nc on the
fixtures tests/gbt/fixtures/bccFe_super_q050 and _q033 (every type the Fe potential, moment rotated). Precondition,
checked first: the supercell S blocks equal the primitive ones for corresponding bonds (gbt_bonds.out of both).
Bug caught: any operator error at finite q and theta, in either order.
With the optional argument ccor, every deck carries ccor_2c (ccor_strict, screening fitted, Sdot). The supercell side builds its
two-centre block from explicit moment vectors (build_ccor_pair_block_noncollinear), the GBT side by lifting S and Sdot with the
spin link, so they are different constructions. Bug caught there: the link missing on Sdot or on the onsite term at finite theta,
Hcc inside h - h o h. The S-block precondition does not cover Sdot, which gbt_bonds.out does not dump.
usage: test_kspace_commensurate.py <rslmto.x> <scratch directory> <first|second> [ccor]"""
import sys
from pathlib import Path

import numpy as np

from gbt_app_common import make_deck, make_super_deck, mesh_index, read_bonds, read_export, run, tolerance

binary, scratch, order = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
ccor = len(sys.argv) > 4 and sys.argv[4] == "ccor"
tol = tolerance("eig_identity_ry")
FIXTURES = {"bccFe_super_q050": ([0.0, 0.5, 1.0, 1.5], 0.5), "bccFe_super_q033": ([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], 1.0 / 3.0)}
NK_SUPER = (4, 4, 2)
worst = 0.0
for fixture, (z, q_z) in FIXTURES.items():
    length = round(z[-1] + 0.5)                     # supercell length in alat; M = q_z L = 1
    winding = round(q_z * length)
    nk_prim = (4, 4, 2 * length)
    prim = make_deck(scratch / f"{fixture}_prim_bonds", "bccFe_sc2", "gbt_single_q", (2, 2, 2), ccor=ccor)
    sup = make_super_deck(scratch / f"{fixture}_super_bonds", fixture, z, q_z, 0.0, "first", (2, 2, 2), gbt=True, ccor=ccor)
    run(binary, prim)
    run(binary, sup)
    p_bonds, s_bonds = read_bonds(prim), read_bonds(sup)
    pre = 0.0
    for atom, d, block in s_bonds:                  # type i of the supercell is sublattice 1 (odd i) or 2 (even i)
        match = [b for a, dd, b in p_bonds if a == 2 - atom % 2 and np.abs(dd - d).max() < 1e-9]
        assert len(match) == 1, f"no primitive bond for supercell type {atom}, d = {d}"
        pre = max(pre, np.abs(match[0] - block).max())
    print(f"{fixture}: precondition, {len(s_bonds)} supercell bonds vs primitive bonds, max |S_super - S_prim| = {pre:.3e}")
    if pre >= tol:
        sys.exit("precondition failed: supercell and primitive S blocks differ; stopping")
    for theta in (90.0, 30.0):
        tag = f"{fixture}_th{theta:.0f}_{order}"
        ref = make_super_deck(scratch / f"{tag}_nc", fixture, z, q_z, theta, order, NK_SUPER, ccor=ccor)
        gbt = make_deck(scratch / f"{tag}_gbt", "bccFe_sc2", "gbt_single_q", nk_prim, q_dir=(0.0, 0.0, q_z), order=order, theta_deg=theta, ccor=ccor)
        run(binary, ref)
        run(binary, gbt)
        (k_s, e_s), (k_g, e_g) = read_export(ref), read_export(gbt)
        index, dev = mesh_index(k_g, nk_prim), 0.0
        for ik, (m1, m2, m3) in enumerate(np.rint(k_s * np.array(NK_SUPER)).astype(int)):
            folded = [index[(m1 % 4, m2 % 4, (m3 - winding + 2 * n) % (2 * length))] for n in range(length)]
            dev = max(dev, np.abs(np.sort(e_s[:, ik]) - np.sort(np.concatenate([e_g[:, j] for j in folded]))).max())
        print(f"{tag}: {len(k_s)} supercell k-points, max |E_super - E_gbt(folded)| = {dev:.3e} Ry")
        worst = max(worst, dev)
print(f"max over all cases {worst:.3e} Ry, tolerance {tol:.1e}")
sys.exit(0 if worst < tol else 1)
