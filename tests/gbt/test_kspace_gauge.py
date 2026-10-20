#!/usr/bin/env python3
"""AppGbtKspaceGauge. Relation R2 (Sandratskii, Adv. Phys. 47, 91 (1998)): at theta = 0 the GBT operator is block
diagonal, spin up = H_up(k - q/2) and spin down = H_down(k + q/2) of the collinear operator. Reference: the
periodic_nc k-space route with independently solved spin sectors (consistency of two production routes, plus the
R2 identity). Bug caught: a wrong k -+ q/2 assignment or q-to-direct conversion, a wrong S extraction or mesh mapping.
(An orbital basis change is unitary and invisible to eigenvalues.)
With the optional argument ccor both sides carry ccor_2c (ccor_strict, screening fitted, Sdot); the two routes then differ in how
the two-centre block gets its link: GBT lifts S and Sdot with U_a^dag diag(e^{i(k-q/2).d}, e^{i(k+q/2).d}) U_b, the collinear
route Bloch-sums its own spin sectors. A wrong Sdot link or transposition, or Hcc inside h - h o h at second order, shows up as a
difference. The shared ccd and lambda helpers feed both sides, so a wrong value of those is not visible here.
usage: test_kspace_gauge.py <rslmto.x> <scratch directory> [ccor]"""
import itertools
import sys
from pathlib import Path

import numpy as np

from gbt_app_common import make_deck, mesh_index, read_export, run, tolerance

binary, scratch = sys.argv[1], Path(sys.argv[2])
ccor = len(sys.argv) > 3 and sys.argv[3] == "ccor"
CASES = [((0.0, 0.0, 0.25), (4, 4, 8)), ((0.1, 0.2, 0.3), (20, 10, 20))]     # q direct; mesh with q/2 on it
tol = tolerance("eig_identity_ry")
worst = 0.0
for order, kind, (iq, (q, nk)) in itertools.product(("first", "second"), ("bccFe_prim", "bccFe_sc2"), enumerate(CASES)):
    tag = f"{kind}_q{iq}_{order}"
    ref = make_deck(scratch / f"{tag}_nc", kind, "periodic_nc", nk, order=order, spin_resolved=True, ccor=ccor)
    gbt = make_deck(scratch / f"{tag}_gbt", kind, "gbt_single_q", nk, q_dir=q, order=order, ccor=ccor)
    run(binary, ref)
    run(binary, gbt)
    k_ref, spin = read_export(ref, spin=True)
    k_gbt, e_gbt = read_export(gbt)
    index = mesh_index(k_ref, nk)
    shift = 0.5 * np.array(q) * np.array(nk)
    assert np.allclose(shift, np.rint(shift)), "q/2 is not a mesh translation"
    shift = np.rint(shift).astype(int)
    m = np.rint(k_gbt * np.array(nk)).astype(int)
    dev = 0.0
    for ik in range(len(k_gbt)):
        up = spin[0, :, index[tuple((m[ik] - shift) % nk)]]
        dn = spin[1, :, index[tuple((m[ik] + shift) % nk)]]
        dev = max(dev, np.abs(np.sort(e_gbt[:, ik]) - np.sort(np.concatenate([up, dn]))).max())
    print(f"{tag}: {len(k_gbt)} k-points, max |E_gbt - E_nc(k -+ q/2)| = {dev:.3e} Ry")
    worst = max(worst, dev)
print(f"max over all cases {worst:.3e} Ry, tolerance {tol:.1e}")
sys.exit(0 if worst < tol else 1)
