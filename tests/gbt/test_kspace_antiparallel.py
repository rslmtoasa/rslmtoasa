#!/usr/bin/env python3
"""AppGbtKspaceAntiparallel. Relations R2 and R3 (Sandratskii, Adv. Phys. 47, 91 (1998)): at q = 0, theta = 0 the GBT
operator of an antiparallel two-site cell equals the ordinary non-collinear operator. Reference: periodic_nc on the same
bcc Cr deck (Cr2 = Cr1 potential with mom = -z). Bug caught: an F1-type sign error for a site with mom(3) < 0.
usage: test_kspace_antiparallel.py <rslmto.x> <scratch directory>"""
import sys
from pathlib import Path

import numpy as np

from gbt_app_common import make_deck, read_export, run, tolerance

binary, scratch = sys.argv[1], Path(sys.argv[2])
decks = {rep: make_deck(scratch / rep, "cr_afm", rep, (4, 4, 4)) for rep in ("periodic_nc", "gbt_single_q")}
for deck in decks.values():
    run(binary, deck)
(k_nc, e_nc), (k_gbt, e_gbt) = read_export(decks["periodic_nc"]), read_export(decks["gbt_single_q"])
assert np.allclose(k_nc, k_gbt)
dev = np.abs(np.sort(e_nc, axis=0) - np.sort(e_gbt, axis=0)).max()
print(f"bcc Cr cell, q = 0, theta = 0: {len(k_gbt)} k-points, max |E_gbt - E_nc| = {dev:.3e} Ry, tolerance {tolerance('eig_identity_ry'):.1e}")
sys.exit(0 if dev < tolerance("eig_identity_ry") else 1)
