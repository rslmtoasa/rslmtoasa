#!/usr/bin/env python3
"""AppPeriodicNcRotation. Relation: global spin-rotation invariance of the spin-only Hamiltonian without SOC; rotating all
moments of a ferromagnet leaves the spectrum unchanged. Reference: the unrotated spectrum of the same periodic_nc deck.
Bug caught: single-precision or inconsistent handling of moment components in the onsite operators (the cmplx-without-kind
class in build_obarm/build_enim). One-site bcc Fe, first and second order, full 8x8x8 mesh.
usage: test_periodic_nc_rotation.py <rslmto.x> <scratch directory>"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gbt"))
from gbt_app_common import make_deck, read_export, run, tolerance  # noqa: E402

binary, scratch = sys.argv[1], Path(sys.argv[2])
tol, worst = tolerance("eig_identity_ry"), 0.0
for order in ("first", "second"):
    spectra = {}
    for theta, phi in ((0, 0), (10, 0), (30, 0), (45, 45), (60, 120)):
        deck = make_deck(scratch / f"{order}_{theta}_{phi}", "bccFe_prim", "periodic_nc", (8, 8, 8), order=order)
        t, p = np.radians(theta), np.radians(phi)
        mom = ", ".join(repr(float(x)) for x in (np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)))
        fe = deck / "Fe.nml"
        fe.write_text(re.sub(r"^ mom = .*$", f" mom = {mom}", fe.read_text(), flags=re.M))
        run(binary, deck)
        spectra[(theta, phi)] = np.sort(read_export(deck)[1], axis=0)
        dev = np.abs(spectra[(theta, phi)] - spectra[(0, 0)]).max()
        print(f"{order} order, theta = {theta}, phi = {phi}: max |E - E(unrotated)| = {dev:.3e} Ry")
        worst = max(worst, dev)
print(f"max over all rotations {worst:.3e} Ry, tolerance {tol:.1e}")
sys.exit(0 if worst < tol else 1)
