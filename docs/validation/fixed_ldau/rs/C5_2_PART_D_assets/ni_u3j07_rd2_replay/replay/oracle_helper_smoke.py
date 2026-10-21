#!/usr/bin/env python3
"""No-driver import/call smoke for the reviewed d-shell correction helper."""
from fractions import Fraction
import importlib.util
from pathlib import Path
import numpy as np

HELPER = Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/local_moment_capture_oracle.py')
EXPECTED_SHA256 = '98d14387be47c3cea4a4803ce331a1c92cb8c0d9d5d3ca412b5491634b68645d'
import hashlib
h = hashlib.sha256(HELPER.read_bytes()).hexdigest()
if h != EXPECTED_SHA256:
    raise SystemExit(f'helper SHA mismatch: {h}')
spec = importlib.util.spec_from_file_location('rd2_smoke_oracle', HELPER)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
P = np.diag(np.array([1., 2., 3., 4., 5.], dtype=np.complex128))
M1 = np.diag(np.array([3., 7., 10., 16., 21.], dtype=np.complex128))
V = np.diag(np.array([.1, .2, .3, .4, .5], dtype=np.complex128))
q0, center, q2 = 15., 2., 30.
result = module.local_correction_no_controls(P, M1, V, q0, center, q2)
t = Fraction(11, 2)
expected_center = Fraction(2) - t / Fraction(15)
expected_q2 = Fraction(30) - 2 * Fraction(106, 10) + Fraction(225, 100) - t * t / Fraction(15)
if abs(result['c_local_ry'] - float(expected_center)) > 1e-14:
    raise SystemExit('center formula smoke mismatch')
if abs(result['q2_local_ry2'] - float(expected_q2)) > 1e-14:
    raise SystemExit('q2 formula smoke mismatch')
print(f'PASS helper_sha256={h} t={float(t):.17g} center={result["c_local_ry"]:.17g} q2={result["q2_local_ry2"]:.17g}')
