"""transport_backend='fast_sp' against the matrix-free loops ('legacy') on periodic fcc Pt (12^3 cell, per_type).

Needs a release build of the full executable: RSLMTO_EXE (default build/bin/rslmto.x).
The loops are the reference. Each case compares cond_total.out where the response is not zero by symmetry.
"""
import os
import shutil
import subprocess
from pathlib import Path

import f90nml
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'tests/postproc/cases/conductivity/fccPt'
EXE = Path(os.environ.get('RSLMTO_EXE', ROOT / 'build/bin/rslmto.x'))
X = [1, 0, 0]
TOL = 1e-6

# name, hoh, linear_out, pol_alpha, v_alpha, v_beta, cond_ll, lld, falls back to legacy
CASES = [
    ('spin', False, 'spin', 'z', [0, 1, 0], [1, 0, 0], 30, 60, False),
    ('charge_xx', False, 'charge', 'z', X, X, 30, 60, False),
    ('hoh_charge_xx', True, 'charge', 'z', X, X, 30, 50, False),
    ('hoh_spin', True, 'spin', 'z', [0, 1, 0], [1, 0, 0], 30, 50, False),
    ('hoh_orbital', True, 'orbital', 'z', [0, 1, 0], [1, 0, 0], 30, 50, False),
    ('hoh_spin_torque', True, 'spin_torque', 'x', [0, 1, 0], [1, 0, 0], 30, 50, True),
]


def run(path, backend, hoh, linear_out, pol, v_alpha, v_beta, ll, lld):
    path.mkdir(parents=True)
    for name in ('input.nml', 'Pt.nml'):
        shutil.copy(FIXTURE / name, path / name)
    nml = f90nml.read(path / 'input.nml')
    nml['lattice'].update({'pbc': True, 'b1': True, 'b2': True, 'b3': True, 'n1': 12, 'n2': 12, 'n3': 12,
                           'strux_backend': 'legacy'})
    nml['hamiltonian'].update({'hoh': hoh, 'pol_alpha': pol, 'pol_beta': 'z', 'v_alpha': v_alpha, 'v_beta': v_beta})
    nml['control'].update({'cond_ll': ll, 'lld': lld, 'cheb_backend': 'legacy', 'transport_backend': backend,
                           'linear_in': 'charge', 'linear_out': linear_out})
    nml.write(path / 'input.nml', force=True)
    result = subprocess.run([str(EXE)], cwd=path, text=True, capture_output=True,
                            env={**os.environ, 'OMP_NUM_THREADS': '2'})
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    return np.loadtxt(path / 'cond_total.out'), result.stdout


@pytest.mark.skipif(not EXE.exists(), reason='set RSLMTO_EXE to a built rslmto.x')
@pytest.mark.parametrize('name,hoh,linear_out,pol,v_alpha,v_beta,ll,lld,fallback', CASES, ids=[c[0] for c in CASES])
def test_fast_sp_matches_legacy(tmp_path, name, hoh, linear_out, pol, v_alpha, v_beta, ll, lld, fallback):
    """Catches: a dropped h*O*v term (relative difference 0.46 without it), a wrong q placement in the
    current or a wrong bond sign (O(1)), and a fast_sp request that silently changes an unsupported setup."""
    legacy, _ = run(tmp_path / 'legacy', 'legacy', hoh, linear_out, pol, v_alpha, v_beta, ll, lld)
    fast, log = run(tmp_path / 'fast', 'fast_sp', hoh, linear_out, pol, v_alpha, v_beta, ll, lld)
    assert ('not available for this setup' in log) == fallback
    if fallback:
        assert np.array_equal(legacy, fast)
        return
    signal = np.abs(legacy[:, 1:]).max(axis=1)
    band = signal > 1e-3 * signal.max()
    scale = np.abs(legacy[band, 1:]).max()
    assert scale > 0
    assert np.abs(legacy[band, 1:] - fast[band, 1:]).max() / scale < TOL
    assert np.abs(fast[:, 2]).max() / np.abs(fast[:, 1]).max() < TOL
