#!/usr/bin/env python3
"""Synthetic evaluate-to-JSON regression for the A4.1 report path."""
from __future__ import annotations

import importlib.util
import json
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('a4_w1', ROOT / 'a4_m21_w1_richardson.py')
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

G0 = 6.4984804203277395
HS = (0.0025, 0.005, 0.01)
LAMBDAS = mod.LAMBDAS
MODES = mod.MODES


def capture_data(kind: str, intercept: float = G0) -> dict:
    # Four exact entropy states keep this end-to-end test fast while exercising
    # the same source-bound rectangular table and report serialization path.
    states = [{'k_index': 0, 'band_index': 0,
               'occupation_text': '0.0000000000000000E+00',
               'weight_text': '1.0000000000000000E+00',
               'eigenvalue_text': '0.0000000000000000E+00'}]
    entropy = {'temperature_k_text': '3.0000000000000000E+02',
               'temperature_ry_text': '1.9000000000000000E-03',
               'fermi_level_text': '0.0000000000000000E+00',
               'states': states, 'expected_n_kpoints': 1,
               'expected_band_indices': [0],
               'band_index_semantics': 'zero-based global spin-inclusive row of eigenvalues(:,ik)'}
    # The endpoint occupation has exactly zero entropy, avoiding any extra
    # transcendental or entropy-model behavior in this report-path test.
    ts = Decimal(0)
    ts_text = '0.0000000000000000E+00'
    p = [[[[G0, 0.0] if i == j == 0 else [0.0, 0.0] for j in range(5)]
          for i in range(5)]]
    p_bound = [[1.0e-14 for _ in range(5)] for _ in range(5)]
    if kind == 'unresolved':
        p_bound = [[1.0e-3 for _ in range(5)] for _ in range(5)]

    def dvalue(h: float) -> float:
        if kind == 'dominant':
            return intercept + 10.0*h*h + 32.0*h**4
        if kind == 'unresolved':
            return intercept + h*h + 20000.0*h**4
        if kind == 'identity_fail':
            return intercept + 0.01 + 10.0*h*h + 32.0*h**4
        if kind == 'h4_stop':
            return intercept + h*h + 20000.0*h**4
        raise ValueError(kind)

    captures = []
    for mode in MODES:
        for lam in LAMBDAS:
            for threshold in ('standard', 'tight'):
                ftarget = 0.0 if lam == 0 else lam * dvalue(abs(lam))
                # Local F is ETOT_local + lambda*Tr(WP) - TS; spectral F is
                # ETOT_spectral - TS.  This synthetic keeps both candidates
                # equal at the U=0 external-coupling model.
                etot_spec = ftarget + float(ts)
                etot_local = etot_spec - lam * G0
                captures.append({
                    'mode': mode, 'lambda_ry': lam, 'threshold': threshold,
                    'P_d': p[0], 'P_d_component_error_bound': p_bound,
                    'ETOT_spectral_text': f'{etot_spec:.16E}',
                    'ETOT_local_text': f'{etot_local:.16E}',
                    'TS_text': ts_text, 'lambda_runtime_text': f'{lam:.17E}',
                    'entropy_inputs': entropy,
                })
    w1 = [[[1.0 if i == j else 0.0, 0.0] for j in range(5)] for i in range(5)]
    return {'W1': w1, 'captures': captures}


def run(kind: str, intercept: float = G0) -> dict:
    result = mod.evaluate(capture_data(kind, intercept))
    # Exercise precisely the JSON serialization that failed on production data.
    return json.loads(json.dumps(result, allow_nan=False, sort_keys=True))


def main() -> None:
    passed = run('dominant')
    for mode in MODES:
        row = passed['modes'][mode]
        assert row['identity_pass'] is True, (mode, row['identity_classification'], row['scaling'])
        assert row['scaling']['resolved'] is True, row['scaling']
        assert 'X' in row['scaling'] and 'Y' in row['scaling']

    failed_identity = run('identity_fail')
    assert all(failed_identity['modes'][m]['scaling']['resolved'] for m in MODES)
    assert all(failed_identity['modes'][m]['identity_classification'] == 'FAIL' for m in MODES)

    stopped = run('unresolved')
    assert all(not stopped['modes'][m]['scaling']['resolved'] for m in MODES)
    assert all(stopped['modes'][m]['identity_classification'] == 'INCONCLUSIVE_STOP' for m in MODES)

    h4_stop = run('h4_stop')
    assert all(not h4_stop['modes'][m]['scaling']['resolved'] for m in MODES)
    print('PASS evaluate-to-JSON: dominant identity PASS; resolved identity FAIL; noisy and h4-dominant STOP')


if __name__ == '__main__':
    main()
