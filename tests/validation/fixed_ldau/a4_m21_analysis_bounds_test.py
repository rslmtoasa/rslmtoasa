#!/usr/bin/env python3
"""Deterministic nonproduction checks for A4.1 analyzer bounds."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODULE_PATH = ROOT / 'a4_m21_w1_richardson.py'
spec = importlib.util.spec_from_file_location('a4_w1', MODULE_PATH)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def main() -> None:
    gamma_check = mod.verify_gamma_bounds()
    assert len(gamma_check['checks']) == 7
    for row in gamma_check['checks']:
        assert row['verified_upper_bound']

    hs = (0.0025, 0.005, 0.01)
    eps = mod.EPS
    # Pure h^2 response resolves under the approved leading-term gate.
    pure_h2 = {h: 1.25 + 2.0*h*h for h in hs}
    bound_h2 = {h: 64.0*eps*abs(pure_h2[h]) for h in hs}
    result_h2 = mod.h2_resolution(pure_h2, bound_h2)
    assert result_h2['resolved'], result_h2

    # A resolved, subordinate h^4 term passes the approved dominance rule.
    dominant = {h: 1.25 + 10.0*h*h + 32.0*h**4 for h in hs}
    bound_dom = {h: 64.0*eps*abs(dominant[h]) for h in hs}
    result_dom = mod.h2_resolution(dominant, bound_dom)
    assert result_dom['resolved'], result_dom

    # h^4 dominance must stop.
    h4_dominant = {h: 1.25 + h*h + 20000.0*h**4 for h in hs}
    bound_h4 = {h: 64.0*eps*abs(h4_dominant[h]) for h in hs}
    result_h4 = mod.h2_resolution(h4_dominant, bound_h4)
    assert result_h4['resolved'] is False, result_h4

    # A mathematically dominant signal unresolved by measurement uncertainty stops.
    noisy = {h: 1.25 + 10.0*h*h + 32.0*h**4 for h in hs}
    noisy_bound = {h: 1e-4 for h in hs}
    result_noisy = mod.h2_resolution(noisy, noisy_bound)
    assert result_noisy['resolved'] is False, result_noisy

    # Near cancellation/zero leading coefficient cannot spuriously pass.
    zero_leading = {h: 1.25 + 32.0*h**4 for h in hs}
    zero_bound = {h: 64.0*eps*abs(zero_leading[h]) for h in hs}
    result_zero = mod.h2_resolution(zero_leading, zero_bound)
    assert result_zero['resolved'] is False, result_zero

    # Exact-binary comparison rejects equality and distinguishes adjacent floats.
    resolved, lhs, rhs = mod.h2_dominance_compare(1.0, 0.5, 0.03125, 0.0)
    assert resolved and lhs == rhs and not (lhs > rhs)
    above = math.nextafter(1.0, math.inf)
    _, lhs_above, rhs_above = mod.h2_dominance_compare(above, 0.5, 0.03125, 0.0)
    assert lhs_above > rhs_above
    below = math.nextafter(1.0, -math.inf)
    _, lhs_below, rhs_below = mod.h2_dominance_compare(below, 0.5, 0.03125, 0.0)
    assert lhs_below < rhs_below

    def entropy_spec(k_indices: tuple[int, int]) -> dict:
        states = []
        for k in k_indices:
            for band in (0, 1):
                states.append({'k_index': k, 'band_index': band,
                               'occupation_text': '5.0E-01', 'weight_text': '5.0E-01',
                               'eigenvalue_text': '0.0E+00'})
        return {'temperature_k_text': '3.0E+02', 'temperature_ry_text': '1.9E-03',
                'fermi_level_text': '0.0E+00', 'TS_text': '0.0E+00', 'states': states,
                'expected_n_kpoints': 2, 'expected_band_indices': [0, 1],
                'band_index_semantics': 'zero-based global spin-inclusive row of eigenvalues(:,ik)'}

    mod.entropy_roundoff_bound(entropy_spec((0, 1)))
    try:
        mod.entropy_roundoff_bound(entropy_spec((1, 2)))
    except ValueError as exc:
        assert 'zero-based source range' in str(exc)
    else:
        raise AssertionError('shifted contiguous global-k domain was accepted')
    print('PASS gamma bounds; pure h2 and dominant h2+h4 resolve; h4-dominant, noisy, and zero-leading cases stop')


if __name__ == '__main__':
    main()
