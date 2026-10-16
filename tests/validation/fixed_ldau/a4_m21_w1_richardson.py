#!/usr/bin/env python3
"""Evaluate the W1 identity with source-format and operation-derived bounds.

Consumes only hash-bound captures; it never launches SCF. Input tokens retain
the production ES/e24.16 decimal spelling so formatting uncertainty is derived
from each value's actual exponent rather than a fixed absolute floor.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path

EPS = 2.220446049250313e-16
HS = (0.0025, 0.005, 0.01)
LAMBDAS = (0.0, -0.0025, 0.0025, -0.005, 0.005, -0.01, 0.01)
MODES = ('spectral', 'local')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gamma(n: int) -> float:
    ne = n * EPS
    if ne >= 1.0:
        raise ValueError('operation count too large for gamma_n bound')
    return math.nextafter(ne / (1.0 - ne), math.inf)


def verify_gamma_bounds() -> dict:
    """Check each operation-count bound against its defining exact ratio."""
    checked = (1, 2, 3, 4, 64, 128, 256)
    rows = []
    with localcontext() as ctx:
        ctx.prec = 80
        eps = Decimal.from_float(EPS)
        for n in checked:
            exact = Decimal(n) * eps / (Decimal(1) - Decimal(n) * eps)
            got = Decimal.from_float(gamma(n))
            if got < exact:
                raise ValueError(f'gamma_{n} rounded below its defining upper bound')
            rows.append({'n': n, 'gamma_n': repr(gamma(n)),
                         'exact_ratio_decimal': str(exact), 'verified_upper_bound': True})
    return {'epsilon_binary64': repr(EPS), 'checks': rows}


def upward(value: float) -> float:
    if not math.isfinite(value) or value < 0.0:
        raise ValueError('cannot outward-round a nonfinite/negative bound')
    return math.nextafter(value, math.inf) if value else 0.0


def downward(value: float) -> float:
    if not math.isfinite(value) or value < 0.0:
        raise ValueError('cannot outward-round a nonfinite/negative bound')
    return math.nextafter(value, -math.inf) if value else 0.0


def q_decimal(token: str) -> Decimal:
    """Exact half quantum from the token's printed mantissa precision."""
    text = str(token).strip().upper().replace('D', 'E')
    if not re.fullmatch(r'[+-]?(?:\d+\.?\d*|\.\d+)E[+-]?\d+', text):
        raise ValueError(f'not a serialized scientific-format number: {token!r}')
    mantissa, exponent = text.split('E')
    decimals = len(mantissa.partition('.')[2])
    return Decimal(5).scaleb(int(exponent) - decimals - 1)


def q_decimal_token(token: str) -> float:
    exact = q_decimal(token)
    rounded = float(exact)
    if Decimal.from_float(rounded) < exact:
        rounded = math.nextafter(rounded, math.inf)
    return rounded


def _dec_ln(x: Decimal, precision: int) -> Decimal:
    # Decimal.ln is correctly rounded to the active context precision.
    with localcontext() as ctx:
        ctx.prec = precision
        return x.ln()


def _outward(value: Decimal, precision: int) -> tuple[Decimal, Decimal]:
    """One-ulp enclosure around a correctly-rounded Decimal transcendental."""
    ulp = Decimal(1).scaleb(value.adjusted() - precision + 1) if value else Decimal(1).scaleb(-precision)
    return (_sub(value, ulp, precision, ROUND_FLOOR),
            _add(value, ulp, precision, ROUND_CEILING))


def _sub(a: Decimal, b: Decimal, precision: int, rounding: str) -> Decimal:
    with localcontext() as ctx:
        ctx.prec, ctx.rounding = precision, rounding
        return a - b


def _add(a: Decimal, b: Decimal, precision: int, rounding: str) -> Decimal:
    with localcontext() as ctx:
        ctx.prec, ctx.rounding = precision, rounding
        return a + b


def _mul(a: Decimal, b: Decimal, precision: int, rounding: str) -> Decimal:
    with localcontext() as ctx:
        ctx.prec, ctx.rounding = precision, rounding
        return a * b


def _h_interval_from_f(flo: Decimal, fhi: Decimal, precision: int) -> tuple[Decimal, Decimal]:
    """Outward entropy bounds over an occupation interval [flo,fhi]."""
    with localcontext() as ctx:
        ctx.prec = precision
        zero = Decimal(0)
        flo, fhi = max(zero, flo), min(Decimal(1), fhi)
        if flo > fhi:
            raise ValueError('invalid occupation interval')
        def hpoint(f: Decimal) -> tuple[Decimal, Decimal]:
            if f <= 0 or f >= 1:
                return zero, zero
            lf_lo = _outward(_dec_ln(f, precision), precision)[0]
            lf_hi = _outward(_dec_ln(f, precision), precision)[1]
            one_f = _sub(Decimal(1), f, precision, ROUND_FLOOR)
            l1_lo = _outward(_dec_ln(one_f, precision), precision)[0]
            l1_hi = _outward(_dec_ln(one_f, precision), precision)[1]
            p1_lo = _mul(f, lf_lo, precision, ROUND_FLOOR)
            p1_hi = _mul(f, lf_hi, precision, ROUND_CEILING)
            p2_lo = _mul(one_f, l1_lo, precision, ROUND_FLOOR)
            p2_hi = _mul(one_f, l1_hi, precision, ROUND_CEILING)
            sum_lo = _add(p1_lo, p2_lo, precision, ROUND_FLOOR)
            sum_hi = _add(p1_hi, p2_hi, precision, ROUND_CEILING)
            return (max(zero, _sub(zero, sum_hi, precision, ROUND_FLOOR)),
                    _sub(zero, sum_lo, precision, ROUND_CEILING))
        at_lo, at_hi = hpoint(flo), hpoint(fhi)
        hmin = min(at_lo[0], at_hi[0])
        hmax = max(at_lo[1], at_hi[1])
        if flo <= Decimal('0.5') <= fhi:
            hmax = _outward(_dec_ln(Decimal(2), precision), precision)[1]
        return hmin, hmax


def matrix(raw: list) -> list[list[complex]]:
    if len(raw) != 5 or any(len(row) != 5 for row in raw):
        raise ValueError('expected 5x5 P matrix encoded as [real,imag] pairs')
    out = []
    for row in raw:
        vals = []
        for pair in row:
            if len(pair) != 2 or not all(math.isfinite(float(x)) for x in pair):
                raise ValueError('nonfinite/bad P matrix element')
            vals.append(complex(float(pair[0]), float(pair[1])))
        out.append(vals)
    return out


def matrix_bounds(raw: list) -> list[list[float]]:
    if len(raw) != 5 or any(len(row) != 5 for row in raw):
        raise ValueError('P_d_component_error_bound must be 5x5')
    out = [[float(x) for x in row] for row in raw]
    if any((not math.isfinite(x) or x < 0.0) for row in out for x in row):
        raise ValueError('invalid P matrix component bound')
    return out


def trace_wp(w: list[list[complex]], p: list[list[complex]]) -> float:
    z = sum(w[i][j] * p[j][i] for i in range(5) for j in range(5))
    return z.real


def trace_abs_scale(w: list[list[complex]], p: list[list[complex]]) -> float:
    return sum(abs(w[i][j] * p[j][i]) for i in range(5) for j in range(5))


def trace_input_bound(w: list[list[complex]], p_bounds: list[list[float]]) -> float:
    return sum(abs(w[i][j]) * p_bounds[j][i] for i in range(5) for j in range(5))


def entropy_roundoff_bound(spec: dict) -> tuple[float, dict]:
    """Enclose TS using the actual captured occupations and k weights.

    The trace writer and fermi_dirac_entropy call the same occupation kernel
    on the same eigenvalue/EF/effective-kT table. Thus serialized occupation
    intervals directly bound the f inputs without an independent libm model.
    """
    precision = 100
    for name in ('temperature_k_text', 'temperature_ry_text', 'fermi_level_text', 'TS_text'):
        if not isinstance(spec.get(name), str):
            raise ValueError(f'{name} must preserve its original serialized token')
    tk_text = spec['temperature_k_text']
    qtk = q_decimal(tk_text)
    temp_k = Decimal(tk_text.replace('D','E'))
    if temp_k <= qtk:
        raise ValueError('serialized temperature interval crosses zero')
    states = spec['states']
    if not states:
        raise ValueError('entropy state table is empty')
    keys: set[tuple[int, int]] = set()
    bands_by_k: dict[int, set[int]] = {}
    weights_by_k: dict[int, str] = {}
    for state in states:
        if not isinstance(state.get('k_index'), int) or not isinstance(state.get('band_index'), int):
            raise ValueError('entropy states require integer global k and band indices')
        k, band = state['k_index'], state['band_index']
        if k < 0 or band < 0 or (k, band) in keys:
            raise ValueError('negative or duplicate global k/band entropy row')
        keys.add((k, band))
        bands_by_k.setdefault(k, set()).add(band)
        wt = state.get('weight_text')
        if not isinstance(wt, str):
            raise ValueError('k weight must preserve its raw serialized token')
        if k in weights_by_k and weights_by_k[k] != wt:
            raise ValueError('k-weight tokens disagree across bands at one global k')
        weights_by_k[k] = wt
    k_indices = sorted(bands_by_k)
    if k_indices != list(range(k_indices[0], k_indices[-1] + 1)):
        raise ValueError('global k-point entropy rows have gaps')
    bands = sorted(bands_by_k[k_indices[0]])
    if bands != list(range(bands[0], bands[-1] + 1)) or any(bands_by_k[k] != set(bands) for k in k_indices):
        raise ValueError('entropy rows do not form a complete rectangular k x full-band table')
    if len(keys) != len(k_indices) * len(bands):
        raise ValueError('entropy state table is incomplete')
    expected_nk = spec.get('expected_n_kpoints')
    expected_bands = spec.get('expected_band_indices')
    semantics = spec.get('band_index_semantics')
    if not isinstance(expected_nk, int) or expected_nk < 1:
        raise ValueError('expected_n_kpoints must be source-bound positive integer metadata')
    if not isinstance(expected_bands, list) or not expected_bands or \
            any(not isinstance(x, int) for x in expected_bands):
        raise ValueError('expected_band_indices must be the source-bound full band-index list')
    if expected_bands != list(range(len(expected_bands))):
        raise ValueError('expected band indices must be the exact zero-based global band range')
    if semantics != 'zero-based global spin-inclusive row of eigenvalues(:,ik)':
        raise ValueError('band-index spin-inclusive semantics are not explicitly bound')
    if len(k_indices) != expected_nk:
        raise ValueError(f'captured k count {len(k_indices)} != expected {expected_nk}')
    if k_indices != list(range(expected_nk)):
        raise ValueError('captured global k indices must equal the exact zero-based source range')
    if bands != expected_bands:
        raise ValueError(f'captured band range {bands} != expected {expected_bands}')

    def compute(prec: int) -> tuple[Decimal, Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = prec
            tlo = temp_k - qtk
            thi = temp_k + qtk
            s_center = Decimal(0)
            s_lo = Decimal(0)
            s_hi = Decimal(0)
            for state in states:
                if not all(isinstance(state.get(k), str) for k in
                           ('occupation_text', 'weight_text', 'eigenvalue_text')):
                    raise ValueError('eigenvalue, occupation and k-weight must preserve original serialized tokens')
                ft, wt = state['occupation_text'], state['weight_text']
                f, w = Decimal(ft.replace('D','E')), Decimal(wt.replace('D','E'))
                qf, qw = q_decimal(ft), q_decimal(wt)
                if f < 0 or f > 1 or w-qw < 0:
                    raise ValueError('invalid occupation or k-weight interval')
                flo, fhi = max(Decimal(0), f-qf), min(Decimal(1), f+qf)
                hlo, hhi = _h_interval_from_f(flo, fhi, prec)
                hcenter, _ = _h_interval_from_f(f, f, prec)
                s_center += w*hcenter
                s_lo += (w-qw)*hlo
                s_hi += (w+qw)*hhi
            kb = Decimal.from_float(6.3336814e-6)
            center = temp_k*kb*s_center
            lo = tlo*kb*s_lo
            hi = thi*kb*s_hi
            return center, lo, hi

    ts_center_100, ts_lo_100, ts_hi_100 = compute(precision)
    ts_center, ts_lo, ts_hi = compute(160)
    hp_width = max(abs(ts_center-ts_center_100), abs(ts_lo-ts_lo_100), abs(ts_hi-ts_hi_100))
    hp_ulp = Decimal(1).scaleb(ts_center.adjusted()-150) if ts_center else Decimal(1).scaleb(-150)
    hp_ulp *= Decimal(256 * len(states) + 256)
    ts_bound = max(abs(ts_center-ts_lo), abs(ts_hi-ts_center)) + hp_width + hp_ulp
    ts_text = str(spec['TS_text'])
    ts_logged = Decimal(ts_text.replace('D','E'))
    qts = q_decimal(ts_text)
    center_difference = abs(ts_logged-ts_center)
    total = qts + center_difference + ts_bound
    total_float = float(total)
    if Decimal.from_float(total_float) < total:
        total_float = math.nextafter(total_float, math.inf)
    return total_float, {
        'TS_high_precision_center_ry': str(ts_center),
        'TS_logged_token': ts_text,
        'TS_center_token_difference_ry': str(center_difference),
        'TS_input_interval_bound_ry': str(ts_bound),
        'TS_high_precision_enclosure_ry': str(hp_width+hp_ulp),
        'HP_precision_enclosure_operation_envelope': 256 * len(states) + 256,
        'B_TS_serialization_ry': str(qts),
        'B_TS_total_ry': str(total),
        'entropy_state_count': len(states),
        'entropy_global_k_count': len(k_indices),
        'entropy_full_band_count_per_k': len(bands),
        'expected_n_kpoints': expected_nk,
        'expected_band_indices': expected_bands,
        'band_index_semantics': semantics,
        'spin_convention': 'band_index is the global spin-inclusive eigenvalue row; no extra spin factor',
        'serialized_input_intervals': 'occupation, k-weight, and temperature K half-quantum intervals',
        'same_state_tokens_bound': ['fermi_level_text', 'temperature_ry_text', 'all eigenvalue_text tokens'],
        'production_entropy_reference': 'captured occupation rows are emitted by the same fermi_dirac_occupation kernel called by fermi_dirac_entropy, with the same k weights',
        'production_entropy_error_treatment': 'captured by |TS token - independent high-precision center|; no assumed libm ulp bound',
    }


def candidate(row: dict, mode: str, w: list[list[complex]], p: list[list[complex]],
              p_bounds: list[list[float]]) -> tuple[float, float, dict]:
    ts_text = row['TS_text']
    if not isinstance(ts_text, str):
        raise ValueError('TS must preserve its original serialized token')
    b_ts, ts_parts = entropy_roundoff_bound(row['entropy_inputs'] | {'TS_text': ts_text})
    ts = float(ts_text.replace('D', 'E'))
    tr = trace_wp(w, p)
    if mode == 'spectral':
        e_text = row['ETOT_spectral_text']
        energy = float(e_text.replace('D', 'E'))
        external = 0.0
        b_external = 0.0
    else:
        e_text = row['ETOT_local_text']
        energy = float(e_text.replace('D', 'E'))
        lam_text = row['lambda_runtime_text']
        lam = float(lam_text.replace('D', 'E'))
        q_lam = 0.5 * math.ulp(lam) if lam != 0.0 else 0.0
        external = lam * tr
        trace_scale = trace_abs_scale(w, p)
        trace_input = trace_input_bound(w, p_bounds)
        b_trace_raw = gamma(256) * trace_scale + trace_input
        b_trace = upward(b_trace_raw + gamma(64) * b_trace_raw)
        b_external_raw = abs(lam) * b_trace + q_lam * (abs(tr) + b_trace)
        b_external = upward(b_external_raw + gamma(3) * b_external_raw)
    if not all(isinstance(x, str) for x in (e_text, row.get('lambda_runtime_text', '0.0'))):
        raise ValueError('ETOT and lambda must preserve their original tokens')
    f = energy + external - ts
    q_energy = q_decimal_token(e_text)
    b_arithmetic = gamma(3) * (abs(energy) + abs(ts) + abs(external))
    b_f = upward(q_energy + b_ts + b_arithmetic + b_external +
                 gamma(4) * (q_energy + b_ts + b_arithmetic + b_external))
    if not all(math.isfinite(x) for x in (f, b_f, tr)):
        raise ValueError('nonfinite candidate or bound')
    return f, b_f, {'B_ETOT_serialization_ry': q_energy, **ts_parts,
                    'B_energy_arithmetic_ry': b_arithmetic,
                    'B_lambda_trace_ry': b_external,
                    'B_F_total_ry': b_f,
                    'Tr_WP': tr}


def h2_resolution(d: dict[float, float], ud: dict[float, float]) -> dict:
    """Approved h²-dominance gate with propagated h²/h⁴ coefficient bounds."""
    # Retain the adjacent-scale differences as diagnostics.  These are not
    # the approved resolution criterion; that is the L/H dominance test below.
    x = d[0.01] - d[0.005]
    y = d[0.005] - d[0.0025]
    ux = upward(ud[0.01] + ud[0.005] + gamma(2) * (ud[0.01] + ud[0.005]))
    uy = upward(ud[0.005] + ud[0.0025] + gamma(2) * (ud[0.005] + ud[0.0025]))
    x_arithmetic = upward(gamma(1) * (abs(d[0.01]) + abs(d[0.005])))
    y_arithmetic = upward(gamma(1) * (abs(d[0.005]) + abs(d[0.0025])))
    a = d[0.005] - d[0.0025]
    b = d[0.01] - d[0.005]
    ba = gamma(1) * (abs(d[0.005]) + abs(d[0.0025]))
    bb = gamma(1) * (abs(d[0.01]) + abs(d[0.005]))
    ua_raw = ud[0.005] + ud[0.0025] + ba
    ub_raw = ud[0.01] + ud[0.005] + bb
    ua = upward(ua_raw + gamma(2) * ua_raw)
    ub = upward(ub_raw + gamma(2) * ub_raw)
    l = (16.0*a - b) / 36.0
    h = (b - 4.0*a) / 180.0
    ul_raw = (16.0*ua + ub) / 36.0 + gamma(3)*(16.0*abs(a)+abs(b))/36.0
    uh_raw = (ub + 4.0*ua) / 180.0 + gamma(3)*(abs(b)+4.0*abs(a))/180.0
    ul = upward(ul_raw + gamma(2)*ul_raw)
    uh = upward(uh_raw + gamma(2)*uh_raw)
    leading_resolved, lhs, rhs = h2_dominance_compare(l, ul, h, uh)
    dominant = leading_resolved and lhs > rhs
    return {'X':x,'Y':y,'uX':ux,'uY':uy,
            'X_arithmetic_bound':x_arithmetic,'Y_arithmetic_bound':y_arithmetic,
            'A':a,'B':b,'uA':ua,'uB':ub,'L_h2':l,'H_h4':h,
            'uL_h2':ul,'uH_h4':uh,'L_minus_uL_exact':str(lhs),
            '16_times_H_plus_uH_exact':str(rhs),
            'dominance_margin_exact':str(lhs-rhs),
            'leading_term_resolved':leading_resolved,'resolved':dominant,
            'criterion':'|L|-uL > 16*(|H|+uH), where L=(16A-B)/36 and H=(B-4A)/180'}


def h2_dominance_compare(l: float, ul: float, h: float, uh: float) -> tuple[bool, Fraction, Fraction]:
    """Compare exact rational values of finite binary64 L/H estimates and bounds."""
    if not all(math.isfinite(x) for x in (l, ul, h, uh)) or ul < 0.0 or uh < 0.0:
        raise ValueError('nonfinite coefficient or invalid uncertainty bound')
    lf, ulf = Fraction.from_float(abs(l)), Fraction.from_float(ul)
    hf, uhf = Fraction.from_float(abs(h)), Fraction.from_float(uh)
    return lf > ulf, lf - ulf, 16 * (hf + uhf)


def evaluate(data: dict) -> dict:
    gamma_check = verify_gamma_bounds()
    w = matrix(data['W1'])
    if any(w[i][j] != complex(1.0 if i == j else 0.0, 0.0) for i in range(5) for j in range(5)):
        raise ValueError('W1 must be the exact 5x5 identity')
    expected = {(m, lam, t) for m in MODES for lam in LAMBDAS for t in ('standard', 'tight')}
    seen = {}
    for row in data['captures']:
        key = (row['mode'], float(row['lambda_ry']), row['threshold'])
        if key not in expected or key in seen:
            raise ValueError(f'unexpected/duplicate capture key: {key}')
        p = matrix(row['P_d'])
        pb = matrix_bounds(row['P_d_component_error_bound'])
        f, bf, parts = candidate(row, key[0], w, p, pb)
        seen[key] = {'P': p, 'P_bound': pb, 'F': f, 'B_F': bf, 'B_F_parts': parts}
    if set(seen) != expected:
        raise ValueError(f'incomplete capture set; missing {sorted(expected-set(seen), key=str)}')

    results = {}
    for mode in MODES:
        ps = seen[(mode, 0.0, 'standard')]['P']
        pt = seen[(mode, 0.0, 'tight')]['P']
        psb = seen[(mode, 0.0, 'standard')]['P_bound']
        ptb = seen[(mode, 0.0, 'tight')]['P_bound']
        g0 = trace_wp(w, ps)
        dp = [[pt[i][j] - ps[i][j] for j in range(5)] for i in range(5)]
        ug_state = abs(trace_wp(w, dp))
        fp_scale = trace_abs_scale(w, ps) + trace_abs_scale(w, pt)
        input_fp = trace_input_bound(w, psb) + trace_input_bound(w, ptb)
        ug_fp_raw = gamma(256) * fp_scale + input_fp
        ug_fp = upward(ug_fp_raw + gamma(128)*ug_fp_raw)
        ug0 = upward(ug_state + ug_fp + gamma(2)*(ug_state + ug_fp))
        d, ud, uf, d_analysis_fp = {}, {}, {}, {}
        for h in HS:
            fp, fm = seen[(mode, h, 'standard')]['F'], seen[(mode, -h, 'standard')]['F']
            d[h] = (fp - fm) / (2.0 * h)
            h_exact = Decimal(str(h))
            q_h = abs(Decimal.from_float(h) - h_exact)
            d_arithmetic = gamma(2) * (abs(fp) + abs(fm)) / (2.0 * h) + \
                abs(d[h]) * float(q_h / h_exact)
            d_analysis_fp[h] = upward(d_arithmetic)
            for sign, lam in (('+', h), ('-', -h)):
                st, tt = seen[(mode, lam, 'standard')], seen[(mode, lam, 'tight')]
                uf_raw = abs(tt['F'] - st['F']) + tt['B_F'] + st['B_F'] + \
                    gamma(1) * (abs(tt['F']) + abs(st['F']))
                uf[f'{sign}{h:g}'] = upward(uf_raw + gamma(3)*(abs(uf_raw) + tt['B_F'] + st['B_F']))
            u_pair = uf[f'+{h:g}'] + uf[f'-{h:g}']
            u_div = u_pair / (2.0 * h)
            u_pair_fp = gamma(2) * u_pair / (2.0 * h)
            ud[h] = upward(u_div + d_analysis_fp[h] + u_pair_fp)
        rf = (4.0 * d[0.0025] - d[0.005]) / 3.0
        rc = (4.0 * d[0.005] - d[0.01]) / 3.0
        rfp = gamma(3) * (4.0 * abs(d[0.0025]) + abs(d[0.005])) / 3.0
        rcp = gamma(3) * (4.0 * abs(d[0.005]) + abs(d[0.01])) / 3.0
        urf_raw = (4.0 * ud[0.0025] + ud[0.005]) / 3.0
        urc_raw = (4.0 * ud[0.005] + ud[0.01]) / 3.0
        urf = upward(urf_raw + rfp + gamma(2) * urf_raw)
        urc = upward(urc_raw + rcp + gamma(2) * urc_raw)
        btrunc_raw = abs(rf - rc) / 15.0 + (urf + urc) / 15.0
        btrunc_fp = gamma(4) * (abs(rf) + abs(rc) + urf + urc) / 15.0
        btrunc = upward(btrunc_raw + btrunc_fp)
        nraw = urf + btrunc + ug0
        n_unc = upward(nraw + gamma(2) * (urf + btrunc + ug0))
        scaling = h2_resolution(d, ud)
        resolved = scaling['resolved']
        bias_decimal = Decimal.from_float(rf) - Decimal.from_float(g0)
        bias = float(bias_decimal)
        bias_pass = abs(bias_decimal) <= Decimal.from_float(n_unc)
        results[mode] = {
            'g0_dimensionless': g0, 'ug0_state_sensitivity_dimensionless': ug_state,
            'ug0_trace_and_input_fp_dimensionless': ug_fp, 'ug0_total_dimensionless': ug0,
            'g0_standard_P': [[[z.real, z.imag] for z in row] for row in ps],
            'g0_tight_P': [[[z.real, z.imag] for z in row] for row in pt],
            'g0_standard_P_abs_error_bound': psb, 'g0_tight_P_abs_error_bound': ptb,
            'B_F_components': {str(l): {t: seen[(mode,l,t)]['B_F_parts'] for t in ('standard','tight')} for l in LAMBDAS},
            'F_standard_ry': {str(l): seen[(mode,l,'standard')]['F'] for l in LAMBDAS},
            'D_dimensionless': {str(h): d[h] for h in HS},
            'D_analysis_roundoff_bound_dimensionless': {str(h): d_analysis_fp[h] for h in HS},
            'uF_ry': uf, 'uD_dimensionless': {str(h): ud[h] for h in HS},
            'R_fine_dimensionless': rf, 'R_coarse_dimensionless': rc,
            'R_fine_analysis_roundoff_bound_dimensionless': rfp,
            'R_coarse_analysis_roundoff_bound_dimensionless': rcp,
            'uR_fine_dimensionless': urf, 'uR_coarse_dimensionless': urc,
            'B_trunc_dimensionless': btrunc, 'B_trunc_analysis_roundoff_bound': btrunc_fp,
            'N_unc_dimensionless': n_unc,
            'B_W_dimensionless': bias,
            'scaling': {**scaling, 'abs_X_gt_uX':abs(scaling['X'])-scaling['X_arithmetic_bound']>scaling['uX'],
                        'abs_Y_gt_uY':abs(scaling['Y'])-scaling['Y_arithmetic_bound']>scaling['uY']},
            'identity_pass': bool(resolved and bias_pass),
            'identity_classification': ('PASS' if resolved and bias_pass else
                                        'FAIL' if resolved else 'INCONCLUSIVE_STOP')}
    return {'schema':'rslmto-a4-m21-w1-richardson-analysis-v5',
            'input_sha256':data.get('_input_sha256'),'criteria_source':data.get('criteria_source'),
            'gamma_bound_verification': gamma_check,
            'serialized_format':'Fortran ES/e24.16; per-token half-quantum derived from printed exponent/precision',
            'lambda_grid_ry':list(LAMBDAS),'status':'ANALYSIS_ONLY_NO_SCFS_RUN_BY_THIS_TOOL','modes':results}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise SystemExit(f'refusing to overwrite {args.output}')
    data = json.loads(args.input.read_text())
    data['_input_sha256'] = sha(args.input)
    result = evaluate(data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),
                      'classifications':{m:v['identity_classification'] for m,v in result['modes'].items()}},sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
