#!/usr/bin/env python3
"""Append the six preregistered RS A2 shifts to its captured A1 series.

The base plus eight A1 responses are read from the immutable accepted RS
capture. The six new responses use the existing production calculate_pl /
NEWRHO fixture driver on the captured frozen radial operands.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import accepted_m03_replay as replay  # noqa: E402
import m0_a11_scaling_oracle as a1  # noqa: E402
import m0_a22_scaling_oracle as a22  # noqa: E402

CAPTURE = Path('/private/tmp/rslmto-localmom-a2-rs-floor3e8')
DRIVER = Path('/private/tmp/rslmto-localmom-a2-u2-d91-bin/accepted_m03_replay_driver')
SHIFTS = (-0.12, -0.08, -0.04, 0.04, 0.08, 0.12)
ALL_SHIFTS = (-0.12, -0.08, -0.04, -0.02, -0.01, -0.005, -0.0025,
              0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.12)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def fit_rs(series: dict[float, dict[str, float]], q0: float, d0: float) -> dict:
    a1.V_MAGNITUDES = a22.V_MAGNITUDES
    a1.SIGNED_V = a22.SIGNED_V
    fit = a1._series_report(series, q0, d0)
    branches = a22._branch_signs(fit['r_ry_by_signed_v'])
    signs_ok = all(x['constant_nonzero_sign_within_branch'] for x in branches.values())
    zero = any(x['zero_is_unassessed_without_pseudocount'] for x in branches.values())
    order_ok = fit['r_loglog_ols_order_all_eight'] is not None and fit['r_loglog_ols_order_all_eight'] >= 2.0
    plateau_ok = fit['plateau_intercept_pass']
    monotone_ok = all(fit['r_abs_strictly_decreases_as_v_halves_by_sign'].values())
    fit['a2_branch_sign_diagnostics'] = branches
    fit['a2_shape_gates'] = {
        'no_sign_reversal_within_either_signed_v_branch': signs_ok,
        'exact_zero_residual_unassessed_without_pseudocount': zero,
        'all_point_loglog_order_at_least_two': order_ok,
        'cubic_constant_within_10D0': plateau_ok,
        'absolute_residual_strictly_decreases_toward_zero_on_both_branches': monotone_ok,
        'pass': signs_ok and not zero and order_ok and plateau_ok and monotone_ok,
    }
    fit['r_loglog_ols_order_all_fourteen'] = fit.pop('r_loglog_ols_order_all_eight')
    fit['r_loglog_intercept_all_fourteen'] = fit.pop('r_loglog_intercept_all_eight')
    fit['status'] = 'UNASSESSED_LOG_ZERO' if zero else ('PASS' if fit['a2_shape_gates']['pass'] else 'FAIL_A2_2')
    return fit


def prepare_captured_raw(rows: list[dict[str, str]], iteration: int) -> tuple[dict, dict]:
    """Rehydrate the complete accepted A1 RAW_PL input snapshot, not LIVE_STATE."""
    stage = 'RAW_PL'
    for key in ('m0_replay_rofi', 'm0_replay_v_input', 'm0_replay_rho_input',
                'm0_replay_ec_seed', 'm0_replay_ev_seed', 'm0_replay_pl_input',
                'm0_replay_ql_input', 'm0_replay_pl_return', 'm0_replay_control'):
        if any(int(r['iteration']) != iteration for r in rows if r['kind'] == key and r['stage'] == stage):
            raise ValueError(f'{key}/{stage} contains a different selected iteration')
    nr = len([r for r in rows if r['kind'] == 'm0_replay_rofi' and r['stage'] == stage])
    ncore = len([r for r in rows if r['kind'] == 'm0_replay_ec_seed' and r['stage'] == stage])
    nval = len([r for r in rows if r['kind'] == 'm0_replay_ev_seed' and r['stage'] == stage])
    if (nr, ncore, nval) != (495, 10, 6):
        raise ValueError(f'unsupported captured RAW_PL dimensions {(nr, ncore, nval)}')
    d = {'rofi': replay.ordered(rows, 'm0_replay_rofi', stage, 'scalar', nr),
         'ec_seed': replay.ordered(rows, 'm0_replay_ec_seed', stage, 'scalar', ncore),
         'ev_seed': replay.ordered(rows, 'm0_replay_ev_seed', stage, 'scalar', nval),
         'params': rows}
    d['v'] = [[float(replay.one(rows, 'm0_replay_v_input', stage, spin=s+1, row=i)['scalar'])
               for s in (0, 1)] for i in range(nr)]
    d['rho_in'] = [[float(replay.one(rows, 'm0_replay_rho_input', stage, spin=s+1, row=i)['scalar'])
                    for s in (0, 1)] for i in range(nr)]
    arrays = {name: np.zeros((3, 2)) for name in ('gc', 'c', 'qp', 'sd', 'pl_seed', 'pl_live')}
    ql = np.zeros((3, 3, 2))
    for spin in (1, 2):
        for ell in range(3):
            for order, name in enumerate(('gc', 'c', 'qp', 'sd', 'pl_seed')):
                arrays[name][ell, spin-1] = float(replay.one(rows, 'm0_replay_pl_input', stage,
                                                              spin=spin, ell=ell, order=order)['scalar'])
            arrays['pl_live'][ell, spin-1] = float(replay.one(rows, 'm0_replay_pl_return', stage,
                                                               spin=spin, ell=ell, order=0)['scalar'])
            for moment in range(3):
                ql[moment, ell, spin-1] = float(replay.one(rows, 'm0_replay_ql_input', stage,
                                                            spin=spin, ell=ell, order=moment)['scalar'])
    d.update(arrays)
    control = {int(r['order']): r for r in rows if r['kind'] == 'm0_replay_control' and r['stage'] == stage}
    if set(control) != set(range(6)):
        raise ValueError('RAW_PL control snapshot is incomplete')
    # The production trace records exact A/Z, B/tolerance, NR/EC/EV counts,
    # VHRMAX and B_fsm in the order-indexed controls. Only f_core/txc are
    # static integer runtime settings not repeated in that tuple.
    params = {
        'radial_A': control[0]['energy'], 'radial_Z': control[0]['scalar'],
        'radial_B': control[1]['energy'], 'radial_tolerance': control[1]['scalar'],
        'radial_NR': control[2]['energy'], 'radial_f_core': '0', 'radial_txc': '1',
        'radial_VHRMAX': control[4]['energy'], 'radial_B_fsm': control[5]['energy'],
        'radial_hyperfine': '0',
    }
    d['raw_ql'] = ql
    # A1 captured this complete matrix. All RAW inputs remain exact; each old
    # shift changes only the d-shell center tuple and its separately recorded
    # shift scalar. This proves the proposed new cases use the same time level.
    frozen_kinds = ('m0_replay_rofi', 'm0_replay_v_input', 'm0_replay_rho_input',
                    'm0_replay_ec_seed', 'm0_replay_ev_seed', 'm0_replay_ql_input',
                    'm0_replay_radial_parameter')
    key_fields = ('site', 'type', 'l', 'spin', 'row', 'col', 'band', 'order')
    value_fields = ('energy', 'scalar', 're', 'im')
    for kind in frozen_kinds:
        raw = {tuple(r[k] for k in key_fields): tuple(r[k].strip() for k in value_fields)
               for r in rows if r['kind'] == kind and r['stage'] == 'RAW_PL'}
        for idx in range(1, 9):
            other = {tuple(r[k] for k in key_fields): tuple(r[k].strip() for k in value_fields)
                     for r in rows if r['kind'] == kind and r['stage'] == f'A1_SHIFT_{idx:02d}'}
            if raw != other:
                raise ValueError(f'{kind} frozen inputs differ between RAW_PL and A1_SHIFT_{idx:02d}')
    raw_pl = {tuple(r[k] for k in key_fields): tuple(r[k].strip() for k in value_fields)
              for r in rows if r['kind'] == 'm0_replay_pl_input' and r['stage'] == 'RAW_PL'}
    for idx in range(1, 9):
        stage_name = f'A1_SHIFT_{idx:02d}'
        other = {tuple(r[k] for k in key_fields): tuple(r[k].strip() for k in value_fields)
                 for r in rows if r['kind'] == 'm0_replay_pl_input' and r['stage'] == stage_name}
        shifts = [r for r in rows if r['kind'] == 'm0_replay_control' and r['stage'] == stage_name and int(r['order']) == 3]
        if len(shifts) != 1 or set(raw_pl) != set(other):
            raise ValueError(f'{stage_name} PL snapshot keys or center shift are incomplete')
        shift = float(shifts[0]['energy'])
        for key in raw_pl:
            expected = list(raw_pl[key])
            observed = other[key]
            if int(key[2]) == 2 and int(key[3]) in (1, 2) and int(key[7]) == 0:
                expected[1] = f'{float(expected[1]) - shift:.16E}'
                if abs(float(observed[1]) - float(expected[1])) > 2e-15:
                    raise ValueError(f'{stage_name} d-shell center is not the registered shift')
                observed_rest = (observed[0], observed[2], observed[3])
                if observed_rest != (expected[0], expected[2], expected[3]):
                    raise ValueError(f'{stage_name} changed non-center d-shell PL inputs')
            elif tuple(expected) != observed:
                raise ValueError(f'{stage_name} changed a non-center PL input at {key}')
        raw_controls = {int(r['order']): tuple(r[k].strip() for k in value_fields)
                        for r in rows if r['kind'] == 'm0_replay_control' and r['stage'] == 'RAW_PL'}
        shift_controls = {int(r['order']): tuple(r[k].strip() for k in value_fields)
                          for r in rows if r['kind'] == 'm0_replay_control' and r['stage'] == stage_name}
        if set(raw_controls) != set(range(6)) or set(shift_controls) != set(range(6)):
            raise ValueError(f'{stage_name} control snapshot is incomplete')
        if any(raw_controls[o] != shift_controls[o] for o in range(6) if o != 3):
            raise ValueError(f'{stage_name} changed a non-shift control')
        if raw_controls[3][1:] != shift_controls[3][1:] or abs(float(shift_controls[3][0]) - shift) > 2e-15:
            raise ValueError(f'{stage_name} shift control does not match d-shell center input')
    return d, params


def text_input_from_captured_raw(d: dict, params: dict, nr: int, ncore: int, nval: int) -> str:
    a, b, z = float(params['radial_A']), float(params['radial_B']), float(params['radial_Z'])
    tol, vhrmax, bfsm = float(params['radial_tolerance']), float(params['radial_VHRMAX']), float(params['radial_B_fsm'])
    fcore, txc = int(float(params['radial_f_core'])), int(float(params['radial_txc']))
    if float(params['radial_hyperfine']) != 0.0:
        raise ValueError('captured hyperfine=true is unsupported by this replay seam')
    lines = [f"2 {nr} {ncore} {nval} {len(d['cases'])} {fcore} {txc}",
             f"{z:.17g} {a:.17g} {b:.17g} {tol:.17g} {vhrmax:.17g} {bfsm:.17g}",
             ' '.join(d['rofi'])]
    lines += [' '.join(f'{x:.17g}' for x in (*d['v'][i], *d['rho_in'][i])) for i in range(nr)]
    lines += [' '.join(d['ec_seed']), ' '.join(d['ev_seed'])]
    for spin in range(2):
        for ell in range(3):
            lines.append(' '.join(f'{x:.17g}' for x in (d['gc'][ell, spin], d['c'][ell, spin],
                d['qp'][ell, spin], d['sd'][ell, spin], d['pl_seed'][ell, spin], d['pl_live'][ell, spin])))
    for _, shift, calculate, ql in d['cases']:
        lines.append(f"{shift[0]:.17g} {shift[1]:.17g} {'.true.' if calculate else '.false.'}")
        for spin in range(2):
            for ell in range(3):
                for moment in range(3):
                    lines.append(f'{ql[moment, ell, spin]:.17g}')
    return '\n'.join(lines) + '\n'


def raw_anchor_check(rows: list[dict[str, str]], case: dict) -> dict:
    """Check the zero-shift replay against the same captured RAW_PL time level."""
    checks = {}
    refs = [r for r in rows if r['kind'] == 'm0_replay_pl_return' and r['stage'] == 'RAW_PL']
    pl_ref = {(int(r['l']), int(r['spin'])): float(r['scalar']) for r in refs}
    pl_error = max(abs(pl_ref[k] - case['pl'][k]) for k in pl_ref)
    checks['PL'] = {'max_abs': pl_error, 'tolerance': 1e-10, 'pass': pl_error <= 1e-10}
    refs = [r for r in rows if r['kind'] == 'm0_replay_rho' and r['stage'] == 'RAW_PL']
    rho_ref = {(int(r['row']), int(r['spin'])): float(r['scalar']) for r in refs}
    rho_error = max(abs(rho_ref[k] - case['rho'][k]) for k in rho_ref)
    checks['rho'] = {'max_abs': rho_error, 'tolerance': 1e-10, 'pass': rho_error <= 1e-10}
    refs = [r for r in rows if r['kind'] == 'm0_replay_energy' and r['stage'] == 'RAW_PL']
    ec_ref = {int(r['spin']): float(r['energy']) for r in refs}
    ev_ref = {int(r['spin']): float(r['scalar']) for r in refs}
    ec_error = max(abs(ec_ref[s] - case['sumec'][s]) for s in ec_ref)
    ev_error = max(abs(ev_ref[s] - case['sumev'][s]) for s in ev_ref)
    checks['SUMEC'] = {'max_abs_ry': ec_error, 'tolerance_ry': 1e-8, 'pass': ec_error <= 1e-8}
    checks['SUMEV'] = {'max_abs_ry': ev_error, 'tolerance_ry': 1e-8, 'pass': ev_error <= 1e-8}
    checks['pass'] = all(v['pass'] for v in checks.values())
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prereg-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--anchor-only', action='store_true', help='run and report the RAW_PL zero anchor, then stop before shifts')
    args = parser.parse_args()
    repo = HERE.parents[2]
    prereg_path = repo / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_preregistration_v3.json'
    prereg = json.loads(prereg_path.read_text())
    check(sha(prereg_path) == args.prereg_sha256, 'RS scalar prereg hash mismatch')
    check(prereg['status'] == 'PREREGISTERED_BEFORE_SIX_NEW_RS_SHIFT_RESPONSES', 'RS scalar prereg not locked')
    cap_manifest_path = CAPTURE / 'manifest.json'
    cap_manifest = json.loads(cap_manifest_path.read_text())
    cap_csv = CAPTURE / 'records.csv'
    driver_source = HERE / 'accepted_m03_replay_driver.f90'
    adapter_source = HERE / 'accepted_m03_replay.py'
    a1_source = HERE / 'm0_a11_scaling_oracle.py'
    a22_source = HERE / 'm0_a22_scaling_oracle.py'
    scan_source = HERE / 'a3_rs_scalar_scan.py'
    check(sha(cap_csv) == prereg['accepted_no_u_state']['csv_sha256'], 'accepted RS CSV changed')
    check(sha(cap_manifest_path) == prereg['accepted_no_u_state']['manifest_sha256'], 'accepted RS manifest changed')
    txc_log = CAPTURE / 'testrun.log'
    check(sha(txc_log) == prereg['accepted_no_u_state']['runtime_log_sha256'], 'accepted RS runtime log changed')
    if 'XC selector: TXC=1' not in txc_log.read_text(errors='replace'):
        raise ValueError('accepted RS runtime log does not bind TXC=1')
    check(sha(DRIVER) == prereg['production_replay_driver']['binary_sha256'], 'production replay driver changed')
    check(sha(driver_source) == prereg['production_replay_driver']['source_sha256'], 'production replay driver source changed')
    bands_source = Path(prereg['production_replay_driver']['bands_source_path'])
    check(sha(bands_source) == prereg['production_replay_driver']['bands_source_sha256'], 'immutable calculate_pl source snapshot changed')
    for path, key in ((adapter_source, 'adapter_sha256'), (a1_source, 'a1_oracle_sha256'),
                       (a22_source, 'a22_oracle_sha256'), (scan_source, 'scan_script_sha256')):
        check(sha(path) == prereg['production_replay_driver'][key], f'{key} changed')

    with cap_csv.open(newline='') as f:
        rows = list(csv.DictReader(f))
    replay._BY_KIND_STAGE = {}
    for row in rows:
        replay._BY_KIND_STAGE.setdefault((row['kind'], row['stage']), []).append(row)
    iteration = int(cap_manifest['selected_iteration'])
    check(iteration == 2, 'RS no-U accepted iteration mismatch')
    d, raw_params = prepare_captured_raw(rows, iteration)
    nr = int(float(raw_params['radial_NR']))
    ncore = len(d['ec_seed'])
    nval = len(d['ev_seed'])

    # Verify frozen q0/q2 and raw SUMEV against the existing A1 capture before
    # launching the six additional response cases.
    hist = {}
    for row in rows:
        if row['kind'] != 'm0_a11_case' or int(row['iteration']) != iteration or int(row['spin']) not in (1, 2):
            continue
        s, v = int(row['spin']), float(row['weight'])
        hist.setdefault(s, {})[v] = {
            'raw_sumev': float(row['energy']), 'case_sumev': float(row['scalar']),
            'q0': float(row['re']), 'q2': float(row['im'])}
    q0_by_spin, q2_by_spin = {}, {}
    for s in (1, 2):
        check(set(hist[s]) == {0.0, *a1.SIGNED_V}, f'missing or duplicate historical RS A1 shifts for spin {s}')
        ql = d['raw_ql']
        ql_q0, ql_q2 = float(ql[0, 2, s-1]), float(ql[2, 2, s-1])
        q0_by_spin[s], q2_by_spin[s] = ql_q0, ql_q2
        check(ql_q0 == hist[s][0.0]['q0'], f'RS A1 q0 differs from captured RAW_PL QL spin {s}')
        check(abs(ql_q2-hist[s][0.0]['q2']) <= 1e-12, f'RS A1 q2 differs from accepted frozen QL spin {s}')
        check(all(x['raw_sumev'] == hist[s][0.0]['raw_sumev'] and x['q0'] == ql_q0 and x['q2'] == hist[s][0.0]['q2']
                  for x in hist[s].values()), f'RS A1 responses are not frozen-input cases spin {s}')

    raw_ql = d['raw_ql'].copy()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    import subprocess
    def run_fixture(prefix: str, cases: list[tuple]) -> tuple[dict[int, dict], dict]:
        d['cases'] = cases
        payload = text_input_from_captured_raw(d, raw_params, nr, ncore, nval)
        input_path = args.output_dir / f'{prefix}.input.txt'
        output_path = args.output_dir / f'{prefix}.output.csv'
        log_path = args.output_dir / f'{prefix}.stdout.log'
        input_path.write_text(payload)
        proc = subprocess.run([str(DRIVER), str(input_path), str(output_path)], capture_output=True, text=True)
        log_path.write_text(proc.stdout + proc.stderr)
        check(proc.returncode == 0, f'production fixture driver returned {proc.returncode}; see {log_path}')
        parsed = replay.output_cases(output_path.read_bytes())
        return parsed, {'input_path': str(input_path), 'input_sha256': sha(input_path),
                        'output_path': str(output_path), 'output_sha256': sha(output_path),
                        'log_path': str(log_path), 'log_sha256': sha(log_path)}
    anchor_cases = [('RS_RAW_ANCHOR', np.zeros(2), True, raw_ql.copy())]
    anchor_results, anchor_artifacts = run_fixture('raw_anchor', anchor_cases)
    check(set(anchor_results) == {1}, 'RAW_PL anchor driver returned unexpected cases')
    baseline_check = raw_anchor_check(rows, anchor_results[1])
    check(baseline_check['pass'], 'zero-shift production replay did not reproduce captured RAW_PL state')
    if args.anchor_only:
        anchor_report = {
            'status': 'RAW_PL_ZERO_ANCHOR_PASS_NO_SHIFTS_RUN',
            'prereg_path': str(prereg_path), 'prereg_sha256': sha(prereg_path),
            'anchor_check': baseline_check, 'anchor_artifacts': anchor_artifacts,
            'radial_f_core': raw_params['radial_f_core'], 'radial_txc': raw_params['radial_txc'],
        }
        (args.output_dir / 'anchor_report.json').write_text(json.dumps(anchor_report, indent=2, sort_keys=True) + '\n')
        print(json.dumps(anchor_report, indent=2, sort_keys=True))
        return
    # Responses remain unexecuted until the independent zero-shift production
    # replay has passed every captured RAW_PL PL/rho/SUMEC/SUMEV guard.
    response_cases = [(f'RS_A2_SHIFT_{i+1:02d}', np.array([v, v]), True, raw_ql.copy())
                      for i, v in enumerate(SHIFTS)]
    new_cases, response_artifacts = run_fixture('six_shift_responses', response_cases)
    check(set(new_cases) == set(range(1, 7)), 'A2 response driver did not return exactly six shift cases')

    d0 = float(prereg['accepted_no_u_state']['D0_RS_ry'])
    spin_reports = {}
    for spin in (1, 2):
        q0 = q0_by_spin[spin]
        series = {0.0: {'delta_sumev_ry': 0.0, 'r_ry': 0.0}}
        for v in a1.SIGNED_V:
            point = hist[spin][v]
            delta = point['raw_sumev'] - point['case_sumev']
            series[v] = {'delta_sumev_ry': delta, 'r_ry': delta-v*q0}
        new_values = {}
        for i, v in enumerate(SHIFTS, 1):
            delta = hist[spin][0.0]['raw_sumev'] - new_cases[i]['sumev'][spin]
            new_values[str(v)] = {'delta_sumev_ry': delta, 'r_ry': delta-v*q0}
            series[v] = new_values[str(v)]
        fit = fit_rs(series, q0, d0)
        spin_reports[str(spin)] = {'q0': q0, 'raw_sumev_ry': hist[spin][0.0]['raw_sumev'],
                                    'existing_a1_values': {str(v): series[v] for v in a1.SIGNED_V},
                                    'new_a2_values': new_values, 'fit': fit}

    report = {
        'schema':'rslmto-a3-rs-a22-scalar-fit-response-v1',
        'status':'MEASURED_RS_A22_SCALAR_RESPONSE',
        'prereg_path':str(prereg_path),'prereg_sha256':sha(prereg_path),
        'capture_csv_sha256':sha(cap_csv),'capture_manifest_sha256':sha(cap_manifest_path),
        'runtime_log_sha256':sha(txc_log),'runtime_xc_selector':'TXC=1 (Barth-Hedin), read from bound accepted-run log',
        'driver_source_sha256':sha(driver_source),'driver_binary_sha256':sha(DRIVER),
        'immutable_bands_source_path':str(bands_source),
        'immutable_bands_source_sha256':sha(bands_source),
        'live_checkout_bands_source_sha256':sha(repo/'source/bands.f90'),
        'anchor_artifacts':anchor_artifacts,
        'six_shift_response_artifacts':response_artifacts,
        'frozen_state_iteration':iteration,'frozen_array_checks':'Rehydrated complete captured RAW_PL stage arrays. Exact stage provenance checks verify rofi, v/rho, EC/EV seeds, QL, radial parameters, and all non-center PL inputs are byte-value identical across RAW_PL and each A1 shift; only d-shell center and its recorded shift differ. q0 and q2 are exact to captured A1 raw rows.',
        'anchor_reference':'Same captured M0 RAW_PL stage outputs: m0_replay_pl_return, m0_replay_rho and m0_replay_energy; generic LIVE_STATE radial_replay operands are not used for RAW_PL.',
        'raw_anchor_reproduction':baseline_check,
        'D0_RS_ry':d0,'u_RS_policy':'reported separately; not added to allowance',
        'E_num_ry':float(prereg['locked_models_and_checks']['E_num_ry']),
        'per_spin':spin_reports,
    }
    out = args.output_dir / 'report.json'
    out.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'report':str(out),'report_sha256':sha(out),'fixture_output_sha256':sha(output_path)},sort_keys=True))


if __name__ == '__main__':
    main()
