#!/usr/bin/env python3
"""Run only the pre-registered A3 sign and spin-provenance controls.

The frozen d91 capture, A2.2 fit, A2 LIVE baseline, and existing production
fixture driver are inputs. Outputs are written to a new private scratch folder;
no A2 evidence is changed.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import accepted_m03_replay as a2  # noqa: E402

CAPTURE = Path('/private/tmp/rslmto-localmom-a2-ks-u2-capture-d91')
DRIVER = Path('/private/tmp/rslmto-localmom-a2-u2-d91-bin/accepted_m03_replay_driver')
FIT = Path('/private/tmp/a2-accepted-iter23-a22-fit.json')
PREREG = Path(__file__).resolve().parents[3] / 'docs/validation/fixed_ldau/evidence/a3_m03_control_preregistration.json'
PREREG_SHA = '1f26829dfebffd2586ac891100f485bb3eda3ec7f6655fd100d56c6267dc177a'
FIT_SHA = '0ad6a7703aa26eb132602e7a43773d7619f70651c6b90817e4f945ff259cb993'
SOURCE_SHA = 'd91a468a4493da21710fc7f60e745b5f1bef8f20d5d667ae939c31d5ba3061ef'
DRIVER_SHA = '0c1c3c12e8501f4642af87168c611d261fb5a34cc4747aaef9c8420c2c6449c4'
A2_CASES = CAPTURE / 'm03_cases.json'
SCRATCH = Path('/private/tmp/a3-m03-sign-spin-response-v3')


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if SCRATCH.exists():
        raise SystemExit(f'refusing to overwrite existing output directory {SCRATCH}')
    prereg = json.loads(PREREG.read_text())
    if digest(PREREG) != PREREG_SHA or prereg.get('schema') != 'rslmto-a3-ks-m03-control-prereg-v1':
        raise ValueError('reviewed A3 preregistration hash/schema mismatch')
    if prereg.get('status') != 'PREDICTIONS_PREREGISTERED; SIGN_AND_SPIN_CONTROLS_RELEASED; ITERATION_CONTROL_APPLICABILITY_PENDING':
        raise ValueError('A3 control preregistration is not locked before responses')
    if digest(FIT) != FIT_SHA or digest(DRIVER) != DRIVER_SHA:
        raise ValueError('locked fit or production fixture driver hash mismatch')
    source = Path(__file__).resolve().parents[3] / 'source/self.f90'
    if digest(source) != SOURCE_SHA:
        raise ValueError('working source differs from frozen d91 source')
    rows, manifest, iteration, nr, ncore, nval, csv_path, manifest_path = a2.export_rows(CAPTURE)
    if iteration != 32 or digest(csv_path) != prereg['accepted_capture']['csv_sha256']:
        raise ValueError('accepted iteration-32 capture does not match A3 preregistration')
    if digest(manifest_path) != prereg['accepted_capture']['manifest_sha256']:
        raise ValueError('accepted capture manifest hash mismatch')
    if manifest.get('source_self_sha256') != SOURCE_SHA:
        raise ValueError('accepted capture source mismatch')
    if digest(FIT) != prereg['locked_A2_fit_sha256']:
        raise ValueError('accepted capture is not paired with the locked A2.2 fit')
    baseline = json.loads((CAPTURE / 'live_baseline.json').read_text())
    if baseline.get('status') != 'PASS_LIVE_BASELINE':
        raise ValueError('accepted capture lacks its passing LIVE_STATE baseline')
    prior_cases = json.loads(A2_CASES.read_text())
    if prior_cases.get('status') != 'MEASURED_REPLAY_DIAGNOSTIC':
        raise ValueError('accepted A2 RAW/LOCAL case report missing')
    if digest(Path(prior_cases['case_output_path'])) != prior_cases['case_output_sha256']:
        raise ValueError('accepted A2 case CSV hash mismatch')

    d = a2.prepare(rows, nr, ncore, nval, 'cases', manifest)
    target = d['target']
    raw_ql = d['cases'][0][3].copy()
    local_ql = d['cases'][1][3].copy()
    sign = prereg['sign_control']
    sign_ql = local_ql.copy()  # preregistered: center only; correct local q2 retained.
    sign_shift = np.array([float(sign['wrong_shift_v_ry'][str(s)]) for s in (1, 2)])
    spin = prereg['spin_provenance_control']
    spin_ql = raw_ql.copy()
    spin_shift = np.zeros(2)
    for s in (1, 2):
        item = spin['per_spin'][str(s)]
        spin_shift[s - 1] = float(item['wrong_shift_ry'])
        spin_ql[2, 2, s - 1] = float(item['q2_local_wrong_ry2'])

    d['cases'] = [
        ('RAW_PL_A3_ANCHOR', np.zeros(2), True, raw_ql),
        ('SIGN_CENTER_WRONG', sign_shift, True, sign_ql),
        ('OPPOSITE_SPIN_VIN', spin_shift, True, spin_ql),
    ]
    payload = a2.text_input(d, rows, nr, ncore, nval).encode()
    SCRATCH.mkdir(parents=True)
    input_path = SCRATCH / 'fixture_input.txt'
    output_path = SCRATCH / 'fixture_output.csv'
    stdout_path = SCRATCH / 'driver.stdout.log'
    input_path.write_bytes(payload)
    completed = subprocess.run([str(DRIVER.resolve()), str(input_path), str(output_path)],
                               capture_output=True, text=True)
    stdout_path.write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise RuntimeError(f'production frozen-state driver failed; see {stdout_path}')
    case_rows = a2.output_cases(output_path.read_bytes())
    if set(case_rows) != {1, 2, 3}:
        raise ValueError(f'unexpected response case indices {sorted(case_rows)}')

    raw, sign_case, spin_case = (case_rows[i] for i in (1, 2, 3))
    true_t = {s: float(np.trace(target['vin'][s - 1] @ target['p'][s - 1]).real) for s in (1, 2)}
    wrong_sign_t = {s: -true_t[s] for s in (1, 2)}
    wrong_spin_t = {s: float(np.trace(target['vin'][2 - s] @ target['p'][s - 1]).real) for s in (1, 2)}
    sign_delta_s = {s: raw['sumev'][s] - sign_case['sumev'][s] for s in (1, 2)}
    spin_delta_s = {s: raw['sumev'][s] - spin_case['sumev'][s] for s in (1, 2)}
    radial_a = float(a2.one(rows, 'radial_A', 'RADIAL_REPLAY')['scalar'])
    radial_b = float(a2.one(rows, 'radial_B', 'RADIAL_REPLAY')['scalar'])

    def case_details(case):
        rofi = [float(x) for x in d['rofi']]
        charge = {s: a2.radial_integral([case['rho'][(i, s)] for i in range(nr)], rofi, radial_a, radial_b)
                  for s in (1, 2)}
        dcharge = {s: a2.radial_integral([case['d_rho'][(i, s)] for i in range(nr)], rofi, radial_a, radial_b)
                   for s in (1, 2)}
        labels = {1: 'SUMEV', 2: 'SUMEC', 3: 'RHOVH', 4: 'RHOMU', 5: 'VNUCL',
                  6: 'ZVNUCL', 7: 'EKIN', 8: 'UTOT', 9: 'RHOEPS', 10: 'ETOT'}
        return {
            'sumec_by_spin_ry': case['sumec'], 'sumev_by_spin_ry': case['sumev'],
            'sphere_charge_by_spin': charge, 'sphere_spin_moment': charge[1] - charge[2],
            'actual_l2_radial_charge_by_spin': dcharge,
            'frozen_input_rho_components_ry': {labels[i]: values[0] for i, values in case['components'].items()},
            'output_rho_diagnostic_components_ry': {labels[i]: values[1] for i, values in case['components'].items()},
            'rvh_input_by_spin_ry': case['rvh_in'], 'rvh_output_by_spin_ry': case['rvh_out'],
        }
    a_sign = float(sign['A_wrong_ry'])
    a_true = float(prereg['A_true_ry'])
    e_num = float(prereg['numerical_guard']['E_num_ry'])
    per_spin_results = {}
    for s in (1, 2):
        item = spin['per_spin'][str(s)]
        predicted_delta = float(item['predicted_DeltaSUMEV_raw_minus_wrong_ry'])
        per_spin_results[str(s)] = {
            'true_target_t_ry': true_t[s],
            'wrong_cross_spin_t_ry': wrong_spin_t[s],
            'measured_delta_sumev_ry': spin_delta_s[s],
            'measured_ratio_to_true_t': spin_delta_s[s] / true_t[s],
            'predicted_delta_sumev_raw_minus_wrong_ry': predicted_delta,
            'measured_error_vs_predicted_delta_sumev_ry': spin_delta_s[s] - predicted_delta,
            'A_wrong_ry': float(item['A_wrong_spin_conservative_ry']),
            'E_num_ry': e_num,
            'consistent_with_wrong_target': abs(spin_delta_s[s] - predicted_delta) <=
                float(item['A_wrong_spin_conservative_ry']) + e_num,
            'deviates_from_true_target': abs(spin_delta_s[s] - true_t[s]) >
                float(item['A_true_spin_conservative_ry']) + e_num,
        }
    sign_delta_total = sum(sign_delta_s.values())
    sign_predicted_delta = float(sign['predicted_DeltaSUMEV_raw_minus_wrong_ry'])
    sign_wrong_t_total = sum(wrong_sign_t.values())
    spin_delta_total = sum(spin_delta_s.values())
    spin_predicted_delta = sum(float(spin['per_spin'][str(s)]['predicted_DeltaSUMEV_raw_minus_wrong_ry']) for s in (1, 2))
    spin_wrong_t_total = sum(wrong_spin_t.values())
    summary = {
        'schema': 'rslmto-a3-ks-m03-response-v1',
        'status': 'MEASURED_A3_SIGN_SPIN_RESPONSES',
        'capture_csv_sha256': digest(csv_path),
        'capture_manifest_sha256': digest(manifest_path),
        'source_self_sha256': digest(source),
        'fit_sha256': digest(FIT),
        'a3_preregistration_path': str(PREREG),
        'a3_preregistration_sha256': digest(PREREG),
        'a3_numerical_guard_ry': e_num,
        'driver_path': str(DRIVER),
        'driver_sha256': digest(DRIVER),
        'fixture_input_path': str(input_path),
        'fixture_input_sha256': digest(input_path),
        'fixture_output_path': str(output_path),
        'fixture_output_sha256': digest(output_path),
        'driver_log_path': str(stdout_path),
        'driver_log_sha256': digest(stdout_path),
        'iteration': iteration,
        'case_names': ['RAW_PL_A3_ANCHOR', 'SIGN_CENTER_WRONG', 'OPPOSITE_SPIN_VIN'],
        'raw_sumEV_by_spin_ry': raw['sumev'],
        'case_details': {
            'RAW_PL_A3_ANCHOR': case_details(raw),
            'SIGN_CENTER_WRONG': case_details(sign_case),
            'OPPOSITE_SPIN_VIN': case_details(spin_case),
        },
        'controls': {
            'SIGN_CENTER_WRONG': {
                'definition': 'c_wrong=c+t/q0; correct local q2 held fixed',
                'q2_local_ry2': [float(sign_ql[2, 2, 0]), float(sign_ql[2, 2, 1])],
                'shift_v_ry': sign_shift.tolist(),
                'delta_sumev_by_spin_ry': sign_delta_s,
                'delta_sumev_total_ry': sign_delta_total,
                'predicted_delta_sumev_raw_minus_wrong_ry': sign_predicted_delta,
                'wrong_target_sum_t_ry': sign_wrong_t_total,
                'ratio_to_true_sum_t': sign_delta_total / sum(true_t.values()),
                'error_vs_predicted_delta_sumev_ry': sign_delta_total - sign_predicted_delta,
                'A_wrong_ry': a_sign,
                'A_true_ry': a_true,
                'within_wrong_target_allowance': abs(sign_delta_total - sign_predicted_delta) <= a_sign + e_num,
                'outside_true_target_allowance': abs(sign_delta_total - sum(true_t.values())) > a_true + e_num,
                'sphere_spin_moment': case_details(sign_case)['sphere_spin_moment'],
            },
            'OPPOSITE_SPIN_VIN': {
                'definition': 't_wrong,s=Tr(V_in,other_spin P_out,s); recompute q2_local with other-spin V',
                'q2_local_wrong_ry2': [float(spin_ql[2, 2, 0]), float(spin_ql[2, 2, 1])],
                'shift_v_ry': spin_shift.tolist(),
                'per_spin': per_spin_results,
                'delta_sumev_total_ry': spin_delta_total,
                'predicted_delta_sumev_raw_minus_wrong_ry': spin_predicted_delta,
                'wrong_target_sum_t_ry': spin_wrong_t_total,
                'error_vs_predicted_delta_sumev_total_ry': spin_delta_total - spin_predicted_delta,
                'A_wrong_total_ry': sum(float(spin['per_spin'][str(s)]['A_wrong_spin_conservative_ry']) for s in (1, 2)),
                'within_wrong_target_allowance_per_spin': all(v['consistent_with_wrong_target'] for v in per_spin_results.values()),
                'outside_true_target_allowance_per_spin': all(v['deviates_from_true_target'] for v in per_spin_results.values()),
                'sphere_spin_moment': case_details(spin_case)['sphere_spin_moment'],
            },
        },
        'limitations': 'Private frozen-state production calculate_pl/newrho responses only; no ETOT/M1/M2 claims.',
    }
    report_path = SCRATCH / 'report.json'
    report_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'report': str(report_path), 'report_sha256': digest(report_path),
                      'fixture_output_sha256': digest(output_path)}, sort_keys=True))


if __name__ == '__main__':
    main()
