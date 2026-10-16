#!/usr/bin/env python3
"""Run the locked A3 nonconverged-iteration control on the accepted KS state."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import accepted_m03_replay as a2  # noqa: E402

CAPTURE = Path('/private/tmp/rslmto-localmom-a2-ks-u2-capture-d91')
ITER1 = Path('/private/tmp/rslmto-a3-u2-iter1-same-run')
DRIVER = Path('/private/tmp/rslmto-localmom-a2-u2-d91-bin/accepted_m03_replay_driver')
FIT = Path('/private/tmp/a2-accepted-iter23-a22-fit.json')
LOCK = Path('/private/tmp/a3-iteration-control-prereg-locked.json')
LOCK_SHA = '3417f716158a951b68c7d5a707a82ea8dfb701416a87be9a8b8011ec532b8158'
FIT_SHA = '0ad6a7703aa26eb132602e7a43773d7619f70651c6b90817e4f945ff259cb993'
SOURCE_SHA = 'd91a468a4493da21710fc7f60e745b5f1bef8f20d5d667ae939c31d5ba3061ef'
DRIVER_SHA = '0c1c3c12e8501f4642af87168c611d261fb5a34cc4747aaef9c8420c2c6449c4'
OUT = Path('/private/tmp/a3-m03-iteration-response')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_matrix(item: dict) -> np.ndarray:
    return np.asarray(item['real'], dtype=float) + 1j*np.asarray(item['imag'], dtype=float)


def captured_vin(path: Path) -> list[np.ndarray]:
    full = np.zeros((18, 18), dtype=np.complex128)
    seen = set()
    with (path/'records.csv').open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['kind'] != 'v_hub' or row['stage'] != 'Pin_H' or int(row['type']) != 1:
                continue
            i, j = int(row['row']), int(row['col'])
            key = (i, j)
            if key in seen:
                raise ValueError(f'duplicate V_in matrix row {key}')
            seen.add(key)
            full[i, j] = complex(float(row['re']), float(row['im']))
    if len(seen) != 18*18:
        raise ValueError(f'incomplete iteration-1 V_in matrix: {len(seen)} entries')
    return [full[4+9*s:9+9*s, 4+9*s:9+9*s] for s in (0, 1)]


def normalized_trajectory(path: Path) -> list[str]:
    rx = re.compile(r'^\[[^]]+\]\s+\[[^]]+\]\s+.*?(FIXED_U_RESIDUAL.*|Charge transfer at atom.*|Excess charge.*|Not converged! Diff=.*|Converged!.*)$')
    result = []
    for line in (path/'testrun.log').read_text(errors='replace').splitlines():
        match = rx.match(line)
        if match:
            result.append(match.group(1))
    return result


def case_details(case: dict, rofi: list[str], a: float, b: float, nr: int) -> dict:
    mesh = [float(x) for x in rofi]
    charge = {s: a2.radial_integral([case['rho'][(i, s)] for i in range(nr)], mesh, a, b)
              for s in (1, 2)}
    dcharge = {s: a2.radial_integral([case['d_rho'][(i, s)] for i in range(nr)], mesh, a, b)
               for s in (1, 2)}
    labels = {1: 'SUMEV', 2: 'SUMEC', 3: 'RHOVH', 4: 'RHOMU', 5: 'VNUCL',
              6: 'ZVNUCL', 7: 'EKIN', 8: 'UTOT', 9: 'RHOEPS', 10: 'ETOT'}
    return {
        'sumec_by_spin_ry': case['sumec'], 'sumev_by_spin_ry': case['sumev'],
        'sphere_charge_by_spin': charge, 'sphere_spin_moment': charge[1]-charge[2],
        'actual_l2_radial_charge_by_spin': dcharge,
        'frozen_input_rho_components_ry': {labels[i]: v[0] for i, v in case['components'].items()},
        'output_rho_diagnostic_components_ry': {labels[i]: v[1] for i, v in case['components'].items()},
        'rvh_input_by_spin_ry': case['rvh_in'], 'rvh_output_by_spin_ry': case['rvh_out'],
    }


def main() -> None:
    if OUT.exists():
        raise SystemExit(f'refusing to overwrite {OUT}')
    lock = json.loads(LOCK.read_text())
    if sha(LOCK) != LOCK_SHA or lock.get('status') != 'PREREGISTERED_BEFORE_ITERATION_CONTROL_RESPONSE; ADMISSIBLE':
        raise ValueError('A3 iteration-control lock hash/status mismatch')
    if sha(FIT) != FIT_SHA or sha(DRIVER) != DRIVER_SHA:
        raise ValueError('locked A2.2 fit or production fixture driver mismatch')
    source = Path(__file__).resolve().parents[3]/'source/self.f90'
    if sha(source) != SOURCE_SHA:
        raise ValueError('source differs from accepted d91 capture')

    iter1_csv = ITER1/'records.csv'
    if sha(iter1_csv) != lock['provenance']['iter1_csv_sha256']:
        raise ValueError('iteration-1 capture differs from the preregistered replay')
    if sha(ITER1/'Fe_out.nml') != lock['provenance']['accepted_endpoint_Fe_out_sha256']:
        raise ValueError('iteration-1 replay final restart differs from accepted endpoint')
    replay_history, accepted_history = normalized_trajectory(ITER1), normalized_trajectory(CAPTURE)
    if len(replay_history) != 128 or replay_history != accepted_history:
        raise ValueError('same-trajectory charge/R05 log history mismatch')

    rows, manifest, iteration, nr, ncore, nval, csv_path, manifest_path = a2.export_rows(CAPTURE)
    if (iteration != 32 or sha(csv_path) != lock['provenance']['accepted_csv_sha256'] or
            sha(manifest_path) != lock['provenance']['accepted_manifest_sha256']):
        raise ValueError('accepted iteration-32 capture or manifest differs from A3 lock')
    if manifest.get('source_self_sha256') != SOURCE_SHA:
        raise ValueError('accepted capture source mismatch')

    d = a2.prepare(rows, nr, ncore, nval, 'cases', manifest)
    target = d['target']
    wrong_v = captured_vin(ITER1)
    locked = lock['spins']
    shifts = np.zeros(2)
    ql_iter1 = d['cases'][0][3].copy()
    predicted = {}
    for idx, spin in enumerate((1, 2)):
        item = locked[str(spin)]
        if np.max(np.abs(wrong_v[idx]-as_matrix(item['V_in_iter1_dblock_Ry']))) > 1.e-14:
            raise ValueError(f'V1 d-block differs from locked matrix for spin {spin}')
        if np.max(np.abs(target['p'][idx]-as_matrix(item['P_OUT_accepted32']))) > 1.e-14:
            raise ValueError(f'accepted P_OUT differs from locked matrix for spin {spin}')
        p, m1 = target['p'][idx], target['m1'][idx]
        q0, center, q2 = target['q0'][idx], target['center'][idx], target['q2'][idx]
        v = wrong_v[idx]
        t = float(np.trace(v@p).real)
        c = center-t/q0
        q2local = q2-2*float(np.trace(v@(m1-center*p)).real)+float(np.trace(v@v@p).real)-t*t/q0
        if q2local < 0:
            raise ValueError(f'iteration-1 q2_local negative for spin {spin}: {q2local}')
        shifts[idx] = t/q0
        ql_iter1[2, 2, idx] = q2local
        predicted[str(spin)] = {'t_wrong_ry': t, 'shift_v_ry': shifts[idx],
            'q2_local_ry2': q2local, 'center_ry': c,
            'locked_t_plus_rfit_ry': item['predicted_wrong_carrier_t_plus_rfit_ry'],
            'locked_true_t_plus_rfit_ry': item['predicted_true_carrier_t_plus_rfit_ry'],
            'A_wrong_ry': item['A_k_s_ry'], 'A_true_ry': item['A_true_s_ry']}

    raw_ql = d['cases'][0][3].copy()
    d['cases'] = [('RAW_PL_A3_C_ANCHOR', np.zeros(2), True, raw_ql),
                  ('ITER1_VIN_CONTROL', shifts, True, ql_iter1)]
    payload = a2.text_input(d, rows, nr, ncore, nval).encode()
    OUT.mkdir(parents=True)
    inp, output, log = OUT/'fixture_input.txt', OUT/'fixture_output.csv', OUT/'driver.stdout.log'
    inp.write_bytes(payload)
    proc = subprocess.run([str(DRIVER.resolve()), str(inp), str(output)], capture_output=True, text=True)
    log.write_text(proc.stdout+proc.stderr)
    if proc.returncode:
        raise RuntimeError(f'production frozen-state driver failed; see {log}')
    cases = a2.output_cases(output.read_bytes())
    if set(cases) != {1, 2}:
        raise ValueError(f'unexpected response cases {sorted(cases)}')
    raw, wrong = cases[1], cases[2]
    prereg = lock['spins']
    atrue = {}; awrong = {}; measured = {}; errors_wrong = {}; errors_true = {}
    for spin in (1, 2):
        item = prereg[str(spin)]
        measured[str(spin)] = raw['sumev'][spin]-wrong['sumev'][spin]
        errors_wrong[str(spin)] = measured[str(spin)]-float(item['predicted_wrong_carrier_t_plus_rfit_ry'])
        errors_true[str(spin)] = measured[str(spin)]-float(item['predicted_true_carrier_t_plus_rfit_ry'])
        awrong[str(spin)] = float(item['A_k_s_ry'])
        atrue[str(spin)] = float(item['A_true_s_ry'])
    e_num = float(lock['allowance_model']['E_num_ry'])
    radial_a = float(a2.one(rows, 'radial_A', 'RADIAL_REPLAY')['scalar'])
    radial_b = float(a2.one(rows, 'radial_B', 'RADIAL_REPLAY')['scalar'])
    result = {
        'schema': 'rslmto-a3-ks-m03-iteration-control-response-v1',
        'status': 'MEASURED_A3_ITERATION_CONTROL_RESPONSE',
        'lock_path': str(LOCK), 'lock_sha256': sha(LOCK),
        'iter1_csv_sha256': sha(iter1_csv), 'iter1_log_sha256': sha(ITER1/'testrun.log'),
        'accepted_csv_sha256': sha(csv_path), 'accepted_manifest_sha256': sha(manifest_path),
        'accepted_endpoint_Fe_out_sha256': sha(CAPTURE/'Fe_out.nml'),
        'source_sha256': sha(source), 'driver_sha256': sha(DRIVER),
        'fixture_input_path': str(inp), 'fixture_input_sha256': sha(inp),
        'fixture_output_path': str(output), 'fixture_output_sha256': sha(output),
        'driver_log_sha256': sha(log), 'same_trajectory_history_rows': len(replay_history),
        'same_trajectory_history_equal': True,
        'case_names': ['RAW_PL_A3_C_ANCHOR','ITER1_VIN_CONTROL'],
        'shifts_v_ry': shifts.tolist(), 'locked_predictions_by_spin': predicted,
        'measured_delta_sumev_by_spin_ry': measured,
        'error_vs_locked_wrong_carrier_by_spin_ry': errors_wrong,
        'error_vs_locked_true_carrier_by_spin_ry': errors_true,
        'A_wrong_by_spin_ry': awrong, 'A_true_by_spin_ry': atrue,
        'E_num_ry': e_num,
        'within_wrong_allowance_by_spin': {s: abs(errors_wrong[s]) <= awrong[s]+e_num for s in ('1','2')},
        'outside_true_allowance_by_spin': {s: abs(errors_true[s]) > atrue[s]+e_num for s in ('1','2')},
        'q2_local_nonnegative': all(v['q2_local_ry2'] >= 0 for v in predicted.values()),
        'case_details': {
            'RAW_PL_A3_C_ANCHOR': case_details(raw, d['rofi'], radial_a, radial_b, nr),
            'ITER1_VIN_CONTROL': case_details(wrong, d['rofi'], radial_a, radial_b, nr),
        },
        'limitations': 'Frozen-state production calculate_pl/newrho only; no ETOT correction or M1/M2 claim.'
    }
    report = OUT/'report.json'; report.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'report_path':str(report),'report_sha256':sha(report),
                      'fixture_output_sha256':sha(output),'measured_delta_sumev_by_spin_ry':measured,
                      'error_vs_locked_wrong_carrier_by_spin_ry':errors_wrong},sort_keys=True))


if __name__ == '__main__':
    main()
