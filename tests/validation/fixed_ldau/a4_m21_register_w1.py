#!/usr/bin/env python3
"""Build hash-bound W1 response PRE decks after reviewed lambda-zero anchors.

This tool creates registration material only; it never launches the executable.
The anchor report contract is intentionally strict and requires per-mode
standard/tight reconstructed P matrices and immutable run-manifest bindings.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

ROOT = Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
BASE = Path('/private/tmp/a4-m21-ks-lambda0-prep')
OUTPUT = Path('/private/tmp/a4-m21-ks-w1-response-prep')
EVIDENCE = ROOT / 'docs/validation/fixed_ldau/evidence'
LAMBDA0_PREREG = EVIDENCE / 'a4_m21_ks_lambda0_calibration_preregistration.json'
W1_TEMPLATE = EVIDENCE / 'a4_m21_w1_response_preregistration.template.json'
ANALYZER = ROOT / 'tests/validation/fixed_ldau/a4_m21_w1_richardson.py'
SOURCE_MANIFEST = EVIDENCE / 'a4_source_ready_manifest_v3.json'
BINARY = Path('/private/tmp/a4-m2-serial/bin/rslmto.x')
SOURCE_SELF_SHA = '2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010'
BINARY_SHA = '536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4'
SOURCE_MANIFEST_SHA = 'db487ea375fb20dd025ed4d62b4ae880940b865f5fa16c3cefa95220985137ed'
LAMBDA0_PREREG_SHA = '50f24c646af2e7c1eedc050d6b0b3c14b7cb55b5d6134e02cfa081913e1ba80f'
W1_TEMPLATE_SHA = '82ad7759f2840d6288d5f71b7a1f1149844c03f065bade1981e762b996f74434'
ANALYZER_SHA = 'be1f30cabc0d8fe36362282c697b8f477d6bb2b46b658350a9629d04c3650330'
ASSETS = ('clust', 'map', 'sbar', 'view.sbar', 'str.out', 'ves.out', 'mad.mat')
HS = (0.0025, 0.005, 0.01)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_matrix(raw: list) -> list[list[complex]]:
    if len(raw) != 5 or any(len(row) != 5 for row in raw):
        raise ValueError('anchor P_d must be 5x5 [real,imag] pairs')
    out = []
    for row in raw:
        vals = []
        for pair in row:
            if len(pair) != 2 or not all(math.isfinite(float(x)) for x in pair):
                raise ValueError('nonfinite/bad P_d matrix element')
            vals.append(complex(float(pair[0]), float(pair[1])))
        out.append(vals)
    return out


def encode_matrix(m: list[list[complex]]) -> list[list[list[float]]]:
    return [[[z.real, z.imag] for z in row] for row in m]


def read_bounds(raw: list) -> list[list[float]]:
    if len(raw) != 5 or any(len(row) != 5 for row in raw):
        raise ValueError('P_d_component_error_bound must be 5x5')
    out = [[float(x) for x in row] for row in raw]
    if any(not math.isfinite(x) or x < 0.0 for row in out for x in row):
        raise ValueError('invalid P_d_component_error_bound')
    return out


def trace(w: list[list[complex]], p: list[list[complex]]) -> float:
    return sum(w[i][j] * p[j][i] for i in range(5) for j in range(5)).real


def trace_scale(w: list[list[complex]], p: list[list[complex]]) -> float:
    return sum(abs(w[i][j] * p[j][i]) for i in range(5) for j in range(5))


def w1() -> list[list[complex]]:
    return [[complex(i == j, 0.0) for j in range(5)] for i in range(5)]


def anchor_metadata(anchor: dict, mode: str, threshold: str) -> tuple[dict, list[list[complex]], list[list[float]]]:
    key = f'{mode}_{threshold}'
    row = anchor['captures'].get(key)
    if not row:
        raise ValueError(f'missing reviewed anchor {key}')
    if row.get('converged') is not True or row.get('r05_applicability') != 'NOT_APPLICABLE_CONFIRMED':
        raise ValueError(f'anchor {key} is not converged with explicit R05 N/A review')
    mpath = Path(row['run_manifest_path'])
    if not mpath.is_file() or sha(mpath) != row['run_manifest_sha256']:
        raise ValueError(f'anchor run manifest hash mismatch: {key}')
    run = json.loads(mpath.read_text())
    if run.get('status') not in ('PROCESS_EXITED', 'COMPLETE_CONVERGED') or run.get('returncode', 0) != 0:
        raise ValueError(f'anchor process did not exit successfully: {key}')
    bindings = run.get('bindings', {})
    if bindings.get('binary_sha256') != BINARY_SHA:
        raise ValueError(f'anchor binary is not the frozen A4 executable: {key}')
    if not row.get('trace_sha256') or not row.get('final_state_sha256'):
        raise ValueError(f'anchor {key} missing trace/final-state hash')
    post = run.get('post_files', {})
    if post.get('records.csv', {}).get('sha256') != row['trace_sha256']:
        raise ValueError(f'anchor trace hash is not bound by run manifest: {key}')
    if post.get('Fe_out.nml', {}).get('sha256') != row['final_state_sha256']:
        raise ValueError(f'anchor endpoint hash is not bound by run manifest: {key}')
    p = read_matrix(row['P_d'])
    p_bound = read_bounds(row['P_d_component_error_bound'])
    return row, p, p_bound


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('reviewed_anchor_report', type=Path)
    ap.add_argument('--output-root', type=Path, default=OUTPUT)
    ap.add_argument('--registration', type=Path, default=EVIDENCE / 'a4_m21_w1_response_preregistration.json')
    args = ap.parse_args()
    if args.output_root.exists() or args.registration.exists():
        raise SystemExit('refusing to overwrite prior PRE material or response registration')
    for path, expected in ((ROOT / 'source/self.f90', SOURCE_SELF_SHA), (BINARY, BINARY_SHA)):
        if not path.is_file() or sha(path) != expected:
            raise SystemExit(f'frozen source/binary hash mismatch: {path}')
    for path, expected in ((SOURCE_MANIFEST, SOURCE_MANIFEST_SHA),
                           (LAMBDA0_PREREG, LAMBDA0_PREREG_SHA),
                           (W1_TEMPLATE, W1_TEMPLATE_SHA),
                           (ANALYZER, ANALYZER_SHA)):
        if not path.is_file() or sha(path) != expected:
            raise SystemExit(f'pre-response registration dependency hash mismatch: {path}')
    anchor_path = args.reviewed_anchor_report
    anchor = json.loads(anchor_path.read_text())
    if anchor.get('review_status') != 'PASS' or anchor.get('status') != 'LAMBDA0_CALIBRATION_REVIEWED':
        raise SystemExit('anchor report must be explicitly reviewed and classified PASS')
    anchor_sha = sha(anchor_path)
    w = w1()
    eps = 2.220446049250313e-16
    gamma = 256.0 * eps / (1.0 - 256.0 * eps)
    anchors = {}
    anchor_numeric = {}
    for mode in ('spectral', 'local'):
        rstd, pstd, pstd_bound = anchor_metadata(anchor, mode, 'standard')
        rtight, ptight, ptight_bound = anchor_metadata(anchor, mode, 'tight')
        delta = [[ptight[i][j] - pstd[i][j] for j in range(5)] for i in range(5)]
        g0 = trace(w, pstd)
        ug_state = abs(trace(w, delta))
        fp_scale = sum(abs(w[i][j]) * (abs(pstd[j][i]) + abs(ptight[j][i]))
                       for i in range(5) for j in range(5))
        input_bound = sum(abs(w[i][j]) * (pstd_bound[j][i] + ptight_bound[j][i])
                          for i in range(5) for j in range(5))
        ug_fp = gamma * fp_scale + input_bound
        anchors[mode] = {
            'g0_dimensionless': g0,
            'ug0_state_sensitivity_dimensionless': ug_state,
            'ug0_trace_fp_dimensionless': ug_fp,
            'ug0_total_dimensionless': ug_state + ug_fp,
            'standard': {'run_manifest_path': rstd['run_manifest_path'],
                         'run_manifest_sha256': rstd['run_manifest_sha256'],
                         'trace_sha256': rstd['trace_sha256'],
                         'final_state_sha256': rstd['final_state_sha256'],
                         'P_d': encode_matrix(pstd),
                         'P_d_component_error_bound': pstd_bound},
            'tight': {'run_manifest_path': rtight['run_manifest_path'],
                      'run_manifest_sha256': rtight['run_manifest_sha256'],
                      'trace_sha256': rtight['trace_sha256'],
                      'final_state_sha256': rtight['final_state_sha256'],
                      'P_d': encode_matrix(ptight),
                      'P_d_component_error_bound': ptight_bound},
        }
        anchor_numeric[mode] = {'g0': g0, 'ug0': ug_state + ug_fp}

    args.output_root.mkdir(parents=True)
    cases = []
    for mode in ('spectral', 'local'):
        for h in HS:
            for lam in (-h, h):
                for threshold in ('standard', 'tight'):
                    name = f'w1_{mode}_{lam:+.4f}_{threshold}'.replace('+', 'p').replace('-', 'm')
                    basedir = BASE / f'w1_{mode}_{threshold}'
                    outdir = args.output_root / name
                    outdir.mkdir()
                    deck = (basedir / 'input.nml.pre').read_text()
                    pattern = re.compile(r'(?im)^(\s*hubbard_test_lambda\s*=\s*)0\.0d0(\s*(?:!.*)?)$')
                    fortran_lam = f'{lam:+.10e}'.replace('e', 'd')
                    deck, n = pattern.subn(rf'\g<1>{fortran_lam}\g<2>', deck)
                    if n != 1:
                        raise ValueError(f'{name}: expected one lambda=0 field, got {n}')
                    (outdir / 'input.nml.pre').write_text(deck)
                    shutil.copyfile(outdir / 'input.nml.pre', outdir / 'input.nml')
                    shutil.copyfile(basedir / 'Fe.nml.pre', outdir / 'Fe.nml.pre')
                    shutil.copyfile(basedir / 'Fe.nml.pre', outdir / 'Fe.nml')
                    for asset in ASSETS:
                        shutil.copyfile(basedir / asset, outdir / asset)
                    env = {
                        'RSLMTO_FIXED_U_TRACE': str(outdir / 'records.csv'),
                        'RSLMTO_FIXED_U_RUN_ID': name,
                        'RSLMTO_FIXED_U_LAMBDA': f'{lam:.17g}',
                        'RSLMTO_FIXED_U_TRACE_ITERATION': '0',
                        'RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY': '1',
                        'RSLMTO_FIXED_U_A2_SCAN': '0',
                        'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                        'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
                    }
                    cases.append({
                        'run_id': name, 'operator': 'w1', 'mode': mode, 'lambda_ry': lam,
                        'threshold': threshold, 'conv_thr': '1e-9' if threshold == 'standard' else '1e-10',
                        'pre_dir': str(outdir), 'input_nml_pre_sha256': sha(outdir / 'input.nml.pre'),
                        'seed_sha256': sha(outdir / 'Fe.nml.pre'),
                        'static_assets': {a: sha(outdir / a) for a in ASSETS},
                        'command': [str(BINARY)], 'environment': env,
                        'outputs': {'stdout': str(outdir / 'stdout.log'), 'stderr': str(outdir / 'stderr.log'),
                                    'full_scf_log': str(outdir / 'stdout.log'),
                                    'trace': str(outdir / 'records.csv'), 'final_state': str(outdir / 'Fe_out.nml')},
                    })
    reg = {
        'schema': 'rslmto-a4-m21-w1-response-prereg-v1',
        'status': 'PRE_RESPONSE_LOCKED_NO_NONZERO_RUNS',
        'protocol_sha256': '812bcafae105eabaea4ae0614f0b419e0d96638df37888342a8708a973f83980',
        'lambda0_preregistration_sha256': sha(LAMBDA0_PREREG),
        'lambda0_reviewed_anchor_report': {'path': str(anchor_path), 'sha256': anchor_sha},
        'source_manifest_sha256': sha(SOURCE_MANIFEST), 'source_self_sha256': SOURCE_SELF_SHA,
        'w1_analyzer_path': str(ANALYZER), 'w1_analyzer_sha256': ANALYZER_SHA,
        'binary_path': str(BINARY), 'binary_sha256': BINARY_SHA,
        'W1': [[1.0 if i == j else 0.0 for j in range(5)] for i in range(5)],
        'shared_start_seed_sha256': 'aab3f56d5f74048d00d92c458abe65fbe099a4a705ba932134d89cfea1b234df',
        'start_seed_rule': 'Every standard/tight and signed-lambda response copies identical accepted no-U Fe_out bytes; no continuation between response cases.',
        'mode_anchors_and_g0': anchors,
        'criteria': {
            'F_spectral': 'ETOT_spectral - TS',
            'F_local': 'ETOT_local + lambda*Tr(W1*P) - TS',
            'uF_ry': 'abs(F_tight-F_standard)+B_F_tight+B_F_standard+gamma_1*(|F_tight|+|F_standard|); no fixed absolute comparator',
            'B_F': 'Q(ETOT)+B_TS+gamma_3*(|ETOT|+|TS|+|lambda*Tr(WP)|)+B_lambda_trace; Q is derived per source output token; B_TS=Q(TS)+|TS_logged-HP_center|+outward interval radius from raw occupation/k-weight/temperature tokens plus the checked 100/160-digit Decimal enclosure; the same occupation kernel feeds capture rows and entropy; B_trace=gamma_256*sum|W_ij P_ji|+sum|W_ij|*P_component_bound_ji; B_lambda_trace=|lambda|*B_trace+Q_lambda*(|Tr(WP)|+B_trace)',
            'analysis_arithmetic': 'Propagate gamma bounds through centered D differences/divisions, Richardson combinations, scaling differences and residual, truncation estimate, and allowance additions; include decimal-grid-to-binary h discrepancy. Round all positive bounds upward; compare B_W via exact Decimal subtraction of binary64 R_fine and g0.',
            'ug0': 'abs(Tr(W1*(P_tight-P_standard))) + gamma_256*sum_ij abs(W1_ij)*(|P_standard_ji|+|P_tight_ji|) + sum_ij abs(W1_ij)*(Qcomplex(P_standard_ji)+Qcomplex(P_tight_ji)); dynamic scales/input error for both matrices',
            'D_h': '(F_standard(+h)-F_standard(-h))/(2h)',
            'uD_h': '(uF(+h)+uF(-h))/(2h)',
            'h_ry': list(HS),
            'richardson_fine': '(4D(.0025)-D(.005))/3',
            'richardson_coarse': '(4D(.005)-D(.01))/3',
            'uR_fine': '(4uD(.0025)+uD(.005))/3',
            'uR_coarse': '(4uD(.005)+uD(.01))/3',
            'B_trunc': 'abs(R_fine-R_coarse)/15+(uR_fine+uR_coarse)/15',
            'N_unc': 'uR_fine+B_trunc+ug0; exclude measured B_W',
            'B_W': 'R_fine-g0; pass only when h^2 scaling resolved and abs(B_W)<=N_unc',
            'h2_resolution': 'Approved h2-dominance rule: A=D(.005)-D(.0025), B=D(.01)-D(.005), L=(16A-B)/36, H=(B-4A)/180. Propagate uA/uB from the frozen uD and arithmetic bounds, then uL/uH through these linear combinations. Require |L|>uL and |L|-uL>16*(|H|+uH), with comparison performed on exact Decimal representations of the computed binary64 values. This replaces only the former strict 4:1 resolution check; g0, ug0, N_unc, Richardson and identity bias gate are unchanged.',
            'R05': 'At U=J=0, require literal NOT APPLICABLE; tighten conv_thr by 10x for tight pairs and do not invent R05 thresholds.',
            'forbidden': ['No extra lambda spacing', 'No bias inside allowance', 'No arbitrary serialization tolerance', 'No A4.2 transfer'],
        },
        'response_cases': cases,
        'execution_status': 'PREPARED_ONLY_NO_PROCESS_LAUNCHED',
    }
    reg_path = args.registration
    reg_path.parent.mkdir(parents=True, exist_ok=True)
    reg_path.write_text(json.dumps(reg, indent=2, sort_keys=True) + '\n')
    summary = {'registration_path': str(reg_path), 'registration_sha256': sha(reg_path),
               'prep_root': str(args.output_root), 'case_count': len(cases),
               'anchor_g0_ug0': anchor_numeric, 'status': reg['status']}
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
