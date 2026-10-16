#!/usr/bin/env python3
"""Prepare, but never execute, the eight A4.1 KS lambda=0 calibration cases."""
from __future__ import annotations
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
BASE = Path('/private/tmp/rslmto-localmom-m02-release/ks')
PREP = Path('/private/tmp/a4-m21-ks-lambda0-prep')
BINARY = Path('/private/tmp/a4-m2-serial/bin/rslmto.x')
SOURCE_MANIFEST = ROOT / 'docs/validation/fixed_ldau/evidence/a4_source_ready_manifest_v3.json'
PROTOCOL_SHA = '812bcafae105eabaea4ae0614f0b419e0d96638df37888342a8708a973f83980'
SOURCE_SELF_SHA = '2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010'
BINARY_SHA = '536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4'
BASE_INPUT_SHA = 'ab09e896ea9cebcae601547142c4d2474090901b3f988bd198f2798e8ee2d72e'
SEED_SHA = 'aab3f56d5f74048d00d92c458abe65fbe099a4a705ba932134d89cfea1b234df'
ASSETS = ('clust', 'map', 'sbar', 'view.sbar', 'str.out', 'ves.out', 'mad.mat')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    assert sha(BINARY) == BINARY_SHA
    assert sha(SOURCE_MANIFEST) == 'db487ea375fb20dd025ed4d62b4ae880940b865f5fa16c3cefa95220985137ed'
    assert sha(ROOT / 'source/self.f90') == SOURCE_SELF_SHA
    assert sha(BASE / 'input.nml.pre') == BASE_INPUT_SHA
    assert sha(BASE / 'Fe_out.nml') == SEED_SHA
    if PREP.exists():
        raise SystemExit(f'refusing to overwrite existing preparation directory: {PREP}')
    PREP.mkdir(parents=True)
    cases = []
    for operator in ('w1', 'w2'):
        for mode in ('spectral', 'local'):
            for threshold in ('standard', 'tight'):
                name = f'{operator}_{mode}_{threshold}'
                d = PREP / name
                d.mkdir()
                deck = (BASE / 'input.nml.pre').read_text()
                if threshold == 'tight':
                    deck, n = deck.replace('conv_thr = 1d-9', 'conv_thr = 1d-10'), deck.count('conv_thr = 1d-9')
                    if n != 1:
                        raise ValueError(f'{name}: expected exactly one base conv_thr, got {n}')
                block = (
                    '\n&hamiltonian\n'
                    'hubbard_u_general(1,3) = 0.0d0\n'
                    'hubbard_j_general(1,3) = 0.0d0\n'
                    f"hubbard_test_operator = '{operator}'\n"
                    'hubbard_test_lambda = 0.0d0\n'
                    f"hubbard_radial_moments = '{mode}'\n"
                    'fixed_u_test_tighten_r05 = .true.\n'
                    '/\n'
                )
                deck += block
                pre = d / 'input.nml.pre'
                pre.write_text(deck)
                shutil.copyfile(pre, d / 'input.nml')
                seed = d / 'Fe.nml.pre'
                shutil.copyfile(BASE / 'Fe_out.nml', seed)
                shutil.copyfile(seed, d / 'Fe.nml')
                for asset in ASSETS:
                    shutil.copyfile(BASE / asset, d / asset)
                case = {
                    'run_id': name,
                    'operator': operator,
                    'radial_moments': mode,
                    'lambda_ry': 0.0,
                    'threshold': threshold,
                    'conv_thr': '1e-9' if threshold == 'standard' else '1e-10',
                    'r05_tightening_flag': True,
                    'r05_applicability': 'NOT_APPLICABLE at numeric U=J=0 external-only test hook; confirm literal runtime log and no 1e-8 R05 targets',
                    'working_directory': str(d),
                    'input_nml_pre_path': str(pre),
                    'input_nml_pre_sha256': sha(pre),
                    'seed_source_path': str(BASE / 'Fe_out.nml'),
                    'seed_source_sha256': SEED_SHA,
                    'Fe_nml_pre_path': str(seed),
                    'Fe_nml_pre_sha256': sha(seed),
                    'static_assets': {a: {'path': str(d / a), 'sha256': sha(d / a)} for a in ASSETS},
                    'command': [str(BINARY)],
                    'environment': {
                        'RSLMTO_FIXED_U_TRACE': str(d / 'records.csv'),
                        'RSLMTO_FIXED_U_RUN_ID': name,
                        'RSLMTO_FIXED_U_LAMBDA': '0.0',
                        'RSLMTO_FIXED_U_TRACE_ITERATION': '0',
                        'RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY': '1',
                        'RSLMTO_FIXED_U_A2_SCAN': '0',
                        'OMP_NUM_THREADS': '1',
                        'OPENBLAS_NUM_THREADS': '1',
                        'MKL_NUM_THREADS': '1',
                        'VECLIB_MAXIMUM_THREADS': '1',
                    },
                    'outputs': {
                        'full_scf_log': str(d / 'testrun.log'),
                        'rolling_last_iteration_trace': str(d / 'records.csv'),
                        'final_state': str(d / 'Fe_out.nml'),
                    },
                }
                cases.append(case)
    prereg = {
        'schema': 'rslmto-a4-m21-ks-lambda0-calibration-prereg-v1',
        'status': 'PREPARED_PRE_RESPONSE_NOT_EXECUTED',
        'purpose': 'Eight zero-coupling calibration captures only. This is not an A4.1 identity response or PASS.',
        'governance': {
            'protocol_sha256': PROTOCOL_SHA,
            'source_ready_manifest_path': str(SOURCE_MANIFEST),
            'source_ready_manifest_sha256': sha(SOURCE_MANIFEST),
            'source_self_sha256': SOURCE_SELF_SHA,
            'binary_path': str(BINARY),
            'binary_sha256': BINARY_SHA,
        },
        'common_state': {
            'backend': 'KS',
            'base_deck_path': str(BASE / 'input.nml.pre'),
            'base_deck_sha256': BASE_INPUT_SHA,
            'base_seed_path': str(BASE / 'Fe_out.nml'),
            'base_seed_sha256': SEED_SHA,
            'base_run_manifest_path': str(BASE / 'run_manifest.json'),
            'base_run_manifest_sha256': sha(BASE / 'run_manifest.json'),
            'base_status': 'COMPLETE_CONVERGED',
            'accepted_iterations': 23,
            'physics': {'nk': [12, 12, 12], 'temperature_K': 300.0, 'nsp': 1, 'lmax': 2,
                        'ct': 4.0, 'r2': 16.0, 'energy_window_ry': [-1.0, 0.5], 'channels_ldos': 2500,
                        'U_eV': 0.0, 'J_eV': 0.0, 'lambda_ry': 0.0},
            'shared_warm_start': 'Every case receives the exact same accepted Fe_out bytes as Fe.nml. No case restarts from another case output.',
        },
        'design': {
            'operators': ['w1', 'w2'], 'modes': ['spectral', 'local'],
            'threshold_pair': {'standard': '1e-9', 'tight': '1e-10'},
            'only_threshold_difference_within_pair': 'self%conv_thr is divided by ten; the test-only R05 flag is true in both, but numeric U=J=0 makes fixed-U R05 inapplicable.',
            'other_changes_from_base': ['explicit U=J=0', 'operator selection at lambda=0', 'radial moment selector', 'fixed_u_test_tighten_r05=.true.'],
            'trace': 'capture-only full SCF; TRACE_ITERATION=0 retains the last completed iteration CSV; full iteration history remains in testrun.log.',
        },
        'pre_response_checks': [
            'Verify every run reaches ordinary SCF convergence from its PRE seed and CSV final iteration matches final endpoint.',
            'Require literal log text for test operator and lambda=0; require literal R05 NOT APPLICABLE and absence of tight R05 targets for all U=0 cases.',
            'Require complete eigenpair, occupation, k-weight, basis-map, Hamiltonian, M1, PL, and radial/energy capture rows.',
            'Reconstruct per-spin d P from captured eigenvectors, occupations and k weights; compute g0=Re sum_s Tr(W P_s) independently for each operator/mode/threshold.',
            'Compute ug0=abs(Tr(W(P_t-P_s))) + gamma_256*sum_ij abs(W_ij P_s,ji), separately per operator/mode; gamma_256=256 eps/(1-256 eps).',
            'At lambda=0 verify the inserted v_hub d-block is exactly zero and compare P/eigenvalue/occupation outputs across operator and mode selectors at fixed threshold for exact identity; report entropy and energy candidates as calibration only.',
            'Do not begin any nonzero-lambda run until these anchors, g0/ug0, and source/input/output provenance are locked in a separate response registration.',
        ],
        'cases': cases,
    }
    out = ROOT / 'docs/validation/fixed_ldau/evidence/a4_m21_ks_lambda0_calibration_preregistration.json'
    if out.exists():
        raise SystemExit(f'refusing to overwrite preregistration: {out}')
    out.write_text(json.dumps(prereg, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'preregistration': str(out), 'prereg_sha256': sha(out),
                      'case_count': len(cases), 'prep_root': str(PREP)}, sort_keys=True))

if __name__ == '__main__':
    main()
