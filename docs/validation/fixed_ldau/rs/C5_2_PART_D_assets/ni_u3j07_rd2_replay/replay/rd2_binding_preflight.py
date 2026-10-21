#!/usr/bin/env python3
"""Structural/hash preflight for the C5 RD2 Ni frozen replay.

This creates a binding record and exact LOCAL_HALL Vin operand only. It does
not run the replay driver or calculate RD2 response gates.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from collections import Counter
from pathlib import Path

REPO = Path('/Users/andersb/Jobb/rslmto_fable_v4b')
CAP = Path('/private/tmp/c52-rd2-ni-u3j07-capture')
OUT = Path('/private/tmp/c52-rd2-private')
RUN_ID = 'rd2_ni_u3j07_endpoint23_fullgreen_i1'
ITERATION = 1
LOCAL_SITE = 1
LOCAL_TYPE = 2
BUILD_ID = 1
NORB = 18


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text())


def need(ok: bool, text: str) -> None:
    if not ok:
        raise ValueError(text)


def import_c2_parser():
    path = REPO / 'tests/validation/fixed_ldau/c2_hubbard_operator_audit.py'
    spec = importlib.util.spec_from_file_location('rd2_c2_audit_parser', path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod, path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    preflight_path = CAP / 'provenance/preflight.json'
    erratum_path = CAP / 'analysis/runtime_mapping_erratum.json'
    trace_path = CAP / 'run/records.csv'
    audit_path = CAP / 'run/c2_hubbard_audit.dat'
    input_path = CAP / 'immutable/input_capture.nml'
    endpoint_path = Path('/private/tmp/c52-rd-ni-impurity-u3j07/accepted_checkpoint_23_window4_23/window_state_0023_Ni-imp_scf.nml')
    driver = Path('/private/tmp/c52-c-private/helper-build/accepted_m03_replay_driver')
    helper_record_path = Path('/private/tmp/c52-c-private/helper-build/helper_build_record.json')
    helper_manifest_path = Path('/private/tmp/c52-c-private/helper-build/source_build_manifest.json')
    driver_source = REPO / 'tests/validation/fixed_ldau/accepted_m03_replay_driver.f90'
    runner_path = REPO / 'tests/validation/fixed_ldau/rs_a3_m03_runner.py'
    fit_path = REPO / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_fit_results.json'
    fit_prereg_path = REPO / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_preregistration_v3.json'
    moment_math_path = REPO / 'source/math.f90'

    preflight, erratum = load_json(preflight_path), load_json(erratum_path)
    need(erratum['original_preflight_sha256'] == sha(preflight_path), 'erratum does not bind original preflight')
    need(erratum['capture_audit']['sha256'] == sha(audit_path), 'erratum does not bind LOCAL_HALL audit')
    need(erratum['runtime_source_of_truth']['str_out_sha256'] == sha(CAP / 'run/str.out'), 'erratum str.out hash mismatch')
    need(erratum['runtime_source_of_truth']['outnewclu_sha256'] == sha(CAP / 'run/outnewclu'), 'erratum outnewclu hash mismatch')
    need(erratum['runtime_recursion_map'] if 'runtime_recursion_map' in erratum else
         erratum['runtime_source_of_truth']['runtime_recursion_map'], 'missing runtime recursion map')
    # The old preflight text is historical and intentionally retained. The
    # runtime erratum is authoritative: irec(1)=cluster 1 -> iz(1)=type 2.
    mapdata = erratum['runtime_source_of_truth']['runtime_recursion_map']
    need(mapdata == {'recursion_index_1': 1, 'irec_1_cluster_site': 1, 'iz_cluster_site_1_type': 2},
         f'unexpected runtime map: {mapdata}')
    need(erratum['capture_audit']['recursion_site_rows']['META'].startswith('META LOCAL_HALL 1 1 2 Ni-imp'),
         'erratum does not identify the recursion-site META row')
    need(sha(input_path) == preflight['input']['capture_input_sha256'], 'capture input hash mismatch')
    endpoint = preflight['accepted_endpoint']
    need(sha(endpoint_path) == endpoint['checkpoint_sha256'], 'accepted endpoint hash mismatch')
    need(sha(driver) == '9e3f44379b8933d21eafb78e8008cbf0d44dd9b65a77140af2920c26713ac9d7', 'driver hash mismatch')
    helper_record = load_json(helper_record_path)
    helper_manifest = load_json(helper_manifest_path)
    need(helper_record['helper_binary_sha256'] == sha(driver), 'helper build record binary mismatch')
    need(helper_record['archive_sha256'] == '36725c2ae3ba7c07c74c17940d9047bf8a314431146782bdd4e09ecc0e5dcead', 'helper is not linked to audited B1 archive')
    need(helper_manifest['head'] == '8a107432e3ca2e00c815c217d52e7aa2a8e30736', 'helper source build head mismatch')
    need(helper_manifest['all_compiled_tracked_sources_match_head'] is True, 'helper source manifest is not complete')
    need(helper_record['helper_source_sha256'] == sha(driver_source), 'replay driver source mismatch')
    need(sha(fit_path) == '0f0b543f81f53bc1c7fdfc46d3b404adef1ca455b8366a90e87154afd2b0987c', 'A3/RB3 fit changed')
    need(sha(fit_prereg_path) == '46ca7042678dc5b74ad15a7c2bc5f66786ca9191f70d2395405e78c69b5b1525', 'A3/RB3 fit preregistration changed')
    fit_data, fit_prereg = load_json(fit_path), load_json(fit_prereg_path)
    need(fit_data.get('prereg_sha256') == sha(fit_prereg_path), 'fit does not bind its registered preregistration')
    fit_spin_map = {}
    for spin in ('1', '2'):
        model = fit_data['per_spin'][spin]['fit']
        fit_spin_map[spin] = {
            'source_fit_path_key': f'per_spin.{spin}.fit',
            'rb3_generic_fit_key_only': f'site=0,type=0,spin={spin}',
            'r_cubic_coefficients_ascending': model['r_cubic_coefficients_ascending'],
            'e_even_a0_ry': model['e_even_fit_coefficients_ascending'][0],
        }
    d0_ry = float(fit_prereg['accepted_no_u_state']['D0_RS_ry'])
    scalar_fit_response_at_plus_001 = {
        spin: sum(float(c) * (0.01 ** j) for j, c in enumerate(fit_spin_map[spin]['r_cubic_coefficients_ascending']))
        for spin in ('1', '2')
    }
    scalar_fit_allowance_ry = sum(abs(x) for x in scalar_fit_response_at_plus_001.values()) + d0_ry
    need(sha(runner_path) == '96c51a8e1e0b0809b35b7875a5c94b9be1e70e82180cb6693c13209fe9a0b515', 'shared runner changed')

    # Validate the captured CSV is one finite, complete full-Green iteration.
    with trace_path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        header = reader.fieldnames
    need(rows and header, 'empty capture CSV')
    need(len(rows) == 552196, f'unexpected captured row count {len(rows)}')
    need({r['run_id'] for r in rows} == {RUN_ID}, 'run id mismatch')
    need({int(r['iteration']) for r in rows} == {ITERATION}, 'iteration mismatch')
    for number, row in enumerate(rows, 1):
        for key in ('weight', 'energy', 'scalar', 're', 'im'):
            need(math.isfinite(float(row[key])), f'nonfinite {key} at CSV row {number}')
    counts = Counter((r['kind'], r['stage'], int(r['site']), int(r['type']), int(r['spin'])) for r in rows)
    basis = [r for r in rows if r['kind'] == 'basis_map' and r['stage'] == 'RS_GREEN' and int(r['site']) == 0]
    need(len(basis) == NORB, f'site0 basis map extent is {len(basis)}, expected 18')
    need({int(r['type']) for r in basis} == {LOCAL_TYPE}, 'site0 Green basis is not mapped to type 2')
    green = [r for r in rows if r['kind'] == 'green_function' and r['stage'] == 'RS_GREEN' and int(r['site']) == 0]
    energy_indices = sorted({int(r['k']) for r in green})
    need(len(energy_indices) == 1567, f'unexpected site0 Green energy count {len(energy_indices)}')
    green_counts = Counter(int(r['k']) for r in green)
    need(all(green_counts[k] == NORB*NORB for k in energy_indices), 'incomplete site0 Green matrix')
    for stage in ('Pin_H', 'P_OUT'):
        for ell in range(3):
            for spin in (1, 2):
                selected = [r for r in rows if r['kind'] == 'ldm' and r['stage'] == stage and int(r['site']) == 0
                            and int(r['l']) == ell and int(r['spin']) == spin]
                need(len(selected) == (2*ell+1)**2, f'incomplete {stage} P rows l={ell} spin={spin}')
    for ell in range(3):
        for spin in (1, 2):
            selected = [r for r in rows if r['kind'] == 'ldm_m1' and r['stage'] == 'M1_OUT' and int(r['site']) == 0
                        and int(r['l']) == ell and int(r['spin']) == spin]
            need(len(selected) == (2*ell+1)**2, f'incomplete M1 rows l={ell} spin={spin}')

    # Select the additive operand recorded at the exact LOCAL_HALL insertion,
    # keyed by runtime cluster site/type, never by a bulk type selector.
    c2, c2_path = import_c2_parser()
    parsed = c2.read_records(audit_path)
    rec = c2.choose_record(parsed, 'LOCAL_HALL', LOCAL_SITE, BUILD_ID)
    need(int(rec['ptype']) == LOCAL_TYPE and rec['symbol'] == 'Ni-imp', 'LOCAL_HALL runtime target mismatch')
    need(int(rec['count']) == 1 and len(rec['block']) == 50, 'LOCAL_HALL insertion is incomplete')
    vin_full = [[0j for _ in range(NORB)] for _ in range(NORB)]
    for (spin, i, j), value in rec['block'].items():
        need(value['count'] == 1, 'LOCAL_HALL block has repeated insertion')
        # The audit's exact additive field is the actual operator operand;
        # place each 5x5 spin block at the d-shell basis slice 4:9 / 13:18.
        offset = 4 + 9 * (spin - 1)
        vin_full[offset+i-1][offset+j-1] = value['correction']
    vhub = {}
    for spin in (1, 2):
        block = [[None for _ in range(5)] for _ in range(5)]
        seen = set()
        for row in rows:
            if row['kind'] == 'v_hub' and row['stage'] == 'Pin_H' and int(row['type']) == LOCAL_TYPE:
                i, j = int(row['row']), int(row['col'])
                if i in range(4+9*(spin-1), 9+9*(spin-1)) and j in range(4+9*(spin-1), 9+9*(spin-1)):
                    block[i-(4+9*(spin-1))][j-(4+9*(spin-1))] = complex(float(row['re']), float(row['im']))
                    seen.add((i,j))
        need(len(seen) == 25 and all(x is not None for line in block for x in line), f'incomplete v_hub type2 Pin_H block spin={spin}')
        vhub[spin] = block
    errors = {}
    for spin in (1, 2):
        off = 4+9*(spin-1)
        errors[str(spin)] = max(abs(vin_full[off+i][off+j]-vhub[spin][i][j]) for i in range(5) for j in range(5))
    need(max(errors.values()) <= 5e-14, f'audit insertion differs from trace type2 Pin_H d blocks: {errors}')

    # Save exact decimal tokens from the audited LOCAL_HALL correction. This
    # is the only permitted RD2 Vin binding; no replay/gate cases are run here.
    local_lines = [line.split() for line in audit_path.read_text().splitlines()
                   if line.startswith('BLOCK LOCAL_HALL 1 1 2 ')]
    need(len(local_lines) == 50, f'expected exactly 50 LOCAL_HALL site1/type2 blocks, got {len(local_lines)}')
    operand = {
        'schema': 'c5-rd2-local-hall-vin-binding-v1',
        'status': 'STRUCTURAL_BINDING_ONLY_NO_RD2_GATE_OR_REPLAY_EXECUTED',
        'source': {'path': str(audit_path), 'sha256': sha(audit_path), 'record_parser': str(c2_path), 'parser_sha256': sha(c2_path)},
        'runtime_map': {'recursion_ordinal': 1, 'irec_cluster_site': 1, 'nlim': LOCAL_SITE,
                        'iz_type': LOCAL_TYPE, 'symbol': rec['symbol'], 'host_sites': '2..19, type1'},
        'iteration_and_capture': {'run_id': RUN_ID, 'iteration': ITERATION, 'csv_sha256': sha(trace_path),
                                  'rows': len(rows), 'site0_is_cluster_site1_source': 'self.f90 writes site=ia-1'},
        'operator': {'assembler': 'LOCAL_HALL', 'build_id': BUILD_ID, 'site': LOCAL_SITE, 'type': LOCAL_TYPE,
                     'u_ry': rec['u_ry'], 'j_ry': rec['j_ry'], 'insert_count': rec['count'],
                     'd_shell_basis_slices_zero_based': {'spin1': [4, 9], 'spin2': [13, 18]},
                     'correction_tokens': [line[9:13] for line in local_lines],
                     'matrix_18x18_complex': [[[float(vin_full[i][j].real), float(vin_full[i][j].imag)] for j in range(NORB)] for i in range(NORB)],
                     'max_abs_delta_from_trace_v_hub_type2_Pin_H_d_blocks_by_spin': errors},
        'host_separation': {'LOCAL_HALL_host_types': sorted({int(v['ptype']) for (a,b,c),v in parsed.items() if a=='LOCAL_HALL' and b==BUILD_ID and c!=LOCAL_SITE}),
                            'host_records_are_not_candidates_for_recursive_Vin': True,
                            'generic_RB3_hub_matrix_default_type1_must_not_select_RD2_Vin': True},
        'capture_completeness': {'backend': 'RS', 'green_site': 0, 'green_type': LOCAL_TYPE,
                                 'energy_points': len(energy_indices), 'matrix_dimension': NORB,
                                 'full_matrices_complete': True, 'P_OUT_P_M1_shells_present': True},
        'driver': {'path': str(driver), 'sha256': sha(driver), 'source_path': str(driver_source),
                   'source_sha256': sha(driver_source), 'build_record': str(helper_record_path),
                   'build_record_sha256': sha(helper_record_path), 'compiled_source_manifest': str(helper_manifest_path),
                   'compiled_source_manifest_sha256': sha(helper_manifest_path),
                   'b1_head': helper_manifest['head'], 'b1_archive_sha256': helper_record['archive_sha256'],
                   'role': 'frozen radial replay driver; it has no impurity-site selector'},
        'math_binding': {'production_path': str(moment_math_path), 'production_sha256': sha(moment_math_path),
                         'python_formula_expression': 't=Re Tr(VP); c\'=c−t/q0; q2\'=q2−2 Re Tr(V(M1−cP))+Re Tr(V²P)−t²/q0',
                         'no_gate_arithmetic_executed': True},
        'fit_control': {'path': str(fit_path), 'sha256': sha(fit_path),
                        'preregistration_path': str(fit_prereg_path),
                        'preregistration_sha256': sha(fit_prereg_path),
                        'fit_source_prereg_sha256_match': True,
                        'spin_fit_mapping': fit_spin_map,
                        'D0_RS_ry': d0_ry,
                        'scalar_plus_0_01_fit_response_by_spin_ry': scalar_fit_response_at_plus_001,
                        'scalar_plus_0_01_registered_allowance_ry': scalar_fit_allowance_ry,
                        'mapping_rule': 'copy only per-spin fitted coefficient vectors from historical fit fields into RB3 generic per-spin slots; historical site/type labels do not select current physical site/type or Vin',
                        'use': 'registered RB3 scalar-shift fitted allowance only; it does not select physical site/type or Vin'},
        'explicit_limitations': ['No radial replay driver was invoked.', 'No carrier ratio, q2-prime gate, scalar-shift response, sign control, or spin-swap response was computed.',
                                 'This record is not an RD2 numerical result or GATE verdict.']
    }
    out = OUT / 'vin_binding.json'
    out.write_text(json.dumps(operand, indent=2, sort_keys=True) + '\n')
    print(f'PASS structural RD2 binding only; no gate/replay executed')
    print(f'Vin binding: {out} sha256={sha(out)}')
    print(f'type2 Pin_H cross-check max abs by spin: {errors}')


if __name__ == '__main__':
    main()
