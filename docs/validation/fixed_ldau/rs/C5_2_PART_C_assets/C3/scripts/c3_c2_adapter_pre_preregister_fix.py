#!/usr/bin/env python3
"""C5.2 C3-only frozen local-versus-spectral replay adapter.

Private preparation harness. It never runs the solver. The RB3 adapter is
left untouched; this profile binds a C2-accepted Chebyshev endpoint and
reports only the C3 carrier-ratio and q2-prime scientific gates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

REPO = Path('/Users/andersb/Jobb/rslmto_fable_v4b')
VALIDATION = REPO / 'tests/validation/fixed_ldau'
sys.path.insert(0, str(VALIDATION))
import rb3_capture_adapter as rb3  # shared CSV schema and generic hash helpers only
import rs_a3_m03_runner as runner
import f90nml
from fortran_namelist_schema import source_schema, validate_input

LOCK_SCHEMA = 'rslmto-c5-c3-c2-accepted-endpoint-v1'
LOCK_STATUS = 'C2_ACCEPTED_CHEBYSHEV_U2_LOCAL_RS'
TRACE_SCHEMA = 'rslmto-fixed-u-trace-v1'
REQUIRED_KINDS = {
    'basis_map', 'green_function', 'radial_replay_mesh', 'radial_replay_v',
    'radial_replay_rho_in', 'radial_replay_ec_seed', 'radial_replay_eigen_seed',
    'radial_replay_ql', 'radial_replay_pl', 'pl_map_input', 'm0_carrier_input',
    'm0_carrier_local', 'm0_carrier_control', 'ldm', 'ldm_m1', 'v_hub',
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text())
    require(isinstance(obj, dict), f'expected JSON object: {path}')
    return obj


def validate_c2_lock(lock: dict[str, Any], lock_path: Path) -> dict[str, Any]:
    require(lock.get('schema') == LOCK_SCHEMA, 'C3 requires a C2 accepted-endpoint lock')
    require(lock.get('status') == LOCK_STATUS, 'C2 endpoint has not been accepted for C3')
    system = lock.get('system', {})
    require(system.get('backend') == 'RS' and system.get('calctype') == 'B', 'C3 endpoint must be RS bulk')
    require(system.get('recurrence') == 'chebyshev', 'C3 endpoint must use Chebyshev recursion')
    require(system.get('radial_moments') == 'local', 'C3 endpoint must be local mode')
    require(int(system.get('nsp', 0)) == 1 and int(system.get('lmax', 0)) == 2,
            'C3 endpoint must be collinear d-shell Fe')
    require(float(system.get('u_ev', float('nan'))) == 2.0 and float(system.get('j_ev', float('nan'))) == 0.0,
            'C3 endpoint must be U=2 eV, J=0')

    acceptance = lock.get('acceptance', {})
    comparison_status = acceptance.get('c2_gate_status')
    require(isinstance(comparison_status, str) and comparison_status.strip(),
            'C2 lock must record the block/Chebyshev comparison status descriptively')
    require(acceptance.get('c1_rs_status') == 'ACCEPTED', 'C2 endpoint lacks accepted RS convergence evidence')
    conv = acceptance.get('convergence', {})
    require(float(conv.get('outer_diff', math.inf)) <= 3.0e-8, 'accepted endpoint misses RS outer-Diff criterion')
    occ = conv.get('last20_d_occupation_by_spin', {})
    require(set(occ) == {'1', '2'} and all(len(occ[s]) == 20 for s in ('1', '2')),
            'accepted endpoint must bind 20 d-occupation samples per spin')
    spreads = {s: max(map(float, occ[s])) - min(map(float, occ[s])) for s in ('1', '2')}
    require(max(spreads.values()) <= 5.0e-7, 'accepted endpoint misses C1 d-occupation spread criterion')
    moments = conv.get('last20_sphere_moment_muB', [])
    require(len(moments) == 20, 'accepted endpoint must bind 20 physical sphere-moment samples')
    moment_spread = max(map(float, moments)) - min(map(float, moments))
    require(moment_spread <= 1.0e-6, 'accepted endpoint misses C1 sphere-moment spread criterion')

    endpoint = lock.get('accepted_endpoint', {})
    require(endpoint.get('run_id') and int(endpoint.get('accepted_iteration', 0)) > 0,
            'C2 lock lacks accepted run/iteration')
    for key in ('deck_path', 'checkpoint_path'):
        path = Path(endpoint.get(key, '')).resolve()
        require(path.is_file(), f'missing accepted endpoint {key}')
        require(endpoint.get(key.replace('_path', '_sha256')) == sha(path),
                f'accepted endpoint {key} hash mismatch')
        if key == 'deck_path':
            validate_input(path, source_schema(REPO))
            deck = f90nml.read(str(path))
            require('self' in deck and deck['self'].get('use_kspace') is False,
                    'accepted C2 deck must set use_kspace=.false. in &self')
            control = deck.get('control', {})
            require(str(control.get('calctype', '')).strip().lower() == 'b',
                    'accepted C2 deck must set calctype=B in &control')
            require(int(control.get('nsp', 0)) == 1,
                    'accepted C2 deck must set nsp=1 in &control')
            require(str(control.get('recur', '')).strip().lower() in {'cheb', 'chebyshev'},
                    'accepted C2 deck must select Chebyshev recurrence in &control')
            hamiltonian = deck.get('hamiltonian', {})
            require(str(hamiltonian.get('hubbard_radial_moments', '')).strip().lower() == 'local',
                    'accepted C2 deck must set local radial moments in &hamiltonian')
    executable = lock.get('executable', {})
    for key in ('binary_path', 'source_path'):
        path = Path(executable.get(key, '')).resolve()
        require(path.is_file(), f'missing C2 executable identity {key}')
        require(executable.get(key.replace('_path', '_sha256')) == sha(path),
                f'C2 executable identity hash mismatch for {key}')
    return {
        'lock_path': str(lock_path.resolve()), 'lock_sha256': sha(lock_path),
        'endpoint': endpoint, 'system': system, 'executable': executable,
        'acceptance': {'c1_occ_spreads': spreads, 'sphere_moment_spread_muB': moment_spread,
                       'outer_diff': float(conv['outer_diff']),
                       'block_chebyshev_comparison_status': comparison_status},
    }


def validate_full_green_rows(rows: list[dict[str, str]]) -> dict[str, Any]:
    """Require a complete square g0 matrix at every traced site/energy point."""
    maps: dict[int, dict[int, tuple[int, int, int]]] = {}
    groups: dict[tuple[int, int], dict[tuple[int, int], int]] = {}
    for row in rows:
        if row['kind'] == 'basis_map' and row['stage'] == 'RS_GREEN':
            site, idx = int(row['site']), int(row['row'])
            require(idx not in maps.setdefault(site, {}), f'duplicate basis-map row site={site} index={idx}')
            maps[site][idx] = (int(row['l']), int(row['col']), int(row['spin']))
        elif row['kind'] == 'green_function' and row['stage'] == 'RS_GREEN':
            site, energy = int(row['site']), int(row['k'])
            pair = (int(row['row']), int(row['col']))
            bucket = groups.setdefault((site, energy), {})
            require(pair not in bucket, f'duplicate Green matrix row site={site} energy={energy} pair={pair}')
            bucket[pair] = 1
    require(maps and groups, 'C3 capture must include full RS Green basis map and g0 rows')
    sizes = {site: len(mapping) for site, mapping in maps.items()}
    for site, mapping in maps.items():
        n = len(mapping)
        require(set(mapping) == set(range(n)), f'non-contiguous basis map for site={site}')
        lmax = max(label[0] for label in mapping.values())
        expected_labels = {(ell, m, spin) for ell in range(lmax + 1)
                           for m in range(-ell, ell + 1) for spin in (0, 1)}
        actual_labels = list(mapping.values())
        require(len(set(actual_labels)) == n,
                f'duplicate (l,m,spin) basis label for site={site}')
        require(set(actual_labels) == expected_labels,
                f'incomplete (l,m,spin) basis labels for site={site}')
    for (site, energy), pairs in groups.items():
        require(site in sizes, f'Green rows have no basis map for site={site}')
        n = sizes[site]
        expected = {(i, j) for i in range(n) for j in range(n)}
        require(set(pairs) == expected,
                f'incomplete Green matrix site={site} energy={energy}: got {len(pairs)} expected {n*n}')
    sites = set(maps)
    energies = {site: {energy for site0, energy in groups if site0 == site} for site in sites}
    require(all(energies[site] for site in sites), 'basis-mapped site has no Green energy grid')
    require(len({tuple(sorted(v)) for v in energies.values()}) == 1,
            'Green energy grids differ across traced sites')
    return {'basis_dimension_by_site': sizes,
            'energy_grid_points_by_site': {str(site): len(energies[site]) for site in sorted(sites)},
            'full_green_matrix_count': len(groups),
            'all_green_matrices_complete': True}


def validate_capture(capture_dir: Path, lock_binding: dict[str, Any], source: Path,
                     binary: Path, driver: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    capture_dir = capture_dir.resolve()
    csv_path = capture_dir / 'records.csv'
    manifest_path = capture_dir / 'manifest.json'
    if not manifest_path.is_file():
        manifest_path = capture_dir / 'manifest.final.json'
    require(csv_path.is_file() and manifest_path.is_file(), 'capture requires records.csv and manifest')
    manifest = load_json(manifest_path)
    csv_sha, manifest_sha = sha(csv_path), sha(manifest_path)
    require(manifest.get('schema') == TRACE_SCHEMA, 'unexpected fixed-U trace schema')
    require(manifest.get('backend') == 'RS', 'C3 capture backend must be RS')
    require(manifest.get('csv_sha256') == csv_sha, 'trace manifest does not bind CSV')
    require(manifest.get('executable_sha256') == sha(binary), 'trace executable hash mismatch')
    require(manifest.get('source_sha') == sha(source), 'trace source hash mismatch')
    capture = manifest.get('capture_execution', {})
    # Historical trace manifests store this data in the run fields instead.
    run_id = capture.get('run_id', manifest.get('run_id'))
    iteration = int(capture.get('iteration', manifest.get('iteration', -1)))
    require(run_id and iteration > 0, 'capture manifest lacks run id/iteration')
    input_hashes = manifest.get('input_sha256', {})
    endpoint = lock_binding['endpoint']
    require(endpoint['checkpoint_sha256'] in input_hashes.values(),
            'capture manifest is not bound to the accepted C2 checkpoint')
    capture_input_sha = capture.get('input_sha256')
    if capture_input_sha is not None:
        require(capture_input_sha in input_hashes.values(), 'capture manifest does not bind continuation deck')

    with csv_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == rb3.CSV_FIELDS, 'unexpected C3 trace CSV header')
        rows = list(reader)
    require(rows, 'empty C3 capture')
    counts: Counter[tuple[str, str, int, int, int]] = Counter()
    runs, iterations = set(), set()
    kinds = Counter()
    for n, row in enumerate(rows, 1):
        runs.add(row['run_id']); iterations.add(int(row['iteration'])); kinds[row['kind']] += 1
        for key in ('lambda', 'weight', 'energy', 'scalar', 're', 'im'):
            require(math.isfinite(float(row[key])), f'nonfinite {key} in row {n}')
        counts[(row['kind'], row['stage'], int(row['site']), int(row['type']), int(row['spin']))] += 1
    require(runs == {run_id} and iterations == {iteration}, 'mixed run/iteration in C3 trace')
    require(REQUIRED_KINDS <= set(kinds), f'missing trace kinds: {sorted(REQUIRED_KINDS-set(kinds))}')
    green_shape = validate_full_green_rows(rows)
    require(sum(1 for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0) == 495,
            'unsupported C3 radial mesh extent')
    require(sum(1 for r in rows if r['kind']=='radial_replay_ec_seed' and r['stage']=='RADIAL_REPLAY') == 10,
            'unsupported C3 core seed extent')
    require(sum(1 for r in rows if r['kind']=='radial_replay_eigen_seed' and r['stage']=='RADIAL_REPLAY') == 6,
            'unsupported C3 valence seed extent')
    for spin in (1, 2):
        for kind, stage, expected in (
            ('m0_carrier_input','FROZEN_REPLAY',1),
            ('m0_carrier_local','FROZEN_REPLAY',1),
            ('m0_carrier_control','VOUT_CONTROL',1),
            ('m0_carrier_control','P_TRANSPOSE_CONTROL',1),
        ):
            got = counts[(kind, stage, 0, 0, spin)]
            require(got == expected, f'expected {expected} {kind}/{stage} rows for spin {spin}; got {got}')
    require(int(manifest.get('iteration', iteration)) == iteration, 'CSV/manifest iteration mismatch')
    return rows, {
        'capture_dir': str(capture_dir), 'csv_path': str(csv_path), 'csv_sha256': csv_sha,
        'manifest_path': str(manifest_path), 'manifest_sha256': manifest_sha,
        'run_id': run_id, 'iteration': iteration, 'row_count': len(rows),
        'kind_counts': dict(sorted(kinds.items())), 'green_shape': green_shape, 'source_sha256': sha(source),
        'binary_sha256': sha(binary), 'driver_sha256': sha(driver),
        'accepted_checkpoint_sha256': endpoint['checkpoint_sha256'],
    }


def _index_rows(rows: list[dict[str, str]]) -> None:
    runner._BY_KIND_STAGE = {}
    for row in rows:
        runner._BY_KIND_STAGE.setdefault((row['kind'], row['stage']), []).append(row)


def prepare_c3(rows: list[dict[str, str]], nr: int, ncore: int, nval: int) -> dict[str, Any]:
    """Create exactly raw/local replay cases without A2.2 fit or A2.1 checks."""
    data = runner.prepare(rows, nr, ncore, nval, 'baseline')
    ql_raw = data['cases'][0][3].copy()
    vin, p, m1, q0, center, q2, carrier, q2prime = [], [], [], [], [], [], [], []
    observed_q2prime = []
    for spin in (1, 2):
        vhub = runner.hub_matrix(rows, 'Pin_H')
        offset = 4 + 9 * (spin - 1)
        v = vhub[offset:offset + 5, offset:offset + 5]
        density = runner.matrix(rows, 'ldm', 'P_OUT', spin, 2, 5)
        first = runner.matrix(rows, 'ldm_m1', 'M1_OUT', spin, 2, 5)
        operand = runner.one(rows, 'm0_carrier_input', 'FROZEN_REPLAY', spin=spin, ell=2)
        q = float(operand['scalar']); c = float(operand['energy']); q2_raw = float(operand['re'])
        t = float((v @ density).trace().real)
        q2p = q2_raw - 2.0 * float((v @ (first - c * density)).trace().real) \
              + float((v @ v @ density).trace().real) - t * t / q
        local = runner.one(rows, 'm0_carrier_local', 'FROZEN_REPLAY', spin=spin, ell=2)
        # Capture's precomputed M0 local value is retained as an independent RECORD;
        # no tolerance is imposed between independently rounded text operands.
        observed_q2prime.append(float(local['re']))
        vin.append(v); p.append(density); m1.append(first); q0.append(q)
        center.append(c); q2.append(q2_raw); carrier.append(t); q2prime.append(q2p)
        ql_raw[2, 2, spin - 1] = q2p
    shifts = [carrier[i] / q0[i] if q0[i] != 0.0 else float('nan') for i in range(2)]
    require(all(math.isfinite(x) for x in shifts), 'zero/nonfinite q0 prevents C3 local replay')
    ql_local = ql_raw.copy()
    cases = [
        ('RAW_PL', [0.0, 0.0], True, ql_raw.copy()),
        ('LOCAL_PL', shifts, True, ql_local),
    ]
    data['cases'] = cases
    data['target'] = {'vin': vin, 'p': p, 'm1': m1, 'q0': q0, 'center': center, 'q2': q2,
                      'carrier_by_spin_ry': carrier, 'q2prime_by_spin_ry2': q2prime,
                      'capture_q2prime_by_spin_ry2': observed_q2prime,
                      'local_shift_by_spin_ry': shifts}
    return data


def summarize_c3(data: dict[str, Any], cases: dict[int, dict[str, Any]]) -> dict[str, Any]:
    require(set(cases) == {1, 2}, 'C3 replay must return exactly RAW_PL and LOCAL_PL cases')
    per_spin = {}
    for spin in (1, 2):
        delta = float(cases[1]['sumev'][spin] - cases[2]['sumev'][spin])
        carrier = float(data['target']['carrier_by_spin_ry'][spin - 1])
        ratio = delta / carrier if carrier != 0.0 else None
        ratio_pass = bool(ratio is not None and math.isfinite(ratio) and 0.99 <= ratio <= 1.01)
        per_spin[str(spin)] = {'measured_raw_minus_local_ry': delta,
                               'captured_trace_vin_p_ry': carrier,
                               'ratio': ratio, 'gate_pass': ratio_pass}
    q2 = {str(spin): float(data['target']['q2prime_by_spin_ry2'][spin - 1]) for spin in (1, 2)}
    ratio_pass = all(v['gate_pass'] for v in per_spin.values())
    q2_pass = all(math.isfinite(v) and v >= 0.0 for v in q2.values())
    gates = {'carrier_ratio_each_spin_within_1_percent': ratio_pass,
             'q2prime_each_spin_nonnegative': q2_pass}
    gates['all_c3_gates_pass'] = all(gates.values())
    return {'carrier_ratio_by_spin': per_spin,
            'q2prime_by_spin_ry2': q2,
            'capture_reported_q2prime_by_spin_ry2': {
                str(spin): float(data['target']['capture_q2prime_by_spin_ry2'][spin - 1]) for spin in (1, 2)},
            'gates': gates}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capture', type=Path, required=True)
    ap.add_argument('--lock', type=Path, required=True, help='accepted C2 C3 endpoint lock')
    ap.add_argument('--source', type=Path, required=True, help='locked source/self.f90')
    ap.add_argument('--binary', type=Path, required=True, help='frozen C2 production executable')
    ap.add_argument('--driver', type=Path, required=True, help='source-connected replay driver')
    ap.add_argument('--mode', choices=('baseline','preregister','cases'), required=True)
    ap.add_argument('--baseline-report', type=Path)
    ap.add_argument('--preregistration', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    lock_path = args.lock.resolve(); source = args.source.resolve(); binary = args.binary.resolve(); driver = args.driver.resolve()
    lock_binding = validate_c2_lock(load_json(lock_path), lock_path)
    require(lock_binding['executable']['binary_sha256'] == sha(binary), 'binary does not match accepted C2 production lock')
    require(lock_binding['executable']['source_sha256'] == sha(source), 'source does not match accepted C2 production lock')
    rows, capture = validate_capture(args.capture, lock_binding, source, binary, driver)
    _index_rows(rows)
    nr = sum(1 for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0)
    ncore = sum(1 for r in rows if r['kind']=='radial_replay_ec_seed' and r['stage']=='RADIAL_REPLAY')
    nval = sum(1 for r in rows if r['kind']=='radial_replay_eigen_seed' and r['stage']=='RADIAL_REPLAY')
    require((nr, ncore, nval) == (495, 10, 6), f'unsupported radial replay dimensions {(nr,ncore,nval)}')
    output_dir = args.output.resolve().parent; output_dir.mkdir(parents=True, exist_ok=True)
    baseline = None
    if args.mode in ('preregister','cases'):
        require(args.baseline_report is not None, 'preregister/cases require --baseline-report')
        baseline = load_json(args.baseline_report)
        require(baseline.get('status') == 'PASS_LIVE_STATE_BASELINE' and
                baseline.get('capture_csv_sha256') == capture['csv_sha256'] and
                baseline.get('manifest_sha256') == capture['manifest_sha256'] and
                baseline.get('driver_sha256') == capture['driver_sha256'],
                'LIVE_STATE baseline is not passing or not bound to this capture/driver')
    if args.mode == 'baseline':
        data = runner.prepare(rows, nr, ncore, nval, 'baseline')
        raw, input_sha, output_sha = runner.run_persistent(driver, runner.text_input(data, rows, nr, ncore, nval), output_dir)
        result = runner.output_cases(raw)
        require(set(result) == {1}, 'baseline must contain one LIVE_STATE case')
        checks = runner.baseline_check(rows, result[1])
        raw_path = args.output.with_suffix('.csv'); raw_path.write_bytes(raw)
        report = {'schema':'rslmto-c5-c3-live-baseline-v1',
                  'status':'PASS_LIVE_STATE_BASELINE' if checks['pass'] else 'FAIL_LIVE_STATE_BASELINE',
                  'capture':capture,'c2_lock_sha256':lock_binding['lock_sha256'],
                  'block_chebyshev_comparison_status':lock_binding['acceptance']['block_chebyshev_comparison_status'],
                  'driver_sha256':sha(driver),'baseline_replay_integrity_checks':checks,
                  'driver_input_sha256':input_sha,'driver_output_path':str(raw_path),
                  'driver_output_sha256':sha(raw_path),'driver_output_file_sha256':output_sha}
    elif args.mode == 'preregister':
        data = prepare_c3(rows, nr, ncore, nval)
        report = {'schema':'rslmto-c5-c3-c2-preregistration-v1','status':'LOCKED_BEFORE_C3_CASES',
                  'capture':capture,'c2_lock_sha256':lock_binding['lock_sha256'],
                  'block_chebyshev_comparison_status':lock_binding['acceptance']['block_chebyshev_comparison_status'],
                  'baseline_report_sha256':sha(args.baseline_report),'driver_sha256':sha(driver),
                  'case_names':[x[0] for x in data['cases']],
                  'target_operands':data['target'],
                  'scientific_gates':['carrier_ratio_each_spin_within_1_percent','q2prime_each_spin_nonnegative'],
                  'notice':'Only the two C5.2 C3 scientific gates are preregistered; no A2.1/A2.2 fitted controls are imported.'}
    else:
        require(args.preregistration is not None, 'cases requires --preregistration')
        prereg = load_json(args.preregistration)
        require(prereg.get('status') == 'LOCKED_BEFORE_C3_CASES' and
                prereg.get('capture',{}).get('csv_sha256') == capture['csv_sha256'] and
                prereg.get('capture',{}).get('manifest_sha256') == capture['manifest_sha256'] and
                prereg.get('baseline_report_sha256') == sha(args.baseline_report) and
                prereg.get('driver_sha256') == sha(driver),
                'C3 preregistration does not bind capture, baseline, and driver')
        data = prepare_c3(rows, nr, ncore, nval)
        raw, input_sha, output_sha = runner.run_persistent(driver, runner.text_input(data, rows, nr, ncore, nval), output_dir)
        result = runner.output_cases(raw)
        measurements = summarize_c3(data, result)
        raw_path = args.output.with_suffix('.csv'); raw_path.write_bytes(raw)
        report = {'schema':'rslmto-c5-c3-measured-cases-v1','status':'MEASURED_C3_CASES',
                  'capture':capture,'c2_lock_sha256':lock_binding['lock_sha256'],
                  'block_chebyshev_comparison_status':lock_binding['acceptance']['block_chebyshev_comparison_status'],
                  'baseline_report_sha256':sha(args.baseline_report),
                  'preregistration_sha256':sha(args.preregistration),'driver_sha256':sha(driver),
                  'measurements':measurements,'driver_input_sha256':input_sha,
                  'driver_output_path':str(raw_path),'driver_output_sha256':sha(raw_path),
                  'driver_output_file_sha256':output_sha}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, default=rb3.json_default) + '\n')
    print(json.dumps({'mode':args.mode,'status':report['status'],'output':str(args.output),
                      'output_sha256':sha(args.output),'capture_csv_sha256':capture['csv_sha256'],
                      'driver_sha256':sha(driver)},sort_keys=True))

if __name__ == '__main__':
    main()
