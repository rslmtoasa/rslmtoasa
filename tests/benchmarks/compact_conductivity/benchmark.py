#!/usr/bin/env python3
"""Compare frozen, fully periodic Pt conductivity with unchanged main.

Requires f90nml, NumPy and GNU time. Paths are resolved before each run.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess

import f90nml
import numpy as np

PBC = dict.fromkeys(('pbc', 'b1', 'b2', 'b3'), True)
TOL = 5e-6
RAW_ZERO_SCALE = 1e-10
RAW_ZERO_TOL = 1e-12
VARIANTS = ('spin', 'orbital', 'charge', 'hoh_spin', 'hoh_orbital', 'hoh_charge',
            'accumulation', 'orbital_pair', 'random_spin', 'multitype', 'lehmann_charge', 'dyson_charge')
ENV = dict(OMP_DYNAMIC='FALSE', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           BLIS_NUM_THREADS='1', RSLMTO_KPM_RANDOM_SEED='917', LC_ALL='C')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(fixture, directory, size, order, channels, variant):
    directory.mkdir(parents=True)
    if variant == 'multitype':
        potentials = fixture.parents[3] / 'scf/cases/impurity/B2FeCo'
        for symbol in ('Fe', 'Co'):
            shutil.copyfile(potentials / f'{symbol}.nml', directory / f'{symbol}.nml')
    else:
        shutil.copyfile(fixture / 'Pt.nml', directory / 'Pt.nml')
    nml = f90nml.read(fixture / 'input.nml')
    nml['lattice'].update(PBC, n1=size, n2=size, n3=size, strux_backend='legacy')
    nml['energy']['channels_ldos'] = channels
    observable = variant.removeprefix('hoh_').removeprefix('random_')
    if variant in ('lehmann_charge', 'dyson_charge'):
        route = variant.removesuffix('_charge')
        observable = 'charge'
        nml['calculation']['gf_route'] = route
        nml['reciprocal'] = dict(nk1=size, nk2=size, nk3=size, use_symmetry_reduction=False)
    observable = {'accumulation': 'orbital_accumulation', 'orbital_pair': 'orbital_accumulation',
                  'multitype': 'charge'}.get(observable, observable)
    nml['control'].update(cond_ll=order, cheb_backend='legacy', transport_backend='legacy',
                          linear_in='charge', linear_out=observable, cond_calctype='per_type')
    if variant == 'orbital_pair':
        nml['control']['linear_in'] = 'orbital_accumulation'
    nml['hamiltonian']['hoh'] = variant.startswith('hoh_')
    if observable == 'charge':
        nml['hamiltonian']['v_alpha'] = [1, 0, 0]
    if variant == 'random_spin':
        nml['control'].update(cond_calctype='random_vec', random_vec_num=3)
    if variant == 'multitype':
        # B2 FeCo with the repository's frozen Fe and Co starting potentials.
        nml['lattice'].update(crystal_sym='b2', ntype=2, alat=2.859, wav=1.406216,
                              ct=[3.0, 3.0], r2=9.0, rc=20)
        nml['atoms']['label'] = ['Fe', 'Co']
        nml['energy'].update(fermi=-0.041385, energy_min=-2.5, energy_max=1.5)
    elif size <= 6:
        nml['lattice'].update(rc=20, ct=[3.0], r2=9.0)
    nml.write(directory / 'input.nml')


def run(binary, fixture, directory, args, order, variant):
    prepare(fixture, directory, args.size, order, args.channels, variant)
    inp = directory / 'input.nml'
    if any(f90nml.read(inp)['lattice'].get(key) is not True for key in PBC):
        raise ValueError(f'Open boundary in {inp}')
    env = {**ENV, 'OMP_NUM_THREADS': str(args.threads)}
    command = [args.time, '-v', '-o', 'time.txt', str(binary)]
    manifest = dict(command=command, env=env, input_sha256=sha(inp),
                    potential_sha256={pot.name: sha(pot) for pot in directory.glob('*.nml') if pot != inp},
                    binary_sha256=sha(binary))
    (directory / 'run.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with (directory / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=directory, env={**os.environ, **env},
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0 or not (directory / 'cond_total.out').is_file():
        raise RuntimeError(f'Run failed ({result.returncode}): {directory / "run.log"}')
    log = (directory / 'run.log').read_text(errors='replace')
    if variant in ('lehmann_charge', 'dyson_charge'):
        if f"gf_route={variant.removesuffix('_charge')} -- filling" not in log:
            raise RuntimeError(f'Reciprocal moment backend not confirmed: {directory / "run.log"}')
    elif not re.search(r'KPM_PROFILE.*moment_backend=\s*cpu_legacy', log):
        raise RuntimeError(f'Legacy moment backend not confirmed: {directory / "run.log"}')
    timings = (directory / 'time.txt').read_text()
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', timings)[1])
    elapsed = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', timings)[1]
    seconds = sum(float(value) * 60**i for i, value in enumerate(reversed(elapsed.split(':'))))
    manifest['geometry_sha256'] = {name: sha(directory / name) for name in ('fort.805', 'map')
                                   if (directory / name).is_file()}
    with (directory / 'clust').open() as cluster:
        manifest['atom_count'] = int(cluster.readline().split('=')[1])
    expected_atoms = args.size**3 * (2 if variant == 'multitype' else 1)
    if manifest['atom_count'] != expected_atoms:
        raise ValueError(f'Expected {expected_atoms} sites in {directory}, found {manifest["atom_count"]}')
    manifest['kpm_profile'] = next((line.strip() for line in log.splitlines() if 'KPM_PROFILE ' in line), '')
    (directory / 'run.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return rss, seconds


def outputs(directory):
    return {p.name for p in directory.glob('*.out')
            if p.name.startswith(('cond_', 'response_')) or '_cond' in p.name or p.name.endswith('_response.out')} | {'fort.123'}


def compare(base, feature):
    left, right = (json.loads((path / 'run.json').read_text()) for path in (base, feature))
    for field in ('input_sha256', 'potential_sha256', 'geometry_sha256', 'atom_count', 'env'):
        if left[field] != right[field]:
            raise ValueError(f'Paired runs differ in {field}: {base} and {feature}')
    if outputs(base) != outputs(feature):
        raise ValueError(f'Output sets differ: {base} and {feature}')
    if sha(base / 'response_units.txt') != sha(feature / 'response_units.txt'):
        raise ValueError(f'Response unit metadata differs: {base} and {feature}')
    meta = (base / 'response_units.txt').read_text()
    conversion = float(re.search(r'^output conversion =\s*(\S+)', meta, re.MULTILINE)[1])
    weights = np.fromstring(re.search(r'^trace weights =\s*(.+)', meta, re.MULTILINE)[1], sep=' ')
    labels = f90nml.read(base / 'input.nml')['atoms']['label']
    symbols = [f90nml.read(base / f'{label}.nml')['element']['symbol'] for label in labels]
    type_weight = dict(zip(symbols, weights))
    for name in sorted(outputs(base)):
        before, after = (np.loadtxt(path / name, ndmin=2) for path in (base, feature))
        if before.shape != after.shape or not np.array_equal(before[:, 0], after[:, 0]):
            raise ValueError(f'Energy grid or array shape differs: {name}')
        if not np.isfinite(before).all() or not np.isfinite(after).all():
            raise ValueError(f'Nonfinite response: {name}')
        # Identify symmetry-null channels in raw units so that a large unit
        # conversion alone (e.g. ~10^18 for orbital density) cannot make
        # roundoff look like a nonzero response. Keep the printed-output
        # tolerance for nonzero channels and use an FP64 absolute test at zero.
        unit_scale = 1.0
        if name.startswith('response_'):
            unit_scale = abs(conversion) * np.sum(np.abs(weights))
        elif name.endswith('_response.out'):
            unit_scale = abs(conversion * type_weight[name.removesuffix('_response.out')])
        diff = after[:, 1:] - before[:, 1:]
        fermi = np.argmin(np.abs(before[:, 0]))
        orbital_imag = name.endswith(('_orb_im.out', '_imag.out'))
        if orbital_imag:
            real_name = name.replace('_orb_im.out', '_orb_real.out').replace('_imag.out', '_real.out')
            baseline_real, optimized_real = (np.max(np.abs(np.loadtxt(path / real_name)[:, 1:]), axis=0)
                                             for path in (base, feature))
        for j in range(diff.shape[1]):
            absolute = np.max(np.abs(diff[:, j]))
            scale = np.max(np.abs(before[:, j + 1]))
            is_imag = orbital_imag or (before.shape[1] == 3 and j == 1)
            if is_imag and before.shape[1] == 3:
                real_scale = max(np.max(np.abs(before[:, 1])), np.finfo(float).tiny)
                optimized_real_scale = max(np.max(np.abs(after[:, 1])), np.finfo(float).tiny)
            elif orbital_imag:
                real_scale, optimized_real_scale = baseline_real[j], optimized_real[j]
            else:
                real_scale = optimized_real_scale = None
            raw_scale = scale / unit_scale
            near_zero = raw_scale <= RAW_ZERO_SCALE
            normalized = absolute / unit_scale if near_zero else absolute / max(1, scale)
            threshold = RAW_ZERO_TOL if near_zero else TOL
            yield dict(file=name, column=j + 2, max_abs=absolute, reference_scale=scale,
                       raw_scale=raw_scale, normalized=normalized, tolerance=threshold,
                       comparison='absolute_raw' if near_zero else 'scaled_output',
                       scaled_difference=absolute / max(1, scale), rms=np.sqrt(np.mean(diff[:, j]**2)),
                       fermi_energy=before[fermi, 0], baseline_fermi=before[fermi, j + 1],
                       optimized_fermi=after[fermi, j + 1],
                       baseline_max_imag=scale if is_imag else '',
                       optimized_max_imag=np.max(np.abs(after[:, j + 1])) if is_imag else '',
                       baseline_imag_ratio=scale / real_scale if real_scale else '',
                       optimized_imag_ratio=np.max(np.abs(after[:, j + 1])) / optimized_real_scale if optimized_real_scale else '',
                       passed=normalized <= threshold)


def check_traces(directory, variant):
    raw_total = np.loadtxt(directory / 'cond_total.out')
    physical_total = np.loadtxt(directory / 'response_total.out')
    expected_types = {'multitype': {'Fe', 'Co'}, 'random_spin': set()}.get(variant, {'Pt'})
    for ending in ('_cond.out', '_cond_orb_real.out', '_cond_orb_im.out', '_response.out'):
        names = {file.name.removesuffix(ending) for file in directory.glob(f'*{ending}')}
        if names != expected_types:
            raise ValueError(f'Expected type outputs for {expected_types}, got {names}: {directory}')
    for total, pattern, divisor in (
        (raw_total, '*_cond.out', max(1, len(expected_types))),
        (physical_total, '*_response.out', 1),
    ):
        traces = [np.loadtxt(file)[:, 1:] for file in sorted(directory.glob(pattern))]
        if traces:
            residual = np.max(np.abs(np.sum(traces, axis=0) / divisor - total[:, 1:]))
            if residual > TOL * max(1, np.max(np.abs(total[:, 1:]))):
                raise ValueError(f'Per-type trace sum disagrees with total: {directory}')
    for total, stem in ((raw_total, 'cond_total_orb'), (physical_total, 'response_orbitals')):
        suffixes = ('real', 'im' if stem.startswith('cond_') else 'imag')
        for column, suffix in enumerate(suffixes, start=1):
            orbital = np.loadtxt(directory / f'{stem}_{suffix}.out')[:, 1:]
            residual = np.max(np.abs(np.sum(orbital, axis=1) - total[:, column]))
            scale = max(1, np.max(np.abs(orbital)), np.max(np.abs(total[:, column])))
            if residual > TOL * scale:
                raise ValueError(f'Orbital trace sum disagrees with total: {directory}')
    # The fully traced, per-type Pt charge and Hall conductivity is real.
    # Seeded finite stochastic estimators and symmetry-forbidden local
    # accumulations still get the ordinary numerical comparison.
    if variant in ('spin', 'orbital', 'charge', 'hoh_spin', 'hoh_orbital', 'hoh_charge'):
        real = np.max(np.abs(raw_total[:, 1]))
        imag = np.max(np.abs(raw_total[:, 2]))
        if real <= 1e-10 or imag / real > 1e-12:
            raise ValueError(f'Unexpected closed-Pt conductivity: {directory}, Re={real:g}, Im={imag:g}')


def write_csv(path, rows):
    with path.open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--optimized', type=Path)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--orders', nargs='+', type=int, default=[50, 100, 250])
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=['spin', 'orbital', 'charge'])
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--size', type=int, default=12)
    parser.add_argument('--channels', type=int, default=2500)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--time', default=shutil.which('time') or '/usr/bin/time')
    parser.add_argument('--memory-gib', type=float, required=True, help='RAM limit available to the unmodified baseline in GiB')
    args = parser.parse_args()
    if min(args.orders) < 2 or min(args.size, args.channels, args.threads, args.repetitions) < 1:
        parser.error('positive sizes, thread counts and repetitions and cond_ll >= 2 are required')
    args.output.mkdir(parents=True, exist_ok=True)
    binaries = {'baseline': args.baseline.resolve()}
    if args.optimized:
        binaries['optimized'] = args.optimized.resolve()
    performance, parity = [], []
    for order in args.orders:
        for variant in args.variants:
            ntrace = {'multitype': 2, 'random_spin': 3}.get(variant, 1)
            site_count = args.size**3 * (2 if variant == 'multitype' else 1)
            gamma_phase = 16 * order**2 * (args.channels + 10 + 18**2*ntrace)
            moment_phase = 16 * 18**2 * order * (site_count + order*ntrace)
            required = max(gamma_phase, moment_phase)
            if required >= args.memory_gib * 2**30:
                row = dict(order=order, variant=variant, baseline_rss_kib='', optimized_rss_kib='',
                           baseline_wall_s='', optimized_wall_s='', saving_pct='', saved_gib='', time_ratio='',
                           status='NOT RUN (OOM risk)', baseline_lower_bound_gib=required / 2**30)
                performance.append(row)
                write_csv(args.output / 'benchmark.csv', performance)
                print(f'{variant} M={order}: NOT RUN (OOM risk), allocation lower bound '
                      f'{required / 2**30:.3f} GiB', flush=True)
                continue
            samples = {label: [] for label in binaries}
            for rep in range(args.repetitions):
                case = args.output / f'{variant}_M{order}_r{rep + 1}'
                sequence = list(binaries) if rep % 2 == 0 else list(reversed(binaries))
                for label in sequence:
                    print(f'{case.name} {label}', flush=True)
                    samples[label].append(run(binaries[label], args.fixture.resolve(), case / label, args, order, variant))
                    check_traces(case / label, variant)
                    if label == 'optimized' and variant not in ('lehmann_charge', 'dyson_charge'):
                        profile = json.loads((case / label / 'run.json').read_text())['kpm_profile']
                        for key in ('bytes_gamma', 'bytes_mu_pack', 'host_full_moment_bytes'):
                            if int(re.search(rf'\b{key}=\s*(\d+)', profile)[1]) != 0:
                                raise ValueError(f'Unexpected full workspace ({key}): {case}')
                        for key in ('gamma_factor_workspace_bytes', 'host_compact_moment_bytes'):
                            if int(re.search(rf'\b{key}=\s*(\d+)', profile)[1]) <= 0:
                                raise ValueError(f'Missing compact workspace ({key}): {case}')
                if args.optimized:
                    rows = [dict(order=order, variant=variant, repetition=rep + 1, **row)
                            for row in compare(case / 'baseline', case / 'optimized')]
                    parity.extend(rows)
                    write_csv(args.output / 'numerical_parity.csv', parity)
                    if not all(row['passed'] for row in rows):
                        raise RuntimeError(f'Numerical gate failed: {case}')
            before_rss, before_time = [statistics.median(x) for x in zip(*samples['baseline'])]
            after_rss, after_time = ([statistics.median(x) for x in zip(*samples['optimized'])]
                                     if args.optimized else (None, None))
            performance.append(dict(order=order, variant=variant, baseline_rss_kib=before_rss,
                                    optimized_rss_kib=after_rss, baseline_wall_s=before_time,
                                    optimized_wall_s=after_time,
                                    saving_pct=100 * (1 - after_rss / before_rss) if after_rss else '',
                                    saved_gib=(before_rss - after_rss) / 2**20 if after_rss else '',
                                    time_ratio=after_time / before_time if after_time else '',
                                    status='PASS' if args.optimized else 'BASELINE',
                                    baseline_lower_bound_gib=required / 2**30))
            write_csv(args.output / 'benchmark.csv', performance)


if __name__ == '__main__':
    main()
