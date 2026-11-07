#!/usr/bin/env python3
"""Copy the two-site Si fixture into independent 150-moment RS and Gaussian KS lessons."""
from pathlib import Path
import argparse
import shutil
import f90nml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('runs/alternative-si'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source = root / 'tests/scf/cases/bulk/diamondSi'
    output = args.output.resolve()
    if output.exists():
        parser.error(f'{output} already exists; choose a new directory to preserve results')
    output.mkdir(parents=True)
    for route in ('rs', 'ks'):
        target = output / route
        shutil.copytree(source, target, ignore=shutil.ignore_patterns("README.md"))
        changes = {
            'lattice': {'pbc': True, 'morton_sfc': False, 'n1': 12, 'n2': 12, 'n3': 12,
                        'strux_backend': 'strux_lib'},
            'self': {'nstep': 25, 'conv_thr': 5e-8, 'cold': True,
                     'use_kspace': route == 'ks', 'force_nonmagnetic': [True, True]},
            'control': {'nsp': 1, 'lmax': 1, 'lld': 150, 'recur': 'chebyshev',
                        'cheb_backend': 'legacy'},
            'energy': {'energy_min': -1.5, 'energy_max': 1.0, 'channels_ldos': 1000},
            'mix': {'beta': 0.1, 'mixtype': 'broyden'},
            'hamiltonian': {'hoh': False},
        }
        if route == 'ks':
            changes['reciprocal'] = {
                'nk1': 8, 'nk2': 8, 'nk3': 8, 'use_symmetry_reduction': False,
                'use_time_reversal': True, 'use_shift': False, 'n_energy_points': 1000,
                'dos_method': 'gaussian', 'gaussian_sigma': 0.02,
                'kspace_ham_order': 'first',
            }
        patched = target / 'input.patched.nml'
        f90nml.patch(target / 'input.nml', changes, patched)
        patched.replace(target / 'input.nml')
        shutil.copy2(target / 'input.nml', target / 'input.initial.nml')
        print(f'Prepared {target}')
    print('Teaching sp-basis inputs; check SCF and numerical convergence separately.')


if __name__ == '__main__':
    main()
