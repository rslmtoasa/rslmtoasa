#!/usr/bin/env python3
"""Plot the seven-column jij.out written by source/exchange.f90; retain directed rows."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_exchange(path):
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith(('#', '!')):
            continue
        fields = line.replace('D', 'E').replace('d', 'e').split()
        if len(fields) != 7:
            raise ValueError(f'{path}:{number}: expected 7 columns from jij.out')
        rows.append([float(value) for value in fields])
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or not np.isfinite(data).all():
        raise ValueError(f'{path}: need finite exchange rows')
    kinds = data[:, :2]
    if np.any(kinds < 1) or np.any(kinds != np.floor(kinds)):
        raise ValueError(f'{path}: endpoint types must be positive integers')
    if np.any(data[:, 6] < 0) or not np.allclose(
            np.linalg.norm(data[:, 2:5], axis=1), data[:, 6], rtol=0, atol=2e-6):
        raise ValueError(f'{path}: distance disagrees with displacement (tolerance 2e-6 a)')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('jij', type=Path)
    parser.add_argument('--types', type=int, nargs='+',
                        help='Keep rows with BOTH endpoints in these verified symbolic types')
    parser.add_argument('--output', type=Path, default=Path('exchange.png'))
    parser.add_argument('--csv', type=Path, help='Save selected original rows with labelled units')
    args = parser.parse_args()
    if args.types is not None and min(args.types) < 1:
        parser.error('--types must be positive symbolic-type IDs')
    paths = [args.output] + ([args.csv] if args.csv else [])
    if len({path.resolve() for path in paths}) != len(paths):
        parser.error('figure and CSV need different paths')
    if any(path.exists() for path in paths):
        parser.error('output already exists; choose fresh figure/CSV paths')
    try:
        data = read_exchange(args.jij)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    selected = data if args.types is None else data[
        np.isin(data[:, 0], args.types) & np.isin(data[:, 1], args.types)]
    print(f'Selected {len(selected)} of {len(data)} directed rows; excluded {len(data)-len(selected)}.')
    if not len(selected):
        parser.error('no rows remain; check symbolic types against the electronic input')
    fig, ax = plt.subplots(figsize=(7, 4))
    for first, second in np.unique(selected[:, :2], axis=0):
        group = selected[(selected[:, 0] == first) & (selected[:, 1] == second)]
        ax.scatter(group[:, 6], group[:, 5] ,
                   label=f'type {int(first)} → {int(second)}', alpha=0.7)
        print(f'Types {int(first)} -> {int(second)}: {len(group)} rows; '
              f'J = [{group[:, 5].min():.8g}, {group[:, 5].max():.8g}] mRy')
    ax.axhline(0, color='black', linewidth=0.6)
    ax.set(xlabel='Pair distance / a', ylabel='J (mRy)', title='Selected directed exchange rows')
    ax.legend()
    fig.tight_layout()
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    if args.csv:
        np.savetxt(args.csv, selected, delimiter=',', fmt='%.12g', comments='',
                   header='type_i,type_j,dx_over_a,dy_over_a,dz_over_a,J_mRy,distance_over_a')
    print(f'Saved {args.output}. No shell averaging, symmetry expansion, or bond deduplication.')
    print('Types are not cluster-site IDs; verify site mapping and Hamiltonian conventions separately.')


if __name__ == '__main__':
    main()
