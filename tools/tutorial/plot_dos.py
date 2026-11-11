#!/usr/bin/env python3
"""Plot all tutorial DOS files in a run directory (default: current directory)."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RY_TO_EV = 13.605693122994


def read_dos(path):
    text = path.read_text().replace('D', 'E').replace('d', 'e')
    from io import StringIO
    data = np.loadtxt(StringIO(text), comments=('#', '!'), ndmin=2)
    if len(data) < 2 or not np.isfinite(data).all():
        raise ValueError('need at least two finite rows')
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError('energy grid must increase strictly')
    return data


def curves(path, data):
    """Return panels of (label, signed DOS); down-spin values are mirrored."""
    ncol = data.shape[1]
    if path.name == 'totaldos.out' and ncol == 2:
        return [[('total', data[:, 1])]]
    if path.name in ('magneticdos.out', 'magneticdos') and ncol == 7:
        return [[('up', data[:, 1]), ('down (mirrored)', -data[:, 2])]]
    if path.name.endswith('_orbital_dos.out') and ncol in (9, 19, 33):
        norb = (ncol - 1) // 2
        # basis.f90: norb=(lmax+1)**2; bands.f90 writes up block then down block.
        return [[(f'{shell} orbital {m-l*l+1} ↑', data[:, 1+m]),
                 (f'{shell} orbital {m-l*l+1} ↓', -data[:, 1+norb+m])]
                for l, shell in enumerate('spdf') if (l+1)**2 <= norb
                for m in range(l*l, (l+1)**2)]
    if path.name.endswith('_dos.out') and ncol == 4:
        return [[('up', data[:, 2]), ('down (mirrored)', -data[:, 3])]]
    raise ValueError(f'unsupported format ({ncol} columns); expected total=2, magnetic=7, local=4, orbital=9/19/33')


def plot_file(path, output):
    data = read_dos(path)
    panels = curves(path, data)
    # A separate panel per orbital avoids hiding small contributions under large ones.
    fig, axes = plt.subplots(len(panels), 1, sharex=True, squeeze=False,
                             figsize=(8, max(4, 1.8*len(panels))))
    for ax, lines in zip(axes[:, 0], panels):
        for index, (label, dos) in enumerate(lines):
            ax.plot(data[:, 0]*RY_TO_EV, dos/RY_TO_EV,
                    label=label, linestyle='--' if index == 1 else '-')
        ax.axhline(0, color='black', linewidth=0.5)
        ax.axvline(0, color='grey', linewidth=0.5)
        ax.set_ylabel('states/eV')
        ax.legend(loc='upper right', fontsize='small')
    axes[-1, 0].set_xlabel('Written E − E_F (eV)')
    fig.suptitle(path.name + ' · system DOS' if path.name in
                 ('totaldos.out', 'magneticdos.out', 'magneticdos') else path.name + ' · local DOS')
    fig.tight_layout()
    output.mkdir(parents=True, exist_ok=True)
    target = output / (path.name + '.png')
    fig.savefig(target, dpi=160)
    plt.close(fig)
    print(f'Saved {target}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, nargs='?', default=Path('.'))
    args = parser.parse_args()
    paths = sorted({p for pattern in ('totaldos.out', 'magneticdos.out', 'magneticdos', '*_dos.out')
                    for p in args.directory.glob(pattern) if p.is_file()})
    if not paths:
        parser.error(f'{args.directory}: no supported DOS files found')
    failed = False
    for path in paths:
        try:
            plot_file(path, args.directory / 'plots')
        except (OSError, ValueError) as error:
            print(f'{path}: {error}')
            failed = True
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
