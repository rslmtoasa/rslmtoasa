#!/usr/bin/env python3
"""Plot completed SCF residuals from a solver log; default: run.log."""
import argparse
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def read_residuals(path):
    text = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', path.read_text(errors='replace'))
    # self.f90 writes these once per completed step, including the converged step.
    pattern = r'(Not converged!\s*Diff=|Converged!)\s*(\S+)'
    records = re.findall(pattern, text)
    if not records:
        raise ValueError('no completed SCF residuals found; check the log path and solver progress')
    residuals = np.array([float(value.replace('D', 'E').replace('d', 'e')) for _, value in records])
    if not np.isfinite(residuals).all() or np.any(residuals < 0):
        raise ValueError('SCF residuals must be finite and nonnegative')
    return residuals, records[-1][0] == 'Converged!'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', type=Path, nargs='?', default=Path('run.log'))
    args = parser.parse_args()
    try:
        residuals, converged = read_residuals(args.log)
    except (OSError, ValueError) as error:
        parser.error(f'{args.log}: {error}')
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(np.arange(1, len(residuals)+1), residuals, 'o-', markersize=3)
    # f12.10 residual output rounds to zero below 5e-11; keep zeros visible.
    ax.set_yscale('symlog', linthresh=5e-11)
    ax.set_ylim(0, max(float(residuals.max())*2, 1e-10))
    ax.set_xlim(0.5, len(residuals)+0.5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    if len(residuals) == 1:
        ax.set_xticks([1])
    ax.set(xlabel='Completed SCF step (log order)', ylabel='Printed SCF mixing residual',
           title=f'{args.log.name}: ' + ('last step reports convergence' if converged else 'no convergence reported at last step'))
    ax.grid(alpha=0.3)
    fig.tight_layout()
    target = args.log.with_name(args.log.stem + '-convergence.png')
    fig.savefig(target, dpi=160)
    plt.close(fig)
    print(f'Saved {target}; {len(residuals)} completed steps; last printed residual {residuals[-1]:.10g}')
    print('Linear axis near zero (|residual| < 5e-11); printed zero is limited by log precision.')
    print('Log order joins appended restarts; SCF convergence does not establish numerical convergence.')


if __name__ == '__main__':
    main()
