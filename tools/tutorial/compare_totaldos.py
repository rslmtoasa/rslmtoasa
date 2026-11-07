#!/usr/bin/env python3
"""Overlay two RS-LMTO totaldos.out files without fitting or normalizing them."""
from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RY_TO_EV = 13.605693122994


def read_dos(path: Path) -> tuple[np.ndarray, np.ndarray, float | None]:
    text = path.read_text()
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(('#', '!')):
            continue
        fields = line.replace('D', 'E').replace('d', 'e').split()
        if len(fields) < 2:
            raise ValueError(f'{path}: expected at least two columns')
        rows.append((float(fields[0]), float(fields[1])))
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[0] < 2 or not np.isfinite(data).all():
        raise ValueError(f'{path}: need at least two finite data rows')
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f'{path}: energy grid must be strictly increasing')
    match = re.search(r'^#\s*Fermi level:\s*([-+0-9.EeDd]+)', text, re.M)
    fermi = float(match[1].replace('D', 'E').replace('d', 'e')) if match else None
    return data[:, 0], data[:, 1], fermi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rs', type=Path)
    parser.add_argument('ks', type=Path)
    parser.add_argument('--rs-fermi', type=float, help='Ry; value used to shift the written RS DOS')
    parser.add_argument('--ks-fermi', type=float, help='Ry; overrides a k-space DOS header')
    parser.add_argument('--rs-energy-min', type=float, help='Ry; known absolute first energy-grid point (alternative to Fermi)')
    parser.add_argument('--ks-energy-min', type=float, help='Ry; known absolute first energy-grid point (alternative to Fermi)')
    parser.add_argument('--alignment', choices=['absolute', 'fermi'], default='absolute')
    parser.add_argument('--sites-per-cell', type=int, default=2,
                        help='Same number of physical sites in both DOS totals (Si primitive cell: 2)')
    parser.add_argument('--unit', choices=['eV', 'Ry'], default='eV')
    parser.add_argument('--output', type=Path, default=Path('si-dos-comparison.png'))
    args = parser.parse_args()
    if args.sites_per_cell < 1:
        parser.error('--sites-per-cell must be positive')
    er, dr, header_r = read_dos(args.rs)
    ek, dk, header_k = read_dos(args.ks)
    fr = args.rs_fermi if args.rs_fermi is not None else header_r
    fk = args.ks_fermi if args.ks_fermi is not None else header_k
    if args.rs_fermi is not None and args.rs_energy_min is not None:
        parser.error('use either --rs-fermi or --rs-energy-min, not both')
    if args.ks_fermi is not None and args.ks_energy_min is not None:
        parser.error('use either --ks-fermi or --ks-energy-min, not both')
    if args.alignment == 'absolute':
        if args.rs_energy_min is not None:
            fr = args.rs_energy_min - er[0]
        if args.ks_energy_min is not None:
            fk = args.ks_energy_min - ek[0]
        if fr is None or fk is None:
            parser.error('absolute alignment needs the Fermi shift or known absolute first grid point for each file (Ry)')
        er, ek = er + fr, ek + fk
        xlabel = 'Absolute energy'
    else:
        xlabel = 'E - E_F (separate Fermi zero for each route)'
    integrate = getattr(np, 'trapezoid', None)
    if integrate is None:
        integrate = np.trapz
    for label, energy, dos in [('RS', er, dr), ('KS', ek, dk)]:
        print(f'{label}: sampled-window integral = {integrate(dos, energy):.8g} states/cell; '
              f'minimum DOS = {dos.min():.6g} states/(Ry cell)')
    lo, hi = max(er[0], ek[0]), min(er[-1], ek[-1])
    if lo >= hi:
        parser.error('energy windows do not overlap with the selected alignment')
    common = np.linspace(lo, hi, 2001)
    rs, ks = np.interp(common, er, dr), np.interp(common, ek, dk)
    rms = np.sqrt(np.mean((rs - ks)**2))
    scale = np.sqrt(np.mean(ks**2))
    print(f'Common-window RMS difference = {rms:.8g} states/(Ry cell)')
    print(f'RMS relative to KS RMS DOS = {rms / scale:.8g}' if scale > 0 else 'KS RMS DOS is zero')
    print('No renormalization, smoothing, fitted shifts, or pass/fail tolerance applied.')
    factor = RY_TO_EV if args.unit == 'eV' else 1.0
    fig, (ax, diff) = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
                                  gridspec_kw={'height_ratios': [3, 1]})
    for energy, dos, label in [(er, dr, 'RS Chebyshev'), (ek, dk, 'k-space Gaussian')]:
        ax.plot(energy * factor, dos / (factor * args.sites_per_cell), label=label)
    ax.set_ylabel(f'DOS (states/({args.unit} atom))')
    ax.legend()
    diff.plot(common * factor, (rs - ks) / (factor * args.sites_per_cell))
    diff.axhline(0, color='black', linewidth=0.6)
    diff.set_ylabel('RS - KS')
    diff.set_xlabel(f'{xlabel} ({args.unit})')
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
