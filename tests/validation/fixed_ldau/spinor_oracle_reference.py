#!/usr/bin/env python3
"""Dump R08 oracle values (spinor_functional_oracle.py, unchanged) for the Fortran unit test.

Usage: spinor_oracle_reference.py OUTPUT. Each case: header, P (re, im), E_int E_DC E_Hub, V (re, im);
matrices are written row-major as 'row col re im'.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spinor_functional_oracle as r08

CASES = ((2.3, 0.0, 71), (2.3, 0.41, 72), (1.4, 0.22, 73))   # U, J, seed


def write_matrix(lines: list[str], m: np.ndarray) -> None:
    for i in range(r08.NDIM):
        for k in range(r08.NDIM):
            lines.append(f"{i+1} {k+1} {m[i, k].real:.17e} {m[i, k].imag:.17e}")


def main() -> None:
    lines = [f"{len(CASES)}"]
    for u, j, seed in CASES:
        f2 = 14*j/1.625
        t = r08.coulomb_tensor(u, f2, 0.625*f2)
        p = r08.random_density(np.random.default_rng(seed))
        eint, edc, ehub, v = r08.energy_potential(p, t, u, j)
        lines.append(f"{u!r} {j!r}")
        write_matrix(lines, p)
        lines.append(f"{eint:.17e} {edc:.17e} {ehub:.17e}")
        write_matrix(lines, v)
    Path(sys.argv[1]).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
