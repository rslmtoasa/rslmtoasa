#!/usr/bin/env python3
"""Independent SymPy reference for the Phase 3 d-shell energy fixture.

This uses SymPy Wigner-3j symbols and normalized complex spherical harmonics;
it does not import production Fortran or the Fortran tensor reference script.
Values are in the code's Ry energy units for U=1.4, Jbar=0.22.
"""

from math import pi, sqrt

import numpy as np
from sympy.physics.wigner import wigner_3j


def gaunt(l, m, k, q, mp):
    return ((-1) ** m * sqrt((2 * l + 1) * (2 * k + 1) * (2 * l + 1) / (4 * pi))
            * float(wigner_3j(l, k, l, 0, 0, 0))
            * float(wigner_3j(l, k, l, -m, q, mp)))


def tensor_d(u, jbar):
    l = 2
    f2 = 14 * jbar / 1.625
    fk = {0: u, 2: f2, 4: 0.625 * f2}
    mvals = range(-l, l + 1)
    t = np.zeros((5, 5, 5, 5), dtype=complex)
    for a, ma in enumerate(mvals):
        for b, mb in enumerate(mvals):
            for c, mc in enumerate(mvals):
                for d, md in enumerate(mvals):
                    for k, value in fk.items():
                        for q in range(-k, k + 1):
                            t[a, b, c, d] += value * 4 * pi / (2 * k + 1) * gaunt(l, ma, k, q, mb) * gaunt(l, md, k, q, mc)
    return t


def fixture_density():
    p = np.zeros((5, 5), dtype=complex)
    q = np.zeros((5, 5), dtype=complex)
    for a in range(5):
        p[a, a] = 0.13 + 0.09 * (a + 1)
        q[a, a] = 0.08 + 0.065 * (a + 1)
        for b in range(a + 1, 5):
            p[a, b] = 0.004 * (a + b + 2) + 0.003j * (a - b)
            p[b, a] = p[a, b].conjugate()
            q[a, b] = 0.003 * (a + b + 2) - 0.002j * (a - b)
            q[b, a] = q[a, b].conjugate()
    return p, q


def evaluate(t, p, q, u, jbar):
    densities = (p, q)
    direct = 0.0
    exchange = 0.0
    potential = [np.zeros((5, 5), dtype=complex) for _ in range(2)]
    for s in range(2):
        for sp in range(2):
            direct += np.einsum("abcd,ba,dc->", t, densities[s], densities[sp]).real
        exchange += np.einsum("abcd,da,bc->", t, densities[s], densities[s]).real
    nspin = [np.trace(x).real for x in densities]
    ntot = sum(nspin)
    eint = 0.5 * (direct - exchange)
    edc = 0.5 * u * ntot * (ntot - 1) - 0.5 * jbar * sum(n * (n - 1) for n in nspin)
    for s in range(2):
        other = 1 - s
        for a in range(5):
            for b in range(5):
                value = sum(t[a, b, c, d] * densities[other][d, c]
                            + (t[a, b, c, d] - t[a, d, c, b]) * densities[s][d, c]
                            for c in range(5) for d in range(5))
                if a == b:
                    value += -u * (ntot - 0.5) + jbar * (nspin[s] - 0.5)
                potential[s][a, b] = value
    trvp = sum(np.trace(v @ density).real for v, density in zip(potential, densities))
    return eint, edc, eint - edc, trvp, potential


if __name__ == "__main__":
    u, jbar = 1.4, 0.22
    p, q = fixture_density()
    values = evaluate(tensor_d(u, jbar), p, q, u, jbar)
    print("Eint={:.12f} EDC={:.12f} EHub={:.12f} TrVP={:.12f}".format(*values[:4]))
    print("Vup(1,2)={:.12f}{:+.12f}i".format(values[4][0][0, 1].real, values[4][0][0, 1].imag))
