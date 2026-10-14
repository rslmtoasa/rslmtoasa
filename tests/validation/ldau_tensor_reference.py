#!/usr/bin/env python3
"""Independent complex-harmonic Coulomb tensor reference for Phase 1.

The reference uses SymPy Wigner 3j symbols and normalized complex spherical
harmonics. It does not import or call production math_mod routines.
"""
from __future__ import annotations

import argparse
import cmath
import math
import subprocess
import numpy as np
from sympy.physics.wigner import wigner_3j


def gaunt(l1: int, m1: int, k: int, q: int, l2: int, m2: int) -> float:
    if abs(m1) > l1 or abs(q) > k or abs(m2) > l2:
        return 0.0
    return (-1) ** m1 * ordinary_gaunt(l1, k, l2, -m1, q, m2)


def ordinary_gaunt(l1: int, l2: int, l3: int, m1: int, m2: int, m3: int) -> float:
    if m1 + m2 + m3 != 0:
        return 0.0
    value = math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * math.pi))
    return value * float(wigner_3j(l1, l2, l3, 0, 0, 0)) * float(wigner_3j(l1, l2, l3, m1, m2, m3))


def tensor(l: int, a: int, b: int, c: int, d: int, slater: dict[int, float]) -> complex:
    """T_abcd = integral Y*_a Y_b(1) v Y*_c Y_d(2) = <a c|b d>."""
    total = 0.0
    for k, fk in slater.items():
        if k > 2 * l or k % 2:
            continue
        for q in range(-k, k + 1):
            total += 4 * math.pi * fk / (2 * k + 1) * gaunt(l, a, k, q, l, b) * gaunt(l, d, k, q, l, c)
    return complex(total)


def legacy_gaunt_ord(l1: int, l2: int, l3: int, u1: int, u2: int, u3: int) -> float:
    root2 = 1 / math.sqrt(2)
    if u1 > 0 and u2 > 0 and u3 > 0:
        if u2 + u3 == u1:
            return root2 * (-1) ** u1 * ordinary_gaunt(l1, l2, l3, -u1, u2, u3)
        if u1 + u3 == u2:
            return root2 * (-1) ** u2 * ordinary_gaunt(l1, l2, l3, u1, -u2, u3)
        if u1 + u2 == u3:
            return root2 * (-1) ** u3 * ordinary_gaunt(l1, l2, l3, u1, u2, -u3)
        return 0.0
    if u1 > 0 and u2 < 0 and u3 < 0:
        if u1 + u3 == u2:
            return root2 * (-1) ** u3 * ordinary_gaunt(l1, l2, l3, u1, -u2, u3)
        if u1 + u2 == u3:
            return root2 * (-1) ** u2 * ordinary_gaunt(l1, l2, l3, u1, u2, -u3)
        if u1 + u2 + u3 == 0:
            return root2 * (-1) ** (u1 + 1) * ordinary_gaunt(l1, l2, l3, u1, u2, u3)
        return 0.0
    if u3 == 0 and u1 == u2 and u1 > 0:
        return (-1) ** u1 * ordinary_gaunt(l1, l2, l3, u1, -u1, 0)
    if u1 == 0 and u2 == u3 and u3 < 0:
        return (-1) ** u2 * ordinary_gaunt(l1, l2, l3, 0, u2, -u2)
    if u1 == 0 and u2 == 0 and u3 == 0:
        return ordinary_gaunt(l1, l2, l3, 0, 0, 0)
    return 0.0


def real_gaunt(l1: int, l2: int, l3: int, u1: int, u2: int, u3: int) -> float:
    if u1 >= u2 and u2 >= u3:
        return legacy_gaunt_ord(l1, l2, l3, u1, u2, u3)
    if u1 >= u3 and u3 >= u2:
        return legacy_gaunt_ord(l1, l3, l2, u1, u3, u2)
    if u2 >= u1 and u1 >= u3:
        return legacy_gaunt_ord(l2, l1, l3, u2, u1, u3)
    if u2 >= u3 and u3 >= u1:
        return legacy_gaunt_ord(l2, l3, l1, u2, u3, u1)
    if u3 >= u1 and u1 >= u2:
        return legacy_gaunt_ord(l3, l1, l2, u3, u1, u2)
    return legacy_gaunt_ord(l3, l2, l1, u3, u2, u1)


def legacy_real_tensor(l: int, a: int, b: int, c: int, d: int, slater: dict[int, float]) -> float:
    total = 0.0
    for k, fk in slater.items():
        if k > 2 * l or k % 2:
            continue
        for q in range(-k, k + 1):
            total += 4 * math.pi * fk / (2 * k + 1) * real_gaunt(l, k, l, a, q, b) * real_gaunt(l, k, l, c, q, d)
    return total


def hcpx_cart_to_spherical(l: int) -> np.ndarray:
    """The l block of `hcpx('cart2sph')`, transcribed from its coefficients."""
    q = 1 / math.sqrt(2)
    if l == 1:
        # cart rows (x,y,z), spherical columns (m=-1,0,+1)
        return np.array([[q, 0, -q], [1j*q, 0, 1j*q], [0, 1, 0]], dtype=complex)
    if l == 2:
        # cart rows (xy,yz,zx,x2-y2,3z2-r2), spherical columns m=-2..+2
        return np.array([
            [1j*q, 0, 0, 0, -1j*q],
            [0, 1j*q, 0, 1j*q, 0],
            [0, q, 0, -q, 0],
            [q, 0, 0, 0, q],
            [0, 0, 1, 0, 0],
        ], dtype=complex)
    raise ValueError(l)


def real_gaunt_to_cart_order(l: int) -> np.ndarray:
    """Reorder realgaunt's m=-l..l harmonics into hcpx Cartesian order."""
    if l == 1:
        # real-gaunt m order (-1,0,+1) is (y,z,x); hcpx rows are (x,y,z)
        return np.array([1, 2, 0])
    # real-gaunt m order (-2,-1,0,+1,+2) is (xy,yz,z2,zx,x2-y2)
    return np.array([0, 1, 4, 2, 3])


def tensor_arrays(l: int, slater: dict[int, float]) -> tuple[np.ndarray, np.ndarray]:
    mvals = range(-l, l + 1)
    shape = (2*l + 1,) * 4
    old = np.empty(shape, dtype=complex)
    reference = np.empty(shape, dtype=complex)
    for ia, a in enumerate(mvals):
        for ib, b in enumerate(mvals):
            for ic, c in enumerate(mvals):
                for id_, d in enumerate(mvals):
                    old[ia, ib, ic, id_] = legacy_real_tensor(l, a, b, c, d, slater)
                    reference[ia, ib, ic, id_] = tensor(l, a, b, c, d, slater)
    return old, reference


def transform_tensor_to_spherical(real_tensor: np.ndarray, l: int) -> np.ndarray:
    # Source realGaunt order -> hcpx Cartesian order -> spherical basis.
    v = hcpx_cart_to_spherical(l)
    permutation = real_gaunt_to_cart_order(l)
    v = v[permutation, :]
    return np.einsum("ia,jb,kc,ld,ijkl->abcd", v.conj(), v, v.conj(), v, real_tensor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", help="optional report path")
    parser.add_argument("--fortran-executable", help="UnitLdaUTensor path; compare its complete p/d dump")
    args = parser.parse_args()
    records: list[str] = ["Independent normalized-complex-harmonic tensor reference (SymPy 1.14)"]
    slaters = {1: {0: 1.0, 2: 0.37, 4: 0.19}, 2: {0: 1.0, 2: 0.37, 4: 0.19}}
    for l, fs in slaters.items():
        old_array, reference = tensor_arrays(l, fs)
        max_mismatch = 0.0
        max_element = None
        for a in range(-l, l + 1):
            for b in range(-l, l + 1):
                for c in range(-l, l + 1):
                    for d in range(-l, l + 1):
                        old = legacy_real_tensor(l, a, b, c, d, fs)
                        new = tensor(l, a, b, c, d, fs).real
                        delta = abs(old - new)
                        if delta > max_mismatch:
                            max_mismatch, max_element = delta, (a, b, c, d, old, new)
        transformed = transform_tensor_to_spherical(old_array, l)
        transform_error = float(np.max(np.abs(transformed - reference)))
        records.append(f"l={l}: direct same-index legacy-real vs complex-spherical maximum difference {max_mismatch:.12g} at {max_element}")
        records.append(f"l={l}: legacy tensor transformed by hcpx orbital map vs independent complex reference max residual {transform_error:.12g}")
    # Elements that distinguish the magnetic selection rule and conjugation.
    records.append(f"d T[-2,-1,-2,-1] = {tensor(2, -2, -1, -2, -1, slaters[2]).real:.12g} (must vanish: -2+-2 != -1+-1)")
    records.append(f"d T[-2,-1,-1,-2] = {tensor(2, -2, -1, -1, -2, slaters[2]).real:.12g} (allowed: -2+-1 = -1+-2)")
    if args.fortran_executable:
        completed = subprocess.run([args.fortran_executable, "--dump"], check=True, capture_output=True, text=True)
        max_residual = 0.0
        count = 0
        for line in completed.stdout.splitlines():
            parts = line.split()
            if not parts or parts[0] != "TENSOR":
                continue
            channel, l, a, b, c, d = map(int, parts[1:7])
            observed = complex(float(parts[7]), float(parts[8]))
            channel_fs = {channel: slaters[l][channel]}
            expected = tensor(l, a, b, c, d, channel_fs)
            max_residual = max(max_residual, abs(observed - expected))
            count += 1
        if count != 2118:
            raise RuntimeError(f"Expected 2118 Fortran tensor elements (706 per Slater channel); received {count}")
        if max_residual > 2e-12:
            raise RuntimeError(f"Fortran vs SymPy tensor residual {max_residual:.4e} exceeds 2e-12")
        records.append(f"Fortran vs independent SymPy p/d F0/F2/F4 tensor grids: {count} elements, max residual {max_residual:.12g} (tol 2e-12)")
    content = "\n".join(records) + "\n"
    print(content, end="")
    if args.report:
        with open(args.report, "w", encoding="utf-8") as stream:
            stream.write(content)


if __name__ == "__main__":
    main()
