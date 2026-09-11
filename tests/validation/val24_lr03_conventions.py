#!/usr/bin/env python3
"""Material-independent LR-03 convention algebra checks.

This is deliberately independent of RS-LMTO response production code and
ground-state/material data.  It checks only the matrix identities documented
in docs/LR_TDDFT_CONVENTIONS.md.
"""

from __future__ import annotations

import cmath
import math


Matrix = tuple[tuple[complex, ...], ...]


def mat_add(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(x + y for x, y in zip(row_a, row_b)) for row_a, row_b in zip(a, b))


def mat_scale(c: complex, a: Matrix) -> Matrix:
    return tuple(tuple(c * x for x in row) for row in a)


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0])))
        for i in range(len(a))
    )


def mat_dagger(a: Matrix) -> Matrix:
    return tuple(tuple(a[j][i].conjugate() for j in range(len(a))) for i in range(len(a)))


def max_error(a: Matrix, b: Matrix) -> float:
    return max(abs(x - y) for row_a, row_b in zip(a, b) for x, y in zip(row_a, row_b))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    zero = 0.0 + 0.0j
    one = 1.0 + 0.0j
    i = 0.0 + 1.0j
    identity: Matrix = ((one, zero), (zero, one))
    sx: Matrix = ((zero, one), (one, zero))
    sy: Matrix = ((zero, -i), (i, zero))
    sz: Matrix = ((one, zero), (zero, -one))
    sp: Matrix = mat_scale(0.5, mat_add(sx, mat_scale(i, sy)))
    sm: Matrix = mat_scale(0.5, mat_add(sx, mat_scale(-i, sy)))
    sigma_plus_expected: Matrix = ((zero, one), (zero, zero))
    sigma_minus_expected: Matrix = ((zero, zero), (one, zero))

    commutator = mat_add(mat_mul(sx, sy), mat_scale(-1.0, mat_mul(sy, sx)))
    require(max_error(commutator, mat_scale(2.0j, sz)) < 1.0e-15,
            "[sigma_x,sigma_y] != 2 i sigma_z")
    require(max_error(sp, sigma_plus_expected) < 1.0e-15,
            "sigma-plus matrix is not the documented halved ladder matrix")
    require(max_error(sm, sigma_minus_expected) < 1.0e-15,
            "sigma-minus matrix is not the documented halved ladder matrix")

    sigma_cap_plus = mat_add(sx, mat_scale(i, sy))
    sigma_cap_minus = mat_add(sx, mat_scale(-i, sy))
    require(max_error(sigma_cap_plus, mat_scale(2.0, sp)) < 1.0e-15,
            "Sigma-plus != 2 sigma-plus")
    require(max_error(sigma_cap_minus, mat_scale(2.0, sm)) < 1.0e-15,
            "Sigma-minus != 2 sigma-minus")

    spin_reversal: Matrix = mat_scale(1j, sy)
    plus_reversed = mat_mul(mat_mul(spin_reversal, sp), mat_dagger(spin_reversal))
    minus_reversed = mat_mul(mat_mul(spin_reversal, sm), mat_dagger(spin_reversal))
    require(max_error(plus_reversed, mat_scale(-1.0, sm)) < 1.0e-15,
            "spin reversal does not map sigma-plus to -sigma-minus")
    require(max_error(minus_reversed, mat_scale(-1.0, sp)) < 1.0e-15,
            "spin reversal does not map sigma-minus to -sigma-plus")
    require(max_error(mat_mul(spin_reversal, mat_mul(sz, mat_dagger(spin_reversal))),
                      mat_scale(-1.0, sz)) < 1.0e-15,
            "spin reversal does not reverse sigma-z")

    # Under the same reversal, xx is unchanged while xy changes sign.  This
    # is the Cartesian statement that the circular response channels exchange.
    chi_xx = 0.7 - 0.2j
    chi_xy = -0.15 + 0.45j
    chi_minus_original = chi_xx + 1j * chi_xy
    chi_plus_reversed = chi_xx - 1j * (-chi_xy)
    require(abs(chi_plus_reversed - chi_minus_original) < 1.0e-15,
            "spin reversal does not exchange the circular susceptibility channels")

    v_up = 1.375
    v_down = -0.625
    v0 = 0.5 * (v_up + v_down)
    b_sigma = 0.5 * (v_up - v_down)
    reconstructed = mat_add(mat_scale(v0, identity), mat_scale(b_sigma, sz))
    target = ((complex(v_up), zero), (zero, complex(v_down)))
    require(max_error(reconstructed, target) < 1.0e-15,
            "scalar/Pauli XC decomposition does not reconstruct both channels")

    # Keep the import visibly standard-library-only and make accidental unused
    # replacement of the exact identity less likely during future edits.
    require(math.isclose(abs(cmath.exp(1j * math.pi)), 1.0, rel_tol=0.0, abs_tol=1.0e-15),
            "standard-library complex arithmetic sanity check failed")
    print("val24_lr03_conventions: PASS (algebraic conventions only)")


if __name__ == "__main__":
    main()
