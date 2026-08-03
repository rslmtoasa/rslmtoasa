#!/usr/bin/env python3
"""Independent algebraic acceptance tests for GBT WP1 / G1.

This file intentionally imports no production GBT or Pauli helper.  The
``pauli_pair_candidate`` route is a local transcription of the elementary
pair algebra; every result is checked against a separately assembled dense
spin-orbital matrix expression.  The convention under test is the blueprint's
physical-displacement convention:

    alpha = q_cart . (R + tau_b - tau_a)

where q_cart includes the 2*pi reciprocal factor.  Sublattice reference
azimuths live in the reference moments, never in alpha.
"""

from __future__ import annotations

import math
import sys

import numpy as np


TOL = 1.0e-12
I2 = np.eye(2, dtype=np.complex128)
SX = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SY = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
SZ = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    """Infinity-norm relative error with a stable zero-scale convention."""
    return float(
        np.linalg.norm(actual - expected, ord=np.inf)
        / max(1.0, float(np.linalg.norm(expected, ord=np.inf)))
    )


def rotz(alpha: float, moment: np.ndarray) -> np.ndarray:
    c, s = math.cos(alpha), math.sin(alpha)
    return np.array(
        [c * moment[0] - s * moment[1], s * moment[0] + c * moment[1], moment[2]],
        dtype=np.float64,
    )


def dspin(alpha: float) -> np.ndarray:
    """D(alpha) = exp(-i sigma_z alpha / 2), built explicitly as dense data."""
    return np.diag(
        [np.exp(-0.5j * alpha), np.exp(0.5j * alpha)]
    ).astype(np.complex128)


def spin_matrix(moment: np.ndarray) -> np.ndarray:
    """Explicit dense sigma . m; this is the independent oracle building block."""
    return moment[0] * SX + moment[1] * SY + moment[2] * SZ


def endpoint_matrix(w0: np.ndarray, w1: np.ndarray, moment: np.ndarray) -> np.ndarray:
    """Dense block-diagonal W = w0 I + w1 sigma.m for one endpoint species."""
    norb = len(w0)
    out = np.zeros((2 * norb, 2 * norb), dtype=np.complex128)
    sm = spin_matrix(moment)
    for orbital in range(norb):
        row = slice(2 * orbital, 2 * orbital + 2)
        out[row, row] = w0[orbital] * I2 + w1[orbital] * sm
    return out


def dense_pair_oracle(
    structure: np.ndarray,
    w0_i: np.ndarray,
    w1_i: np.ndarray,
    w0_j: np.ndarray,
    w1_j: np.ndarray,
    moment_i: np.ndarray,
    moment_j_ref: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Dense W_i [S tensor D(alpha)] W_j reference-frame oracle."""
    return (
        endpoint_matrix(w0_i, w1_i, moment_i)
        @ np.kron(structure, dspin(alpha))
        @ endpoint_matrix(w0_j, w1_j, moment_j_ref)
    )


def pauli_pair_candidate(
    structure: np.ndarray,
    w0_i: np.ndarray,
    w1_i: np.ndarray,
    w0_j: np.ndarray,
    w1_j: np.ndarray,
    moment_i: np.ndarray,
    moment_j_ref: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Independent local Pauli-component route used only as the tested side.

    It uses the elementary dot/cross pair expression and an explicit final
    right multiplication.  It deliberately does not call a production helper.
    """
    moment_j = rotz(alpha, moment_j_ref)
    dot = float(np.dot(moment_i, moment_j))
    cross = np.cross(moment_i, moment_j)
    norb = structure.shape[0]
    raw = np.zeros((2 * norb, 2 * norb), dtype=np.complex128)

    for left in range(norb):
        for right in range(norb):
            s = structure[left, right]
            h0 = w0_i[left] * s * w0_j[right] + w1_i[left] * s * w1_j[right] * dot
            hv = (
                w1_i[left] * s * w0_j[right] * moment_i
                + w0_i[left] * s * w1_j[right] * moment_j
                + 1.0j * w1_i[left] * s * w1_j[right] * cross
            )
            row = slice(2 * left, 2 * left + 2)
            col = slice(2 * right, 2 * right + 2)
            raw[row, col] = np.array(
                [[h0 + hv[2], hv[0] - 1.0j * hv[1]], [hv[0] + 1.0j * hv[1], h0 - hv[2]]],
                dtype=np.complex128,
            )
    return raw @ np.kron(np.eye(norb, dtype=np.complex128), dspin(alpha))


def reference_moment(theta: float, azimuth: float) -> np.ndarray:
    return np.array(
        [math.sin(theta) * math.cos(azimuth), math.sin(theta) * math.sin(azimuth), math.cos(theta)],
        dtype=np.float64,
    )


def test_cell_phase_equivalence() -> float:
    """Cartesian and fractional phases agree for non-equivalent cell metrics."""
    cells = {
        "cubic": np.diag([3.1, 3.1, 3.1]),
        "hexagonal": np.array([[2.7, -1.35, 0.0], [0.0, 2.7 * math.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, 4.4]]),
        "triclinic": np.array([[3.2, 0.7, -0.4], [0.2, 2.9, 0.6], [-0.1, 0.5, 4.1]]),
    }
    q_fracs = [np.array([0.17, -0.31, 0.23]), np.array([-0.42, 0.11, 0.37])]
    d_fracs = [np.array([1.0, -2.0, 3.0]), np.array([-0.35, 0.625, 0.125])]
    max_error = 0.0
    for cell in cells.values():
        for q_frac in q_fracs:
            # A maps direct fractional coordinates to Cartesian coordinates.
            q_cart = 2.0 * math.pi * np.linalg.solve(cell.T, q_frac)
            for d_frac in d_fracs:
                phase_frac = 2.0 * math.pi * float(q_frac @ d_frac)
                phase_cart = float(q_cart @ (cell @ d_frac))
                # Compare modulo 2*pi to make the convention explicit.
                max_error = max(max_error, float(abs(np.exp(1.0j * phase_cart) - np.exp(1.0j * phase_frac))))
    return max_error


def test_stoner_shifted_k_identity() -> float:
    """A one-orbital Stoner model has exactly epsilon(k-/+q/2) diagonals."""
    translations = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    hoppings = np.array([-0.43, -0.43, 0.19, 0.19, 0.07, 0.07], dtype=np.float64)
    k = np.array([0.61, -0.37, 0.0])
    q = np.array([0.29, 0.17, 0.0])
    theta, exchange, onsite = 0.73, 0.41, -0.12

    direct = onsite * I2 + exchange * (math.cos(theta) * SZ + math.sin(theta) * SX)
    for displacement, hopping in zip(translations, hoppings):
        alpha = float(q @ displacement)
        direct += hopping * np.exp(1.0j * float(k @ displacement)) * dspin(alpha)

    def epsilon(argument: np.ndarray) -> complex:
        return onsite + sum(
            hopping * np.exp(1.0j * float(argument @ displacement))
            for displacement, hopping in zip(translations, hoppings)
        )

    expected = np.array(
        [
            [epsilon(k - q / 2.0) + exchange * math.cos(theta), exchange * math.sin(theta)],
            [exchange * math.sin(theta), epsilon(k + q / 2.0) - exchange * math.cos(theta)],
        ],
        dtype=np.complex128,
    )
    return relerr(direct, expected)


def random_species_pair_cases() -> tuple[float, float]:
    """Random dense unequal-endpoint pair oracle and its reverse-bond check."""
    rng = np.random.default_rng(20260803)
    max_dense_error = 0.0
    max_reverse_error = 0.0
    for norb in (2, 3):
        for _ in range(12):
            structure = rng.normal(size=(norb, norb)) + 1.0j * rng.normal(size=(norb, norb))
            w0_a = rng.uniform(0.15, 0.85, size=norb)
            w1_a = rng.uniform(-0.55, 0.55, size=norb)
            # Deliberately distinct endpoint species: every orbital product is exposed.
            w0_b = rng.uniform(0.25, 1.05, size=norb)
            w1_b = rng.uniform(-0.75, 0.45, size=norb)
            ma = reference_moment(rng.uniform(0.2, 2.6), rng.uniform(-math.pi, math.pi))
            mb = reference_moment(rng.uniform(0.2, 2.6), rng.uniform(-math.pi, math.pi))
            alpha = float(rng.uniform(-2.8, 2.8))

            candidate = pauli_pair_candidate(structure, w0_a, w1_a, w0_b, w1_b, ma, mb, alpha)
            oracle = dense_pair_oracle(structure, w0_a, w1_a, w0_b, w1_b, ma, mb, alpha)
            max_dense_error = max(max_dense_error, relerr(candidate, oracle))

            reverse = pauli_pair_candidate(structure.conj().T, w0_b, w1_b, w0_a, w1_a, mb, ma, -alpha)
            max_reverse_error = max(max_reverse_error, relerr(reverse, candidate.conj().T))
    return max_dense_error, max_reverse_error


def test_multisublattice_phase_ownership() -> tuple[float, float]:
    """Physical-displacement convention rejects a basis/reference double phase."""
    # R, tau_a, and tau_b are fractional direct coordinates of a triclinic cell.
    cell = np.array([[3.3, 0.6, -0.3], [0.1, 2.8, 0.4], [0.2, 0.5, 4.0]])
    q_frac = np.array([0.23, -0.19, 0.31])
    q_cart = 2.0 * math.pi * np.linalg.solve(cell.T, q_frac)
    r_cell = np.array([1.0, -1.0, 0.0])
    tau_a = np.array([0.15, 0.20, 0.05])
    tau_b = np.array([0.62, 0.37, 0.41])
    basis_phase = float(q_cart @ (cell @ (tau_b - tau_a)))
    alpha = float(q_cart @ (cell @ (r_cell + tau_b - tau_a)))

    # Non-geometric reference azimuths belong only in m_a0/m_b0.
    phi_a, phi_b = 0.41, -0.83
    ma = reference_moment(0.91, phi_a)
    mb = reference_moment(1.27, phi_b)
    structure = np.array([[0.7 + 0.2j, -0.3 + 0.4j], [0.1 - 0.5j, 0.8 - 0.1j]])
    w0_a, w1_a = np.array([0.31, 0.52]), np.array([0.48, -0.21])
    w0_b, w1_b = np.array([0.88, 0.27]), np.array([-0.36, 0.44])

    correct = pauli_pair_candidate(structure, w0_a, w1_a, w0_b, w1_b, ma, mb, alpha)
    oracle = dense_pair_oracle(structure, w0_a, w1_a, w0_b, w1_b, ma, mb, alpha)
    correct_error = relerr(correct, oracle)

    # This deliberately mixes conventions: alpha already contains tau_b-tau_a,
    # while the extra basis and reference azimuth differences count them again.
    double_alpha = alpha + basis_phase + (phi_b - phi_a)
    double_counted = pauli_pair_candidate(structure, w0_a, w1_a, w0_b, w1_b, ma, mb, double_alpha)
    separation = relerr(double_counted, oracle)
    return correct_error, separation


def main() -> int:
    phase_error = test_cell_phase_equivalence()
    stoner_error = test_stoner_shifted_k_identity()
    dense_error, reverse_error = random_species_pair_cases()
    multisublattice_error, double_phase_separation = test_multisublattice_phase_ownership()

    print("GBT WP1 conventions: alpha=q_cart.(R+tau_b-tau_a), q_cart=2*pi*A^-T*q_frac; reference azimuths live only in endpoint moments.")
    print(f"cell Cartesian/fractional phase maximum relative error : {phase_error:.3e}")
    print(f"one-orbital Stoner k-/+q/2 maximum relative error      : {stoner_error:.3e}")
    print(f"dense unequal-species pair maximum relative error       : {dense_error:.3e}")
    print(f"reverse directed-bond maximum relative error            : {reverse_error:.3e}")
    print(f"multi-sublattice correct-convention relative error      : {multisublattice_error:.3e}")
    print(f"double-counted basis/reference phase separation         : {double_phase_separation:.3e}")

    errors = [phase_error, stoner_error, dense_error, reverse_error, multisublattice_error]
    if max(errors) >= TOL:
        print(f"G1 FAIL: maximum algebraic relative error {max(errors):.3e} is not < {TOL:.1e}")
        return 1
    if double_phase_separation <= 1.0e-3:
        print("G1 FAIL: the multi-sublattice fixture does not detect phase double counting")
        return 1
    print(f"G1 PASS: maximum algebraic relative error {max(errors):.3e} < {TOL:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
