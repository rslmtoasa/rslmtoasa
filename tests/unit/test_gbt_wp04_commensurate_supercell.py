#!/usr/bin/env python3
"""WP04 exact matched-operator commensurate-supercell oracle.

The period-2 and period-3 cases use one immutable finite directed primitive
operator. No production GBT helper is imported: the primitive side
transcribes the endpoint-link equation and the explicit side is the lab-frame
supercell assembled from the same table.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
from gbt_wp04_matched_operator import (  # noqa: E402
    DirectedBond,
    MatchedPrimitiveOperator,
    block_diagonal,
    su2_frame,
)


TOL = 2.0e-12
NEGATIVE_TOL = 1.0e-3
ALAT = 3.7
THETA = 0.71
PHI = 0.37


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(
        np.linalg.norm(actual - expected, ord="fro")
        / max(1.0, float(np.linalg.norm(expected, ord="fro")))
    )


def make_fixture() -> MatchedPrimitiveOperator:
    """Return the canonical finite one-sublattice primitive bond table."""
    onsite = np.array(
        [[0.42 + 0.0j, 0.09 - 0.04j], [0.09 + 0.04j, -0.17 + 0.0j]],
        dtype=np.complex128,
    )
    plus_one = np.array(
        [[0.31 + 0.07j, -0.12 + 0.16j], [0.08 - 0.05j, 0.21 - 0.03j]],
        dtype=np.complex128,
    )
    plus_two = np.array(
        [[-0.13 + 0.02j, 0.06 + 0.11j], [-0.04 + 0.09j, 0.07 - 0.08j]],
        dtype=np.complex128,
    )

    def bond(delta: int, structure: np.ndarray) -> DirectedBond:
        return DirectedBond(
            source=0,
            target=0,
            cell_displacement=(delta, 0, 0),
            physical_displacement=(ALAT * delta, 0.0, 0.0),
            structure=structure,
        )

    bonds = (
        bond(0, onsite),
        bond(1, plus_one),
        bond(-1, plus_one.conj().T),
        bond(2, plus_two),
        bond(-2, plus_two.conj().T),
    )
    w0 = (np.array([1.07, 0.76], dtype=np.complex128),)
    w1 = (np.array([0.24, -0.19], dtype=np.complex128),)
    c0 = (np.array([-0.11, 0.18], dtype=np.complex128),)
    c1 = (np.array([0.29, -0.13], dtype=np.complex128),)
    frames = (su2_frame(THETA, PHI),)
    return MatchedPrimitiveOperator(ALAT, bonds, w0, w1, c0, c1, frames)


def run_case(operator: MatchedPrimitiveOperator, nperiod: int, qx: float, kx: float) -> dict[str, float]:
    q = np.array([qx, 0.0, 0.0], dtype=float)
    k_super = np.array([kx, 0.0, 0.0], dtype=float)
    supercell = operator.build_explicit_supercell(nperiod, k_super, q)
    folded = operator.folded_kpoints(nperiod, k_super)
    primitive = block_diagonal(operator.build_primitive_gbt(k, q) for k in folded)
    mapping = operator.folding_unitary(nperiod, k_super, q)

    unitary_error = relerr(mapping.conj().T @ mapping, np.eye(supercell.shape[0], dtype=np.complex128))
    matrix_error = relerr(supercell @ mapping, mapping @ primitive)
    spectrum_error = float(
        np.max(np.abs(np.linalg.eigvalsh(supercell) - np.linalg.eigvalsh(primitive)))
    )
    trace_error = float(abs(np.trace(supercell) - np.trace(primitive)))
    frobenius_error = float(abs(np.linalg.norm(supercell, ord="fro") - np.linalg.norm(primitive, ord="fro")))
    hermiticity_error = float(np.linalg.norm(supercell - supercell.conj().T, ord="fro"))

    # Deliberately omit the spinor half-period shift. This is a negative
    # control for the folding map, not an alternate acceptance criterion.
    wrong_fold = block_diagonal(
        operator.build_primitive_gbt(k_super + np.array([ell / nperiod, 0.0, 0.0]), q)
        for ell in range(nperiod)
    )
    wrong_fold_error = float(
        np.max(np.abs(np.linalg.eigvalsh(supercell) - np.linalg.eigvalsh(wrong_fold)))
    )
    return {
        "unitary": unitary_error,
        "matrix": matrix_error,
        "spectrum": spectrum_error,
        "trace": trace_error,
        "frobenius": frobenius_error,
        "hermiticity": hermiticity_error,
        "wrong_fold": wrong_fold_error,
    }


def main() -> int:
    operator = make_fixture()
    # The single table object is passed to both builders in every case.
    assert len(operator.canonical_bond_table()) == 5
    assert operator.basis_slice(0) == slice(0, operator.block_size)
    assert operator.basis_index(0, 1, operator.norb - 1) == operator.block_size - 1
    cases = ((2, 0.5, 0.019), (3, 1.0 / 3.0, -0.021))
    all_results: dict[int, dict[str, float]] = {}
    for nperiod, qx, kx in cases:
        result = run_case(operator, nperiod, qx, kx)
        all_results[nperiod] = result
        print(
            f"WP04 period-{nperiod} q=({qx:.16g},0,0) cone(theta,phi)=({THETA:.6f},{PHI:.6f}) "
            f"matrix={result['matrix']:.3e} spectrum={result['spectrum']:.3e} "
            f"trace={result['trace']:.3e} frobenius={result['frobenius']:.3e} "
            f"unitary={result['unitary']:.3e} hermiticity={result['hermiticity']:.3e} "
            f"wrong-fold={result['wrong_fold']:.3e}"
        )

    exact_errors = [
        value
        for result in all_results.values()
        for key, value in result.items()
        if key != "wrong_fold"
    ]
    if max(exact_errors) >= TOL:
        print(f"WP04 FAIL: exact residual {max(exact_errors):.3e} is not below {TOL:.1e}")
        return 1
    if min(result["wrong_fold"] for result in all_results.values()) <= NEGATIVE_TOL:
        print("WP04 FAIL: wrong folding map was not detected")
        return 1
    print(f"WP04 PASS: exact residual {max(exact_errors):.3e} < {TOL:.1e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
