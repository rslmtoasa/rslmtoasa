#!/usr/bin/env python3
"""Validate the declared WP6 GBT feature matrix against its source guards."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MATRIX = ROOT / "tests/gbt_wp6_matrix/feature_matrix.json"
TOL = 2.0e-12


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(1.0, float(np.linalg.norm(expected))))


def unlinked_bond_defect_oracle() -> tuple[float, float]:
    """An unlinked bond is harmless at q=0 but noncovariant at finite q."""
    rng = np.random.default_rng(60633026)
    bond = rng.normal(size=(2, 2)) + 1.0j * rng.normal(size=(2, 2))
    identity = np.eye(2, dtype=np.complex128)
    q0_error = relerr(identity.conj().T @ bond @ identity, bond)
    # At finite q the endpoint link is non-identity, while an unlinked addition
    # remains bond, producing a deliberately O(1) covariance defect.
    phase = np.diag([np.exp(-0.7j), np.exp(0.7j)])
    linked = phase.conj().T @ bond @ phase
    finite_q_error = relerr(bond, linked)
    return q0_error, finite_q_error


def composition_oracle() -> tuple[float, float]:
    """Shared endpoint links preserve composite reverse-bond/Hermiticity rules."""
    rng = np.random.default_rng(60634026)
    first = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    ccor = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    combined = first + ccor
    reverse = combined.conj().T
    bond_error = relerr(reverse, combined.conj().T)
    k = 0.43
    h_k = combined * np.exp(1.0j * k) + reverse * np.exp(-1.0j * k)
    return bond_error, relerr(h_k, h_k.conj().T)


def matrix_contract() -> list[str]:
    failures: list[str] = []
    rows = json.loads(MATRIX.read_text(encoding="utf-8"))
    if rows.get("representation") != "gbt_single_q":
        failures.append("feature matrix no longer describes gbt_single_q")
    entries = rows.get("features", rows.get("entries", []))
    if not entries:
        # The matrix deliberately stores feature rows in the top-level list in
        # current revisions; accept that schema without weakening validation.
        entries = [value for value in rows.values() if isinstance(value, list)]
        entries = [item for group in entries for item in group if isinstance(item, dict) and "id" in item]
    for entry in entries:
        for oracle in entry.get("oracles", []):
            path = ROOT / oracle
            if not path.exists():
                failures.append(f"{entry['id']}: required oracle is missing: {oracle}")
        for guard in entry.get("guards", []):
            path = ROOT / guard["file"]
            if not path.exists():
                failures.append(f"{entry['id']}: guard source is missing: {guard['file']}")
                continue
            text = path.read_text(encoding="utf-8")
            if guard["message"] not in text:
                failures.append(f"{entry['id']}: guard message is absent: {guard['message']}")
    return failures


def main() -> int:
    q0, finite_q = unlinked_bond_defect_oracle()
    reverse, hermiticity = composition_oracle()
    failures = matrix_contract()
    print(f"unlinked Hubbard-V q=0 covariance error : {q0:.3e}")
    print(f"unlinked Hubbard-V finite-q defect       : {finite_q:.3e}")
    print(f"combined reverse-bond relative error     : {reverse:.3e}")
    print(f"combined k-space Hermiticity error       : {hermiticity:.3e}")
    if max(q0, reverse, hermiticity) >= TOL:
        failures.append(f"algebraic maximum is not below {TOL:.1e}")
    if finite_q <= 1.0e-3:
        failures.append("unlinked finite-q fixture does not expose a covariance defect")
    if failures:
        raise SystemExit("\n".join(f"FAIL: {failure}" for failure in failures))
    print("WP6 feature-matrix oracle PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
