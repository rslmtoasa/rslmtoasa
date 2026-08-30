#!/usr/bin/env python3
"""Dependency-free unit checks for the WP11 LKAG/J(q) analysis path."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "analyze_gbt_wp11.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("gbt_wp11_lkag_jq", TOOL)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load {TOOL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def record(index: int, species_i: int, species_j: int, bond: tuple[float, float, float], jij_mry: float, multiplicity: float = 1.0):
    return {
        "index": index,
        "species_i": species_i,
        "species_j": species_j,
        "R": bond,
        "J_Ry": jij_mry * 1.0e-3,
        "J_mRy": jij_mry,
        "distance": math.sqrt(sum(value * value for value in bond)),
        "multiplicity": multiplicity,
    }


def bcc_two_shells():
    return [
        record(1, 1, 1, (-0.5, -0.5, -0.5), 1.0, 8.0),
        record(2, 1, 1, (0.0, 0.0, -1.0), 0.5, 6.0),
    ]


def test_scalar_jq_uses_total_shell_multiplicity_and_wp11_prefactor() -> None:
    tool = load_tool()
    records = bcc_two_shells()
    jq0 = tool.construct_jq(records, (0.0, 0.0, 0.0)).real
    jq_half = tool.construct_jq(records, (0.0, 0.0, 0.5)).real
    assert abs(jq0 - 0.011) < 1.0e-14
    # q=(0,0,1/2): the first shell has cos(pi/2)=0 and the second has cos(pi)=-1.
    assert abs(jq_half + 0.003) < 1.0e-14
    omega = 2.0 * (jq0 - jq_half) / 2.0
    assert abs(omega - 0.014) < 1.0e-14


def test_directed_fourier_phase_and_conjugate_reversal() -> None:
    tool = load_tool()
    records = [record(1, 1, 1, (0.25, 0.0, 0.0), 2.0)]
    plus = tool.construct_jq(records, (0.2, 0.0, 0.0), centrosymmetric=False)
    minus = tool.construct_jq(records, (-0.2, 0.0, 0.0), centrosymmetric=False)
    assert abs(plus - minus.conjugate()) < 1.0e-14
    assert abs(plus.imag - 2.0e-3 * math.sin(2.0 * math.pi * 0.05)) < 1.0e-14


def test_multisublattice_matrix_has_goldstone_and_optic_branches() -> None:
    tool = load_tool()
    records = [
        record(1, 1, 1, (1.0, 0.0, 0.0), 1.0, 2.0),
        record(2, 2, 2, (1.0, 0.0, 0.0), 1.0, 2.0),
        record(3, 1, 2, (0.5, 0.0, 0.0), 0.5),
        record(4, 2, 1, (-0.5, 0.0, 0.0), 0.5),
    ]
    matrix = tool.construct_exchange_matrix(records, (0.0, 0.0, 0.0), n_sublattices=2)
    assert abs(matrix[0][1].real - 0.0005) < 1.0e-14
    frequencies = tool.construct_acoustic_frequencies(records, (0.0, 0.0, 0.0), [1.0, 1.0])
    assert abs(frequencies[0]) < 1.0e-12
    assert abs(frequencies[1] - 0.002) < 1.0e-12


def test_range_convergence_is_not_claimed_for_one_exchange_file() -> None:
    tool = load_tool()
    rows = []
    sin2 = math.sin(math.radians(5.0)) ** 2
    for q in (0.0, 0.02, -0.02, 0.04, -0.04):
        omega = 2.0 * (tool.construct_jq(bcc_two_shells(), (0.0, 0.0, 0.0)).real -
                       tool.construct_jq(bcc_two_shells(), (0.0, 0.0, q)).real) / 2.0
        rows.append({
            "source": "synthetic", "line": len(rows) + 1, "mode": "mft",
            "q": (0.0, 0.0, q), "theta_deg": 5.0, "mesh": (16, 16, 16),
            "mtot": 2.0, "sin2theta": sin2,
            "delta_gauge": omega * 2.0 * sin2 / 4.0,
            "omega": omega, "gauge_available": True,
        })
    result = tool.analyze_wp11([{"name": "two_shell", "records": bcc_two_shells()}], rows)
    assert result["acceptance"]["exchange_range"]["status"] == "not_assessed"
    assert result["acceptance"]["q0_goldstone"]["status"] == "within_tolerance"


if __name__ == "__main__":
    test_scalar_jq_uses_total_shell_multiplicity_and_wp11_prefactor()
    test_directed_fourier_phase_and_conjugate_reversal()
    test_multisublattice_matrix_has_goldstone_and_optic_branches()
    test_range_convergence_is_not_claimed_for_one_exchange_file()
    print("GBT WP11 LKAG/J(q) unit checks PASS")
