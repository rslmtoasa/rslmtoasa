#!/usr/bin/env python3
"""Independent checks for the DRESP-09W Fe radial provenance artifact."""

from __future__ import annotations

import math
import sys


def read(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            if "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return values


def number(values: dict[str, str], key: str) -> float:
    value = float(values[key])
    if not math.isfinite(value):
        raise AssertionError(f"{key} is not finite")
    return value


def vector(values: dict[str, str], key: str) -> list[float]:
    result = [float(item) for item in values[key].split()]
    assert len(result) == 3, key
    assert all(math.isfinite(item) for item in result), key
    return result


def main() -> None:
    values = read(sys.argv[1])
    assert values["accepted_snapshot_complete"] == "CLOSED"
    assert values["accepted_snapshot_reproduces_SR1"] == "CLOSED"
    assert values["core_l_resolved_available"] == "T"
    assert values["primary_classification"] == "SR_RADIAL_RECONSTRUCTION_OPEN"
    assert values["verdict"] == "PASS-B"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"
    assert "CROSS-REPRESENTATION DIAGNOSTIC" in values["DRESP09V_0p305325"]
    assert values["origin_handling"].startswith("excluded from L2/max diagnostics")

    assert number(values, "P1_P2_relative") < 1.0e-8
    assert number(values, "SR2_SR3_relative") < 1.0e-8
    assert number(values, "SR1_SR4_relative") < 1.0e-5
    assert abs(number(values, "old_cross_representation_residual") - 0.305325) < 1.0e-12
    assert number(values, "new_SR_to_SR_valence_residual") > 1.0e-3
    assert number(values, "P1_P3_relative") > 1.0e-3
    assert number(values, "P2_P3_relative") > 1.0e-3
    assert number(values, "Pauli_response_vs_P3_relative") > 1.0e-3
    assert number(values, "SR_response_vs_SR2_relative") > 1.0e-3

    for key in (
        "P1_P2_per_l_s_p_d",
        "P1_P3_per_l_s_p_d",
        "P2_P3_per_l_s_p_d",
        "SR2_SR3_per_l_s_p_d",
        "Pauli_response_vs_P1_per_l_s_p_d",
        "Pauli_response_vs_P3_per_l_s_p_d",
        "SR_response_vs_SR2_per_l_s_p_d",
    ):
        assert all(item >= 0.0 for item in vector(values, key)), key

    for key in (
        "core_integrated_spin_moment",
        "core_weighted_L2_norm",
        "core_to_max_valence_radial_amplitude",
        "core_s_like_integrated_spin_moment",
        "core_s_like_weighted_L2_norm",
        "core_s_like_max_radial_amplitude",
    ):
        number(values, key)

    print("Dresp09W Fe artifact: target hierarchies and core provenance close; SR/Pauli response seam remains explicit")

if __name__ == "__main__":
    main()
