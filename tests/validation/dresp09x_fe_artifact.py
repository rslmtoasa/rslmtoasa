#!/usr/bin/env python3
"""Independent checks for the DRESP-09X Fe physical-spin artifact."""

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
    assert math.isfinite(value), key
    return value


def vector(values: dict[str, str], key: str, length: int) -> list[float]:
    result = [float(item) for item in values[key].split()]
    assert len(result) == length, key
    assert all(math.isfinite(item) for item in result), key
    return result


def main() -> None:
    values = read(sys.argv[1])

    assert values["primary_classification"] == "SECOND_ORDER_ENDPOINT_CLOSED_AUGMENTATION_ROTATION_REQUIRED"
    assert values["verdict"] == "PASS-B"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"
    assert values["observable_rotation_contribution"].startswith("OPEN")
    assert values["complete_rigid_response_vs_SRspin2"].startswith("OPEN")
    assert "SUPERSEDED TARGET MISMATCH" in values["historical_DRESP09W_SR_response_vs_SR2"]

    assert max(abs(item) for item in vector(values, "KH_angular_oracle_x_y_z", 3)) < 1.0e-9
    assert number(values, "SRspin2_SRspin3_relative") < 1.0e-10
    assert number(values, "Pauli_six_vs_P3_relative") < 1.0e-10
    assert number(values, "SR_probability_six_vs_SR2_relative") < 1.0e-10
    assert number(values, "SR_physical_spin_six_vs_SRspin2_relative") < 1.0e-10
    assert abs(number(values, "SR2_minus_SRspin2_integrated_difference")) > 1.0e-6
    assert abs(number(values, "SRspin2_minus_SRspin3_integrated_difference")) < 1.0e-10

    endpoint = vector(values, "endpoint_tangent_00_10_01_11_20_02", 6)
    assert max(endpoint) < 3.0e-10
    ratios = vector(values, "finite_angle_theta_ratios", 2)
    assert all(3.0 < ratio < 5.0 for ratio in ratios)

    for key in (
        "old_four_branch_vs_SR2_relative",
        "old_four_branch_vs_SRspin2_relative",
        "six_branch_fixed_observable_vs_SRspin2_relative",
        "Pauli_response_vs_P1_relative",
        "Pauli_response_vs_P3_relative",
        "Pauli_P1_vs_P3_relative",
    ):
        number(values, key)

    for key in ("SR2_minus_SRspin2_per_l_s_p_d", "six_spin_per_l_s_p_d", "fixed_observable_per_l_s_p_d"):
        assert all(item >= 0.0 for item in vector(values, key, 3)), key

    print("Dresp09X Fe artifact: physical-spin target and six-branch closure pass; augmentation tangent remains PASS-B")


if __name__ == "__main__":
    main()
