#!/usr/bin/env python3
"""Independent checks for the DRESP-09Y accepted-Fe artifact."""

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

    assert values["primary_classification"] == "AUGMENTATION_FRAME_TANGENT_CLOSED"
    assert values["verdict"] == "PASS-A"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"
    assert number(values, "augmentation_analytic_deltaA_fd_max") < 4.0e-8
    assert number(values, "augmentation_upper_fd_max") < 4.0e-8
    assert number(values, "augmentation_lower_small_fd_max") < 4.0e-8
    assert number(values, "augmentation_lower_angular_fd_max") < 4.0e-8
    assert number(values, "augmentation_total_fd_max") < 4.0e-8
    assert number(values, "finite_angle_branch_upper_max") < 4.0e-8
    assert number(values, "finite_angle_branch_lower_small_max") < 4.0e-8
    assert number(values, "finite_angle_branch_lower_angular_max") < 4.0e-8
    assert number(values, "finite_angle_branch_total_max") < 4.0e-8
    assert number(values, "circular_delta_O_hermitian_max") < 4.0e-8
    assert number(values, "complete_response_vs_SRspin2_relative") < 3.0e-8
    assert abs(number(values, "complete_minus_SRspin2_integrated")) < 3.0e-12
    assert number(values, "Pauli_complete_vs_P3_relative") < 3.0e-8
    ratios = vector(values, "finite_angle_ratios", 2)
    assert all(ratio > 3.0 for ratio in ratios)
    assert max(abs(item) for item in vector(values, "delta_O_branch_00_10_01_11_20_02", 6)) < 4.0e-8
    for key in (
        "fixed_observable_vs_SRspin2_relative",
        "augmentation_vs_SRspin2_relative",
        "Pauli_fixed_vs_P1_relative",
        "Pauli_fixed_vs_P3_relative",
        "Pauli_complete_vs_P1_relative",
    ):
        number(values, key)
    for key in ("fixed_observable_per_l_s_p_d", "augmentation_per_l_s_p_d", "complete_per_l_s_p_d"):
        vector(values, key, 3)

    print("Dresp09Y Fe artifact: augmentation-frame tangent and complete SR closure pass")


if __name__ == "__main__":
    main()
