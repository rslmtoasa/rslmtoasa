#!/usr/bin/env python3
"""Independent checks for the DRESP-09U production representation artifact."""

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


def main() -> None:
    values = read(sys.argv[1])
    assert values["classification"] == "FIELD_CLOSED_DENSITY_REPRESENTATION_OPEN"
    assert values["verdict"] == "PASS-B"
    assert values["accepted_kpoints"] == "64"
    assert values["predls_analytic_tangent"] == "CLOSED"
    assert values["predls_finite_difference_oracle"] == "CLOSED_IN_UNIT_FIXTURE"
    assert values["explicit_X_available"] == "NO"
    for key in (
        "ee_onsite_cex_tangent_residual",
        "ee_left_endpoint_tangent_residual",
        "ee_right_endpoint_tangent_residual",
        "ee_combined_tangent_residual",
        "obarm_tangent_residual",
        "enim_tangent_residual",
        "first_order_h_weighted_rms",
        "second_order_H_weighted_rms",
    ):
        assert number(values, key) < 1.0e-9, key
    theta = [
        number(values, "finite_rotation_representation_theta_1e-2"),
        number(values, "finite_rotation_representation_theta_5e-3"),
        number(values, "finite_rotation_representation_theta_2p5e-3"),
    ]
    assert theta[0] > theta[1] > theta[2] > 0.0
    assert 3.0 < theta[0] / theta[1] < 5.0
    assert 3.0 < theta[1] / theta[2] < 5.0
    print("Dresp09U Fe artifact: production representation tangent closes; density coefficient map remains explicitly open")


if __name__ == "__main__":
    main()
