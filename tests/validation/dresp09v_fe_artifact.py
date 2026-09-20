#!/usr/bin/env python3
"""Independent checks for the DRESP-09V Fe density-moment artifact."""

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
    assert values["accepted_kpoints"] == "64"
    assert values["classification"] == "DENSITY_MOMENTS_CLOSED_RADIAL_TARGET_OPEN"
    assert values["verdict"] == "PASS-B"
    assert values["explicit_coefficient_X_required"] == "NO"
    assert values["DRESP09T_density_correction"] == "RETIRED_NOT_USED"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"

    for key in (
        "commutator_H_rho_residual",
        "M0_product_rule_vs_Frechet",
        "M1_product_rule_vs_Frechet",
        "M2_product_rule_vs_Frechet",
        "M0_equilibrium_product_residual",
        "rigid_delta_M0_commutator",
        "rigid_delta_M1_commutator",
        "rigid_delta_M2_commutator",
        "endpoint_D1_identity",
        "endpoint_D2_identity",
        "endpoint_D3_identity",
        "endpoint_D4_identity",
        "endpoint_D10_minus_D01",
    ):
        assert number(values, key) < 1.0e-9, key

    finite = [
        number(values, "finite_angle_theta_1e-2"),
        number(values, "finite_angle_theta_5e-3"),
        number(values, "finite_angle_theta_2p5e-3"),
    ]
    assert finite[0] > finite[1] > finite[2] > 0.0
    assert 3.0 < finite[0] / finite[1] < 5.0
    assert 3.0 < finite[1] / finite[2] < 5.0

    assert number(values, "production_spin_density_equilibrium_mapping") < 2.0e-12
    for key in ("production_spin_density_M0_fd", "production_spin_density_M1_fd", "production_spin_density_M2_fd"):
        assert number(values, key) < 2.0e-6, key

    assert abs(number(values, "old_DRESP09S_fixed_H_residual") - 4.08210e-2) < 1.0e-6
    assert number(values, "endpoint_H_contribution_norm") > 0.0
    assert number(values, "complete_DRESP09V_residual") > 1.0e-3
    for key in (
        "integrated_valence_moment_residual",
        "maximum_radial_relative_residual",
        "maximum_radial_absolute_residual",
        "per_l_s",
        "per_l_p",
        "per_l_d",
    ):
        assert number(values, key) >= 0.0, key

    for key in (
        "field_density_fixed_H_pairing_residual",
        "field_density_endpoint_H_pairing_residual",
        "field_density_total_pairing_residual",
    ):
        assert number(values, key) < 1.0e-9, key

    print("Dresp09V Fe artifact: coefficient moment and production spin-density tangent close; radial target remains open")


if __name__ == "__main__":
    main()
