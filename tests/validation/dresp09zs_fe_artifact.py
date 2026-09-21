#!/usr/bin/env python3
"""Independent checks for the DRESP-09ZS Fe compact-span audit artifact."""

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

    assert values["primary_classification"] == (
        "RAW_ARBITRARY_L_VERTEX_CLOSED_PRODUCT_BASIS_EXTENSION_REQUIRED"
    )
    assert values["verdict"] == "PASS-B"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"
    assert values["raw_arbitrary_L_vertex"] == "CLOSED"
    assert values["rank_sensitivity_stable"] == "T"
    assert values["total_old_retained_rank_tau1_tau10_tau100"] == "232 232 232"
    assert values["total_six_retained_rank_tau1_tau10_tau100"] == "348 348 348"
    assert values["total_rank_increase_tau1"] == "116"

    assert max(vector(values, "old_span_outside_six_span maximum RMS", 2)) < 1.0e-8
    assert max(vector(values, "physical_global_max_RMS", 2)) > 1.0e-8
    assert max(vector(values, "delta_O_global_max_RMS", 2)) > 1.0e-8
    assert number(values, "mixed_LM_full_SR_residual") > 1.0e-8
    assert number(values, "double_weighting_negative_residual") > 1.0e-3
    assert number(values, "L0_scalar_projection_regression") < 1.0e-12
    assert number(values, "full_candidate_residual_global") > 1.0e-8

    duality = vector(values, "field_density_duality_single_mixed_circular_real", 4)
    assert max(duality) < 1.0e-12
    assert max(vector(values, "R20 maximum RMS", 2)) > 1.0e-8
    assert max(vector(values, "R02 maximum RMS", 2)) > 1.0e-8
    assert number(values, "R20 median") > 0.0
    assert number(values, "R02 median") > 0.0

    print("Dresp09ZS Fe artifact: six-branch span audit, physical activity, and duality pass")


if __name__ == "__main__":
    main()
