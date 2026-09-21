#!/usr/bin/env python3
"""Independent checks for the DRESP-09ZS-R Fe production-span reconciliation."""

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

    assert values["primary_classification"] == "PRODUCTION_SPAN_STRICT_SUBSPACE_EXTENSION_REQUIRED"
    assert values["nesting"] == "STRICT_SUBSPACE"
    assert values["verdict"] == "PASS-A"
    assert values["recommended_DRESP_10A_migration"] == "EXTEND"
    assert values["ALSDA_Ward"] == "NOT RUN"
    assert values["BES_Halle"] == "OFF"
    assert values["raw_arbitrary_L_vertex"] == "CLOSED"
    assert values["live_production_four_branch_candidate_count"] == "232"
    assert values["second_order_four_branch_candidate_count"] == "232"
    assert values["second_order_six_branch_candidate_count"] == "348"
    assert values["live_production_retained_product_dimension"] == "232"
    assert values["rank_sensitivity_stable"] == "T"
    assert values["total_S4_prod_retained_rank_tau1_tau10_tau100"] == "232 232 232"
    assert values["total_S4_2nd_retained_rank_tau1_tau10_tau100"] == "232 232 232"
    assert values["total_S6_2nd_retained_rank_tau1_tau10_tau100"] == "348 348 348"
    assert values["total_S6_2nd_minus_S4_prod_rank_tau1"] == "116"
    assert values["production_oracle_rank_and_subspace_agreement"] == "T"

    assert max(vector(values, "old_span_outside_six_span maximum RMS", 2)) < 1.0e-8
    assert number(values, "S4_prod_vs_S4_2nd_principal_minimum") < 0.99
    assert number(values, "S4_prod_vs_S4_2nd_prod_outside") > 1.0e-2
    assert number(values, "S4_prod_vs_S4_2nd_second_order_outside") > 1.0e-2
    assert number(values, "S4_prod_vs_S6_2nd_principal_minimum") > 1.0 - 1.0e-8
    assert number(values, "S4_prod_vs_S6_2nd_prod_outside") < 1.0e-6
    assert number(values, "S4_prod_vs_S6_2nd_six_outside") > 1.0e-3
    assert max(vector(values, "physical_global_max_RMS", 2)) > 1.0e-8
    assert max(vector(values, "delta_O_global_max_RMS", 2)) > 1.0e-8
    assert number(values, "mixed_LM_full_SR_residual") > 1.0e-8
    assert number(values, "double_weighting_negative_residual") > 1.0e-3
    assert number(values, "L0_scalar_projection_regression") < 1.0e-12
    assert number(values, "full_candidate_residual_global") > 1.0e-8
    assert number(values, "S4_prod component augmentation delta-O") > 1.0e-8
    assert number(values, "S4_prod component SR total") > 1.0e-8

    duality = vector(values, "field_density_duality_single_mixed_circular_real", 4)
    assert max(duality) < 1.0e-12
    assert max(vector(values, "R20 maximum RMS", 2)) > 1.0e-8
    assert max(vector(values, "R02 maximum RMS", 2)) > 1.0e-8
    assert number(values, "R20 median") > 0.0
    assert number(values, "R02 median") > 0.0
    assert number(values, "mixed_LM_full_SR_residual") > 1.0e-8
    assert max(vector(values, "raw_pairing_abs", 4)) > 0.0
    assert max(vector(values, "projected_pairing_abs", 4)) > 0.0
    assert max(vector(values, "projection_loss", 4)) > 1.0e-14

    print("Dresp09ZS Fe artifact: six-branch span audit, physical activity, and duality pass")


if __name__ == "__main__":
    main()
