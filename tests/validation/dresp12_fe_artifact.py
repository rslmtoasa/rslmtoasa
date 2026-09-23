#!/usr/bin/env python3
"""Independent checks for the DRESP-12 accepted-Fe covariance artifact."""

from __future__ import annotations

import math
import sys


def read(path: str) -> tuple[str, dict[str, str]]:
    values: dict[str, str] = {}
    verdict = ""
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("DRESP-12 verdict:"):
                verdict = stripped.split(":", 1)[1].strip()
            elif "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return verdict, values


def number(values: dict[str, str], key: str) -> float:
    value = float(values[key])
    assert math.isfinite(value), key
    return value


def vector(values: dict[str, str], key: str, size: int = 6) -> list[float]:
    fields = values[key].split()
    assert len(fields) == size
    result = [float(field) for field in fields]
    assert all(math.isfinite(field) for field in result)
    return result


def main() -> None:
    verdict, values = read(sys.argv[1])
    assert verdict == "PASS-A"
    assert values["Pauli compact dimension"] == "348"
    assert values["branches"] == "00,10,01,11,20,02"
    assert values["Goldstone correction"] == "OFF"
    assert values["BES/Halle production"] == "OFF"
    assert values["Dynamics"] == "NOT RUN"
    assert values["covariant_branch_decomposition"] == "CLOSED_BY_COMPLETE_ENDPOINT_TANGENT"
    assert values["arbitrary_L_within_ASA"] == "NONSPHERICAL_GROUND_STATE_RESPONSE_NOT_DEFINED"
    assert values["DRESP-11 regression"] == "LIVE_MATRIX_FREE_DENOMINATOR_ACTION"
    assert values["Model boundary"] == "NONSPHERICAL_RESPONSE_ON_SPHERICAL_ASA_GROUND_STATE"

    for key in (
        "field_identity_max_k_residual",
        "field_identity_weighted_RMS",
        "field_identity_max_matrix_element",
        "DRESP08_native_product_tangent_vs_connection_residual",
        "Frechet_oracle_residual",
        "independent_compact_connection_residual",
        "fixed_basis_goldstone_relative",
        "dm_cov_rho_vs_dm_B_plus_dm_conn",
        "dm_cov_frozen_endpoint_vs_direct_fixed_H",
        "dm_cov_complete_vs_frozen_plus_endpoint_H",
        "covariance_accounting_relative_to_fixed_defect",
        "master_identity_relative",
        "DRESP11_full_space_DmG_reconstruction_residual",
        "dmg_l0_raw_residual",
        "dmg_l0_compact_residual",
        "l0_mode_cosine",
        "l0_parallel_coefficient",
        "l0_orthogonal_fraction",
        "integrated_n_up_minus_down",
        "integrated_m_xc_common_radial_metric",
        "integrated_P3",
        "integrated_physical_SR_spin",
        "m_xc_vs_P3_weighted_profile_difference",
        "m_xc_vs_P3_weighted_L2_difference",
        "m_xc_vs_P3_integrated_difference",
        "m_xc_vs_P3_max_physical_difference",
        "m_xc_vs_P3_max_relative_difference",
        "integrated_core_m_xc",
        "m_xc_core_bookkeeping_residual",
        "constraining_field_ry",
        "constraining_field_max_abs_ry",
        "DmG_full",
        "DmG_norm_reconstruction_residual",
        "DRESP09Y_endpoint_measurement_norm",
        "DRESP12_endpoint_measurement_norm",
        "DRESP12_vs_DRESP09Y_endpoint_residual",
        "independent_cartesian_trace_oracle_vs_DRESP09Y_residual",
        "independent_cartesian_trace_oracle_vs_DRESP12_residual",
        "circular_plus_norm",
        "circular_minus_norm",
        "cartesian_x_reconstruction_residual",
        "cartesian_y_reconstruction_residual",
        "observable_deltaO_plus_minus_to_x_residual",
        "observable_deltaO_y_reconstruction_residual",
        "Embedded_DRESP09Y_Pauli_complete_vs_P3",
        "Current_DRESP12_Pauli_complete_vs_P3",
        "Pauli_measurement_seam_residual",
        "l0_residual",
        "nonspherical_norm",
        "nonspherical_residual",
        "nonspherical_L4_fraction",
        "r_fixed_full",
        "r_fixed_norm_reconstruction_residual",
        "master_full",
        "master_norm_reconstruction_residual",
        "circular_plus_accounting_residual",
        "circular_minus_accounting_residual",
    ):
        number(values, key)

    assert number(values, "field_identity_max_k_residual") < 1.0e-10
    assert number(values, "Frechet_oracle_residual") < 1.0e-10
    assert number(values, "independent_compact_connection_residual") < 1.0e-10
    assert number(values, "fixed_basis_goldstone_relative") < 0.25
    assert number(values, "dm_cov_rho_vs_dm_B_plus_dm_conn") < 1.0e-10
    assert number(values, "dm_cov_frozen_endpoint_vs_direct_fixed_H") < 1.0e-10
    assert number(values, "dm_cov_complete_vs_frozen_plus_endpoint_H") < 1.0e-10
    assert number(values, "DRESP09Y_endpoint_measurement_norm") > 0.0
    assert number(values, "DRESP12_endpoint_measurement_norm") > 0.0
    assert number(values, "DRESP12_vs_DRESP09Y_endpoint_residual") < 1.0e-10
    assert number(values, "independent_cartesian_trace_oracle_vs_DRESP09Y_residual") < 1.0e-10
    assert number(values, "independent_cartesian_trace_oracle_vs_DRESP12_residual") < 1.0e-10
    assert number(values, "circular_plus_norm") > 0.0
    assert number(values, "circular_minus_norm") > 0.0
    assert number(values, "cartesian_x_reconstruction_residual") < 1.0e-10
    assert number(values, "cartesian_y_reconstruction_residual") < 1.0e-10
    assert number(values, "observable_deltaO_plus_minus_to_x_residual") < 1.0e-10
    assert number(values, "observable_deltaO_y_reconstruction_residual") < 1.0e-10
    assert number(values, "Embedded_DRESP09Y_Pauli_complete_vs_P3") < 3.0e-8
    assert number(values, "Current_DRESP12_Pauli_complete_vs_P3") < 3.0e-8
    assert number(values, "Pauli_measurement_seam_residual") < 3.0e-8
    assert values["Measurement seam classification"] == "MEASUREMENT_CONVENTIONS_IDENTICAL"
    assert number(values, "r_fixed_norm_reconstruction_residual") < 1.0e-10
    assert number(values, "master_norm_reconstruction_residual") < 1.0e-10
    assert number(values, "DmG_norm_reconstruction_residual") < 1.0e-10
    for key in (
        "endpoint_Hermitian_branch_swap_residual_00_10_01_11_20_02",
        "radial_spin_direction_branch_swap_residual_00_10_01_11_20_02",
        "DRESP09Y_channel1_branch_norm_00_10_01_11_20_02",
        "DRESP09Y_channel2_branch_norm_00_10_01_11_20_02",
        "DRESP12_branch_norm_00_10_01_11_20_02",
        "correct_circular_reconstruction_branch_norm_00_10_01_11_20_02",
        "branch_circular_reconstruction_residual_00_10_01_11_20_02",
        "branch_DRESP12_vs_correct_residual_00_10_01_11_20_02",
    ):
        vector(values, key)
    assert max(vector(values, "endpoint_Hermitian_branch_swap_residual_00_10_01_11_20_02")) < 1.0e-10
    assert max(vector(values, "radial_spin_direction_branch_swap_residual_00_10_01_11_20_02")) < 1.0e-10
    assert max(vector(values, "branch_circular_reconstruction_residual_00_10_01_11_20_02")) < 1.0e-10
    branch_keys = (
        "endpoint_branch_complete_frozen_endpoint_H_residual_00_10_01_11_20_02",
        "endpoint_H_subtraction_vs_delta_rho_zero_00_10_01_11_20_02",
        "endpoint_complete_vs_authoritative_covariant_00_10_01_11_20_02",
    )
    for key in branch_keys:
        fields = values[key].split()
        assert len(fields) == 6 and all(math.isfinite(float(field)) for field in fields)
        assert max(float(field) for field in fields) < 1.0e-10
    for key in (
        "r_fixed_by_L_0_1_2_3_4",
        "master_by_L_0_1_2_3_4",
        "DmG_by_L_0_1_2_3_4",
    ):
        vector(values, key, size=5)
    assert values["Primary classification"] == "ASA_L0_RIGID_RESPONSE_CLOSED"
    assert abs(number(values, "l0_residual")) < 1.0e-8
    assert number(values, "dmg_l0_raw_residual") < 1.0e-6
    assert number(values, "dmg_l0_compact_residual") >= 0.0
    assert values["dmg_l0_compact_residual_role"] == "PROJECTION_COMPRESSION_DIAGNOSTIC"
    assert values["Ward-focused campaign"] == "CLOSED"
    assert values["NEXT"] == "NATIVE_ROTATION_DYNAMICS"
    assert values["m_xc definition"] == "n_up-n_down from the live VXC0SP spherical density"
    assert values["m_xc_vs_P3_sign_convention"] == "PASS_UP_MINUS_DOWN"
    assert values["core_inclusion"] == "m_xc includes captured frozen core plus valence; P3 is accepted reciprocal valence large-component density"
    assert values["bxc_pauli_definition"] == "0.5*(vxc_up-vxc_down); XC field only; constraining field excluded"
    assert abs(number(values, "constraining_field_max_abs_ry")) < 1.0e-12
    assert abs(number(values, "m_xc_core_bookkeeping_residual")) < 1.0e-12

    ktable = sys.argv[1] + ".kpoints.csv"
    with open(ktable, encoding="utf-8") as stream:
        rows = stream.read().splitlines()
    assert rows[0] == "k_index,R_H,band_action_max,distance_to_EF"
    assert len(rows) == 65
    for row in rows[1:]:
        fields = row.split(",")
        assert len(fields) == 4
        for field in fields[1:]:
            assert math.isfinite(float(field))

    print("Dresp12 Fe artifact: covariance accounting and compact reconstruction gates pass")


if __name__ == "__main__":
    main()
