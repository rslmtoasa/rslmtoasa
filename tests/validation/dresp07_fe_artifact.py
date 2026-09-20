#!/usr/bin/env python3
"""Independent artifact checks for the bounded DRESP-07 Fe campaign."""

from __future__ import annotations

import math
import sys
from pathlib import Path


def read_artifact(path: Path) -> tuple[dict[str, float], list[list[float]], str]:
    values: dict[str, float] = {}
    ladder: list[list[float]] = []
    classification = ""
    in_ladder = False
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("# product_dimension = "):
            values["product_dimension"] = float(line.split("=", 1)[1].strip())
            continue
        if line.startswith("# eta_ladder"):
            in_ladder = True
            continue
        if line.startswith("# ward_ladder"):
            in_ladder = False
            continue
        if in_ladder and line and not line.startswith("#"):
            ladder.append([float(item) for item in line.split()])
            continue
        if line.startswith("classification = "):
            classification = line.split("=", 1)[1].strip()
            continue
        if " = " in line and not line.startswith("#"):
            key, value = line.split(" = ", 1)
            try:
                values[key.strip()] = float(value.strip())
            except ValueError:
                pass
    return values, ladder, classification


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/dresp07_fe_4k.dat")
    if not path.is_file():
        raise SystemExit(f"missing DRESP-07 artifact: {path}")
    values, ladder, classification = read_artifact(path)
    required = {
        "product_dimension",
        "max_eigenpair_residual",
        "max_eigenpair_frobenius_residual",
        "max_hermiticity_residual",
        "max_density_rotation_spectral_relative",
        "response_rotation_vs_spectral_relative",
        "response_exact_field_vs_spectral_relative",
        "accepted_core_to_total_ratio",
        "best_local_spherical_operator_relative",
        "best_local_spherical_operator_action_relative",
        "radial_best_BH_matrix_relative",
        "radial_best_BH_action_relative",
        "radial_map_svd_rank",
        "radial_map_svd_condition",
        "radial_map_svd_residual",
        "radial_map_svd_relative_residual",
        "radial_map_minimum_norm_profile",
        "radial_map_singular_value_1",
        "radial_map_singular_value_2",
        "radial_map_singular_value_3",
        "R_exactH_val",
        "R_exactH_total",
        "R_BKS_val",
        "R_Bxc_val",
        "Bxc_matrix_relative",
        "Bxc_action_relative",
        "BKS_matrix_relative",
        "BKS_action_relative",
        "best_BH_matrix_relative",
        "best_BH_action_relative",
        "gamma_BH_onsite_frobenius",
        "gamma_BH_frobenius",
        "gamma_BH_nonlocal_frobenius",
        "gamma_BH_orbital_offdiagonal_frobenius",
        "gamma_BH_within_l_diagonal_anisotropy_frobenius",
        "gamma_BH_cross_l_frobenius",
        "gamma_BH_intersite_frobenius",
        "BH_k_dependence_relative",
        "BH_global_local_projection_relative",
    }
    missing = sorted(required - values.keys())
    if missing:
        raise SystemExit(f"missing DRESP-07 keys: {', '.join(missing)}")
    if values["product_dimension"] != 232:
        raise SystemExit("DRESP-07 product dimension is not 232")
    for key in (
        "max_eigenpair_residual",
        "max_eigenpair_frobenius_residual",
        "max_hermiticity_residual",
        "max_density_rotation_spectral_relative",
        "response_rotation_vs_spectral_relative",
        "response_exact_field_vs_spectral_relative",
    ):
        if not math.isfinite(values[key]) or values[key] >= 1.0e-9:
            raise SystemExit(f"DRESP-07 exact gate failed: {key}={values[key]:.6e}")
    if values["accepted_core_to_total_ratio"] <= 0.0:
        raise SystemExit("DRESP-07 core bookkeeping channel is empty")
    for key in required - {"product_dimension"}:
        if not math.isfinite(values[key]):
            raise SystemExit(f"DRESP-07 non-finite diagnostic: {key}={values[key]}")
    if values["best_local_spherical_operator_relative"] > values["Bxc_matrix_relative"] + 1.0e-10:
        raise SystemExit("DRESP-07R best local operator is worse than mapped Bxc")
    if values["best_local_spherical_operator_relative"] > values["BKS_matrix_relative"] + 1.0e-10:
        raise SystemExit("DRESP-07R best local operator is worse than mapped BKS")
    if values["radial_map_svd_rank"] != 3 or values["radial_map_svd_relative_residual"] >= 1.0e-10:
        raise SystemExit("DRESP-07R radial inverse did not pass the full-rank SVD oracle")
    if abs(values["radial_best_BH_matrix_relative"] - values["best_local_spherical_operator_relative"]) >= 1.0e-10:
        raise SystemExit("DRESP-07R radial realization does not reproduce the best local operator")
    if values["gamma_BH_orbital_offdiagonal_frobenius"] <= 0.0 or values["gamma_BH_within_l_diagonal_anisotropy_frobenius"] <= 0.0:
        raise SystemExit("DRESP-07R Fe field decomposition lost the non-spherical onsite defects")
    if values["gamma_BH_cross_l_frobenius"] >= 1.0e-10 or values["gamma_BH_intersite_frobenius"] >= 1.0e-10:
        raise SystemExit("DRESP-07R Fe cross-l/intersite decomposition regressed")
    if values["BH_k_dependence_relative"] <= 1.0e-3:
        raise SystemExit("DRESP-07R Fe k-dependence diagnostic is unexpectedly empty")
    if len(ladder) != 4 or [row[0] for row in ladder] != [0.04, 0.02, 0.01, 0.005]:
        raise SystemExit("DRESP-07 eta ladder is not the required 0.04/0.02/0.01/0.005 sequence")
    for row in ladder:
        if len(row) != 5 or not all(math.isfinite(value) for value in row):
            raise SystemExit("DRESP-07 eta ladder contains a non-finite value")
    for column in range(1, 5):
        if any(ladder[i][column] >= ladder[i - 1][column] for i in range(1, len(ladder))):
            raise SystemExit("DRESP-07 eta ladder is not decreasing for every Ward source")
    if classification != "SECOND_ORDER_LMTO_MAPPING_REQUIRED":
        raise SystemExit(f"unexpected DRESP-07 classification: {classification}")
    print("Dresp07FeArtifact: PASS (independent exact-gate, bookkeeping, field, and eta-ladder checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
