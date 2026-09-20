#!/usr/bin/env python3
"""Independent artifact checks for the bounded DRESP-08 Fe campaign."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


def read_artifact(path: Path) -> tuple[dict[str, float], str, str]:
    values: dict[str, float] = {}
    classification = ""
    verdict = ""
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("classification = "):
            classification = line.split("=", 1)[1].strip()
            continue
        if line.startswith("verdict = "):
            verdict = line.split("=", 1)[1].strip()
            continue
        if " = " in line and not line.startswith("#"):
            key, value = line.split(" = ", 1)
            try:
                values[key.strip()] = float(value.strip())
            except ValueError:
                pass
    return values, classification, verdict


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/dresp08_fe_4k.dat")
    if not path.is_file():
        raise SystemExit(f"missing DRESP-08 artifact: {path}")
    values, classification, verdict = read_artifact(path)
    required = {
        "product_dimension",
        "H2_max_element_residual",
        "H2_max_frobenius_residual",
        "H2_max_relative_frobenius_residual",
        "H2_max_hermiticity_residual",
        "H2_fourier_eeo_vs_h_o_residual",
        "potential_parameter_average_identity_residual",
        "obarm_parameter_reconstruction_residual",
        "enim_parameter_reconstruction_residual",
        "spin_resolved_BH_max_residual",
        "spin_resolved_BH_max_relative_residual",
        "first_order_tangent_vs_commutator",
        "overlap_tangent_vs_commutator",
        "Enu_tangent_vs_commutator",
        "second_order_product_tangent_vs_commutator",
        "second_order_product_tangent_max_element",
        "B_E_weighted_norm",
        "B_h_weighted_norm",
        "B_hoh_weighted_norm",
        "B_H_weighted_norm",
        "Rk_B_h",
        "Rk_B_hoh",
        "Rk_B_H",
        "Rk_delta_h",
        "Rk_delta_H_complete",
        "Gamma_BH_local_spherical_projection_relative",
        "Gamma_BH_same_l_orbital_offdiagonal",
        "Gamma_BH_within_l_m_anisotropy",
        "Gamma_BH_cross_l",
        "Gamma_BH_intersite",
        "old_radial_Bxc_vs_native_matrix",
        "old_radial_Bxc_vs_native_action",
        "old_radial_BKS_vs_native_matrix",
        "old_radial_BKS_vs_native_action",
        "native_product_response_vs_exact_static",
        "native_product_response_vs_exact_rotation",
    }
    missing = sorted(required - values.keys())
    if missing:
        raise SystemExit(f"missing DRESP-08 keys: {', '.join(missing)}")
    if values["product_dimension"] != 232:
        raise SystemExit("DRESP-08 product dimension is not 232")
    if classification != "NATIVE_SECOND_ORDER_MAPPING_CLOSED" or verdict != "PASS":
        raise SystemExit(f"DRESP-08 did not close: classification={classification!r}, verdict={verdict!r}")
    text = path.read_text(encoding="utf-8")
    for flag in ("BES_Halle = OFF", "Kxc_tuning = OFF", "Goldstone_correction = OFF"):
        if f"# {flag}" not in text:
            raise SystemExit(f"DRESP-08 forbidden correction flag is not OFF: {flag}")

    exact_keys = [
        key
        for key in required
        if key.endswith("residual") or "tangent" in key or key.startswith("native_product_response")
    ]
    for key in exact_keys:
        if not math.isfinite(values[key]) or values[key] >= 2.0e-10:
            raise SystemExit(f"DRESP-08 exact gate failed: {key}={values[key]:.6e}")
    for key, value in values.items():
        if not math.isfinite(value):
            raise SystemExit(f"DRESP-08 non-finite diagnostic: {key}={value}")

    if values["Rk_B_h"] <= 1.0e-3 or values["Rk_B_hoh"] <= 1.0e-3 or values["Rk_B_H"] <= 1.0e-3:
        raise SystemExit("DRESP-08 Fe k-dependence diagnostic is unexpectedly empty")
    if values["Gamma_BH_cross_l"] >= 1.0e-10 or values["Gamma_BH_intersite"] >= 1.0e-10:
        raise SystemExit("DRESP-08 angular audit lost the spherical/no-SOC cross-l or intersite closure")
    if values["B_h_weighted_norm"] <= values["B_hoh_weighted_norm"]:
        raise SystemExit("DRESP-08 Fe decomposition no longer distinguishes B_h from B_hoh")

    eta = []
    for key, value in values.items():
        match = re.fullmatch(r"eta_(\d+)", key)
        if match:
            eta.append((int(match.group(1)), value))
    eta.sort()
    if [value for _, value in eta] != [0.04, 0.02, 0.01, 0.005]:
        raise SystemExit("DRESP-08 eta ladder is not the required 0.04/0.02/0.01/0.005 sequence")
    for index in range(1, 5):
        key = f"native_vs_exact_eta_relative_{index}"
        if key not in values or values[key] >= 2.0e-10:
            raise SystemExit(f"DRESP-08 eta response gate failed: {key}={values.get(key)}")

    k_residuals = {
        int(match.group(1)): value
        for key, value in values.items()
        if (match := re.fullmatch(r"k_(\d+)_native_vs_commutator_relative", key))
    }
    nk = 64  # bounded Fe fixture is explicitly 4x4x4
    if sorted(k_residuals) != list(range(1, nk + 1)):
        raise SystemExit(f"DRESP-08 per-k tangent audit is incomplete: found {len(k_residuals)} k points")
    if max(k_residuals.values()) >= 2.0e-10:
        raise SystemExit("DRESP-08 per-k native tangent residual exceeded the exact gate")
    k_spherical = {
        int(match.group(1)): value
        for key, value in values.items()
        if (match := re.fullmatch(r"k_(\d+)_BH_local_spherical_projection_relative", key))
    }
    if sorted(k_spherical) != list(range(1, nk + 1)) or not all(math.isfinite(value) for value in k_spherical.values()):
        raise SystemExit("DRESP-08 per-k spherical B_H audit is incomplete")
    print("Dresp08FeArtifact: PASS (independent reconstruction, provenance, angular, tangent, and response checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
