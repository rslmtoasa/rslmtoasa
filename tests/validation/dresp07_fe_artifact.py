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
        "best_BH_relative_residual",
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
        "gamma_BH_nonlocal_frobenius",
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
    if not math.isfinite(values["best_BH_relative_residual"]) or values["best_BH_relative_residual"] <= 1.0e-3:
        raise SystemExit("DRESP-07 best-field representation defect was not retained")
    if len(ladder) != 4 or [row[0] for row in ladder] != [0.04, 0.02, 0.01, 0.005]:
        raise SystemExit("DRESP-07 eta ladder is not the required 0.04/0.02/0.01/0.005 sequence")
    for row in ladder:
        if len(row) != 5 or not all(math.isfinite(value) for value in row):
            raise SystemExit("DRESP-07 eta ladder contains a non-finite value")
    if ladder[-1][1] >= ladder[0][1]:
        raise SystemExit("DRESP-07 exact-H retarded response did not approach the static oracle")
    if classification != "FIELD_REPRESENTATION_FAILURE":
        raise SystemExit(f"unexpected DRESP-07 classification: {classification}")
    print("Dresp07FeArtifact: PASS (independent exact-gate, bookkeeping, field, and eta-ladder checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
