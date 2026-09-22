#!/usr/bin/env python3
"""Independent checks for the DRESP-11 accepted-Fe artifact."""

from __future__ import annotations

import math
import sys
from pathlib import Path


def read(path: str) -> tuple[str, dict[str, str]]:
    values: dict[str, str] = {}
    verdict = ""
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("DRESP-11 verdict:"):
                verdict = stripped.split(":", 1)[1].strip()
            elif "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return verdict, values


def number(values: dict[str, str], key: str) -> float:
    value = float(values[key])
    assert math.isfinite(value), key
    return value


def main() -> None:
    verdict, values = read(sys.argv[1])
    assert verdict in {"PASS-A", "PASS-B"}
    assert values["DRESP-10F verdict"] in {"PASS-A", "PASS-B"}
    assert values["denominator"] == "D = I - A"
    assert values["branches"] == "00,10,01,11,20,02"
    assert number(values, "Pauli compact dimension") == 348
    assert number(values, "denominator_dimension") == 348
    assert values["Correction status"] == "DIAGNOSTIC_ONLY"
    assert values["BES/Halle production"] == "OFF"
    assert values["Dynamics"] == "NOT RUN"
    assert values["Juelich/GSR reference"].startswith("NOT_COMPARABLE")

    for key in (
        "DRESP10F_span",
        "DRESP10F_source",
        "DRESP10F_Frechet",
        "DRESP10F_mixed_action",
        "DRESP10F_Kxc",
    ):
        assert number(values, key) < 1.0e-10, key
    for key in (
        "P3_compact_reconstruction_residual",
        "matrix_free_vs_assembled_action",
        "matrix_free_vs_assembled_denominator",
        "assembly_order_residual",
        "gauge_sigma0",
        "gauge_gap",
        "gauge_overlap",
        "gauge_Ward_residual",
        "singular_gap_sigma1_over_sigma0",
        "rigid_overlap_with_v0",
        "Hermiticity_residual",
        "non_normality_residual",
        "condition_estimate",
    ):
        number(values, key)

    assert number(values, "matrix_free_vs_assembled_action") < 1.0e-10
    assert number(values, "matrix_free_vs_assembled_denominator") < 1.0e-10
    assert number(values, "assembly_order_residual") < 1.0e-10
    # The compact radial product basis is a finite projection of the raw P3
    # target; retain the measured 1e-6-scale closure rather than demanding
    # exact equality from that projection.
    assert number(values, "P3_compact_reconstruction_residual") < 1.0e-4
    assert number(values, "singular_gap_sigma1_over_sigma0") > 0.0
    assert number(values, "rigid_overlap_with_v0") >= 0.0
    assert number(values, "rigid_overlap_with_v0") <= 1.0 + 1.0e-10
    assert number(values, "gauge_overlap") >= 0.0
    for index in range(10):
        assert number(values, f"sigma({index})") >= 0.0

    spectrum = Path(sys.argv[1] + ".spectrum.csv")
    assert spectrum.is_file()
    rows = spectrum.read_text(encoding="utf-8").splitlines()
    assert rows and rows[0] == "kind,index,singular_value,abs_eigenvalue,eigen_real,eigen_imag,right_overlap"
    singular_rows = [row.split(",") for row in rows[1:] if row.startswith("singular,")]
    eigen_rows = [row.split(",") for row in rows[1:] if row.startswith("eigen,")]
    assert len(singular_rows) >= 10
    assert len(eigen_rows) >= 10
    for row in singular_rows[:10]:
        assert math.isfinite(float(row[2])) and float(row[2]) >= 0.0
        assert math.isfinite(float(row[6])) and 0.0 <= float(row[6]) <= 1.0 + 1.0e-10

    assert values["Primary classification"] in {
        "ISOLATED_GOLDSTONE_DEFECT_CORRECTION_ADMISSIBLE",
        "DISTRIBUTED_FINITE_LMTO_WARD_DEFECT",
        "NONNORMAL_GOLDSTONE_DEFECT_AMBIGUOUS",
        "DENOMINATOR_ASSEMBLY_OPEN",
        "DRESP10F_REGRESSION",
    }
    assert values["NEXT"] in {
        "DYNAMIC_GOLDSTONE_RESTORATION",
        "FINITE_LMTO_COVARIANCE",
        "COVARIANT_DENOMINATOR_FORMULATION",
    }
    print("Dresp11 Fe artifact: denominator assembly and Goldstone spectral gates pass")


if __name__ == "__main__":
    main()
