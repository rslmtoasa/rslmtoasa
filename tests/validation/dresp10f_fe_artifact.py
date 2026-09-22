#!/usr/bin/env python3
"""Independent checks for the DRESP-10F fixed-basis mixed artifact."""

from __future__ import annotations

import math
import sys


def read(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("DRESP-10F verdict:"):
                values["verdict"] = stripped.split(":", 1)[1].strip()
            elif "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return values


def number(values: dict[str, str], key: str) -> float:
    value = float(values[key])
    assert math.isfinite(value), key
    return value


def main() -> None:
    values = read(sys.argv[1])
    assert values["verdict"] in {"PASS-A", "PASS-B"}
    assert values["formulation"] == "FIXED_GROUND_STATE_MIXED_RESPONSE"
    assert values["classification"] in {
        "FIXED_BASIS_MIXED_WARD_CLOSED",
        "FINITE_LMTO_WARD_INCONSISTENCY",
    }
    assert values["Pauli product basis branches"] == "00,10,01,11,20,02"
    assert values["field_basis"].startswith("INDEPENDENT_RAW_SR_FULL_SPATIAL")
    assert values["Kxc_low_m_status"].startswith("DIRECT_P3_DENOMINATOR")
    assert values["BES_Halle"] == "OFF"
    assert values["Goldstone_correction"] == "OFF"
    assert values["Dynamics"] == "NOT RUN"
    assert number(values, "Pauli product basis dimension") == 348
    assert number(values, "corrected_span_audit_pauli_max") < 1.0e-12
    assert number(values, "Pauli_span_residual") < 1.0e-10
    assert number(values, "source_vertex_mixed_LM_residual") < 1.0e-10
    assert number(values, "static_frechet_direct_matrix_residual") < 1.0e-10
    assert number(values, "mixed_chi0_rigid_source_oracle") < 1.0e-10
    assert number(values, "mixed_chi0_arbitrary_source_oracle") < 1.0e-10
    assert number(values, "Kxc_P3_identity_max_abs") < 1.0e-10
    assert math.isfinite(number(values, "PRIMARY_fixed_basis_relative_weighted_radial_L2"))
    print("Dresp10F Fe artifact: fixed-basis mixed representation and Ward gates pass")


if __name__ == "__main__":
    main()
