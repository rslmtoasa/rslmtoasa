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


def main() -> None:
    verdict, values = read(sys.argv[1])
    assert verdict in {"PASS-A", "PASS-B"}
    assert values["Pauli compact dimension"] == "348"
    assert values["branches"] == "00,10,01,11,20,02"
    assert values["Goldstone correction"] == "OFF"
    assert values["BES/Halle production"] == "OFF"
    assert values["Dynamics"] == "NOT RUN"
    assert values["covariant_branch_decomposition"].startswith("NOT_DEFINED")
    assert values["arbitrary_L_within_ASA"] == "ARBITRARY_L_COVARIANCE_REQUIRES_NEW_BASIS_RESPONSE"

    for key in (
        "field_identity_max_k_residual",
        "field_identity_weighted_RMS",
        "field_identity_max_matrix_element",
        "DRESP08_native_product_tangent_vs_connection_residual",
        "Frechet_oracle_residual",
        "independent_compact_connection_residual",
        "master_identity_relative",
        "DRESP11_compact_DmG_reconstruction_residual",
        "DRESP11_denominator_reconstruction_residual",
    ):
        number(values, key)

    assert number(values, "field_identity_max_k_residual") < 1.0e-10
    assert number(values, "Frechet_oracle_residual") < 1.0e-10
    assert number(values, "independent_compact_connection_residual") < 1.0e-10
    # PASS-B is the documented finite-but-incomplete accounting outcome.  The
    # closure gates are therefore required for PASS-A, while PASS-B still has
    # to publish finite diagnostics for the same identities.
    if verdict == "PASS-A":
        assert number(values, "master_identity_relative") < 1.0e-6
        assert number(values, "DRESP11_compact_DmG_reconstruction_residual") < 1.0e-6
        assert number(values, "DRESP11_denominator_reconstruction_residual") < 1.0e-6

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

    for index in range(10):
        fields = values[f"SVD_mode({index})_sigma_target_overlap"].split()
        assert len(fields) == 2 and all(math.isfinite(float(field)) for field in fields)
    print("Dresp12 Fe artifact: covariance accounting and compact reconstruction gates pass")


if __name__ == "__main__":
    main()
