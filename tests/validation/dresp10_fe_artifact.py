#!/usr/bin/env python3
"""Independent checks for the DRESP-10 RAW full-spatial artifact."""

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


def main() -> None:
    values = read(sys.argv[1])
    assert values["source_space_complete"].lower() == "f"
    assert values["source_probe_nonzero"].lower() == "t"
    assert values["target_primary"] == "P3_PAULI_MOMENT_DENSITY"
    assert values["mixed_K_contact_operator"] == "NOT_ASSEMBLED"
    assert values["dynamic_spectra"] == "NOT_RUN"
    assert values["raw_denominator"].startswith("NOT_ASSEMBLED")
    assert values["verdict"] == "BLOCKED"
    assert values["classification"] == "SR_SOURCE_SPACE_INCOMPLETE"
    assert number(values, "response_lmax") == 4
    assert number(values, "product_dimension") == 348
    assert number(values, "Kxc_identity_max_abs") < 1.0e-10
    assert number(values, "static_frechet_action_relative") < 1.0e-10
    assert number(values, "static_frechet_pseudorandom_relative") < 1.0e-10
    assert number(values, "static_frechet_physical_source_relative") < 1.0e-10
    assert number(values, "field_insertion_relative") > 1.0e-3
    assert number(values, "source_span_physical_global") > 1.0e-3
    assert number(values, "source_span_delta_O_global") > 1.0e-3
    assert values["DRESP09Y_contact_regression"] == "CLOSED"
    assert number(values, "contact_angular_rank_recomposition_relative") < 1.0e-10
    assert number(values, "complete_response_vs_SRspin2") < 3.0e-8
    assert number(values, "Pauli_complete_vs_P3") < 3.0e-8
    assert number(values, "legacy_l0_field_insertion_relative") > 1.0e-3
    print("Dresp10 Fe artifact: source-space gate passes; RAW Ward gates correctly block")


if __name__ == "__main__":
    main()
