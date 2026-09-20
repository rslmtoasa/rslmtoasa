#!/usr/bin/env python3
"""Independent expected-failure validation for the DRESP-09R Fe gate."""

import math
import sys


def read_artifact(path):
    values = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if "=" not in line or line.lstrip().startswith("#"):
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def number(values, key):
    value = float(values[key])
    if not math.isfinite(value):
        raise AssertionError(f"{key} is not finite")
    return value


def main():
    values = read_artifact(sys.argv[1])
    assert values["accepted_kpoints"] == "64"
    assert values["product_dimension_plus"] == "232"
    assert values["product_dimension_minus"] == "232"
    assert values["classification"] == "PAULI_PROJECTION_INSUFFICIENT_FOR_NATIVE_TANGENT"
    assert values["verdict"] == "BLOCKED"
    assert values["ward_stage"].startswith("NOT RUN")
    assert values["alsda_stage"].startswith("NOT RUN")
    assert values["full_sr_capability"] == "BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED"
    assert number(values, "pauli_direct_vs_compact_adjoint_max_relative") < 2.0e-10
    assert number(values, "native_tangent_vs_pauli_direct_max_relative_frobenius") > 2.0e-10
    print("DRESP-09R Fe artifact: expected native-gate blocker and duality closure verified")


if __name__ == "__main__":
    main()
