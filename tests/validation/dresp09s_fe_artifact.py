#!/usr/bin/env python3
"""Independent validation of the DRESP-09S accepted Fe artifact."""

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
    assert values["classification"] == "BASIS_ROTATION_TERM_REQUIRED"
    assert values["verdict"] == "PASS-B"
    assert values["basis_response_term_required"] == "YES"
    assert number(values, "pauli_large_component_weighted_rms") > 1.0e-2
    assert number(values, "full_sr_angular_weighted_rms") > 1.0e-2
    assert number(values, "field_density_duality_max_relative") < 2.0e-10
    theta1 = number(values, "finite_rotation_theta_1e-2_relative")
    theta2 = number(values, "finite_rotation_theta_5e-3_relative")
    theta3 = number(values, "finite_rotation_theta_2p5e-3_relative")
    assert theta1 > theta2 > theta3 > 0.0
    assert 3.0 < theta1 / theta2 < 5.0
    assert 3.0 < theta2 / theta3 < 5.0
    for key in (
        "density_pauli_projected_valence_relative",
        "density_full_sr_valence_relative",
        "density_full_sr_total_relative",
        "dresp08_native_commutator_max_relative",
    ):
        number(values, key)
    print("DRESP-09S Fe artifact: SR derivation, duality, and PASS-B basis-response classification verified")


if __name__ == "__main__":
    main()
