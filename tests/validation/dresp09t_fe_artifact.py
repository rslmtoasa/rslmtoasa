#!/usr/bin/env python3
"""Independent checks for the DRESP-09T bounded Fe artifact."""

from __future__ import annotations

import math
import sys


def read(path: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in open(path, encoding="utf-8"):
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def number(values: dict[str, str], key: str) -> float:
    value = float(values[key])
    if not math.isfinite(value):
        raise AssertionError(f"{key} is not finite")
    return value


def main() -> None:
    values = read(sys.argv[1])
    assert values["accepted_kpoints"] == "64"
    assert values["classification"] == "ORTHOGONALIZATION_RESPONSE_REQUIRED"
    assert values["verdict"] == "BLOCKED"
    assert values["orthogonalization_response_required"] == "YES"
    assert number(values, "metric_covariant_residual") > 1.0e-3
    assert number(values, "first_order_h_weighted_residual") > 1.0e-2
    assert number(values, "second_order_native_max_relative_residual") > 1.0e-2
    assert number(values, "density_fixed_basis_relative") > 1.0e-2
    assert number(values, "density_corrected_relative") > number(values, "density_fixed_basis_relative")
    theta = [number(values, f"finite_rotation_theta_{key}_relative") for key in ("1e-2", "5e-3", "2p5e-3")]
    assert theta[0] > theta[1] > theta[2] > 0.0
    assert 3.0 < theta[0] / theta[1] < 5.0
    assert 3.0 < theta[1] / theta[2] < 5.0
    print("Dresp09T Fe artifact: independent radial connection passes; live orthogonalization response remains open")


if __name__ == "__main__":
    main()
