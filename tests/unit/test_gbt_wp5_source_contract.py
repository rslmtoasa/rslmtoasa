#!/usr/bin/env python3
"""Static WP5 guard: reciprocal GBT reconstruction must stay deleted."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source"


def main() -> None:
    failures: list[str] = []
    for path in SOURCE.rglob("*.f90"):
        text = path.read_text(encoding="utf-8")
        if "fourier_transform_gbt" in text:
            failures.append(f"legacy reciprocal GBT transform remains in {path.relative_to(ROOT)}")

        if "gbt_kspace" in text and path.relative_to(ROOT).as_posix() not in {
            "source/hamiltonian.f90",
            "source/hamiltonian_build.f90",
            "source/include_codes/namelists/hamiltonian.f90",
        }:
            failures.append(f"gbt_kspace escaped its deprecated input boundary in {path.relative_to(ROOT)}")

    reciprocal_fourier = (SOURCE / "reciprocal_fourier.f90").read_text(encoding="utf-8")
    hamiltonian_build = (SOURCE / "hamiltonian_build.f90").read_text(encoding="utf-8")
    if "call this%fourier_transform_array(this%hamiltonian%ee" not in reciprocal_fourier:
        failures.append("reciprocal ee assembly does not call the ordinary Fourier transform")
    if "hamiltonian%gbt_kspace" in reciprocal_fourier:
        failures.append("reciprocal assembly still branches on gbt_kspace")
    compatibility_block = hamiltonian_build[
        hamiltonian_build.find("if (this%gbt_kspace)") : hamiltonian_build.find(
            "this%js_alpha = js_alpha"
        )
    ]
    if "this%magnetic_representation =" in compatibility_block:
        failures.append("gbt_kspace still changes magnetic representation")

    if failures:
        raise SystemExit("\n".join(f"FAIL: {item}" for item in failures))
    print("WP5 source contract PASS: legacy transforms absent; gbt_kspace is input-only")


if __name__ == "__main__":
    main()
