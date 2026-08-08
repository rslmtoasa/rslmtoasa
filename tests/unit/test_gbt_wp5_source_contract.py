#!/usr/bin/env python3
"""Static WP5 guard: reciprocal GBT reconstruction must stay deleted.

gbt_kspace itself was removed outright in WP10 (it had carried no physics
role since WP5); its absence is asserted here rather than in a dedicated
WP10 test because this file already owns the "legacy transform is gone"
contract.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source"


def main() -> None:
    failures: list[str] = []
    for path in SOURCE.rglob("*.f90"):
        text = path.read_text(encoding="utf-8")
        if "fourier_transform_gbt" in text:
            failures.append(f"legacy reciprocal GBT transform remains in {path.relative_to(ROOT)}")

        if "gbt_kspace" in text:
            failures.append(f"gbt_kspace was removed in WP10 but still appears in {path.relative_to(ROOT)}")

    reciprocal_fourier = (SOURCE / "reciprocal_fourier.f90").read_text(encoding="utf-8")
    if "call this%fourier_transform_array(this%hamiltonian%ee" not in reciprocal_fourier:
        failures.append("reciprocal ee assembly does not call the ordinary Fourier transform")

    if failures:
        raise SystemExit("\n".join(f"FAIL: {item}" for item in failures))
    print("WP5 source contract PASS: legacy transforms absent; gbt_kspace removed")


if __name__ == "__main__":
    main()
