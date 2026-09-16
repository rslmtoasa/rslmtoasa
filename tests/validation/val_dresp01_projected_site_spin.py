#!/usr/bin/env python3
"""DRESP-01 accepted bcc-Fe material gate.

The persisted accepted state contains the production spin-density band
moments, not the in-memory LMTO eigenvectors or radial Pauli arrays.  This
gate therefore checks the material-level d/spd selector arithmetic against
the accepted band-moment artifact and checks the spd result against the
accepted SCF moment.  The Fortran unit test covers the live radial/product
operator contract independently.
"""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FE_DIR = ROOT / "example" / "susceptibility" / "bccFe"
BAND_MOMENTS = FE_DIR / "band_moments.dat"
KSPACE_STATE = FE_DIR / "kspace_scf_state.dat"
RADIAL_STATE = FE_DIR / "radial_ground_state_Fe_1.dat"


def header_value(path: pathlib.Path, key: str) -> str | None:
    prefix = f"# {key} ="
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            return line[len(prefix) :].strip()
    return None


def main() -> int:
    required = (BAND_MOMENTS, KSPACE_STATE, RADIAL_STATE)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        print("DRESP-01 bcc-Fe material gate: SKIP (missing accepted artifacts)")
        print("  " + "\n  ".join(missing))
        return 77

    if header_value(RADIAL_STATE, "accepted") != "true":
        raise RuntimeError("radial ground-state artifact is not accepted")
    if header_value(KSPACE_STATE, "scf_converged") != "T":
        raise RuntimeError("k-space artifact is not an accepted converged state")

    moments: dict[tuple[int, int], float] = {}
    for line in BAND_MOMENTS.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 6:
            raise RuntimeError(f"unexpected band-moment row: {line}")
        site, orbital, spin = (int(fields[index]) for index in range(3))
        moments[(orbital, spin)] = float(fields[3])
        if site != 1:
            raise RuntimeError("bcc-Fe gate expects one accepted site")

    required_channels = {(orbital, spin) for orbital in (1, 2, 3) for spin in (1, 2)}
    if not required_channels.issubset(moments):
        raise RuntimeError("accepted band-moment artifact lacks s/p/d spin channels")

    d_moment = moments[(3, 1)] - moments[(3, 2)]
    spd_moment = sum(moments[(orbital, 1)] - moments[(orbital, 2)] for orbital in (1, 2, 3))
    total_moment = float(header_value(KSPACE_STATE, "scf_moment_muB"))
    electron_sum = sum(moments[(orbital, spin)] for orbital in (1, 2, 3) for spin in (1, 2))
    spd_error = abs(spd_moment - total_moment)

    # The band-moment artifact is printed to 6 decimals; retain a tolerance
    # wider than that rounding while still rejecting a selector mismatch.
    tolerance = 2.0e-5
    if spd_error > tolerance or abs(electron_sum - 8.0) > tolerance:
        raise RuntimeError(
            f"accepted bcc-Fe material gate failed: spd-total={spd_error:.6e}, "
            f"electron-sum={electron_sum:.6e}"
        )

    print("DRESP-01 bcc-Fe material gate: PASS")
    print(f"  Fe d projected moment (mu_B, valence-only)   = {d_moment:.9f}")
    print(f"  Fe spd projected moment (mu_B, valence-only)  = {spd_moment:.9f}")
    print(f"  Fe accepted total moment (mu_B)               = {total_moment:.9f}")
    print(f"  spd-vs-total residual                         = {spd_error:.6e}")
    print(f"  accepted s/p/d electron sum                   = {electron_sum:.9f}")
    print("  core contribution                             = excluded by DRESP-01 policy")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError, TypeError) as error:
        print(f"DRESP-01 bcc-Fe material gate: FAIL: {error}", file=sys.stderr)
        sys.exit(1)
