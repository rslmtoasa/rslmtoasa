#!/usr/bin/env python3
"""WP6c GBT velocity/operator algebra and source-contract oracle."""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TOL = 2.0e-12


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(1.0, float(np.linalg.norm(expected))))


def velocity_oracle() -> tuple[float, float, float, float]:
    """Check reverse bonds, Fourier Hermiticity, dH/dk, and current operators."""
    rng = np.random.default_rng(60632026)
    n = 4
    raw = rng.normal(size=(n, n)) + 1.0j * rng.normal(size=(n, n))
    onsite = 0.15 * (raw + raw.conj().T)
    forward = 0.08 * (rng.normal(size=(n, n)) + 1.0j * rng.normal(size=(n, n)))
    reverse = forward.conj().T
    # rij = -R, so v(R)=(dir.rij)/i H(R).  For R=+1 this is +i H(+1).
    velocity_forward = 1.0j * forward
    velocity_reverse = -1.0j * reverse
    reverse_error = relerr(velocity_reverse, velocity_forward.conj().T)

    max_hermiticity = 0.0
    max_derivative = 0.0
    max_current = 0.0
    spin = np.diag([1.0, 1.0, -1.0, -1.0]).astype(np.complex128)
    orbital = np.diag([-1.5, -0.5, 0.5, 1.5]).astype(np.complex128)
    dk = 1.0e-5
    for k in np.linspace(-np.pi, np.pi, 31):
        phase = np.exp(1.0j * k)
        h_k = onsite + forward * phase + reverse * phase.conjugate()
        v_k = velocity_forward * phase + velocity_reverse * phase.conjugate()
        central = (
            forward * np.exp(1.0j * (k + dk))
            + reverse * np.exp(-1.0j * (k + dk))
            - forward * np.exp(1.0j * (k - dk))
            - reverse * np.exp(-1.0j * (k - dk))
        ) / (2.0 * dk)
        js_k = 0.5 * (spin @ v_k + v_k @ spin)
        jl_k = 0.5 * (orbital @ v_k + v_k @ orbital)
        max_hermiticity = max(max_hermiticity, relerr(v_k, v_k.conj().T))
        max_derivative = max(max_derivative, relerr(v_k, central))
        max_current = max(max_current, relerr(js_k, js_k.conj().T), relerr(jl_k, jl_k.conj().T))
    return reverse_error, max_hermiticity, max_derivative, max_current


def routine_body(text: str, name: str) -> str:
    start = text.index(f"module subroutine {name}")
    end = text.index(f"end subroutine {name}", start)
    return text[start:end]


def source_contract() -> list[str]:
    failures: list[str] = []
    build = (ROOT / "source/hamiltonian_build.f90").read_text(encoding="utf-8")
    velocity = routine_body(build, "build_realspace_velocity_operators")
    spin_torque = routine_body(build, "build_realspace_spin_torque_operators")
    orbital_torque = routine_body(build, "build_realspace_orbital_torque_operators")
    if "this%ee(:, :, m, ntype)" not in velocity:
        failures.append("velocity operator is not derived from the already linked ee block")
    if "this%obarm(:, :, ji)" not in velocity:
        failures.append("HOH companion velocity is not formed with the endpoint obarm")
    if "q_ss" in velocity or "gbt_" in velocity:
        failures.append("velocity builder adds a q/GBT phase instead of inheriting ee's gauge")
    if "gbt_single_q with cond_type=spin_torque is unsupported" not in spin_torque:
        failures.append("GBT spin-torque response lacks its explicit rejection")
    if "this%ee" not in orbital_torque or "this%lsham" not in orbital_torque:
        failures.append("orbital torque no longer uses the documented ee/lsham operator split")
    return failures


def main() -> int:
    reverse, hermiticity, derivative, current = velocity_oracle()
    failures = source_contract()
    print(f"velocity reverse-bond relative error : {reverse:.3e}")
    print(f"velocity k-space Hermiticity error  : {hermiticity:.3e}")
    print(f"velocity dH/dk central-difference    : {derivative:.3e}")
    print(f"current-operator Hermiticity error   : {current:.3e}")
    if max(reverse, hermiticity, current) >= TOL:
        failures.append(f"exact algebraic maximum is not below {TOL:.1e}")
    if derivative >= 5.0e-9:
        failures.append("central-difference dH/dk error exceeds 5.0e-9")
    if failures:
        raise SystemExit("\n".join(f"FAIL: {failure}" for failure in failures))
    print("WP6c velocity/operator oracle PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
