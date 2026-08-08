#!/usr/bin/env python3
"""Independent WP6b dense and Hermiticity oracles for GBT CCOR."""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TOL = 2.0e-12


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(1.0, float(np.linalg.norm(expected))))


def random_unitary(rng: np.random.Generator) -> np.ndarray:
    raw = rng.normal(size=(2, 2)) + 1.0j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(raw)
    phases = np.diag(r).copy()
    phases /= np.where(np.abs(phases) == 0.0, 1.0, np.abs(phases))
    return q @ np.diag(phases.conj())


def endpoint_link(rng: np.random.Generator, alpha: float) -> np.ndarray:
    ua, ub = random_unitary(rng), random_unitary(rng)
    rz = np.diag([np.exp(-0.5j * alpha), np.exp(0.5j * alpha)])
    return ua.conj().T @ rz @ ub


def endpoint_w(w0: np.ndarray, w1: np.ndarray, sign: float) -> np.ndarray:
    return np.diag(np.concatenate((w0 + sign * w1, w0 - sign * w1))).astype(np.complex128)


def orbital_diag(values: np.ndarray) -> np.ndarray:
    return np.diag(np.tile(values, 2)).astype(np.complex128)


def primitive_d(
    primitive: np.ndarray,
    link: np.ndarray,
    w0_i: np.ndarray,
    w1_i: np.ndarray,
    w0_j: np.ndarray,
    w1_j: np.ndarray,
    sign_i: float,
    sign_j: float,
) -> np.ndarray:
    """Dense spin-major W_i [G_ij tensor S_ij] W_j oracle."""
    return endpoint_w(w0_i, w1_i, sign_i) @ np.kron(link, primitive) @ endpoint_w(
        w0_j, w1_j, sign_j
    )


def scalar_candidate(dmat: np.ndarray, ddot: np.ndarray, a1_i: np.ndarray, a1_j: np.ndarray, lam: float) -> np.ndarray:
    norb = len(a1_i)
    out = np.zeros_like(dmat)
    for si in range(2):
        for sj in range(2):
            for i in range(norb):
                for j in range(norb):
                    row, col = si * norb + i, sj * norb + j
                    out[row, col] = lam * (ddot[row, col] + (a1_i[i] + a1_j[j]) * dmat[row, col])
    return out


def pair_surface_candidate(
    dmat: np.ndarray,
    ddot: np.ndarray,
    a1_i: np.ndarray,
    a1_j: np.ndarray,
    lambda_i: np.ndarray,
    lambda_j: np.ndarray,
) -> np.ndarray:
    norb = len(a1_i)
    out = np.zeros_like(dmat)
    for si in range(2):
        for sj in range(2):
            for i in range(norb):
                for j in range(norb):
                    row, col = si * norb + i, sj * norb + j
                    out[row, col] = 0.5 * (lambda_i[si] + lambda_j[sj]) * ddot[row, col]
                    out[row, col] += (
                        a1_i[i] * lambda_i[si] + a1_j[j] * lambda_j[sj]
                    ) * dmat[row, col]
    return out


def lambda_channels(pair: np.ndarray, sign: float) -> np.ndarray:
    l0, l1 = 0.5 * (pair[0] + pair[1]), 0.5 * (pair[0] - pair[1])
    return np.array([l0 + sign * l1, l0 - sign * l1])


def dense_and_reverse_oracle() -> tuple[float, float, float, float]:
    rng = np.random.default_rng(60622026)
    max_scalar, max_surface = 0.0, 0.0
    max_reverse, max_onsite = 0.0, 0.0
    for norb in (2, 3, 5):
        for _ in range(8):
            s = rng.normal(size=(norb, norb)) + 1.0j * rng.normal(size=(norb, norb))
            sdot = rng.normal(size=(norb, norb)) + 1.0j * rng.normal(size=(norb, norb))
            link = endpoint_link(rng, float(rng.uniform(-2.8, 2.8)))
            sign_i, sign_j = rng.choice([-1.0, 1.0], size=2)
            w0_i, w0_j = rng.uniform(0.2, 1.1, size=(2, norb))
            w1_i, w1_j = rng.uniform(-0.5, 0.5, size=(2, norb))
            a0_i = rng.uniform(0.01, 0.2, size=norb)
            a1_i, a1_j = rng.uniform(-0.3, 0.3, size=(2, norb))
            dmat = primitive_d(s, link, w0_i, w1_i, w0_j, w1_j, sign_i, sign_j)
            ddot = primitive_d(sdot, link, w0_i, w1_i, w0_j, w1_j, sign_i, sign_j)

            lam = float(rng.uniform(-0.9, 0.9))
            scalar = scalar_candidate(dmat, ddot, a1_i, a1_j, lam)
            scalar_dense = lam * (ddot + orbital_diag(a1_i) @ dmat + dmat @ orbital_diag(a1_j))
            max_scalar = max(max_scalar, relerr(scalar, scalar_dense))

            pair = rng.uniform(-1.0, 1.0, size=2)
            li, lj = lambda_channels(pair, sign_i), lambda_channels(pair, sign_j)
            surface = pair_surface_candidate(dmat, ddot, a1_i, a1_j, li, lj)
            li_dense = np.diag(np.repeat(li, norb))
            lj_dense = np.diag(np.repeat(lj, norb))
            surface_dense = 0.5 * (li_dense @ ddot + ddot @ lj_dense)
            surface_dense += li_dense @ orbital_diag(a1_i) @ dmat
            surface_dense += dmat @ orbital_diag(a1_j) @ lj_dense
            max_surface = max(max_surface, relerr(surface, surface_dense))

            # Reverse the primitive directed factors and swap all endpoints.
            d_rev = primitive_d(s.conj().T, link.conj().T, w0_j, w1_j, w0_i, w1_i, sign_j, sign_i)
            ddot_rev = primitive_d(sdot.conj().T, link.conj().T, w0_j, w1_j, w0_i, w1_i, sign_j, sign_i)
            scalar_rev = scalar_candidate(d_rev, ddot_rev, a1_j, a1_i, lam)
            surface_rev = pair_surface_candidate(d_rev, ddot_rev, a1_j, a1_i, lj, li)
            max_reverse = max(max_reverse, relerr(scalar_rev, scalar.conj().T))
            max_reverse = max(max_reverse, relerr(surface_rev, surface.conj().T))

            # Onsite alpha=0/common frame: G=I and the term is diagonal in the
            # shared rotating-frame spin channels.
            wi = endpoint_w(w0_i, w1_i, sign_i)
            onsite_scalar = lam * wi @ orbital_diag(a0_i) @ wi
            onsite_candidate = np.diag(
                lam * np.tile(a0_i, 2) * np.diag(wi) ** 2
            )
            onsite_surface = li_dense @ wi @ orbital_diag(a0_i) @ wi
            onsite_surface_candidate = np.diag(
                np.repeat(li, norb) * np.tile(a0_i, 2) * np.diag(wi) ** 2
            )
            max_onsite = max(max_onsite, relerr(onsite_candidate, onsite_scalar))
            max_onsite = max(max_onsite, relerr(onsite_surface_candidate, onsite_surface))
    return max_scalar, max_surface, max_reverse, max_onsite


def kspace_hermiticity_oracle() -> tuple[float, float]:
    rng = np.random.default_rng(6062026)
    norb = 3
    maxima = [0.0, 0.0]
    for mode in range(2):
        raw0 = rng.normal(size=(2 * norb, 2 * norb)) + 1.0j * rng.normal(size=(2 * norb, 2 * norb))
        onsite = 0.1 * (raw0 + raw0.conj().T)
        forward = 0.08 * (rng.normal(size=raw0.shape) + 1.0j * rng.normal(size=raw0.shape))
        reverse = forward.conj().T
        for k in np.linspace(-np.pi, np.pi, 41):
            hk = onsite + forward * np.exp(1.0j * k) + reverse * np.exp(-1.0j * k)
            maxima[mode] = max(maxima[mode], relerr(hk, hk.conj().T))
    return maxima[0], maxima[1]


def source_contract() -> list[str]:
    failures: list[str] = []
    ccor = (ROOT / "source/hamiltonian_ccor.f90").read_text(encoding="utf-8")
    build = (ROOT / "source/hamiltonian_build.f90").read_text(encoding="utf-8")
    lifecycle = (ROOT / "source/reciprocal_lifecycle.f90").read_text(encoding="utf-8")
    parent = (ROOT / "source/hamiltonian.f90").read_text(encoding="utf-8")
    forbidden = ("alpha_ss", "kx_ss", "ky_ss", "hx_ss", "hy_ss", "ccor_apply_spin_spiral")
    if any(token in ccor or token in parent for token in forbidden):
        failures.append("a deleted full-angle/completed-object CCOR rotation remains")
    start = ccor.index("module subroutine build_ccor_pair_block_gbt")
    end = ccor.index("end subroutine build_ccor_pair_block_gbt", start)
    gbt = ccor[start:end]
    if gbt.count("call gbt_lift_orbital_block") != 2:
        failures.append("GBT CCOR must link primitive S and Sdot exactly once each")
    if "raw_orb, raw_sdot, ccor_spin" not in build:
        failures.append("common main-Hamiltonian endpoint link is not passed to CCOR")
    if "call ccor_select_endpoint_moments" not in ccor or "this%texture_moments" not in ccor:
        failures.append("ordinary explicit_texture CCOR does not use site endpoint moments")
    if "contracted CCOR pair builder" not in ccor:
        failures.append("GBT is not rejected from the ordinary completed-object CCOR path")
    if "with CCOR is guarded" in build or "GBT with CCOR is unsupported" in lifecycle:
        failures.append("obsolete blanket GBT+CCOR guard remains after all enabled-term oracles")
    if "Invalid ccor_vmt_mode" not in build or "ccor_2c requires lattice%sdot" not in ccor:
        failures.append("unsupported/missing CCOR inputs do not fail early")
    if "gbt_single_q with SOC is unsupported" not in build:
        failures.append("mixed GBT+CCOR+SOC is not covered by the early SOC rejection")
    return failures


def main() -> int:
    scalar, surface, reverse, onsite = dense_and_reverse_oracle()
    kscalar, ksurface = kspace_hermiticity_oracle()
    failures = source_contract()
    print(f"scalar CCOR dense relative error          : {scalar:.3e}")
    print(f"pair-surface CCOR dense relative error    : {surface:.3e}")
    print(f"directed reverse-bond relative error      : {reverse:.3e}")
    print(f"onsite rotating-frame relative error      : {onsite:.3e}")
    print(f"scalar CCOR k-space Hermiticity error     : {kscalar:.3e}")
    print(f"pair-surface k-space Hermiticity error    : {ksurface:.3e}")
    checked = [scalar, surface, reverse, onsite, kscalar, ksurface]
    if max(checked) >= TOL:
        failures.append(f"algebraic maximum {max(checked):.3e} is not below {TOL:.1e}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("WP6b CCOR oracle PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
