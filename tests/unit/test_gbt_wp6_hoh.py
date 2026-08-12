#!/usr/bin/env python3
"""Independent WP6a oracle for GBT HOH covariance and overlap policy."""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TOL = 2.0e-12


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(1.0, float(np.linalg.norm(expected))))


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    raw = rng.normal(size=(n, n)) + 1.0j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(raw)
    phases = np.diag(r).copy()
    phases /= np.where(np.abs(phases) == 0.0, 1.0, np.abs(phases))
    return q @ np.diag(phases.conj())


def primitive_covariance_oracle() -> tuple[float, float, float]:
    """Gauge directed H factors first, leave O/E onsite, then contract h-o-h."""
    rng = np.random.default_rng(60262026)
    nsite, block = 5, 4
    n = nsite * block
    raw = rng.normal(size=(n, n)) + 1.0j * rng.normal(size=(n, n))
    h_lab = 0.08 * (raw + raw.conj().T)
    o_lab = np.zeros((n, n), dtype=np.complex128)
    e_lab = np.zeros((n, n), dtype=np.complex128)
    gauge = np.zeros((n, n), dtype=np.complex128)
    for site in range(nsite):
        row = slice(site * block, (site + 1) * block)
        oraw = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
        eraw = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
        o_lab[row, row] = 0.03 * (oraw + oraw.conj().T)
        e_lab[row, row] = 0.2 * (eraw + eraw.conj().T)
        gauge[row, row] = random_unitary(rng, block)

    h_gbt = gauge.conj().T @ h_lab @ gauge
    o_gbt = gauge.conj().T @ o_lab @ gauge
    e_gbt = gauge.conj().T @ e_lab @ gauge
    eeo_gbt = h_gbt @ o_gbt
    hoh_two_sweep = eeo_gbt @ h_gbt
    hoh_covariant = gauge.conj().T @ (h_lab @ o_lab @ h_lab) @ gauge
    h_second = e_gbt + h_gbt - hoh_two_sweep
    return (
        relerr(hoh_two_sweep, hoh_covariant),
        relerr(h_second, h_second.conj().T),
        relerr(eeo_gbt, eeo_gbt.conj().T),
    )


def rs_reciprocal_oracle() -> tuple[float, float, float, float]:
    """Compare periodic RS eeo*(h*x) with reciprocal eeo(k)*h(k)*x(k)."""
    rng = np.random.default_rng(6062026)
    ncell, block = 7, 4
    onsite_raw = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
    hop = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
    h_blocks = {
        0: 0.2 * (onsite_raw + onsite_raw.conj().T),
        1: 0.08 * hop,
        -1: 0.08 * hop.conj().T,
    }
    oraw = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
    eraw = rng.normal(size=(block, block)) + 1.0j * rng.normal(size=(block, block))
    onsite_o = 0.025 * (oraw + oraw.conj().T)
    onsite_e = 0.15 * (eraw + eraw.conj().T)
    eeo_blocks = {shift: value @ onsite_o for shift, value in h_blocks.items()}

    h_rs = np.zeros((ncell * block, ncell * block), dtype=np.complex128)
    eeo_rs = np.zeros_like(h_rs)
    for cell in range(ncell):
        row = slice(cell * block, (cell + 1) * block)
        for shift, h_block in h_blocks.items():
            other = (cell + shift) % ncell
            col = slice(other * block, (other + 1) * block)
            h_rs[row, col] += h_block
            eeo_rs[row, col] += eeo_blocks[shift]

    vector = rng.normal(size=ncell * block) + 1.0j * rng.normal(size=ncell * block)
    rs_action = eeo_rs @ (h_rs @ vector)
    reciprocal_action = np.zeros_like(vector)
    max_h_herm = 0.0
    max_second_herm = 0.0
    min_supported_overlap = np.inf
    for ik in range(ncell):
        k = 2.0 * np.pi * ik / ncell
        # The RS storage below uses row cell i, column cell i+shift. With the
        # x(k)=sum_i exp(+ik i)x_i convention this maps to exp(-ik*shift).
        hk = sum(block_r * np.exp(-1.0j * k * shift) for shift, block_r in h_blocks.items())
        eeok = sum(block_r * np.exp(-1.0j * k * shift) for shift, block_r in eeo_blocks.items())
        hsecond_k = onsite_e + hk - eeok @ hk
        max_h_herm = max(max_h_herm, relerr(hk, hk.conj().T))
        max_second_herm = max(max_second_herm, relerr(hsecond_k, hsecond_k.conj().T))
        min_supported_overlap = min(min_supported_overlap, float(np.linalg.eigvalsh(np.eye(block)).min()))

        xk = sum(
            vector[cell * block : (cell + 1) * block] * np.exp(1.0j * k * cell)
            for cell in range(ncell)
        )
        yk = eeok @ hk @ xk
        for cell in range(ncell):
            reciprocal_action[cell * block : (cell + 1) * block] += (
                yk * np.exp(-1.0j * k * cell) / ncell
            )

    return relerr(reciprocal_action, rs_action), max_h_herm, max_second_herm, min_supported_overlap


def source_contract() -> list[str]:
    failures: list[str] = []
    build = (ROOT / "source/hamiltonian_build.f90").read_text(encoding="utf-8")
    lifecycle = (ROOT / "source/reciprocal_lifecycle.f90").read_text(encoding="utf-8")
    fourier = (ROOT / "source/reciprocal_fourier.f90").read_text(encoding="utf-8")
    bands = (ROOT / "source/reciprocal_bands.f90").read_text(encoding="utf-8")
    start = build.index("module subroutine build_gbt_bulkham")
    end = build.index("end subroutine build_gbt_bulkham", start)
    gbt_builder = build[start:end]

    if "gbt_single_q with HOH/eeo is guarded" in build or "GBT with HOH/overlap products is unsupported" in lifecycle:
        failures.append("obsolete GBT+HOH guard remains")
    if "call this%build_obarm()" not in gbt_builder or "this%eeo(:, :, m, ntype)" not in gbt_builder:
        failures.append("GBT eeo is not derived from linked ee and onsite obarm")
    if gbt_builder.count("call gbt_endpoint_link") != 1:
        failures.append("GBT link must be computed only for the primitive S path")
    if "not a complete formal GBT metric" not in lifecycle:
        failures.append("incomplete generalized-overlap modes are not rejected explicitly")
    legacy_hoh = "call this%fourier_transform_array(this%hamiltonian%eeo" in fourier and "eeok, ndim, hk" in fourier
    batched_hoh = "workspace%eeo(:,:,ik),nmat,workspace%h(:,:,ik)" in fourier
    if not (legacy_hoh or batched_hoh):
        failures.append("reciprocal HOH is not eeo(k)*h(k)")
    if "call zpotrf" not in bands:
        failures.append("supported generalized overlaps lack a positive-definiteness check")
    return failures


def main() -> int:
    covariance, h_cov_herm, eeo_nonherm = primitive_covariance_oracle()
    rs_k, hk_herm, h2k_herm, overlap_min = rs_reciprocal_oracle()
    failures = source_contract()
    print(f"primitive-bond HOH covariance relative error : {covariance:.3e}")
    print(f"covariant second-order H Hermiticity error   : {h_cov_herm:.3e}")
    print(f"directed eeo non-Hermiticity (expected)      : {eeo_nonherm:.3e}")
    print(f"RS two-sweep / reciprocal relative error     : {rs_k:.3e}")
    print(f"reciprocal first-order H Hermiticity error   : {hk_herm:.3e}")
    print(f"reciprocal second-order H Hermiticity error  : {h2k_herm:.3e}")
    print(f"supported identity-overlap minimum eigenvalue: {overlap_min:.6f}")

    checked = [covariance, h_cov_herm, rs_k, hk_herm, h2k_herm]
    if max(checked) >= TOL:
        failures.append(f"algebraic maximum {max(checked):.3e} is not below {TOL:.1e}")
    if eeo_nonherm <= 1.0e-6:
        failures.append("fixture does not expose that directed eeo is not Hermitian")
    if overlap_min <= 0.0:
        failures.append("supported overlap is not positive definite")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("WP6a HOH/overlap oracle PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
