#!/usr/bin/env python3
"""Independent finite-H reference for the M1.4 local-moment helper.

This module has no production imports. It fixes a deterministic four-case
fixture (two inequivalent sites by two collinear spins), diagonalizes each
complex Hermitian finite-H block, and computes P/M1 by eigenpair sums. The
direct corrected second moment is evaluated from H_LDA acting on the same
eigenvectors, independently of the trace correction formula.
"""
from __future__ import annotations

import json
from pathlib import Path
import argparse
import sys

import numpy as np


SEED = 20261015
N = 8
ND = 5
U_EV = 2.0
RY2EV = 13.605703976
U_RY = U_EV / RY2EV
MU_RY = 0.017
TEMP_RY = 0.061
EPS = np.finfo(np.float64).eps


def hermitian(rng: np.random.Generator, n: int, scale: float) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return scale * (z + z.conj().T) / 2.0


def density(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND))
    q, _ = np.linalg.qr(z)
    occ = np.linspace(0.12, 0.84, ND)
    return (q * occ[None, :]) @ q.conj().T


def gamma(nops: int) -> float:
    x = nops * EPS
    return x / (1.0 - x)


def one_case(rng: np.random.Generator, site: int, spin: int) -> dict:
    h0 = hermitian(rng, N, 0.095)
    h0 += np.diag(np.linspace(-0.31, 0.29, N))
    pin = density(rng)
    vin = U_RY * (0.5 * np.eye(ND) - pin)
    vembed = np.zeros((N, N), dtype=np.complex128)
    vembed[:ND, :ND] = vin
    h = h0 + vembed
    e, u = np.linalg.eigh(h)
    f = 1.0 / (1.0 + np.exp(np.clip((e - MU_RY) / TEMP_RY, -700.0, 700.0)))
    a = u[:ND, :]
    p = (a * f[None, :]) @ a.conj().T
    m1 = (a * (f * e)[None, :]) @ a.conj().T
    q0 = float(np.trace(p).real)
    c = float(np.trace(m1).real / q0)
    q2 = float(np.sum(f * (e - c) ** 2 * np.sum(np.abs(a) ** 2, axis=0)))
    t = float(np.trace(vin @ p).real)
    m1c = m1 - c * p
    vv = vin @ vin
    cprime = c - t / q0
    q2_formula = q2 - 2.0 * float(np.trace(vin @ m1c).real) \
        + float(np.trace(vv @ p).real) - t * t / q0
    h0u = (h0 @ u)[:ND, :]
    direct_vec = h0u - cprime * a
    q2_direct = float(np.sum(f[None, :] * np.abs(direct_vec) ** 2))
    c_direct = float(np.sum(f[None, :] * np.conj(a) * h0u).real / q0)
    vout = U_RY * (0.5 * np.eye(ND) - p)

    def corrected(v: np.ndarray, pp: np.ndarray = p, mm: np.ndarray = m1) -> tuple[float, float]:
        tt = float(np.trace(v @ pp).real)
        cc = c - tt / q0
        mmc = mm - c * pp
        qq = q2 - 2.0 * float(np.trace(v @ mmc).real) \
            + float(np.trace(v @ v @ pp).real) - tt * tt / q0
        return cc, qq

    c_vout, q2_vout = corrected(vout)
    c_pt, q2_pt = corrected(vin, p.T, m1.T)
    scalar_v = 0.137
    c_scalar, q2_scalar = corrected(scalar_v * np.eye(ND, dtype=np.complex128))
    peig = np.linalg.eigvalsh(p)
    p_bad = p - (float(peig.min()) + 0.01) * np.eye(ND)
    p_bad_min_eig = float(np.linalg.eigvalsh(p_bad).min())
    # A deliberately bad input second moment creates q2'=-0.01 Ry^2;
    # this is far outside the precomputed forward-error envelope.
    q2_bad_input = q2 - q2_formula - 0.01

    # Conservative operation-count bound for double-complex traces, products,
    # quotient and direct projected squared norms. The residual term accounts
    # separately for the observed finite-H eigensolver residual.
    resid = h @ u - u * e[None, :]
    resid_local = resid[:ND, :]
    gamma_formula = gamma(16384)
    formula_scale = (abs(q2)
        + 2.0 * float(np.sum(np.abs(vin) * (np.abs(m1.T) + abs(c) * np.abs(p.T))))
        + float(np.sum((np.abs(vin) @ np.abs(vin)) * np.abs(p.T)))
        + abs(t * t / q0)
        + float(np.sum(np.abs(vin) * np.abs(p.T))) * (1.0 + abs(t / q0)))
    direct_intermediate = np.abs(h0) @ np.abs(u)
    direct_intermediate = direct_intermediate[:ND, :] + abs(cprime) * np.abs(a)
    direct_scale = float(np.sum(f[None, :] * direct_intermediate ** 2))
    xnorm = np.linalg.norm(direct_vec, axis=0)
    rnorm = np.linalg.norm(resid_local, axis=0)
    residual_q2_bound = float(np.sum(f * (2.0 * xnorm * rnorm + rnorm * rnorm)))
    q2_error_bound = gamma_formula * (formula_scale + direct_scale) + residual_q2_bound
    c_residual_bound = float(np.sum(f * np.linalg.norm(a, axis=0) * rnorm) / q0)
    c_scale = (abs(c) + abs(t / q0) + abs(c_direct)
        + float(np.sum(f[None, :] * np.abs(a) * (np.abs(h0) @ np.abs(u))[:ND, :]) / q0)
        + float(np.sum(np.abs(vin) * np.abs(p.T)) / q0)
        + float(np.sum(np.abs(m1)) / q0))
    c_error_bound = gamma(4096) * c_scale + c_residual_bound

    return {
        "site": site, "spin": spin,
        "H0_re": h0.real.tolist(), "H0_im": h0.imag.tolist(),
        "Vin_re": vin.real.tolist(), "Vin_im": vin.imag.tolist(),
        "Vout_re": vout.real.tolist(), "Vout_im": vout.imag.tolist(),
        "eigenvalues_ry": e.tolist(), "occupations": f.tolist(),
        "eigenvectors_re": u.real.tolist(), "eigenvectors_im": u.imag.tolist(),
        "P_re": p.real.tolist(), "P_im": p.imag.tolist(),
        "M1_re": m1.real.tolist(), "M1_im": m1.imag.tolist(),
        "q0": q0, "c_ry": c, "q2_ry2": q2, "t_ry": t,
        "cprime_expected_ry": cprime, "cprime_direct_ry": c_direct,
        "q2prime_formula_ry2": q2_formula, "q2prime_direct_ry2": q2_direct,
        "c_error_bound_ry": c_error_bound,
        "q2_error_bound_ry2": q2_error_bound,
        "vout_cprime_ry": c_vout, "vout_q2prime_ry2": q2_vout,
        "p_transpose_cprime_ry": c_pt, "p_transpose_q2prime_ry2": q2_pt,
        "scalar_v_ry": scalar_v, "scalar_cprime_ry": c_scalar,
        "scalar_q2prime_ry2": q2_scalar,
        "scalar_center_expected_ry": c - scalar_v,
        "non_psd_min_eigenvalue": p_bad_min_eig,
        "non_psd_target_min_eigenvalue": -0.01,
        "negative_q2prime_target_ry2": -0.01,
        "negative_q2_input_ry2": q2_bad_input,
        "physical_p_min_eigenvalue": float(peig.min()),
        "max_abs_im_p": float(np.max(np.abs(p.imag))),
        "max_abs_im_vin": float(np.max(np.abs(vin.imag))),
        "max_eigensolver_residual": float(np.max(np.abs(resid))),
        "max_orthogonality_residual": float(np.max(np.abs(u.conj().T @ u - np.eye(N)))),
        "formula_scale_ry2": formula_scale,
        "direct_scale_ry2": direct_scale,
        "residual_q2_bound_ry2": residual_q2_bound,
        "corrected_vout_separation_c_ry": abs(c_vout - cprime),
        "corrected_vout_separation_q2_ry2": abs(q2_vout - q2_formula),
        "corrected_transpose_separation_c_ry": abs(c_pt - cprime),
        "corrected_transpose_separation_q2_ry2": abs(q2_pt - q2_formula),
        "scalar_center_error_ry": abs(c_scalar - (c - scalar_v)),
        "scalar_q2_error_ry2": abs(q2_scalar - q2),
    }


def make_fixture() -> dict:
    rng = np.random.default_rng(SEED)
    cases = [one_case(rng, site, spin) for site in (1, 2) for spin in (1, 2)]
    by_key = {(x["site"], x["spin"]): x for x in cases}
    for spin in (1, 2):
        src, dst = by_key[(1, spin)], by_key[(2, spin)]
        src["site_swap_cprime_ry"] = dst["cprime_expected_ry"]
        src["site_swap_q2prime_ry2"] = dst["q2prime_direct_ry2"]
        src["site_swap_separation_c_ry"] = abs(src["cprime_expected_ry"] - dst["cprime_expected_ry"])
        src["site_swap_separation_q2_ry2"] = abs(src["q2prime_direct_ry2"] - dst["q2prime_direct_ry2"])
    return {
        "schema": "rslmto-a3-m14-independent-finite-h-reference-v1",
        "seed": SEED, "dimension": N, "correlated_dimension": ND,
        "u_eV": U_EV, "ry2ev": RY2EV, "u_ry": U_RY,
        "fermi_mu_ry": MU_RY, "temperature_ry": TEMP_RY,
        "occupations": "Fermi-Dirac weights from the stated finite-H eigenvalues",
        "cases": cases,
        "direct_identity": "q2prime=sum_n f_n ||(H0 U_n)_d-cprime (U_n)_d||^2",
        "error_model": {
            "epsilon": EPS,
            "gamma_formula_ops": 16384,
            "gamma_center_ops": 4096,
            "q2_bound": "gamma_16384*(S_formula+S_direct)+sum_n f_n*(2*||x_n||*||r_n,d||+||r_n,d||^2). S_formula bounds absolute intermediates using |V|(|M1|+|c||P|), |V||V||P|, |V||P| and scalar terms; S_direct=sum_n f_n || |H0||U_n|+|cprime||U_n,d| ||^2. r=H U-U diag(e), x=(H0 U)_d-cprime U_d.",
            "center_bound": "gamma_4096*S_center+sum_n f_n*||U_n,d||*||r_n,d||/q0. S_center includes absolute trace/intermediate scales for c, t/q0 and the direct H0 expectation; q0 is the shared serialized trace input.",
            "operation_count_basis": "16384 real floating operations is a conservative upper bound covering four 5x5 complex matmuls/traces, the 8x8 complex H0*U action (512 complex multiply-accumulates), scalar products/divisions, and 8 weighted five-component squared norms; gamma_n=n*eps/(1-n*eps). The center path is bounded with gamma_4096. Absolute-value intermediate scales avoid relying on cancellation. The explicit eigenpair residual term bounds finite-H identity deviation; orthogonality residual is also recorded as a diagnostic (the algebra uses the same supplied columns and does not assume exact unitarity).",
            "separation_decision": "wrong-control predictions must separate from correct prediction by more than the sum of their case-specific bounds; no additional magnitude cutoff",
        },
    }


def write_fortran_input(obj: dict, path: Path) -> None:
    """Serialize matrices in Fortran column-major order for the unit driver.

    File grammar: header `A3M14 1`, case count, then per case a scalar row
    (site spin q0 c q2 cprime_direct q2prime_direct Ecenter Eq2), followed by
    P, M1, Vin, Vout. Each 5x5 matrix is 25 real/imag pairs with first index
    varying fastest. This avoids any JSON dependency in the Fortran test.
    """
    out = ["A3M14 1", str(len(obj["cases"]))]
    for case in obj["cases"]:
        vals = [case[k] for k in ("site", "spin", "q0", "c_ry", "q2_ry2",
                                   "cprime_direct_ry", "q2prime_direct_ry2",
                                   "c_error_bound_ry", "q2_error_bound_ry2")]
        out.append(" ".join(f"{v:.17g}" for v in vals))
        for label in ("P", "M1", "Vin", "Vout"):
            re = np.asarray(case[f"{label}_re"], dtype=float)
            im = np.asarray(case[f"{label}_im"], dtype=float)
            for j in range(ND):
                for i in range(ND):
                    out.append(f"{re[i,j]:.17g} {im[i,j]:.17g}")
    path.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_json", nargs="?", type=Path)
    parser.add_argument("--fortran-input", type=Path)
    args = parser.parse_args()
    obj = make_fixture()
    payload = json.dumps(obj, indent=2, sort_keys=True) + "\n"
    if args.output_json is None:
        print(payload, end="")
    else:
        args.output_json.write_text(payload)
    if args.fortran_input is not None:
        write_fortran_input(obj, args.fortran_input)
