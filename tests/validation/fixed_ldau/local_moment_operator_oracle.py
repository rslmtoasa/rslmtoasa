#!/usr/bin/env python3
"""Independent finite-H coefficient-space oracle for local moments.

The candidate formula is checked against direct projection of H_LDA u onto
the correlated subspace, not against a second evaluation of the same trace
identity.  This is a deterministic mathematical fixture, not material data.
"""

from __future__ import annotations

import json

import numpy as np


SEED = 20261015
N = 8
NC = 5
WRONG_OPERATOR_FLOOR = 1.0e-4  # Ry^2, declared before evaluating controls


def hermitian(rng: np.random.Generator, n: int, scale: float) -> np.ndarray:
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return scale * (raw + raw.conj().T) / 2.0


def spectrum(h0: np.ndarray, v: np.ndarray, mu: float = 0.02, temp: float = 0.055):
    embed = np.zeros_like(h0)
    embed[:NC, :NC] = v
    h = h0 + embed
    e, u = np.linalg.eigh(h)
    f = 1.0 / (1.0 + np.exp(np.clip((e - mu) / temp, -700.0, 700.0)))
    a = u[:NC, :]
    p = (a * f[None, :]) @ a.conj().T
    m1 = (a * (f * e)[None, :]) @ a.conj().T
    q0 = float(np.trace(p).real)
    c = float(np.trace(m1).real / q0)
    q2 = float(np.sum(f * (e - c) ** 2 * np.sum(np.abs(a) ** 2, axis=0)))
    if np.max(np.abs(h @ u - u * e[None, :])) > 2e-14:
        raise AssertionError("fixture eigensystem residual")
    return e, u, f, p, m1, q0, c, q2


def formula(p: np.ndarray, m1: np.ndarray, v: np.ndarray, q0: float, c: float, q2: float):
    t = float(np.trace(v @ p).real)
    c_local = c - t / q0
    m1c = m1 - c * p
    q2_local = q2 - 2.0 * float(np.trace(v @ m1c).real) \
        + float(np.trace(v @ v @ p).real) - t * t / q0
    return t, c_local, q2_local


def coefficient_oracle(h0: np.ndarray, u: np.ndarray, f: np.ndarray, c_local: float) -> float:
    # Directly use the unperturbed operator on the full eigenvectors. Since
    # H=H0+V and Hu=E u, this evaluates ((E-c')a - V a) without reusing P/M1.
    local_action = (h0 @ u)[:NC, :] - c_local * u[:NC, :]
    return float(np.sum(f[None, :] * np.abs(local_action) ** 2))


def run_case(h0: np.ndarray, v: np.ndarray, *, name: str) -> dict[str, float | str]:
    e, u, f, p, m1, q0, c, q2 = spectrum(h0, v)
    t, c_local, q2_local = formula(p, m1, v, q0, c, q2)
    direct = coefficient_oracle(h0, u, f, c_local)
    p_eigs = np.linalg.eigvalsh(p)
    if p_eigs.min() < -1e-12 or p_eigs.max() > 1.0 + 1e-12:
        raise AssertionError(f"physical P is not PSD/contractive in {name}: {p_eigs}")
    if q2_local < -1e-12:
        raise AssertionError(f"local second moment is negative in {name}: {q2_local}")
    direct_center = float(sum(f[b] * np.vdot(u[:NC, b], (h0 @ u)[:NC, b]).real
                              for b in range(N)) / q0)
    if abs(c_local - direct_center) > 2e-13:
        raise AssertionError(f"first-moment coefficient identity failed in {name}")
    err = abs(q2_local - direct)
    if err > 3e-13:
        raise AssertionError(f"independent direct coefficient oracle failed in {name}: {err:.3e}")
    return {
        "case": name,
        "q0": q0,
        "center_ry": c,
        "tr_vp_ry": t,
        "center_local_ry": c_local,
        "q2_legacy_ry2": q2,
        "q2_local_formula_ry2": q2_local,
        "q2_local_direct_h0u_ry2": direct,
        "coefficient_oracle_abs_error_ry2": err,
        "p_min_eigenvalue": float(p_eigs.min()),
        "p_max_eigenvalue": float(p_eigs.max()),
        "max_abs_v_imag_ry": float(np.max(np.abs(v.imag))),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    h0 = hermitian(rng, N, 0.18)
    h0 += np.diag(np.linspace(-0.22, 0.25, N))
    vraw = hermitian(rng, NC, 0.075)
    vraw += np.diag(np.linspace(-0.04, 0.06, NC))
    assert np.max(np.abs(vraw.imag)) > 1e-3

    primary = run_case(h0, vraw, name="complex_anisotropic_vin")

    # Exact no-U no-op: V=0 makes c'=c and q2'=q2 without cancellation.
    zero = np.zeros((NC, NC), dtype=np.complex128)
    zero_case = run_case(h0, zero, name="V_zero")
    if zero_case["center_local_ry"] != zero_case["center_ry"] or \
       zero_case["q2_local_formula_ry2"] != zero_case["q2_legacy_ry2"]:
        raise AssertionError("V=0 path is not exactly unchanged")

    # Scalar shift control, with a fresh exact eigensystem for H0+vI_d.
    shift = 0.137
    scalar_case = run_case(h0, shift * np.eye(NC, dtype=np.complex128), name="scalar_v_identity")
    scalar_err = abs(float(scalar_case["q2_local_formula_ry2"]) -
                     float(scalar_case["q2_legacy_ry2"]))
    if scalar_err > 3e-13:
        raise AssertionError(f"scalar-shift q2 invariance failed: {scalar_err:.3e}")
    if abs(float(scalar_case["center_local_ry"]) -
           (float(scalar_case["center_ry"]) - shift)) > 3e-13:
        raise AssertionError("scalar-shift center did not move by -v")

    # Wrong-input controls: a stale/output potential and a transposed P must
    # disagree with the direct H0-u coefficient result by a predeclared amount.
    e, u, f, p, m1, q0, c, q2 = spectrum(h0, vraw)
    ref_c = float(primary["center_local_ry"])
    ref_q2 = float(primary["q2_local_direct_h0u_ry2"])
    vout = vraw + hermitian(rng, NC, 0.045)
    _, c_vout, q2_vout = formula(p, m1, vout, q0, c, q2)
    _, c_pt, q2_pt = formula(p.T, m1, vraw, q0, c, q2)
    wrong_vout_center = abs(c_vout - ref_c)
    wrong_vout_q2 = abs(q2_vout - ref_q2)
    wrong_pt_center = abs(c_pt - ref_c)
    wrong_pt_q2 = abs(q2_pt - ref_q2)
    if wrong_vout_q2 < WRONG_OPERATOR_FLOOR or wrong_pt_q2 < WRONG_OPERATOR_FLOOR:
        raise AssertionError(f"wrong-input q2 control too weak: Vout={wrong_vout_q2:.3e}, "
                             f"Ptranspose={wrong_pt_q2:.3e} Ry^2")

    # Deliberately non-PSD P must be rejected by the acceptance predicate.
    p_bad = p.copy()
    p_bad[0, 0] -= 2.0
    bad_eigs = np.linalg.eigvalsh(0.5 * (p_bad + p_bad.conj().T))
    non_psd_rejected = bool(bad_eigs.min() < -1e-3)
    if not non_psd_rejected:
        raise AssertionError("deliberately non-PSD density was not detected")

    print(json.dumps({
        "evidence_class": "A mathematical finite-H oracle; no material/SCF claim",
        "seed": SEED,
        "fixture": {"full_dimension": N, "correlated_dimension": NC,
                    "occupation": "finite-temperature Fermi weights from exact Hermitian eigensystem",
                    "V": "complex Hermitian anisotropic d-block, embedded in H=H0+V"},
        "predeclared_wrong_operator_floor_ry2": WRONG_OPERATOR_FLOOR,
        "cases": [primary, zero_case, scalar_case],
        "scalar_shift_q2_error_ry2": scalar_err,
        "wrong_vout_center_mismatch_ry": wrong_vout_center,
        "wrong_vout_q2_mismatch_ry2": wrong_vout_q2,
        "wrong_transposed_p_center_mismatch_ry": wrong_pt_center,
        "wrong_transposed_p_q2_mismatch_ry2": wrong_pt_q2,
        "non_psd_control_min_eigenvalue": float(bad_eigs.min()),
        "non_psd_control_rejected": non_psd_rejected,
        "direct_identity": "q2_local = sum_n f_n ||[(H0 u_n)_d - c_local a_n]||^2",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
