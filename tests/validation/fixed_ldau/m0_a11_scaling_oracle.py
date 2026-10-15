#!/usr/bin/env python3
"""Independent A1.1 signed scalar-shift scaling analysis.

Consumes a hash-bound KS fixed-U trace after the source-connected replay has
recorded the eight signed centre shifts.  It derives all response values from
raw and shifted SUMEV rows, checks that only the d-shell centre changed, and
compares the reproduced no-U state with the accepted M0.2 capture.  No
production moment or fit result is trusted as input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

import local_moment_capture_oracle as m0


V_MAGNITUDES = (0.0025, 0.005, 0.01, 0.02)
SIGNED_V = (-0.02, -0.01, -0.005, -0.0025, 0.0025, 0.005, 0.01, 0.02)
STAGE_BY_V = {
    -0.02: "A1_SHIFT_01", -0.01: "A1_SHIFT_02", -0.005: "A1_SHIFT_03", -0.0025: "A1_SHIFT_04",
    0.0025: "A1_SHIFT_05", 0.005: "A1_SHIFT_06", 0.01: "A1_SHIFT_07", 0.02: "A1_SHIFT_08",
}
# The zero-shift row reuses the existing exact RAW_PL replay snapshot.  Its
# distinct kind (`m0_a11_case`) prevents ambiguity while avoiding a second
# serialization of the baseline arrays.
RAW_STAGE = "RAW_PL"
SLOPE_INTERCEPT_TOL = 1.0e-6
M1_STATE_TOL_RY = 5.0e-12
MOMENT_STATE_TOL = 1.0e-10
STATE_VALUE_TOL = 1.0e-10


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _one(rows: list[dict[str, str]], *, kind: str, stage: str, iteration: int,
         site: int, type_id: int, ell: int, spin: int, order: int | None = None) -> dict[str, str]:
    found = [r for r in rows if r["kind"] == kind and r["stage"] == stage
             and m0._i(r, "iteration") == iteration and m0._i(r, "site") == site
             and m0._i(r, "type") == type_id and m0._i(r, "l") == ell
             and m0._i(r, "spin") == spin
             and (order is None or m0._i(r, "order") == order)]
    if len(found) != 1:
        raise ValueError(f"expected one {kind}/{stage} row for site,type,l,spin,order="
                         f"{site},{type_id},{ell},{spin},{order}; found {len(found)}")
    return found[0]


def _row_map(rows: list[dict[str, str]], kind: str, stage: str, iteration: int,
             key_fields: tuple[str, ...], value_field: str = "scalar") -> dict[tuple[int, ...], float]:
    out: dict[tuple[int, ...], float] = {}
    for r in rows:
        if r["kind"] != kind or r["stage"] != stage or m0._i(r, "iteration") != iteration:
            continue
        key = tuple(m0._i(r, field) for field in key_fields)
        if key in out:
            raise ValueError(f"duplicate {kind}/{stage} key {key}")
        out[key] = m0._f(r, value_field)
    if not out:
        raise ValueError(f"missing {kind}/{stage} rows")
    return out


def _assert_same_map(label: str, expected: dict[tuple[int, ...], float],
                     actual: dict[tuple[int, ...], float], *, exact: bool = True,
                     tol: float = 0.0) -> float:
    if expected.keys() != actual.keys():
        raise ValueError(f"{label}: row keys differ ({len(expected)} vs {len(actual)})")
    maximum = max(abs(expected[k] - actual[k]) for k in expected)
    if (exact and maximum != 0.0) or (not exact and maximum > tol):
        raise ValueError(f"{label}: frozen values differ by {maximum:.6e}")
    return maximum


def _manifest_hash(manifest: dict[str, Any], basename: str) -> str:
    hits = [str(digest) for name, digest in manifest.get("input_sha256", {}).items()
            if Path(name).name == basename]
    if len(hits) != 1:
        raise ValueError(f"manifest must contain exactly one {basename} hash")
    return hits[0]


def _check_reference_state(candidate_dir: Path, candidate_manifest: dict[str, Any],
                           candidate_rows: list[dict[str, str]], reference_dir: Path,
                           reference_manifest: dict[str, Any], reference_rows: list[dict[str, str]]) -> dict[str, Any]:
    if str(candidate_manifest.get("backend", "")).upper() != "KS":
        raise ValueError("A1.1 oracle is KS-only")
    if str(reference_manifest.get("backend", "")).upper() != "KS":
        raise ValueError("accepted reference capture must be KS")
    for name in ("input.nml.pre", "Fe.nml.pre"):
        if _manifest_hash(candidate_manifest, name) != _manifest_hash(reference_manifest, name):
            raise ValueError(f"candidate {name} does not match the accepted M0.2 pre-run input")
    for field in ("nsp", "lmax", "metric", "kpoint_count", "temperature_ry", "electron_count"):
        if candidate_manifest.get(field) != reference_manifest.get(field):
            raise ValueError(f"candidate/reference manifest mismatch in {field}")
    iteration = int(candidate_manifest["iteration"])
    ref_iteration = int(reference_manifest["iteration"])
    if iteration != ref_iteration:
        raise ValueError(f"A1 candidate iteration {iteration} is not accepted reference iteration {ref_iteration}")

    # Independent eigenpair reconstruction provides P because this no-U
    # capture intentionally marks the live LDM unavailable.
    p_candidate, _ = m0.independent_spectral_moments(candidate_manifest, candidate_rows)
    p_reference, _ = m0.independent_spectral_moments(reference_manifest, reference_rows)
    if p_candidate.keys() != p_reference.keys():
        raise ValueError("candidate/reference independently reconstructed P keys differ")
    max_p = max(float(np.max(np.abs(p_candidate[k][0] - p_reference[k][0]))) for k in p_candidate)
    if max_p > M1_STATE_TOL_RY:
        raise ValueError(f"accepted-state P difference {max_p:.3e} exceeds {M1_STATE_TOL_RY:.3e}")

    max_m1 = 0.0
    m1_keys = sorted({(m0._i(r, "site"), m0._i(r, "l"), m0._i(r, "spin"))
                      for r in candidate_rows if r["kind"] == "ldm_m1" and r["stage"] == "M1_OUT"})
    ref_m1_keys = sorted({(m0._i(r, "site"), m0._i(r, "l"), m0._i(r, "spin"))
                          for r in reference_rows if r["kind"] == "ldm_m1" and r["stage"] == "M1_OUT"})
    if not m1_keys or m1_keys != ref_m1_keys:
        raise ValueError("candidate/reference M1_OUT shell keys differ or are absent")
    for site, ell, spin in m1_keys:
        dim = 2 * ell + 1
        a = m0.matrix(candidate_rows, kind="ldm_m1", stage="M1_OUT", iteration=iteration,
                      site=site, ell=ell, spin=spin, dim=dim)
        b = m0.matrix(reference_rows, kind="ldm_m1", stage="M1_OUT", iteration=ref_iteration,
                      site=site, ell=ell, spin=spin, dim=dim)
        max_m1 = max(max_m1, float(np.max(np.abs(a - b))))
    if max_m1 > M1_STATE_TOL_RY:
        raise ValueError(f"accepted-state M1 difference {max_m1:.3e} Ry exceeds {M1_STATE_TOL_RY:.3e}")

    max_moment = 0.0
    # q0, q0*c, and q2 are captured from the same P_OUT radial moments.
    key_fields = ("site", "type", "l", "spin", "order")
    a = _row_map(candidate_rows, "radial_moment", "P_OUT", iteration, key_fields)
    b = _row_map(reference_rows, "radial_moment", "P_OUT", ref_iteration, key_fields)
    max_moment = _assert_same_map("accepted q0/c/q2", a, b, exact=False, tol=MOMENT_STATE_TOL)

    max_pl = 0.0
    a = _row_map(candidate_rows, "pl_map_input", "CALCULATE_PL", iteration,
                 ("site", "type", "l", "spin", "order"))
    b = _row_map(reference_rows, "pl_map_input", "CALCULATE_PL", ref_iteration,
                 ("site", "type", "l", "spin", "order"))
    max_pl = _assert_same_map("accepted calculate_pl inputs/outputs", a, b, exact=False, tol=STATE_VALUE_TOL)

    frozen_kinds = ("m0_replay_rofi", "m0_replay_v_input", "m0_replay_rho_input",
                    "m0_replay_ec_seed", "m0_replay_ev_seed", "m0_replay_ql_input",
                    "m0_replay_pl_input", "m0_replay_pl_return")
    max_frozen = 0.0
    for kind in frozen_kinds:
        a_stage = RAW_STAGE
        b_stage = "RAW_PL"
        a = _row_map(candidate_rows, kind, a_stage, iteration,
                     ("site", "type", "l", "spin", "row", "col", "order"))
        b = _row_map(reference_rows, kind, b_stage, ref_iteration,
                     ("site", "type", "l", "spin", "row", "col", "order"))
        max_frozen = max(max_frozen, _assert_same_map(f"accepted frozen {kind}", a, b, exact=True))

    return {
        "accepted_state": "PASS",
        "reference_iteration": ref_iteration,
        "candidate_iteration": iteration,
        "max_abs_P_difference": max_p,
        "max_abs_M1_difference_ry": max_m1,
        "max_abs_radial_moment_difference": max_moment,
        "max_abs_calculate_pl_difference": max_pl,
        "max_abs_frozen_replay_input_difference": max_frozen,
        "tolerances": {
            "P_and_M1": M1_STATE_TOL_RY,
            "q0_c_q2": MOMENT_STATE_TOL,
            "calculate_pl": STATE_VALUE_TOL,
            "frozen_arrays_seeds_and_replay_inputs": "bitwise equal as serialized doubles",
        },
    }


def _fit(x: np.ndarray, y: np.ndarray, degree: int) -> np.ndarray:
    design = np.column_stack([x**power for power in range(degree + 1)])
    coeff, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coeff


def _eval(coeff: np.ndarray, x: float) -> float:
    return float(sum(coeff[p] * x**p for p in range(len(coeff))))


def _series_report(values: dict[float, dict[str, float]], q0: float, floor: float) -> dict[str, Any]:
    deltas = {v: item["delta_sumev_ry"] for v, item in values.items()}
    residual = {v: item["r_ry"] for v, item in values.items()}
    magnitudes = np.asarray(V_MAGNITUDES, dtype=float)
    slopes: dict[float, float] = {}
    e: dict[float, float] = {}
    for v in V_MAGNITUDES:
        slopes[v] = (deltas[v] - deltas[-v]) / (2.0 * v)
        e[v] = slopes[v] / q0 - 1.0
    e_values = np.asarray([e[v] for v in V_MAGNITUDES], dtype=float)
    e_coeff = _fit(magnitudes**2, e_values, 2)
    e_linear = _fit(magnitudes**2, e_values, 1)
    e_loo = []
    for drop in range(len(magnitudes)):
        keep = np.arange(len(magnitudes)) != drop
        e_loo.append(_fit(magnitudes[keep]**2, e_values[keep], 2)[0])
    u_e_report = max(abs(e_coeff[0] - e_linear[0]), *(abs(e_coeff[0] - x) for x in e_loo))
    e_monotonic = all(abs(e[magnitudes[i]]) < abs(e[magnitudes[i + 1]])
                      for i in range(len(magnitudes) - 1))
    # Coefficients are stored ascending by power.
    signed_v = np.asarray(SIGNED_V, dtype=float)
    signed_r = np.asarray([residual[float(v)] for v in SIGNED_V], dtype=float)
    r_coeff = _fit(signed_v, signed_r, 3)
    r_quadratic = _fit(signed_v, signed_r, 2)
    log_defined = bool(np.all(np.isfinite(signed_r)) and np.all(np.abs(signed_r) > 0.0))
    if log_defined:
        log_order, log_intercept = np.polyfit(np.log(np.abs(signed_v)), np.log(np.abs(signed_r)), 1)
        log_order = float(log_order)
        log_intercept = float(log_intercept)
        sign_orders = {}
        pair_orders = {}
        for sign, sign_name in ((-1, "negative"), (1, "positive")):
            branch_v = np.asarray([sign * v for v in V_MAGNITUDES], dtype=float)
            branch_r = np.asarray([residual[float(v)] for v in branch_v], dtype=float)
            sign_orders[sign_name] = float(np.polyfit(np.log(np.abs(branch_v)),
                                                     np.log(np.abs(branch_r)), 1)[0])
            ordered = sorted(zip(np.abs(branch_v), np.abs(branch_r)))
            pair_orders[sign_name] = [
                float(math.log(y2 / y1) / math.log(x2 / x1))
                for (x1, y1), (x2, y2) in zip(ordered[:-1], ordered[1:])
            ]
    else:
        log_order = None
        log_intercept = None
        sign_orders = {}
        pair_orders = {}
    r_loo = []
    for drop in range(len(V_MAGNITUDES)):
        keep = np.asarray([abs(v) != V_MAGNITUDES[drop] for v in SIGNED_V])
        r_loo.append(_fit(signed_v[keep], signed_r[keep], 3))

    def ur(x: float) -> float:
        changes = [abs(_eval(r_coeff, x) - _eval(r_quadratic, x))]
        changes.extend(abs(_eval(r_coeff, x) - _eval(c, x)) for c in r_loo)
        return max(changes)

    r_fit_values = {float(v): _eval(r_coeff, float(v)) for v in SIGNED_V}
    r_monotonic: dict[str, bool] = {}
    for sign, sign_name in ((-1, "negative"), (1, "positive")):
        ordered = [abs(residual[sign * v]) for v in V_MAGNITUDES]
        # With magnitudes ascending, |r| must increase strictly; equivalently
        # it decreases strictly as |v| is halved toward zero.
        r_monotonic[sign_name] = all(ordered[i] > ordered[i - 1] for i in range(1, len(ordered)))
    plateau_limit = 10.0 * floor
    plateau_intercept = float(r_coeff[0])
    plateau_pass = abs(plateau_intercept) <= plateau_limit
    e_fit_residuals = e_values - np.asarray([_eval(e_coeff, v*v) for v in magnitudes])
    r_fit_residuals = signed_r - np.asarray([_eval(r_coeff, float(v)) for v in signed_v])
    slope_pass = abs(float(e_coeff[0])) < SLOPE_INTERCEPT_TOL
    order_pass = log_order is not None and log_order >= 2.0
    r_shape_pass = all(r_monotonic.values()) and order_pass and plateau_pass
    unassessed = not log_defined
    return {
        "q0": q0,
        "delta_sumev_ry_by_signed_v": {str(v): deltas[v] for v in SIGNED_V},
        "r_ry_by_signed_v": {str(v): residual[v] for v in SIGNED_V},
        "symmetric_slope_ry_by_abs_v": {str(v): slopes[v] for v in V_MAGNITUDES},
        "e_slope_error_by_abs_v": {str(v): e[v] for v in V_MAGNITUDES},
        "e_even_fit_coefficients_ascending": [float(c) for c in e_coeff],
        "e_linear_sensitivity_intercept": float(e_linear[0]),
        "e_leave_one_magnitude_out_intercepts": [float(x) for x in e_loo],
        "e_intercept_sensitivity_report_only": float(u_e_report),
        "e_fit_residuals": [float(x) for x in e_fit_residuals],
        "slope_extrapolation_pass": slope_pass,
        "e_error_strictly_decreases_as_v_halves": e_monotonic,
        "r_cubic_coefficients_ascending": [float(c) for c in r_coeff],
        "r_quadratic_coefficients_ascending": [float(c) for c in r_quadratic],
        "r_fit_residuals": [float(x) for x in r_fit_residuals],
        "r_fit_sensitivity_diagnostic_by_signed_v": {
            str(v): ur(float(v)) for v in SIGNED_V
        },
        "r_fitted_at_measured_points_ry": {str(v): r_fit_values[float(v)] for v in SIGNED_V},
        "wrong_delta_sign_negative_control_max_abs_residual_ry": max(
            abs(-deltas[v] - v*q0) for v in SIGNED_V
        ),
        "r_loglog_ols_order_all_eight": log_order,
        "r_loglog_intercept_all_eight": log_intercept,
        "r_loglog_orders_by_sign_diagnostic": sign_orders,
        "r_pairwise_log_orders_by_sign_diagnostic": pair_orders,
        "r_abs_strictly_decreases_as_v_halves_by_sign": r_monotonic,
        "plateau_intercept_ry": plateau_intercept,
        "plateau_limit_ry": plateau_limit,
        "observed_smallest_shift_r_ry": {
            "negative": residual[-0.0025], "positive": residual[0.0025]
        },
        "plateau_intercept_pass": plateau_pass,
        "r_log_fit_defined": log_defined,
        "status": "UNASSESSED_LOG_ZERO" if unassessed else
                  ("PASS" if slope_pass and e_monotonic and r_shape_pass else "FAIL_A1_1"),
    }


def analyse_rows(manifest: dict[str, Any], rows: list[dict[str, str]], floor_json: dict[str, Any]) -> dict[str, Any]:
    if str(manifest.get("backend", "")).upper() != "KS":
        raise ValueError("A1.1 is registered for KS only")
    if int(manifest.get("nsp", -1)) != 1 or int(manifest.get("lmax", -1)) != 2:
        raise ValueError("A1.1 requires scalar collinear spd (nsp=1,lmax=2)")
    if manifest.get("metric") != "identity":
        raise ValueError("A1.1 accepted-state oracle requires the captured identity metric")
    if str(floor_json.get("backend", "")).upper() != "KS":
        raise ValueError("floor JSON backend must be KS")
    floor = float(floor_json["max_sum_ev_floor_ry"])
    if not (math.isfinite(floor) and floor > 0.0):
        raise ValueError("KS no-U floor must be finite and positive")
    iteration = int(manifest["iteration"])
    case_rows = [r for r in rows if r["kind"] == "m0_a11_case" and m0._i(r, "iteration") == iteration]
    if not case_rows:
        raise ValueError("capture has no A1 case rows")
    case_keys = {(m0._i(r, "site"), m0._i(r, "type"), m0._i(r, "spin")) for r in case_rows}
    site_types = sorted({(site, type_id) for site, type_id, _spin in case_keys})
    for site, type_id in site_types:
        if {(site, type_id, spin) for spin in (-1, 1, 2)} - case_keys:
            raise ValueError(f"A1 cases lack required spin/aggregate rows at site/type {(site, type_id)}")
    expected_stages = {RAW_STAGE, *STAGE_BY_V.values()}
    actual_stages = {r["stage"] for r in case_rows}
    if actual_stages != expected_stages:
        raise ValueError(f"A1 stage set mismatch: expected {sorted(expected_stages)}, got {sorted(actual_stages)}")
    expected_by_v = {0.0: RAW_STAGE, **{v: stage for v, stage in STAGE_BY_V.items()}}
    grouped: dict[tuple[int, int, int], dict[float, dict[str, float]]] = defaultdict(dict)
    for r in case_rows:
        key = (m0._i(r, "site"), m0._i(r, "type"), m0._i(r, "spin"))
        v = m0._f(r, "weight")
        expected_stage = expected_by_v.get(v)
        if expected_stage is None or r["stage"] != expected_stage:
            raise ValueError(f"A1 signed v/stage mismatch: v={v}, stage={r['stage']}")
        if v in grouped[key]:
            raise ValueError(f"duplicate A1 case at {key}, v={v}")
        grouped[key][v] = {
            "raw_sumev": m0._f(r, "energy"),
            "case_sumev": m0._f(r, "scalar"),
            "q0": m0._f(r, "re"),
            "q2": m0._f(r, "im"),
        }
    for key, series in grouped.items():
        if set(series) != set(expected_by_v):
            raise ValueError(f"A1 case {key} lacks one or more signed shifts")
        raw = series[0.0]
        grouped[key][0.0]["delta_sumev_ry"] = 0.0
        grouped[key][0.0]["r_ry"] = 0.0
        for v in SIGNED_V:
            row = series[v]
            if row["raw_sumev"] != raw["raw_sumev"]:
                raise ValueError(f"raw SUMEV changes across shifts at {key}")
            if row["q0"] != raw["q0"] or row["q2"] != raw["q2"]:
                raise ValueError(f"A1 q0/q2 are not exactly frozen at {key}, v={v}")
            delta = raw["raw_sumev"] - row["case_sumev"]
            grouped[key][v]["delta_sumev_ry"] = delta
            grouped[key][v]["r_ry"] = delta - v * raw["q0"]

    # Verify center and returned-PL rows against the fixed signed convention.
    center_rows = [r for r in rows if r["kind"] == "m0_a11_center" and m0._i(r, "iteration") == iteration]
    center_map = {}
    for r in center_rows:
        key = (m0._i(r, "site"), m0._i(r, "type"), m0._i(r, "spin"), m0._f(r, "weight"))
        if key in center_map:
            raise ValueError(f"duplicate m0_a11_center row {key}")
        center_map[key] = r
    for key, series in grouped.items():
        site, type_id, spin = key
        # Aggregate spin=-1 rows exist for SUMEV/q0/q2 accounting, but the
        # production replay records calculate_pl inputs/outputs per physical
        # spin channel only.  Aggregate cases are independently checked below
        # against the sum of the two spin rows.
        if spin not in (1, 2):
            continue
        for v in (0.0, *SIGNED_V):
            row = center_map.get((site, type_id, spin, v))
            if row is None:
                raise ValueError(f"missing/duplicate m0_a11_center at {key}, v={v}")
            if abs((m0._f(row, "scalar") - m0._f(row, "energy")) + v) > STATE_VALUE_TOL:
                raise ValueError(f"center shift is not exactly -v at {key}, v={v}")
            if abs(m0._f(row, "re") - m0._f(row, "im")) > 0.0 and v == 0.0:
                raise ValueError(f"raw PL differs from returned PL at zero shift for {key}")

    stages = [(0.0, RAW_STAGE), *[(v, STAGE_BY_V[v]) for v in SIGNED_V]]
    # Bind the compact q0/q2 carrier fields to the actual QL array inputs and
    # the captured calculate_pl input/return records.
    for (site, type_id, spin), series in grouped.items():
        if spin not in (1, 2):
            continue
        for v, stage in stages:
            ql0 = _one(rows, kind="m0_replay_ql_input", stage=stage, iteration=iteration,
                       site=site, type_id=type_id, ell=2, spin=spin, order=0)
            ql2 = _one(rows, kind="m0_replay_ql_input", stage=stage, iteration=iteration,
                       site=site, type_id=type_id, ell=2, spin=spin, order=2)
            if m0._f(ql0, "scalar") != series[v]["q0"] or m0._f(ql2, "scalar") != series[v]["q2"]:
                raise ValueError(f"A1 q0/q2 carrier values do not match QL inputs at {(site,type_id,spin,v)}")
            pl_input = _one(rows, kind="m0_replay_pl_input", stage=stage, iteration=iteration,
                            site=site, type_id=type_id, ell=2, spin=spin, order=0)
            center = center_map[(site, type_id, spin, v)]
            if abs(m0._f(pl_input, "scalar") - m0._f(center, "scalar")) > STATE_VALUE_TOL:
                raise ValueError(f"A1 recorded center differs from calculate_pl center input at {(site,type_id,spin,v)}")
            pl_return = _one(rows, kind="m0_replay_pl_return", stage=stage, iteration=iteration,
                             site=site, type_id=type_id, ell=2, spin=spin, order=0)
            expected_pl = m0._f(center, "re") if v == 0.0 else m0._f(center, "im")
            if abs(m0._f(pl_return, "scalar") - expected_pl) > STATE_VALUE_TOL:
                raise ValueError(f"A1 recorded PL return differs from calculate_pl output at {(site,type_id,spin,v)}")

    # Input rows are serialized for every clone.  All maps must have identical
    # keys/values except the d-shell centre (PL-input order zero).
    frozen_specs = (
        ("m0_replay_rofi", ("site", "type", "l", "spin", "row", "col", "order")),
        ("m0_replay_v_input", ("site", "type", "l", "spin", "row", "col", "order")),
        ("m0_replay_rho_input", ("site", "type", "l", "spin", "row", "col", "order")),
        ("m0_replay_ec_seed", ("site", "type", "l", "spin", "row", "col", "order")),
        ("m0_replay_ev_seed", ("site", "type", "l", "spin", "row", "col", "order")),
        ("m0_replay_ql_input", ("site", "type", "l", "spin", "row", "col", "order")),
    )
    frozen_diffs = {}
    for kind, key_fields in frozen_specs:
        base = _row_map(rows, kind, RAW_STAGE, iteration, key_fields)
        max_diff = 0.0
        for v, stage in stages[1:]:
            other = _row_map(rows, kind, stage, iteration, key_fields)
            max_diff = max(max_diff, _assert_same_map(f"frozen {kind} v={v}", base, other, exact=True))
        frozen_diffs[kind] = max_diff
    pl_key_fields = ("site", "type", "l", "spin", "row", "col", "order")
    pl_base = _row_map(rows, "m0_replay_pl_input", RAW_STAGE, iteration, pl_key_fields)
    pl_diff = 0.0
    for v, stage in stages[1:]:
        pl_other = _row_map(rows, "m0_replay_pl_input", stage, iteration, pl_key_fields)
        if pl_base.keys() != pl_other.keys():
            raise ValueError(f"calculate_pl input keys differ at v={v}")
        for key, value in pl_base.items():
            expected = value - v if key[2] == 2 and key[3] in (1, 2) and key[-1] == 0 else value
            diff = abs(pl_other[key] - expected)
            pl_diff = max(pl_diff, diff)
            if diff > STATE_VALUE_TOL:
                raise ValueError(f"calculate_pl input changed unexpectedly at v={v}, key={key}, diff={diff:.3e}")

    # Cross-check aggregate rows against the two spin channels to prevent a
    # cancellation or spin-index misbinding from passing per-channel fits.
    for site, type_id in site_types:
        for v in (0.0, *SIGNED_V):
            up, down, total = (grouped[(site, type_id, s)][v] for s in (1, 2, -1))
            for field in ("raw_sumev", "case_sumev", "q0", "q2"):
                if abs(total[field] - up[field] - down[field]) > STATE_VALUE_TOL:
                    raise ValueError(f"spin aggregate mismatch {field} at site/type/v {(site, type_id, v)}")

    reports = {}
    aggregate_reports = {}
    for key, series in grouped.items():
        if key[2] == -1:
            # The summed channel is a required independent accounting check;
            # it is reported separately and does not replace the two physical
            # spin-channel A1 gates.
            q0_total = series[0.0]["q0"]
            if q0_total == 0.0:
                raise ValueError(f"zero aggregate q0 prevents normalized slope at {key}")
            aggregate_reports["site=%d,type=%d,spin=%d" % key] = _series_report(series, q0_total, floor)
            continue
        if key[2] not in (1, 2):
            continue
        q0 = series[0.0]["q0"]
        if q0 == 0.0:
            raise ValueError(f"zero q0 prevents normalized slope at {key}")
        reports["site=%d,type=%d,spin=%d" % key] = _series_report(series, q0, floor)

    all_statuses = [x["status"] for x in reports.values()]
    overall_status = "PASS" if all_statuses and all(s == "PASS" for s in all_statuses) else "FAIL_A1_1"
    if "UNASSESSED_LOG_ZERO" in all_statuses:
        overall_status = "UNASSESSED_LOG_ZERO"
    return {
        "schema": "rslmto-ldau-a1-1-scaling-v1",
        "verdict": overall_status,
        "scope": "A1.1 no-U KS scalar-shift scaling only; no M0.3/U-state or energy-closure claim",
        "capture_run_id": manifest["run_id"],
        "iteration": iteration,
        "capture_csv_sha256": manifest["csv_sha256"],
        "capture_executable_sha256": manifest.get("executable_sha256", manifest.get("binary_sha256")),
        "capture_source_sha": manifest.get("source_sha", manifest.get("self_f90_sha256")),
        "input_sha256": manifest.get("input_sha256"),
        "accepted_floor_ry": floor,
        "fit_registration": {
            "signed_v_ry": list(SIGNED_V),
            "delta_definition": "SUMEV_RAW - SUMEV_case",
            "r_definition": "delta - v*q0",
            "slope_definition": "[delta(+v)-delta(-v)]/(2v)",
            "symmetric_slope_fit": "unweighted a0+a2*v^2+a4*v^4 over four positive magnitudes; gate abs(a0)<1e-6",
            "r_fit": "unweighted cubic in signed v over all eight points; unweighted signed quadratic sensitivity",
            "primary_log_fit": "OLS log(abs(r)) versus log(abs(v)) across all eight signed points, order >=2",
            "per_sign_and_pairwise_log_orders": "diagnostics only",
            "plateau": "report cubic constant and observed r values against 10*D0; gate constant magnitude <=10*D0",
            "fit_residuals_and_sensitivity": "reported, not additional gates",
        },
        "accepted_state_comparison": None,
        "frozen_replay_input_checks": {"verdict": "PASS", "max_differences": frozen_diffs,
                                        "max_calculate_pl_input_difference_after_expected_shift": pl_diff},
        "per_site_spin": reports,
        "aggregate_spin_diagnostics": aggregate_reports,
        "floor_definition": floor_json.get("definition"),
    }


def _self_test() -> None:
    """Exercise signed conventions and fit gates on a known analytic series."""
    q0 = 3.25
    values = {}
    for v in (0.0, *SIGNED_V):
        # The tiny odd term makes the symmetric-slope error nonzero and
        # strictly scale down, while keeping the all-sign log order at 2.
        r = 0.4*v*v + 1.0e-7*v**3
        delta = q0*v + r
        values[v] = {"delta_sumev_ry": delta, "r_ry": delta - q0*v}
    report = _series_report(values, q0, floor=1.0e-5)
    if report["status"] != "PASS":
        raise AssertionError(f"analytic A1 series should pass: {report['status']}")
    if abs(report["r_loglog_ols_order_all_eight"] - 2.0) > 0.02:
        raise AssertionError("analytic quadratic residual should have log order near two")
    if report["wrong_delta_sign_negative_control_max_abs_residual_ry"] < 0.01:
        raise AssertionError("wrong-sign control did not detect the sign convention")

    # A genuine v-independent residual above 10*D0 must be caught as a
    # plateau; this checks the measured-shape and fitted-constant gates.
    plateau = {}
    for v in (0.0, *SIGNED_V):
        r = 2.0e-4
        delta = q0*v + r
        plateau[v] = {"delta_sumev_ry": delta, "r_ry": delta - q0*v}
    plateau_report = _series_report(plateau, q0, floor=1.0e-5)
    if plateau_report["status"] != "FAIL_A1_1" or plateau_report["plateau_intercept_pass"]:
        raise AssertionError("above-10*D0 flat plateau was not rejected")

    # A smooth first-order residual must fail the mandated log-order gate.
    linear = {}
    for v in (0.0, *SIGNED_V):
        r = 0.01*v
        delta = q0*v + r
        linear[v] = {"delta_sumev_ry": delta, "r_ry": delta - q0*v}
    linear_report = _series_report(linear, q0, floor=1.0e-5)
    if linear_report["r_loglog_ols_order_all_eight"] is None or linear_report["r_loglog_ols_order_all_eight"] >= 2.0:
        raise AssertionError("linear residual did not fail the log-order gate")
    if linear_report["status"] != "FAIL_A1_1":
        raise AssertionError("linear residual was not stopped for A1.3 review")

    # Exact zero is intentionally not regularized with a pseudocount.
    zero = {v: {"delta_sumev_ry": q0*v, "r_ry": 0.0} for v in (0.0, *SIGNED_V)}
    zero_report = _series_report(zero, q0, floor=1.0e-5)
    if zero_report["status"] != "UNASSESSED_LOG_ZERO":
        raise AssertionError("exact-zero log residual must be marked unassessed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture_dir", type=Path, nargs="?",
                        help="A1.1 trace directory containing manifest.json/records.csv")
    parser.add_argument("--reference-dir", type=Path,
                        help="accepted iteration-23 M0.2 KS trace directory")
    parser.add_argument("--floor-json", type=Path,
                        help="hash-bound KS M0.2 floor JSON")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true", help="run synthetic math/sign checks only")
    args = parser.parse_args()

    if args.self_test:
        _self_test()
        print("A1.1 signed fit self-test PASS")
        return
    if args.capture_dir is None or args.reference_dir is None or args.floor_json is None or args.output is None:
        parser.error("capture_dir, --reference-dir, --floor-json, and --output are required")

    manifest, rows = m0.load_capture(args.capture_dir)
    ref_manifest, ref_rows = m0.load_capture(args.reference_dir)
    accepted_state = _check_reference_state(
        args.capture_dir, manifest, rows, args.reference_dir, ref_manifest, ref_rows)
    result = analyse_rows(manifest, rows, json.loads(args.floor_json.read_text(encoding="utf-8")))
    result["accepted_state_comparison"] = accepted_state
    result["reference_capture_csv_sha256"] = ref_manifest["csv_sha256"]
    result["reference_capture_source_sha"] = ref_manifest.get("source_sha")
    result["floor_json_sha256"] = _sha(args.floor_json)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": result["verdict"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
