#!/usr/bin/env python3
"""Independent A2.2 extension of the preregistered A1.1 KS shift analysis.

The original A1 analyzer and its eight-shift outputs remain unchanged. This
wrapper reuses its capture/state checks and numerical fit implementations with
the six preregistered A2 shifts added. The only new shape check is constant
residual sign within each signed-v branch; exact zero remains unassessed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import m0_a11_scaling_oracle as a1
import local_moment_capture_oracle as m0


V_MAGNITUDES = (0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.12)
SIGNED_V = (-0.12, -0.08, -0.04, -0.02, -0.01, -0.005, -0.0025,
            0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.12)
STAGE_BY_V = {
    -0.12: "A2_SHIFT_11", -0.08: "A2_SHIFT_10", -0.04: "A2_SHIFT_09",
    -0.02: "A1_SHIFT_01", -0.01: "A1_SHIFT_02", -0.005: "A1_SHIFT_03",
    -0.0025: "A1_SHIFT_04", 0.0025: "A1_SHIFT_05", 0.005: "A1_SHIFT_06",
    0.01: "A1_SHIFT_07", 0.02: "A1_SHIFT_08", 0.04: "A2_SHIFT_12",
    0.08: "A2_SHIFT_13", 0.12: "A2_SHIFT_14",
}
A1_SIGNED_V = (-0.02, -0.01, -0.005, -0.0025, 0.0025, 0.005, 0.01, 0.02)


def _bind_preregistered_scan() -> None:
    """Bind A1's unchanged fit/check functions to the A2.2 point set."""
    a1.V_MAGNITUDES = V_MAGNITUDES
    a1.SIGNED_V = SIGNED_V
    a1.STAGE_BY_V = STAGE_BY_V


def _branch_signs(residuals: dict[str, float]) -> dict[str, dict[str, Any]]:
    report: dict[str, dict[str, Any]] = {}
    for sign, name in ((-1, "negative_v"), (1, "positive_v")):
        values = [float(residuals[str(sign * mag)]) for mag in V_MAGNITUDES]
        signs = [0 if x == 0.0 else (1 if x > 0.0 else -1) for x in values]
        report[name] = {
            "r_by_increasing_abs_v_ry": values,
            "sign_by_increasing_abs_v": signs,
            "constant_nonzero_sign_within_branch": bool(signs and signs[0] != 0
                                                           and all(x == signs[0] for x in signs)),
            "zero_is_unassessed_without_pseudocount": any(x == 0 for x in signs),
        }
    return report


def a21_intercept_bound(a1_result: dict[str, Any], t_by_site_spin: dict[str, float],
                        d0_floor_ry: float) -> dict[str, Any]:
    """Propagate the immutable A1 a0 fits to the accepted U-state t values."""
    per_spin = a1_result.get("per_site_spin", {})
    if not per_spin or set(per_spin) != set(t_by_site_spin):
        raise ValueError("A2.1 requires one t value for each physical A1 site/type/spin fit")
    terms = {}
    for key, t_value in t_by_site_spin.items():
        if ",spin=-1" in key:
            raise ValueError("A2.1 uses physical spin rows, not aggregate diagnostics")
        a0 = float(per_spin[key]["e_even_fit_coefficients_ascending"][0])
        terms[key] = {"a0": a0, "t_s_ry": float(t_value),
                      "absolute_error_bound_ry": abs(a0) * abs(float(t_value))}
    bound = sum(v["absolute_error_bound_ry"] for v in terms.values())
    return {"terms": terms, "sum_abs_a0_t_ry": bound, "D0_ry": float(d0_floor_ry),
            "pass": bound <= float(d0_floor_ry)}


def a22_m03_allowance(a22_result: dict[str, Any], q0_by_site_spin: dict[str, float],
                      t_by_site_spin: dict[str, float], d0_floor_ry: float) -> dict[str, Any]:
    """Evaluate the preregistered cubic at actual shifts and add D0 once."""
    per_spin = a22_result.get("per_site_spin", {})
    if not per_spin or set(per_spin) != set(q0_by_site_spin) or set(per_spin) != set(t_by_site_spin):
        raise ValueError("A2.2 allowance needs matching q0/t entries for every physical spin fit")
    terms = {}
    for key in per_spin:
        q0 = float(q0_by_site_spin[key])
        t_value = float(t_by_site_spin[key])
        if q0 <= 0.0:
            raise ValueError(f"A2.2 requires positive q0 for {key}")
        shift = t_value / q0
        if abs(shift) > 0.12:
            raise ValueError(f"A2.2 shift for {key} is outside scanned [-0.12,0.12] Ry: {shift:.9g}")
        coeff = per_spin[key]["r_cubic_coefficients_ascending"]
        fitted = sum(float(c) * shift**i for i, c in enumerate(coeff))
        terms[key] = {"q0_s": q0, "t_s_ry": t_value, "shift_ry": shift,
                      "r_fit_ry": fitted, "absolute_r_fit_ry": abs(fitted)}
    model_error = sum(v["absolute_r_fit_ry"] for v in terms.values())
    return {"terms": terms, "sum_abs_r_fit_ry": model_error, "D0_ry": float(d0_floor_ry),
            "allowance_ry": model_error + float(d0_floor_ry),
            "floor_added_once": True}


def analyse_rows(manifest: dict[str, Any], rows: list[dict[str, str]],
                 floor_json: dict[str, Any]) -> dict[str, Any]:
    _bind_preregistered_scan()
    result = a1.analyse_rows(manifest, rows, floor_json)
    if result["capture_csv_sha256"] != manifest.get("csv_sha256"):
        raise ValueError("A2.2 result CSV hash does not match its manifest")

    statuses = []
    for section_name, section in (("per_site_spin", result["per_site_spin"]),
                                  ("aggregate_spin_diagnostics", result["aggregate_spin_diagnostics"])):
        for item in section.values():
            signs = _branch_signs(item["r_ry_by_signed_v"])
            sign_ok = all(v["constant_nonzero_sign_within_branch"] for v in signs.values())
            exact_zero = any(v["zero_is_unassessed_without_pseudocount"] for v in signs.values())
            order_ok = item["r_loglog_ols_order_all_eight"] is not None and \
                item["r_loglog_ols_order_all_eight"] >= 2.0
            plateau_ok = item["plateau_intercept_pass"]
            monotonic_ok = all(item["r_abs_strictly_decreases_as_v_halves_by_sign"].values())
            shape_ok = order_ok and plateau_ok and monotonic_ok and sign_ok
            item["a2_branch_sign_diagnostics"] = signs
            item["a2_shape_gates"] = {
                "all_signed_log_order_at_least_two": order_ok,
                "cubic_constant_within_existing_10D0_plateau_limit": plateau_ok,
                "absolute_residual_strictly_decreases_toward_zero_on_both_branches": monotonic_ok,
                "no_sign_reversal_within_either_signed_v_branch": sign_ok,
                "exact_zero_log_residual_unassessed": exact_zero,
            }
            # A2 supersedes A1.1(i); keep the even-quartic slope fit as a
            # diagnostic, but do not retain the A1 slope-intercept gate.
            item["slope_extrapolation_pass"] = "SUPERSEDED_BY_A2_1_DECISION_TEST"
            item["r_loglog_ols_order_all_fourteen"] = item.pop("r_loglog_ols_order_all_eight")
            item["r_loglog_intercept_all_fourteen"] = item.pop("r_loglog_intercept_all_eight")
            item["status"] = "UNASSESSED_LOG_ZERO" if exact_zero else ("PASS" if shape_ok else "FAIL_A2_2")
            if section_name == "per_site_spin":
                statuses.append(item["status"])

    result["schema"] = "rslmto-ldau-a2-2-scaling-v1"
    result["scope"] = "A2.2 extended no-U KS scalar-shift scan; no M0.3 or energy-closure claim"
    result["fit_registration"] = {
        "signed_v_ry_in_fit_order": list(SIGNED_V),
        "A1_points_preserved": list(A1_SIGNED_V),
        "A2_stage_to_v_ry": {stage: v for v, stage in STAGE_BY_V.items()},
        "delta_definition": "SUMEV_RAW - SUMEV_case",
        "r_definition": "delta - v*q0",
        "slope_definition": "[delta(+v)-delta(-v)]/(2v)",
        "slope_model": "unchanged A1 unweighted even quartic a0+a2*v^2+a4*v^4 on seven positive magnitudes; a0 is diagnostic under A2.1",
        "r_model": "unchanged A1 unweighted signed cubic for fit and signed quadratic sensitivity",
        "primary_log_fit": "unchanged A1 OLS log(abs(r)) versus log(abs(v)) across all fourteen signed points; order >=2",
        "shape_checks": "strict |r| monotonicity toward zero on each sign branch; no sign reversal within either branch; cubic constant <= 10*D0",
        "zero_rule": "exact zero makes the log fit undefined; no pseudocount; report unassessed",
        "shape_error_model": "D0 is the additive radial inversion offset in the M0.3 allowance, not a per-point noise band; A2 shape uses deterministic same-state differences and existing A1 tests without a new cutoff",
    }
    result["verdict"] = "UNASSESSED_LOG_ZERO" if "UNASSESSED_LOG_ZERO" in statuses else \
        ("PASS" if statuses and all(x == "PASS" for x in statuses) else "FAIL_A2_2")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture_dir", type=Path, nargs="?", help="A2.2 trace dir with manifest.json and records.csv")
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--floor-json", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true", help="test A2 fit and branch-shape logic only")
    args = parser.parse_args()
    _bind_preregistered_scan()
    if args.self_test:
        q0 = 3.25
        series = {}
        for v in SIGNED_V:
            r = 0.4 * v * v + 1.0e-7 * v**3
            series[v] = {"delta_sumev_ry": q0 * v + r, "r_ry": r}
        report = a1._series_report(series, q0, 1.0e-5)
        signs = _branch_signs(report["r_ry_by_signed_v"])
        if report["r_log_fit_defined"] is not True:
            raise AssertionError("nonzero analytic A2 residual should define the log fit")
        if report["r_loglog_ols_order_all_eight"] is None or \
                report["r_loglog_ols_order_all_eight"] < 2.0:
            raise AssertionError("analytic quadratic residual should retain order >=2")
        if not all(item["constant_nonzero_sign_within_branch"] for item in signs.values()):
            raise AssertionError("analytic residual should have no within-branch sign reversal")
        if not all(report["r_abs_strictly_decreases_as_v_halves_by_sign"].values()):
            raise AssertionError("analytic residual should be monotone on each branch")
        broken = dict(report["r_ry_by_signed_v"])
        broken["-0.04"] = -abs(float(broken["-0.04"]))
        if _branch_signs(broken)["negative_v"]["constant_nonzero_sign_within_branch"]:
            raise AssertionError("within-branch sign reversal was not detected")
        zero = dict(report["r_ry_by_signed_v"])
        zero["0.04"] = 0.0
        if not _branch_signs(zero)["positive_v"]["zero_is_unassessed_without_pseudocount"]:
            raise AssertionError("exact zero should make its branch unassessed")
        print("A2.2 fit and branch-sign self-test PASS")
        return
    if args.capture_dir is None or args.reference_dir is None or args.floor_json is None or args.output is None:
        parser.error("capture_dir, --reference-dir, --floor-json, and --output are required")
    manifest, rows = m0.load_capture(args.capture_dir)
    ref_manifest, ref_rows = m0.load_capture(args.reference_dir)
    accepted = a1._check_reference_state(args.capture_dir, manifest, rows,
                                         args.reference_dir, ref_manifest, ref_rows)
    floor_json = json.loads(args.floor_json.read_text(encoding="utf-8"))
    result = analyse_rows(manifest, rows, floor_json)
    result["accepted_state_comparison"] = accepted
    result["reference_capture_csv_sha256"] = ref_manifest["csv_sha256"]
    result["reference_capture_source_sha"] = ref_manifest.get("source_sha")
    result["floor_json_sha256"] = hashlib.sha256(args.floor_json.read_bytes()).hexdigest()
    result["a1_analyzer_sha256"] = hashlib.sha256(Path(a1.__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": result["verdict"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
