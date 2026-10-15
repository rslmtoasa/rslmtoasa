#!/usr/bin/env python3
"""Apply the preregistered A1/A2 fits to original A1 rows plus private replay A2 rows.

This adapter does not synthesize a capture CSV. It extracts the immutable eight
A1 measurements from the accepted historical CSV and appends the six measured
SUMEV values from the hash-bound accepted-state replay JSON, then calls the
existing A1 fit implementation with the preregistered A2 point globals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m0_a11_scaling_oracle as a1  # noqa: E402
import m0_a22_scaling_oracle as a22  # noqa: E402

CSV_SHA = "55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4"
A1_RESULT_SHA = "48b0cb79291b12fd5e3c027b92379f0f70b45e8418882348616a1ea1abed27ae"
BASELINE_GATE_SHA = "549edf2183ce13fbec7556f2bcde3d23b4f1c78eb908c04e616641a4d694d179"
EXTRA_SHA = "93236aeaa5c85afcbf8c0075104f78462888c45cdbcdb2a70b51fa066c86fc79"
FIXTURE_SHA = "e88fcff8f1a83831b13834324a6730fa6f8433d7ca3f3b389119d6736b1adbb3"
DRIVER_SHA = "496b550a30f372dc6e623627fb46c2ef5e1e8854228df58b91a4f2255b72a008"
SOURCE_SHA = "577749f1dc11b0f1ccbce559e83d91925b56abf375823ab1cb4bd95705d383a9"
U2_CSV_SHA = "dfffa07b366ae47076b8b7aba3ee610e41e65af5e8bc28bb82c0425a62588600"
U2_MANIFEST_SHA = "5d1a052fce4242b0cbaae8a79b920872a5102d4b037e88ee9b34ff2cfff861ae"
A1_STAGES = {
    -0.02: "A1_SHIFT_01", -0.01: "A1_SHIFT_02", -0.005: "A1_SHIFT_03",
    -0.0025: "A1_SHIFT_04", 0.0025: "A1_SHIFT_05", 0.005: "A1_SHIFT_06",
    0.01: "A1_SHIFT_07", 0.02: "A1_SHIFT_08",
}
A2_STAGES = {
    -0.04: "A2_SHIFT_09", -0.08: "A2_SHIFT_10", -0.12: "A2_SHIFT_11",
    0.04: "A2_SHIFT_12", 0.08: "A2_SHIFT_13", 0.12: "A2_SHIFT_14",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_cases(capture_csv: Path) -> dict[tuple[str, int], dict[str, float]]:
    if sha(capture_csv) != CSV_SHA:
        raise ValueError("historical accepted CSV hash mismatch")
    found: dict[tuple[str, int], dict[str, float]] = {}
    with capture_csv.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["kind"] != "m0_a11_case" or int(row["iteration"]) != 23:
                continue
            spin = int(row["spin"])
            if spin not in (1, 2, -1):
                continue
            key = (row["stage"], spin)
            if key in found:
                raise ValueError(f"duplicate historical A1 case {key}")
            found[key] = {name: float(row[col]) for name, col in
                          (("raw_sumev", "energy"), ("case_sumev", "scalar"),
                           ("q0", "re"), ("q2", "im"), ("v", "weight"))}
    return found


def verify_spin_aggregates(cases: dict[tuple[str, int], dict[str, float]]) -> None:
    """Repeat A1's per-spin versus aggregate accounting check on original rows."""
    tol = a1.STATE_VALUE_TOL
    for stage in ("RAW_PL", *A1_STAGES.values()):
        for field in ("raw_sumev", "case_sumev", "q0", "q2"):
            total = cases[(stage, -1)][field]
            up = cases[(stage, 1)][field]
            down = cases[(stage, 2)][field]
            if abs(total - up - down) > tol:
                raise ValueError(f"historical spin aggregate mismatch {field} at {stage}")


def make_values(cases: dict[tuple[str, int], dict[str, float]], extra: dict[str, Any],
                spin: int) -> tuple[dict[float, dict[str, float]], float]:
    raw = cases[("RAW_PL", spin)]
    q0 = raw["q0"]
    q2 = raw["q2"]
    if q0 <= 0.0:
        raise ValueError(f"nonpositive accepted q0 for spin {spin}")
    values: dict[float, dict[str, float]] = {0.0: {"delta_sumev_ry": 0.0, "r_ry": 0.0}}
    for v, stage in A1_STAGES.items():
        row = cases[(stage, spin)]
        if row["v"] != v or row["raw_sumev"] != raw["raw_sumev"] or row["q0"] != q0 or row["q2"] != q2:
            raise ValueError(f"historical A1 operands changed at {stage}, spin {spin}")
        delta = raw["raw_sumev"] - row["case_sumev"]
        values[v] = {"delta_sumev_ry": delta, "r_ry": delta - v*q0}
    for v, stage in A2_STAGES.items():
        point = extra["cases"][stage]
        if float(point["v"]) != v:
            raise ValueError(f"A2 shift mapping mismatch at {stage}")
        sumev = float(point["sumev"][f"spin_{spin}"])
        delta = raw["raw_sumev"] - sumev
        values[v] = {"delta_sumev_ry": delta, "r_ry": delta - v*q0}
    expected_v = {0.0, *a22.SIGNED_V}
    if set(values) != expected_v:
        raise ValueError("combined A1/A2 series does not contain exactly fourteen shifts plus zero")
    return values, q0


def fit(values: dict[float, dict[str, float]], q0: float, d0: float) -> dict[str, Any]:
    a22._bind_preregistered_scan()
    report = a1._series_report(values, q0, d0)
    branches = a22._branch_signs(report["r_ry_by_signed_v"])
    sign_ok = all(x["constant_nonzero_sign_within_branch"] for x in branches.values())
    exact_zero = any(x["zero_is_unassessed_without_pseudocount"] for x in branches.values())
    order_ok = report["r_loglog_ols_order_all_eight"] is not None and report["r_loglog_ols_order_all_eight"] >= 2.0
    plateau_ok = report["plateau_intercept_pass"]
    monotonic_ok = all(report["r_abs_strictly_decreases_as_v_halves_by_sign"].values())
    shape_ok = order_ok and plateau_ok and monotonic_ok and sign_ok
    report["a2_branch_sign_diagnostics"] = branches
    report["a2_shape_gates"] = {
        "all_signed_log_order_at_least_two": order_ok,
        "cubic_constant_within_existing_10D0_plateau_limit": plateau_ok,
        "absolute_residual_strictly_decreases_toward_zero_on_both_branches": monotonic_ok,
        "no_sign_reversal_within_either_signed_v_branch": sign_ok,
        "exact_zero_log_residual_unassessed": exact_zero,
    }
    report["r_loglog_ols_order_all_fourteen"] = report.pop("r_loglog_ols_order_all_eight")
    report["r_loglog_intercept_all_fourteen"] = report.pop("r_loglog_intercept_all_eight")
    report["slope_extrapolation_pass"] = "SUPERSEDED_BY_A2_1_DECISION_TEST"
    report["status"] = "UNASSESSED_LOG_ZERO" if exact_zero else ("PASS" if shape_ok else "FAIL_A2_2")
    return report


def accepted_u2_targets(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Recompute U-state q0 and Tr(Vin Pout) from hash-bound iteration 32 rows."""
    if sha(path) != U2_CSV_SHA:
        raise ValueError("accepted U2 capture CSV hash mismatch")
    manifest_path = path.with_name("manifest.final.json")
    if sha(manifest_path) != U2_MANIFEST_SHA:
        raise ValueError("accepted U2 final manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("output", {}).get("csv_sha256") != U2_CSV_SHA or
            manifest.get("output", {}).get("captured_iterations") != [32] or
            manifest.get("stage") != "A2.1 eligible state capture" or
            manifest.get("status") != "COMPLETE_CONVERGED"):
        raise ValueError("accepted U2 manifest does not bind the expected KS iteration-32 capture")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if {int(r["iteration"]) for r in rows} != {32}:
        raise ValueError("U2 target capture is not the accepted rolling-final iteration 32")

    def isrow(r: dict[str, str], kind: str, stage: str, spin: int | None = None,
              ell: int | None = None, order: int | None = None) -> bool:
        return (r["kind"] == kind and r["stage"] == stage and int(r["site"]) == 0 and
                (spin is None or int(r["spin"]) == spin) and
                (ell is None or int(r["l"]) == ell) and
                (order is None or int(r["order"]) == order))

    vrows = [r for r in rows if r["kind"] == "v_hub" and r["stage"] == "Pin_H" and
             int(r["iteration"]) == 32 and int(r["site"]) == -1 and int(r["type"]) == 1]
    if len(vrows) != 324:
        raise ValueError(f"expected complete 18x18 Vin matrix for Fe type 1; got {len(vrows)}")
    vin = np.zeros((18, 18), dtype=np.complex128)
    seen_v = set()
    for r in vrows:
        ij = (int(r["row"]), int(r["col"]))
        if ij in seen_v or not all(0 <= x < 18 for x in ij):
            raise ValueError(f"duplicate/out-of-range Vin matrix key {ij}")
        seen_v.add(ij)
        vin[ij] = complex(float(r["re"]), float(r["im"]))

    q0s: dict[str, float] = {}
    ts: dict[str, float] = {}
    for spin in (1, 2):
        prows = [r for r in rows if isrow(r, "ldm", "P_OUT", spin=spin, ell=2) and int(r["type"]) == 1]
        if len(prows) != 25:
            raise ValueError(f"expected complete 5x5 P_OUT d projector spin={spin}; got {len(prows)}")
        p = np.zeros((5, 5), dtype=np.complex128)
        seen_p = set()
        for r in prows:
            ij = (int(r["row"]), int(r["col"]))
            if ij in seen_p or not all(0 <= x < 5 for x in ij):
                raise ValueError(f"duplicate/out-of-range P_OUT key {ij}")
            seen_p.add(ij)
            p[ij] = complex(float(r["re"]), float(r["im"]))
        qrows = [r for r in rows if isrow(r, "radial_moment", "P_OUT", spin=spin, ell=2, order=0)]
        if len(qrows) != 1:
            raise ValueError(f"expected one U-state q0 for spin={spin}; got {len(qrows)}")
        q0 = float(qrows[0]["scalar"])
        block = slice(4 + 9*(spin - 1), 9 + 9*(spin - 1))
        t = float(np.trace(vin[block, block] @ p).real)
        if q0 <= 0.0:
            raise ValueError(f"nonpositive accepted U-state q0 for spin={spin}")
        key = f"site=0,type=0,spin={spin}"
        q0s[key], ts[key] = q0, t
    return q0s, ts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture_csv", type=Path)
    parser.add_argument("a1_result", type=Path)
    parser.add_argument("baseline_gate", type=Path)
    parser.add_argument("extra_six", type=Path)
    parser.add_argument("u2_capture_csv", type=Path)
    parser.add_argument("--d0-ry", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path, expected, label in ((args.a1_result, A1_RESULT_SHA, "historical A1 result"),
                                  (args.baseline_gate, BASELINE_GATE_SHA, "A1 replay gate"),
                                  (args.extra_six, EXTRA_SHA, "six-case replay")):
        if sha(path) != expected:
            raise ValueError(f"{label} SHA256 mismatch")
    a1_result = json.loads(args.a1_result.read_text())
    gate = json.loads(args.baseline_gate.read_text())
    extra = json.loads(args.extra_six.read_text())
    if (a1_result.get("capture_csv_sha256") != CSV_SHA or a1_result.get("verdict") != "FAIL_A1_1" or
            gate.get("status") != "PASS_A1_GATE_STOP_BEFORE_A2" or
            extra.get("status") != "A2_EXTRA_CASES_CAPTURED_AFTER_HASH_BOUND_A1_GATE"):
        raise ValueError("historical/gate/replay status does not match the registered analysis inputs")
    if extra.get("a2_gate", {}).get("report_sha256") != BASELINE_GATE_SHA:
        raise ValueError("six-case output is not bound to the reviewed passing baseline+A1 gate")
    prov = extra.get("provenance", {})
    if (prov.get("capture_csv_sha256") != CSV_SHA or prov.get("fixture_sha256") != FIXTURE_SHA or
            prov.get("replay_driver_sha256") != DRIVER_SHA or prov.get("replay_source_sha256") != SOURCE_SHA or
            prov.get("accepted_iteration") != 23 or prov.get("backend") != "KS"):
        raise ValueError("six-case replay provenance is not the accepted KS iteration-23 fixture")
    stages = extra.get("a2_extra_run", {}).get("stage_mapping", [])
    if stages != [A2_STAGES[v] for v in (-0.04, -0.08, -0.12, 0.04, 0.08, 0.12)]:
        raise ValueError("six-case stage mapping differs from A2.2 preregistration")

    fixture_path = Path(prov["fixture_path"]).resolve()
    if sha(fixture_path) != FIXTURE_SHA:
        raise ValueError("hash-bound accepted-state fixture is missing or changed")
    fixture = json.loads(fixture_path.read_text())
    if fixture.get("schema") != "rslmto-accepted-iter23-private-radial-replay-fixture-v1":
        raise ValueError("unsupported accepted-state fixture schema")

    rows = read_cases(args.capture_csv)
    verify_spin_aggregates(rows)
    d0 = float(args.d0_ry)
    if abs(d0 - 5.028166848219939e-6) > 1e-18:
        raise ValueError("D0 differs from the frozen A1 accepted no-U floor")
    per_spin = {}
    all_series = {}
    combined_by_spin = {}
    for spin in (1, 2):
        values, q0 = make_values(rows, extra, spin)
        fixture_q0 = float(fixture["ql_by_spin_l"][f"spin_{spin}"]["l_2"][0])
        fixture_q2 = float(fixture["ql_by_spin_l"][f"spin_{spin}"]["l_2"][2])
        raw_q2 = rows[("RAW_PL", spin)]["q2"]
        if q0 != fixture_q0 or raw_q2 != fixture_q2:
            raise ValueError(f"fixture QL q0/q2 do not exactly match historical raw operands for spin {spin}")
        # Independently reproduce the accepted eight-point A1 coefficients
        # before binding the same fit code to the fourteen-point A2 set.
        a1_values = {v: values[v] for v in (0.0, *A1_STAGES.keys())}
        a1_report = a1._series_report(a1_values, q0, d0)
        a1_key = f"site=0,type=0,spin={spin}"
        ref_a1 = a1_result["per_site_spin"][a1_key]
        for field in ("e_even_fit_coefficients_ascending", "r_cubic_coefficients_ascending"):
            got, expected = a1_report[field], ref_a1[field]
            if len(got) != len(expected) or any(float(x) != float(y) for x, y in zip(got, expected)):
                raise ValueError(f"A1 fit coefficient reproduction differs for {a1_key}: {field}")
        all_series[f"site=0,type=0,spin={spin}"] = values
        combined_by_spin[spin] = (values, q0)
    # Bind the extended A2 point arrays only after both historical eight-point
    # fits have been checked against the immutable reviewed A1 result.
    for spin in (1, 2):
        values, q0 = combined_by_spin[spin]
        per_spin[f"site=0,type=0,spin={spin}"] = {"q0": q0, **fit(values, q0, d0)}
    # Aggregate diagnostic is the sum of spin-channel deltas and residuals.
    agg_values = {}
    for v in (0.0, *a22.SIGNED_V):
        agg_values[v] = {
            "delta_sumev_ry": sum(all_series[f"site=0,type=0,spin={s}"][v]["delta_sumev_ry"] for s in (1, 2)),
            "r_ry": sum(all_series[f"site=0,type=0,spin={s}"][v]["r_ry"] for s in (1, 2)),
        }
    q0_total = per_spin["site=0,type=0,spin=1"]["q0"] + per_spin["site=0,type=0,spin=2"]["q0"]
    aggregate = {"site=0,type=0,spin=-1": {"q0": q0_total, **fit(agg_values, q0_total, d0)}}
    u2_q0, u2_t = accepted_u2_targets(args.u2_capture_csv)
    a21 = a22.a21_intercept_bound(a1_result, u2_t, d0)
    allowance = a22.a22_m03_allowance({"per_site_spin": per_spin}, u2_q0, u2_t, d0)
    scalar_v = 0.01
    scalar_terms = {}
    for key, item in per_spin.items():
        coeff = item["r_cubic_coefficients_ascending"]
        rfit = sum(float(c) * scalar_v**power for power, c in enumerate(coeff))
        scalar_terms[key] = {"no_u_scan_q0_s": item["q0"], "u_state_q0_s": u2_q0[key],
                             "shift_ry": scalar_v,
                             "r_fit_ry": rfit, "absolute_r_fit_ry": abs(rfit)}
    scalar_model_error = sum(x["absolute_r_fit_ry"] for x in scalar_terms.values())
    scalar_allowance = scalar_model_error + d0
    statuses = [v["status"] for v in per_spin.values()]
    verdict = "UNASSESSED_LOG_ZERO" if "UNASSESSED_LOG_ZERO" in statuses else ("PASS" if statuses and all(x == "PASS" for x in statuses) else "FAIL_A2_2")
    result = {
        "schema": "rslmto-ldau-a2-2-accepted-replay-adapter-v1",
        "verdict": verdict,
        "scope": "A2.2 extension using immutable A1 accepted capture rows and accepted-state private replay; no M0.3 or energy-closure claim",
        "provenance": {
            "historical_a1_capture_csv_sha256": CSV_SHA,
            "historical_a1_result_sha256": sha(args.a1_result),
            "historical_a1_result_verdict": a1_result["verdict"],
            "baseline_a1_gate_sha256": sha(args.baseline_gate),
            "six_case_replay_sha256": sha(args.extra_six),
            "accepted_u2_capture_csv_sha256": sha(args.u2_capture_csv),
            "accepted_u2_manifest_sha256": U2_MANIFEST_SHA,
            "a1_fit_oracle_sha256": sha(Path(a1.__file__)),
            "a2_fit_oracle_sha256": sha(Path(a22.__file__)),
            "adapter_sha256": sha(Path(__file__)),
            "fixture_sha256": prov["fixture_sha256"], "replay_driver_sha256": prov["replay_driver_sha256"],
            "replay_source_sha256": prov["replay_source_sha256"], "accepted_iteration": 23,
            "backend": "KS", "d0_floor_ry": d0,
        },
        "fit_registration": {
            "signed_v_ry_in_fit_order": list(a22.SIGNED_V),
            "A1_stage_to_v_ry": {stage: v for v, stage in A1_STAGES.items()},
            "A2_stage_to_v_ry": {stage: v for v, stage in A2_STAGES.items()},
            "delta_definition": "SUMEV_RAW - SUMEV_case", "r_definition": "delta - v*q0",
            "models": "unchanged A1 _series_report fits: signed cubic r(v), signed quadratic sensitivity, all-signed log|r| versus log|v| OLS, and even quartic symmetric-slope fit on seven positive magnitudes",
            "shape_gates": "existing A1 monotonicity/order>=2/plateau<=10D0 plus A2 within-branch constant nonzero sign; exact zero unassessed",
            "d0_role": "A2 additive M0.3 offset and existing plateau scale, not a per-point noise band",
            "aggregate_note": "historical eight-point spin aggregate rows independently checked; six replay cases report physical spins only, aggregate diagnostics are their arithmetic sum",
        },
        "per_site_spin": per_spin,
        "aggregate_spin_diagnostics": aggregate,
        "a21_target_recheck": a21,
        "m03_allowance_at_accepted_u2_targets": allowance,
        "scalar_shift_control_preregistration": {
            "v_ry": scalar_v,
            "per_spin_terms": scalar_terms,
            "sum_abs_r_fit_ry": scalar_model_error,
            "D0_ry": d0,
            "allowance_ry": scalar_allowance,
            "expected_delta_sumev_ry": scalar_v * sum(u2_q0.values()),
            "criterion": "abs(DeltaSUMEV - v*sum_s(q0_s)) <= sum_s(abs(r_fit,s(v))) + D0; q2_local equals q2 per original protocol",
            "basis": "original protocol M0.3 scalar-shift control with A2.2 allowance substitution; v is inside the extended scan, so fitted r is interpolated",
        },
        "wrong_control_preregistration": {
            "criterion": "evaluate V_out and P-transpose controls against the true V_in target DeltaSUMEV=sum_s(t_s); detectable violation means absolute deviation exceeds the same true-target A2.2 allowance above",
            "allowance_ry": allowance["allowance_ry"],
            "basis": "original protocol requires a detectable violation; A2 replaces the M0.3 allowance and supplies no additional threshold",
            "no_measurement_outcome_included": True,
        },
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verdict": verdict, "output": str(args.output), "output_sha256": sha(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
