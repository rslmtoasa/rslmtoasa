#!/usr/bin/env python3
"""Two-phase A5.4 W2 frozen-state prediction and dual replay.

`predict` first checks a single-iteration KS capture and runs the exact LIVE
radial replay twice. It writes the predicted native-vs-P^T (=P*) carrier
separation and a numerical allowance before any dual-control replay.
`replay` consumes that hash-bound prediction and runs the frozen native and
joint wrong-convention cases through the unchanged production
calculate_pl/NEWRHO fixture driver. No SCF or A2 fit is used here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("a5_rs_replay", HERE / "rs_a3_m03_runner.py")
replay = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(replay)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_for(csv_path: Path):
    rows = list(csv.DictReader(csv_path.open(newline="")))
    replay._BY_KIND_STAGE = {}
    for row in rows:
        replay._BY_KIND_STAGE.setdefault((row["kind"], row["stage"]), []).append(row)
        for key in ("weight", "energy", "scalar", "re", "im"):
            if not math.isfinite(float(row[key])):
                raise ValueError(f"nonfinite {row['kind']}/{row['stage']}/{key}")
    its = {int(r["iteration"]) for r in rows}
    if len(its) != 1:
        raise ValueError(f"capture must contain one selected iteration, got {its}")
    has_projector = (any(r["kind"] == "ldm" and r["stage"] == "P_OUT" for r in rows) or
                     any(r["kind"] == "test_density_matrix" and r["stage"] == "P_OUT" for r in rows))
    required = ("ldm_m1", "radial_moment", "radial_parameter", "ql",
                "pl_map_input", "radial_replay_ql", "radial_replay_pl",
                "radial_replay_rho_in", "radial_replay_v", "radial_replay_mesh")
    missing = [k for k in required if not any(r["kind"] == k for r in rows)]
    if not has_projector:
        missing.insert(0, "ldm or test_density_matrix at P_OUT")
    if missing:
        raise ValueError(f"incomplete W2 frozen-state capture: {missing}")
    return rows, next(iter(its))


def make_cases(rows, lam: float):
    nr, ncore, nval = 495, 10, 6
    base = replay.prepare(rows, nr, ncore, nval, "baseline")
    ql_live = base["cases"][0][3]
    ql_native, ql_wrong = ql_live.copy(), ql_live.copy()
    shifts_native, shifts_wrong = np.zeros(2), np.zeros(2)
    w = np.zeros((5, 5), dtype=np.complex128)
    w[0, 1], w[1, 0] = 1j, -1j
    w[2, 2], w[3, 3] = 0.5, -0.25
    per_spin = []
    calibration_cases = [("LIVE_STATE", np.zeros(2), False, ql_live.copy())]
    calibration_keys = {}
    for spin in (1, 2):
        p_kind = "ldm" if any(r["kind"] == "ldm" and r["stage"] == "P_OUT" for r in rows) else "test_density_matrix"
        p_stage = "P_OUT"
        p = replay.matrix(rows, p_kind, p_stage, spin if p_kind == "ldm" else spin - 1, 2, 5)
        m1 = replay.matrix(rows, "ldm_m1", "M1_OUT", spin, 2, 5)
        q0 = float(replay.one(rows, "radial_moment", "P_OUT", spin=spin, ell=2, order=0)["scalar"])
        c = float(replay.one(rows, "radial_moment", "P_OUT", spin=spin, ell=2, order=1)["energy"])
        # This is the centered q2 in the P_OUT radial parameter contract.
        q2 = float(replay.one(rows, "radial_parameter", "P_OUT", spin=spin, ell=2, order=3)["scalar"])
        ql0 = float(replay.one(rows, "ql", "P_OUT", spin=spin, ell=2, order=0)["scalar"])
        if abs(q0 - ql0) > 2e-13:
            raise ValueError(f"P_OUT q0/QL mismatch spin={spin}: {q0} vs {ql0}")
        antiherm = float(np.max(np.abs(p - p.conj().T)))
        if antiherm > 5e-12:
            raise ValueError(f"captured P is not Hermitian within capture precision: {antiherm}")
        cases = []
        for label, pc, mc in (("native", p, m1), ("joint_Ptranspose_Pconjugate", p.T, m1.T)):
            v = lam * w
            t = float(np.trace(v @ pc).real)
            center_shift = t / q0
            q2_local = (q2 - 2.0 * float(np.trace(v @ (mc - c * pc)).real)
                        + float(np.trace(v @ v @ pc).real) - t * t / q0)
            cases.append({"case": label, "t_ry": t, "center_shift_ry": center_shift,
                          "q2_local_ry2": q2_local})
        shifts_native[spin - 1] = cases[0]["center_shift_ry"]
        shifts_wrong[spin - 1] = cases[1]["center_shift_ry"]
        ql_native[2, 2, spin - 1] = cases[0]["q2_local_ry2"]
        ql_wrong[2, 2, spin - 1] = cases[1]["q2_local_ry2"]
        gc = float(base["gc"][2, spin - 1])
        qpar = float(base["qp"][2, spin - 1])
        srdel = float(base["sd"][2, spin - 1])
        radial_scale = max(abs(c - gc), abs(srdel * srdel / qpar), 1.0e-15)
        q2_live = float(ql_live[2, 2, spin - 1])
        for control_index, label in enumerate(("native", "joint_Ptranspose_Pconjugate")):
            center_shift = cases[control_index]["center_shift_ry"]
            q2_delta = cases[control_index]["q2_local_ry2"] - q2_live
            for fraction in (0.25, 0.5):
                for direction in ("center", "q2", "combined"):
                    name = f"CAL_{direction}_{label}_SPIN{spin}_F{fraction:.2f}"
                    shift_vector = np.zeros(2)
                    ql_case = ql_live.copy()
                    if direction in ("center", "combined"):
                        shift_vector[spin - 1] = fraction * center_shift
                    if direction in ("q2", "combined"):
                        ql_case[2, 2, spin - 1] = q2_live + fraction * q2_delta
                    calibration_cases.append((name, shift_vector, True, ql_case))
                    calibration_keys[(spin, label, fraction, direction)] = len(calibration_cases)
            cases[control_index]["radial_scale_ry"] = radial_scale
            cases[control_index]["center_shift_over_radial_scale"] = abs(center_shift) / radial_scale
            cases[control_index]["q2_delta_over_radial_scale2"] = abs(q2_delta) / max(abs(q2), radial_scale ** 2, 1.0e-30)
        per_spin.append({"spin": spin, "q0": q0, "center": c, "q2_centered": q2,
                         "max_abs_imag_P": float(np.max(np.abs(p.imag))),
                         "max_abs_P_minus_PT": float(np.max(np.abs(p - p.T))),
                         "antihermitian_capture_error": antiherm, "controls": cases})
    base["cases"] = [("NATIVE", shifts_native, True, ql_native),
                     ("JOINT_PTRANSPOSE_PCONJUGATE", shifts_wrong, True, ql_wrong)]
    calibration = base.copy()
    calibration["cases"] = calibration_cases
    return base, calibration, per_spin, calibration_keys


def sum_ev(case):
    return sum(case["sumev"].values())


def prediction(csv_path: Path, driver: Path, lam: float, out: Path):
    rows, iteration = rows_for(csv_path)
    if not driver.is_file():
        raise ValueError(f"replay driver not found: {driver}")
    data, calibration, per_spin, calibration_keys = make_cases(rows, lam)
    # Measure the frozen radial inversion residual and deterministic repeat
    # before any convention control is run.
    baseline = replay.prepare(rows, 495, 10, 6, "baseline")
    baseline["cases"] = [("LIVE_STATE", np.zeros(2), False, baseline["cases"][0][3])]
    payload = replay.text_input(baseline, rows, 495, 10, 6)
    outputs, inp_sha, out_sha = [], [], []
    for n in range(2):
        dest = out.parent / f"{out.stem}.live_repeat{n+1}"
        raw, ih, oh = replay.run_persistent(driver, payload, dest)
        result = replay.output_cases(raw)
        outputs.append(result[1]); inp_sha.append(ih); out_sha.append(oh)
    captured_sumev = {int(r["spin"]) + 1: float(r["scalar"]) for r in rows
                      if r["kind"] == "radial_replay_sumev" and r["stage"] == "RADIAL_REPLAY"}
    replay_residual = {s: outputs[0]["sumev"][s] - captured_sumev[s] for s in (1, 2)}
    repeat = {s: outputs[0]["sumev"][s] - outputs[1]["sumev"][s] for s in (1, 2)}
    cal_payload = replay.text_input(calibration, rows, 495, 10, 6)
    cal_dir = out.parent / f"{out.stem}.shift_calibration"
    cal_raw, cal_in_sha, cal_out_sha = replay.run_persistent(driver, cal_payload, cal_dir)
    cal_outputs = replay.output_cases(cal_raw)
    cal_live = cal_outputs[1]
    calibration_by_spin = {}
    for block in per_spin:
        spin = block["spin"]
        details = {}
        for label, case_index in (("native", 0), ("joint_Ptranspose_Pconjugate", 1)):
            ctrl = block["controls"][case_index]
            shift, target_t = ctrl["center_shift_ry"], ctrl["t_ry"]
            by_fraction = {}
            for fraction in (0.25, 0.5):
                center_case = cal_outputs[calibration_keys[(spin, label, fraction, "center")]]
                q2_case = cal_outputs[calibration_keys[(spin, label, fraction, "q2")]]
                combined_case = cal_outputs[calibration_keys[(spin, label, fraction, "combined")]]
                center_observed = cal_live["sumev"][spin] - center_case["sumev"][spin]
                center_target = fraction * shift * block["q0"]
                center_error = abs(center_observed - center_target)
                q2_envelope = abs(cal_live["sumev"][spin] - q2_case["sumev"][spin])
                mixed = (combined_case["sumev"][spin] - center_case["sumev"][spin]
                         - q2_case["sumev"][spin] + cal_live["sumev"][spin])
                by_fraction[str(fraction)] = {"center_carrier_observed_ry": center_observed,
                    "center_carrier_target_ry": center_target, "center_carrier_error_ry": center_error,
                    "q2_only_response_ry": q2_envelope, "mixed_factorial_term_ry": mixed}
            center_linear = [2.0 * by_fraction["0.5"]["center_carrier_error_ry"],
                             4.0 * by_fraction["0.25"]["center_carrier_error_ry"]]
            center_quadratic = [4.0 * by_fraction["0.5"]["center_carrier_error_ry"],
                                16.0 * by_fraction["0.25"]["center_carrier_error_ry"]]
            e_center = max(center_linear + center_quadratic) + max(
                abs(center_linear[i] - center_quadratic[i]) for i in range(2))
            q2_linear = [2.0 * by_fraction["0.5"]["q2_only_response_ry"],
                         4.0 * by_fraction["0.25"]["q2_only_response_ry"]]
            q2_quadratic = [4.0 * by_fraction["0.5"]["q2_only_response_ry"],
                            16.0 * by_fraction["0.25"]["q2_only_response_ry"]]
            e_q2 = max(q2_linear + q2_quadratic) + max(
                abs(q2_linear[i] - q2_quadratic[i]) for i in range(2))
            mixed_half = by_fraction["0.5"]["mixed_factorial_term_ry"]
            mixed_quarter = by_fraction["0.25"]["mixed_factorial_term_ry"]
            mixed_full_half = 4.0 * mixed_half
            mixed_full_quarter = 16.0 * mixed_quarter
            mixed_allowance = max(abs(mixed_full_half), abs(mixed_full_quarter)) + abs(mixed_full_half - mixed_full_quarter)
            details[label] = {"center_shift_ry": shift,
                "center_calibration_error_extrapolated_ry": e_center,
                "center_linear_extrapolations_ry": center_linear,
                "center_quadratic_extrapolations_ry": center_quadratic,
                "q2_only_response_envelope_extrapolated_ry": e_q2,
                "q2_linear_extrapolations_ry": q2_linear,
                "q2_quadratic_extrapolations_ry": q2_quadratic,
                "mixed_term_extrapolated_half_ry": mixed_full_half,
                "mixed_term_extrapolated_quarter_ry": mixed_full_quarter,
                "mixed_term_allowance_ry": mixed_allowance,
                "fractional_calibrations": by_fraction,
                "radial_scale_ry": ctrl["radial_scale_ry"],
                "center_shift_over_radial_scale": ctrl["center_shift_over_radial_scale"],
                "q2_delta_over_radial_scale2": ctrl["q2_delta_over_radial_scale2"],
                "direct_W_carrier_t_ry": target_t}
        calibration_by_spin[spin] = details
    t_native = sum(x["controls"][0]["t_ry"] for x in per_spin)
    t_wrong = sum(x["controls"][1]["t_ry"] for x in per_spin)
    separation_by_spin = {x["spin"]: abs(x["controls"][0]["t_ry"] - x["controls"][1]["t_ry"])
                          for x in per_spin}
    allowance_by_spin = {}
    for spin in (1, 2):
        scale = max(1.0, abs(captured_sumev[spin]),
                    abs(per_spin[spin-1]["controls"][0]["t_ry"]),
                    abs(per_spin[spin-1]["controls"][1]["t_ry"]))
        roundoff = np.finfo(float).eps * scale * 64
        # ES24.16 stores 16 fractional digits with a 0.1..1 normalized
        # mantissa: half an output ulp is at most 5e-16 relative. Bound 64
        # rounded operands/terms in the 5x5 matrix and moment expressions.
        capture_precision = 64.0 * 5.0e-16 * scale
        shift_budget = sum(calibration_by_spin[spin][label]["center_calibration_error_extrapolated_ry"] +
                           calibration_by_spin[spin][label]["q2_only_response_envelope_extrapolated_ry"] +
                           calibration_by_spin[spin][label]["mixed_term_allowance_ry"]
                           for label in ("native", "joint_Ptranspose_Pconjugate"))
        inversion = 2.0 * abs(replay_residual[spin])
        repeat_scatter = abs(repeat[spin])
        carrier_allowance = inversion + shift_budget
        allowance_by_spin[spin] = {"radial_inversion_pair": inversion,
            "shift_dependent_center_and_q2_calibration": shift_budget,
            "repeat_scatter": repeat_scatter, "roundoff": float(roundoff),
            "capture_precision": float(capture_precision),
            "sum": carrier_allowance + repeat_scatter + roundoff + capture_precision,
            "predicted_separation_ry": separation_by_spin[spin],
            "admissibility": "admissible" if separation_by_spin[spin] > carrier_allowance + repeat_scatter + roundoff + capture_precision else "symmetry-/construction-null"}
    separation = abs(t_native - t_wrong)
    prereg = {
        "schema": "A5.4-W2-prediction-v1", "csv_sha256": sha(csv_path),
        "driver_sha256": sha(driver), "selected_iteration": iteration,
        "lambda_ry": lam, "operator": [[0, 1, "+i"], [1, 0, "-i"], [2, 2, 0.5], [3, 3, -0.25]],
        "basis": "complex spherical d block m=-2..+2", "per_spin": per_spin,
        "native_total_t_ry": t_native, "joint_wrong_total_t_ry": t_wrong,
        "predicted_separation_ry": separation,
        "per_spin_admissibility": allowance_by_spin,
        "matched_shift_calibration": calibration_by_spin,
        "admissibility": "at least one spin admissible" if any(v["admissibility"] == "admissible" for v in allowance_by_spin.values()) else "symmetry-/construction-null",
        "baseline_replay": {"sumEV_residuals_vs_capture_ry": replay_residual,
                            "repeat_delta_ry": repeat, "driver_input_sha256": inp_sha,
                            "driver_output_sha256": out_sha,
                            "shift_calibration_input_sha256": cal_in_sha,
                            "shift_calibration_output_sha256": cal_out_sha},
        "decision": "If admissible, native frozen control should follow the captured P convention and joint P^T=P* should separate by predicted amount. If null, report nondiscriminating; no gate.",
        "slope": "HF slope is diagnostic only and is not gated."}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(prereg, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"prediction": str(out), "sha256": sha(out),
                      "separation_ry": separation, "per_spin_admissibility": allowance_by_spin,
                      "admissibility": prereg["admissibility"]}, indent=2))


def do_replay(csv_path: Path, driver: Path, prereg_path: Path, outdir: Path, lam: float):
    prereg = json.loads(prereg_path.read_text())
    if prereg.get("schema") != "A5.4-W2-prediction-v1" or prereg["csv_sha256"] != sha(csv_path):
        raise ValueError("prediction does not bind this capture")
    if prereg["driver_sha256"] != sha(driver) or float(prereg["lambda_ry"]) != lam:
        raise ValueError("driver/lambda differs from preregistered W2 prediction")
    rows, _ = rows_for(csv_path)
    data, _, _, _ = make_cases(rows, lam)
    payload = replay.text_input(data, rows, 495, 10, 6)
    raw, inp_sha, out_sha = replay.run_persistent(driver, payload, outdir)
    cases = replay.output_cases(raw)
    native, wrong = cases[1], cases[2]
    measured = sum_ev(wrong) - sum_ev(native)
    predicted_native_minus_wrong = prereg["joint_wrong_total_t_ry"] - prereg["native_total_t_ry"]
    report = {"prediction_sha256": sha(prereg_path), "capture_sha256": sha(csv_path),
              "driver_sha256": sha(driver), "measured_native_minus_joint_wrong_ry": -measured,
              "predicted_native_minus_joint_wrong_ry": predicted_native_minus_wrong,
              "absolute_prediction_residual_ry": abs(-measured - predicted_native_minus_wrong),
              "per_spin_admissibility": prereg["per_spin_admissibility"], "admissibility": prereg["admissibility"],
              "native_sumEV_by_spin_ry": native["sumev"], "joint_wrong_sumEV_by_spin_ry": wrong["sumev"],
              "driver_input_sha256": inp_sha, "driver_output_sha256": out_sha,
              "interpretation": "state-level complex-convention diagnostic only; no production gate"}
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "result.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("predict")
    p.add_argument("--csv", type=Path, required=True); p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--lambda-ry", type=float, required=True); p.add_argument("--out", type=Path, required=True)
    r = sub.add_parser("replay")
    r.add_argument("--csv", type=Path, required=True); r.add_argument("--driver", type=Path, required=True)
    r.add_argument("--prediction", type=Path, required=True); r.add_argument("--outdir", type=Path, required=True)
    r.add_argument("--lambda-ry", type=float, required=True)
    a = ap.parse_args()
    if a.cmd == "predict": prediction(a.csv, a.driver, a.lambda_ry, a.out)
    else: do_replay(a.csv, a.driver, a.prediction, a.outdir, a.lambda_ry)


if __name__ == "__main__":
    main()
