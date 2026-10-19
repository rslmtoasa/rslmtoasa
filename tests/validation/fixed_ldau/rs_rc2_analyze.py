#!/usr/bin/env python3
"""Analyze completed C3 RC2 endpoint runs; never launches an SCF.

The analyzer joins each case's selected terminal energy record to the same
iteration's complete P_OUT trace, obtains fixed-P production E_Hub derivatives
from the central local endpoints, then evaluates the registered allowance.
Incomplete campaigns are written as blocked results rather than guessed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CAPTURE = ROOT / "tests/validation/fixed_ldau/rs_rc2_fixed_p_capture.py"
DEFAULT_RC_JSON = ROOT / "docs/validation/fixed_ldau/rs/C3_RC_PROVENANCE.json"
FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")
ENERGY_PREFIX = "HUBBARD_MOMENT_ENERGY_TOTAL"
DIFF_RE = re.compile(r"(?:Not converged!\s*)?Diff=\s*([-+0-9.eEdD]+)")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite_float(value: str, where: str) -> float:
    normalized = value.strip().replace("D", "E").replace("d", "e")
    if not FLOAT_RE.fullmatch(normalized):
        raise ValueError(f"missing/non-numeric value for {where}: {value!r}")
    number = float(normalized)
    if not math.isfinite(number):
        raise ValueError(f"non-finite value for {where}: {value!r}")
    return number


def parse_assignments(line: str) -> dict[str, str]:
    return {key: value.rstrip(",;") for key, value in
            re.findall(r"([A-Za-z][A-Za-z0-9_]*)=\s*([^\s]+)", line)}


def read_energy_rows(log_path: Path) -> list[dict]:
    rows = []
    for line_no, line in enumerate(log_path.read_text(errors="replace").splitlines(), 1):
        if ENERGY_PREFIX not in line:
            continue
        fields = parse_assignments(line)
        row = {"log_line": line_no, "raw_line_sha256": hashlib.sha256(line.encode()).hexdigest()}
        for key in ("physical_total_energy", "candidate_a", "candidate_b", "candidate_c"):
            value = fields.get(key, "UNAVAILABLE")
            if value.upper().startswith("UNAVAILABLE"):
                row[key] = None
            else:
                row[key] = finite_float(value, f"{key} at {log_path}:{line_no}")
        rows.append(row)
    return rows


def selected_phase(case: dict) -> tuple[dict, dict, str]:
    phases = case.get("phases")
    if not isinstance(phases, list) or not phases:
        raise ValueError(f"{case.get('case')}: no phase records")
    accepted = [p for p in phases if p.get("acceptance", {}).get("accepted") is True]
    # Candidate comparisons are preregistered against the earliest C1-
    # accepted phase. Preserve later phases as provenance, but never let an
    # additional recovery continuation silently replace that endpoint.
    phase = accepted[0] if accepted else phases[-1]
    manifest = phase.get("solver_manifest")
    if not isinstance(manifest, dict):
        run_dir = Path(phase.get("run_dir", ""))
        mp = Path(phase.get("trace_manifest_path", run_dir / "TRACE_ARCHIVE_MANIFEST.json"))
        if not mp.is_file():
            raise ValueError(f"{case.get('case')}: missing runner manifest {mp}")
        manifest = json.loads(mp.read_text())
    acceptance = phase.get("acceptance", {})
    status = "C1_ACCEPTED" if acceptance.get("accepted") is True else "PROVISIONAL_OR_UNACCEPTED"
    return phase, manifest, status


def read_case_endpoint(case: dict) -> dict:
    phase, manifest, endpoint_status = selected_phase(case)
    run = Path(phase.get("run_dir", manifest.get("run_dir", ""))).resolve()
    trace = Path(manifest.get("aggregate_path", run / "fixed_u_trace_history.csv")).resolve()
    log = Path(manifest.get("console_path", run / "console.log")).resolve()
    if not trace.is_file() or not log.is_file():
        raise ValueError(f"{case.get('case')}: missing trace or console ({trace}, {log})")
    phase_manifest_path = Path(phase.get("trace_manifest_path", run / "TRACE_ARCHIVE_MANIFEST.json"))
    if phase_manifest_path.is_file() and phase.get("trace_manifest_sha256"):
        if sha(phase_manifest_path) != phase["trace_manifest_sha256"]:
            raise ValueError(f"{case.get('case')}: phase runner-manifest hash mismatch")

    iterations = manifest.get("captured_iterations", [])
    if not iterations or iterations != list(range(1, int(iterations[-1]) + 1)):
        raise ValueError(f"{case.get('case')}: captured iteration IDs are not a contiguous 1-based prefix")
    endpoint_it = int(iterations[-1])
    live_trace = trace
    checkpoint_trace_used = False
    ignored_live_suffix = []
    if endpoint_status == "C1_ACCEPTED":
        stated = phase.get("acceptance", {}).get("window_iterations", [])
        checkpoint = manifest.get("early_c1_checkpoint")
        if manifest.get("intentional_accepted_stop") is True and checkpoint:
            endpoint_it = int(checkpoint.get("iteration", -1))
            checkpoint_dir = Path(checkpoint.get("checkpoint_dir", ""))
            checkpoint_trace = checkpoint_dir / "fixed_u_trace_history.csv"
            expected_hash = checkpoint.get("files_sha256", {}).get("fixed_u_trace_history.csv")
            if (endpoint_it < 20 or not checkpoint_trace.is_file() or not expected_hash or
                    sha(checkpoint_trace) != expected_hash or
                    checkpoint.get("trace_sha256") != expected_hash):
                raise ValueError(f"{case.get('case')}: frozen C1 checkpoint trace is missing or hash-mismatched")
            checkpoint_window = checkpoint.get("window_iterations", [])
            if (not stated or stated != checkpoint_window or int(stated[-1]) != endpoint_it or
                    stated != list(range(endpoint_it - 19, endpoint_it + 1))):
                raise ValueError(f"{case.get('case')}: accepted C1 checkpoint window is not a contiguous final 20")
            # The live rolling writer may complete one more P_OUT after the
            # wrapper freezes its accepted checkpoint. Use the hash-bound
            # checkpoint history and ignore only this post-checkpoint suffix.
            if endpoint_it > int(iterations[-1]):
                raise ValueError(f"{case.get('case')}: accepted checkpoint exceeds the live trace")
            ignored_live_suffix = [int(it) for it in iterations if int(it) > endpoint_it]
            trace = checkpoint_trace
            checkpoint_trace_used = True
        elif not stated or int(stated[-1]) != endpoint_it:
            raise ValueError(f"{case.get('case')}: accepted iteration does not match archived trace suffix")

    by_it: dict[int, list[dict[str, str]]] = {}
    with trace.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("stage") == "P_OUT":
                try:
                    it = int(row["iteration"])
                except (KeyError, ValueError) as exc:
                    raise ValueError(f"{case.get('case')}: malformed P_OUT iteration") from exc
                by_it.setdefault(it, []).append(row)
    trace_iterations = sorted(by_it)
    if (not trace_iterations or trace_iterations != list(range(1, endpoint_it + 1)) or
            endpoint_it not in by_it):
        raise ValueError(f"{case.get('case')}: selected P_OUT history is not a contiguous prefix through accepted iteration")
    records = by_it.get(endpoint_it, [])
    matrices = [r for r in records if r.get("kind") == "test_density_matrix" and r.get("site") == "0" and r.get("l") == "2"]
    ql = [r for r in records if r.get("kind") == "ql" and r.get("site") == "0"]
    wp = [r for r in records if r.get("kind") == "test_trace_wp" and r.get("site") == "0"]
    if len(matrices) != 50 or len(ql) != 6 or len(wp) != 2:
        raise ValueError(f"{case.get('case')}: terminal P_OUT incomplete: matrix={len(matrices)}, ql={len(ql)}, W1={len(wp)}")
    perchecks = manifest.get("per_iteration_matrix_checks", {}).get(str(endpoint_it), {})
    if not perchecks or perchecks.get("pass_1e_11") is not True:
        raise ValueError(f"{case.get('case')}: terminal P_OUT matrix checks fail: {perchecks}")

    energy_rows = read_energy_rows(log)
    diff_values = [finite_float(x, f"Diff in {log}") for x in DIFF_RE.findall(log.read_text(errors="replace"))]
    if len(energy_rows) < endpoint_it or len(diff_values) < endpoint_it:
        raise ValueError(f"{case.get('case')}: energy/Diff log ends before selected P_OUT iteration {endpoint_it}: energies={len(energy_rows)}, diffs={len(diff_values)}")
    energy = energy_rows[endpoint_it - 1]
    selected_diff = diff_values[endpoint_it - 1]
    if endpoint_status == "C1_ACCEPTED":
        expected_diff = float(phase.get("acceptance", {}).get("final_diff",
                            manifest.get("early_c1_checkpoint", {}).get("final_diff", math.nan)))
        if not math.isfinite(expected_diff) or abs(selected_diff - expected_diff) > 1e-15:
            raise ValueError(f"{case.get('case')}: selected Diff does not match frozen C1 checkpoint")
    expected = "candidate_b" if case.get("mode") == "local" else "candidate_c"
    if energy.get(expected) is None:
        raise ValueError(f"{case.get('case')}: {expected} unavailable at selected iteration {endpoint_it}")
    return {
        "case": case["case"], "family": case["family"], "mode": case["mode"],
        "u_eV": float(case["u_eV"]), "j_eV": float(case["j_eV"]),
        "status": endpoint_status,
        "phase": phase.get("phase"), "phase_dir": str(run),
        "phase_selection_policy": "earliest C1-accepted phase; if none, terminal recorded phase is provisional",
        "all_phase_records": [
            {
                "phase": p.get("phase"), "run_dir": p.get("run_dir"),
                "acceptance": p.get("acceptance", {}),
                "runner_returncode": p.get("runner_returncode"),
                "trace_manifest_path": p.get("trace_manifest_path"),
                "trace_manifest_sha256": p.get("trace_manifest_sha256"),
                "solver_manifest_trace_sha256": p.get("solver_manifest", {}).get("aggregate_sha256"),
            }
            for p in case.get("phases", [])
        ],
        "phase_acceptance": phase.get("acceptance", {}),
        "run_manifest_sha256": sha(phase_manifest_path) if phase_manifest_path.is_file() else None,
        "trace_path": str(trace), "trace_sha256": sha(trace),
        "live_trace_path": str(live_trace), "live_trace_sha256": sha(live_trace),
        "checkpoint_trace_used": checkpoint_trace_used,
        "ignored_live_post_checkpoint_iterations": ignored_live_suffix,
        "console_path": str(log), "console_sha256": sha(log),
        "selected_iteration": endpoint_it,
        "energy_log_line": energy["log_line"], "energy_row_sha256": energy["raw_line_sha256"],
        "physical_total_energy_Ry": energy["physical_total_energy"],
        "candidate_a_Ry": energy["candidate_a"], "candidate_b_Ry": energy["candidate_b"],
        "candidate_c_Ry": energy["candidate_c"], "selected_diff": selected_diff,
        "matrix_entry_count": len(matrices), "ql_entry_count": len(ql), "w1_trace_count": len(wp),
    }


def evaluator_values(stdout_path: Path) -> dict[float, dict[str, float]]:
    values: dict[float, dict[str, float]] = {}
    for line in stdout_path.read_text(errors="replace").splitlines():
        if "h_eV=" not in line:
            continue
        fields = parse_assignments(line)
        h = finite_float(fields["h_eV"], "fixed-P h")
        values[h] = {
            "dEHub_dU_Ry_per_eV": finite_float(fields["dEHub_dU_Ry_per_eV"], "fixed-P dEHub/dU"),
            "dEHub_dJ_Ry_per_eV": finite_float(fields["dEHub_dJ_Ry_per_eV"], "fixed-P dEHub/dJ"),
        }
    if set(values) != {0.05, 0.1}:
        raise ValueError(f"fixed-P evaluator must return h=.05 and .10; got {sorted(values)}")
    return values


def run_fixed_p_capture(center: dict, driver: Path, output_dir: Path) -> dict:
    capture_dir = output_dir / "fixed_p" / center["case"]
    result_path = capture_dir / "FIXED_P_RESULT.json"
    if result_path.is_file():
        result = json.loads(result_path.read_text())
        expected_trace_sha = center["trace_sha256"]
        expected_driver_sha = sha(driver)
        if (result.get("trace_sha256") != expected_trace_sha or
                result.get("iteration") != center["selected_iteration"] or
                result.get("driver_sha256") != expected_driver_sha):
            raise ValueError(f"cached fixed-P result does not match current endpoint: {result_path}")
        stdout = Path(result["stdout_path"])
        return {
            "case": center["case"], "status": center["status"],
            "selected_iteration": center["selected_iteration"],
            "matrix_path": result["P_matrix_path"], "matrix_sha256": result["P_matrix_sha256"],
            "matrix_entry_count": result["matrix_entry_count"],
            "result_path": str(result_path), "result_sha256": sha(result_path),
            "capture_script_sha256": sha(CAPTURE),
            "driver_path": str(driver.resolve()), "driver_sha256": expected_driver_sha,
            "evaluator_stdout_path": str(stdout), "evaluator_stdout_sha256": sha(stdout),
            "fixed_P_derivatives": evaluator_values(stdout),
        }
    cmd = [sys.executable, str(CAPTURE), center["trace_path"], str(center["selected_iteration"]),
           f"{center['u_eV']:.16g}", f"{center['j_eV']:.16g}", "--site", "0",
           "--driver", str(driver.resolve()), "--outdir", str(capture_dir)]
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError(f"fixed-P capture failed for {center['case']} (rc={proc.returncode}): {proc.stderr[-2000:]}")
    result = json.loads(result_path.read_text())
    if result.get("returncode") != 0 or result.get("matrix_entry_count") != 50:
        raise ValueError(f"invalid fixed-P capture result {result_path}")
    stdout = Path(result["stdout_path"])
    return {
        "case": center["case"], "status": center["status"],
        "selected_iteration": center["selected_iteration"],
        "matrix_path": result["P_matrix_path"], "matrix_sha256": result["P_matrix_sha256"],
        "matrix_entry_count": result["matrix_entry_count"],
        "result_path": str(result_path), "result_sha256": sha(result_path),
        "capture_script_sha256": sha(CAPTURE),
        "driver_path": str(driver.resolve()), "driver_sha256": sha(driver),
        "evaluator_stdout_path": str(stdout), "evaluator_stdout_sha256": sha(stdout),
        "fixed_P_derivatives": evaluator_values(stdout),
    }


def slope(minus: float, plus: float, h: float) -> float:
    return (plus - minus) / (2.0 * h)


def centered_slopes(index: dict, family: str, mode: str, candidate: str,
                    center: float, fixed: float, variable: str) -> dict[str, float]:
    result = {}
    for h in (0.05, 0.10):
        if variable == "u":
            minus = index[(family, round(center - h, 8), round(fixed, 8), mode)]
            plus = index[(family, round(center + h, 8), round(fixed, 8), mode)]
        else:
            minus = index[(family, round(fixed, 8), round(center - h, 8), mode)]
            plus = index[(family, round(fixed, 8), round(center + h, 8), mode)]
        minus_value = minus[f"candidate_{candidate}_Ry"]
        plus_value = plus[f"candidate_{candidate}_Ry"]
        if minus_value is None or plus_value is None:
            raise ValueError(f"candidate {candidate} unavailable for {family} h={h}")
        result[f"h_{h:.2f}"] = slope(minus_value, plus_value, h)
    return result


def endpoint_index(endpoints: list[dict]) -> dict:
    index = {}
    for item in endpoints:
        key = (item["family"], round(item["u_eV"], 8), round(item["j_eV"], 8), item["mode"])
        if key in index:
            raise ValueError(f"duplicate candidate endpoint key: {key}")
        index[key] = item
    return index


def analyze_u_derivative(index: dict, fixed_p: dict, r_rc1: float) -> dict:
    u_candidates = {
        "a": centered_slopes(index, "U", "spectral", "a", 2.0, 0.0, "u"),
        "b": centered_slopes(index, "U", "local", "b", 2.0, 0.0, "u"),
        "c": centered_slopes(index, "U", "spectral", "c", 2.0, 0.0, "u"),
    }
    u_target = fixed_p["u2_j0"]["fixed_P_derivatives"]
    u_A = 3.0 * abs(r_rc1) * abs(u_target[0.05]["dEHub_dU_Ry_per_eV"]) + abs(
        u_candidates["b"]["h_0.10"] - u_candidates["b"]["h_0.05"])
    u_per_h = {}
    for h in (0.05, 0.10):
        key = f"h_{h:.2f}"
        target = u_target[h]["dEHub_dU_Ry_per_eV"]
        bdist = abs(u_candidates["b"][key] - target)
        cdist = abs(u_candidates["c"][key] - target)
        u_per_h[key] = {
            "target_Ry_per_eV": target,
            "candidate_slopes_Ry_per_eV": {c: u_candidates[c][key] for c in ("a", "b", "c")},
            "b_within_A": bdist <= u_A, "b_abs_error_Ry_per_eV": bdist,
            "c_miss_gt_10A": cdist > 10.0 * u_A, "c_abs_error_Ry_per_eV": cdist,
        }
    return {
        "center_U_eV": 2.0, "fixed_J_eV": 0.0, "fixed_P_target": u_target,
        "candidate_slopes": u_candidates, "common_allowance_A_Ry_per_eV": u_A,
        "per_h": u_per_h,
        "gate": "PASS" if all(x["b_within_A"] and x["c_miss_gt_10A"] for x in u_per_h.values()) else "FAIL",
    }


def analyze_j_derivative(index: dict, fixed_p: dict, r_rc1: float, center_j: float) -> dict:
    family = f"J{center_j:g}"
    b_slopes = centered_slopes(index, family, "local", "b", center_j, 2.0, "j")
    a_slopes = centered_slopes(index, family, "spectral", "a", center_j, 2.0, "j")
    c_slopes = centered_slopes(index, family, "spectral", "c", center_j, 2.0, "j")
    target = fixed_p[f"u2_j{center_j:g}"]["fixed_P_derivatives"]
    normalizer = abs(target[0.05]["dEHub_dU_Ry_per_eV"])
    allowance = 3.0 * abs(r_rc1) * normalizer + abs(b_slopes["h_0.10"] - b_slopes["h_0.05"])
    per_h = {}
    for h in (0.05, 0.10):
        key = f"h_{h:.2f}"
        j_target = target[h]["dEHub_dJ_Ry_per_eV"]
        bdist = abs(b_slopes[key] - j_target)
        cdist = abs(c_slopes[key] - j_target)
        per_h[key] = {
            "target_TJ_Ry_per_eV": j_target,
            "fixed_P_U_normalizer_TU_h005_Ry_per_eV": normalizer,
            "candidate_slopes_Ry_per_eV": {"a": a_slopes[key], "b": b_slopes[key], "c": c_slopes[key]},
            "b_within_A": bdist <= allowance, "b_abs_error_Ry_per_eV": bdist,
            "c_miss_gt_10A": cdist > 10.0 * allowance, "c_abs_error_Ry_per_eV": cdist,
        }
    return {
        "center_U_eV": 2.0, "center_J_eV": center_j, "fixed_P_target": target,
        "candidate_slopes": {"a": a_slopes, "b": b_slopes, "c": c_slopes},
        "common_allowance_A_Ry_per_eV": allowance, "per_h": per_h,
        "gate": "PASS" if all(x["b_within_A"] and x["c_miss_gt_10A"] for x in per_h.values()) else "FAIL",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign", type=Path, required=True)
    ap.add_argument("--driver", type=Path, required=True, help="pinned fixed-P production energy evaluator")
    ap.add_argument("--rc-json", type=Path, default=DEFAULT_RC_JSON)
    ap.add_argument("--out", type=Path, required=True, help="new result JSON path")
    ap.add_argument("--partial", action="store_true",
                    help="analyze only families with all required completed endpoints; no final verdict")
    args = ap.parse_args()
    output = {"schema": "c3-rc2-endpoint-analysis-v1", "campaign_path": str(args.campaign.resolve()),
              "campaign_sha256": sha(args.campaign), "driver_path": str(args.driver.resolve()),
              "driver_sha256": sha(args.driver) if args.driver.is_file() else None,
              "analyzer_path": str(Path(__file__).resolve()), "analyzer_sha256": sha(Path(__file__).resolve()),
              "status": "BLOCKED", "blockers": []}
    try:
        campaign = json.loads(args.campaign.read_text())
        if not args.partial and campaign.get("status") not in {"COMPLETE", "COMPLETE_WITH_ERRORS"}:
            raise ValueError(f"campaign is not complete: {campaign.get('status')!r}")
        if (not isinstance(campaign.get("cases"), list) or len(campaign["cases"]) > 30 or
                (not args.partial and len(campaign["cases"]) != 30)):
            raise ValueError(f"expected 30 completed campaign case records, got {len(campaign.get('cases', []))}")
        if not args.driver.is_file():
            raise FileNotFoundError(args.driver)
        rc_data = json.loads(args.rc_json.read_text())
        prereg = rc_data["rc2_preregistration"]
        gate = rc_data["rc1"]["final_gate"]
        if gate.get("eligible_h005_all_endpoints_c1_accepted") is not True:
            raise ValueError("RC1 h=.005 endpoint is not C1-qualified for the allowance scale")
        r_rc1 = abs(finite_float(str(gate["eligible_h005_absolute_relative_residual"]),
                                 "eligible RC1 h=.005 relative residual"))
        staged_path = Path(campaign.get("staged_manifest_path", ""))
        staged = json.loads(staged_path.read_text()) if staged_path.is_file() else None
        if staged is None or len(staged.get("cases", [])) != 30:
            raise ValueError(f"staged-deck manifest unavailable/incomplete: {staged_path}")
        if campaign.get("staged_manifest_sha256") != sha(staged_path):
            raise ValueError("campaign/staged-deck manifest hash mismatch")
        by_name = {x["case"]: x for x in staged["cases"]}
        if len(by_name) != 30:
            raise ValueError("duplicate case identifiers in staged manifest")
        campaign_names = [x.get("case") for x in campaign["cases"]]
        if len(set(campaign_names)) != len(campaign_names) or not set(campaign_names).issubset(by_name):
            raise ValueError("campaign case set contains duplicates or cases outside the frozen staged manifest")
        if not args.partial and (len(campaign_names) != 30 or set(campaign_names) != set(by_name)):
            raise ValueError("campaign case set differs from the frozen 30-case staged manifest")
        endpoints = []
        case_errors = []
        for case in campaign["cases"]:
            try:
                if case.get("status") == "CAMPAIGN_DRIVER_ERROR":
                    raise ValueError(case.get("error", "campaign driver error"))
                if case.get("case") not in by_name:
                    raise ValueError("case absent from frozen staged manifest")
                staged_case = by_name[case["case"]]
                for field in ("family", "mode", "u_eV", "j_eV"):
                    if field not in case or case[field] != staged_case[field]:
                        raise ValueError(f"campaign/staged {field} mismatch")
                endpoints.append(read_case_endpoint(case))
            except Exception as exc:
                case_errors.append({"case": case.get("case"), "error": str(exc)})
        output.update({"campaign_status": campaign.get("status"), "campaign_completed_unix": campaign.get("completed_unix"),
                       "staged_manifest_path": str(staged_path), "staged_manifest_sha256": sha(staged_path),
                       "rc1_status": gate.get("status"), "rc1_relative_residuals": gate.get("relative_residuals"),
                       "registered_allowance_formulas": prereg.get("allowance"),
                       "registered_candidate_definitions": prereg.get("candidate_definitions"),
                       "r_rc1_used_max_abs": r_rc1, "case_endpoints": endpoints, "case_errors": case_errors,
                       "fixed_P_captures": {}, "derivatives": None})
        if not args.partial and (case_errors or len(endpoints) != 30):
            raise ValueError(f"incomplete/invalid case endpoint set ({len(endpoints)}/30)")
        capture_specs = {
            "U": ("u2_j0", "U", 2.0, 0.0),
            "J0.5": ("u2_j0.5", "J0.5", 2.0, 0.5),
            "J0.9": ("u2_j0.9", "J0.9", 2.0, 0.9),
        }
        endpoint_by_key = {
            (x["family"], round(x["u_eV"], 8), round(x["j_eV"], 8), x["mode"]): x
            for x in endpoints
        }
        endpoint_by_name = {x["case"]: x for x in endpoints}
        errors_by_name = {x["case"]: x["error"] for x in case_errors}
        available_fixed_p = {}
        family_results = {}
        family_progress = {}
        assets = args.out.parent / (args.out.stem + "_assets")
        for family, (capture_key, _, center_u, center_j) in capture_specs.items():
            required = [x["case"] for x in staged["cases"] if x["family"] == family]
            present = [name for name in required if name in endpoint_by_name]
            failed = [name for name in required if name in errors_by_name]
            missing = [name for name in required if name not in endpoint_by_name and name not in errors_by_name]
            complete = len(required) == 10 and len(present) == 10 and not failed
            family_progress[family] = {
                "status": "READY_FOR_DIAGNOSTIC" if complete else "NOT_YET_EVALUABLE",
                "gate": "NOT_YET_EVALUABLE",
                "required_case_count": len(required), "valid_endpoint_count": len(present),
                "missing_cases": missing, "invalid_cases": [{"case": n, "error": errors_by_name[n]} for n in failed],
            }
            if not complete:
                continue
            center_key = (family, round(center_u, 8), round(center_j, 8), "local")
            if center_key not in endpoint_by_key:
                family_progress[family]["status"] = "NOT_YET_EVALUABLE"
                family_progress[family]["missing_central_local_endpoint"] = list(center_key)
                continue
            center = endpoint_by_key[center_key]
            capture = run_fixed_p_capture(center, args.driver, assets)
            available_fixed_p[capture_key] = capture
            result = (analyze_u_derivative(endpoint_index(endpoints), {capture_key: available_fixed_p[capture_key]}, r_rc1)
                      if family == "U" else
                      analyze_j_derivative(endpoint_index(endpoints), {capture_key: available_fixed_p[capture_key]}, r_rc1, center_j))
            all_family_accepted = all(endpoint_by_name[n]["status"] == "C1_ACCEPTED" for n in required)
            family_results[family] = {
                **result,
                "endpoint_status": "ALL_C1_ACCEPTED" if all_family_accepted else "HAS_PROVISIONAL_OR_UNACCEPTED_ENDPOINTS",
                "diagnostic_gate": result["gate"],
                "final_gate_verdict": ("NOT_FINALIZED_PARTIAL_MODE" if args.partial else result["gate"]),
            }
            family_progress[family]["status"] = "DIAGNOSTIC_COMPLETE"

        all_campaign_cases = (campaign.get("status") in {"COMPLETE", "COMPLETE_WITH_ERRORS"} and
                              len(campaign_names) == 30 and set(campaign_names) == set(by_name))
        output["fixed_P_captures"] = available_fixed_p
        output["derivative_families"] = family_results
        output["family_progress"] = family_progress
        output["r_rc1_used_eligible_h005_abs"] = r_rc1
        output["endpoint_status"] = ("ALL_C1_ACCEPTED" if endpoints and
                                     all(x["status"] == "C1_ACCEPTED" for x in endpoints)
                                     else "HAS_PROVISIONAL_OR_UNACCEPTED_ENDPOINTS")
        if args.partial or not all_campaign_cases:
            output["status"] = "PARTIAL_PROGRESS"
            output["final_gate_verdict"] = "NOT_YET_EVALUABLE"
        elif case_errors or len(endpoints) != 30 or len(family_results) != 3:
            output["status"] = "BLOCKED"
            output["final_gate_verdict"] = "NOT_YET_EVALUABLE"
        else:
            all_accepted = all(x["status"] == "C1_ACCEPTED" for x in endpoints)
            output["status"] = "ANALYZED" if all_accepted else "ANALYZED_WITH_PROVISIONAL_ENDPOINTS"
            output["final_gate_verdict"] = (
                "PASS" if all(v["gate"] == "PASS" for v in family_results.values()) else "FAIL"
            ) if all_accepted else "NOT_FINALIZED_PROVISIONAL_ENDPOINTS"
    except Exception as exc:
        output["blockers"].append(str(exc))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"], "out": str(args.out.resolve()), "blockers": output["blockers"]}, indent=2))
    if output["status"] == "BLOCKED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
