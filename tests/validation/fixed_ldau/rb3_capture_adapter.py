#!/usr/bin/env python3
"""Validate, replay, and analyze one accepted RB2 full-trace capture for RB3.

The ordered modes establish a LIVE_STATE baseline, preregister capture-derived
operands and allowances, then run the registered response/control cases. Every
mode accepts only a hash-bound C1-accepted RB2 state and full CAPTURE_ONLY trace.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

FIT_DEFAULT = (Path(__file__).resolve().parents[3] /
               "docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_fit_results.json")
FIT_SHA256 = "0f0b543f81f53bc1c7fdfc46d3b404adef1ca455b8366a90e87154afd2b0987c"
CSV_FIELDS = ("kind", "run_id", "lambda", "iteration", "stage", "site", "type", "l", "spin",
              "k", "band", "row", "col", "order", "weight", "energy", "scalar", "re", "im")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def json_default(value: Any) -> Any:
    if isinstance(value, complex):
        return {"re": value.real, "im": value.imag}
    if hasattr(value, "tolist"):
        return value.tolist()
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def validate(args: argparse.Namespace) -> dict[str, Any]:
    lock_path = args.lock.resolve()
    capture_dir = args.capture.resolve()
    lock = load_json(lock_path)
    require(lock.get("schema") == "rslmto-rb3-capture-lock-v2", "unsupported RB3 capture lock")
    require(lock.get("status") == "LOCKED_BEFORE_RB3_RESPONSE", "RB3 lock must be frozen before any response run")
    require(lock.get("backend") == "RS", "RB3 requires the accepted RS RB2 endpoint")
    require(lock.get("capture_only") is True and lock.get("ql_only") is False,
            "capture lock must attest full CAPTURE_ONLY=1 with QL_ONLY disabled")
    convergence = lock.get("convergence", {})
    require(convergence.get("outer_diff") is not None and float(convergence["outer_diff"]) <= 3e-8,
            "RB2 endpoint does not satisfy C1 outer Diff <= 3e-8")
    occ = convergence.get("last20_d_occupation_by_spin", {})
    require(set(occ) == {"1", "2"} and all(len(occ[s]) == 20 for s in ("1", "2")),
            "C1 lock must record 20 d-occupation samples for each spin")
    occ_spreads = {s: max(map(float, occ[s])) - min(map(float, occ[s])) for s in ("1", "2")}
    require(max(occ_spreads.values()) <= 5e-7, "RB2 endpoint fails C1 last-20 d-occupation spread")
    moments = convergence.get("last20_sphere_moment_muB", [])
    require(len(moments) == 20, "C1 lock must record 20 sphere-moment samples")
    moment_spread = max(map(float, moments)) - min(map(float, moments))
    require(moment_spread <= 1e-6, "RB2 endpoint fails C1 last-20 sphere-moment spread")

    fit = args.fit.resolve()
    fit_sha = sha256(fit)
    require(fit_sha == FIT_SHA256, "locked RS A2.2 fit SHA mismatch")
    require(lock.get("fit_sha256") == FIT_SHA256,
            "pre-response lock does not bind the locked RS A2.2 fit")

    # The accepted production state and the full diagnostic capture have
    # different executable/source identities. The capture build changes only
    # bounds-checked trace helpers; it is never represented as the RB2 solver.
    rb2_binary = Path(args.rb2_binary).resolve()
    capture_source, capture_binary = Path(args.source).resolve(), Path(args.binary).resolve()
    rb2_binary_sha = sha256(rb2_binary)
    source_sha, binary_sha = sha256(capture_source), sha256(capture_binary)
    production = lock.get("rb2_production", {})
    capture_build = lock.get("capture_build", {})
    production_manifest_path = Path(production.get("manifest_path", "")).resolve()
    production_manifest = load_json(production_manifest_path)
    production_case = production_manifest.get("cases", {}).get(production.get("case_key"), {})
    require(production.get("manifest_sha256") == sha256(production_manifest_path),
            "frozen RB2 endpoint manifest hash mismatch")
    require(production_manifest.get("status") == "C1_ACCEPTED_RB2_ENDPOINTS_FROZEN",
            "RB2 source-of-truth manifest is not frozen/accepted")
    rb2_source_sha = production_manifest.get("source_self_sha256")
    require(production.get("source_sha256") == rb2_source_sha, "RB2 production source hash mismatch")
    require(production.get("executable_sha256") == rb2_binary_sha, "RB2 production executable hash mismatch")
    require(production_case.get("binary_sha256") == rb2_binary_sha, "RB2 case manifest binary mismatch")
    require(production_case.get("initial_input_sha256") == lock.get("initial_input_sha256"),
            "RB2 case manifest initial input mismatch")
    require(production_case.get("initial_seed_sha256") == lock.get("initial_seed_sha256"),
            "RB2 case manifest initial seed mismatch")
    require(capture_build.get("source_sha256") == source_sha, "capture diagnostic source hash mismatch")
    require(capture_build.get("executable_sha256") == binary_sha, "capture diagnostic executable hash mismatch")
    require(capture_build.get("diagnostic_only") is True,
            "capture source delta must be explicitly classified as diagnostic-only")
    driver = args.driver.resolve()
    require(lock.get("driver_sha256") == sha256(driver), "frozen replay-driver hash differs from RB3 lock")
    adapter_sha = sha256(Path(__file__).resolve())
    require(lock.get("adapter_sha256") == adapter_sha, "RB3 adapter hash differs from pre-response lock")

    input_path = Path(lock["initial_input_path"]).resolve()
    seed_path = Path(lock["initial_seed_path"]).resolve()
    input_sha, seed_sha = sha256(input_path), sha256(seed_path)
    require(lock.get("initial_input_sha256") == input_sha, "RB2 initial deck hash mismatch")
    require(lock.get("initial_seed_sha256") == seed_sha, "RB2 initial seed hash mismatch")
    capture_execution = lock.get("capture_execution", {})
    capture_input = Path(capture_execution.get("input_path", "")).resolve()
    checkpoint = Path(capture_execution.get("checkpoint_path", "")).resolve()
    capture_input_sha, checkpoint_sha = sha256(capture_input), sha256(checkpoint)
    require(capture_execution.get("input_sha256") == capture_input_sha,
            "one-step capture deck hash mismatch")
    require(capture_execution.get("accepted_checkpoint_sha256") == checkpoint_sha,
            "one-step capture checkpoint differs from frozen accepted RB2 state")
    require(isinstance(capture_execution.get("run_id"), str) and capture_execution.get("run_id"),
            "capture execution requires its own run_id")
    require(int(capture_execution.get("iteration", -1)) == 1,
            "capture must be the registered single continuation iteration")

    csv_path = capture_dir / "records.csv"
    manifest_path = capture_dir / "manifest.json"
    if not manifest_path.exists():
        manifest_path = capture_dir / "manifest.final.json"
    require(csv_path.is_file() and manifest_path.is_file(), "capture requires records.csv and manifest")
    manifest = load_json(manifest_path)
    csv_sha, manifest_sha = sha256(csv_path), sha256(manifest_path)
    require(manifest.get("schema") == "rslmto-fixed-u-trace-v1", "unexpected trace manifest schema")
    require(manifest.get("backend") == "RS", "trace manifest backend is not RS")
    require(manifest.get("executable_sha256") == binary_sha, "trace executable differs from diagnostic capture build")
    require(manifest.get("csv_sha256") == csv_sha, "trace manifest does not bind records.csv")
    inputs = manifest.get("input_sha256", {})
    require(capture_input_sha in inputs.values(), "trace manifest does not bind the capture deck")
    require(checkpoint_sha in inputs.values(), "trace manifest does not bind the accepted checkpoint")
    require(int(manifest.get("iteration", -1)) == int(capture_execution["iteration"]),
            "manifest iteration differs from registered capture iteration")
    require(manifest.get("run_id") == capture_execution["run_id"], "manifest run_id differs from capture execution")
    require(manifest.get("source_sha") == capture_build.get("source_sha256"),
            "trace source identifier differs from diagnostic capture source")

    counts: Counter[tuple[str, str, int, int, int]] = Counter()
    run_ids: set[str] = set()
    iterations: set[int] = set()
    row_count = 0
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == CSV_FIELDS, "unexpected capture CSV header")
        for row in reader:
            row_count += 1
            run_ids.add(row["run_id"])
            iterations.add(int(row["iteration"]))
            for key in ("lambda", "weight", "energy", "scalar", "re", "im"):
                require(math.isfinite(float(row[key])), f"non-finite {key} in trace row {row_count}")
            counts[(row["kind"], row["stage"], int(row["site"]), int(row["type"]), int(row["spin"]))] += 1
    require(row_count > 0 and run_ids == {capture_execution["run_id"]}, "capture has empty/mixed/unexpected run IDs")
    require(iterations == {int(capture_execution["iteration"])}, "capture must contain exactly the selected continuation iteration")
    require(int(manifest.get("iteration", -1)) == next(iter(iterations)), "CSV/manifest iteration mismatch")

    # Check the rows needed by the M0.3 carrier and its two registered controls.
    def exactly(kind: str, stage: str, site: int, typ: int, spin: int, n: int = 1) -> None:
        got = counts[(kind, stage, site, typ, spin)]
        require(got == n, f"expected {n} {kind}/{stage} site={site} type={typ} spin={spin} rows, got {got}")

    for spin in (1, 2):
        exactly("m0_carrier_input", "FROZEN_REPLAY", 0, 0, spin)
        exactly("m0_carrier_control", "VOUT_CONTROL", 0, 0, spin)
        exactly("m0_carrier_control", "P_TRANSPOSE_CONTROL", 0, 0, spin)
    kinds = Counter()
    with csv_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            kinds[row["kind"]] += 1
    required = {"radial_replay_mesh", "radial_replay_v", "radial_replay_rho_in", "radial_replay_ec_seed",
                "radial_replay_eigen_seed", "radial_replay_ql", "radial_replay_pl", "pl_map_input",
                "m0_carrier_input", "m0_carrier_control"}
    require(required <= set(kinds), f"capture misses required full replay rows: {sorted(required-set(kinds))}")

    return {
        "schema": "rslmto-rb3-capture-binding-v1",
        "status": "VALIDATED_CAPTURE_BINDING_ONLY_NO_RESPONSE_RUN",
        "rb2": {"run_id": lock["run_id"], "accepted_iteration": int(lock["accepted_iteration"]),
                "backend": "RS", "initial_input_path": str(input_path), "initial_input_sha256": input_sha,
                "initial_seed_path": str(seed_path), "initial_seed_sha256": seed_sha,
                "source_sha256": rb2_source_sha, "source_id": lock["source_id"],
                "executable_sha256": rb2_binary_sha, "convergence": lock["convergence"],
                "c1_acceptance_recomputed": {"outer_diff_le_3e-8": True,
                  "last20_d_occupation_spread_by_spin": occ_spreads,
                  "last20_sphere_moment_spread_muB": moment_spread,
                  "r05_pv_residuals_are_logged_only": True,
                  "r05_residual_record_if_present": convergence.get("fixed_u_residual_tail")}},
        "capture": {"directory": str(capture_dir), "csv_path": str(csv_path), "csv_sha256": csv_sha,
                    "manifest_path": str(manifest_path), "manifest_sha256": manifest_sha,
                    "row_count": row_count, "iteration": int(manifest["iteration"]),
                    "run_id": manifest["run_id"], "trace_source_id": manifest["source_sha"],
                    "kind_counts": dict(sorted(kinds.items())),
                    "capture_build_source_sha256": source_sha, "capture_build_executable_sha256": binary_sha,
                    "diagnostic_only_source_delta_sha256": capture_build.get("diagnostic_only_source_delta_sha256"),
                    "input_sha256": capture_input_sha, "accepted_checkpoint_sha256": checkpoint_sha,
                    "joined_c1_acceptance": lock.get("joined_c1_acceptance")},
        "capture_mode": {"capture_only": True, "ql_only": False},
        "locked_rs_a22_fit": {"path": str(fit), "sha256": fit_sha},
        "lock_path": str(lock_path), "lock_sha256": sha256(lock_path),
        "adapter_path": str(Path(__file__).resolve()), "adapter_sha256": adapter_sha,
        "driver_path": str(driver), "driver_sha256": sha256(driver),
    }


def read_rows(capture: Path) -> list[dict[str, str]]:
    csv_path = capture / "records.csv"
    with csv_path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("capture", type=Path, help="accepted capture directory containing records.csv and manifest")
    ap.add_argument("--lock", type=Path, required=True, help="immutable pre-response accepted RB2 lock JSON")
    ap.add_argument("--source", type=Path, required=True, help="exact source file whose hash is locked (self.f90)")
    ap.add_argument("--binary", type=Path, required=True, help="exact executable whose hash is locked")
    ap.add_argument("--rb2-binary", type=Path, required=True, help="frozen production RB2 executable")
    ap.add_argument("--driver", type=Path, required=True, help="frozen accepted_m03_replay_driver executable")
    ap.add_argument("--fit", type=Path, default=FIT_DEFAULT, help="locked RS A2.2 fit JSON")
    ap.add_argument("--mode", choices=("baseline", "preregister", "cases"), required=True)
    ap.add_argument("--baseline-report", type=Path, help="passing LIVE_STATE baseline report for preregister/cases")
    ap.add_argument("--preregistration", type=Path, help="locked preregistration JSON for cases")
    ap.add_argument("--output", type=Path, required=True, help="new report JSON output")
    args = ap.parse_args()
    binding = validate(args)
    capture = args.capture.resolve()
    rows = read_rows(capture)
    runner_path = Path(__file__).resolve().with_name("rs_a3_m03_runner.py")
    sys.path.insert(0, str(runner_path.parent))
    import rs_a3_m03_runner as runner  # noqa: E402

    runner.FIT_PATH = args.fit.resolve()
    runner.FIT_SHA = FIT_SHA256
    runner._BY_KIND_STAGE = {}
    for row in rows:
        runner._BY_KIND_STAGE.setdefault((row["kind"], row["stage"]), []).append(row)
    nr = sum(1 for r in rows if r["kind"] == "radial_replay_mesh" and int(r["spin"]) == 0)
    ncore = sum(1 for r in rows if r["kind"] == "radial_replay_ec_seed" and r["stage"] == "RADIAL_REPLAY")
    nval = sum(1 for r in rows if r["kind"] == "radial_replay_eigen_seed" and r["stage"] == "RADIAL_REPLAY")
    require((nr, ncore, nval) == (495, 10, 6), f"unsupported replay dimensions {(nr, ncore, nval)}")
    output_dir = args.output.resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)
    # The locked A2.2 evidence stores the fit under per_spin.{fit}; the older
    # RS replay helper expects per_site_spin. Adapt only the JSON shape in a
    # derived file while retaining and validating the immutable A2.2 SHA above.
    fit_data = load_json(args.fit.resolve())
    if "per_site_spin" not in fit_data:
        prereg_path = args.fit.resolve().with_name("a3_rs_a22_scalar_preregistration_v3.json")
        prereg = load_json(prereg_path)
        require(sha256(prereg_path) == fit_data.get("prereg_sha256"),
                "A2.2 fit does not bind the adjacent RS preregistration")
        d0 = float(prereg["accepted_no_u_state"]["D0_RS_ry"])
        per_site = {}
        for spin in ("1", "2"):
            model = fit_data["per_spin"][spin]["fit"]
            per_site[f"site=0,type=0,spin={spin}"] = {
                "r_cubic_coefficients_ascending": model["r_cubic_coefficients_ascending"],
                "a0_ry": model["e_even_fit_coefficients_ascending"][0],
            }
        fit_data = {"schema": "rslmto-rb3-fit-shape-adapter-v1",
                    "source_fit_sha256": FIT_SHA256, "source_prereg_sha256": sha256(prereg_path),
                    "D0_RS_ry": d0, "per_site_spin": per_site}
        fit_adapter_path = output_dir / "a2_2_fit_shape_adapter.json"
        fit_adapter_path.write_text(json.dumps(fit_data, indent=2, sort_keys=True) + "\n")
        runner.FIT_PATH = fit_adapter_path
        runner.FIT_SHA = sha256(fit_adapter_path)
    fit_adapter_sha = sha256(runner.FIT_PATH)
    capture_manifest = binding["capture"]["manifest_sha256"]
    capture_csv = binding["capture"]["csv_sha256"]
    driver = args.driver.resolve()

    if args.mode in ("preregister", "cases"):
        require(args.baseline_report is not None, "preregister/cases requires --baseline-report")
        baseline = load_json(args.baseline_report)
        require(baseline.get("status") == "PASS_LIVE_BASELINE" and
                baseline.get("capture_csv_sha256") == capture_csv and
                baseline.get("capture_manifest_sha256") == capture_manifest and
                baseline.get("driver_sha256") == sha256(driver),
                "baseline is not passing or is not bound to this capture/driver")
    if args.mode == "cases":
        require(args.preregistration is not None, "cases requires --preregistration")
        prereg = load_json(args.preregistration)
        require(prereg.get("status") == "LOCKED_BEFORE_RB3_CASES" and
                prereg.get("capture_csv_sha256") == capture_csv and
                prereg.get("capture_manifest_sha256") == capture_manifest and
                prereg.get("baseline_report_sha256") == sha256(args.baseline_report) and
                prereg.get("driver_sha256") == sha256(driver) and
                prereg.get("adapter_sha256") == sha256(Path(__file__).resolve()) and
                prereg.get("fit_sha256") == FIT_SHA256,
                "RB3 case preregistration does not bind this capture/baseline/driver/adapter/fit")

    if args.mode == "baseline":
        data = runner.prepare(rows, nr, ncore, nval, "baseline")
        raw, input_hash, driver_output_hash = runner.run_persistent(
            driver, runner.text_input(data, rows, nr, ncore, nval), output_dir)
        cases = runner.output_cases(raw)
        require(set(cases) == {1}, "LIVE baseline driver returned unexpected case set")
        checks = runner.baseline_check(rows, cases[1])
        csvout = args.output.with_suffix(".csv")
        csvout.write_bytes(raw)
        report = {"schema": "rslmto-rb3-live-baseline-v1",
                  "status": "PASS_LIVE_BASELINE" if checks["pass"] else "FAIL_LIVE_BASELINE",
                  "capture_csv_sha256": capture_csv, "capture_manifest_sha256": capture_manifest,
                  "driver_sha256": sha256(driver), "adapter_sha256": sha256(Path(__file__).resolve()),
                  "binding_sha256": binding["lock_sha256"], "source_sha256": binding["rb2"]["source_sha256"],
                  "binary_sha256": binding["rb2"]["executable_sha256"], "iteration": binding["capture"]["iteration"],
                  "fit_shape_adapter_sha256": fit_adapter_sha,
                  "checks": checks, "driver_input_sha256": input_hash,
                  "driver_output_path": str(csvout), "driver_output_sha256": sha256(csvout),
                  "driver_output_file_sha256": driver_output_hash}
    elif args.mode == "preregister":
        data = runner.prepare_rs(rows, nr, ncore, nval)
        report = {"schema": "rslmto-rb3-case-preregistration-v1", "status": "LOCKED_BEFORE_RB3_CASES",
                  "capture_csv_sha256": capture_csv, "capture_manifest_sha256": capture_manifest,
                  "baseline_report_sha256": sha256(args.baseline_report), "driver_sha256": sha256(driver),
                  "adapter_sha256": sha256(Path(__file__).resolve()), "fit_sha256": FIT_SHA256,
                  "fit_shape_adapter_sha256": fit_adapter_sha,
                  "source_sha256": binding["rb2"]["source_sha256"], "binary_sha256": binding["rb2"]["executable_sha256"],
                  "iteration": binding["capture"]["iteration"], "case_names": [x[0] for x in data["cases"]],
                  "pre_case_gate": data["preregistration"], "target_operands": data["target"],
                  "locked_rs_a22_fit_path": str(args.fit.resolve()),
                  "notice": "This response preregistration is generated after a passing LIVE baseline and before cases are run."}
    else:
        data = runner.prepare_rs(rows, nr, ncore, nval)
        raw, input_hash, driver_output_hash = runner.run_persistent(
            driver, runner.text_input(data, rows, nr, ncore, nval), output_dir)
        cases = runner.output_cases(raw)
        measurements = runner.cases_summary_rs(rows, data, cases)
        legacy_diagnostics = measurements["gates"]
        ratios = measurements["carrier"]["signed_ratio_by_spin"]
        carrier_ratio_pass = all(0.99 <= float(ratios[s]) <= 1.01 for s in (1, 2))
        q2_pass = float(measurements["q2_local_min_ry2"]) >= 0.0
        scalar = measurements["scalar_plus_001"]
        scalar_pass = float(scalar["absolute_deviation_from_v_sum_q0_ry"]) <= float(scalar["allowance_ry"])
        carrier_error = float(measurements["carrier"]["absolute_error_ry"])
        controls = measurements["wrong_controls"]
        control_separation = {
            name: float(controls[name]["absolute_separation_from_true_carrier_ry"])
            for name in ("SIGN_WRONG", "OPPOSITE_SPIN")
        }
        controls_pass = all(value > 10.0 * carrier_error for value in control_separation.values())
        formal_gates = {
            "carrier_ratio_each_spin_within_1_percent": carrier_ratio_pass,
            "q2_prime_nonnegative": q2_pass,
            "scalar_control_within_locked_a22_allowance": scalar_pass,
            "sign_and_spin_swap_separation_exceeds_10x_measured_carrier_error": controls_pass,
            "control_separation_ry": control_separation,
            "measured_carrier_error_ry": carrier_error,
            "legacy_runner_diagnostic_gates_not_used_as_formal_c1_gates": legacy_diagnostics,
        }
        formal_gates["all_c1_rb3_gates_pass"] = all(
            formal_gates[k] for k in ("carrier_ratio_each_spin_within_1_percent", "q2_prime_nonnegative",
                                      "scalar_control_within_locked_a22_allowance",
                                      "sign_and_spin_swap_separation_exceeds_10x_measured_carrier_error"))
        csvout = args.output.with_suffix(".csv")
        csvout.write_bytes(raw)
        report = {"schema": "rslmto-rb3-measured-cases-v1", "status": "MEASURED_RB3_CASES",
                  "capture_csv_sha256": capture_csv, "capture_manifest_sha256": capture_manifest,
                  "baseline_report_sha256": sha256(args.baseline_report),
                  "preregistration_sha256": sha256(args.preregistration),
                  "driver_sha256": sha256(driver), "adapter_sha256": sha256(Path(__file__).resolve()),
                  "fit_sha256": FIT_SHA256, "source_sha256": binding["rb2"]["source_sha256"],
                  "fit_shape_adapter_sha256": fit_adapter_sha,
                  "binary_sha256": binding["rb2"]["executable_sha256"], "iteration": binding["capture"]["iteration"],
                  "pre_case_gate": data["preregistration"], "measurements": measurements,
                  "c1_rb3_gates": formal_gates,
                  "driver_input_sha256": input_hash, "driver_output_path": str(csvout),
                  "driver_output_sha256": sha256(csvout), "driver_output_file_sha256": driver_output_hash}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps({"mode": args.mode, "status": report["status"], "output": str(args.output),
                      "output_sha256": sha256(args.output), "capture_csv_sha256": capture_csv,
                      "driver_sha256": sha256(driver)}, sort_keys=True))


if __name__ == "__main__":
    main()
