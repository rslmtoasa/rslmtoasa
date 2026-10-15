#!/usr/bin/env python3
"""Hash-bound accepted-iteration PL/NEWRHO clone replay.

The zero-shift output must reproduce archived RAW_PL results before any A1
case is run. The default workflow stops after the repeated baseline and all
eight historical A1 rows pass. Six A2 shifts require a separate explicit
invocation bound to the hash of that passing A1 gate report.
The numerical kernels are the real bands%calculate_pl and unchanged NEWRHO
inside the private Fortran fixture seam.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

CSV_SHA = "55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4"
FIXTURE_SHA = "e88fcff8f1a83831b13834324a6730fa6f8433d7ca3f3b389119d6736b1adbb3"
PL_VALUE_TOL = 1.0e-10
RAW_RHO_TOL = 1.0e-10
RADIAL_ENERGY_TOL_RY = 1.0e-8
A1_SHIFTS = (-2.0e-2, -1.0e-2, -5.0e-3, -2.5e-3,
             2.5e-3, 5.0e-3, 1.0e-2, 2.0e-2)
A2_SHIFTS = (-4.0e-2, -8.0e-2, -1.2e-1, 4.0e-2, 8.0e-2, 1.2e-1)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def f(value: Any) -> float:
    return float(value)


def max_abs(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise ValueError(f"array cardinality mismatch: {len(left)} vs {len(right)}")
    if not left:
        return 0.0
    return max(abs(a - b) for a, b in zip(left, right))


def canonical_sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def verify_artifact(fixture: dict[str, Any], fixture_path: Path) -> tuple[Path, dict[str, Any]]:
    if fixture.get("schema") != "rslmto-accepted-iter23-private-radial-replay-fixture-v1":
        raise ValueError("unsupported accepted replay fixture schema")
    prov = fixture["provenance"]
    if sha(fixture_path) != FIXTURE_SHA:
        raise ValueError("fixture hash differs from the reviewed accepted-input export")
    if int(prov["accepted_iteration"]) != 23 or prov["capture_csv_sha256"] != CSV_SHA:
        raise ValueError("fixture does not bind the accepted iteration-23 CSV")
    capture = Path(prov["capture_dir"]).resolve()
    filehashes = {
        "records.csv": prov["capture_csv_sha256"],
        "manifest.json": prov["capture_manifest_sha256"],
        "Fe.nml": prov["fe_nml_sha256"],
        "Fe_out.nml": prov["accepted_fe_out_sha256"],
        "input.nml": prov["input_nml_postrun_sha256"],
    }
    for name, expected in filehashes.items():
        path = capture / name
        if not path.is_file() or sha(path) != expected:
            raise ValueError(f"hash-bound accepted input changed or missing: {name}")
    manifest = json.loads((capture / "manifest.json").read_text())
    if manifest.get("post_run_sha256", {}).get("records.csv") != CSV_SHA:
        raise ValueError("accepted raw manifest does not bind the approved CSV")
    if (manifest.get("self_f90_sha256") != prov["capture_self_sha256"] or
            manifest.get("binary_sha256") != prov["release_binary_sha256"] or
            manifest.get("expected_accepted_iteration") != 23 or manifest.get("backend") != "KS"):
        raise ValueError("fixture source/binary/backend provenance differs from the raw capture manifest")
    input_hashes = manifest.get("input_sha256", {})
    if (input_hashes.get("input.nml.pre") != prov["input_pre_sha256"] or
            input_hashes.get("Fe.nml.pre") != prov["fe_pre_sha256"]):
        raise ValueError("fixture pre-run input hashes differ from the raw manifest")
    return capture, manifest


def exact_input_text(fixture: dict[str, Any], shifts: tuple[float, ...]) -> str:
    dims, params, radial = fixture["dimensions"], fixture["parameters"], fixture["radial"]
    lmax, nr = int(dims["lmax"]), int(dims["nr"])
    if (lmax, nr, int(params["ncore"]), int(params["nval"])) != (2, 495, 10, 6):
        raise ValueError("captured input dimensions do not match the frozen Fe schema")
    lines = [f"{lmax} {nr} {params['ncore']} {params['nval']} {len(shifts)} {fixture['static_atom']['f_core']}",
             f"{params['z']} {params['a']} {params['b']} {params['tol']}",
             " ".join(format(v, ".17g") for v in shifts),
             " ".join(radial["rofi"])]
    v = radial["v_by_row_spin"]
    for ir in range(nr):
        lines.append(f"{v[2*ir]} {v[2*ir+1]}")
    lines.append(" ".join(radial["ec_seed"]))
    lines.append(" ".join(radial["ev_seed"]))
    for spin in (1, 2):
        for ell in range(lmax + 1):
            operands = fixture["pl_input_by_spin_l"][f"spin_{spin}"][f"l_{ell}"]
            lines.append(" ".join(operands[str(i)] for i in range(5)))
    for spin in (1, 2):
        for ell in range(lmax + 1):
            for moment in fixture["ql_by_spin_l"][f"spin_{spin}"][f"l_{ell}"]:
                lines.append(moment)
    return "\n".join(lines) + "\n"


def run_driver(driver: Path, fixture: dict[str, Any], shifts: tuple[float, ...]) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="accepted-pl-newrho-") as temp:
        work = Path(temp)
        inp, out = work / "captured_state.in", work / "replay.csv"
        inp.write_text(exact_input_text(fixture, shifts))
        proc = subprocess.run([str(driver), str(inp), str(out)], check=False,
                              text=True, capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(f"private replay driver failed ({proc.returncode}): {proc.stderr[-1000:]}")
        records: dict[int, dict[str, Any]] = {}
        with out.open(newline="") as stream:
            for row in csv.reader(stream):
                kind, iv = row[0], int(row[1])
                case = records.setdefault(iv, {"pl": {}, "rho": {}, "sumec": {}, "sumev": {}, "ec": {}, "ev": {}})
                if kind == "PL":
                    case["pl"][(int(row[2]), int(row[3]))] = float(row[4])
                elif kind == "RHO":
                    case["rho"][(int(row[2])-1, int(row[3]))] = float(row[4])
                elif kind in ("SUMEC", "SUMEV"):
                    case[kind.lower()][int(row[2])] = float(row[3])
                elif kind in ("EC", "EV"):
                    case[kind.lower()][int(row[2])-1] = float(row[3])
                else:
                    raise ValueError(f"unexpected private replay output row kind {kind}")
        if set(records) != set(range(1, len(shifts) + 1)):
            raise ValueError("private replay output has missing or extra cases")
        for iv, case in records.items():
            if (len(case["pl"]) != 6 or len(case["rho"]) != 990 or len(case["sumec"]) != 2 or
                    len(case["sumev"]) != 2 or len(case["ec"]) != 10 or len(case["ev"]) != 6):
                raise ValueError(f"private replay output cardinality failure for case {iv}")
        return {"cases": records, "input_sha256": sha(inp), "output_sha256": sha(out)}


def compare_baseline(case: dict[str, Any], fixture: dict[str, Any]) -> dict[str, Any]:
    raw_pl = fixture["raw_pl_return_by_spin"]
    ref_pl = [f(raw_pl[f"spin_{s}"][ell]) for ell in range(3) for s in (1, 2)]
    got_pl = [case["pl"][(ell, s)] for ell in range(3) for s in (1, 2)]
    ref_rho = [f(x) for x in fixture["raw_rho_output_by_row_spin"]]
    got_rho = [case["rho"][(ir, s)] for ir in range(495) for s in (1, 2)]
    raw_energy = fixture["raw_neWRHO_reference"]
    ref_sumec = [f(raw_energy["ec_output_by_spin"][f"spin_{s}"]) for s in (1, 2)]
    ref_sumev = [f(raw_energy["ev_output_by_spin"][f"spin_{s}"]) for s in (1, 2)]
    got_sumec = [case["sumec"][s] for s in (1, 2)]
    got_sumev = [case["sumev"][s] for s in (1, 2)]
    raw_ev = sorted(fixture["raw_valence_ev_reference_by_state"], key=lambda r: int(r["band"]))
    ref_ev = [f(r["energy"]) for r in raw_ev]
    got_ev = [case["ev"][i] for i in range(6)]
    checks = {
        "calculate_pl_dimensionless": {"max_abs": max_abs(got_pl, ref_pl), "tol": PL_VALUE_TOL},
        "rho": {"max_abs": max_abs(got_rho, ref_rho), "tol": RAW_RHO_TOL},
        "SUMEC": {"max_abs": max_abs(got_sumec, ref_sumec), "tol": RADIAL_ENERGY_TOL_RY},
        "SUMEV": {"max_abs": max_abs(got_sumev, ref_sumev), "tol": RADIAL_ENERGY_TOL_RY},
        "EV_vector": {"max_abs": max_abs(got_ev, ref_ev), "tol": RADIAL_ENERGY_TOL_RY},
    }
    checks["pass"] = all(x["max_abs"] <= x["tol"] for x in checks.values())
    return checks


def compare_a1(cases: dict[int, dict[str, Any]], fixture: dict[str, Any]) -> dict[str, Any]:
    summaries: dict[str, Any] = {}
    historical = fixture["historical_a1_case_rows_by_stage_spin"]
    for index, shift in enumerate(A1_SHIFTS, start=1):
        iv = index + 1  # case 1 in this invocation is the repeated zero baseline
        stage = f"A1_SHIFT_{index:02d}"
        ref_rows = historical[stage]
        spin_errors = {}
        pl_errors = {}
        spin_actual = {}
        for spin in (1, 2):
            expected = f(ref_rows[spin-1]["scalar"])
            actual = cases[iv]["sumev"][spin]
            spin_errors[f"spin_{spin}"] = abs(actual - expected)
            spin_actual[f"spin_{spin}"] = actual
            ref_pl = fixture["historical_a1_pl_return_by_stage_spin"][stage][f"spin_{spin}"]
            pl_errors[f"spin_{spin}"] = max_abs(
                [cases[iv]["pl"][(ell, spin)] for ell in range(3)], [f(x) for x in ref_pl])
        summaries[stage] = {"v": shift, "spin_max_abs_sumev_error_ry": spin_errors,
                            "max_abs_calculate_pl_error": pl_errors,
                            "calculate_pl_tol_dimensionless": PL_VALUE_TOL,
                            "actual_sumev_by_spin_ry": spin_actual,
                            "max_abs_sumev_error_ry": max(spin_errors.values()),
                            "max_abs_calculate_pl_error_dimensionless": max(pl_errors.values()),
                            "sumev_tol_ry": RADIAL_ENERGY_TOL_RY,
                            "pass": (max(spin_errors.values()) <= RADIAL_ENERGY_TOL_RY and
                                    max(pl_errors.values()) <= PL_VALUE_TOL)}
    return {"cases": summaries, "pass": all(c["pass"] for c in summaries.values())}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fixture", type=Path, required=True)
    ap.add_argument("--driver", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--through-a1", action="store_true",
                      help="run baseline and historical A1 reproduction only; never run A2 shifts")
    mode.add_argument("--a2-only", action="store_true",
                      help="run six A2 shifts only after validating a hash-bound passing A1 report")
    ap.add_argument("--gate-report", type=Path,
                    help="required with --a2-only; output of a passing --through-a1 invocation")
    args = ap.parse_args()
    if args.a2_only and args.gate_report is None:
        ap.error("--a2-only requires --gate-report")
    if args.through_a1 and args.gate_report is not None:
        ap.error("--gate-report is only valid with --a2-only")
    fixture = json.loads(args.fixture.read_text())
    capture, raw_manifest = verify_artifact(fixture, args.fixture)
    driver = args.driver.resolve()
    if not driver.is_file():
        raise ValueError("private Fortran replay driver is missing")

    gate_report_sha = None
    if args.a2_only:
        gate_report = json.loads(args.gate_report.read_text())
        gate_provenance = gate_report.get("provenance", {})
        if (gate_report.get("schema") != "rslmto-accepted-iter23-private-clone-replay-v1" or
                gate_report.get("status") != "PASS_A1_GATE_STOP_BEFORE_A2" or
                gate_report.get("baseline", {}).get("checks", {}).get("pass") is not True or
                gate_report.get("a1_validation", {}).get("checks", {}).get("pass") is not True or
                gate_provenance.get("fixture_sha256") != sha(args.fixture) or
                gate_provenance.get("replay_driver_sha256") != sha(driver) or
                gate_provenance.get("replay_source_sha256") != sha(Path(__file__).resolve().parents[3] / "source" / "self.f90")):
            raise ValueError("A2 gate report is not a passing A1 report for these exact fixture/driver/source hashes")
        gate_report_sha = sha(args.gate_report)

    result: dict[str, Any] = {
        "schema": "rslmto-accepted-iter23-private-clone-replay-v1",
        "diagnostic": "A2.2 accepted-state replay; no SCF run",
        "provenance": {"fixture_path": str(args.fixture.resolve()), "fixture_sha256": sha(args.fixture),
            "capture_path": str(capture), "capture_csv_sha256": fixture["provenance"]["capture_csv_sha256"],
            "capture_manifest_sha256": fixture["provenance"]["capture_manifest_sha256"],
            "capture_source_sha256": fixture["provenance"]["capture_self_sha256"],
            "capture_binary_sha256": fixture["provenance"]["release_binary_sha256"],
            "replay_driver_sha256": sha(driver),
            "replay_source_sha256": sha(Path(__file__).resolve().parents[3] / "source" / "self.f90"),
            "accepted_iteration": 23, "backend": raw_manifest.get("backend")},
        "predeclared_tolerances": {"calculate_pl_dimensionless": PL_VALUE_TOL,
            "SUMEC_SUMEV_EV_vector_ry": RADIAL_ENERGY_TOL_RY,
            "rho_max_abs": RAW_RHO_TOL,
            "basis": {"calculate_pl_dimensionless": "existing m0_a11_scaling_oracle.py STATE_VALUE_TOL",
                "radial_energies_and_density": "existing local_moment_capture_oracle.py production replay limits"}},
        "historical_reference_limitations": {
            "individual_returned_EC_vector": "not serialized; exact EC seeds and per-spin SUMEC are available",
            "same_call_neWRHO_errors": fixture.get("raw_live_neWRHO_verification_error", {})},
        "frozen_radial_operand_sha256": canonical_sha({
            "radial": fixture["radial"], "ql": fixture["ql_by_spin_l"],
            "parameters": fixture["parameters"], "static_atom": fixture["static_atom"]}),
        "cases": {},
    }

    if args.a2_only:
        # The previously reviewed gate report is the sole authorization to
        # expose the six additional outcomes. Do not rerun baseline/A1 here.
        extra_run = run_driver(driver, fixture, A2_SHIFTS)
        for i, shift in enumerate(A2_SHIFTS, start=9):
            stage = f"A2_SHIFT_{i:02d}"
            case = extra_run["cases"][i - 8]
            result["cases"][stage] = {"v": shift, "pl": {f"l{ell}_spin{s}": case["pl"][(ell, s)]
                for ell in range(3) for s in (1, 2)}, "sumec": {f"spin_{s}": case["sumec"][s] for s in (1, 2)},
                "sumev": {f"spin_{s}": case["sumev"][s] for s in (1, 2)},
                "rho": [[case["rho"][(ir, s)] for s in (1, 2)] for ir in range(495)],
                "ec_sha256": hashlib.sha256("\n".join(format(case["ec"][j], ".17g") for j in range(10)).encode()).hexdigest(),
                "ec_vector": [case["ec"][j] for j in range(10)], "ev": [case["ev"][j] for j in range(6)]}
        result["a2_gate"] = {"report_path": str(args.gate_report.resolve()),
                              "report_sha256": gate_report_sha,
                              "status": "PASS_A1_GATE_STOP_BEFORE_A2"}
        result["a2_extra_run"] = {"shift_vector_v": list(A2_SHIFTS),
            "stage_mapping": [f"A2_SHIFT_{i:02d}" for i in range(9, 15)],
            "input_sha256": extra_run["input_sha256"], "output_sha256": extra_run["output_sha256"]}
        result["status"] = "A2_EXTRA_CASES_CAPTURED_AFTER_HASH_BOUND_A1_GATE"
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"status": result["status"], "output": str(args.output.resolve()),
                          "output_sha256": sha(args.output), "A1_gate_report_sha256": gate_report_sha}, sort_keys=True))
        return 0

    baseline_run = run_driver(driver, fixture, (0.0,))
    baseline = baseline_run["cases"][1]
    baseline_checks = compare_baseline(baseline, fixture)
    result["baseline"] = {"v": 0.0, "checks": baseline_checks,
                          "input_sha256": baseline_run["input_sha256"],
                          "output_sha256": baseline_run["output_sha256"],
                          "ec_output_vector": [baseline["ec"][i] for i in range(10)],
                          "ec_vector_output_uncompared": True,
                          "ev_output_vector": [baseline["ev"][i] for i in range(6)]}
    if not baseline_checks["pass"]:
        result["status"] = "FAIL_RAW_BASELINE_STOP"
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"status": result["status"], "output": str(args.output.resolve()),
                          "output_sha256": sha(args.output)}, sort_keys=True))
        return 1

    old_run = run_driver(driver, fixture, (0.0, *A1_SHIFTS))
    repeated_baseline = compare_baseline(old_run["cases"][1], fixture)
    old_checks = compare_a1(old_run["cases"], fixture)
    result["a1_validation"] = {"repeated_zero_baseline": repeated_baseline,
        "stage_shift_mapping": {f"A1_SHIFT_{i:02d}": A1_SHIFTS[i-1] for i in range(1, 9)},
        "checks": old_checks, "input_sha256": old_run["input_sha256"],
        "output_sha256": old_run["output_sha256"]}
    if not repeated_baseline["pass"] or not old_checks["pass"]:
        result["status"] = "FAIL_A1_REPRODUCTION_STOP_A2"
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"status": result["status"], "output": str(args.output.resolve()),
                          "output_sha256": sha(args.output)}, sort_keys=True))
        return 1
    result["status"] = "PASS_A1_GATE_STOP_BEFORE_A2"
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "output": str(args.output.resolve()),
                      "output_sha256": sha(args.output), "A1_validation": "PASS", "A2": "NOT_RUN"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
