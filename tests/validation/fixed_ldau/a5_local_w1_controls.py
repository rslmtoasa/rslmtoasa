#!/usr/bin/env python3
"""Run the preregistered A5.2 local-W1 KS controls from one common seed."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
LAMBDAS = (-0.005, 0.0, 0.005)
CONTROLS = {
    "base_300K_12cubed": {"hoh": ".false.", "temperature": 300.0, "mesh": 12},
    "hoh_300K_12cubed": {"hoh": ".true.", "temperature": 300.0, "mesh": 12,
                         "fixed_u_test_allow_hoh_w1_u0": ".true."},
    "temperature_150K_12cubed": {"hoh": ".false.", "temperature": 150.0, "mesh": 12},
    "temperature_600K_12cubed": {"hoh": ".false.", "temperature": 600.0, "mesh": 12},
    "mesh_20cubed_300K": {"hoh": ".false.", "temperature": 300.0, "mesh": 20},
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def edit_group(deck: str, group: str, assignments: dict[str, str]) -> str:
    lines = deck.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip().lower() == f"&{group}"), None)
    if start is None:
        lines.extend([f"&{group}", "/"])
        start = len(lines) - 2
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == "/"), None)
    if end is None:
        raise ValueError(f"unterminated &{group} namelist")
    remaining = dict(assignments)
    kept = []
    for line in lines[start + 1:end]:
        stripped = line.lstrip()
        key = stripped.split("=", 1)[0].strip().lower() if "=" in stripped and not stripped.startswith("!") else ""
        if key in remaining:
            kept.append(f"{key} = {remaining.pop(key)}")
        else:
            kept.append(line)
    kept.extend(f"{key} = {value}" for key, value in remaining.items())
    lines[start + 1:end] = kept
    return "\n".join(lines) + "\n"


def parse_last(log: str, label: str) -> float:
    matches = re.findall(rf"{re.escape(label)}=\s*([^\s]+)", log)
    if not matches or matches[-1] in {"UNAVAILABLE", "NaN", "nan"}:
        raise ValueError(f"missing finite {label} in final runtime log")
    value = float(matches[-1])
    if not (value == value and abs(value) != float("inf")):
        raise ValueError(f"nonfinite {label}")
    return value


def run_case(binary: Path, base_input: Path, base_seed: Path, assets: Path,
             output: Path, control_name: str, cfg: dict, lam: float) -> dict:
    run_id = f"{control_name}_lambda_{lam:+.3f}"
    case_dir = output / run_id
    case_dir.mkdir(parents=True, exist_ok=False)
    deck = edit_group(base_input.read_text(), "self", {})
    deck = edit_group(deck, "hamiltonian", {
        "hoh": cfg["hoh"], "hubbard_u_general(1,3)": "0.0", "hubbard_j_general(1,3)": "0.0",
        "hubbard_u_potential_form": "'liechtenstein'", "hubbard_radial_moments": "'local'",
        "hubbard_test_operator": "'w1'", "hubbard_test_lambda": f"{lam:.17g}",
        "hubbard_test_spin": "0",
    } | ({"fixed_u_test_allow_hoh_w1_u0": cfg["fixed_u_test_allow_hoh_w1_u0"]}
         if "fixed_u_test_allow_hoh_w1_u0" in cfg else {}))
    deck = edit_group(deck, "reciprocal", {
        "temperature": f"{cfg['temperature']:.17g}", "nk1": str(cfg["mesh"]),
        "nk2": str(cfg["mesh"]), "nk3": str(cfg["mesh"]),
    })
    input_path, seed_path = case_dir / "input.nml", case_dir / "Fe.nml"
    input_path.write_text(deck)
    shutil.copyfile(base_seed, seed_path)
    input_sha, seed_sha = sha(input_path), sha(seed_path)
    for asset in ASSETS:
        shutil.copyfile(assets / asset, case_dir / asset)
    trace = case_dir / "moments_trace.csv"
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "RSLMTO_FIXED_U_TRACE": str(trace),
                "RSLMTO_FIXED_U_TRACE_ITERATION": "0", "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
                "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1", "RSLMTO_FIXED_U_LAMBDA": f"{lam:.17g}",
                "RSLMTO_FIXED_U_RUN_ID": run_id})
    stdout = case_dir / "stdout.log"
    # Stream output to disk while the SCF child is running so an interrupted
    # campaign still preserves the full residual history.
    with stdout.open("wb") as stdout_file:
        result = subprocess.run([str(binary)], cwd=case_dir, env=env, stdin=subprocess.DEVNULL,
                                stdout=stdout_file, stderr=subprocess.STDOUT, check=False)
    log = stdout.read_text(errors="replace")
    if (case_dir / "testrun.log").is_file():
        log += "\n" + (case_dir / "testrun.log").read_text(errors="replace")
    if result.returncode != 0 or "Converged!" not in log or not trace.is_file():
        failure = {"run_id": run_id, "control": control_name, "lambda_ry": lam,
            "returncode": result.returncode, "converged_marker": "Converged!" in log,
            "input_pre_sha256": input_sha, "seed_pre_sha256": seed_sha,
            "runtime_Fe_nml_sha256": sha(seed_path) if seed_path.is_file() else None,
            "runtime_restart_sha256": {name: sha(case_dir / name) for name in
                ("Fe_out.nml", "Fe_scf.nml", "Fe_scf_beta.nml") if (case_dir / name).is_file()},
            "trace_sha256": sha(trace) if trace.is_file() else None,
            "stdout_sha256": sha(stdout), "trace_bytes": trace.stat().st_size if trace.is_file() else None,
            "failure_reason": "returncode nonzero, no Converged! marker, or missing trace; preserve this raw case and stop campaign"}
        (case_dir / "failure.json").write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
        raise RuntimeError(f"{run_id}: SCF did not exit successfully and converge; see {case_dir}")
    with trace.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"test_density_matrix", "test_q0_trace", "test_trace_wp"}
    kinds = {row["kind"] for row in rows}
    if not required <= kinds:
        raise RuntimeError(f"{run_id}: missing trace rows {sorted(required-kinds)}")
    f_local = parse_last(log, "test_F_local")
    return {"run_id": run_id, "control": control_name, "lambda_ry": lam,
            "temperature_k": cfg["temperature"], "mesh": cfg["mesh"], "hoh": cfg["hoh"],
            "returncode": result.returncode, "converged": True, "F_local_ry": f_local,
            "P_OUT_test_trace_wp_by_spin": {
                spin: sum(float(row["scalar"]) for row in rows
                          if row["kind"] == "test_trace_wp" and row["stage"] == "P_OUT" and row["spin"] == spin)
                for spin in ("0", "1")},
            "input_sha256": input_sha, "seed_input_sha256": seed_sha,
            "runtime_Fe_nml_sha256": sha(seed_path),
            "runtime_Fe_out_nml_sha256": sha(case_dir / "Fe_out.nml") if (case_dir / "Fe_out.nml").is_file() else None,
            "trace_output_sha256": sha(trace), "stdout_sha256": sha(stdout),
            "trace_bytes": trace.stat().st_size}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.binary, args.input, args.seed, *(args.assets / asset for asset in ASSETS)):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=False)
    results = []
    for name, cfg in CONTROLS.items():
        for lam in LAMBDAS:
            results.append(run_case(args.binary, args.input, args.seed, args.assets,
                                    args.output, name, cfg, lam))
    grouped = {name: {str(row["lambda_ry"]): row for row in results if row["control"] == name}
               for name in CONTROLS}
    for name, cases in grouped.items():
        plus, zero, minus = (cases[str(x)] for x in (0.005, 0.0, -0.005))
        fd = (plus["F_local_ry"] - minus["F_local_ry"]) / 0.01
        direct = sum(zero["P_OUT_test_trace_wp_by_spin"].values())
        residual = fd - direct
        cases["summary"] = {"central_dF_dlambda": fd, "direct_Tr_WP_at_zero": direct,
                            "signed_residual": residual,
                            "relative_residual": residual / direct if direct else None}
    report = {"status": "MEASURED_A5_2_LOCAL_W1_CONTROLS", "initial_state_label": "shared immutable near-converged restart; not an accepted/converged endpoint",
              "binary_sha256": sha(args.binary),
              "base_input_sha256": sha(args.input), "base_seed_sha256": sha(args.seed),
              "runs": results, "controls": grouped}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "runs": len(results),
                      "summary": str(args.output / "summary.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
