#!/usr/bin/env python3
"""Run and reduce the preregistered spin-selective W1 response probe."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
SOURCES = (
    "source/hamiltonian.f90",
    "source/hamiltonian_build.f90",
    "source/hamiltonian_hubbard.f90",
    "source/include_codes/namelists/hamiltonian.f90",
    "source/self.f90",
)
EPS = 2.220446049250313e-16
FLOAT_TOKEN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    kept: list[str] = []
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


def number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def has_converged_marker(log: str) -> bool:
    return re.search(r"(?<!Not )Converged!", log, re.IGNORECASE) is not None


def parse_logs(case_dir: Path, stdout: str, require_test_energy: bool = True) -> tuple[float | None, float | None, bool]:
    log = stdout
    testrun = case_dir / "testrun.log"
    if testrun.exists():
        log += "\n" + testrun.read_text(errors="replace")
    f_values = re.findall(r"test_F_local\s*=\s*(" + FLOAT_TOKEN.pattern + r")", log)
    e_values = re.findall(r"test_lambda_trace_wp\s*=\s*(" + FLOAT_TOKEN.pattern + r")", log)
    if require_test_energy and (not f_values or not e_values):
        raise RuntimeError(f"missing local free energy or external trace in {case_dir}")
    converged = has_converged_marker(log)
    return (number(f_values[-1]) if f_values else None,
            number(e_values[-1]) if e_values else None, converged)


def parse_final_thermodynamics(case_dir: Path, stdout: str) -> tuple[float, float]:
    """Read final physical ETOT and TS for reconstructing external-field F."""
    log = stdout
    testrun = case_dir / "testrun.log"
    if testrun.exists():
        log += "\n" + testrun.read_text(errors="replace")
    etot_tokens = re.findall(r"physical_total_energy\s*=\s*(" + FLOAT_TOKEN.pattern + r")", log)
    ts_tokens = re.findall(r"\bTS\s*=\s*(" + FLOAT_TOKEN.pattern + r")", log)
    if not etot_tokens or not ts_tokens:
        raise RuntimeError(f"missing finite ETOT or TS in {case_dir}")
    return number(etot_tokens[-1]), number(ts_tokens[-1])


def reconstruct_local_free_energy(etot: float, lam: float, tr_wp: float, ts: float) -> float:
    """Use ETOT + lambda*Tr(WP) - TS, all taken from this same final state."""
    return etot + lam * tr_wp - ts


def parse_trace(path: Path) -> dict:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows = [r for r in rows if r["stage"] == "P_OUT" and int(r["l"]) == 2]
    qrows = [r for r in rows if r["kind"] == "test_q0_trace"]
    wrows = [r for r in rows if r["kind"] == "test_trace_wp"]
    prows = [r for r in rows if r["kind"] == "test_density_matrix"]
    if not qrows or not wrows or not prows:
        raise RuntimeError(f"missing P_OUT spin moments in {path}")
    iteration = max(int(r["iteration"]) for r in qrows)
    qrows = [r for r in qrows if int(r["iteration"]) == iteration]
    wrows = [r for r in wrows if int(r["iteration"]) == iteration]
    prows = [r for r in prows if int(r["iteration"]) == iteration]
    q = [0.0, 0.0]
    w = [0.0, 0.0]
    p: list[dict[tuple[int, int, int], complex]] = [{}, {}]
    for row in qrows:
        spin = int(row["spin"])
        q[spin] += number(row["scalar"])
    for row in wrows:
        spin = int(row["spin"])
        w[spin] += number(row["scalar"])
    for row in prows:
        spin, i, j = int(row["spin"]), int(row["row"]), int(row["col"])
        p[spin][(int(row["site"]), i, j)] = complex(number(row["re"]), number(row["im"]))
    p_trace = [sum(value.real for (site, i, j), value in block.items() if i == j) for block in p]
    return {"iteration": iteration, "q": q, "w": w, "p_trace": p_trace,
            "p": [{f"{site}:{i}:{j}": [value.real, value.imag] for (site, i, j), value in block.items()}
                  for block in p]}


def case_specs() -> list[dict]:
    cases = [{"name": "base_h0", "lambda": 0.0, "spin": 0, "conv_thr": 1.0e-9}]
    for spin, label in ((1, "up"), (2, "down")):
        for h in (0.005, 0.0025):
            for sign in (-1, 1):
                cases.append({"name": f"{label}_{sign:+d}_{h:g}", "lambda": sign*h,
                              "spin": spin, "conv_thr": 1.0e-9})
    for spin, label in ((1, "up"), (2, "down")):
        for sign in (-1, 1):
            cases.append({"name": f"{label}_{sign:+d}_0.005_tight", "lambda": sign*0.005,
                          "spin": spin, "conv_thr": 1.0e-10})
    cases.append({"name": "base_h0_tight", "lambda": 0.0, "spin": 0, "conv_thr": 1.0e-10})
    return cases


def configure_case_deck(base_text: str, spec: dict) -> str:
    self_assignments = {"use_kspace": ".true."}
    if spec["conv_thr"] != 1.0e-9:
        self_assignments["conv_thr"] = "1.0d-10"
    deck = edit_group(base_text, "self", self_assignments)
    ham = {
        "hoh": ".false.",
        "hubbard_u_general(1,3)": "0.0",
        "hubbard_j_general(1,3)": "0.0",
        "hubbard_u_potential_form": "'liechtenstein'",
        "hubbard_radial_moments": "'local'",
        "hubbard_test_operator": "'w1'",
        "hubbard_test_lambda": f"{spec['lambda']:.16g}",
        "hubbard_test_spin": str(spec["spin"]),
    }
    return edit_group(deck, "hamiltonian", ham)


def run_case(spec: dict, binary: Path, base_input: Path, base_seed: Path,
             assets_dir: Path, out_dir: Path, source_root: Path) -> dict:
    case_dir = out_dir / spec["name"]
    case_dir.mkdir()
    # The historical archived input defaults to recursion. Explicitly select
    # reciprocal KS in each generated case rather than inheriting that default.
    deck = configure_case_deck(base_input.read_text(), spec)
    input_path = case_dir / "input.nml"
    # The accepted KS restart is a prior result; the executable loads it from
    # Fe_out.nml. Keep this pre-run input hash separate from the mutated
    # post-run restart hash.
    seed_path = case_dir / "Fe_out.nml"
    input_path.write_text(deck)
    shutil.copyfile(base_seed, seed_path)
    input_sha, seed_input_sha = sha(input_path), sha(seed_path)
    for asset in ASSETS:
        shutil.copyfile(assets_dir / asset, case_dir / asset)

    trace_path = case_dir / "moments.csv"
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "RSLMTO_FIXED_U_TRACE": str(trace_path),
        "RSLMTO_FIXED_U_TRACE_ITERATION": "0",
        "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1",
        "RSLMTO_FIXED_U_LAMBDA": f"{spec['lambda']:.16g}",
        "RSLMTO_FIXED_U_RUN_ID": spec["name"],
    })
    stdout_path = case_dir / "stdout.log"
    # Stream output directly to disk so a long SCF attempt remains auditable
    # even if the Python runner is interrupted while the child is running.
    with stdout_path.open("wb") as stdout_file:
        result = subprocess.run([str(binary)], cwd=case_dir, env=env, stdin=subprocess.DEVNULL,
                                stdout=stdout_file, stderr=subprocess.STDOUT, check=False)
    stdout = stdout_path.read_text(errors="replace")
    log = stdout
    testrun = case_dir / "testrun.log"
    if testrun.exists():
        log += "\n" + testrun.read_text(errors="replace")
    ks_branch = "run_recursion: use_kspace=.true., skipping recursion solver stage" in log
    thermo_pairs = re.findall(r"entropy=\s*([^\s]+).*?TS=\s*([^\s]+)", log)
    number_pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
    finite_thermo = bool(thermo_pairs) and all(
        re.fullmatch(number_pattern, value) for value in thermo_pairs[-1])
    if not ks_branch or not finite_thermo:
        failure = {"case": spec, "returncode": result.returncode,
                   "ks_branch": ks_branch, "finite_entropy_TS": finite_thermo,
                   "input_sha256": input_sha, "seed_input_sha256": seed_input_sha,
                   "stdout_sha256": sha(stdout_path),
                   "trace_sha256": sha(trace_path) if trace_path.is_file() else None}
        (case_dir / "setup_failure.json").write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
        raise RuntimeError(f"{spec['name']} lacks native-KS/finite-thermo validation; raw case preserved")
    local_f_logged, external_logged, converged = parse_logs(
        case_dir, stdout, require_test_energy=False)
    trace = parse_trace(trace_path)
    if result.returncode != 0 or not converged:
        raise RuntimeError(f"{spec['name']} did not converge cleanly (rc={result.returncode})")
    lam = spec["lambda"]
    active_spin = spec["spin"] - 1
    expected_external = lam * trace["w"][active_spin] if spec["spin"] else lam * sum(trace["w"])
    etot, ts = parse_final_thermodynamics(case_dir, stdout)
    local_f = reconstruct_local_free_energy(etot, lam,
                                            trace["w"][active_spin] if spec["spin"] else sum(trace["w"]), ts)
    row = {
        **spec, "returncode": result.returncode, "converged_marker": converged,
        "ks_branch": ks_branch, "finite_entropy_TS": finite_thermo,
        "input_sha256": input_sha, "seed_input_sha256": seed_input_sha,
        "stdout_sha256": sha(stdout_path), "trace_sha256": sha(trace_path),
        "postrun_Fe_out_nml_sha256": sha(seed_path),
        "ETOT_Ry": etot, "TS_Ry": ts,
        "F_local_Ry": local_f, "F_local_from_runtime_logger_Ry": local_f_logged,
        "F_local_logger_minus_reconstructed_Ry": (local_f_logged - local_f) if local_f_logged is not None else None,
        "external_lambda_trace_logged_Ry": external_logged,
        "external_lambda_trace_from_P_Ry": expected_external,
        "external_trace_residual_Ry": external - expected_external if external is not None else None,
        **trace,
    }
    (case_dir / "result.json").write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    return row


def invert2(a: list[list[float]]) -> list[list[float]]:
    det = a[0][0]*a[1][1] - a[0][1]*a[1][0]
    if det == 0.0:
        raise ArithmeticError("exactly singular response matrix")
    return [[a[1][1]/det, -a[0][1]/det], [-a[1][0]/det, a[0][0]/det]]


def mv(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j]*x[j] for j in range(len(x))) for i in range(len(a))]


def derive(rows: list[dict], u_response_path: Path | None) -> dict:
    cases = {r["name"]: r for r in rows}
    q0 = cases["base_h0"]["q"]
    q0_tight = cases["base_h0_tight"]["q"]
    b_by_h: dict[float, list[float]] = {}
    r_by_h: dict[float, list[list[float]]] = {}
    for h in (0.005, 0.0025):
        b = []
        rmat = [[0.0, 0.0], [0.0, 0.0]]
        for col, label in enumerate(("up", "down")):
            plus, minus = cases[f"{label}_+1_{h:g}"], cases[f"{label}_-1_{h:g}"]
            d_f = (plus["F_local_Ry"] - minus["F_local_Ry"])/(2*h)
            b.append(d_f - q0[col])
            for s in range(2):
                rmat[s][col] = (plus["q"][s] - minus["q"][s])/(2*h)
        b_by_h[h], r_by_h[h] = b, rmat
    b_primary, r_primary = b_by_h[0.0025], r_by_h[0.0025]
    b_step = [abs(b_by_h[0.005][i] - b_primary[i]) for i in range(2)]
    r_step = [[abs(r_by_h[0.005][i][j] - r_primary[i][j]) for j in range(2)] for i in range(2)]

    b_repeat = [0.0, 0.0]
    r_repeat = [[0.0, 0.0], [0.0, 0.0]]
    for col, label in enumerate(("up", "down")):
        plus, minus = cases[f"{label}_+1_0.005"], cases[f"{label}_-1_0.005"]
        plus_t, minus_t = cases[f"{label}_+1_0.005_tight"], cases[f"{label}_-1_0.005_tight"]
        d_f_std = (plus["F_local_Ry"] - minus["F_local_Ry"])/0.01
        d_f_tight = (plus_t["F_local_Ry"] - minus_t["F_local_Ry"])/0.01
        b_repeat[col] = abs(d_f_tight - d_f_std) + abs(q0_tight[col] - q0[col])
        for s in range(2):
            r_std = (plus["q"][s] - minus["q"][s])/0.01
            r_tight = (plus_t["q"][s] - minus_t["q"][s])/0.01
            r_repeat[s][col] = abs(r_tight - r_std)
    b_round = []
    for i, label in enumerate(("up", "down")):
        plus = cases[f"{label}_+1_0.0025"]["F_local_Ry"]
        minus = cases[f"{label}_-1_0.0025"]["F_local_Ry"]
        b_round.append(16*EPS*(abs(plus) + abs(minus))/(2*0.0025) +
                       16*EPS*abs(q0[i]))
    # Conservative sum of step sensitivity, tightened-repeat difference, and
    # output-rounding scale. The value uncertainty is not used as a gate.
    u_b = [b_step[i] + b_repeat[i] + b_round[i] for i in range(2)]
    u_r = [[r_step[i][j] + r_repeat[i][j] +
            16*EPS*(abs(r_primary[i][j]) + 1.0) for j in range(2)] for i in range(2)]

    result = {
        "q0_up_down": q0, "q0_tight_up_down": q0_tight,
        "B_by_h_Ry_per_Ry": {str(h): b_by_h[h] for h in b_by_h},
        "R_by_h_per_Ry": {str(h): r_by_h[h] for h in r_by_h},
        "B_primary": b_primary, "R_primary": r_primary,
        "B_uncertainty": u_b, "R_uncertainty": u_r,
        "B_step_sensitivity": b_step, "R_step_sensitivity": r_step,
        "B_repeat_sensitivity": b_repeat, "R_repeat_sensitivity": r_repeat,
        "B_sum_primary": sum(b_primary),
        "B_sum_uncertainty": sum(u_b),
        "B_sum_target_A5_1_local": -0.009214878893257783,
        "B_sum_A5_1_uncertainty_bound": 5.0e-5,
        "B_sum_difference_from_A5_1": sum(b_primary) + 0.009214878893257783,
        "B_sum_A5_1_comparison_role": "historical definition diagnostic only; no agreement gate",
        "external_trace_residual_max_Ry": max(abs(r["external_trace_residual_Ry"]) for r in rows),
    }
    # Exact singularity is reported, with no arbitrary condition-number gate.
    try:
        a = [[r_primary[0][0], r_primary[1][0]], [r_primary[0][1], r_primary[1][1]]]
        a_inv = invert2(a)
        g = mv(a_inv, b_primary)
        u_g = []
        for i in range(2):
            u_g.append(sum(abs(a_inv[i][j])*(u_b[j] + sum(u_r[k][j]*abs(g[k]) for k in range(2)))
                           for j in range(2)))
        result.update({"response_matrix_transpose": a, "response_matrix_determinant": a[0][0]*a[1][1]-a[0][1]*a[1][0],
                       "g_Ry_per_electron": g, "g_uncertainty": u_g})
        if u_response_path:
            u = json.loads(u_response_path.read_text())
            dqdu = [float(v) for v in u["dq_dU_electron_per_eV"]]
            udqdu = [float(v) for v in u["dq_dU_uncertainty"]]
            budget = sum(abs(g[i]*dqdu[i]) for i in range(2))
            u_budget = sum(abs(dqdu[i])*u_g[i] + abs(g[i])*udqdu[i] for i in range(2))
            result.update({"dq_dU_electron_per_eV": dqdu, "dq_dU_uncertainty": udqdu,
                           "N_ASA_Ry_per_eV": budget, "N_ASA_uncertainty_Ry_per_eV": u_budget,
                           "N_ASA_interval_Ry_per_eV": [max(0.0, budget-u_budget), budget+u_budget],
                           "u_response_input_sha256": sha(u_response_path)})
    except ArithmeticError as exc:
        result["gradient_status"] = str(exc)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True, help="common base input deck")
    parser.add_argument("--seed", type=Path, required=True,
                        help="immutable common RS-origin restart seed; not a KS-accepted endpoint")
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--u-response", type=Path)
    args = parser.parse_args()
    for path in (args.binary, args.input, args.seed, *(args.assets / f for f in ASSETS),
                 *(args.source_root / f for f in SOURCES)):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=False)
    rows = [run_case(spec, args.binary, args.input, args.seed, args.assets, args.output, args.source_root)
            for spec in case_specs()]
    summary = {
        "status": "CAPTURED",
        "binary_sha256": sha(args.binary),
        "common_input_sha256": sha(args.input),
        "common_seed_sha256": sha(args.seed),
        "source_sha256": {name: sha(args.source_root / name) for name in SOURCES},
        "asset_sha256": {name: sha(args.assets / name) for name in ASSETS},
        "cases": rows,
        "analysis": derive(rows, args.u_response),
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": summary["status"], "case_count": len(rows),
                      "summary": str(summary_path)}, sort_keys=True))


if __name__ == "__main__":
    main()
