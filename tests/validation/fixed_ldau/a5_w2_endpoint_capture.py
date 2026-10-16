#!/usr/bin/env python3
"""Converge and capture the two A5.4 W2 endpoint states.

Each endpoint first converges from an independent copy of the common initial
restart with lightweight moments-only output. A second nstep=1 run restarts
from that endpoint and captures the complete selected-iteration radial state
needed by a5_w2_frozen_replay.py. The full trace is therefore written once per
endpoint, not during SCF convergence.
"""
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

from a5_local_w1_controls import ASSETS, edit_group, sha

LAMBDAS = (-0.01, 0.01)


def has_converged_marker(log: str) -> bool:
    """Match the positive terminal marker, excluding `Not converged!` lines."""
    return re.search(r"(?<!Not )Converged!", log, re.IGNORECASE) is not None


def configure(deck: str, lam: float, nstep: int | None = None) -> str:
    deck = edit_group(deck, "self", {"use_kspace": ".true."})
    deck = edit_group(deck, "hamiltonian", {
        "hoh": ".false.", "hubbard_u_general(1,3)": "0.0", "hubbard_j_general(1,3)": "0.0",
        "hubbard_u_potential_form": "'liechtenstein'", "hubbard_radial_moments": "'local'",
        "hubbard_test_operator": "'w2'", "hubbard_test_lambda": f"{lam:.17g}",
        "hubbard_test_spin": "0"})
    if nstep is not None:
        deck = edit_group(deck, "self", {"nstep": str(nstep)})
    return deck


def stage_assets(assets: Path, dest: Path) -> None:
    for name in ASSETS:
        shutil.copyfile(assets / name, dest / name)


def execute(binary: Path, case: Path, trace: Path, lam: float, capture: bool,
            moments_only: bool | None = None) -> tuple[int, str]:
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1", "RSLMTO_FIXED_U_TRACE": str(trace),
                "RSLMTO_FIXED_U_TRACE_ITERATION": "1" if capture else "0",
                "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
                "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1" if (moments_only if moments_only is not None else not capture) else "0",
                "RSLMTO_FIXED_U_LAMBDA": f"{lam:.17g}",
                "RSLMTO_FIXED_U_RUN_ID": f"a5_w2_lambda_{lam:+.3f}_{'capture' if capture else 'converge'}"})
    stdout_path = case / "stdout.log"
    with stdout_path.open("wb") as stdout_file:
        proc = subprocess.run([str(binary)], cwd=case, env=env, stdin=subprocess.DEVNULL,
                              stdout=stdout_file, stderr=subprocess.STDOUT, check=False)
    log = stdout_path.read_text(errors="replace")
    return proc.returncode, log


def require_ks_thermo(log: str, case: Path) -> None:
    if "run_recursion: use_kspace=.true., skipping recursion solver stage" not in log:
        raise RuntimeError(f"reciprocal KS branch marker missing in {case}; do not interpret this run")
    pairs = re.findall(r"entropy=\s*([^\s]+).*?TS=\s*([^\s]+)", log)
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    if not pairs or not all(re.fullmatch(number, value) for value in pairs[-1]):
        raise RuntimeError(f"finite entropy/TS marker missing in {case}; do not interpret this run")


def converge(args) -> None:
    args.output.mkdir(parents=True, exist_ok=False)
    receipts = []
    for lam in LAMBDAS:
        case = args.output / f"w2_lambda_{lam:+.3f}_converge"
        case.mkdir()
        input_path, seed_path = case / "input.nml", case / "Fe.nml"
        input_path.write_text(configure(args.input.read_text(), lam))
        shutil.copyfile(args.seed, seed_path)
        pre_input_sha, pre_seed_sha = sha(input_path), sha(seed_path)
        stage_assets(args.assets, case)
        trace = case / "moments_trace.csv"
        rc, log = execute(args.binary, case, trace, lam, capture=False)
        require_ks_thermo(log, case)
        if rc != 0 or not has_converged_marker(log) or not (case / "Fe_out.nml").is_file():
            raise RuntimeError(f"W2 lambda={lam:+.3f} did not converge under the existing threshold; see {case}")
        receipts.append({"lambda_ry": lam, "status": "CONVERGED", "returncode": rc,
            "initial_seed_label": "shared near-converged initial restart; not an accepted endpoint",
            "input_pre_sha256": pre_input_sha, "seed_pre_sha256": pre_seed_sha,
            "runtime_Fe_nml_sha256": sha(seed_path), "endpoint_Fe_out_sha256": sha(case / "Fe_out.nml"),
            "trace_sha256": sha(trace), "stdout_sha256": sha(case / "stdout.log"),
            "endpoint_dir": str(case)})
    report = {"status": "A5.4_W2_ENDPOINTS_CONVERGED", "binary_sha256": sha(args.binary),
              "base_input_sha256": sha(args.input), "base_initial_seed_sha256": sha(args.seed),
              "runs": receipts}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "summary": str(args.output / "summary.json")}, indent=2))


def capture(args) -> None:
    source = args.endpoint / "Fe_out.nml"
    if not source.is_file():
        raise FileNotFoundError(f"converged W2 endpoint restart missing: {source}")
    args.outdir.mkdir(parents=True, exist_ok=False)
    input_path = args.outdir / "input.nml"
    input_path.write_text(configure(args.input.read_text(), args.lambda_ry, nstep=1))
    restart_path = args.outdir / "Fe_out.nml"
    shutil.copyfile(source, restart_path)
    restart_sha = sha(restart_path)
    stage_assets(args.assets, args.outdir)
    trace = args.outdir / "records.csv"
    moments_only = args.mode == "moments"
    rc, log = execute(args.binary, args.outdir, trace, args.lambda_ry, capture=True,
                      moments_only=moments_only)
    require_ks_thermo(log, args.outdir)
    if rc != 0 or not trace.is_file():
        raise RuntimeError(f"one-step W2 endpoint capture failed; see {args.outdir}")
    with trace.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    iters = {int(r["iteration"]) for r in rows}
    if iters != {1}:
        raise ValueError(f"capture must contain only selected iteration 1, got {iters}")
    kinds = {r["kind"] for r in rows}
    required = ({"test_density_matrix", "test_q0_trace", "test_trace_wp", "radial_moment",
                 "radial_parameter", "ql"} if moments_only else
                {"ldm_m1", "radial_moment", "radial_parameter", "ql", "radial_replay_ql",
                 "radial_replay_pl", "radial_replay_mesh", "radial_replay_v",
                 "radial_replay_rho_in", "radial_replay_ec_seed", "radial_replay_eigen_seed"})
    missing = required - kinds
    if missing:
        raise ValueError(f"{args.mode} one-step W2 capture is missing {sorted(missing)}")
    # Check complete 5x5 complex P/M1 blocks where supplied by this mode.
    for kind, stage, spins in (("test_density_matrix", "P_OUT", (0, 1)),
                               ("ldm_m1", "M1_OUT", (1, 2))):
        if kind not in kinds:
            continue
        for spin in spins:
            block = [r for r in rows if r["kind"] == kind and r["stage"] == stage
                     and int(r["spin"]) == spin and int(r["l"]) == 2]
            if len(block) != 25:
                raise ValueError(f"incomplete {kind}/{stage} spin={spin}: {len(block)} rows")
    receipt = {"status": "A5.4_W2_SELECTED_STATE_CAPTURED", "lambda_ry": args.lambda_ry,
        "binary_sha256": sha(args.binary), "converged_endpoint_restart_sha256": sha(source),
        "capture_restart_pre_sha256": restart_sha, "input_pre_sha256": sha(input_path),
        "capture_runtime_Fe_out_sha256": sha(restart_path),
        "records_sha256": sha(trace), "records_bytes": trace.stat().st_size,
        "selected_iteration": 1, "row_count": len(rows), "trace_kinds": sorted(kinds),
        "all_rows_one_iteration": True, "stdout_sha256": sha(args.outdir / "stdout.log"),
        "full_trace_capture_only": True, "moments_only": moments_only, "capture_mode": args.mode,
        "capture_log_has_converged_marker": has_converged_marker(log),
        "capture_dir": str(args.outdir)}
    (args.outdir / "manifest.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("converge")
    c.add_argument("--binary", type=Path, required=True); c.add_argument("--input", type=Path, required=True)
    c.add_argument("--seed", type=Path, required=True); c.add_argument("--assets", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True); c.set_defaults(func=converge)
    s = sub.add_parser("capture")
    s.add_argument("--binary", type=Path, required=True); s.add_argument("--input", type=Path, required=True)
    s.add_argument("--endpoint", type=Path, required=True); s.add_argument("--assets", type=Path, required=True)
    s.add_argument("--lambda-ry", type=float, choices=(-0.01, 0.01), required=True)
    s.add_argument("--mode", choices=("full", "moments"), default="full")
    s.add_argument("--outdir", type=Path, required=True); s.set_defaults(func=capture)
    args = ap.parse_args()
    for pth in (args.binary, args.input, args.assets, *(args.assets / n for n in ASSETS)):
        if not pth.exists(): raise FileNotFoundError(pth)
    args.func(args)


if __name__ == "__main__":
    main()
