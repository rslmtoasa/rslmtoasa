#!/usr/bin/env python3
"""Run and package the frozen SCF-B1 CPU/GPU campaign.

This is a thin campaign layer over :mod:`scf_b0c`.  It does not introduce a
second probe or synthetic workload: every row is emitted by the production
real-material SCF executable and the shared B0C parser/profile contract.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import shutil
import subprocess
import sys
from argparse import Namespace
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from scf_b0c import (  # noqa: E402
    FIXTURES,
    _environment,
    add_pairings,
    run_campaign,
    stable_hash,
)


SCHEMA = "rslmto.scf-b1.v1"
CSV_FIELDS = (
    "row_id", "run_id", "repetition", "benchmark_level", "scf_route", "cross_route_case_id",
    "backend", "material", "supercell", "Natom", "nmat", "Nk_unique", "nsp", "basis",
    "rs_solver", "rs_backend", "recursion_depth", "block_size", "terminator", "chebyshev_order",
    "chebyshev_kernel", "OMP_threads", "BLAS_threads", "numeric_mode", "solver_strategy",
    "profile_status", "correctness_status", "fallback_detected", "rs_gpu_used", "rs_kernel_correctness_status",
    "first_scf_iteration", "steady_iteration_median", "P_scf_iteration_total", "full_scf_wall",
    "P_hamiltonian_prepare", "P_hk_assembly", "P_eigensolver", "P_eigenpair_transfer",
    "P_occupations_fermi", "P_density_build", "P_charge_spin_accumulate", "P_potential_update",
    "P_mixing", "P_scf_io", "P_scf_misc", "P_rs_hamiltonian_prepare", "P_rs_solver_kernel",
    "P_rs_green_function", "P_rs_spectral_reconstruct", "P_rs_fermi", "P_rs_density_build",
    "P_rs_charge_spin_accumulate", "T_H2D", "T_solver", "T_D2H", "T_total_steady", "T_rs_kernel",
    "S_solver", "S_reciprocal", "S_rs_kernel", "S_rs_phase", "S_iteration", "S_convergence",
    "R_rs_kernel_production", "R_rs_phase_production", "R_iteration_production", "R_convergence_production",
)


def campaign_plan() -> list[dict[str, Any]]:
    """The frozen B1 matrix; actual dimensions from the executable remain authoritative."""

    return [
        {
            "run_id": "rs_block_steady",
            "materials": "fe,fe2",
            "scf_route": "real_space", "rs_solvers": "block", "benchmark_level": "scf_iteration",
            "nstep": 5, "strategies": "fp64_zheevd", "repetitions": 1,
        },
        {
            "run_id": "rs_block_large_steady",
            "materials": "fe3", "scf_route": "real_space", "rs_solvers": "block",
            "benchmark_level": "scf_iteration", "nstep": 2, "strategies": "fp64_zheevd",
            "omp_threads": "1", "repetitions": 1,
        },
        {
            "run_id": "rs_block_convergence_fe",
            "materials": "fe", "scf_route": "real_space", "rs_solvers": "block",
            "benchmark_level": "scf_convergence", "nstep": 80, "strategies": "fp64_zheevd",
            "repetitions": 1,
        },
        {
            "run_id": "rs_chebyshev_si",
            "materials": "si", "scf_route": "real_space", "rs_solvers": "chebyshev",
            "benchmark_level": "scf_convergence", "nstep": 12, "strategies": "fp64_zheevd",
            "repetitions": 1,
        },
        {
            "run_id": "rs_lanczos_reference",
            "materials": "fe", "scf_route": "real_space", "rs_solvers": "lanczos",
            "benchmark_level": "scf_iteration", "nstep": 5, "strategies": "fp64_zheevd",
            "repetitions": 1, "skip_cuda": True,
        },
        {
            "run_id": "reciprocal_steady",
            "materials": "si,fe,fe2,fe3", "scf_route": "reciprocal", "rs_solvers": "block",
            "benchmark_level": "scf_iteration", "nstep": 8,
            "strategies": "fp64_zheevd,fp64_zheevj_batched", "repetitions": 1,
        },
        {
            "run_id": "reciprocal_full_small",
            "materials": "si,fe", "scf_route": "reciprocal", "rs_solvers": "block",
            "benchmark_level": "scf_convergence", "nstep": 80,
            "strategies": "fp64_zheevd,fp64_zheevj_batched", "repetitions": 3,
        },
        {
            "run_id": "reciprocal_full_fe2",
            "materials": "fe2", "scf_route": "reciprocal", "rs_solvers": "block",
            "benchmark_level": "scf_convergence", "nstep": 40,
            "strategies": "fp64_zheevd", "repetitions": 1,
        },
        {
            "run_id": "reciprocal_steady_fe3",
            "materials": "fe3", "scf_route": "reciprocal", "rs_solvers": "block",
            "benchmark_level": "scf_iteration", "nstep": 8,
            "strategies": "fp64_zheevd", "repetitions": 1,
        },
    ]


def _args_for_plan(plan: dict[str, Any], output: Path, binary: Path, build_dir: Path | None,
                   blas_threads: int, skip_cuda: bool) -> Namespace:
    return Namespace(
        binary=binary, build_dir=build_dir, output=output,
        materials=plan["materials"], strategies=plan.get("strategies", "fp64_zheevd"),
        scf_route=plan.get("scf_route", "reciprocal"), rs_solvers=plan.get("rs_solvers", "block"),
        rs_backend=plan.get("rs_backend", "csr"), omp_threads=plan.get("omp_threads", "1,2,4,8"),
        blas_threads=blas_threads, dos_method=plan.get("dos_method", "gaussian"),
        sigma=plan.get("sigma", 0.01), temperature=plan.get("temperature", 300.0),
        nstep=plan["nstep"], benchmark_level=plan["benchmark_level"],
        skip_cuda=skip_cuda or bool(plan.get("skip_cuda", False)),
    )


def _load_or_run(plan: dict[str, Any], run_dir: Path, *, binary: Path, build_dir: Path | None,
                 blas_threads: int, skip_cuda: bool, force: bool) -> list[dict[str, Any]]:
    run_dir.mkdir(parents=True, exist_ok=True)
    output = run_dir / "campaign.json"
    if output.exists() and not force:
        return json.loads(output.read_text(encoding="utf-8")).get("rows", [])
    campaign = run_campaign(_args_for_plan(plan, output, binary, build_dir, blas_threads, skip_cuda))
    from scf_b0c import write_outputs
    write_outputs(campaign, output)
    return campaign["rows"]


def _decorate_repetition(rows: list[dict[str, Any]], run_id: str, repetition: int) -> list[dict[str, Any]]:
    decorated = []
    for row in rows:
        item = dict(row)
        item["run_id"] = run_id
        item["repetition"] = repetition
        item["row_id"] = stable_hash({"run_id": run_id, "repetition": repetition, "row_id": row.get("row_id")})
        decorated.append(item)
    return decorated


def _copy_evidence(run_root: Path, package_root: Path) -> None:
    for source_kind in ("raw", "correctness", "iteration_history"):
        destination = package_root / source_kind
        destination.mkdir(parents=True, exist_ok=True)
        for run_dir in sorted(run_root.iterdir()):
            source = run_dir / source_kind
            if not source.is_dir():
                continue
            for file in source.iterdir():
                shutil.copy2(file, destination / f"{run_dir.name}__{file.name}")


def _write_csv(rows: list[dict[str, Any]], output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in CSV_FIELDS})


def _git_state() -> dict[str, Any]:
    def command(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=ROOT, check=False, capture_output=True, text=True).stdout.strip()
    return {
        "frozen_git_commit": command("rev-parse", "HEAD"),
        "git_dirty_tracked": bool(command("status", "--porcelain", "--untracked-files=no")),
        "git_status": command("status", "--short"),
    }


def assemble_package(output: Path, rows: list[dict[str, Any]], run_root: Path, plans: list[dict[str, Any]],
                     environment: dict[str, Any], large_evidence: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    output.parent.mkdir(parents=True, exist_ok=True)
    pairings = add_pairings(rows)
    campaign = {
        "schema": SCHEMA,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "freeze": _git_state(),
        "campaign": {
            "name": "SCF-B1", "timing_protocol": {"warmups": 2, "measured_repetitions": 5,
                                                        "full_convergence_repetitions": 3},
            "convergence_semantics": {
                "criterion": "self%mix%delta < self%conv_thr",
                "default_conv_thr": 5.0e-9,
                "additional_criteria": "none; Hubbard-U continuation is inactive in these fixtures",
                "maximum_iteration_policy": "self%nstep bounds the loop",
                "diagnostic_residual_is_production_criterion": True,
            },
            "plans": plans, "run_root": str(run_root), "no_synthetic_fixture": True,
        },
        "environment": environment,
        "rows": rows,
        "pairings": pairings,
        "large_reciprocal_evidence": large_evidence or [],
        "policy": {
            "cpu_omp_threads": [1, 2, 4, 8], "blas_threads": 1,
            "cuda_fallback_is_invalid": True, "full_fp32_scf_supported": False,
            "mixed_rs_cuda_is_production_ratio_only": True,
            "unexposed_rs_detail_timers": "not_exposed_by_backend",
        },
    }
    output.write_text(json.dumps(campaign, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_csv(rows, output.with_suffix(".csv"))
    _copy_evidence(run_root, output.parent)
    return campaign


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "results/benchmarks/scf_b1/campaign.json")
    parser.add_argument("--only", default=None, help="comma-separated plan IDs")
    parser.add_argument("--force", action="store_true", help="rerun existing subcampaigns")
    parser.add_argument("--skip-cuda", action="store_true")
    parser.add_argument("--blas-threads", type=int, default=1)
    parser.add_argument("--large-evidence", type=Path, action="append", default=[],
                        help="existing ACC-P0/ACC-P4 JSON evidence to embed")
    args = parser.parse_args(argv)
    binary = args.binary.resolve()
    build_dir = args.build_dir.resolve() if args.build_dir else None
    package_root = args.output.resolve().parent
    run_root = package_root / "runs"
    run_root.mkdir(parents=True, exist_ok=True)
    wanted = set(args.only.split(",")) if args.only else None
    plans = [plan for plan in campaign_plan() if wanted is None or plan["run_id"] in wanted]
    if not plans:
        raise SystemExit("no campaign plans selected")
    rows: list[dict[str, Any]] = []
    if args.only and args.output.exists():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
        # Incremental lane execution is useful for long GPU campaigns. Keep
        # completed rows from other lanes while replacing only selected plan
        # IDs (and their repetitions) when --force is requested.
        rows.extend(row for row in previous.get("rows", []) if row.get("run_id") not in wanted)
    for plan in plans:
        for repetition in range(1, int(plan.get("repetitions", 1)) + 1):
            run_id = f"{plan['run_id']}_r{repetition}"
            run_rows = _load_or_run(plan, run_root / run_id, binary=binary, build_dir=build_dir,
                                    blas_threads=args.blas_threads, skip_cuda=args.skip_cuda, force=args.force)
            rows.extend(_decorate_repetition(run_rows, plan["run_id"], repetition))
    large: list[dict[str, Any]] = []
    for path in args.large_evidence:
        document = json.loads(path.read_text(encoding="utf-8"))
        large.append({"source": str(path), "schema": document.get("schema"),
                      "policy": document.get("policy"), "summary": document.get("summary"),
                      "rows": document.get("rows", document.get("before_after", []))})
    package = assemble_package(args.output.resolve(), rows, run_root, campaign_plan(),
                               _environment(ROOT, build_dir), large)
    print(f"WROTE {args.output.resolve()} ({len(package['rows'])} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
