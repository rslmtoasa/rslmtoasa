#!/usr/bin/env python3
"""Canonical SCF-B0C CPU/GPU benchmark harness.

The executable is the existing real-material SCF probe.  This driver stages
fixtures, controls CPU/GPU choices, archives raw output, checks the emitted
exclusive profile, and derives JSON/CSV/Markdown from one row dataset for
both reciprocal and real-space routes.  It never substitutes a synthetic
Hamiltonian or silently changes a production input.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "rslmto.scf-b0c.v2"
BENCHMARK_LEVELS = (
    "eigensolver", "reciprocal_phase", "rs_kernel", "rs_electronic_structure",
    "scf_iteration", "scf_convergence",
)
PROFILE_PHASES = (
    "P_hamiltonian_prepare", "P_hk_assembly", "P_eigensolver", "P_eigenpair_transfer",
    "P_occupations_fermi", "P_density_build", "P_charge_spin_accumulate",
    "P_potential_update", "P_mixing", "P_scf_io", "P_scf_misc",
    "P_rs_hamiltonian_prepare", "P_rs_solver_kernel", "P_rs_green_function",
    "P_rs_spectral_reconstruct", "P_rs_energy_integration", "P_rs_fermi",
    "P_rs_density_build", "P_rs_charge_spin_accumulate",
)
DETAIL_TIMERS = (
    "T_H2D", "T_solver", "T_D2H", "T_total_steady",
    "T_rs_H2D", "T_rs_kernel", "T_rs_D2H", "T_rs_sync", "T_rs_setup",
)
BACKEND_DETAIL_FIELDS = (
    "gpu_device", "T_host_conversion", "T_host_staging", "T_host_widen", "T_sync",
    "h2d_bytes", "d2h_values_bytes", "d2h_vectors_bytes", "workspace_queries",
    "workspace_reuses", "pinned_host_active", "gpu_free_bytes", "gpu_total_bytes",
)
CPU_THREAD_SWEEP = (1, 2, 4, 8)
FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?$")
ROW_KEY_FIELDS = (
    "scf_route", "cross_route_case_id", "material", "fixture_id", "fixture_revision", "supercell", "Natom", "nmat", "nsp", "SOC",
    "basis", "lmax", "strux_backend", "Nk_nominal", "Nk_unique", "k_mesh", "smearing",
    "temperature", "electron_count", "fermi_policy", "mixing_method", "mixing_parameters",
    "convergence_threshold", "starting_state", "potential_identity", "feature_flags",
    "eigenvectors", "hamiltonian_precision", "eigensolver_precision", "eigenvector_precision",
    "density_accumulation_precision", "scf_canonical_precision", "numeric_mode",
    "rs_solver", "rs_backend", "recursion_depth", "block_size", "terminator",
    "chebyshev_order", "chebyshev_kernel", "spectral_bounds_policy",
)
PRECISION_KEY_FIELDS = {
    "hamiltonian_precision", "eigensolver_precision", "eigenvector_precision",
    "density_accumulation_precision", "scf_canonical_precision", "numeric_mode",
}

CSV_FIELDS = (
    "row_id", "benchmark_level", "scf_route", "cross_route_case_id", "backend", "material", "supercell", "Natom", "nmat", "nnz", "nsp", "SOC", "basis",
    "Nk_nominal", "Nk_unique", "starting_state", "numeric_mode", "solver_strategy", "vectors",
    "rs_solver", "rs_backend", "recursion_depth", "block_size", "terminator", "chebyshev_order", "chebyshev_kernel", "spectral_bounds_policy",
    "OMP_threads", "BLAS_threads", "GPU_device", "profile_status", "correctness_status", "rs_detail_timers_status",
    "P_hamiltonian_prepare", "P_hk_assembly", "P_eigensolver", "P_occupations_fermi", "P_density_build", "P_potential_update",
    "P_mixing", "P_scf_io", "P_scf_misc", "P_scf_iteration_total", "steady_iteration_median",
    "P_eigenpair_transfer", "P_charge_spin_accumulate", "T_H2D", "T_solver", "T_D2H", "T_total_steady",
    "P_rs_hamiltonian_prepare", "P_rs_solver_kernel", "P_rs_green_function", "P_rs_spectral_reconstruct",
    "P_rs_energy_integration", "P_rs_fermi", "P_rs_density_build", "P_rs_charge_spin_accumulate",
    "T_rs_H2D", "T_rs_kernel", "T_rs_D2H", "T_rs_sync", "T_rs_setup",
    "T_host_conversion", "T_host_staging", "T_host_widen", "T_sync", "h2d_bytes", "d2h_values_bytes",
    "d2h_vectors_bytes", "workspace_queries", "workspace_reuses", "pinned_host_active", "gpu_free_bytes", "gpu_total_bytes",
    "first_scf_iteration", "cold_process_wall", "n_scf_iterations", "full_scf_wall",
    "S_solver", "S_reciprocal", "S_rs_kernel", "S_rs_phase", "S_iteration", "S_convergence",
    "R_rs_kernel_production", "R_rs_phase_production", "R_iteration_production", "R_convergence_production",
    "headline_speedup_eligible", "best_cpu_row_id", "best_cpu_convergence_row_id",
)


def _number(value: str) -> int | float | str:
    if not FLOAT_RE.match(value):
        return value
    number = float(value.replace("D", "E").replace("d", "e"))
    return int(number) if number.is_integer() else number


def _parse_tokens(text: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    tokens = text.split()
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if "=" not in token:
            index += 1
            continue
        key, value = token.split("=", 1)
        if not value and index + 1 < len(tokens) and "=" not in tokens[index + 1]:
            value = tokens[index + 1]
            index += 1
        result[key] = _number(value)
        index += 1
    return result


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def stable_json(value: Any) -> str:
    return json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def stable_hash(value: Any, length: int = 16) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()[:length]


def numeric_mode(hamiltonian: Any, solver: Any, eigenvector: Any, density: Any, canonical: Any) -> str:
    values = {str(value).lower() for value in (hamiltonian, solver, eigenvector, density, canonical)}
    if values == {"fp64"}:
        return "fp64"
    if values == {"fp32"}:
        return "fp32"
    return "mixed"


def build_comparison_key(row: dict[str, Any]) -> dict[str, Any]:
    """Return the strict physics/start-state/numeric fingerprint for an SCF row."""

    defaults = {
        "fixture_revision": "fixture-v1",
        "Nk_nominal": row.get("k_mesh"),
        "k_mesh": row.get("Nk_nominal"),
        "smearing": {"method": row.get("dos_method", "gaussian"), "sigma": row.get("sigma")},
        "temperature": row.get("temperature"),
        "electron_count": row.get("electron_count"),
        "fermi_policy": "auto_from_eigenvalues",
        "mixing_method": row.get("mixing_method", "namelist"),
        "mixing_parameters": row.get("mixing_parameters", {}),
        "convergence_threshold": row.get("convergence_threshold"),
        "starting_state": row.get("starting_state", "normal_initial"),
        "potential_identity": row.get("potential_identity", "fixture-potential"),
        "feature_flags": row.get("feature_flags", {}),
        "eigenvectors": row.get("eigenvectors", "required"),
    }
    key = {field: row.get(field, defaults.get(field)) for field in ROW_KEY_FIELDS}
    # Electron count is a physical input invariant here, while the emitted
    # charge is a floating-point diagnostic. The production mixed route can
    # differ from the FP64 oracle by a few 1e-8 in large supercells; retain
    # that diagnostic in final-state correctness but do not let it split the
    # same physical pairing.
    if key.get("electron_count") is not None:
        key["electron_count"] = round(float(key["electron_count"]), 4)
    return key


def comparison_key_fingerprint(row: dict[str, Any], *, without_mode: bool = False) -> str:
    key = build_comparison_key(row)
    if without_mode:
        # Pairing a mixed GPU route with its FP64 CPU oracle compares the
        # same physical workload, while numeric precision remains a separate
        # headline eligibility gate.
        for field in PRECISION_KEY_FIELDS:
            key.pop(field, None)
    return stable_hash(key)


def _final_metrics(row: dict[str, Any]) -> dict[str, Any]:
    return row.get("final_state", {}) or {}


def validate_scf_profile(iterations: list[dict[str, Any]], closure_tolerance: float = 0.03) -> dict[str, Any]:
    if not iterations:
        return {"status": "FAIL", "reason": "missing_scf_iteration_profile", "profile_closure_error": None}
    errors: list[float] = []
    misc_fractions: list[float] = []
    for item in iterations:
        total = float(item.get("P_scf_iteration_total", 0.0) or 0.0)
        exclusive = sum(float(item.get(name, 0.0) or 0.0) for name in PROFILE_PHASES)
        error = abs(exclusive - total) / max(abs(total), 1.0e-12)
        errors.append(error)
        misc_fractions.append(float(item.get("P_scf_misc", 0.0) or 0.0) / max(total, 1.0e-12))
    max_error = max(errors)
    steady = misc_fractions[1:] if len(misc_fractions) > 1 else misc_fractions
    max_misc = max(steady)
    first_misc = misc_fractions[0]
    passed = max_error <= closure_tolerance and max_misc <= 0.05
    return {
        "status": "PASS" if passed else "FAIL",
        "reason": None if passed else ("profile_closure_failed" if max_error > closure_tolerance else "scf_misc_exceeds_5_percent"),
        "profile_closure_error": max_error,
        "max_misc_fraction": max_misc,
        "first_iteration_misc_fraction": first_misc,
        "iterations": len(iterations),
        "misc_explanation": "unprofiled production boundary retained in P_scf_misc; investigate before headline use",
    }


def validate_pairing(cpu: dict[str, Any], gpu: dict[str, Any], *, convergence: bool = True) -> dict[str, Any]:
    reasons: list[str] = []
    cpu_key = build_comparison_key(cpu)
    gpu_key = build_comparison_key(gpu)
    cpu_mode = cpu.get("numeric_mode", "mixed")
    gpu_mode = gpu.get("numeric_mode", "mixed")
    for field in PRECISION_KEY_FIELDS:
        cpu_key.pop(field, None)
        gpu_key.pop(field, None)
    if cpu.get("scf_route", "reciprocal") != gpu.get("scf_route", "reciprocal"):
        reasons.append("scf_route_mismatch")
    if cpu.get("scf_route", "reciprocal") == "real_space" and gpu.get("scf_route", "reciprocal") == "real_space":
        for field, reason in (
            ("rs_solver", "rs_solver_mismatch"),
            ("recursion_depth", "recursion_depth_mismatch"),
            ("block_size", "block_size_mismatch"),
            ("terminator", "terminator_mismatch"),
            ("chebyshev_order", "chebyshev_order_mismatch"),
            ("chebyshev_kernel", "chebyshev_kernel_mismatch"),
        ):
            if cpu.get(field) != gpu.get(field):
                reasons.append(reason)
    if cpu_key != gpu_key:
        reasons.append("physics_mismatch")
    if cpu_mode != gpu_mode:
        reasons.append("numeric_mode_mismatch")
    if cpu_mode not in {"fp32", "fp64"} or gpu_mode not in {"fp32", "fp64"}:
        reasons.append("mixed_precision_not_headline")
    if gpu.get("backend") != "cuda":
        reasons.append("gpu_backend_not_cuda")
    if gpu.get("fallback_detected"):
        reasons.append("silent_fallback_detected")
    for row in (cpu, gpu):
        if row.get("profile_status") != "PASS":
            reasons.append("profile_failed")
        if row.get("scf_route", "reciprocal") == "real_space":
            kernel_status = row.get("rs_kernel_correctness_status")
            if kernel_status is None:
                kernel_status = row.get("final_state", {}).get("rs_kernel_correctness_status")
            if kernel_status != "PASS":
                reasons.append("rs_kernel_correctness_failed" if kernel_status == "FAIL" else "rs_kernel_correctness_missing")
        status = (_final_metrics(row).get("status") or row.get("correctness_status") or
                  row.get("correctness", {}).get("status"))
        if status == "FAIL":
            reasons.append("correctness_failed")
        elif status != "PASS":
            reasons.append("correctness_missing")
    if convergence and cpu.get("starting_state") != gpu.get("starting_state"):
        reasons.append("starting_state_mismatch")
    unique = list(dict.fromkeys(reasons))
    return {"eligible": not unique, "reasons": unique, "comparison_key_equal": not ("physics_mismatch" in unique)}


def _median(values: Iterable[Any]) -> float | None:
    numbers = [float(value) for value in values if value is not None]
    return statistics.median(numbers) if numbers else None


def _profile_summary(iterations: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {name: _median(item.get(name) for item in iterations) for name in PROFILE_PHASES}
    summary.update({name: _median(item.get(name) for item in iterations) for name in DETAIL_TIMERS})
    summary["P_scf_iteration_total"] = _median(item.get("P_scf_iteration_total") for item in iterations)
    summary["first_scf_iteration"] = float(iterations[0].get("P_scf_iteration_total", 0.0))
    summary["first_reciprocal_solve"] = float(iterations[0].get("T_total_steady", 0.0))
    summary["first_rs_solver_kernel"] = float(iterations[0].get("P_rs_solver_kernel", 0.0) or 0.0)
    steady = iterations[1:] if len(iterations) > 1 else iterations
    summary["steady_iteration_median"] = _median(item.get("P_scf_iteration_total") for item in steady)
    summary["n_scf_iterations"] = len(iterations)
    return summary


def _solver_time(row: dict[str, Any]) -> float | None:
    """Return the comparable route kernel time for this backend."""

    if row.get("scf_route", "reciprocal") == "real_space":
        value = row.get("T_rs_kernel") or row.get("P_rs_solver_kernel")
        return None if value is None else float(value)

    if row.get("backend") == "cuda" and row.get("T_solver") is not None:
        # CUDA device solve work is nested in T_solver. P_eigensolver is the
        # exclusive outer stage and intentionally excludes that nested work.
        return float(row["T_solver"])
    value = row.get("P_eigensolver")
    return None if value is None else float(value)


def _reciprocal_time(row: dict[str, Any]) -> float | None:
    if row.get("scf_route", "reciprocal") != "reciprocal":
        return None
    assembly = row.get("P_hk_assembly")
    solver = _solver_time(row)
    if assembly is None or solver is None:
        return None
    return float(assembly) + solver


def _rs_phase_time(row: dict[str, Any]) -> float | None:
    if row.get("scf_route", "reciprocal") != "real_space":
        return None
    values = [row.get(name) for name in (
        "P_rs_hamiltonian_prepare", "P_rs_solver_kernel", "P_rs_green_function",
        "P_rs_spectral_reconstruct", "P_rs_energy_integration", "P_rs_fermi",
        "P_rs_density_build", "P_rs_charge_spin_accumulate",
    )]
    if not any(value is not None for value in values):
        return None
    return sum(float(value or 0.0) for value in values)


def _environment(repo: Path, build_dir: Path | None = None) -> dict[str, Any]:
    # Reuse the established null-safe provenance collector.  Import lazily so
    # this module remains runnable as a standalone script from any directory.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from benchmark_harness import capture_environment  # type: ignore
    return capture_environment(repo, build_dir)


FIXTURES: dict[str, dict[str, Any]] = {
    "si": {
        "material": "diamondSi", "source": ROOT / "tests/scf/cases/bulk/diamondSi", "supercell": "1x1x1",
        "fixture_id": "diamondSi_reciprocal_scf", "nsp": 1, "soc": "off", "basis": "sp",
    },
    "fe": {
        "material": "bccFe", "source": ROOT / "tests/scf/cases/bulk/bccFe", "supercell": "1x1x1",
        "fixture_id": "bccFe_reciprocal_scf", "nsp": 2, "soc": "off", "basis": "spd",
    },
    "fe2": {
        "material": "bccFe", "source": ROOT / "example/bulk/supercellFe/2x2x2", "supercell": "2x2x2",
        "fixture_id": "bccFe_reciprocal_scf_supercell", "nsp": 2, "soc": "off", "basis": "spd",
    },
    "fe3": {
        "material": "bccFe", "source": ROOT / "example/bulk/supercellFe/3x3x3", "supercell": "3x3x3",
        "fixture_id": "bccFe_reciprocal_scf_supercell", "nsp": 2, "soc": "off", "basis": "spd",
    },
}


def stage_material_fixture(name: str, target: Path) -> dict[str, Any]:
    if name not in FIXTURES:
        raise ValueError(f"unknown material fixture {name!r}; choose from {sorted(FIXTURES)}")
    spec = FIXTURES[name]
    source = Path(spec["source"])
    if not (source / "input.nml").exists():
        raise FileNotFoundError(f"fixture is incomplete: {source}")
    shutil.copytree(source, target)
    return dict(spec)


def _controlled_environment(omp_threads: int | None, blas_threads: int | None) -> dict[str, str]:
    environment = os.environ.copy()
    if omp_threads is not None:
        environment.update({"OMP_NUM_THREADS": str(omp_threads), "OMP_PROC_BIND": "true", "OMP_PLACES": "cores"})
    if blas_threads is not None:
        value = str(blas_threads)
        environment.update({"BLAS_NUM_THREADS": value, "MKL_NUM_THREADS": value, "OPENBLAS_NUM_THREADS": value})
    environment["RSLMTO_SCF_B0C_PROFILE"] = "1"
    return environment


def parse_probe_output(text: str) -> dict[str, Any]:
    iterations: list[dict[str, Any]] = []
    result: dict[str, Any] = {}
    unsupported: str | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("SCF_B0C_ITER "):
            fields = stripped.split()
            item = {"iteration": int(fields[1])}
            item.update(_parse_tokens(" ".join(fields[2:])))
            iterations.append(item)
        elif stripped.startswith("SCF_B0C_RESULT "):
            result = _parse_tokens(stripped[len("SCF_B0C_RESULT "):])
        elif stripped.startswith("SCF_B0C status=UNSUPPORTED") or stripped.startswith("ACCP1B_SCF status=UNSUPPORTED"):
            prefix = "SCF_B0C " if stripped.startswith("SCF_B0C ") else "ACCP1B_SCF "
            unsupported = str(_parse_tokens(stripped[len(prefix):]).get("reason", "unsupported"))
    profile = validate_scf_profile(iterations)
    if unsupported is not None:
        profile = {"status": "UNSUPPORTED", "reason": unsupported}
    return {"iterations": iterations, "result": result, "profile": profile, "unsupported_reason": unsupported}


def _correctness(reference: dict[str, Any] | None, candidate: dict[str, Any], mode: str, *, require_convergence: bool = True) -> dict[str, Any]:
    if candidate.get("status") == "UNSUPPORTED":
        return {"status": "UNSUPPORTED", "reason": candidate.get("unsupported_reason")}
    if not reference:
        return {"status": "NOT_APPLICABLE", "reason": "reference_missing"}
    if require_convergence and candidate.get("converged") not in (True, "true", "TRUE"):
        return {"status": "FAIL", "reason": "candidate_not_converged"}
    if not require_convergence:
        # Kernel/phase/steady-iteration rows deliberately stop before the
        # production convergence flag.  Their validity is the route-specific
        # profile and invariant contract, not a claim of full SCF convergence.
        return {
            "status": "PASS",
            "reason": "full_convergence_not_required_for_benchmark_level",
            "convergence_required": False,
            "reference_converged": reference.get("converged") in (True, "true", "TRUE"),
            "candidate_converged": candidate.get("converged") in (True, "true", "TRUE"),
        }
    tolerances = {"fp64": 5.0e-6, "mixed": 5.0e-4, "fp32": 5.0e-4}
    tolerance = tolerances.get(mode, 5.0e-4)
    fields = ("final_total_energy", "fermi_energy", "total_charge", "site_charge_transfer", "site_moment", "final_residual")
    differences: dict[str, float | None] = {}
    failures: list[str] = []
    for field in fields:
        left, right = reference.get(field), candidate.get(field)
        if left is None or right is None:
            differences[field] = None
            failures.append(field)
            continue
        difference = abs(float(right) - float(left))
        differences[field] = difference
        scale = max(1.0, abs(float(left)))
        if difference > tolerance * scale:
            failures.append(field)
    return {
        "status": "PASS" if not failures else "FAIL", "tolerance": tolerance,
        "differences": differences, "failed_fields": failures,
        "reference_iteration_count": reference.get("scf_iterations"),
        "candidate_iteration_count": candidate.get("scf_iterations"),
        "iteration_count_difference_allowed": True,
    }


def run_probe(*, binary: Path, fixture: Path, backend: str, strategy: str, scf_route: str,
              rs_solver: str, rs_backend: str, dos_method: str,
              sigma: float, temperature: float, nstep: int, benchmark_level: str, omp_threads: int | None,
              blas_threads: int | None, raw_path: Path) -> dict[str, Any]:
    command = [str(binary), "--input", "input.nml", "--backend", backend,
               "--solver-strategy", strategy, "--dos-method", dos_method,
               "--scf-route", scf_route, "--rs-solver", rs_solver, "--rs-backend", rs_backend,
               "--sigma", str(sigma), "--temperature", str(temperature),
               "--nstep", str(nstep), "--benchmark-level", benchmark_level, "--profile"]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=fixture, env=_controlled_environment(omp_threads, blas_threads),
                               capture_output=True, text=True)
    wall = time.perf_counter() - started
    text = (completed.stdout or "") + (completed.stderr or "")
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text(text, encoding="utf-8")
    parsed = parse_probe_output(text)
    result = parsed["result"]
    unsupported = parsed["unsupported_reason"]
    status = "UNSUPPORTED" if unsupported else ("PASS" if completed.returncode == 0 and result else "FAIL")
    summary = _profile_summary(parsed["iterations"]) if parsed["iterations"] else {}
    row = {
        "backend": backend, "solver_strategy": strategy, "dos_method": dos_method,
        "OMP_threads": omp_threads, "BLAS_threads": blas_threads,
        # The physical probe is a fresh process per row. It does not yet
        # expose separate CUDA context/backend-init counters, so preserve the
        # complete cold process wall and keep it distinct from steady timers.
        "cold_process_wall": wall, "full_scf_wall": wall,
        "returncode": completed.returncode, "profile_status": parsed["profile"].get("status", "FAIL"),
        "profile": parsed["profile"], "status": status, "unsupported_reason": unsupported,
        "iterations": parsed["iterations"], "final_state": result, "final_state_raw": result,
        **summary,
    }
    row["fallback_detected"] = backend == "cuda" and not result and not unsupported
    return row


def _precision_fields(strategy: str, backend: str, scf_route: str, rs_solver: str) -> dict[str, str]:
    if scf_route == "real_space":
        if backend != "cuda":
            return {name: "fp64" for name in ("hamiltonian_precision", "eigensolver_precision", "eigenvector_precision", "density_accumulation_precision", "scf_canonical_precision")}
        # The production RS CUDA plugin has FP32 working paths and FP64
        # storage/reconstruction.  Until a route-specific end-to-end precision
        # probe proves otherwise, retain the honest mixed classification.
        return {
            "hamiltonian_precision": "fp64", "eigensolver_precision": "fp32", "eigenvector_precision": "fp32",
            "density_accumulation_precision": "fp64", "scf_canonical_precision": "fp64",
        }
    if backend != "cuda" or not strategy.startswith("fp32_"):
        return {name: "fp64" for name in ("hamiltonian_precision", "eigensolver_precision", "eigenvector_precision", "density_accumulation_precision", "scf_canonical_precision")}
    return {
        "hamiltonian_precision": "fp64", "eigensolver_precision": "fp32", "eigenvector_precision": "fp32",
        "density_accumulation_precision": "fp64", "scf_canonical_precision": "fp64",
    }


def _decorate_row(raw: dict[str, Any], spec: dict[str, Any], *, name: str, nstep: int,
                  sigma: float, temperature: float, benchmark_level: str, potential_identity: str,
                  scf_route: str, rs_backend: str) -> dict[str, Any]:
    result = raw.get("final_state", {})
    route = str(result.get("scf_route", scf_route))
    rs_solver = result.get("rs_solver") if route == "real_space" else None
    if route == "real_space" and not rs_solver:
        rs_solver = raw.get("rs_solver")
    precision = _precision_fields(str(raw["solver_strategy"]), str(raw["backend"]), route, str(rs_solver or ""))
    def _nullable_number(field: str) -> Any:
        value = result.get(field)
        return None if value in (None, 0, "0") else value
    row = {
        **raw, **precision,
        "benchmark_level": benchmark_level, "scf_route": route,
        "material": spec["material"], "fixture_id": spec["fixture_id"],
        "fixture_revision": "current-head-production-fixture-v1", "supercell": spec["supercell"],
        "Natom": result.get("natom"), "nmat": result.get("nmat"),
        "Nk_unique": result.get("nk_unique") if route == "reciprocal" else None,
        "Nk_nominal": result.get("nk_mesh") if route == "reciprocal" else None,
        "k_mesh": result.get("nk_mesh") if route == "reciprocal" else None, "nsp": spec["nsp"],
        "SOC": spec["soc"], "basis": spec["basis"], "lmax": 2 if spec["basis"] == "spd" else 1,
        "strux_backend": "strux_lib" if name == "si" else "legacy_or_fixture",
        "starting_state": "normal_initial", "potential_identity": potential_identity,
        "eigenvectors": "required", "sigma": sigma, "temperature": temperature,
        "electron_count": result.get("total_charge"), "convergence_threshold": "input.nml",
        "mixing_method": "input.nml", "mixing_parameters": "input.nml",
        "feature_flags": {"SOC": spec["soc"], "Hubbard": False, "GBT": False, "CCOR": False},
        "vectors": True, "n_scf_iterations": raw.get("n_scf_iterations"),
        "steady_iteration_median": raw.get("steady_iteration_median"),
        "cross_route_case_id": f"{spec['material']}:{spec['supercell']}:{spec['basis']}:{spec['soc']}",
        "rs_solver": rs_solver,
        "rs_backend": result.get("rs_backend") if route == "real_space" else None,
        "recursion_depth": _nullable_number("recursion_depth"),
        "block_size": _nullable_number("block_size"),
        "terminator": _nullable_number("terminator"),
        "chebyshev_order": _nullable_number("chebyshev_order"),
        "chebyshev_kernel": result.get("chebyshev_kernel") if route == "real_space" and rs_solver == "chebyshev" else None,
        "spectral_bounds_policy": result.get("spectral_bounds_policy") if route == "real_space" and rs_solver == "chebyshev" else None,
        "rs_kernel_correctness_status": result.get("rs_kernel_correctness_status") if route == "real_space" else None,
        "rs_gpu_used": result.get("rs_gpu_used") if route == "real_space" else None,
        "rs_detail_timers_status": result.get("rs_detail_timers_status") if route == "real_space" else None,
        "P_hk_assembly": raw.get("P_hk_assembly"), "P_eigensolver": raw.get("P_eigensolver"),
        "P_eigenpair_transfer": raw.get("P_eigenpair_transfer"),
        "P_occupations_fermi": raw.get("P_occupations_fermi"), "P_density_build": raw.get("P_density_build"),
        "P_charge_spin_accumulate": raw.get("P_charge_spin_accumulate"),
        "P_potential_update": raw.get("P_potential_update"), "P_mixing": raw.get("P_mixing"),
        "P_scf_io": raw.get("P_scf_io"), "P_scf_misc": raw.get("P_scf_misc"),
        "P_scf_iteration_total": raw.get("P_scf_iteration_total"),
        **{name: raw.get(name) for name in PROFILE_PHASES if name.startswith("P_rs_")},
        "T_H2D": raw.get("T_H2D"), "T_solver": raw.get("T_solver"),
        "T_D2H": raw.get("T_D2H"), "T_total_steady": raw.get("T_total_steady"),
        **{name: raw.get(name) for name in DETAIL_TIMERS if name.startswith("T_rs_")},
        "first_scf_iteration": raw.get("first_scf_iteration"),
        **{field: (result.get(field) if raw["backend"] == "cuda" else None)
           for field in BACKEND_DETAIL_FIELDS},
        "correctness_status": "NOT_APPLICABLE", "correctness": {"status": "NOT_APPLICABLE"},
    }
    row["numeric_mode"] = numeric_mode(*(row[field] for field in (
        "hamiltonian_precision", "eigensolver_precision", "eigenvector_precision",
        "density_accumulation_precision", "scf_canonical_precision")))
    row["row_id"] = stable_hash({"name": name, "backend": raw["backend"], "strategy": raw["solver_strategy"],
                                  "omp": raw["OMP_threads"], "dos": raw["dos_method"], "level": benchmark_level,
                                  "route": route, "rs_solver": rs_solver, "rs_backend": rs_backend})
    row["comparison_key"] = build_comparison_key(row)
    row["comparison_key_fingerprint"] = comparison_key_fingerprint(row)
    row["starting_state_id"] = "normal_initial"
    row["nstep_requested"] = nstep
    return row


def add_pairings(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    pairings: list[dict[str, Any]] = []
    for gpu in (row for row in rows if row.get("backend") == "cuda"):
        if gpu.get("status") == "UNSUPPORTED":
            gpu["headline_speedup_eligible"] = False
            gpu["headline_rejection_reasons"] = ["unsupported_gpu_route"]
            continue
        candidates = [row for row in rows if row.get("backend") == "lapack" and comparison_key_fingerprint(row) == gpu.get("comparison_key_fingerprint")]
        # A mixed-precision GPU row is still physically checked against the
        # FP64 CPU oracle, but cannot become an equal-precision headline.
        if not candidates:
            candidates = [row for row in rows if row.get("backend") == "lapack" and comparison_key_fingerprint(row, without_mode=True) == comparison_key_fingerprint(gpu, without_mode=True)]
        cpu = min(candidates, key=lambda row: float(row.get("steady_iteration_median") or math.inf), default=None)
        if cpu is None:
            gpu["headline_speedup_eligible"] = False
            gpu["headline_rejection_reasons"] = ["best_cpu_missing"]
            continue
        decision = validate_pairing(cpu, gpu, convergence=gpu.get("benchmark_level") == "scf_convergence")
        gpu["best_cpu_row_id"] = cpu["row_id"]
        gpu["headline_speedup_eligible"] = decision["eligible"]
        gpu["headline_rejection_reasons"] = decision["reasons"]
        convergence_candidates = [row for row in candidates
                                  if row.get("profile_status") == "PASS"
                                  and row.get("correctness", {}).get("status") == "PASS"
                                  and _final_metrics(row).get("converged") in (True, "true", "TRUE")]
        best_convergence_cpu = min(convergence_candidates,
                                   key=lambda row: float(row.get("full_scf_wall") or math.inf),
                                   default=None)
        if best_convergence_cpu is not None:
            gpu["best_cpu_convergence_row_id"] = best_convergence_cpu["row_id"]
        cpu_iteration = float(cpu.get("steady_iteration_median") or 0.0)
        gpu_iteration = float(gpu.get("steady_iteration_median") or 0.0)
        if cpu_iteration > 0 and gpu_iteration > 0:
            gpu["S_iteration"] = cpu_iteration / gpu_iteration
        cpu_solver = _solver_time(cpu)
        gpu_solver = _solver_time(gpu)
        if cpu_solver and gpu_solver:
            gpu["S_solver"] = cpu_solver / gpu_solver
        cpu_reciprocal = _reciprocal_time(cpu)
        gpu_reciprocal = _reciprocal_time(gpu)
        if cpu_reciprocal and gpu_reciprocal:
            gpu["S_reciprocal"] = cpu_reciprocal / gpu_reciprocal
        cpu_rs_kernel = _solver_time(cpu) if cpu.get("scf_route", "reciprocal") == "real_space" else None
        gpu_rs_kernel = _solver_time(gpu) if gpu.get("scf_route", "reciprocal") == "real_space" else None
        if cpu_rs_kernel and gpu_rs_kernel:
            gpu["S_rs_kernel"] = cpu_rs_kernel / gpu_rs_kernel
        cpu_rs_phase = _rs_phase_time(cpu)
        gpu_rs_phase = _rs_phase_time(gpu)
        if cpu_rs_phase and gpu_rs_phase:
            gpu["S_rs_phase"] = cpu_rs_phase / gpu_rs_phase
        if (gpu.get("benchmark_level") == "scf_convergence" and best_convergence_cpu is not None
                and best_convergence_cpu.get("full_scf_wall") and gpu.get("full_scf_wall")):
            gpu["S_convergence"] = float(best_convergence_cpu["full_scf_wall"]) / float(gpu["full_scf_wall"])
        # Mixed RS CUDA is a valid production-route comparison when the
        # physical/profile/invariant contract passes, but it is not an
        # equal-precision speedup.  Keep those ratios under explicit R_* names.
        if cpu.get("numeric_mode") != gpu.get("numeric_mode"):
            disallowed = {
                "scf_route_mismatch", "physics_mismatch", "rs_solver_mismatch",
                "recursion_depth_mismatch", "block_size_mismatch", "terminator_mismatch",
                "chebyshev_order_mismatch", "chebyshev_kernel_mismatch", "gpu_backend_not_cuda",
                "silent_fallback_detected", "profile_failed", "correctness_failed",
                "correctness_missing", "starting_state_mismatch", "rs_kernel_correctness_failed",
                "rs_kernel_correctness_missing",
            }
            gpu["production_comparison_eligible"] = not any(reason in disallowed for reason in decision["reasons"])
            if cpu_rs_kernel and gpu_rs_kernel:
                gpu["R_rs_kernel_production"] = cpu_rs_kernel / gpu_rs_kernel
            if cpu_rs_phase and gpu_rs_phase:
                gpu["R_rs_phase_production"] = cpu_rs_phase / gpu_rs_phase
            if cpu_iteration > 0 and gpu_iteration > 0:
                gpu["R_iteration_production"] = cpu_iteration / gpu_iteration
            if (gpu.get("benchmark_level") == "scf_convergence" and best_convergence_cpu is not None
                    and best_convergence_cpu.get("full_scf_wall") and gpu.get("full_scf_wall")):
                gpu["R_convergence_production"] = float(best_convergence_cpu["full_scf_wall"]) / float(gpu["full_scf_wall"])
            for field in ("S_solver", "S_reciprocal", "S_rs_kernel", "S_rs_phase", "S_iteration", "S_convergence"):
                gpu.pop(field, None)
        else:
            gpu["production_comparison_eligible"] = False
        pairings.append({"cpu_row_id": cpu["row_id"], "gpu_row_id": gpu["row_id"], **decision})
    return pairings


def _csv_value(row: dict[str, Any], field: str) -> Any:
    if field == "profile_status":
        return row.get("profile_status")
    if field == "correctness_status":
        return row.get("correctness", {}).get("status", row.get("correctness_status"))
    if field == "GPU_device":
        return row.get("gpu_device") if row.get("backend") == "cuda" else None
    if field == "full_scf_wall":
        return row.get("full_scf_wall")
    return row.get(field)


def write_outputs(campaign: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_jsonable(campaign), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    correctness_root = output.parent / "correctness"
    history_root = output.parent / "iteration_history"
    correctness_root.mkdir(parents=True, exist_ok=True)
    history_root.mkdir(parents=True, exist_ok=True)
    for row in campaign.get("rows", []):
        row_id = str(row.get("row_id"))
        (correctness_root / f"{row_id}.json").write_text(
            json.dumps(_jsonable({"row_id": row_id, "correctness": row.get("correctness"),
                                  "profile": row.get("profile"), "final_state": row.get("final_state")}),
                       indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (history_root / f"{row_id}.json").write_text(
            json.dumps(_jsonable({"row_id": row_id, "iterations": row.get("iterations", [])}),
                       indent=2, sort_keys=True) + "\n", encoding="utf-8")
    csv_path = output.with_suffix(".csv")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in campaign.get("rows", []):
            writer.writerow({field: _csv_value(row, field) for field in CSV_FIELDS})
    rows = campaign.get("rows", [])
    eligible = [row for row in rows if row.get("headline_speedup_eligible")]
    lines = [
        "# SCF-B0C benchmark closure", "", "Harness validation only; these rows are not SCF-B1 application conclusions.", "",
        "## Environment", "", "```json", json.dumps(_jsonable(campaign.get("environment", {})), indent=2, sort_keys=True), "```", "",
        "## Profile closure", "", "| row | route | solver | level | backend | mode | profile | correctness | iteration (s) | convergence wall (s) |", "|---|---|---|---|---|---|---|---|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row.get('row_id')} | {row.get('scf_route', 'reciprocal')} | {row.get('rs_solver') or '-'} | {row.get('benchmark_level')} | {row.get('backend')} | {row.get('numeric_mode')} | {row.get('profile_status')} | {row.get('correctness', {}).get('status')} | {row.get('steady_iteration_median') or '-'} | {row.get('full_scf_wall') or '-'} |")
    lines += ["", "## Headline eligibility", "", f"Eligible GPU rows: {len(eligible)}", "", "| GPU row | route | best CPU | S_solver | S_reciprocal | S_rs_kernel | S_rs_phase | S_iteration | S_convergence | reasons |", "|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for row in (item for item in rows if item.get("backend") == "cuda"):
        lines.append(f"| {row.get('row_id')} | {row.get('scf_route', 'reciprocal')} | {row.get('best_cpu_row_id', '-')} | {row.get('S_solver', '-')} | {row.get('S_reciprocal', '-')} | {row.get('S_rs_kernel', '-')} | {row.get('S_rs_phase', '-')} | {row.get('S_iteration', '-')} | {row.get('S_convergence', '-')} | {','.join(row.get('headline_rejection_reasons', [])) or 'PASS'} |")
    lines += ["", "## Definitions", "", "- `S_solver`, `S_reciprocal`, `S_rs_kernel`, `S_rs_phase`, `S_iteration`, and `S_convergence` are separate ratios.", "- Reciprocal solver ratios use nested `T_solver`; RS kernel ratios use the production `P_rs_solver_kernel` boundary.", "- GPU SCF headline rows require matching route/physics/start state, valid kernel and final correctness, profile closure, and equal numeric mode.", "- Mixed precision rows remain visible and are checked against the FP64 CPU oracle, but are not equal-precision headlines.", "- CUDA unsupported rows are retained with a reason and are never treated as CPU timings."]
    output.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_campaign(args: argparse.Namespace) -> dict[str, Any]:
    args.binary = args.binary.resolve()
    materials = [item.strip() for item in args.materials.split(",") if item.strip()]
    for material in materials:
        if material not in FIXTURES:
            raise ValueError(f"unknown material fixture {material!r}")
    strategies = [item.strip() for item in args.strategies.split(",") if item.strip()]
    threads = [int(item) for item in args.omp_threads.split(",") if item.strip()]
    scf_route = getattr(args, "scf_route", "reciprocal")
    rs_solvers = [item.strip() for item in getattr(args, "rs_solvers", "block").split(",") if item.strip()]
    rs_backend = getattr(args, "rs_backend", "csr")
    if scf_route not in {"reciprocal", "real_space"}:
        raise ValueError(f"unknown SCF route {scf_route!r}")
    if scf_route == "real_space" and any(item not in {"block", "chebyshev", "lanczos"} for item in rs_solvers):
        raise ValueError("RS solvers must be block, chebyshev, or lanczos")
    output = args.output
    raw_root = output.parent / "raw"
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    environment = _environment(ROOT, args.build_dir)
    with tempfile.TemporaryDirectory(prefix="rslmto-scf-b0c-") as temporary:
        scratch = Path(temporary)
        for material in materials:
            spec = FIXTURES[material]
            route_solvers = rs_solvers if scf_route == "real_space" else [""]
            for route_solver in route_solvers:
                solver_label = route_solver or "reciprocal"
                # CPU OMP 1/2/4/8 rows are retained independently.  The same
                # real-material fixture and physical controls are used for
                # every backend/precision route.
                cpu_rows: list[dict[str, Any]] = []
                for omp in threads:
                    fixture = scratch / f"{material}_{solver_label}_cpu_{omp}"
                    stage_material_fixture(material, fixture)
                    raw = run_probe(binary=args.binary, fixture=fixture, backend="lapack",
                                    strategy="fp64_zheevd", scf_route=scf_route,
                                    rs_solver=route_solver, rs_backend=rs_backend,
                                    dos_method=args.dos_method, sigma=args.sigma, temperature=args.temperature,
                                    nstep=args.nstep, benchmark_level=args.benchmark_level, omp_threads=omp,
                                    blas_threads=args.blas_threads,
                                    raw_path=raw_root / f"{material}_{solver_label}_cpu_{omp}.log")
                    row = _decorate_row(raw, spec, name=material, nstep=args.nstep, sigma=args.sigma,
                                        temperature=args.temperature, benchmark_level=args.benchmark_level,
                                        potential_identity=material, scf_route=scf_route, rs_backend=rs_backend)
                    row["environment"] = {**environment, "omp_threads": omp, "blas_threads": args.blas_threads}
                    cpu_rows.append(row)
                    rows.append(row)
                reference = next((row for row in cpu_rows if row.get("status") == "PASS"), None)
                gpu_specs = (() if args.skip_cuda else (strategies if scf_route == "reciprocal" else [route_solver]))
                for strategy in gpu_specs:
                    backend = "cuda"
                    fixture = scratch / f"{material}_{solver_label}_cuda_{strategy}"
                    stage_material_fixture(material, fixture)
                    raw = run_probe(binary=args.binary, fixture=fixture, backend=backend, strategy=strategy,
                                    scf_route=scf_route, rs_solver=route_solver, rs_backend=rs_backend,
                                    dos_method=args.dos_method, sigma=args.sigma, temperature=args.temperature,
                                    nstep=args.nstep, benchmark_level=args.benchmark_level, omp_threads=1,
                                    blas_threads=args.blas_threads,
                                    raw_path=raw_root / f"{material}_{solver_label}_cuda_{strategy}.log")
                    row = _decorate_row(raw, spec, name=material, nstep=args.nstep, sigma=args.sigma,
                                        temperature=args.temperature, benchmark_level=args.benchmark_level,
                                        potential_identity=material, scf_route=scf_route, rs_backend=rs_backend)
                    row["environment"] = {**environment, "omp_threads": 1, "blas_threads": args.blas_threads}
                    if raw.get("status") == "UNSUPPORTED":
                        skipped.append({"row_id": row["row_id"], "status": "UNSUPPORTED", "reason": raw.get("unsupported_reason"), "row": row})
                    else:
                        row["correctness"] = _correctness(
                            _final_metrics(reference) if reference else None,
                            _final_metrics(row),
                            row["numeric_mode"],
                            require_convergence=args.benchmark_level == "scf_convergence",
                        )
                        row["correctness_status"] = row["correctness"].get("status")
                    rows.append(row)
    # CPU self/reference rows are correctness PASS only when the production
    # run converged; their final state remains the oracle for candidate rows.
    for row in rows:
        if row.get("backend") == "lapack":
            final = _final_metrics(row)
            converged = final.get("converged") in (True, "true", "TRUE")
            row["correctness"] = {
                "status": "PASS" if converged or row.get("benchmark_level") != "scf_convergence" else "NOT_CONVERGED",
                "oracle": True,
                "convergence_required": row.get("benchmark_level") == "scf_convergence",
            }
            row["correctness_status"] = row["correctness"]["status"]
    pairings = add_pairings(rows)
    return {
        "schema": SCHEMA, "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "campaign": {"name": "SCF-B0C", "benchmark_level": args.benchmark_level, "scf_route": scf_route,
                      "rs_solvers": rs_solvers, "rs_backend": rs_backend, "dos_method": args.dos_method,
                      "strategies": strategies, "omp_threads": threads, "nstep": args.nstep,
                      "production_executable": str(args.binary), "no_synthetic_fixture": True},
        "environment": environment, "rows": rows, "skipped_rows": skipped, "pairings": pairings,
        "policy": {"cpu_oracle": True, "cuda_fallback_is_invalid": True, "full_fp32_scf_supported": False,
                   "eigenvectors_required": True, "closure_timing_only": True},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="?", default="run")
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--materials", default="si,fe")
    parser.add_argument("--strategies", default="fp64_zheevd")
    parser.add_argument("--scf-route", choices=("reciprocal", "real_space"), default="reciprocal")
    parser.add_argument("--rs-solvers", default="block")
    parser.add_argument("--rs-backend", choices=("csr", "bsr", "fft", "conv"), default="csr")
    parser.add_argument("--omp-threads", default="1,2,4,8")
    parser.add_argument("--blas-threads", type=int, default=1)
    parser.add_argument("--dos-method", choices=("gaussian", "tetrahedron", "blochl"), default="gaussian")
    parser.add_argument("--sigma", type=float, default=0.01)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--nstep", type=int, default=12)
    parser.add_argument("--benchmark-level", choices=BENCHMARK_LEVELS, default="scf_convergence")
    parser.add_argument("--skip-cuda", action="store_true")
    args = parser.parse_args(argv)
    campaign = run_campaign(args)
    write_outputs(campaign, args.output)
    print(f"WROTE {args.output}")
    print(f"WROTE {args.output.with_suffix('.csv')}")
    print(f"WROTE {args.output.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
