#!/usr/bin/env python3
"""ACC-00 benchmark recorder and comparison tool.

The harness deliberately treats timing as evidence, not as a correctness
assertion.  It runs an existing production command, records repeated wall
times and optional phase records, and writes a small JSON document that can be
compared across runs.

Examples:

  python3 tests/benchmarks/benchmark_harness.py run \
    --name reciprocal_cpu_profile --class component \
    --labels performance microbenchmark reciprocal eigensolver \
    --command build/bin/UnitTddftCpuProfile \
    --output results/benchmarks/reciprocal_cpu_profile.json

  python3 tests/benchmarks/benchmark_harness.py compare \
    results/benchmarks/cpu.json results/benchmarks/gpu.json
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "rslmto.accelerator-benchmark.v1"
REQUIRED_METADATA = (
    "git_commit",
    "compiler",
    "build_type",
    "blas_lapack",
    "omp_threads",
    "mpi_ranks",
    "cuda_toolkit",
    "cusolver",
    "gpu_model",
    "cpu_model",
)

KPM_PROFILE_METADATA = (
    "backend",
    "moment_backend",
    "moment_precision",
    "reconstruction_backend",
    "reconstruction_precision",
    "precision",
    "estimator",
    "N",
    "nnz",
    "M",
    "lld",
    "Ntrace",
    "trace_block_width",
    "trace_batches",
    "OMP_NUM_THREADS",
    "BLAS_NUM_THREADS",
    "omp_threads",
    "blas_threads",
    "random_seed",
    "clock_source",
    "output_mode",
)

NUMERIC_MODES = frozenset(("fp64", "fp32", "mixed"))
KPM_TOLERANCE_SETS = {
    # These envelopes are the established G1/G1.3/VAL-09 closure contract:
    # FP64 comparisons are effectively reference-equal, while the FP32
    # envelope covers the largest validated same-precision Pt spin case.
    # Production files are emitted with ES16.6 formatting.  The FP64
    # comparison envelope is therefore two printed ulps, while the numerical
    # FP64 route remains otherwise reference precision.
    "fp64": {"max_abs": 2.0e-6, "rel_l2": 5.0e-6, "integrated_rel": 5.0e-6},
    "fp32": {"max_abs": 5.0e-4, "rel_l2": 5.0e-4, "integrated_rel": 5.0e-4},
    "mixed": {"max_abs": 5.0e-4, "rel_l2": 5.0e-4, "integrated_rel": 5.0e-4},
}

COMPARISON_KEY_FIELDS = (
    "material", "fixture_id", "fixture_revision", "replication", "N", "nnz",
    "nsp", "soc_state", "cond_type", "current_operator", "cond_calctype",
    "Ntrace", "projector_count", "random_seed", "vector_contract", "M", "lld",
    "NE", "channels", "rc", "energy_min", "energy_max", "fermi", "kernel",
    "chebyshev_scaling", "numeric_mode", "canonical_output_precision",
)

KPM_PHASES = (
    "P_operator_setup",
    "P_trace_setup",
    "P_moments_total",
    "P_gamma_basis_setup",
    "P_gamma_generation",
    "P_reconstruction_total",
    "P_result_unpack",
    "P_energy_integration",
    "P_tensor_postprocess",
    "P_output_prepare",
    "P_output_io",
    "P_stack_setup",
    "P_moment_finalize",
    "P_misc",
)
KPM_PHASES_LEGACY = (
    "P_operator",
    "P_trace_setup",
    "P_moments_total",
    "P_gamma",
    "P_reconstruction_total",
    "P_energy_integration",
    "P_output_io",
    "P_other",
)
KPM_MOMENT_CHILDREN = (
    "D_moment_H2D",
    "D_moment_GPU_kernel",
    "D_moment_D2H",
    "D_conversion",
    "D_gpu_mu_pack",
)
KPM_RECONSTRUCTION_CHILDREN = (
    "D_mu_pack",
    "D_reconstruction_BLAS",
    "D_reconstruction_D2H",
)
KPM_GAMMA_CHILDREN = ("D_gamma_basis", "D_gamma_fill")

PROFILE_DIMENSIONS = re.compile(
    r"^PROFILE_DIMENSIONS\s+(?P<label>\S+)\s+"
    r"sites=(?P<sites>\d+)\s+spinor_basis=(?P<matrix_dimension>\d+)\s+"
    r"nk=(?P<k_points>\d+)\s+mesh=\s*(?P<mesh>.*?)\s+"
    r"nw=(?P<energy_points>\d+)\s*$"
)
PROFILE_RECIPROCAL = re.compile(
    r"^PROFILE_RECIPROCAL\s+(?P<label>\S+)\s+"
    r"fourier_assembly=\s*(?P<fourier_assembly_s>\S+)\s+"
    r"k_eigensolution=\s*(?P<eigensolver_s>\S+)\s+"
    r"arbitrary_kq_assembly_eigensolution=\s*(?P<arbitrary_k_eigensolver_s>\S+)\s+"
    r"pair_operator_construction=\s*(?P<pair_operator_s>\S+)\s*$"
)
PROFILE_MEMORY = re.compile(
    r"^PROFILE_MEMORY_MIB\s+(?P<label>\S+)\s+"
    r"(?P<rest>.*)$"
)
KPM_PROFILE = re.compile(r"^KPM_PROFILE\s+(?P<rest>.*)$")
ACC06_DIMENSIONS = re.compile(
    r"^ACC06_DIMENSIONS\s+fixture=(?P<fixture>\S+)\s+"
    r"backend=(?P<backend>\S+)\s+strategy=(?P<strategy>\S+)\s+"
    r"sites=(?P<sites>\d+)\s+matrix_dimension=(?P<matrix_dimension>\d+)\s+"
    r"nk=(?P<k_points>\d+)\s+tile_size=(?P<tile_size>\d+)\s+"
    r"eigenvectors=(?P<eigenvectors>\d+)\s+lmax=(?P<lmax>\d+)"
)
ACC06_TIMING = re.compile(
    r"^ACC06_TIMING\s+fixture=(?P<fixture>\S+)\s+"
    r"assembly_s=\s*(?P<assembly_s>\S+)\s+solve_s=\s*(?P<solve_s>\S+)\s+"
    r"total_s=\s*(?P<total_s>\S+)"
)
ACCP0_DIMENSIONS = re.compile(r"^ACCP0_DIMENSIONS\s+(?P<rest>.*)$")
ACCP0_TIMING = re.compile(r"^ACCP0_TIMING\s+(?P<rest>.*)$")


def _number(value: str) -> int | float | str:
    try:
        number = float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        return value
    return int(number) if number.is_integer() else number


def _parse_key_values(text: str) -> dict[str, int | float | str]:
    result: dict[str, int | float | str] = {}
    # Fortran list-directed output may render either ``key= value`` or
    # ``key=value`` depending on the format descriptor and compiler.
    tokens = text.split()
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if "=" in token:
            key, value = token.split("=", 1)
            if value:
                result[key] = _number(value)
            elif index + 1 < len(tokens):
                result[key] = _number(tokens[index + 1])
                index += 1
        index += 1
    return result


def numeric_mode(moment_precision: Any, reconstruction_precision: Any) -> str:
    """Classify the end-to-end arithmetic contract, never by backend name."""

    moment = str(moment_precision).lower()
    reconstruction = str(reconstruction_precision).lower()
    if moment == reconstruction == "fp64":
        return "fp64"
    if moment == reconstruction == "fp32":
        return "fp32"
    return "mixed"


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


def stable_hash(value: Any, *, length: int = 16) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()[:length]


def build_comparison_key(row: dict[str, Any]) -> dict[str, Any]:
    """Build the physics/numerics fingerprint shared by CPU/GPU rows.

    Performance-strategy fields (OMP count, BLAS count, backend spelling,
    CUDA tiles and stochastic block width) are intentionally absent.
    """

    mode = row.get("numeric_mode") or numeric_mode(
        row.get("moment_precision"), row.get("reconstruction_precision")
    )
    operator = row.get("current_operator")
    if operator is None:
        operator = {
            "va": row.get("va"),
            "vb": row.get("vb"),
            "spin_orbital_selector": row.get("spin_orbital_selector"),
        }
    key = {
        "material": row.get("material"),
        "fixture_id": row.get("fixture_id"),
        "fixture_revision": row.get("fixture_revision"),
        "replication": row.get("replication"),
        "N": row.get("N"),
        "nnz": row.get("nnz"),
        "nsp": row.get("nsp"),
        "soc_state": row.get("soc_state"),
        "cond_type": row.get("cond_type"),
        "current_operator": operator,
        "cond_calctype": row.get("cond_calctype"),
        "Ntrace": row.get("Ntrace", row.get("random_vec_num")),
        "projector_count": row.get("projector_count", row.get("Ntrace")),
        "random_seed": row.get("random_seed") if row.get("cond_calctype") == "random_vec" else
                       row.get("projector_identifier", "per_type_projectors"),
        "vector_contract": row.get(
            "vector_contract",
            "uniform_unit_phase_normalized_by_sqrt_sites" if row.get("cond_calctype") == "random_vec"
            else "deterministic_type_projectors",
        ),
        "M": row.get("M"),
        "lld": row.get("lld"),
        "NE": row.get("NE"),
        "channels": row.get("channels", row.get("NE")),
        "rc": row.get("rc"),
        "energy_min": row.get("energy_min"),
        "energy_max": row.get("energy_max"),
        "fermi": row.get("fermi"),
        "kernel": row.get("kernel", {"name": "Lorentz", "alpha": 6.0}),
        "chebyshev_scaling": row.get(
            "chebyshev_scaling",
            {"contract": "(E-b)/a", "window": [row.get("energy_min"), row.get("energy_max")]},
        ),
        "numeric_mode": mode,
        "canonical_output_precision": row.get("canonical_output_precision", "fp64"),
    }
    return {field: key[field] for field in COMPARISON_KEY_FIELDS}


def comparison_key_fingerprint(row: dict[str, Any]) -> tuple[dict[str, Any], str]:
    key = build_comparison_key(row)
    return key, stable_hash(key)


def _comparison_key_without_mode(row: dict[str, Any]) -> dict[str, Any]:
    key = build_comparison_key(row)
    key.pop("numeric_mode", None)
    return key


def empty_correctness(reason: str = "correctness_missing") -> dict[str, Any]:
    return {
        "status": "NOT_APPLICABLE",
        "reference_row_id": None,
        "tolerance_set": None,
        "moment": {"available": False, "max_abs": None, "rel_l2": None},
        "conductivity_spectrum": {"available": False, "max_abs": None, "rel_l2": None},
        "integrated_or_tensor": {"available": False, "max_abs": None, "rel": None},
        "validation_evidence_id": None,
        "reason": reason,
    }


def _numeric_file_values(path: Path) -> list[float]:
    values: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        for token in line.split():
            try:
                value = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if math.isfinite(value):
                values.append(value)
    return values


def _file_metrics(reference: list[float], candidate: list[float]) -> dict[str, Any]:
    if len(reference) != len(candidate) or not reference:
        return {"available": False, "max_abs": None, "rel_l2": None, "reason": "shape_mismatch_or_empty"}
    differences = [abs(left - right) for left, right in zip(reference, candidate)]
    scale = math.sqrt(sum(value * value for value in reference))
    diff_norm = math.sqrt(sum(value * value for value in (left - right for left, right in zip(reference, candidate))))
    return {
        "available": True,
        "max_abs": max(differences),
        "rel_l2": diff_norm / max(scale, 1.0e-30),
    }


def compare_production_outputs(reference_dir: Path, candidate_dir: Path, *, mode: str,
                               reference_row_id: str | None = None,
                               evidence_dir: Path | None = None) -> dict[str, Any]:
    """Compare canonical production conductivity files without reimplementing physics."""

    if mode not in KPM_TOLERANCE_SETS:
        result = empty_correctness("mixed_rows_are_not_equal_precision_evidence")
        result["reference_row_id"] = reference_row_id
        return result
    reference_files = {path.name: path for path in reference_dir.glob("*cond*.out")}
    candidate_files = {path.name: path for path in candidate_dir.glob("*cond*.out")}
    common = sorted(set(reference_files) & set(candidate_files))
    if not common:
        result = empty_correctness("production_conductivity_outputs_missing")
        result["reference_row_id"] = reference_row_id
        return result
    file_results = {
        name: _file_metrics(_numeric_file_values(reference_files[name]), _numeric_file_values(candidate_files[name]))
        for name in common
    }
    available = [item for item in file_results.values() if item.get("available")]
    if not available:
        result = empty_correctness("production_conductivity_output_shape_mismatch")
        result["reference_row_id"] = reference_row_id
        result["conductivity_spectrum"]["files"] = file_results
        return result
    max_abs = max(float(item["max_abs"]) for item in available)
    rel_l2 = max(float(item["rel_l2"]) for item in available)
    # The first numeric row nearest the Fermi edge is the established VAL-09
    # tensor/integrated-value anchor; all energy samples remain in the spectrum
    # metric above.
    tensor_files = [name for name in common if "cond_total" in name]
    tensor_values = [file_results[name] for name in tensor_files if file_results[name].get("available")]
    tensor_abs = max((float(item["max_abs"]) for item in tensor_values), default=max_abs)
    tensor_rel = max((float(item["rel_l2"]) for item in tensor_values), default=rel_l2)
    tolerance = KPM_TOLERANCE_SETS[mode]
    passed = max_abs <= tolerance["max_abs"] and rel_l2 <= tolerance["rel_l2"] and tensor_rel <= tolerance["integrated_rel"]
    evidence_id = stable_hash({"reference": str(reference_dir), "candidate": str(candidate_dir),
                               "files": common, "mode": mode})
    result = {
        "status": "PASS" if passed else "FAIL",
        "reference_row_id": reference_row_id,
        "tolerance_set": {"name": mode, **tolerance},
        "moment": {"available": False, "max_abs": None, "rel_l2": None,
                   "evidence": "optimized production path intentionally keeps full moments off host",
                   "evidence_id": evidence_id},
        "conductivity_spectrum": {"available": True, "max_abs": max_abs, "rel_l2": rel_l2,
                                   "files": file_results},
        "integrated_or_tensor": {"available": True, "max_abs": tensor_abs, "rel": tensor_rel},
        "validation_evidence_id": evidence_id,
        "compared_files": common,
    }
    if evidence_dir is not None:
        evidence_dir.mkdir(parents=True, exist_ok=True)
        (evidence_dir / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def validate_pairing(cpu: dict[str, Any], gpu: dict[str, Any]) -> dict[str, Any]:
    """Return an auditable decision; never select a nearby row implicitly."""

    reasons: list[str] = []
    cpu_mode = cpu.get("numeric_mode") or numeric_mode(cpu.get("moment_precision"), cpu.get("reconstruction_precision"))
    gpu_mode = gpu.get("numeric_mode") or numeric_mode(gpu.get("moment_precision"), gpu.get("reconstruction_precision"))
    if _comparison_key_without_mode(cpu) != _comparison_key_without_mode(gpu):
        reasons.append("physics_mismatch")
    if cpu_mode != gpu_mode:
        reasons.append("precision_mismatch")
    if cpu_mode not in {"fp32", "fp64"} or gpu_mode not in {"fp32", "fp64"}:
        reasons.append("mixed_precision_not_headline")
    if cpu.get("cond_calctype") == "random_vec" or gpu.get("cond_calctype") == "random_vec":
        if cpu.get("random_seed") != gpu.get("random_seed"):
            reasons.append("seed_mismatch")
        if cpu.get("Ntrace") != gpu.get("Ntrace"):
            reasons.append("Ntrace_mismatch")
        if cpu.get("vector_contract") != gpu.get("vector_contract"):
            reasons.append("vector_contract_mismatch")
    for row in (cpu, gpu):
        correctness = row.get("correctness") or {}
        if correctness.get("status") == "FAIL":
            reasons.append("correctness_failed")
        elif correctness.get("status") != "PASS":
            reasons.append("correctness_missing")
        if row.get("profile_status") != "PASS":
            reasons.append("profile_failed")
        if row.get("output_mode", "production") == "benchmark_no_write":
            reasons.append("benchmark_no_write")
    unique_reasons = list(dict.fromkeys(reasons))
    return {
        "eligible": not unique_reasons,
        "reasons": unique_reasons,
        "cpu_numeric_mode": cpu_mode,
        "gpu_numeric_mode": gpu_mode,
        "comparison_key_equal": _comparison_key_without_mode(cpu) == _comparison_key_without_mode(gpu),
    }


def validate_kpm_profile(
    profile: dict[str, Any], *, closure_tolerance: float = 0.03,
    child_tolerance: float | None = None,
) -> dict[str, Any]:
    """Validate the fixed KPM phase/detail accounting contract.

    Legacy ``T_*``-only records remain parseable as historical evidence, but
    they are deliberately not valid G1.2 profiles because their scopes cannot
    establish closure.
    """

    metrics = profile.get("metrics", {})
    phase_names = KPM_PHASES if all(name in metrics for name in KPM_PHASES) else KPM_PHASES_LEGACY
    if not all(name in metrics for name in phase_names):
        return {
            "valid": False,
            "status": "LEGACY_UNVALIDATED",
            "reason": "missing exclusive P_* phase namespace",
        }
    total = float(metrics.get("T_transport_total", 0.0))
    exclusive_sum = sum(float(metrics[name]) for name in phase_names)
    calculated_error = abs(exclusive_sum - total) / max(abs(total), 1.0e-12)
    tolerance = child_tolerance if child_tolerance is not None else max(0.05 * max(total, 1.0e-12), 0.01)

    def child_sum(names: tuple[str, ...]) -> float:
        return sum(float(metrics.get(name, 0.0)) for name in names)

    moment_children = child_sum(KPM_MOMENT_CHILDREN)
    reconstruction_children = child_sum(KPM_RECONSTRUCTION_CHILDREN)
    gamma_children = child_sum(KPM_GAMMA_CHILDREN)
    parents = {
        "moments": float(metrics["P_moments_total"]),
        "reconstruction": float(metrics["P_reconstruction_total"]),
        "gamma": (
            float(metrics["P_gamma_basis_setup"]) + float(metrics["P_gamma_generation"])
            if "P_gamma_basis_setup" in metrics
            else float(metrics["P_gamma"])
        ),
    }
    children_ok = (
        moment_children <= parents["moments"] + tolerance
        and reconstruction_children <= parents["reconstruction"] + tolerance
        and gamma_children <= parents["gamma"] + tolerance
    )
    other_name = "P_misc" if "P_misc" in metrics else "P_other"
    other_ok = float(metrics[other_name]) >= -tolerance
    reported_status = str(metrics.get("PROFILE_STATUS", "UNKNOWN"))
    valid = calculated_error <= closure_tolerance and children_ok and other_ok and reported_status != "FAIL"
    return {
        "valid": valid,
        "status": "PASS" if valid else "FAIL",
        "reported_status": reported_status,
        "exclusive_sum_s": exclusive_sum,
        "transport_total_s": total,
        "profile_closure_error": calculated_error,
        "moment_children_s": moment_children,
        "reconstruction_children_s": reconstruction_children,
        "gamma_children_s": gamma_children,
        "child_tolerance_s": tolerance,
        "children_ok": children_ok,
        "other_ok": other_ok,
    }


def parse_profile_output(output: str) -> list[dict[str, Any]]:
    """Extract phase records emitted by the existing CPU profile executable."""

    dimensions: dict[str, dict[str, Any]] = {}
    records: dict[str, dict[str, Any]] = {}
    acc06_label: str | None = None
    accp0_label: str | None = None
    for line in output.splitlines():
        match = KPM_PROFILE.match(line.strip())
        if match:
            values = _parse_key_values(match.group("rest"))
            record = records.setdefault("kpm_transport", {"name": "kpm_transport"})
            record["metadata"] = {
                key: values.pop(key)
                for key in KPM_PROFILE_METADATA
                if key in values
            }
            record["metrics"] = {
                key: value for key, value in values.items()
                if key.startswith(("P_", "D_", "T_", "bytes_", "profile_")) or
                key in {
                    "PROFILE_STATUS", "gpu_energy_block", "resident_device_moment_bytes",
                    "host_full_moment_bytes", "full_moment_d2h_bytes",
                    "reduced_result_d2h_bytes",
                }
            }
            record["validation"] = validate_kpm_profile(record)
            record["class"] = "component"
            record["labels"] = ["performance", "component", "rs", "kpm", "transport"]
            continue
        match = ACCP0_DIMENSIONS.match(line.strip())
        if match:
            values = _parse_key_values(match.group("rest"))
            fixture = str(values.get("fixture", "unknown"))
            backend = str(values.get("backend", "unknown"))
            workload = str(values.get("workload", "crossover"))
            accp0_label = "accp0_" + "_".join(
                str(values.get(key, "unknown"))
                for key in ("fixture", "backend", "L", "nominal_mesh", "tile", "eigenvectors")
            )
            dimensions[accp0_label] = {
                "fixture": fixture,
                "backend": backend,
                "workload": workload,
                **values,
            }
            continue
        match = ACCP0_TIMING.match(line.strip())
        if match:
            values = _parse_key_values(match.group("rest"))
            fixture = str(values.get("fixture", "unknown"))
            label = accp0_label
            if label is None or dimensions.get(label, {}).get("fixture") != fixture:
                label = next(
                    (key for key, item in dimensions.items() if item.get("fixture") == fixture),
                    f"accp0_{fixture}",
                )
            record = records.setdefault(label, {"name": label})
            record["metrics"] = {
                key: value for key, value in values.items() if key not in {"fixture", "backend"}
            }
            continue
        match = ACC06_DIMENSIONS.match(line.strip())
        if match:
            values = match.groupdict()
            fixture = values.pop("fixture")
            backend = values.pop("backend")
            strategy = values.pop("strategy")
            acc06_label = f"acc06_{fixture}_{backend}_{strategy}"
            dimensions[acc06_label] = {
                "fixture": fixture,
                "backend": backend,
                "strategy": strategy,
                **{key: _number(value) for key, value in values.items()},
            }
            continue
        match = ACC06_TIMING.match(line.strip())
        if match:
            values = match.groupdict()
            fixture = values.pop("fixture")
            label = acc06_label
            if label is None or dimensions.get(label, {}).get("fixture") != fixture:
                label = next((key for key, item in dimensions.items() if item.get("fixture") == fixture),
                             f"acc06_{fixture}")
            record = records.setdefault(label, {"name": label})
            record["metrics"] = {key: float(value) for key, value in values.items()}
            continue
        match = PROFILE_DIMENSIONS.match(line.strip())
        if match:
            values = match.groupdict()
            label = values.pop("label")
            dimensions[label] = {
                key: _number(value) for key, value in values.items()
            }
            continue
        match = PROFILE_RECIPROCAL.match(line.strip())
        if match:
            values = match.groupdict()
            label = values.pop("label")
            record = records.setdefault(label, {"name": label})
            record["metrics"] = {
                key: float(value) for key, value in values.items()
            }
            continue
        match = PROFILE_MEMORY.match(line.strip())
        if match:
            label = match.group("label")
            record = records.setdefault(label, {"name": label})
            record.setdefault("metrics", {}).update(
                {f"{key}_mib": value for key, value in _parse_key_values(match.group("rest")).items()}
            )

    for label, record in records.items():
        if label in dimensions:
            record["metadata"] = dimensions[label]
        else:
            record.setdefault("metadata", {})
        record.setdefault("class", "component")
        record.setdefault("labels", ["performance", "reciprocal"])
        if record.get("name") == "kpm_transport":
            record["validation"] = validate_kpm_profile(record)
    return list(records.values())


def _read_cache(build_dir: Path | None) -> dict[str, str]:
    if build_dir is None:
        return {}
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return {}
    values: dict[str, str] = {}
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.split(":", 1)[0]
        values[key] = value
    return values


def _command_version(command: str) -> str | None:
    try:
        completed = subprocess.run(
            [command, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    line = (completed.stdout or completed.stderr).splitlines()
    return line[0].strip() if line else None


def _git_commit(repo: Path) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() or None


def _git_dirty(repo: Path) -> bool | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
            check=False, capture_output=True, text=True, timeout=3,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode:
        return None
    if completed.stdout.strip():
        return True
    try:
        untracked = subprocess.run(
            ["git", "-C", str(repo), "ls-files", "--others", "--exclude-standard", "--directory", "-z"],
            check=False, capture_output=True, timeout=5,
        )
    except subprocess.TimeoutExpired:
        # A repository with a very large build tree is conservatively dirty;
        # provenance must never claim a clean tree after an incomplete scan.
        return True
    except (OSError, subprocess.SubprocessError):
        return None
    if untracked.returncode:
        return None
    return bool(untracked.stdout)


def _command_text(command: list[str], timeout: float = 2.0) -> str | None:
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode:
        return None
    return completed.stdout.strip() or None


def _selected_gpu_metadata() -> dict[str, Any]:
    selected = os.environ.get("RSLMTO_GPU_DEVICE") or os.environ.get("CUDA_DEVICE")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    selected_index = selected if selected is not None else (visible.split(",")[0].strip() if visible else "0")
    query = _command_text([
        "nvidia-smi", f"--id={selected_index}",
        "--query-gpu=index,name,memory.total,compute_cap,driver_version",
        "--format=csv,noheader,nounits",
    ])
    if not query:
        return {
            "selected_gpu_index": selected_index if visible or selected else None,
            "gpu_model": None, "gpu_vram_mib": None, "gpu_compute_capability": None,
            "cuda_driver": None,
        }
    fields = [part.strip() for part in query.splitlines()[0].split(",")]
    return {
        "selected_gpu_index": fields[0] if len(fields) > 0 else selected_index,
        "gpu_model": fields[1] if len(fields) > 1 else None,
        "gpu_vram_mib": int(float(fields[2])) if len(fields) > 2 and fields[2].replace(".", "", 1).isdigit() else None,
        "gpu_compute_capability": fields[3] if len(fields) > 3 else None,
        "cuda_driver": fields[4] if len(fields) > 4 else None,
    }


def _cpu_model() -> str | None:
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name") and ":" in line:
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine() or None


def _memory_mib() -> int | None:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) // 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def _cpu_counts() -> tuple[int | None, int | None]:
    logical = os.cpu_count()
    physical: int | None = None
    cores: str | None = None
    sockets: str | None = None
    text = _command_text(["lscpu"])
    if text:
        for line in text.splitlines():
            if line.startswith("Core(s) per socket:"):
                cores = line.split(":", 1)[1].strip()
            elif line.startswith("Socket(s):"):
                sockets = line.split(":", 1)[1].strip()
            else:
                continue
            if cores is None or sockets is None:
                continue
    if cores is not None and sockets is not None:
        try:
            physical = int(cores) * int(sockets)
        except ValueError:
            physical = None
    return physical, logical


def capture_environment(
    repo: Path,
    build_dir: Path | None = None,
    *,
    omp_threads: int | None = None,
    blas_threads: int | None = None,
    mpi_ranks: int | None = None,
    selected_gpu_index: str | None = None,
) -> dict[str, Any]:
    """Capture stable environment fields, retaining nulls for unavailable GPU data."""

    cache = _read_cache(build_dir)
    compiler = cache.get("CMAKE_Fortran_COMPILER") or cache.get("CMAKE_Fortran_COMPILER_ID")
    compiler_version = _command_version(compiler) if compiler else None
    if compiler_version:
        compiler = f"{compiler} ({compiler_version})"
    library_paths = [
        value for key, value in cache.items()
        if (key.startswith("BLAS_") or key.startswith("LAPACK_"))
        and "LIBRARY" in key and value and "NOTFOUND" not in value
    ]
    library_text = ";".join(sorted(set(library_paths)))
    if "mkl" in library_text.lower():
        blas = "Intel oneMKL"
    elif "openblas" in library_text.lower():
        blas = "OpenBLAS"
    else:
        blas = cache.get("BLA_VENDOR") or library_text or None
    omp = omp_threads if omp_threads is not None else os.environ.get("OMP_NUM_THREADS")
    if blas_threads is not None:
        blas_thread_value = str(blas_threads)
    else:
        blas_thread_value = (
            os.environ.get("BLAS_NUM_THREADS")
            or os.environ.get("MKL_NUM_THREADS")
            or os.environ.get("OPENBLAS_NUM_THREADS")
        )
    ranks = mpi_ranks if mpi_ranks is not None else os.environ.get("RSLMTO_MPI_RANKS", "1")
    gpu = _selected_gpu_metadata()
    if selected_gpu_index is not None:
        gpu["selected_gpu_index"] = selected_gpu_index
    cuda_version = cache.get("CUDAToolkit_VERSION") or cache.get("CMAKE_CUDA_COMPILER_VERSION")
    if cuda_version is None:
        nvcc_version = _command_text(["nvcc", "--version"])
        if nvcc_version:
            match = re.search(r"release\s+([0-9.]+)", nvcc_version)
            cuda_version = match.group(1) if match else nvcc_version.splitlines()[-1]
    flags = {
        key: value for key, value in cache.items()
        if key in {"CMAKE_Fortran_FLAGS", "CMAKE_Fortran_FLAGS_RELEASE", "CMAKE_C_FLAGS", "CMAKE_CUDA_FLAGS"}
    }
    physical_cpu_count, logical_cpu_count = _cpu_counts()
    unavailable: dict[str, str] = {}
    for key in ("cpu_governor", "numa_topology"):
        unavailable[key] = "not collected by default" if not os.environ.get("RSLMTO_CAPTURE_OPTIONAL_HW") else "unavailable"
    return {
        "git_commit": _git_commit(repo),
        "git_dirty": _git_dirty(repo),
        "compiler": compiler,
        "compiler_version": compiler_version,
        "build_type": cache.get("CMAKE_BUILD_TYPE") or "unspecified",
        "optimization_flags": flags or None,
        "blas_lapack": blas,
        "omp_threads": int(omp) if str(omp).isdigit() else omp,
        "blas_threads": int(blas_thread_value) if str(blas_thread_value).isdigit() else blas_thread_value,
        "omp_proc_bind": os.environ.get("OMP_PROC_BIND"),
        "omp_places": os.environ.get("OMP_PLACES"),
        "mpi_ranks": int(ranks) if str(ranks).isdigit() else ranks,
        "cuda_toolkit": cuda_version,
        "cusolver": cache.get("CUDAToolkit_cusolver_VERSION"),
        "cuda_driver": gpu["cuda_driver"],
        "gpu_model": gpu["gpu_model"],
        "gpu_vram_mib": gpu["gpu_vram_mib"],
        "gpu_compute_capability": gpu["gpu_compute_capability"],
        "selected_gpu_index": gpu["selected_gpu_index"],
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "cpu_model": _cpu_model(),
        "physical_cpu_count": physical_cpu_count,
        "logical_cpu_count": logical_cpu_count,
        "ram_mib": _memory_mib(),
        "os": platform.platform(),
        "kernel": platform.release(),
        "optional_metadata_unavailable": unavailable,
    }


def _normalise_command(command: Iterable[str]) -> list[str]:
    values = list(command)
    if not values:
        raise ValueError("benchmark command must not be empty")
    return values


def run_command(
    command: Iterable[str],
    *,
    warmups: int,
    repetitions: int,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    persistent: bool = False,
) -> tuple[list[float], list[dict[str, Any]], str]:
    command = _normalise_command(command)
    if warmups < 0 or repetitions < 1:
        raise ValueError("warmups must be non-negative and repetitions must be positive")
    last_output = ""
    if persistent:
        start = time.perf_counter()
        completed = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        last_output = (completed.stdout or "") + (completed.stderr or "")
        if completed.returncode:
            raise RuntimeError(f"persistent benchmark command failed ({completed.returncode}):\n{last_output}")
        return [elapsed], [{item["name"]: item for item in parse_profile_output(last_output)}], last_output
    for _ in range(warmups):
        completed = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
        if completed.returncode:
            raise RuntimeError(f"warm-up command failed ({completed.returncode}):\n{completed.stdout}\n{completed.stderr}")

    wall_times: list[float] = []
    profile_samples: list[dict[str, Any]] = []
    for _ in range(repetitions):
        start = time.perf_counter()
        completed = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        last_output = (completed.stdout or "") + (completed.stderr or "")
        if completed.returncode:
            raise RuntimeError(f"benchmark command failed ({completed.returncode}):\n{last_output}")
        wall_times.append(elapsed)
        profile_samples.append({item["name"]: item for item in parse_profile_output(last_output)})
    return wall_times, profile_samples, last_output


def benchmark_environment(omp_threads: int | None, blas_threads: int | None) -> dict[str, str] | None:
    """Return a controlled child environment when thread counts were given."""

    if omp_threads is None and blas_threads is None:
        return None
    environment = os.environ.copy()
    if omp_threads is not None:
        environment["OMP_NUM_THREADS"] = str(omp_threads)
    if blas_threads is not None:
        value = str(blas_threads)
        environment["BLAS_NUM_THREADS"] = value
        environment["MKL_NUM_THREADS"] = value
        environment["OPENBLAS_NUM_THREADS"] = value
    return environment


def _metadata_with_overrides(metadata: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    overrides = {
        "matrix_dimension": args.matrix_dimension,
        "k_points": args.k_points,
        "tile_size": args.tile_size,
        "eigenvectors": args.eigenvectors,
        "problem_type": args.problem_type,
        "transfer_policy": args.transfer_policy,
    }
    return {key: value for key, value in overrides.items() if value is not None} | metadata


def make_document(
    *,
    name: str,
    benchmark_class: str,
    labels: list[str],
    command: list[str],
    metadata: dict[str, Any],
    wall_times: list[float],
    profile_samples: list[dict[str, Any]],
    last_output: str,
    warmups: int,
    persistent: bool = False,
) -> dict[str, Any]:
    workload_metadata = {
        "matrix_dimension": None,
        "k_points": None,
        "tile_size": None,
        "eigenvectors": None,
        "problem_type": None,
        "transfer_policy": None,
    }
    workload_metadata.update(metadata)
    samples = [{"wall_time_s": value} for value in wall_times]
    records: list[dict[str, Any]] = []
    for sample_index, profiles in enumerate(profile_samples):
        for profile_name, profile in profiles.items():
            validation = profile.get("validation")
            if validation and validation.get("status") == "FAIL":
                raise ValueError(
                    f"KPM profile failed stage accounting: {validation.get('reason', validation)}"
                )
            if sample_index == 0:
                records.append(
                    {
                        "name": f"{name}.{profile_name}",
                        "class": profile.get("class", "component"),
                        "labels": profile.get("labels", labels),
                        "metadata": workload_metadata | profile.get("metadata", {}),
                        "validation": profile.get("validation"),
                        "samples": [],
                    }
                )
            record = records[-1] if records and records[-1]["name"] == f"{name}.{profile_name}" else next(
                item for item in records if item["name"] == f"{name}.{profile_name}"
            )
            record["samples"].append(profile.get("metrics", {}))

    document = {
        "schema": SCHEMA,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "benchmark": {
            "name": name,
            "class": benchmark_class,
            "labels": labels,
            "command": command,
            "metadata": workload_metadata,
            "policy": {
                "warmups": warmups,
                "repetitions": len(wall_times),
                "persistent_process": persistent,
            },
            "samples": samples,
        },
        "profile_records": records,
        "last_output": last_output[-4000:],
    }
    return document


def _load_document(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != SCHEMA:
        raise ValueError(f"{path}: unsupported benchmark schema")
    return document


def _sample_values(document: dict[str, Any]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    benchmark = document.get("benchmark", {})
    for sample in benchmark.get("samples", []):
        for key, value in sample.items():
            if isinstance(value, (int, float)):
                result.setdefault(key, []).append(float(value))
    for record in document.get("profile_records", []):
        for sample in record.get("samples", []):
            for key, value in sample.items():
                if isinstance(value, (int, float)):
                    result.setdefault(f"{record['name']}:{key}", []).append(float(value))
    return result


def _summary(values: list[float]) -> dict[str, float]:
    median = statistics.median(values)
    deviations = [abs(value - median) for value in values]
    sorted_values = sorted(values)
    return {
        "median": median,
        "minimum": min(values),
        "maximum": max(values),
        "spread": max(values) - min(values),
        "mad": statistics.median(deviations),
        "iqr": statistics.quantiles(sorted_values, n=4, method="inclusive")[2] -
               statistics.quantiles(sorted_values, n=4, method="inclusive")[0]
               if len(sorted_values) >= 2 else 0.0,
    }


def compare_documents(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    base_metadata = baseline.get("benchmark", {}).get("metadata", {})
    candidate_metadata = candidate.get("benchmark", {}).get("metadata", {})
    warnings = [
        f"{key}: baseline={base_metadata.get(key)!r}, candidate={candidate_metadata.get(key)!r}"
        for key in REQUIRED_METADATA
        if base_metadata.get(key) != candidate_metadata.get(key)
    ]
    base_values = _sample_values(baseline)
    candidate_values = _sample_values(candidate)
    rows: list[dict[str, Any]] = []
    for name in sorted(set(base_values) | set(candidate_values)):
        row: dict[str, Any] = {"metric": name}
        if name in base_values:
            row["baseline"] = _summary(base_values[name])
        if name in candidate_values:
            row["candidate"] = _summary(candidate_values[name])
        if name in base_values and name in candidate_values and candidate_values[name]:
            row["speedup"] = statistics.median(base_values[name]) / statistics.median(candidate_values[name])
        rows.append(row)
    return {
        "schema": SCHEMA,
        "baseline": baseline.get("benchmark", {}).get("name"),
        "candidate": candidate.get("benchmark", {}).get("name"),
        "environment_mismatch_warnings": warnings,
        "metrics": rows,
    }


def _add_common_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--name", required=True)
    parser.add_argument("--class", dest="benchmark_class", default="component", choices=["microbenchmark", "component", "end_to_end"])
    parser.add_argument("--labels", nargs="+", default=["performance"])
    parser.add_argument("--command", nargs=argparse.REMAINDER, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--cwd", type=Path)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--omp-threads", type=int)
    parser.add_argument("--blas-threads", type=int)
    parser.add_argument("--mpi-ranks", type=int)
    parser.add_argument("--matrix-dimension", type=int)
    parser.add_argument("--k-points", type=int)
    parser.add_argument("--tile-size", type=int)
    parser.add_argument("--eigenvectors", action="store_true")
    parser.add_argument("--problem-type")
    parser.add_argument("--transfer-policy")
    parser.add_argument(
        "--persistent",
        action="store_true",
        help="run the command once; the command owns warm-ups and repetitions inside one process",
    )


def _expand_manifest_command(command: list[str], values: dict[str, str]) -> list[str]:
    expanded: list[str] = []
    for token in command:
        try:
            value = token.format(**values)
        except KeyError as error:
            raise ValueError(f"manifest command uses unknown placeholder {{{error.args[0]}}}") from error
        if value:
            expanded.append(value)
    return expanded


def run_manifest(args: argparse.Namespace) -> int:
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("schema") != "rslmto.accelerator-benchmark-manifest.v1":
        raise ValueError(f"{args.manifest}: unsupported benchmark manifest schema")
    selected = set(args.name or [])
    repo = args.repo.resolve()
    values = {
        "repo": str(repo),
        "binary": str(args.binary.resolve()),
        "profile_binary": str(args.profile_binary.resolve()) if args.profile_binary else "",
        "gpu_flag": "--gpu-plugin" if args.gpu_plugin else "",
        "python": sys.executable,
        "scratch_root": str(args.scratch_root.resolve()),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for entry in manifest["benchmarks"]:
        if selected and entry["name"] not in selected:
            continue
        needs_profile_binary = any(
            "{profile_binary}" in token for token in entry.get("command", [])
        )
        if needs_profile_binary and not args.profile_binary:
            message = (
                f"benchmark {entry['name']!r} requires --profile-binary; "
                "skipping optional profile entry"
            )
            if entry["name"] in selected:
                raise ValueError(message)
            print(f"SKIP {message}")
            continue
        command = _expand_manifest_command(entry["command"], values)
        metadata = capture_environment(
            repo,
            args.build_dir.resolve() if args.build_dir else None,
            omp_threads=args.omp_threads,
            blas_threads=args.blas_threads,
            mpi_ranks=args.mpi_ranks,
        )
        metadata.update(entry.get("metadata", {}))
        metadata["gpu_plugin"] = bool(args.gpu_plugin and "{gpu_flag}" in " ".join(entry.get("command", [])))
        wall_times, profile_samples, last_output = run_command(
            command,
            warmups=args.warmups,
            repetitions=args.repetitions,
            cwd=repo,
            env=benchmark_environment(args.omp_threads, args.blas_threads),
            persistent=False,
        )
        document = make_document(
            name=entry["name"],
            benchmark_class=entry["class"],
            labels=entry["labels"],
            command=command,
            metadata=metadata,
            wall_times=wall_times,
            profile_samples=profile_samples,
            last_output=last_output,
            warmups=args.warmups,
            persistent=False,
        )
        output = args.output_dir / f"{entry['name']}.json"
        output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"WROTE {output}: median={statistics.median(wall_times):.6f}s minimum={min(wall_times):.6f}s")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    run_parser = subparsers.add_parser("run", help="run one benchmark command and write JSON")
    _add_common_run_arguments(run_parser)
    compare_parser = subparsers.add_parser("compare", help="compare two JSON benchmark records")
    compare_parser.add_argument("baseline", type=Path)
    compare_parser.add_argument("candidate", type=Path)
    compare_parser.add_argument("--output", type=Path)
    manifest_parser = subparsers.add_parser("run-manifest", help="run the production benchmark inventory")
    manifest_parser.add_argument("--manifest", type=Path, required=True)
    manifest_parser.add_argument("--binary", type=Path, required=True)
    manifest_parser.add_argument("--profile-binary", type=Path)
    manifest_parser.add_argument("--gpu-plugin", action="store_true",
                                 help="Pass --gpu-plugin to manifest commands that opt into CUDA recursion")
    manifest_parser.add_argument("--output-dir", type=Path, required=True)
    manifest_parser.add_argument("--scratch-root", type=Path, default=Path("/tmp/rslmto-acc00"))
    manifest_parser.add_argument("--name", action="append")
    manifest_parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    manifest_parser.add_argument("--build-dir", type=Path)
    manifest_parser.add_argument("--warmups", type=int, default=1)
    manifest_parser.add_argument("--repetitions", type=int, default=3)
    manifest_parser.add_argument("--omp-threads", type=int)
    manifest_parser.add_argument("--blas-threads", type=int)
    manifest_parser.add_argument("--mpi-ranks", type=int)
    args = parser.parse_args(argv)

    if args.action == "compare":
        report = compare_documents(_load_document(args.baseline), _load_document(args.candidate))
        text = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(text, encoding="utf-8")
        else:
            print(text, end="")
        for warning in report["environment_mismatch_warnings"]:
            print(f"WARNING: environment mismatch: {warning}", file=sys.stderr)
        return 0

    if args.action == "run-manifest":
        return run_manifest(args)

    command = _normalise_command(args.command)
    metadata = capture_environment(
        args.repo.resolve(),
        args.build_dir.resolve() if args.build_dir else None,
        omp_threads=args.omp_threads,
        blas_threads=args.blas_threads,
        mpi_ranks=args.mpi_ranks,
    )
    metadata = _metadata_with_overrides(metadata, args)
    wall_times, profile_samples, last_output = run_command(
        command,
        warmups=args.warmups,
        repetitions=args.repetitions,
        cwd=args.cwd.resolve() if args.cwd else None,
        env=benchmark_environment(args.omp_threads, args.blas_threads),
        persistent=args.persistent,
    )
    document = make_document(
        name=args.name,
        benchmark_class=args.benchmark_class,
        labels=args.labels,
        command=command,
        metadata=metadata,
        wall_times=wall_times,
        profile_samples=profile_samples,
        last_output=last_output,
        warmups=args.warmups,
        persistent=args.persistent,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE {args.output}: median={statistics.median(wall_times):.6f}s minimum={min(wall_times):.6f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
