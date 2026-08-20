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
import datetime as dt
import json
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
    "OMP_NUM_THREADS",
    "BLAS_NUM_THREADS",
    "omp_threads",
    "blas_threads",
    "clock_source",
)

KPM_PHASES = (
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


def validate_kpm_profile(
    profile: dict[str, Any], *, closure_tolerance: float = 0.05,
    child_tolerance: float | None = None,
) -> dict[str, Any]:
    """Validate the fixed KPM phase/detail accounting contract.

    Legacy ``T_*``-only records remain parseable as historical evidence, but
    they are deliberately not valid G1.2 profiles because their scopes cannot
    establish closure.
    """

    metrics = profile.get("metrics", {})
    if not all(name in metrics for name in KPM_PHASES):
        return {
            "valid": False,
            "status": "LEGACY_UNVALIDATED",
            "reason": "missing exclusive P_* phase namespace",
        }
    total = float(metrics.get("T_transport_total", 0.0))
    exclusive_sum = sum(float(metrics[name]) for name in KPM_PHASES)
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
        "gamma": float(metrics["P_gamma"]),
    }
    children_ok = (
        moment_children <= parents["moments"] + tolerance
        and reconstruction_children <= parents["reconstruction"] + tolerance
        and gamma_children <= parents["gamma"] + tolerance
    )
    other_ok = float(metrics["P_other"]) >= -tolerance
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
                key in {"PROFILE_STATUS", "gpu_energy_block"}
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


def _gpu_model() -> str | None:
    try:
        completed = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode:
        return None
    names = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    return ", ".join(names) if names else None


def capture_environment(
    repo: Path,
    build_dir: Path | None = None,
    *,
    omp_threads: int | None = None,
    blas_threads: int | None = None,
    mpi_ranks: int | None = None,
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
    return {
        "git_commit": _git_commit(repo),
        "compiler": compiler,
        "build_type": cache.get("CMAKE_BUILD_TYPE") or "unspecified",
        "blas_lapack": blas,
        "omp_threads": int(omp) if str(omp).isdigit() else omp,
        "blas_threads": int(blas_thread_value) if str(blas_thread_value).isdigit() else blas_thread_value,
        "omp_proc_bind": os.environ.get("OMP_PROC_BIND"),
        "omp_places": os.environ.get("OMP_PLACES"),
        "mpi_ranks": int(ranks) if str(ranks).isdigit() else ranks,
        "cuda_toolkit": cache.get("CUDAToolkit_VERSION") or cache.get("CMAKE_CUDA_COMPILER_VERSION"),
        "cusolver": cache.get("CUDAToolkit_cusolver_VERSION"),
        "gpu_model": _gpu_model(),
        "cpu_model": platform.processor() or platform.machine(),
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
