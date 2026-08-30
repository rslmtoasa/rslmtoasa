#!/usr/bin/env python3
"""Resume-safe SCF-B1R campaign driver.

The driver deliberately keeps orchestration and evidence collection here.  It
does not silently turn a failed build or an unsupported accelerator case into a
timing row.  Long campaign execution is user-controlled; ``--dry-run`` and
``--preflight-only`` are intended for the short Phase-A checks.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import scf_b0c  # noqa: E402


SCHEMA = "rslmto.scf-b1r.v1"
SENTINEL = "B1R_RUN_COMPLETE"
DEFAULT_OUTPUT = ROOT / "results" / "benchmarks" / "scf_b1r"
DEFAULT_BUILD = ROOT / "build-b1r-release-cuda"
BUILD_BINARIES = {
    "scf": "rslmto.x",
    "p0": "ReciprocalAccP0Benchmark",
    "p1b": "ReciprocalAccP1bPhysicalSCF",
}

OPTIMIZATION_RE = re.compile(r"(?<![A-Za-z0-9_])(?:-O(?:0|1|2|3|s|fast|g|z)|/O(?:d|[0-3]))(?![A-Za-z0-9_])")
SOURCE_RE = re.compile(r"(?:^|[/\\])([A-Za-z0-9_+.-]+\.(?:f90|F90|f|F|cu|cpp|c))(?:\s|$)")
REPRESENTATIVE_SOURCES = {
    "scf_driver": ("main.f90", "self.f90", "calculation.f90"),
    "rs_recursion": ("green_block.f90", "recursion_core.f90", "recursion.f90"),
    "hamiltonian": ("hamiltonian_build.f90", "hamiltonian.f90"),
    "reciprocal": ("reciprocal.f90", "reciprocal_bands.f90", "self_reciprocal.f90"),
}

CSV_FIELDS = [
    "case_id",
    "route",
    "material",
    "physical_replication",
    "reciprocal_basis_replication",
    "backend",
    "solver_strategy",
    "numeric_mode",
    "rs_solver",
    "rs_backend",
    "omp_threads",
    "blas_threads",
    "nmat",
    "nk_nominal",
    "nk_unique",
    "vectors",
    "repetitions",
    "warmups",
    "P_eigensolver",
    "T_solver",
    "T_rs_kernel",
    "T_iteration",
    "T_convergence",
    "status",
    "status_reason",
    "correctness_status",
    "correctness_reason",
    "equal_precision_eligible",
    "production_comparison_eligible",
    "ineligibility_reason",
    "potential_identity",
    "starting_state_id",
    "row_id",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=_json_default) + "\n", encoding="utf-8")


def read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return default


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(*args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def command_text(command: Iterable[str]) -> str:
    return shlex.join([str(item) for item in command])


def ensure_output_tree(output: Path) -> None:
    for name in ("PRE", "RS1", "RS2", "K1", "K2", "X1", "raw", "correctness", "iteration_history", "cross_route", "plots"):
        (output / name).mkdir(parents=True, exist_ok=True)


def parse_cmake_cache(build_dir: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    cache = build_dir / "CMakeCache.txt"
    if not cache.is_file():
        return result
    for line in cache.read_text(errors="replace").splitlines():
        if not line or line.startswith("//") or line.startswith("#") or ":" not in line or "=" not in line:
            continue
        left, value = line.split("=", 1)
        key = left.split(":", 1)[0]
        result[key] = value
    return result


def _compile_commands_from_json(build_dir: Path) -> list[str]:
    path = build_dir / "compile_commands.json"
    data = read_json(path, [])
    commands: list[str] = []
    if not isinstance(data, list):
        return commands
    for entry in data:
        if not isinstance(entry, dict):
            continue
        command = entry.get("command")
        if not command and isinstance(entry.get("arguments"), list):
            command = command_text(entry["arguments"])
        if command:
            commands.append(str(command))
    return commands


def _compile_commands_from_ninja(build_dir: Path) -> list[str]:
    try:
        result = subprocess.run(
            ["ninja", "-C", str(build_dir), "-t", "commands"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return []
    if result.returncode != 0:
        return []
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _compile_commands_from_make(build_dir: Path) -> list[str]:
    commands: list[str] = []
    for path in build_dir.glob("**/flags.make"):
        text = path.read_text(errors="replace")
        for key in ("C_FLAGS", "Fortran_FLAGS", "FFLAGS", "CUDA_FLAGS"):
            match = re.search(rf"^{key}\s*=\s*(.*)$", text, flags=re.MULTILINE)
            if match:
                compiler = "gfortran" if "Fortran" in key or key == "FFLAGS" else "cc"
                commands.append(f"{compiler} {match.group(1)} -c {path}")
    return commands


def collect_compile_commands(build_dir: Path) -> list[str]:
    commands = _compile_commands_from_json(build_dir)
    if commands:
        return commands
    commands = _compile_commands_from_ninja(build_dir)
    if commands:
        return commands
    return _compile_commands_from_make(build_dir)


def _source_names(command: str) -> set[str]:
    return {match.group(1).lower() for match in SOURCE_RE.finditer(command)}


def _last_optimization(command: str) -> str | None:
    matches = OPTIMIZATION_RE.findall(command)
    return matches[-1] if matches else None


def build_preflight(build_dir: Path) -> dict[str, Any]:
    cache = parse_cmake_cache(build_dir)
    all_commands = collect_compile_commands(build_dir)
    representative: dict[str, list[str]] = {}
    conclusions: dict[str, str] = {}
    failures: list[str] = []
    for label, sources in REPRESENTATIVE_SOURCES.items():
        selected = []
        for command in all_commands:
            names = _source_names(command)
            if "-c" in command and any(source.lower() in names for source in sources):
                selected.append(command)
        representative[label] = selected[:4]
        if not selected:
            failures.append(f"no actual compile command found for {label} ({', '.join(sources)})")
            continue
        effective = [_last_optimization(command) for command in selected]
        conclusions[label] = ", ".join(sorted({item or "missing" for item in effective}))
        if any(item in {"-O0", "/Od"} for item in effective):
            failures.append(f"{label} has an effective trailing -O0/debug optimization flag")
        if any(item is None for item in effective):
            failures.append(f"{label} has no detectable optimization flag")

    build_type = cache.get("CMAKE_BUILD_TYPE", "unknown")
    compiler = cache.get("CMAKE_Fortran_COMPILER")
    if not compiler:
        for command in all_commands:
            match = re.match(r"\s*([^\s]+)", command)
            if match:
                compiler = match.group(1)
                break
    if not build_dir.is_dir():
        failures.append("build directory does not exist")
    status = "PASS" if not failures else "FAIL"
    if status == "PASS":
        conclusion = "representative actual compile commands have optimized effective flags; no trailing -O0"
    else:
        conclusion = "; ".join(failures)
    return {
        "schema": "rslmto.scf-b1r.build-preflight.v1",
        "timestamp": utc_now(),
        "build_dir": str(build_dir),
        "compiler": compiler,
        "build_type": build_type,
        "generator": cache.get("CMAKE_GENERATOR"),
        "cache_optimization_flags_non_authoritative": cache.get("CMAKE_Fortran_FLAGS_RELEASE"),
        "representative_compile_commands": representative,
        "effective_optimization_by_component": conclusions,
        "effective_optimization_conclusion": conclusion,
        "status": status,
        "reason": None if status == "PASS" else "; ".join(failures),
        "timing_reuse_from_B1": "forbidden",
        "timing_reuse_reason": "B1 Release timings were produced before the B1R -O0 integrity fix and must be rerun.",
        "source_commit": git_value("rev-parse", "HEAD"),
        "source_cmake_flag_fix": "cmake/SetFortranFlags.cmake no longer appends -O0 to Release flags",
    }


def _cache_bool(cache: dict[str, str], key: str, default: bool = False) -> bool:
    value = cache.get(key)
    if value is None:
        return default
    return value.upper() in {"ON", "TRUE", "YES", "1"}


def configure_args_from_cache(cache: dict[str, str]) -> list[str]:
    values = {
        "CMAKE_BUILD_TYPE": cache.get("CMAKE_BUILD_TYPE", "Release"),
        "CMAKE_Fortran_COMPILER": cache.get("CMAKE_Fortran_COMPILER", "gfortran"),
        "ENABLE_CUDA_PLUGIN": _cache_bool(cache, "ENABLE_CUDA_PLUGIN", True),
        "ENABLE_CUDA_RECIPROCAL": _cache_bool(cache, "ENABLE_CUDA_RECIPROCAL", True),
        "ENABLE_FUSED_RECIPROCAL": _cache_bool(cache, "ENABLE_FUSED_RECIPROCAL", True),
        "ENABLE_MPI": _cache_bool(cache, "ENABLE_MPI", False),
        "ENABLE_OPENMP": _cache_bool(cache, "ENABLE_OPENMP", True),
        "ENABLE_SPGLIB": _cache_bool(cache, "ENABLE_SPGLIB", True),
        "ENABLE_MARCH_NATIVE": _cache_bool(cache, "ENABLE_MARCH_NATIVE", True),
        # The production benchmark executables are registered in the same
        # CMake block as the standalone unit targets.  Keep the build target
        # narrow below, but enable this block so the benchmark targets exist.
        "RUN_UNIT_TESTS": True,
        "RUN_REG_TESTS": False,
        "RUN_EXAMPLE_TESTS": False,
    }
    args = [f"-D{key}={'ON' if value is True else 'OFF' if value is False else value}" for key, value in values.items()]
    # Make a reconfigure of an accidentally reused build harmless: the old
    # cache may contain the B1-era trailing -O0 even though the source fix is
    # now correct.  The post-build preflight remains authoritative.
    args.extend(["-DCMAKE_Fortran_FLAGS_RELEASE=-O3", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"])
    for key in ("CMAKE_CUDA_COMPILER", "BLA_VENDOR"):
        if cache.get(key):
            args.append(f"-D{key}={cache[key]}")
    return args


def run_logged_command(
    command: list[str],
    evidence_dir: Path,
    *,
    cwd: Path = ROOT,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    evidence_dir.mkdir(parents=True, exist_ok=True)
    (evidence_dir / "command.txt").write_text(command_text(command) + "\n", encoding="utf-8")
    started = time.time()
    start_iso = utc_now()
    try:
        # Production Fortran output is not guaranteed to be UTF-8 (some
        # legacy diagnostics use an extended single-byte character).  Capture
        # bytes first so logging can never abort before a case is classified.
        process = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=False, check=False)
        stdout_bytes = process.stdout or b""
        stderr_bytes = process.stderr or b""
        stdout = stdout_bytes.decode("utf-8", errors="replace")
        stderr = stderr_bytes.decode("utf-8", errors="replace")
        returncode = process.returncode
    except OSError as exc:
        stdout_bytes = b""
        stderr_bytes = b""
        stdout = ""
        stderr = f"{type(exc).__name__}: {exc}\n"
        returncode = 127
    (evidence_dir / "stdout.bin").write_bytes(stdout_bytes)
    (evidence_dir / "stderr.bin").write_bytes(stderr_bytes)
    (evidence_dir / "stdout.log").write_text(stdout, encoding="utf-8", errors="replace")
    (evidence_dir / "stderr.log").write_text(stderr, encoding="utf-8", errors="replace")
    (evidence_dir / "combined.log").write_text(stdout + ("\n" if stdout and stderr else "") + stderr, encoding="utf-8", errors="replace")
    return {
        "command": [str(item) for item in command],
        "command_text": command_text(command),
        "started": start_iso,
        "finished": utc_now(),
        "wall_seconds": time.time() - started,
        "returncode": returncode,
        "stdout": stdout,
        "stderr": stderr,
    }


def ensure_build(output: Path, requested: Path | None, dry_run: bool) -> tuple[Path, dict[str, Any]]:
    candidate = requested
    if candidate is None:
        if DEFAULT_BUILD.exists():
            candidate = DEFAULT_BUILD
        else:
            candidate = ROOT / "build-b1-frozen-cuda"
    candidate = candidate.resolve()
    initial = build_preflight(candidate)
    binaries_present = all((candidate / "bin" / name).is_file() for name in BUILD_BINARIES.values())
    if initial["status"] == "PASS" and binaries_present:
        return candidate, initial
    if dry_run:
        initial["dry_run"] = True
        initial["reason"] = initial.get("reason") or "dry-run: build was not executed"
        return candidate, initial

    build_dir = candidate if requested is not None else DEFAULT_BUILD
    cache_source = parse_cmake_cache(candidate)
    configure = [
        "cmake",
        "-S",
        str(ROOT),
        "-B",
        str(build_dir),
        "-G",
        cache_source.get("CMAKE_GENERATOR", "Ninja"),
        *configure_args_from_cache(cache_source),
    ]
    build_log = output / "raw" / "build"
    configure_result = run_logged_command(configure, build_log / "configure")
    if configure_result["returncode"] != 0:
        initial["status"] = "FAIL"
        initial["reason"] = "CMake configure failed; see raw/build/configure"
        write_json(output / "build_preflight.json", initial)
        raise RuntimeError(initial["reason"])
    build_command = [
        "cmake",
        "--build",
        str(build_dir),
        "--target",
        *BUILD_BINARIES.values(),
        "--parallel",
    ]
    build_result = run_logged_command(build_command, build_log / "compile")
    if build_result["returncode"] != 0:
        initial["status"] = "FAIL"
        initial["reason"] = "Release build failed; see raw/build/compile"
        write_json(output / "build_preflight.json", initial)
        raise RuntimeError(initial["reason"])
    final = build_preflight(build_dir)
    final["configure_command"] = configure
    final["build_command"] = build_command
    final["previous_candidate"] = str(candidate)
    if not all((build_dir / "bin" / name).is_file() for name in BUILD_BINARIES.values()):
        final["status"] = "FAIL"
        final["reason"] = "one or more required benchmark binaries are missing after build"
    write_json(output / "build_preflight.json", final)
    if final["status"] != "PASS":
        raise RuntimeError(final.get("reason") or "build preflight failed")
    return build_dir, final


def _case_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for omp in (1, 2, 4, 8):
        specs.append({
            "case_id": f"rs1_fe3_cpu_omp{omp}", "kind": "scf", "route": "real_space",
            "material": "fe3", "backend": "cpu", "solver_strategy": "lapack",
            "rs_solver": "block", "rs_backend": "csr", "omp_threads": omp,
            "blas_threads": 1, "nstep": 5, "benchmark_level": "scf_iteration",
            "warmups": 2, "repetitions": 3, "vectors": False,
        })
    specs.append({
        "case_id": "rs1_fe3_cuda_mixed", "kind": "scf", "route": "real_space",
        "material": "fe3", "backend": "cuda", "solver_strategy": "fp64_zheevd",
        "rs_solver": "block", "rs_backend": "csr", "omp_threads": 1,
        "blas_threads": 1, "nstep": 5, "benchmark_level": "scf_iteration",
        "warmups": 2, "repetitions": 3, "vectors": False,
    })
    for length in (1, 2, 3):
        for backend in ("cpu", "cuda"):
            specs.append({
                "case_id": f"k1_fe{length}_{backend}", "kind": "scf", "route": "reciprocal",
                "material": f"fe{length}", "backend": backend, "solver_strategy": "fp64_zheevd",
                "rs_solver": "block", "rs_backend": "csr", "omp_threads": 1,
                "blas_threads": 1, "nstep": 3, "benchmark_level": "scf_iteration",
                "warmups": 1, "repetitions": 1, "vectors": True,
            })
    specs.append({
        "case_id": "rs2_reference_fe_cpu", "kind": "reference", "route": "real_space",
        "material": "fe", "backend": "cpu", "solver_strategy": "lapack",
        "rs_solver": "block", "rs_backend": "csr", "omp_threads": 1,
        "blas_threads": 1, "nstep": 80, "benchmark_level": "scf_convergence",
        "warmups": 0, "repetitions": 1, "vectors": False,
    })
    for backend in ("cpu", "cuda"):
        specs.append({
            "case_id": f"rs2_common_fe_{'cpu_fp64' if backend == 'cpu' else 'cuda_mixed'}",
            "kind": "common_state", "route": "real_space", "material": "fe", "backend": backend,
            "solver_strategy": "lapack" if backend == "cpu" else "fp64_zheevd", "rs_solver": "block",
            "rs_backend": "csr", "omp_threads": 1, "blas_threads": 1, "nstep": 3,
            "benchmark_level": "scf_iteration", "warmups": 1, "repetitions": 1, "vectors": False,
        })
    for length in (3, 4, 5):
        specs.append({
            "case_id": f"k2_fe{length}", "kind": "large_reciprocal", "route": "reciprocal",
            "material": f"fe{length}", "backend": "mixed", "solver_strategy": "fp64_zheevd",
            "warmups": 2, "repetitions": 3, "vectors": True,
        })
    for route in ("real_space", "reciprocal"):
        specs.append({
            "case_id": f"x1_fe_cpu_{'rs' if route == 'real_space' else 'kspace'}",
            "kind": "cross_route", "route": route, "material": "fe", "backend": "cpu",
            "solver_strategy": "lapack", "rs_solver": "block", "rs_backend": "csr",
            "omp_threads": 1, "blas_threads": 1, "nstep": 5, "benchmark_level": "scf_iteration",
            "warmups": 1, "repetitions": 1, "vectors": route == "reciprocal",
        })
    matrix_by_kind = {
        "scf": "RS1",
        "reference": "RS2",
        "common_state": "RS2",
        "large_reciprocal": "K2",
        "cross_route": "X1",
    }
    for spec in specs:
        spec["matrix"] = matrix_by_kind[spec["kind"]]
        if spec["case_id"].startswith("k1_"):
            spec["matrix"] = "K1"
    return specs


def _spec_by_id() -> dict[str, dict[str, Any]]:
    return {spec["case_id"]: spec for spec in _case_specs()}


def _material_fixture(material: str) -> dict[str, Any]:
    if material in scf_b0c.FIXTURES:
        return dict(scf_b0c.FIXTURES[material])
    if material.startswith("fe") and material[2:].isdigit():
        length = int(material[2:])
        if length == 1:
            return dict(scf_b0c.FIXTURES["fe"])
        return {
            "material": "bccFe",
            "source": ROOT / "example" / "bulk" / "supercellFe" / f"{length}x{length}x{length}",
            "supercell": f"{length}x{length}x{length}",
            "fixture_id": "bccFe_reciprocal_scf_supercell",
            "nsp": 2, "soc": "off", "basis": "spd",
        }
    raise ValueError(f"unknown B1R material fixture {material!r}")


def _restart_copy(source: Path | None, fixture_dir: Path) -> str | None:
    if source is None or not source.is_file():
        return None
    target = fixture_dir / "Fe_out.nml"
    shutil.copy2(source, target)
    return sha256_file(source)


def _row_status(parsed: dict[str, Any], returncode: int, level: str) -> tuple[str, str | None]:
    if parsed.get("unsupported_reason"):
        return "UNSUPPORTED", parsed.get("unsupported_reason") or "benchmark reported unsupported"
    if returncode != 0:
        return "FAIL", f"command exited with status {returncode}"
    if not parsed.get("result"):
        return "FAIL", "no SCF_B0C_RESULT record was found"
    if parsed.get("profile", {}).get("status") != "PASS":
        return "FAIL", parsed.get("profile", {}).get("reason") or "SCF profile closure failed"
    if level == "scf_convergence" and (parsed.get("result") or {}).get("converged") not in (True, "true", "TRUE"):
        return "NOT_CONVERGED", "reference run did not reach the requested convergence criterion"
    return "PASS", None


def _make_probe_command(binary: Path, fixture_dir: Path, spec: dict[str, Any]) -> list[str]:
    return [
        str(binary),
        "--input", str(fixture_dir / "input.nml"),
        "--backend", "cuda" if spec["backend"] == "cuda" else "lapack",
        "--solver-strategy", spec.get("solver_strategy", "lapack"),
        "--dos-method", "tetrahedron", "--scf-route", spec["route"],
        "--rs-solver", spec.get("rs_solver", "block"),
        "--rs-backend", spec.get("rs_backend", "csr"),
        "--sigma", "0.02", "--temperature", "300", "--nstep", str(spec["nstep"]),
        "--benchmark-level", spec["benchmark_level"], "--profile",
    ]


def _run_one_scf(
    binary: Path,
    spec: dict[str, Any],
    run_dir: Path,
    restart: Path | None = None,
    potential_identity: str | None = None,
    starting_state_id: str | None = None,
    reference_final: dict[str, Any] | None = None,
) -> dict[str, Any]:
    fixture_dir = run_dir / "fixture"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    fixture = _material_fixture(spec["material"])
    shutil.copytree(Path(fixture["source"]), fixture_dir, dirs_exist_ok=True)
    restart_hash = _restart_copy(restart, fixture_dir)
    env = scf_b0c._controlled_environment(spec.get("omp_threads", 1), spec.get("blas_threads", 1))
    env["RSLMTO_SCF_B0C_PROFILE"] = "1"
    command = _make_probe_command(binary, fixture_dir, spec)
    execution = run_logged_command(command, run_dir / "evidence", cwd=fixture_dir, env=env)
    parsed = scf_b0c.parse_probe_output(execution["stdout"] + execution["stderr"])
    result = parsed.get("result") or {}
    status, reason = _row_status(parsed, execution["returncode"], spec["benchmark_level"])
    raw = {
        "backend": "cuda" if spec["backend"] == "cuda" else "lapack",
        "solver_strategy": spec.get("solver_strategy"),
        "dos_method": "tetrahedron",
        "OMP_threads": spec.get("omp_threads"),
        "BLAS_threads": spec.get("blas_threads"),
        "cold_process_wall": execution["wall_seconds"],
        "full_scf_wall": execution["wall_seconds"],
        "returncode": execution["returncode"],
        "profile_status": parsed.get("profile", {}).get("status", "FAIL"),
        "profile": parsed.get("profile", {}),
        "iterations": parsed.get("iterations", []),
        "final_state": result,
        "final_state_raw": result,
        "status": status,
        "unsupported_reason": parsed.get("unsupported_reason"),
        "rs_solver": spec.get("rs_solver"),
    }
    if parsed.get("iterations"):
        raw.update(scf_b0c._profile_summary(parsed["iterations"]))
    row = scf_b0c._decorate_row(
        raw,
        fixture,
        name=spec["material"],
        nstep=spec["nstep"],
        sigma=0.02,
        temperature=300.0,
        benchmark_level=spec["benchmark_level"],
        potential_identity=potential_identity or "normal_initial",
        scf_route=spec["route"],
        rs_backend=spec.get("rs_backend", "csr"),
    )
    row.update({
        "case_id": spec["case_id"],
        "status": status,
        "status_reason": reason,
        "failure_ineligibility_reason": reason,
        "physical_replication": fixture.get("supercell"),
        "reciprocal_basis_replication": "runtime_authoritative",
        "reciprocal_basis_replication_source": "SCF_B0C_RESULT:nmat",
        "vectors": bool(spec.get("vectors")),
        "omp_threads": spec.get("omp_threads"),
        "blas_threads": spec.get("blas_threads"),
        "restart_identity": restart_hash or potential_identity,
        "starting_state_id": starting_state_id,
        "equal_precision_eligible": False,
        "production_comparison_eligible": False,
    })
    if status == "PASS" and reference_final:
        candidate_for_check = dict(result)
        # The fixed-state RS anchor intentionally measures only a short
        # production segment; compare observables without requiring another
        # full convergence cycle.
        candidate_for_check["converged"] = True
        correctness = scf_b0c._correctness(reference_final, candidate_for_check, row["numeric_mode"], require_convergence=True)
        correctness["reference_source"] = "rs2_reference_fe_cpu"
        if not correctness.get("reason"):
            if correctness.get("status") == "PASS":
                correctness["reason"] = "reference observable comparison passed"
            else:
                failed_fields = correctness.get("failed_fields") or []
                correctness["reason"] = ", ".join(str(field) for field in failed_fields) or str(correctness.get("status", "comparison unavailable"))
    else:
        correctness = {
            "status": "NOT_APPLICABLE" if status == "PASS" else status,
            "reason": "common-state or route comparison is required" if status == "PASS" else reason,
        }
    row["correctness"] = correctness
    row["correctness_status"] = correctness["status"]
    row["correctness_reason"] = correctness.get("reason", correctness.get("status", "comparison unavailable"))
    row["timer_validation"] = {
        "status": "PENDING",
        "nmat_source": "runtime SCF_B0C_RESULT",
        "nk_unique_source": "runtime SCF_B0C_RESULT",
        "cpu_timer": row.get("P_eigensolver") if row.get("backend") == "lapack" else None,
        "gpu_timer": row.get("T_solver") if row.get("backend") == "cuda" else None,
    }
    try:
        row["evidence_dir"] = str(run_dir.relative_to(ROOT))
    except ValueError:
        row["evidence_dir"] = str(run_dir)
    row["row_id"] = hashlib.sha256(json.dumps(row, sort_keys=True, default=_json_default).encode()).hexdigest()[:16]
    return row


def _aggregate_rows(rows: list[dict[str, Any]], spec: dict[str, Any]) -> dict[str, Any]:
    if not rows:
        return {
            "case_id": spec["case_id"], "status": "SKIPPED",
            "status_reason": "no measured repetitions completed", "repetitions": 0,
            "warmups": spec.get("warmups", 0),
        }
    successful = [row for row in rows if row.get("status") == "PASS"]
    pool = successful or rows
    timing_fields = ("T_iteration", "T_convergence", "T_rs_kernel", "P_eigensolver", "T_solver")
    representative = copy.deepcopy(pool[len(pool) // 2])
    for field in timing_fields:
        values = [row[field] for row in pool if isinstance(row.get(field), (int, float))]
        if values:
            representative[field] = statistics.median(values)
    representative["repetitions"] = len(rows)
    representative["warmups"] = spec.get("warmups", 0)
    representative["repetition_rows"] = rows
    if len(successful) != len(rows):
        statuses = {row.get("status") for row in rows}
        if "FAIL" in statuses:
            aggregate_status = "FAIL"
        elif "NOT_CONVERGED" in statuses:
            aggregate_status = "NOT_CONVERGED"
        elif "UNSUPPORTED" in statuses:
            aggregate_status = "UNSUPPORTED"
        else:
            aggregate_status = "SKIPPED"
        representative["status"] = aggregate_status
        representative["status_reason"] = "; ".join(sorted({row.get("status_reason") or row.get("status", "") for row in rows if row.get("status") != "PASS"}))
        representative["failure_ineligibility_reason"] = representative["status_reason"]
    representative["row_id"] = hashlib.sha256(json.dumps(representative, sort_keys=True, default=_json_default).encode()).hexdigest()[:16]
    return representative


def _case_state_path(output: Path, case_id: str) -> Path:
    return output / "raw" / case_id / "case.json"


def _state_is_resumable(state: Any) -> bool:
    if not isinstance(state, dict) or not state.get("completed"):
        return False
    aggregate = state.get("aggregate", {})
    status = aggregate.get("status")
    if status in {"PASS", "UNSUPPORTED", "SKIPPED"}:
        return True
    # A completed probe can be scientifically ineligible (for example, the
    # strict profile-misc gate) while still being valid evidence.  Do not
    # rerun that row forever; distinguish it from an old command/launch
    # failure by requiring successful process exits and a parsed final state.
    if status in {"FAIL", "NOT_CONVERGED"}:
        repetitions = aggregate.get("repetition_rows", [])
        return bool(repetitions) and all(
            row.get("returncode") == 0 and isinstance(row.get("final_state"), dict) and row.get("final_state")
            for row in repetitions
        )
    return False


def _load_completed_case(output: Path, case_id: str, force: bool) -> dict[str, Any] | None:
    if force:
        return None
    state = read_json(_case_state_path(output, case_id))
    # Failed/non-converged evidence must be retried after a runner fix.  Only
    # successful or intentional terminal outcomes are resumable by default.
    if _state_is_resumable(state):
        return state
    return None


def run_scf_case(
    output: Path,
    binary: Path,
    spec: dict[str, Any],
    *,
    force: bool = False,
    restart: Path | None = None,
    potential_identity: str | None = None,
    starting_state_id: str | None = None,
    reference_final: dict[str, Any] | None = None,
) -> dict[str, Any]:
    prior_state = read_json(_case_state_path(output, spec["case_id"]))
    if not force and isinstance(prior_state, dict) and prior_state.get("completed") and not _state_is_resumable(prior_state):
        # Do not combine repaired runner logic with stale per-repetition rows.
        force = True
    cached = _load_completed_case(output, spec["case_id"], force)
    if cached:
        return cached
    case_dir = output / "raw" / spec["case_id"]
    case_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for index in range(int(spec.get("warmups", 0))):
        run_dir = case_dir / f"warmup-{index + 1:02d}"
        if not (run_dir / "evidence" / "combined.log").exists() or force:
            _run_one_scf(binary, spec, run_dir, restart, potential_identity, starting_state_id, reference_final)
    for index in range(int(spec.get("repetitions", 1))):
        run_dir = case_dir / f"rep-{index + 1:02d}"
        row_path = run_dir / "row.json"
        if row_path.is_file() and not force:
            row = read_json(row_path, {})
        else:
            row = _run_one_scf(binary, spec, run_dir, restart, potential_identity, starting_state_id, reference_final)
            write_json(row_path, row)
        if row:
            rows.append(row)
    aggregate = _aggregate_rows(rows, spec)
    state = {
        "schema": "rslmto.scf-b1r.case.v1",
        "case_id": spec["case_id"],
        "completed": True,
        "finished": utc_now(),
        "spec": spec,
        "aggregate": aggregate,
    }
    write_json(_case_state_path(output, spec["case_id"]), state)
    return state


def run_reference(output: Path, binary: Path, spec: dict[str, Any], force: bool) -> tuple[dict[str, Any], Path | None]:
    state = run_scf_case(output, binary, spec, force=force)
    reference_fixture = output / "raw" / spec["case_id"] / "rep-01" / "fixture"
    candidates = [reference_fixture / "Fe_scf.nml", reference_fixture / "Fe_out.nml"]
    restart = next((path for path in candidates if path.is_file()), None)
    return state, restart


def reference_final_state(state: dict[str, Any]) -> dict[str, Any] | None:
    aggregate = state.get("aggregate", {})
    final = aggregate.get("final_state")
    if isinstance(final, dict):
        return final
    return None


def run_large_case(output: Path, binary: Path, spec: dict[str, Any], force: bool) -> dict[str, Any]:
    prior_state = read_json(_case_state_path(output, spec["case_id"]))
    if not force and isinstance(prior_state, dict) and prior_state.get("completed") and not _state_is_resumable(prior_state):
        force = True
    cached = _load_completed_case(output, spec["case_id"], force)
    if cached:
        return cached
    case_dir = output / "raw" / spec["case_id"]
    p0_output = case_dir / "p0"
    length = int(spec["material"][2:])
    command = [
        sys.executable,
        str(ROOT / "tests" / "benchmarks" / "accp0_real_material.py"),
        "--binary", str(binary), "--build-dir", str(binary.parent.parent),
        "--output-dir", str(p0_output), "--fe-lengths", str(length), "--meshes", "1",
        "--tiles", "1", "--vectors", "--warmups", "2", "--repetitions", "3",
        "--cuda-strategies", "fp64_zheevd",
    ]
    execution = run_logged_command(command, case_dir / "evidence")
    result_file = p0_output / "accp0_results.json"
    payload = read_json(result_file)
    rows: list[dict[str, Any]] = []
    if isinstance(payload, dict):
        all_rows = payload.get("rows", []) if isinstance(payload.get("rows"), list) else []
        # Keep the canonical B1R case focused on the requested Fe length and
        # vectors=yes production workload; the child P0 driver also emits Si,
        # values-only, and matched-density rows for its own validation.
        rows = [
            row for row in all_rows
            if row.get("fixture") == "bccFe"
            and int(row.get("L", -1)) == length
            and int(row.get("vectors", 0)) == 1
            and row.get("workload") == "crossover"
        ]
    if rows:
        for row in rows:
            row["case_id"] = spec["case_id"]
            row["status"] = row.get("status", "PASS")
            row["status_reason"] = row.get("status_reason")
            row["physical_replication"] = f"{length}x{length}x{length}"
            row["reciprocal_basis_replication"] = "runtime_authoritative"
            row["reciprocal_basis_replication_source"] = "P0 runtime metadata"
            row["vectors"] = True
            row["equal_precision_eligible"] = row.get("backend") in {"cpu", "cuda"}
            row["production_comparison_eligible"] = row.get("backend") in {"cpu", "cuda"}
        status = "PASS"
        reason = None
    else:
        status = "SKIPPED" if length == 5 or execution["returncode"] != 0 else "FAIL"
        tail = (execution["stderr"] or execution["stdout"])[-1000:].strip()
        reason = f"P0 matrix produced no rows (exit {execution['returncode']}): {tail}"
    aggregate = {
        "case_id": spec["case_id"], "completed": True, "status": status,
        "status_reason": reason, "rows": rows, "command": execution["command"],
        "repetitions": 3, "warmups": 2,
    }
    state = {"schema": "rslmto.scf-b1r.case.v1", "case_id": spec["case_id"], "completed": True, "finished": utc_now(), "spec": spec, "aggregate": aggregate}
    write_json(_case_state_path(output, spec["case_id"]), state)
    return state


def _flatten_rows(states: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for state in states:
        aggregate = state.get("aggregate", {})
        if isinstance(aggregate.get("rows"), list):
            rows.extend(aggregate["rows"])
        elif aggregate:
            rows.append(aggregate)
    return rows


def _annotate_pairing(rows: list[dict[str, Any]]) -> None:
    for row in rows:
        # ACC-P0 rows use a slightly different field vocabulary. Normalize
        # them before writing the unified B1R campaign so K2 remains visible
        # in JSON/CSV/Markdown without pretending it is a full SCF row.
        if str(row.get("case_id", "")).startswith("k2_"):
            row.setdefault("scf_route", "reciprocal")
            row.setdefault("material", row.get("fixture", "bccFe"))
            row.setdefault("supercell", f"{row.get('L')}x{row.get('L')}x{row.get('L')}")
            row.setdefault("physical_replication", row.get("supercell"))
            row.setdefault("Nk_unique", row.get("actual_unique_nk"))
            row.setdefault("numeric_mode", "fp64")
            row["reciprocal_basis_replication"] = f"runtime:nmat={row.get('nmat')}"
        elif row.get("scf_route") == "reciprocal":
            row["reciprocal_basis_replication"] = f"runtime:nmat={row.get('nmat')}"
        elif row.get("scf_route") == "real_space":
            row["reciprocal_basis_replication"] = "not_applicable"
    try:
        scf_b0c.add_pairings(rows)
    except (KeyError, TypeError, ValueError):
        # A partially completed campaign must remain reportable.  The report
        # generator will show the explicit reason rather than inventing a pair.
        pass
    for row in rows:
        route = row.get("scf_route")
        backend = row.get("backend")
        # Equal precision is a property of the actual route and numeric mode,
        # not of whether the implementation happens to run on CUDA. RS CUDA
        # remains mixed precision; reciprocal FP64 CPU/GPU rows are comparable
        # at this gate.
        equal_precision = (
            row.get("status") == "PASS"
            and row.get("numeric_mode") == "fp64"
            and backend in {"lapack", "cpu", "cuda"}
            and not (route == "real_space" and backend == "cuda")
        )
        if row.get("status") != "PASS":
            row["equal_precision_eligible"] = False
            row["production_comparison_eligible"] = False
            row.setdefault("failure_ineligibility_reason", row.get("status_reason"))
        elif not row.get("equal_precision_eligible"):
            row["equal_precision_eligible"] = equal_precision
            if backend in {"lapack", "cpu"}:
                row["production_comparison_eligible"] = True
        row.setdefault("equal_precision_eligible", False)
        row.setdefault("production_comparison_eligible", False)
        row.setdefault("correctness_status", row.get("correctness", {}).get("status", "NOT_APPLICABLE"))
        row.setdefault("correctness_reason", row.get("correctness", {}).get("reason"))
        rejection = row.get("headline_rejection_reasons") or []
        row["ineligibility_reason"] = ", ".join(str(value) for value in rejection) or row.get("failure_ineligibility_reason") or row.get("correctness_reason")
        if row.get("headline_speedup_eligible") is not True:
            for field in ("S_solver", "S_reciprocal", "S_rs_kernel", "S_rs_phase", "S_iteration", "S_convergence"):
                row.pop(field, None)
        if row.get("production_comparison_eligible") is not True:
            for field in ("R_rs_kernel_production", "R_rs_phase_production", "R_iteration_production", "R_convergence_production"):
                row.pop(field, None)


def _large_reciprocal_evidence(states: list[dict[str, Any]]) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    for state in states:
        aggregate = state.get("aggregate", {})
        rows = aggregate.get("rows")
        if not isinstance(rows, list) or not rows:
            continue
        spec = state.get("spec", {})
        if spec.get("kind") != "large_reciprocal":
            continue
        evidence.append({
            "case_id": state.get("case_id"),
            "schema": "rslmto.scf-b1r.large-reciprocal.v1",
            "status": aggregate.get("status"),
            "reason": aggregate.get("status_reason"),
            "rows": rows,
        })
    return evidence


def validate_reciprocal_timers(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Check that CPU and CUDA reciprocal timers are present at every size.

    This is a semantic instrumentation check, not a claim that CPU and GPU
    timer names are interchangeable: CPU uses the outer ``P_eigensolver``
    boundary and CUDA uses nested device ``T_solver`` work.
    """

    checks: list[dict[str, Any]] = []
    for row in rows:
        if row.get("scf_route") != "reciprocal" or not str(row.get("case_id", "")).startswith("k1_"):
            continue
        field = "T_solver" if row.get("backend") == "cuda" else "P_eigensolver"
        timer = row.get(field)
        status = "PASS" if row.get("status") == "PASS" and isinstance(timer, (int, float)) and timer > 0 else "INCONCLUSIVE"
        reason = None if status == "PASS" else f"{field} missing/non-positive or run status is {row.get('status')}"
        row.setdefault("timer_validation", {}).update({
            "status": status,
            "authoritative_timer": field,
            "reason": reason,
            "runtime_nmat": row.get("nmat"),
            "runtime_nk_unique": row.get("Nk_unique"),
        })
        checks.append({
            "case_id": row.get("case_id"), "backend": row.get("backend"),
            "nmat": row.get("nmat"), "Nk_unique": row.get("Nk_unique"),
            "timer": field, "value": timer, "status": status, "reason": reason,
        })
    return checks


def write_campaign_outputs(output: Path, states: list[dict[str, Any]], preflight: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    rows = _flatten_rows(states)
    _annotate_pairing(rows)
    timer_validation = validate_reciprocal_timers(rows)
    campaign = {
        "schema": SCHEMA,
        "created": manifest.get("created"),
        "updated": utc_now(),
        "build_preflight": preflight,
        "provenance": {
            "commit": git_value("rev-parse", "HEAD"),
            "dirty": bool(git_value("status", "--porcelain")),
            "host": os.uname().nodename if hasattr(os, "uname") else None,
        },
        "policy": {
            "performance_metric": "steady production SCF iteration wall time",
            "cycle_count_role": "correctness and stability only",
            "ineligible_speedups_are_suppressed": True,
        },
        "reciprocal_timer_validation": timer_validation,
        "large_reciprocal_evidence": _large_reciprocal_evidence(states),
        "manifest": manifest,
        "rows": rows,
    }
    write_json(output / "campaign.json", campaign)
    for row in rows:
        row_id = str(row.get("row_id") or row.get("case_id"))
        write_json(output / "correctness" / f"{row_id}.json", {
            "row_id": row_id,
            "status": row.get("status"),
            "reason": row.get("status_reason"),
            "correctness": row.get("correctness"),
            "timer_validation": row.get("timer_validation"),
        })
        write_json(output / "iteration_history" / f"{row_id}.json", {
            "row_id": row_id,
            "iterations": row.get("iterations", []),
            "repetition_rows": row.get("repetition_rows", []),
        })
    fields = CSV_FIELDS
    with (output / "campaign.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            values = {field: row.get(field) for field in fields}
            values["route"] = row.get("route") or row.get("scf_route") or ("reciprocal" if str(row.get("case_id", "")).startswith("k2_") else None)
            values["material"] = row.get("material") or row.get("fixture")
            values["physical_replication"] = row.get("physical_replication") or row.get("supercell")
            values["reciprocal_basis_replication"] = row.get("reciprocal_basis_replication") or (f"runtime:nmat={row.get('nmat')}" if values["route"] == "reciprocal" else "not_applicable")
            values["nk_unique"] = row.get("nk_unique") if row.get("nk_unique") is not None else row.get("Nk_unique", row.get("actual_unique_nk"))
            writer.writerow(values)
    lines = [
        "# SCF-B1R campaign evidence",
        "",
        f"Schema: `{SCHEMA}`  ",
        f"Build preflight: **{preflight.get('status')}** — {preflight.get('effective_optimization_conclusion')}  ",
        "",
        "| case | route | material | backend/solver | nmat | Nk_unique | status | equal precision | reason |",
        "|---|---|---|---|---:|---:|---|---|---|",
    ]
    for row in rows:
        solver = f"{row.get('backend', '')}/{row.get('solver_strategy', '')}"
        row_reason = (
            row.get("status_reason")
            or row.get("failure_ineligibility_reason")
            or ", ".join(str(value) for value in row.get("headline_rejection_reasons", []))
            or row.get("correctness_reason")
            or row.get("correctness", {}).get("reason")
        )
        lines.append(
            f"| {row.get('case_id', '')} | {row.get('scf_route', '') or ('reciprocal' if str(row.get('case_id', '')).startswith('k2_') else '')} | {row.get('material', '') or row.get('fixture', '')} | {solver} | "
            f"{row.get('nmat', '')} | {row.get('Nk_unique', row.get('actual_unique_nk', ''))} | {row.get('status', '')} | "
            f"{row.get('equal_precision_eligible', False)} | {row_reason or ''} |"
        )
    (output / "campaign.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return campaign


def update_manifest(output: Path, manifest: dict[str, Any], states: list[dict[str, Any]]) -> None:
    by_id = {state.get("case_id"): state for state in states}
    for entry in manifest["cases"]:
        state = by_id.get(entry["case_id"])
        if state:
            aggregate = state.get("aggregate", {})
            entry.update({
                "completed": bool(state.get("completed")),
                "status": aggregate.get("status"),
                "status_reason": aggregate.get("status_reason"),
            })
    manifest["updated"] = utc_now()
    write_json(output / "manifest.json", manifest)


def collect_completed_states(output: Path, specs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    states: list[dict[str, Any]] = []
    for spec in specs:
        state = read_json(_case_state_path(output, spec["case_id"]))
        if isinstance(state, dict) and state.get("completed"):
            states.append(state)
    return states


def cross_route_summary(output: Path, rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_case = {row.get("case_id"): row for row in rows}
    rs = by_case.get("x1_fe_cpu_rs")
    reciprocal = by_case.get("x1_fe_cpu_kspace")
    result: dict[str, Any] = {
        "schema": "rslmto.scf-b1r.cross-route.v1",
        "status": "INCONCLUSIVE",
        "reason": "common-potential route comparison requires both successful fixed-state rows",
        "timing_eligible": False,
    }
    if rs and reciprocal and rs.get("status") == reciprocal.get("status") == "PASS":
        metrics = {}
        aliases = {
            "fermi": "fermi_energy", "charge": "total_charge",
            "moment": "site_moment", "total_energy": "final_total_energy",
        }
        left_state = rs.get("final_state", {})
        right_state = reciprocal.get("final_state", {})
        natom = max(float(rs.get("Natom") or reciprocal.get("Natom") or 1), 1.0)
        for key, source_key in aliases.items():
            left = left_state.get(source_key)
            right = right_state.get(source_key)
            if isinstance(left, (int, float)) and isinstance(right, (int, float)):
                if key in {"charge", "total_energy"}:
                    left /= natom
                    right /= natom
                metrics[key] = {"real_space": left, "reciprocal": right, "absolute_difference": abs(left - right)}
        limits = {"total_energy": 1.0e-4, "fermi": 5.0e-5, "moment": 5.0e-4, "charge": 5.0e-5}
        failed = [key for key, limit in limits.items() if key in metrics and metrics[key]["absolute_difference"] > limit]
        missing = sorted(set(limits) - set(metrics))
        result.update({"metrics": metrics, "limits": limits})
        if failed:
            result.update({"status": "FAIL", "reason": f"accuracy tolerance exceeded: {', '.join(failed)}"})
        elif missing:
            result["reason"] = f"missing accuracy metrics: {', '.join(missing)}"
        else:
            result.update({"status": "PASS", "reason": "common-potential CPU FP64 route comparison passed"})
            result["timing_eligible"] = True
    write_json(output / "cross_route" / "x1.json", result)
    return result


def planned_command(spec: dict[str, Any], binary: Path, output: Path) -> list[str]:
    if spec["kind"] == "large_reciprocal":
        return [sys.executable, str(ROOT / "tests" / "benchmarks" / "accp0_real_material.py"), "--binary", str(binary), "--output-dir", str(output / "raw" / spec["case_id"] / "p0")]
    return _make_probe_command(binary, output / "raw" / spec["case_id"] / "fixture", spec)


def run_campaign(args: argparse.Namespace) -> int:
    output = Path(args.output).resolve()
    ensure_output_tree(output)
    all_specs = _case_specs()
    specs = all_specs
    if args.case:
        selected = set(args.case)
        specs = [spec for spec in specs if spec["case_id"] in selected]
        missing = selected - {spec["case_id"] for spec in specs}
        if missing:
            raise ValueError(f"unknown B1R case(s): {', '.join(sorted(missing))}")
    planned_build = Path(args.build_dir).resolve() if args.build_dir else DEFAULT_BUILD
    planned_binary = Path(args.binary).resolve() if args.binary else planned_build / "bin" / BUILD_BINARIES["p1b"]
    manifest = read_json(output / "manifest.json", {})
    if not isinstance(manifest, dict) or manifest.get("schema") != "rslmto.scf-b1r.manifest.v1":
        manifest = {
            "schema": "rslmto.scf-b1r.manifest.v1", "created": utc_now(), "updated": utc_now(),
            "repo_root": str(ROOT), "command": [sys.executable, *sys.argv],
            "matrices": ["PRE", "RS1", "RS2", "K1", "K2", "X1"],
            "cases": [{"case_id": spec["case_id"], "matrix": spec["matrix"], "kind": spec["kind"], "planned_command": planned_command(spec, planned_binary, output)} for spec in specs],
            "sentinel": SENTINEL,
        }
    else:
        known = {entry["case_id"] for entry in manifest.get("cases", [])}
        for spec in specs:
            if spec["case_id"] not in known:
                manifest.setdefault("cases", []).append({"case_id": spec["case_id"], "matrix": spec["matrix"], "kind": spec["kind"], "planned_command": planned_command(spec, planned_binary, output)})
    write_json(output / "manifest.json", manifest)

    if args.aggregate_only:
        preflight = read_json(output / "build_preflight.json", {})
        states = collect_completed_states(output, all_specs)
        missing = sorted(set(spec["case_id"] for spec in all_specs) - {state.get("case_id") for state in states})
        if missing:
            raise RuntimeError(f"cannot aggregate incomplete B1R campaign; missing: {', '.join(missing)}")
        campaign = write_campaign_outputs(output, states, preflight, manifest)
        campaign["cross_route"] = cross_route_summary(output, campaign["rows"])
        write_json(output / "campaign.json", campaign)
        manifest["completed"] = True
        manifest["completed_at"] = manifest.get("completed_at") or utc_now()
        update_manifest(output, manifest, states)
        (output / SENTINEL).write_text(f"{utc_now()} {SENTINEL}\n", encoding="utf-8")
        print(f"B1R aggregate output: {output}")
        print(SENTINEL)
        return 0

    build_dir, preflight = ensure_build(output, Path(args.build_dir).resolve() if args.build_dir else None, args.dry_run)
    write_json(output / "build_preflight.json", preflight)
    write_json(output / "PRE" / "build_preflight.json", preflight)
    # The B1R SCF probe CLI is implemented by the physical SCF benchmark
    # executable, not by the ordinary production driver.
    binary = Path(args.binary).resolve() if args.binary else build_dir / "bin" / BUILD_BINARIES["p1b"]
    p0_binary = build_dir / "bin" / BUILD_BINARIES["p0"]
    if args.dry_run:
        print(f"B1R dry-run output: {output}")
        print(f"Build preflight: {preflight.get('status')} — {preflight.get('reason') or preflight.get('effective_optimization_conclusion')}")
        for spec in specs:
            print(f"PLAN {spec['case_id']}: {command_text(planned_command(spec, binary, output))}")
        return 0
    if preflight.get("status") != "PASS":
        raise RuntimeError(f"build preflight failed: {preflight.get('reason')}")
    if not binary.is_file():
        raise RuntimeError(f"SCF binary does not exist: {binary}")

    states: list[dict[str, Any]] = []
    reference_restart: Path | None = None
    reference_identity: str | None = None
    reference_final: dict[str, Any] | None = None
    rerun = args.force or not args.resume
    for spec in specs:
        print(f"B1R case {spec['case_id']}", flush=True)
        if spec["kind"] == "reference":
            state, reference_restart = run_reference(output, binary, spec, rerun)
            reference_final = reference_final_state(state)
            if reference_restart:
                reference_identity = sha256_file(reference_restart)
        elif spec["kind"] == "common_state":
            state = run_scf_case(output, binary, spec, force=rerun, restart=reference_restart, potential_identity=reference_identity, starting_state_id="rs2_reference_fe_cpu", reference_final=reference_final)
        elif spec["kind"] == "cross_route":
            state = run_scf_case(output, binary, spec, force=rerun, restart=reference_restart, potential_identity=reference_identity, starting_state_id="rs2_reference_fe_cpu", reference_final=reference_final)
        elif spec["kind"] == "large_reciprocal":
            if not p0_binary.is_file():
                state = {"case_id": spec["case_id"], "completed": True, "aggregate": {"case_id": spec["case_id"], "status": "SKIPPED", "status_reason": f"missing binary: {p0_binary}"}}
                write_json(_case_state_path(output, spec["case_id"]), state)
            else:
                state = run_large_case(output, p0_binary, spec, rerun)
        else:
            state = run_scf_case(output, binary, spec, force=rerun)
        states.append(state)
        persisted_states = collect_completed_states(output, all_specs)
        update_manifest(output, manifest, persisted_states)
        write_campaign_outputs(output, persisted_states, preflight, manifest)

    states = collect_completed_states(output, all_specs)
    campaign = write_campaign_outputs(output, states, preflight, manifest)
    campaign["cross_route"] = cross_route_summary(output, campaign["rows"])
    write_json(output / "campaign.json", campaign)
    manifest["completed"] = True
    manifest["completed_at"] = utc_now()
    update_manifest(output, manifest, states)
    (output / SENTINEL).write_text(f"{utc_now()} {SENTINEL}\n", encoding="utf-8")
    print(f"B1R output: {output}")
    print(SENTINEL)
    return 0


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="run the complete B1R matrix")
    parser.add_argument("--case", action="append", help="run only a named case; may be repeated")
    parser.add_argument("--binary", help="SCF binary; defaults to the validated build")
    parser.add_argument("--build-dir", help="build directory to preflight/use")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--force", action="store_true", help="rerun completed cases")
    resume = parser.add_mutually_exclusive_group()
    resume.add_argument("--resume", dest="resume", action="store_true", default=True, help="reuse valid completed case state (default)")
    resume.add_argument("--no-resume", dest="resume", action="store_false", help="rerun every selected case")
    parser.add_argument("--dry-run", action="store_true", help="print the plan without configuring, building, or timing")
    parser.add_argument("--preflight-only", action="store_true", help="run build preflight and stop")
    parser.add_argument("--aggregate-only", action="store_true", help="rebuild canonical outputs from completed evidence without running benchmarks")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    if not args.all and not args.case and not args.dry_run and not args.preflight_only and not args.aggregate_only:
        args.all = True
    if args.preflight_only:
        output = Path(args.output).resolve()
        ensure_output_tree(output)
        build_dir, preflight = ensure_build(output, Path(args.build_dir).resolve() if args.build_dir else None, True)
        write_json(output / "build_preflight.json", preflight)
        write_json(output / "PRE" / "build_preflight.json", preflight)
        print(json.dumps(preflight, indent=2, sort_keys=True))
        return 0 if preflight.get("status") == "PASS" else 1
    try:
        return run_campaign(args)
    except (RuntimeError, ValueError) as exc:
        print(f"B1R_ABORT: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
