#!/usr/bin/env python3
"""Resume-safe Phase-A/Phase-B driver for the SCF-B1R2 campaign.

The driver is intentionally a campaign layer over the production
``ReciprocalAccP1bPhysicalSCF`` probe and the shared SCF-B0C parser.  It does
not change the production algorithm or silently replace a failed CUDA route
with a CPU measurement.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
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


SCHEMA = "rslmto.scf-b1r2.v2"
MANIFEST_SCHEMA = "rslmto.scf-b1r2.manifest.v2"
SENTINEL = "B1R2_RUN_COMPLETE"
DEFAULT_OUTPUT = ROOT / "results" / "benchmarks" / "scf_b1r2"
DEFAULT_BUILD = ROOT / "build-b1r2-release-cuda"
DEFAULT_BINARY_NAME = "ReciprocalAccP1bPhysicalSCF"
CPU_THREADS = (1, 2, 4, 8)
SI_SIZES = (1, 2, 3, 4)
CHEB_WARMUPS = 2
CHEB_REPETITIONS = 5
RECIP_WARMUPS = 1
RECIP_REPETITIONS = 3
CHEB_TIMING_NSTEP = 1
RECIP_TIMING_NSTEP = 1
CHEB_ORDER = 200
TIER_NAMES = ("lean", "full")
TIER_CONFIG = {
    "lean": {
        "si_sizes": (1, 2),
        "cpu_threads": {1: (1, 4, 8), 2: (8,)},
        "cheb_warmups": 0,
        "cheb_repetitions": 1,
        "recip_warmups": 0,
        "recip_repetitions": 1,
    },
    "full": {
        "si_sizes": SI_SIZES,
        "cpu_threads": {size: CPU_THREADS for size in SI_SIZES},
        "cheb_warmups": CHEB_WARMUPS,
        "cheb_repetitions": CHEB_REPETITIONS,
        "recip_warmups": RECIP_WARMUPS,
        "recip_repetitions": RECIP_REPETITIONS,
    },
}
DOS_TOLERANCES = {"relative_l2": 0.20, "max_abs": 12.5, "integrated": 5.0e-3}
OBSERVABLE_TOLERANCES = {
    "fermi_energy": 5.0e-4,
    "total_charge": 5.0e-4,
    "site_charge_transfer": 5.0e-4,
    "site_moment": 5.0e-4,
    "energy_per_atom": 5.0e-4,
    "density_residual": 5.0e-4,
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
    "case_id", "lane", "kind", "size", "supercell", "backend", "omp_threads", "blas_threads",
    "warmups", "repetitions", "Natom", "nmat", "nnz", "chebyshev_order", "chebyshev_kernel",
    "spectral_bounds_policy", "M", "P_rs_solver_kernel", "P_rs_phase", "P_scf_iteration_total",
    "steady_iteration_median", "T_rs_kernel", "T_rs_H2D", "T_rs_D2H", "T_rs_sync", "T_rs_setup",
    "numeric_mode", "hamiltonian_precision", "recurrence_precision", "moment_precision",
    "spectral_reconstruction_precision", "density_precision", "scf_canonical_precision",
    "status", "profile_status", "correctness_status", "common_state_id", "equal_precision_eligible",
    "production_comparison_eligible", "S_iteration", "S_cheb_kernel", "S_cheb_phase",
    "R_iteration_production", "R_cheb_kernel_production", "R_cheb_phase_production", "reason",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return default


def command_text(command: Iterable[str]) -> str:
    import shlex
    return shlex.join([str(item) for item in command])


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_hash(value: Any, length: int = 16) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:length]


def git_value(*args: str) -> str:
    result = subprocess.run(["git", *args], cwd=ROOT, check=False, capture_output=True, text=True)
    return result.stdout.strip()


def git_state() -> dict[str, Any]:
    tracked = subprocess.run(["git", "diff", "--quiet"], cwd=ROOT, check=False).returncode != 0
    staged = subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=ROOT, check=False).returncode != 0
    status = git_value("status", "--short")
    return {
        "commit": git_value("rev-parse", "HEAD"),
        "git_dirty_tracked": tracked or staged,
        "git_untracked_present": bool(git_value("status", "--porcelain", "--untracked-files=all")),
        "git_status": status,
    }


def ensure_output_tree(output: Path) -> None:
    for name in ("raw", "correctness", "iteration_history", "chebyshev", "reciprocal_common_state", "plots"):
        (output / name).mkdir(parents=True, exist_ok=True)


def parse_cache(build_dir: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    path = build_dir / "CMakeCache.txt"
    if not path.is_file():
        return result
    for line in path.read_text(errors="replace").splitlines():
        if ":" not in line or "=" not in line or line.startswith(("#", "//")):
            continue
        left, value = line.split("=", 1)
        result[left.split(":", 1)[0]] = value
    return result


def compile_commands(build_dir: Path) -> list[str]:
    path = build_dir / "compile_commands.json"
    data = read_json(path, [])
    if isinstance(data, list):
        values = []
        for item in data:
            if not isinstance(item, dict):
                continue
            command = item.get("command")
            if not command and isinstance(item.get("arguments"), list):
                command = command_text(item["arguments"])
            if command:
                values.append(str(command))
        if values:
            return values
    result = subprocess.run(["ninja", "-C", str(build_dir), "-t", "commands"], cwd=ROOT, check=False, capture_output=True, text=True)
    return [line.strip() for line in result.stdout.splitlines() if line.strip()] if result.returncode == 0 else []


def build_preflight(build_dir: Path) -> dict[str, Any]:
    commands = compile_commands(build_dir)
    representative: dict[str, list[str]] = {}
    failures: list[str] = []
    for label, names in REPRESENTATIVE_SOURCES.items():
        matches = []
        for command in commands:
            source_names = {match.group(1).lower() for match in SOURCE_RE.finditer(command)}
            if "-c" in command and any(name.lower() in source_names for name in names):
                matches.append(command)
        representative[label] = matches[:4]
        if not matches:
            failures.append(f"no actual compile command found for {label}")
        elif any((OPTIMIZATION_RE.findall(command)[-1] if OPTIMIZATION_RE.findall(command) else None) in {"-O0", "/Od"} for command in matches):
            failures.append(f"{label} has a trailing debug optimization flag")
    cache = parse_cache(build_dir)
    if cache.get("CMAKE_BUILD_TYPE", "").lower() != "release":
        failures.append("build type is not Release")
    return {
        "schema": "rslmto.scf-b1r2.build-preflight.v1",
        "timestamp": utc_now(),
        "build_dir": str(build_dir),
        "build_type": cache.get("CMAKE_BUILD_TYPE", "unknown"),
        "compiler": cache.get("CMAKE_Fortran_COMPILER"),
        "effective_optimization": "-O3" if commands and not failures else None,
        "representative_compile_commands": representative,
        "status": "PASS" if not failures else "FAIL",
        "reason": "; ".join(failures) if failures else None,
    }


def run_logged(command: list[str], evidence: Path, *, cwd: Path = ROOT, env: dict[str, str] | None = None) -> dict[str, Any]:
    evidence.mkdir(parents=True, exist_ok=True)
    (evidence / "command.txt").write_text(command_text(command) + "\n", encoding="utf-8")
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=False, check=False)
    stdout_bytes = completed.stdout or b""
    stderr_bytes = completed.stderr or b""
    stdout = stdout_bytes.decode("utf-8", errors="replace")
    stderr = stderr_bytes.decode("utf-8", errors="replace")
    (evidence / "stdout.bin").write_bytes(stdout_bytes)
    (evidence / "stderr.bin").write_bytes(stderr_bytes)
    (evidence / "stdout.log").write_text(stdout, encoding="utf-8")
    (evidence / "stderr.log").write_text(stderr, encoding="utf-8")
    (evidence / "combined.log").write_text(stdout + ("\n" if stdout and stderr else "") + stderr, encoding="utf-8")
    return {"command": command, "returncode": completed.returncode, "stdout": stdout, "stderr": stderr, "wall_seconds": time.perf_counter() - started}


def ensure_build(output: Path, build_dir: Path, *, dry_run: bool) -> tuple[Path, dict[str, Any]]:
    build_dir = build_dir.resolve()
    binary = build_dir / "bin" / DEFAULT_BINARY_NAME
    preflight = build_preflight(build_dir)
    if binary.is_file() and preflight["status"] == "PASS":
        return binary, preflight
    if dry_run:
        preflight["dry_run"] = True
        preflight["reason"] = preflight.get("reason") or "dry-run: build not executed"
        return binary, preflight
    cache = parse_cache(build_dir)
    configure = [
        "cmake", "-S", str(ROOT), "-B", str(build_dir), "-G", cache.get("CMAKE_GENERATOR", "Ninja"),
        "-DCMAKE_BUILD_TYPE=Release", "-DENABLE_OPENMP=ON", "-DENABLE_CUDA_PLUGIN=ON",
        "-DENABLE_CUDA_RECIPROCAL=ON", "-DENABLE_SPGLIB=ON", "-DRUN_UNIT_TESTS=ON",
        "-DCMAKE_Fortran_FLAGS_RELEASE=-O3",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
    ]
    if cache.get("CMAKE_CUDA_COMPILER"):
        configure.append(f"-DCMAKE_CUDA_COMPILER={cache['CMAKE_CUDA_COMPILER']}")
    if cache.get("BLA_VENDOR"):
        configure.append(f"-DBLA_VENDOR={cache['BLA_VENDOR']}")
    configure_result = run_logged(configure, output / "raw" / "build" / "configure")
    if configure_result["returncode"] != 0:
        raise RuntimeError("CMake configure failed; see results/benchmarks/scf_b1r2/raw/build/configure")
    build_command = ["cmake", "--build", str(build_dir), "--target", DEFAULT_BINARY_NAME, "--parallel"]
    build_result = run_logged(build_command, output / "raw" / "build" / "compile")
    if build_result["returncode"] != 0:
        raise RuntimeError("Release build failed; see results/benchmarks/scf_b1r2/raw/build/compile")
    preflight = build_preflight(build_dir)
    preflight.update({"configure_command": configure, "build_command": build_command})
    write_json(output / "build_preflight.json", preflight)
    if preflight["status"] != "PASS" or not binary.is_file():
        raise RuntimeError(preflight.get("reason") or "clean Release build preflight failed")
    return binary, preflight


def _replace_symbol(text: str, old: str, new: str) -> str:
    return re.sub(rf"(symbol\s*=\s*')\s*{re.escape(old)}(\s*')", rf"\g<1>{new}\g<2>", text, count=1, flags=re.IGNORECASE)


def _si_lattice(size: int) -> str:
    primitive_vectors = ((0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    # The expanded basis is periodic in the replicated cell.  Keeping the
    # primitive lattice vectors here would place all translated sites in one
    # primitive cell and make the neighbour map unphysically dense.
    vectors = tuple(tuple(size * value for value in vector) for vector in primitive_vectors)
    bases = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
    sites = []
    for i in range(size):
        for j in range(size):
            for k in range(size):
                shift = tuple(i * primitive_vectors[0][axis] + j * primitive_vectors[1][axis] + k * primitive_vectors[2][axis] for axis in range(3))
                sites.extend(tuple(shift[axis] + base[axis] for axis in range(3)) for base in bases)
    nsite = len(sites)
    lines = [
        "&lattice", f"  nbulk_bulk = {nsite}", f"  ntot = {nsite}", f"  nbas = {nsite}", f"  nrec = {nsite}", "",
    ]
    for col, vector in enumerate(vectors, 1):
        lines.append(f"  a(:, {col}) = {vector[0]:.10f}, {vector[1]:.10f}, {vector[2]:.10f}")
    lines.append("")
    for index, coord in enumerate(sites, 1):
        lines.append(f"  crd(:, {index}) = {coord[0]:.10f}, {coord[1]:.10f}, {coord[2]:.10f}")
    lines.append("")
    for field in ("izp", "no", "iu", "ib", "irec"):
        lines.extend(f"  {field}({index}) = {index}" for index in range(1, nsite + 1))
        lines.append("")
    lines.append("/")
    return "\n".join(lines) + "\n"


def stage_si_fixture(size: int, target: Path) -> dict[str, Any]:
    source = ROOT / "tests" / "scf" / "cases" / "bulk" / "diamondSi"
    shutil.copytree(source, target)
    if size != 1:
        keep = {"input.nml", "lattice.nml", "Si1.nml", "Si2.nml", "README.md"}
        for path in target.iterdir():
            if path.name not in keep:
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        input_path = target / "input.nml"
        input_text = input_path.read_text(encoding="utf-8")
        input_text = re.sub(r"(ntype\s*=\s*)\d+", rf"\g<1>{2 * size**3}", input_text, count=1)
        atoms = ["&atoms", "  database = './'"] + [f"  label({index}) = 'Si{index}'" for index in range(1, 2 * size**3 + 1)] + ["/"]
        input_text = re.sub(r"(?ms)^&atoms\b.*?^/\s*", "\n".join(atoms) + "\n", input_text, count=1, flags=re.IGNORECASE)
        input_path.write_text(input_text, encoding="utf-8")
        (target / "lattice.nml").write_text(_si_lattice(size), encoding="utf-8")
        si1 = (target / "Si1.nml").read_text(encoding="utf-8")
        si2 = (target / "Si2.nml").read_text(encoding="utf-8")
        for index in range(1, 2 * size**3 + 1):
            template = si1 if index % 2 else si2
            old = "Si1" if index % 2 else "Si2"
            (target / f"Si{index}.nml").write_text(_replace_symbol(template, old, f"Si{index}"), encoding="utf-8")
    digest = hashlib.sha256()
    for path in sorted(target.glob("*.nml")):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return {
        "material": "diamondSi", "source": str(source), "supercell": f"{size}x{size}x{size}",
        "fixture_id": "diamondSi_chebyshev_supercell", "fixture_revision": "b1r2-generated-periodic-si-v1",
        "nsp": 1, "soc": "off", "basis": "sp", "size": size,
        "common_state_id": f"si{size}-{digest.hexdigest()[:16]}", "source_hash": digest.hexdigest(),
    }


def stage_fixture(spec: dict[str, Any], target: Path) -> dict[str, Any]:
    if spec["lane"] == "chebyshev":
        return stage_si_fixture(int(spec["size"]), target)
    source = scf_b0c.FIXTURES[spec["material"]]["source"]
    shutil.copytree(source, target)
    fixture = dict(scf_b0c.FIXTURES[spec["material"]])
    fixture.update({"common_state_id": f"{spec['material']}-normal-initial-v1", "size": None})
    return fixture


def case_specs(tier: str = "full") -> list[dict[str, Any]]:
    if tier not in TIER_CONFIG:
        raise ValueError(f"unknown B1R2 tier {tier!r}; choose from {TIER_NAMES}")
    config = TIER_CONFIG[tier]
    specs: list[dict[str, Any]] = []
    for size in config["si_sizes"]:
        for omp in config["cpu_threads"][size]:
            specs.append({
                "case_id": f"cheb_si{size}_cpu_omp{omp}", "label": f"Chebyshev Si{size} CPU OMP{omp}",
                "lane": "chebyshev", "kind": "timing", "size": size, "material": "si", "backend": "lapack",
                "omp_threads": omp, "warmups": config["cheb_warmups"], "repetitions": config["cheb_repetitions"],
                "nstep": CHEB_TIMING_NSTEP, "benchmark_level": "scf_iteration", "rs_solver": "chebyshev",
            })
        specs.append({
            "case_id": f"cheb_si{size}_gpu", "label": f"Chebyshev Si{size} GPU", "lane": "chebyshev",
            "kind": "timing", "size": size, "material": "si", "backend": "cuda", "omp_threads": 1,
            "warmups": config["cheb_warmups"], "repetitions": config["cheb_repetitions"], "nstep": CHEB_TIMING_NSTEP,
            "benchmark_level": "scf_iteration", "rs_solver": "chebyshev",
        })
        for backend, label, binary_backend, omp in (("cpu", "CPU", "lapack", 1), ("gpu", "GPU", "cuda", 1)):
            specs.append({
                "case_id": f"cheb_si{size}_common_{backend}", "label": f"Chebyshev Si{size} {label} common-state",
                "lane": "chebyshev", "kind": "common_state", "size": size, "material": "si",
                "backend": binary_backend, "omp_threads": omp, "warmups": 0, "repetitions": 1, "nstep": 1,
                "benchmark_level": "rs_electronic_structure", "rs_solver": "chebyshev",
            })
    for length in (2, 3):
        material = f"fe{length}"
        for backend, label in (("lapack", "CPU"), ("cuda", "GPU")):
            specs.append({
                "case_id": f"recip_fe{length}_{'cpu' if backend == 'lapack' else 'gpu'}",
                "label": f"Reciprocal Fe{length} {label}", "lane": "reciprocal", "kind": "timing",
                "size": length, "material": material, "backend": backend, "omp_threads": 1,
                "warmups": config["recip_warmups"], "repetitions": config["recip_repetitions"], "nstep": RECIP_TIMING_NSTEP,
                "benchmark_level": "scf_iteration", "rs_solver": "block",
            })
            specs.append({
                "case_id": f"recip_fe{length}_common_{'cpu' if backend == 'lapack' else 'gpu'}",
                "label": f"Reciprocal Fe{length} {label} common-state", "lane": "reciprocal",
                "kind": "common_state", "size": length, "material": material, "backend": backend,
                "omp_threads": 1, "warmups": 0, "repetitions": 1, "nstep": 1,
                "benchmark_level": "reciprocal_phase", "rs_solver": "block",
            })
    for spec in specs:
        spec["tier"] = tier
    return specs


def _spec_metadata(spec: dict[str, Any], fixture: dict[str, Any]) -> dict[str, Any]:
    if spec["lane"] == "chebyshev":
        return {
            "material": "diamondSi", "fixture_id": fixture["fixture_id"], "fixture_revision": fixture["fixture_revision"],
            "supercell": fixture["supercell"], "nsp": 1, "soc": "off", "basis": "sp",
        }
    base = scf_b0c.FIXTURES[spec["material"]]
    return {"material": base["material"], "fixture_id": base["fixture_id"], "fixture_revision": "current-head-production-fixture-v1",
            "supercell": base["supercell"], "nsp": base["nsp"], "soc": base["soc"], "basis": base["basis"]}


def run_one(spec: dict[str, Any], binary: Path, case_dir: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    all_runs = [("warmup", index) for index in range(1, spec["warmups"] + 1)] + [("measured", index) for index in range(1, spec["repetitions"] + 1)]
    for phase, index in all_runs:
        run_dir = case_dir / f"{phase}-{index:02d}"
        fixture_dir = run_dir / "fixture"
        fixture = stage_fixture(spec, fixture_dir)
        env = scf_b0c._controlled_environment(spec["omp_threads"], 1)
        raw = scf_b0c.run_probe(
            binary=binary, fixture=fixture_dir, backend=spec["backend"], strategy="fp64_zheevd",
            scf_route="real_space" if spec["lane"] == "chebyshev" else "reciprocal",
            rs_solver=spec["rs_solver"], rs_backend="csr", dos_method="gaussian", sigma=0.01,
            temperature=300.0, nstep=spec["nstep"], benchmark_level=spec["benchmark_level"],
            omp_threads=spec["omp_threads"], blas_threads=1,
            raw_path=run_dir / "evidence" / "probe.log",
        )
        # run_probe controls its own environment; keep the explicit contract
        # visible in the case record for reproducibility.
        del env
        metadata = _spec_metadata(spec, fixture)
        row = scf_b0c._decorate_row(
            raw, metadata, name="si" if spec["lane"] == "chebyshev" else spec["material"], nstep=spec["nstep"],
            sigma=0.01, temperature=300.0, benchmark_level=spec["benchmark_level"],
            potential_identity=fixture["common_state_id"], scf_route="real_space" if spec["lane"] == "chebyshev" else "reciprocal",
            rs_backend="csr",
        )
        try:
            evidence_dir = str(run_dir.relative_to(ROOT))
        except ValueError:
            evidence_dir = str(run_dir)
        row.update({
            "case_id": spec["case_id"], "phase": phase, "repetition": index, "common_state_id": fixture["common_state_id"],
            "fixture_hash": fixture.get("source_hash"), "lane": spec["lane"], "kind": spec["kind"],
            "size": spec["size"], "warmups": spec["warmups"], "repetitions": spec["repetitions"],
            "evidence_dir": evidence_dir,
        })
        if spec["lane"] == "chebyshev":
            row.update({
                "M": row.get("chebyshev_order"), "recurrence_precision": "fp32" if spec["backend"] == "cuda" else "fp32-working/CPU route",
                "moment_precision": "fp32" if spec["backend"] == "cuda" else "fp64 storage with fast CPU working path",
                "spectral_reconstruction_precision": "fp64", "density_precision": "fp64", "scf_canonical_precision": "fp64",
            })
            row["rs_detail_timers_status"] = row.get("rs_detail_timers_status") or "not_exposed_by_backend"
        rows.append(row)
    measured = [row for row in rows if row["phase"] == "measured"]
    status_values = {row.get("status") for row in measured}
    if status_values == {"UNSUPPORTED"}:
        status = "UNSUPPORTED"
    elif any(value == "FAIL" for value in status_values):
        status = "FAIL"
    elif any(value == "UNSUPPORTED" for value in status_values):
        status = "UNSUPPORTED"
    else:
        status = "PASS"
    aggregate = aggregate_rows(measured, spec)
    aggregate.update({
        "case_id": spec["case_id"], "lane": spec["lane"], "kind": spec["kind"], "size": spec["size"],
        "status": status, "status_reason": None if status == "PASS" else "; ".join(sorted(str(row.get("status_reason") or row.get("unsupported_reason") or row.get("status")) for row in measured)),
        "common_state_id": measured[0].get("common_state_id") if measured else None,
        "repetition_rows": rows, "completed": True,
    })
    if spec["kind"] == "common_state" and measured:
        aggregate["dos"] = read_dos(case_dir / "measured-01" / "fixture") if spec["lane"] == "chebyshev" else None
    return {"schema": SCHEMA, "case_id": spec["case_id"], "spec": spec, "aggregate": aggregate, "completed": True}


def aggregate_rows(rows: list[dict[str, Any]], spec: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    numeric_fields = set()
    for row in rows:
        numeric_fields.update(key for key, value in row.items() if isinstance(value, (int, float)) and not isinstance(value, bool))
    for field in sorted(numeric_fields):
        values = [float(row[field]) for row in rows if isinstance(row.get(field), (int, float)) and math.isfinite(float(row[field]))]
        if values:
            result[field] = statistics.median(values)
            result[f"{field}_min"] = min(values)
            result[f"{field}_max"] = max(values)
            result[f"{field}_mad"] = statistics.median([abs(value - statistics.median(values)) for value in values])
    for field in ("backend", "solver_strategy", "profile_status", "numeric_mode", "rs_solver", "rs_backend", "chebyshev_kernel", "spectral_bounds_policy",
                  "chebyshev_order", "Natom", "nmat", "nnz", "hamiltonian_precision", "recurrence_precision", "moment_precision",
                  "spectral_reconstruction_precision", "density_precision", "scf_canonical_precision", "rs_detail_timers_status"):
        values = [row.get(field) for row in rows if row.get(field) is not None]
        if values:
            result[field] = values[0]
    # Runtime metadata uses canonical uppercase names; expose lowercase
    # aliases too so the CSV/Markdown presentation remains readable.
    if "OMP_threads" in result:
        result["omp_threads"] = result["OMP_threads"]
    if "BLAS_threads" in result:
        result["blas_threads"] = result["BLAS_threads"]
    result["steady_iteration_median"] = result.get("steady_iteration_median", result.get("P_scf_iteration_total"))
    result["final_state"] = rows[0].get("final_state", {}) if rows else {}
    if spec["lane"] == "chebyshev":
        phases = ("P_rs_hamiltonian_prepare", "P_rs_solver_kernel", "P_rs_green_function", "P_rs_spectral_reconstruct",
                  "P_rs_energy_integration", "P_rs_fermi", "P_rs_density_build", "P_rs_charge_spin_accumulate")
        result["P_rs_phase"] = sum(float(result.get(field, 0.0) or 0.0) for field in phases)
        result["P_rs_phase_min"] = None
        result["P_rs_phase_max"] = None
    result["profile_status"] = "PASS" if rows and all(row.get("profile_status") == "PASS" for row in rows) else "FAIL"
    return result


def read_dos(fixture: Path) -> dict[str, Any] | None:
    for name in ("totaldos.out", "density_of_states.dat", "density_of_states.out"):
        path = fixture / name
        if not path.is_file():
            continue
        values: list[tuple[float, float]] = []
        for line in path.read_text(errors="replace").splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                values.append((float(parts[0].replace("D", "E").replace("d", "e")), float(parts[1].replace("D", "E").replace("d", "e"))))
            except ValueError:
                continue
        if values:
            return {"path": str(path), "grid": [item[0] for item in values], "dos": [item[1] for item in values], "normalization": "production-file-first-two-columns"}
    return None


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _energy(state: dict[str, Any], natom: Any) -> float | None:
    for field in ("final_total_energy", "band_energy"):
        if _finite(state.get(field)):
            return float(state[field]) / max(float(natom or 1), 1.0)
    return None


def compare_common(cpu: dict[str, Any], gpu: dict[str, Any], *, require_dos: bool) -> dict[str, Any]:
    if cpu.get("status") != "PASS" or gpu.get("status") == "UNSUPPORTED":
        return {"status": "UNSUPPORTED" if gpu.get("status") == "UNSUPPORTED" else "INCONCLUSIVE", "reason": "common-state run unavailable"}
    if gpu.get("status") != "PASS":
        return {"status": "INCONCLUSIVE", "reason": "GPU common-state run failed"}
    left = cpu.get("final_state", {})
    right = gpu.get("final_state", {})
    differences: dict[str, float | None] = {}
    failures: list[str] = []
    for field in ("fermi_energy", "total_charge", "site_charge_transfer", "site_moment", "final_residual"):
        if _finite(left.get(field)) and _finite(right.get(field)):
            difference = abs(float(left[field]) - float(right[field]))
            differences[field] = difference
            key = "density_residual" if field == "final_residual" else field
            if difference > OBSERVABLE_TOLERANCES[key]:
                failures.append(field)
        else:
            differences[field] = None
            failures.append(field)
    energy_left = _energy(left, cpu.get("Natom"))
    energy_right = _energy(right, gpu.get("Natom"))
    if energy_left is None or energy_right is None:
        differences["energy_per_atom"] = None
        failures.append("energy_per_atom")
    else:
        differences["energy_per_atom"] = abs(energy_left - energy_right)
        if differences["energy_per_atom"] > OBSERVABLE_TOLERANCES["energy_per_atom"]:
            failures.append("energy_per_atom")
    dos_result: dict[str, Any] = {"status": "NOT_REQUESTED"}
    if require_dos:
        left_dos, right_dos = cpu.get("dos"), gpu.get("dos")
        if not left_dos or not right_dos:
            dos_result = {"status": "INCONCLUSIVE", "reason": "production DOS file missing"}
            failures.append("dos_missing")
        elif len(left_dos["grid"]) != len(right_dos["grid"]):
            dos_result = {"status": "FAIL", "reason": "DOS grid length mismatch"}
            failures.append("dos_grid")
        else:
            grid_delta = max(abs(a - b) for a, b in zip(left_dos["grid"], right_dos["grid"]))
            if grid_delta > 1.0e-10:
                dos_result = {"status": "FAIL", "reason": f"DOS energy grid mismatch ({grid_delta:g})"}
                failures.append("dos_grid")
            else:
                delta = [a - b for a, b in zip(left_dos["dos"], right_dos["dos"])]
                norm = math.sqrt(sum(value * value for value in left_dos["dos"]))
                relative_l2 = math.sqrt(sum(value * value for value in delta)) / max(norm, 1.0e-15)
                max_abs = max(abs(value) for value in delta)
                integrated_left = sum((x2 - x1) * (y1 + y2) * 0.5 for x1, x2, y1, y2 in zip(left_dos["grid"], left_dos["grid"][1:], left_dos["dos"], left_dos["dos"][1:]))
                integrated_right = sum((x2 - x1) * (y1 + y2) * 0.5 for x1, x2, y1, y2 in zip(right_dos["grid"], right_dos["grid"][1:], right_dos["dos"], right_dos["dos"][1:]))
                integrated = abs(integrated_left - integrated_right)
                dos_result = {"status": "PASS" if relative_l2 <= DOS_TOLERANCES["relative_l2"] and max_abs <= DOS_TOLERANCES["max_abs"] and integrated <= DOS_TOLERANCES["integrated"] else "FAIL",
                              "relative_l2": relative_l2, "max_abs": max_abs, "integrated_difference": integrated,
                              "grid_delta": grid_delta, "normalization": "same grid and production broadening"}
                if dos_result["status"] != "PASS":
                    failures.append("dos_tolerance")
    status = "PASS" if not failures else ("INCONCLUSIVE" if "dos_missing" in failures else "FAIL")
    return {"status": status, "reason": None if status == "PASS" else ", ".join(failures), "differences": differences, "dos": dos_result,
            "tolerances": {**OBSERVABLE_TOLERANCES, **DOS_TOLERANCES}}


def annotate_correctness(states: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id = {state["case_id"]: state for state in states}
    correctness: list[dict[str, Any]] = []
    for size in SI_SIZES:
        cpu = by_id.get(f"cheb_si{size}_common_cpu", {}).get("aggregate")
        gpu = by_id.get(f"cheb_si{size}_common_gpu", {}).get("aggregate")
        if cpu and gpu:
            result = compare_common(cpu, gpu, require_dos=True)
            record = {"lane": "chebyshev", "size": size, "common_state_id": cpu.get("common_state_id"), **result}
            correctness.append(record)
            for aggregate in (cpu, gpu):
                aggregate["correctness_status"] = result["status"]
                aggregate["correctness"] = result
    for length in (2, 3):
        cpu = by_id.get(f"recip_fe{length}_common_cpu", {}).get("aggregate")
        gpu = by_id.get(f"recip_fe{length}_common_gpu", {}).get("aggregate")
        if cpu and gpu:
            result = compare_common(cpu, gpu, require_dos=False)
            record = {"lane": "reciprocal", "size": length, "common_state_id": cpu.get("common_state_id"), **result}
            correctness.append(record)
            for aggregate in (cpu, gpu):
                aggregate["correctness_status"] = result["status"]
                aggregate["correctness"] = result
    return correctness


def annotate_ratios(rows: list[dict[str, Any]], correctness: list[dict[str, Any]]) -> None:
    corr = {(item["lane"], item["size"]): item for item in correctness}
    for size in SI_SIZES:
        cpu_rows = [row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("size") == size and row.get("backend") == "lapack" and row.get("status") == "PASS"]
        gpu = next((row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("size") == size and row.get("backend") == "cuda"), None)
        if not cpu_rows or not gpu or gpu.get("status") != "PASS":
            continue
        best_kernel = min(cpu_rows, key=lambda row: float(row.get("P_rs_solver_kernel") or math.inf))
        best_phase = min(cpu_rows, key=lambda row: float(row.get("P_rs_phase") or math.inf))
        best_iteration = min(cpu_rows, key=lambda row: float(row.get("steady_iteration_median") or math.inf))
        eligible = corr.get(("chebyshev", size), {}).get("status") == "PASS"
        gpu["best_cpu_kernel_case_id"] = best_kernel.get("case_id")
        gpu["best_cpu_phase_case_id"] = best_phase.get("case_id")
        gpu["best_cpu_iteration_case_id"] = best_iteration.get("case_id")
        gpu["equal_precision_eligible"] = False
        gpu["production_comparison_eligible"] = eligible
        if eligible:
            gpu["R_cheb_kernel_production"] = float(best_kernel["P_rs_solver_kernel"]) / float(gpu["P_rs_solver_kernel"])
            gpu["R_cheb_phase_production"] = float(best_phase["P_rs_phase"]) / float(gpu["P_rs_phase"])
            gpu["R_iteration_production"] = float(best_iteration["steady_iteration_median"]) / float(gpu["steady_iteration_median"])
    for length in (2, 3):
        cpu = next((row for row in rows if row.get("lane") == "reciprocal" and row.get("kind") == "timing" and row.get("size") == length and row.get("backend") == "lapack"), None)
        gpu = next((row for row in rows if row.get("lane") == "reciprocal" and row.get("kind") == "timing" and row.get("size") == length and row.get("backend") == "cuda"), None)
        eligible = corr.get(("reciprocal", length), {}).get("status") == "PASS"
        if cpu and gpu:
            gpu["equal_precision_eligible"] = eligible
            gpu["production_comparison_eligible"] = eligible
            if eligible and float(gpu.get("steady_iteration_median") or 0.0) > 0:
                gpu["S_iteration"] = float(cpu["steady_iteration_median"]) / float(gpu["steady_iteration_median"])
            gpu["solver_component_timer_status"] = "INVALID_FOR_K1"
            gpu["solver_component_timer_reason"] = "K1 CPU P_eigensolver is suppressed; use complete iteration timing and K2 for solver evidence"


def load_preserved_evidence() -> dict[str, Any]:
    source = ROOT / "results" / "benchmarks" / "scf_b1r" / "campaign.json"
    campaign = read_json(source, {})
    return {
        "source": str(source), "source_commit": campaign.get("provenance", {}).get("commit"),
        "large_reciprocal_evidence": campaign.get("large_reciprocal_evidence", []),
        "block_recursion_rows": [row for row in campaign.get("rows", []) if row.get("rs_solver") == "block"],
        "reuse_policy": "preserved; no rerun unless clean B1R2 source changes numerical/solver code",
    }


def row_reason(row: dict[str, Any]) -> str:
    return str(row.get("status_reason") or row.get("solver_component_timer_reason") or row.get("correctness", {}).get("reason") or "-")


def write_campaign_outputs(output: Path, states: list[dict[str, Any]], manifest: dict[str, Any], preflight: dict[str, Any], provenance: dict[str, Any]) -> dict[str, Any]:
    correctness = annotate_correctness(states)
    rows = [state["aggregate"] for state in states]
    for row in rows:
        if "OMP_threads" in row:
            row["omp_threads"] = row["OMP_threads"]
        if "BLAS_threads" in row:
            row["blas_threads"] = row["BLAS_threads"]
    annotate_ratios(rows, correctness)
    campaign = {
        "schema": SCHEMA, "created": manifest.get("created"), "updated": utc_now(), "tier": manifest.get("tier"), "manifest": manifest,
        "build_preflight": preflight, "provenance": provenance, "rows": rows, "correctness": correctness,
        "chebyshev_configuration": {
            "order_M": CHEB_ORDER, "kernel": "fast", "spectral_bounds_policy": "resolved_runtime",
            "energy_window": "fixture input.nml: -1.5..1.0", "energy_grid": "fixture input.nml: channels_ldos=2000",
            "smearing": {"method": "gaussian", "sigma": 0.01}, "temperature": 300.0,
        },
        "precision_policy": {
            "cpu": "FP64 canonical production route; fast CPU working path recorded by runtime metadata",
            "gpu": "FP32 working recurrence with FP64 canonical state; numeric_mode=mixed",
            "rs_detail_timers_status": "not_exposed_by_backend",
        },
        "k1_timer_audit": {
            "solver_component_timer_status": "INVALID_FOR_K1", "action": "suppressed",
            "reason": "Existing B1R CPU P_eigensolver remains approximately 1e-5 s for nmat=18,144,486 and is not a physical eigensolver component timer.",
            "authoritative_k1_metric": "complete reciprocal electronic-structure / SCF iteration timing",
            "authoritative_solver_metric": "preserved K2 frozen-potential FP64 evidence",
        },
        "preserved_b1r_evidence": load_preserved_evidence(),
        "policy": {
            "performance_metric": "steady production SCF iteration wall time",
            "cycle_count_role": "correctness and stability diagnostic only",
            "equal_precision_speedups_require_matching_numeric_mode": True,
            "rs_vs_kspace_formulation_comparison": "INCONCLUSIVE and deferred",
        },
    }
    write_json(output / "campaign.json", campaign)
    with (output / "campaign.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            values = {field: row.get(field) for field in CSV_FIELDS}
            values["reason"] = row_reason(row)
            writer.writerow(values)
    tier = manifest.get("tier", "unknown")
    lines = ["# SCF-B1R2 campaign evidence", "", f"Schema: `{SCHEMA}`", f"Tier: `{tier}`", "", "## Table C1 — CPU Chebyshev scaling", "", "| size | Natom | nmat | M | OMP | kernel s | phase s | iteration s | status |", "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for row in rows:
        if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("backend") == "lapack":
            omp = row.get("omp_threads", row.get("OMP_threads", "-"))
            if isinstance(omp, float) and omp.is_integer():
                omp = int(omp)
            lines.append(f"| Si{row.get('size')} | {row.get('Natom','-')} | {row.get('nmat','-')} | {row.get('chebyshev_order','-')} | {omp} | {row.get('P_rs_solver_kernel','-')} | {row.get('P_rs_phase','-')} | {row.get('steady_iteration_median','-')} | {row.get('status','-')} |")
    lines += ["", "## Table C2 — CPU/GPU Chebyshev production comparison", "", "| size | best CPU OMP (iteration) | CPU kernel | GPU kernel | kernel ratio | CPU phase | GPU phase | phase ratio | CPU iteration | GPU iteration | iteration ratio | correctness | mode |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|"]
    for size in SI_SIZES:
        cpus = [row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("size") == size and row.get("backend") == "lapack" and row.get("status") == "PASS"]
        gpu = next((row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("size") == size and row.get("backend") == "cuda"), None)
        if cpus and gpu:
            best = min(cpus, key=lambda row: float(row.get("steady_iteration_median") or math.inf))
            common = next((item for item in correctness if item.get("lane") == "chebyshev" and item.get("size") == size), {})
            omp = best.get("omp_threads", best.get("OMP_threads", "-"))
            if isinstance(omp, float) and omp.is_integer():
                omp = int(omp)
            lines.append(f"| Si{size} | OMP{omp} | {best.get('P_rs_solver_kernel','-')} | {gpu.get('P_rs_solver_kernel','-')} | {gpu.get('R_cheb_kernel_production','-')} | {best.get('P_rs_phase','-')} | {gpu.get('P_rs_phase','-')} | {gpu.get('R_cheb_phase_production','-')} | {best.get('steady_iteration_median','-')} | {gpu.get('steady_iteration_median','-')} | {gpu.get('R_iteration_production','-')} | {common.get('status','-')} | {gpu.get('numeric_mode','-')} |")
    lines += ["", "## Table C3 — Chebyshev common-state correctness", "", "| size | common state | Fermi diff | charge diff | energy/atom diff | DOS rel L2 | DOS max abs | integrated DOS diff | status | reason |", "|---|---|---:|---:|---:|---:|---:|---:|---|---|"]
    for item in correctness:
        if item.get("lane") != "chebyshev":
            continue
        differences = item.get("differences", {})
        dos = item.get("dos", {})
        lines.append(f"| Si{item.get('size')} | {item.get('common_state_id','-')} | {differences.get('fermi_energy','-')} | {differences.get('total_charge','-')} | {differences.get('energy_per_atom','-')} | {dos.get('relative_l2','-')} | {dos.get('max_abs','-')} | {dos.get('integrated_difference','-')} | {item.get('status','-')} | {item.get('reason') or '-'} |")
    lines += ["", "## Reciprocal Fe2/Fe3 common-state gate", "", "| size | status | iteration ratio | reason |", "|---|---|---:|---|"]
    for length in (2, 3):
        row = next((item for item in rows if item.get("lane") == "reciprocal" and item.get("kind") == "timing" and item.get("size") == length and item.get("backend") == "cuda"), {})
        item = next((item for item in correctness if item.get("lane") == "reciprocal" and item.get("size") == length), {})
        lines.append(f"| Fe{length} | {item.get('status','-')} | {row.get('S_iteration','-')} | {item.get('reason') or '-'} |")
    lines += ["", "K1 component solver timing is explicitly suppressed as `INVALID_FOR_K1`; K2 frozen-potential evidence remains authoritative.", "", "RS-vs-k-space formulation equivalence remains `INCONCLUSIVE` and deferred."]
    (output / "campaign.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for row in rows:
        row_id = str(row.get("case_id"))
        write_json(output / "correctness" / f"{row_id}.json", {"case_id": row_id, "status": row.get("correctness_status"), "correctness": row.get("correctness")})
        write_json(output / "iteration_history" / f"{row_id}.json", {"case_id": row_id, "repetition_rows": row.get("repetition_rows", [])})
    write_json(output / "correctness" / "chebyshev_common_state.json", [item for item in correctness if item.get("lane") == "chebyshev"])
    write_json(output / "reciprocal_common_state" / "fe2_fe3.json", [item for item in correctness if item.get("lane") == "reciprocal"])
    return campaign


def write_plots(output: Path, campaign: dict[str, Any]) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return
    rows = campaign.get("rows", [])
    plot_dir = output / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    cpu = [row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("backend") == "lapack" and row.get("status") == "PASS"]
    gpu = [row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "timing" and row.get("backend") == "cuda" and row.get("status") == "PASS"]
    best_kernel = {size: min((row for row in cpu if row.get("size") == size), key=lambda row: float(row.get("P_rs_solver_kernel") or math.inf), default=None) for size in SI_SIZES}
    best_iteration = {size: min((row for row in cpu if row.get("size") == size), key=lambda row: float(row.get("steady_iteration_median") or math.inf), default=None) for size in SI_SIZES}
    for field, name, ylabel in (("P_rs_solver_kernel", "C1_chebyshev_kernel_vs_size.png", "Chebyshev kernel wall time (s)"), ("steady_iteration_median", "C2_chebyshev_iteration_vs_size.png", "Chebyshev SCF iteration wall time (s)")):
        plt.figure()
        selected = best_kernel if field == "P_rs_solver_kernel" else best_iteration
        xs = [size for size, row in selected.items() if row and row.get(field) is not None]
        plt.plot(xs, [selected[size][field] for size in xs], "o-", label="CPU best threaded")
        gx = [row.get("size") for row in gpu if row.get(field) is not None]
        plt.plot(gx, [row[field] for row in gpu if row.get(field) is not None], "o-", label="GPU mixed")
        plt.xlabel("Si replication")
        plt.ylabel(ylabel)
        plt.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / name, dpi=140)
        plt.close()
    plt.figure()
    ratios = [row for row in gpu if row.get("R_iteration_production") is not None]
    plt.axhline(1.0, color="black", linewidth=0.8)
    plt.plot([row["size"] for row in ratios], [row["R_iteration_production"] for row in ratios], "o-", label="CPU FP64 / GPU mixed")
    plt.xlabel("Si replication")
    plt.ylabel("production ratio (parity = 1)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / "C3_chebyshev_production_ratio_vs_size.png", dpi=140)
    plt.close()
    representative = next((item for item in campaign.get("correctness", []) if item.get("lane") == "chebyshev" and item.get("status") == "PASS"), None)
    if representative:
        size = representative["size"]
        state_rows = [row for row in rows if row.get("lane") == "chebyshev" and row.get("kind") == "common_state" and row.get("size") == size]
        plt.figure()
        for row in state_rows:
            dos = row.get("dos")
            if dos:
                plt.plot(dos["grid"], dos["dos"], label=row.get("backend"))
        plt.xlabel("energy")
        plt.ylabel("DOS")
        plt.title(f"Si{size} common-state DOS")
        plt.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / "C4_chebyshev_dos_common_state.png", dpi=140)
        plt.close()


def planned_command(spec: dict[str, Any], binary: Path) -> list[str]:
    return [str(binary), "--input", "input.nml", "--backend", spec["backend"], "--solver-strategy", "fp64_zheevd",
            "--dos-method", "gaussian", "--scf-route", "real_space" if spec["lane"] == "chebyshev" else "reciprocal",
            "--rs-solver", spec["rs_solver"], "--rs-backend", "csr", "--sigma", "0.01", "--temperature", "300",
            "--nstep", str(spec["nstep"]), "--benchmark-level", spec["benchmark_level"], "--profile"]


def update_manifest(output: Path, manifest: dict[str, Any], states: list[dict[str, Any]]) -> None:
    by_id = {state.get("case_id"): state for state in states}
    for entry in manifest["cases"]:
        state = by_id.get(entry["case_id"])
        if state:
            aggregate = state.get("aggregate", {})
            entry.update({"completed": bool(state.get("completed")), "status": aggregate.get("status"), "reason": aggregate.get("status_reason")})
    manifest["updated"] = utc_now()
    write_json(output / "manifest.json", manifest)


def completed_states(output: Path, specs: list[dict[str, Any]], *, include_failed: bool = False) -> list[dict[str, Any]]:
    states = []
    for spec in specs:
        state = read_json(output / "raw" / spec["case_id"] / "case.json")
        if not isinstance(state, dict) or not state.get("completed"):
            continue
        if state.get("schema") != SCHEMA or state.get("spec") != spec:
            continue
        if not include_failed and state.get("aggregate", {}).get("status") == "FAIL":
            continue
        states.append(state)
    return states


def make_manifest(output: Path, specs: list[dict[str, Any]], binary: Path, tier: str) -> dict[str, Any]:
    current = read_json(output / "manifest.json", {})
    if (isinstance(current, dict) and current.get("schema") == MANIFEST_SCHEMA
            and current.get("tier") == tier):
        return current
    return {"schema": MANIFEST_SCHEMA, "tier": tier, "created": utc_now(), "updated": utc_now(), "repo_root": str(ROOT), "sentinel": SENTINEL,
            "cases": [{"case_id": spec["case_id"], "label": spec["label"], "lane": spec["lane"], "kind": spec["kind"], "planned_command": planned_command(spec, binary), "completed": False} for spec in specs]}


def run_campaign(args: argparse.Namespace) -> int:
    output = Path(args.output).resolve() if args.output else (DEFAULT_OUTPUT / args.tier).resolve()
    ensure_output_tree(output)
    specs = case_specs(args.tier)
    build_dir = Path(args.build_dir).resolve() if args.build_dir else DEFAULT_BUILD
    binary, preflight = ensure_build(output, build_dir, dry_run=args.dry_run)
    provenance = git_state()
    write_json(output / "build_preflight.json", preflight)
    write_json(output / "PRE" / "provenance.json", {"provenance": provenance, "build_preflight": preflight})
    manifest = make_manifest(output, specs, binary, args.tier)
    write_json(output / "manifest.json", manifest)
    if args.dry_run:
        print(f"B1R2 {args.tier} dry-run output: {output}")
        print(f"tracked source clean: {not provenance['git_dirty_tracked']}")
        print(f"build preflight: {preflight.get('status')} ({preflight.get('reason') or 'ready'})")
        for index, spec in enumerate(specs, 1):
            print(f"[B1R2 {index:02d}/{len(specs)}] {spec['label']}: {command_text(planned_command(spec, binary))}")
        return 0
    if provenance["git_dirty_tracked"]:
        raise RuntimeError("tracked source is dirty; commit the B1R2 source/build/harness changes before benchmarking")
    if preflight.get("status") != "PASS" or not binary.is_file():
        raise RuntimeError(f"clean Release build preflight failed: {preflight.get('reason') or binary}")
    if args.aggregate_only:
        states = completed_states(output, specs, include_failed=True)
        campaign = write_campaign_outputs(output, states, manifest, preflight, provenance)
        write_plots(output, campaign)
        print(f"B1R2 aggregate output: {output}")
        return 0
    for index, spec in enumerate(specs, 1):
        state_path = output / "raw" / spec["case_id"] / "case.json"
        previous = read_json(state_path)
        if (args.resume and not args.force and isinstance(previous, dict) and previous.get("schema") == SCHEMA
                and previous.get("spec") == spec and previous.get("completed")
                and previous.get("aggregate", {}).get("status") in {"PASS", "UNSUPPORTED", "INCONCLUSIVE"}):
            print(f"[B1R2 {index:02d}/{len(specs)}] {spec['label']} (resume: {previous['aggregate'].get('status')})", flush=True)
            continue
        print(f"[B1R2 {index:02d}/{len(specs)}] {spec['label']}", flush=True)
        try:
            state = run_one(spec, binary, output / "raw" / spec["case_id"])
        except Exception as exc:
            state = {"schema": SCHEMA, "case_id": spec["case_id"], "spec": spec, "completed": True,
                     "aggregate": {"case_id": spec["case_id"], "status": "FAIL", "status_reason": f"driver exception: {type(exc).__name__}: {exc}", "repetition_rows": []}}
        write_json(state_path, state)
        states = completed_states(output, specs, include_failed=True)
        update_manifest(output, manifest, states)
        campaign = write_campaign_outputs(output, states, manifest, preflight, provenance)
        write_plots(output, campaign)
    states = completed_states(output, specs)
    if len(states) != len(specs):
        missing = sorted(set(spec["case_id"] for spec in specs) - {state["case_id"] for state in states})
        print(f"B1R2_RUN_INCOMPLETE missing={','.join(missing)}")
        return 2
    campaign = write_campaign_outputs(output, states, manifest, preflight, provenance)
    write_plots(output, campaign)
    manifest["completed"] = True
    manifest["completed_at"] = utc_now()
    update_manifest(output, manifest, states)
    (output / SENTINEL).write_text(f"{utc_now()} {SENTINEL}\n", encoding="utf-8")
    print(f"B1R2 output: {output}")
    print(SENTINEL)
    return 0


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="run the complete B1R2 matrix")
    parser.add_argument("--tier", choices=TIER_NAMES, default="lean", help="lean desktop smoke or full HPC campaign (default: lean)")
    parser.add_argument("--output", default=None, help="result directory (default: results/benchmarks/scf_b1r2/<tier>)")
    parser.add_argument("--build-dir", default=str(DEFAULT_BUILD))
    parser.add_argument("--force", action="store_true", help="rerun completed cases")
    parser.add_argument("--resume", action="store_true", default=True, help="reuse completed PASS/UNSUPPORTED/INCONCLUSIVE cases")
    parser.add_argument("--dry-run", action="store_true", help="print the full plan without building or timing")
    parser.add_argument("--aggregate-only", action="store_true", help="rebuild canonical outputs from completed evidence")
    return parser


def main(argv: list[str] | None = None) -> int:
    return run_campaign(make_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
