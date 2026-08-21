#!/usr/bin/env python3
"""Canonical KPM-B0C CPU/GPU transport evidence pipeline.

The historical filename is retained so existing G1.2/G2 commands continue to
work. The campaign itself is now material-general, precision-explicit, and
strict about pairing. JSON is canonical; CSV and Markdown are derived from
the same row objects.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from benchmark_harness import (
    KPM_TOLERANCE_SETS,
    build_comparison_key,
    capture_environment,
    compare_production_outputs,
    comparison_key_fingerprint,
    empty_correctness,
    numeric_mode,
    stable_hash,
    validate_kpm_profile,
    validate_pairing,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from val09_kubo_bastin_transport import patch_case  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tests/run_binary.sh"


@dataclass(frozen=True)
class MaterialSpec:
    name: str
    base: Path
    fixture_id: str
    soc_state: str
    nsp: int
    default_cond_type: str
    va: tuple[int, int, int]
    vb: tuple[int, int, int]
    default_rc: int


MATERIALS = {
    "pt": MaterialSpec(
        "fccPt_SOC", ROOT / "tests/postproc/cases/conductivity/fccPt",
        "fccPt_SOC_val09", "SOC", 2, "spin", (0, 1, 0), (1, 0, 0), 20,
    ),
    "fe": MaterialSpec(
        "bccFe_magnetic", ROOT / "tests/regression/triad_bccFe_conductivity",
        "bccFe_magnetic_triad", "magnetic_nsp2", 2, "charge", (1, 0, 0), (1, 0, 0), 80,
    ),
}

PROFILE_STAGE_NAMES = (
    "P_operator_setup", "P_trace_setup", "P_moments_total", "P_gamma_basis_setup",
    "P_gamma_generation", "P_reconstruction_total", "P_result_unpack",
    "P_energy_integration", "P_tensor_postprocess", "P_output_prepare",
    "P_output_io", "P_stack_setup", "P_moment_finalize", "P_misc",
    "T_transport_total",
)

CSV_COLUMNS = (
    "row_id", "material", "size", "N", "nnz", "M", "NE", "lld", "cond_type",
    "cond_calctype", "Ntrace", "seed", "moment_backend", "moment_precision",
    "reconstruction_backend", "reconstruction_precision", "numeric_mode",
    "canonical_output_precision",
    "OMP_threads", "BLAS_threads", "GPU_strategy", "block_width", "profile_status",
    "correctness_status", "transport_median", "whole_wall_median", "moment_median",
    "Gamma_median", "reconstruction_median", "output_median", "S_moments",
    "S_transport", "S_whole", "headline_speedup_eligible", "best_cpu_row_id",
)


def summary(values: list[float]) -> dict[str, float]:
    if not values:
        return {"median": 0.0, "minimum": 0.0, "maximum": 0.0, "mad": 0.0, "iqr": 0.0}
    median = statistics.median(values)
    quartiles = statistics.quantiles(values, n=4, method="inclusive") if len(values) >= 2 else [median] * 3
    return {
        "median": median, "minimum": min(values), "maximum": max(values),
        "mad": statistics.median([abs(value - median) for value in values]),
        "iqr": quartiles[2] - quartiles[0],
    }


def environment(omp_threads: int, blas_threads: int) -> dict[str, str]:
    """Create one controlled, non-oversubscribed child environment."""

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    # run_binary.sh uses this explicit control for serial launches. Keep it
    # synchronized with the requested sweep instead of silently pinning all
    # benchmark rows to OMP=1.
    env["RSLMTO_OMP_THREADS_SERIAL"] = str(omp_threads)
    env["OMP_PROC_BIND"] = env.get("OMP_PROC_BIND", "close")
    env["OMP_PLACES"] = env.get("OMP_PLACES", "cores")
    env["BLAS_NUM_THREADS"] = str(blas_threads)
    env["MKL_NUM_THREADS"] = str(blas_threads)
    env["OPENBLAS_NUM_THREADS"] = str(blas_threads)
    return env


def fixture_revision(spec: MaterialSpec) -> str:
    digest = hashlib.sha256()
    for path in sorted(spec.base.glob("*.nml")):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()[:16]


def _copy_output_evidence(scratch: Path, destination: Path) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=True)
    checksums: dict[str, str] = {}
    for path in sorted(scratch.glob("*cond*.out")):
        target = destination / path.name
        shutil.copy2(path, target)
        checksums[path.name] = hashlib.sha256(target.read_bytes()).hexdigest()
    return checksums


def _profile_metric(profile: dict[str, Any], name: str) -> float:
    try:
        return float(profile.get("metrics", {}).get(name, 0.0))
    except (TypeError, ValueError):
        return 0.0


def _row_config(material: MaterialSpec, **values: Any) -> dict[str, Any]:
    return {
        "material": material.name, "fixture_id": material.fixture_id,
        "fixture_revision": fixture_revision(material), **values,
    }


def run_one(
    binary: Path, scratch: Path, *, material: MaterialSpec | str = MATERIALS["pt"],
    replication: int, cond_type: str, cheb_backend: str, gpu_plugin: bool,
    gpu_precision: str, omp_threads: int, blas_threads: int, cond_ll: int,
    lld: int, channels: int, rc: int, warmups: int, repetitions: int,
    cpu_reconstruction_precision: str = "fp64", raw_root: Path | None = None,
    no_write: bool = False, cond_calctype: str = "per_type", random_vec_num: int = 1,
    gpu_stochastic_block: int = 1, fermi: float | None = None,
    energy_min: float = -2.5, energy_max: float = 1.2,
) -> dict[str, Any]:
    if isinstance(material, str):
        material = MATERIALS[material]
    va = [1, 0, 0] if cond_type == "charge" else list(material.va)
    vb = list(material.vb)
    if fermi is None:
        fermi = -0.069099 if material.name.startswith("bccFe") else -0.085837
    config = _row_config(
        material, replication=replication, cond_type=cond_type, cheb_backend=cheb_backend,
        gpu_plugin=gpu_plugin, gpu_precision=gpu_precision,
        cpu_reconstruction_precision=cpu_reconstruction_precision, omp_threads=omp_threads,
        blas_threads=blas_threads, M=cond_ll, lld=lld, NE=channels + 10, rc=rc,
        cond_calctype=cond_calctype, Ntrace=random_vec_num,
        random_seed="271828" if cond_calctype == "random_vec" else "per_type_projectors",
        block=gpu_stochastic_block, output_mode="benchmark_no_write" if no_write else "production",
    )
    row_id = f"{material.name}_{stable_hash(config)}"
    row_raw = (raw_root / row_id) if raw_root is not None else scratch / "raw" / row_id
    row_raw.mkdir(parents=True, exist_ok=True)
    patch_case(
        material.base, scratch, cond_type=cond_type, va=va, vb=vb,
        replication=replication, cond_ll=cond_ll, lld=lld, channels=channels, rc=rc,
        fermi=fermi, energy_min=energy_min, energy_max=energy_max, gpu_plugin=gpu_plugin,
        gpu_backend="csr", cheb_backend=cheb_backend, gpu_precision=gpu_precision,
        cond_calctype=cond_calctype, random_vec_num=random_vec_num,
        gpu_stochastic_block=gpu_stochastic_block,
        cpu_reconstruction_precision=cpu_reconstruction_precision,
    )
    env = environment(omp_threads, blas_threads)
    if cond_calctype == "random_vec":
        env["RSLMTO_KPM_RANDOM_SEED"] = "271828"
    if no_write:
        env["RSLMTO_KPM_BENCHMARK_NO_WRITE"] = "1"

    raw_logs: dict[str, str] = {}
    for warmup in range(warmups):
        result = subprocess.run(["/bin/bash", str(RUNNER), str(binary)], cwd=scratch, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        log = (scratch / "testrun.log").read_text(errors="replace") if (scratch / "testrun.log").exists() else result.stdout
        log_path = row_raw / f"warmup_{warmup + 1:02d}.log"
        log_path.write_text(log, encoding="utf-8")
        raw_logs[f"warmup_{warmup + 1:02d}"] = str(log_path)
        if result.returncode:
            raise RuntimeError(f"KPM-B0C warmup failed for {row_id}:\n{log[-8000:]}")

    samples: list[dict[str, Any]] = []
    for repetition in range(repetitions):
        started = time.perf_counter()
        result = subprocess.run(["/bin/bash", str(RUNNER), str(binary)], cwd=scratch, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        wall_time = time.perf_counter() - started
        log = (scratch / "testrun.log").read_text(errors="replace") if (scratch / "testrun.log").exists() else result.stdout
        log_path = row_raw / f"sample_{repetition + 1:02d}.log"
        log_path.write_text(log, encoding="utf-8")
        raw_logs[f"sample_{repetition + 1:02d}"] = str(log_path)
        if result.returncode:
            raise RuntimeError(f"KPM-B0C measurement failed for {row_id}:\n{log[-8000:]}")
        from benchmark_harness import parse_profile_output
        profiles = [item for item in parse_profile_output(log) if item.get("name") == "kpm_transport"]
        if len(profiles) != 1:
            raise RuntimeError(f"{row_id}: expected one KPM profile, found {len(profiles)}")
        profile = profiles[0]
        validation = validate_kpm_profile(profile)
        if not validation["valid"]:
            raise RuntimeError(f"{row_id}: invalid KPM profile: {validation}")
        samples.append({"sample_id": f"{row_id}:sample:{repetition + 1:02d}",
                        "repetition": repetition + 1, "wall_time_s": wall_time,
                        "profile": profile, "raw_log": raw_logs[f"sample_{repetition + 1:02d}"]})

    first = samples[0]["profile"]
    metadata = first["metadata"]
    validation = first.get("validation") or validate_kpm_profile(first)
    moment_precision = str(metadata.get("moment_precision", "unknown"))
    reconstruction_precision = str(metadata.get("reconstruction_precision", "unknown"))
    trace_count = int(metadata.get("Ntrace") or (random_vec_num if cond_calctype == "random_vec" else 1))
    row: dict[str, Any] = {
        "row_id": row_id, "material": material.name, "size": f"r{replication}",
        "fixture_id": material.fixture_id, "fixture_revision": fixture_revision(material),
        "replication": replication, "N": metadata.get("N"), "nnz": metadata.get("nnz"),
        "M": cond_ll, "NE": channels + 10, "channels": channels, "lld": lld,
        "nsp": material.nsp, "soc_state": material.soc_state, "rc": rc,
        "energy_min": energy_min, "energy_max": energy_max, "fermi": fermi,
        "cond_type": cond_type, "cond_calctype": cond_calctype, "Ntrace": trace_count,
        "projector_count": trace_count,
        "random_seed": metadata.get("random_seed") or ("271828" if cond_calctype == "random_vec" else "per_type_projectors"),
        "vector_contract": "uniform_unit_phase_normalized_by_sqrt_sites" if cond_calctype == "random_vec" else "deterministic_type_projectors",
        "current_operator": {"va": va, "vb": vb, "spin_orbital_selector": "S_z" if cond_type == "spin" else ("L_z" if cond_type == "orbital" else None)},
        "kernel": {"name": "Lorentz", "alpha": 6.0},
        "chebyshev_scaling": {"contract": "(E-b)/a", "window": [energy_min, energy_max], "a_denominator": 1.7},
        "moment_backend": metadata.get("moment_backend"), "moment_precision": moment_precision,
        "reconstruction_backend": metadata.get("reconstruction_backend"), "reconstruction_precision": reconstruction_precision,
        "canonical_output_precision": "fp64",
        "numeric_mode": numeric_mode(moment_precision, reconstruction_precision),
        "backend": metadata.get("backend"), "gpu_plugin": gpu_plugin,
        "moment_strategy": "cuda_resident" if gpu_plugin else cheb_backend,
        "GPU_strategy": "resident_tiled" if gpu_plugin else None,
        "gpu_precision_request": gpu_precision if gpu_plugin else None,
        "block_width": int(metadata.get("trace_block_width") or gpu_stochastic_block),
        "gpu_stochastic_block": gpu_stochastic_block, "OMP_threads": omp_threads,
        "BLAS_threads": blas_threads, "OMP_NUM_THREADS": omp_threads, "BLAS_NUM_THREADS": blas_threads,
        "environment": {"OMP_PROC_BIND": env.get("OMP_PROC_BIND"), "OMP_PLACES": env.get("OMP_PLACES")},
        "output_mode": "benchmark_no_write" if no_write else "production",
        "profile_status": validation.get("status", "FAIL"), "profile_validation": validation,
        "warmups": warmups, "repetitions": repetitions, "profile_metadata": metadata,
        "samples": samples, "raw_logs": raw_logs,
        "correctness": empty_correctness("pending_pairing" if not no_write else "benchmark_no_write"),
    }
    key, fingerprint = comparison_key_fingerprint(row)
    row["comparison_key"] = key
    row["comparison_key_hash"] = fingerprint
    row["pairing_fingerprint"] = fingerprint
    output_dir = row_raw / "outputs"
    row["output_provenance"] = {"directory": str(output_dir),
                                "files": _copy_output_evidence(scratch, output_dir) if not no_write else {}}
    row["output_directory"] = str(output_dir)
    wall_values = [float(sample["wall_time_s"]) for sample in samples]
    metric_values = {name: [_profile_metric(sample["profile"], name) for sample in samples]
                     for name in PROFILE_STAGE_NAMES}
    row["statistics"] = {"wall_time_s": summary(wall_values)}
    row["statistics"].update({name: summary(values) for name, values in metric_values.items()})
    moment_per_trace = [value / max(trace_count, 1) for value in metric_values["P_moments_total"]]
    row["derived"] = {"moment_time_per_trace_s": summary(moment_per_trace),
                       "traces_per_second": summary([trace_count / max(value, 1.0e-30) for value in metric_values["P_moments_total"]])}
    row["gamma_bytes"] = first["metrics"].get("bytes_gamma")
    row["mu_packed_bytes"] = first["metrics"].get("bytes_mu_pack")
    return row


def memory_limited_row(*, material: MaterialSpec, replication: int, cond_type: str,
                       gpu_precision: str, omp_threads: int, cond_ll: int, lld: int,
                       random_vec_num: int, gpu_stochastic_block: int, error: str) -> dict[str, Any]:
    reason = next((line.strip() for line in reversed(error.splitlines()) if "stochastic workspace does not fit" in line),
                  "stochastic workspace does not fit")
    config = {"material": material.name, "replication": replication, "cond_type": cond_type,
              "M": cond_ll, "lld": lld, "random_vec_num": random_vec_num, "block": gpu_stochastic_block,
              "precision": gpu_precision}
    row_id = f"{material.name}_{stable_hash(config)}"
    row = {"row_id": row_id, "material": material.name, "size": f"r{replication}",
           "replication": replication, "cond_type": cond_type, "cond_calctype": "random_vec",
           "M": cond_ll, "lld": lld, "Ntrace": random_vec_num, "random_seed": "271828",
           "gpu_plugin": True, "gpu_precision_request": gpu_precision, "block_width": gpu_stochastic_block,
           "canonical_output_precision": "fp64",
           "OMP_threads": omp_threads, "status": "skipped_memory_limit", "reason": reason,
           "correctness": empty_correctness("memory_limited"), "headline_speedup_eligible": False,
           "ineligible_reasons": ["memory_limited"]}
    key, fingerprint = comparison_key_fingerprint(row)
    row["comparison_key"] = key
    row["comparison_key_hash"] = fingerprint
    return row


def attach_correctness(rows: list[dict[str, Any]], correctness_root: Path) -> None:
    """Attach production-output evidence outside timed samples."""

    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("status", "").startswith("skipped"):
            continue
        groups.setdefault(stable_hash({"physics": row.get("comparison_key", {}), "mode": row.get("numeric_mode")}), []).append(row)
    for candidates in groups.values():
        cpu = [row for row in candidates if not row.get("gpu_plugin")]
        for row in candidates:
            if row.get("output_mode") == "benchmark_no_write":
                row["correctness"] = empty_correctness("benchmark_no_write")
                continue
            if not row.get("output_provenance", {}).get("files"):
                row["correctness"] = empty_correctness("production_outputs_missing")
                continue
            reference = cpu[0] if cpu else row
            mode = row.get("numeric_mode", "mixed")
            evidence_id = stable_hash({"reference": reference["row_id"], "candidate": row["row_id"]})
            row["correctness"] = compare_production_outputs(
                Path(reference["output_directory"]), Path(row["output_directory"]), mode=mode,
                reference_row_id=reference["row_id"], evidence_dir=correctness_root / evidence_id)
            row["correctness"]["validation_evidence_id"] = evidence_id


def add_speedups(rows: list[dict[str, Any]]) -> None:
    """Compute speedups only after the strict pairing/correctness contract."""

    for row in rows:
        row.setdefault("headline_speedup_eligible", False)
        row.setdefault("ineligible_reasons", [])
        row.setdefault("speedups", {})
        row.setdefault("best_cpu_row_id", None)
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("status", "").startswith("skipped"):
            continue
        key = dict(row.get("comparison_key") or build_comparison_key(row))
        key.pop("numeric_mode", None)
        groups.setdefault(stable_hash(key), []).append(row)
    for candidates in groups.values():
        cpus = [row for row in candidates if not row.get("gpu_plugin")]
        for gpu in [row for row in candidates if row.get("gpu_plugin")]:
            possible = [cpu for cpu in cpus if cpu.get("numeric_mode") == gpu.get("numeric_mode")]
            valid = [cpu for cpu in possible if cpu.get("OMP_threads") in {1, 2, 4, 8} and validate_pairing(cpu, gpu)["eligible"]]
            if not valid:
                comparisons = [validate_pairing(cpu, gpu) for cpu in cpus]
                gpu["ineligible_reasons"] = list(dict.fromkeys(comparisons[0]["reasons"] if comparisons else ["no_cpu_reference"]))
                continue
            best = min(valid, key=lambda item: item["statistics"]["T_transport_total"]["median"])
            gpu["best_cpu_row_id"] = best["row_id"]
            gpu["best_cpu_reference"] = {"row_id": best["row_id"], "moment_backend": best.get("moment_backend"),
                                          "reconstruction_backend": best.get("reconstruction_backend"), "OMP_threads": best.get("OMP_threads"),
                                          "BLAS_threads": best.get("BLAS_threads"),
                                          "transport_median": best["statistics"]["T_transport_total"]["median"],
                                          "wall_median": best["statistics"]["wall_time_s"]["median"]}
            base, current = best["statistics"], gpu["statistics"]
            gpu["speedups"] = {
                "S_moments": base["P_moments_total"]["median"] / max(current["P_moments_total"]["median"], 1.0e-30),
                "S_transport": base["T_transport_total"]["median"] / max(current["T_transport_total"]["median"], 1.0e-30),
                "S_whole": base["wall_time_s"]["median"] / max(current["wall_time_s"]["median"], 1.0e-30),
            }
            gpu["headline_speedup_eligible"] = True
            gpu["ineligible_reasons"] = []
        for cpu in cpus:
            cpu.setdefault("ineligible_reasons", ["cpu_reference_row"])


def _csv_value(row: dict[str, Any], column: str) -> Any:
    mapping = {
        "OMP_threads": row.get("OMP_threads"), "BLAS_threads": row.get("BLAS_threads"), "GPU_strategy": row.get("GPU_strategy"),
        "block_width": row.get("block_width"), "profile_status": row.get("profile_status"),
        "correctness_status": row.get("correctness", {}).get("status"),
        "transport_median": row.get("statistics", {}).get("T_transport_total", {}).get("median"),
        "whole_wall_median": row.get("statistics", {}).get("wall_time_s", {}).get("median"),
        "moment_median": row.get("statistics", {}).get("P_moments_total", {}).get("median"),
        "Gamma_median": row.get("statistics", {}).get("P_gamma_basis_setup", {}).get("median", 0.0) +
                         row.get("statistics", {}).get("P_gamma_generation", {}).get("median", 0.0),
        "reconstruction_median": row.get("statistics", {}).get("P_reconstruction_total", {}).get("median"),
        "output_median": row.get("statistics", {}).get("P_output_io", {}).get("median"),
        "S_moments": row.get("speedups", {}).get("S_moments"), "S_transport": row.get("speedups", {}).get("S_transport"),
        "S_whole": row.get("speedups", {}).get("S_whole"), "headline_speedup_eligible": row.get("headline_speedup_eligible", False),
        "best_cpu_row_id": row.get("best_cpu_row_id"), "seed": row.get("random_seed"), "numeric_mode": row.get("numeric_mode"),
    }
    return mapping[column] if column in mapping else row.get(column)


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _csv_value(row, column) for column in CSV_COLUMNS})


def write_markdown(report: dict[str, Any], path: Path) -> None:
    rows = report.get("rows", [])
    env = report.get("environment", {})
    lines = ["# KPM-B0C closure campaign", "",
             "Strict fair CPU/GPU evidence package. JSON is canonical; this summary is derived from its rows.", "",
             "## Environment", "",
             f"- Git commit: `{env.get('git_commit')}`; dirty: `{env.get('git_dirty')}`",
             f"- Compiler: `{env.get('compiler')}`; build type: `{env.get('build_type')}`; BLAS/LAPACK: `{env.get('blas_lapack')}`",
             f"- CPU: `{env.get('cpu_model')}` ({env.get('physical_cpu_count')} physical / {env.get('logical_cpu_count')} logical); RAM: `{env.get('ram_mib')} MiB`",
             f"- CUDA: toolkit `{env.get('cuda_toolkit')}`, driver `{env.get('cuda_driver')}`, GPU `{env.get('gpu_model')}`, selected `{env.get('selected_gpu_index')}`, VRAM `{env.get('gpu_vram_mib')} MiB`, compute `{env.get('gpu_compute_capability')}`",
             "", "## Physics/workload and CPU thread winner by group", "",
             "| material | numeric mode | backend | OMP | transport median (s) | correctness | headline eligible |",
             "|---|---|---|---:|---:|---|---|"]
    for row in rows:
        if not row.get("gpu_plugin") and not row.get("status", "").startswith("skipped"):
            stat = row.get("statistics", {}).get("T_transport_total", {}).get("median", 0.0)
            lines.append(f"| {row.get('material')} | {row.get('numeric_mode', '-')} | {row.get('moment_backend', '-')} | {row.get('OMP_threads', '-')} | {stat:.5g} | {row.get('correctness', {}).get('status', 'NOT_APPLICABLE')} | {'yes' if row.get('headline_speedup_eligible') else 'no'} |")
    for title, mode in (("Equal-FP64 pairs", "fp64"), ("Equal-FP32 pairs", "fp32")):
        lines += ["", f"## {title}", "", "| GPU row | CPU reference | S_moments | S_transport | S_whole |", "|---|---|---:|---:|---:|"]
        for row in rows:
            if row.get("gpu_plugin") and row.get("numeric_mode") == mode:
                speedups = row.get("speedups", {})
                lines.append(f"| {row.get('row_id')} | {row.get('best_cpu_row_id')} | {speedups.get('S_moments', '-')} | {speedups.get('S_transport', '-')} | {speedups.get('S_whole', '-')} |")
    mixed = [row for row in rows if row.get("numeric_mode") == "mixed"]
    lines += ["", "## Mixed practical rows", "", f"{len(mixed)} mixed rows retained; they are never used for equal-precision headline speedups.", "", "## Failed/ineligible and memory-limited rows", ""]
    for row in rows:
        reasons = row.get("ineligible_reasons", [])
        if reasons or row.get("status", "").startswith("skipped"):
            lines.append(f"- `{row.get('row_id')}`: {', '.join(reasons or [row.get('status', 'ineligible')])}")
    counts: dict[str, int] = {}
    for row in rows:
        status = row.get("correctness", {}).get("status", "NOT_APPLICABLE")
        counts[status] = counts.get(status, 0) + 1
    lines += ["", "## Correctness summary", "", ", ".join(f"{key}: {value}" for key, value in sorted(counts.items())) or "no rows",
              "", "## Definitions", "", "- `S_moments = best_same_precision_CPU(P_moments_total) / GPU(P_moments_total)`.",
              "- `S_transport = best_same_precision_CPU(T_transport_total) / GPU(T_transport_total)`.",
              "- `S_whole = best_same_precision_CPU(whole_wall) / GPU(whole_wall)`.",
              "- `T_transport_total` is the internal transport phase; `whole_wall` is the complete process invocation."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(report: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(report["rows"], output.with_suffix(".csv"))
    write_markdown(report, output.with_suffix(".md"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--material", action="append", choices=tuple(MATERIALS), dest="material_one")
    parser.add_argument("--materials", nargs="+", choices=tuple(MATERIALS), default=None)
    parser.add_argument("--replications", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument("--cond-types", nargs="+", choices=("charge", "spin", "orbital"), default=None)
    parser.add_argument("--omp-threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--gpu-omp-threads", type=int, nargs="+", default=None,
                        help="OMP sweep for GPU rows; defaults to --omp-threads")
    parser.add_argument("--blas-threads", type=int, nargs="+", default=[1])
    parser.add_argument("--cond-ll", type=int, default=500)
    parser.add_argument("--lld", type=int, default=150)
    parser.add_argument("--channels", type=int, default=2500)
    parser.add_argument("--rc", type=int)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--cheb-backends", nargs="+", choices=("legacy", "fast", "fast_dp"), default=["legacy", "fast", "fast_dp"])
    parser.add_argument("--cpu-reconstruction-precisions", nargs="+", choices=("fp32", "fp64"), default=["fp64", "fp32"])
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-only", action="store_true", help="run only CUDA rows; requires --gpu")
    parser.add_argument("--gpu-precisions", nargs="+", choices=("fp32", "fp64"), default=["fp32", "fp64"])
    parser.add_argument("--cond-calctype", choices=("per_type", "random_vec"), default="per_type")
    parser.add_argument("--random-vec-num", type=int, default=16)
    parser.add_argument("--gpu-stochastic-block", type=int, nargs="+", default=[1])
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    if args.gpu_only and not args.gpu:
        parser.error("--gpu-only requires --gpu")
    materials = args.material_one or args.materials or ["pt"]
    gpu_omp_threads = args.gpu_omp_threads or args.omp_threads
    binary = args.binary.resolve()
    build_dir = args.build_dir.resolve() if args.build_dir else binary.parent.parent
    scratch_root = args.scratch_root.resolve()
    raw_root = args.output.resolve().parent / "raw"
    correctness_root = args.output.resolve().parent / "correctness"
    scratch_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    skipped_rows: list[dict[str, Any]] = []
    for material_name in materials:
        material = MATERIALS[material_name]
        for replication in args.replications:
            for cond_type in args.cond_types or [material.default_cond_type]:
                if not args.gpu_only:
                    for backend in args.cheb_backends:
                        for reconstruction_precision in args.cpu_reconstruction_precisions:
                            if backend != "fast" and reconstruction_precision != "fp64":
                                continue
                            for omp_threads in args.omp_threads:
                                for blas_threads in args.blas_threads:
                                    rows.append(run_one(
                                        binary, scratch_root / f"{material_name}_r{replication}_{cond_type}_{backend}_cpu{reconstruction_precision}_omp{omp_threads}_blas{blas_threads}",
                                        material=material, replication=replication, cond_type=cond_type,
                                        cheb_backend=backend, gpu_plugin=False, gpu_precision="fp64",
                                        cpu_reconstruction_precision=reconstruction_precision,
                                        omp_threads=omp_threads, blas_threads=blas_threads, cond_ll=args.cond_ll,
                                        lld=args.lld, channels=args.channels, rc=args.rc or material.default_rc,
                                        warmups=args.warmups, repetitions=args.repetitions, raw_root=raw_root,
                                        no_write=args.no_write, cond_calctype=args.cond_calctype,
                                        random_vec_num=args.random_vec_num))
                if args.gpu:
                    for gpu_precision in args.gpu_precisions:
                        widths = args.gpu_stochastic_block if args.cond_calctype == "random_vec" else [1]
                        for block_width in widths:
                            for omp_threads in gpu_omp_threads:
                                try:
                                    rows.append(run_one(
                                        binary, scratch_root / f"{material_name}_r{replication}_{cond_type}_cuda_{gpu_precision}_b{block_width}_omp{omp_threads}",
                                        material=material, replication=replication, cond_type=cond_type,
                                        cheb_backend="legacy", gpu_plugin=True, gpu_precision=gpu_precision,
                                        omp_threads=omp_threads, blas_threads=1, cond_ll=args.cond_ll, lld=args.lld,
                                        channels=args.channels, rc=args.rc or material.default_rc,
                                        warmups=args.warmups, repetitions=args.repetitions, raw_root=raw_root,
                                        no_write=args.no_write, cond_calctype=args.cond_calctype,
                                        random_vec_num=args.random_vec_num, gpu_stochastic_block=block_width))
                                except RuntimeError as exc:
                                    if args.cond_calctype == "random_vec" and "stochastic workspace does not fit" in str(exc):
                                        skipped_rows.append(memory_limited_row(
                                            material=material, replication=replication, cond_type=cond_type,
                                            gpu_precision=gpu_precision, omp_threads=omp_threads, cond_ll=args.cond_ll,
                                            lld=args.lld, random_vec_num=args.random_vec_num,
                                            gpu_stochastic_block=block_width, error=str(exc)))
                                        continue
                                    raise
    # Correctness is deliberately kept in this driver so the raw-output and
    # pair provenance remain in the same canonical campaign tree.
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(stable_hash({"key": row.get("comparison_key"), "mode": row.get("numeric_mode")}), []).append(row)
    for candidates in groups.values():
        cpus = [row for row in candidates if not row.get("gpu_plugin")]
        for row in candidates:
            if row.get("output_mode") == "benchmark_no_write" or not row.get("output_provenance", {}).get("files"):
                row["correctness"] = empty_correctness("benchmark_no_write" if args.no_write else "production_outputs_missing")
                continue
            reference = cpus[0] if cpus else row
            evidence_id = stable_hash({"reference": reference["row_id"], "candidate": row["row_id"]})
            row["correctness"] = compare_production_outputs(
                Path(reference["output_directory"]), Path(row["output_directory"]), mode=row.get("numeric_mode", "mixed"),
                reference_row_id=reference["row_id"], evidence_dir=correctness_root / evidence_id)
            row["correctness"]["validation_evidence_id"] = evidence_id
    add_speedups(rows)
    campaign_environment = capture_environment(ROOT, build_dir, mpi_ranks=1)
    campaign_environment.update({"omp_threads": args.omp_threads[0] if len(args.omp_threads) == 1 else None,
                                 "blas_threads": args.blas_threads[0] if len(args.blas_threads) == 1 else None,
                                 "omp_threads_sweep": args.omp_threads, "blas_threads_sweep": args.blas_threads,
                                 "omp_proc_bind": os.environ.get("OMP_PROC_BIND", "close"),
                                 "omp_places": os.environ.get("OMP_PLACES", "cores")})
    report = {
        "schema": "rslmto.kpm-b0c.v1",
        "scope": "strict fair CPU/GPU KPM/Kubo-Bastin transport benchmark evidence",
        "physics": {"materials": [MATERIALS[name].name for name in materials], "M": args.cond_ll,
                     "NE": args.channels + 10, "lld": args.lld, "estimator": args.cond_calctype,
                     "cond_calctype": args.cond_calctype, "random_vec_num": args.random_vec_num,
                     "gpu_stochastic_block": args.gpu_stochastic_block, "kernel": {"name": "Lorentz", "alpha": 6.0},
                     "paired_rows_share_input": True},
        "precision_contract": {"fp64": "moment_precision=fp64 and reconstruction_precision=fp64",
                               "fp32": "moment_precision=fp32 and reconstruction_precision=fp32",
                               "mixed": "all other supported combinations; canonical host output may still be FP64"},
        "pairing_contract": {"key_fields": list(build_comparison_key({}).keys()),
                             "performance_strategy_exclusions": ["OMP_threads", "BLAS_threads", "GPU_strategy", "block_width", "backend"],
                             "headline_requirements": ["identical comparison key", "same numeric mode", "correctness PASS", "profile PASS", "production output mode"]},
        "correctness_contract": {"source": "production *_cond.out files; no parallel Python Kubo-Bastin estimator",
                                 "tolerances": KPM_TOLERANCE_SETS,
                                 "moment_evidence": "optimized resident GPU path avoids full-moment D2H; diagnostic evidence is referenced when available"},
        "policy": {"warmups": args.warmups, "repetitions": args.repetitions, "persistent_process": False,
                    "omp_sweep": args.omp_threads, "blas_sweep": args.blas_threads, "cpu_rows_included": not args.gpu_only,
                    "gpu_only": args.gpu_only, "output_mode": "benchmark_no_write" if args.no_write else "production",
                    "speedup_definitions": {"S_moments": "best_same_precision_CPU(P_moments_total)/GPU(P_moments_total)",
                                             "S_transport": "best_same_precision_CPU(T_transport_total)/GPU(T_transport_total)",
                                             "S_whole": "best_same_precision_CPU(whole_wall)/GPU(whole_wall)"}},
        "environment": campaign_environment, "rows": rows, "skipped_rows": skipped_rows,
    }
    write_outputs(report, args.output.resolve())
    print(f"WROTE {args.output}: {len(rows)} rows, {len(skipped_rows)} skipped rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
