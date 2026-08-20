#!/usr/bin/env python3
"""KPM-G1.2 fair CPU/GPU transport benchmark campaign.

This driver keeps the Pt physics and production dimensions fixed while
sweeping backend/precision/thread choices.  It records every measured profile
row, rejects a KPM profile with failed exclusive-stage accounting, and writes
median/min/max/MAD/IQR summaries.  The default campaign is intentionally
explicit: it is suitable for a controlled host/GPU machine, not a quick CI
smoke test.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from benchmark_harness import capture_environment, parse_profile_output, validate_kpm_profile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from val09_kubo_bastin_transport import patch_case  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tests/run_binary.sh"
PT_BASE = ROOT / "tests/postproc/cases/conductivity/fccPt"


def summary(values: list[float]) -> dict[str, float]:
    median = statistics.median(values)
    quartiles = statistics.quantiles(values, n=4, method="inclusive") if len(values) >= 2 else [median] * 3
    return {
        "median": median,
        "minimum": min(values),
        "maximum": max(values),
        "mad": statistics.median([abs(value - median) for value in values]),
        "iqr": quartiles[2] - quartiles[0],
    }


def environment(omp_threads: int, blas_threads: int) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["OMP_PROC_BIND"] = env.get("OMP_PROC_BIND", "close")
    env["OMP_PLACES"] = env.get("OMP_PLACES", "cores")
    env["BLAS_NUM_THREADS"] = str(blas_threads)
    env["MKL_NUM_THREADS"] = str(blas_threads)
    env["OPENBLAS_NUM_THREADS"] = str(blas_threads)
    return env


def run_one(
    binary: Path,
    scratch: Path,
    *,
    replication: int,
    cond_type: str,
    cheb_backend: str,
    gpu_plugin: bool,
    gpu_precision: str,
    omp_threads: int,
    blas_threads: int,
    cond_ll: int,
    lld: int,
    channels: int,
    rc: int,
    warmups: int,
    repetitions: int,
) -> dict[str, Any]:
    patch_case(
        PT_BASE,
        scratch,
        cond_type=cond_type,
        va=[0, 1, 0] if cond_type != "charge" else [1, 0, 0],
        vb=[1, 0, 0],
        replication=replication,
        cond_ll=cond_ll,
        lld=lld,
        channels=channels,
        rc=rc,
        gpu_plugin=gpu_plugin,
        gpu_backend="csr",
        cheb_backend=cheb_backend,
        gpu_precision=gpu_precision,
    )
    env = environment(omp_threads, blas_threads)
    for _ in range(warmups):
        result = subprocess.run(
            ["/bin/bash", str(RUNNER), str(binary)], cwd=scratch, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
        )
        if result.returncode:
            raise RuntimeError(f"KPM-G1.2 warmup failed:\n{result.stdout[-8000:]}")

    samples: list[dict[str, Any]] = []
    for repetition in range(repetitions):
        started = time.perf_counter()
        result = subprocess.run(
            ["/bin/bash", str(RUNNER), str(binary)], cwd=scratch, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
        )
        wall_time = time.perf_counter() - started
        log_path = scratch / "testrun.log"
        log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
        if result.returncode:
            raise RuntimeError(f"KPM-G1.2 measurement failed (rep={repetition}):\n{log[-8000:]}")
        profiles = [item for item in parse_profile_output(log) if item.get("name") == "kpm_transport"]
        if len(profiles) != 1:
            raise RuntimeError(f"expected one KPM profile, found {len(profiles)}")
        profile = profiles[0]
        validation = validate_kpm_profile(profile)
        if not validation["valid"]:
            raise RuntimeError(f"invalid KPM profile: {validation}")
        samples.append({
            "repetition": repetition + 1,
            "wall_time_s": wall_time,
            "profile": profile,
        })

    first = samples[0]["profile"]
    metrics = first["metrics"]
    return {
        "size": f"r{replication}",
        "replication": replication,
        "N": first["metadata"].get("N"),
        "nnz": first["metadata"].get("nnz"),
        "M": cond_ll,
        "NE": channels + 10,
        "lld": lld,
        "cond_type": cond_type,
        "cond_calctype": "per_type",
        "moment_backend": first["metadata"].get("moment_backend"),
        "moment_precision": first["metadata"].get("moment_precision"),
        "reconstruction_backend": first["metadata"].get("reconstruction_backend"),
        "reconstruction_precision": first["metadata"].get("reconstruction_precision"),
        "gpu_plugin": gpu_plugin,
        "gpu_precision_request": gpu_precision if gpu_plugin else None,
        "OMP_NUM_THREADS": omp_threads,
        "BLAS_NUM_THREADS": blas_threads,
        "environment": {
            "OMP_PROC_BIND": env.get("OMP_PROC_BIND"),
            "OMP_PLACES": env.get("OMP_PLACES"),
        },
        "warmups": warmups,
        "repetitions": repetitions,
        "profile_metadata": first["metadata"],
        "samples": samples,
        "statistics": {
            "wall_time_s": summary([float(item["wall_time_s"]) for item in samples]),
            **{
                name: summary([float(item["profile"]["metrics"][name]) for item in samples])
                for name in (
                    "P_moments_total", "P_gamma", "P_reconstruction_total",
                    "P_energy_integration", "P_output_io", "P_other", "T_transport_total",
                )
            },
        },
        "gamma_bytes": metrics.get("bytes_gamma"),
        "mu_packed_bytes": metrics.get("bytes_mu_pack"),
        "correctness_status": "attach validated moment/conductivity comparison before speedup claim",
    }


def add_speedups(rows: list[dict[str, Any]]) -> None:
    """Attach equal-precision speedups against the best CPU row in each group."""

    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["size"], row["cond_type"], row["moment_precision"], row["reconstruction_precision"])
        groups.setdefault(key, []).append(row)
    for candidates in groups.values():
        cpu = [row for row in candidates if row["moment_backend"] in ("cpu_legacy", "cpu_fast")]
        if not cpu:
            continue
        best = min(cpu, key=lambda row: row["statistics"]["T_transport_total"]["median"])
        base_metrics = best["statistics"]
        for row in candidates:
            stats = row["statistics"]
            row["best_cpu_reference"] = {
                "moment_backend": best["moment_backend"],
                "OMP_NUM_THREADS": best["OMP_NUM_THREADS"],
                "BLAS_NUM_THREADS": best["BLAS_NUM_THREADS"],
            }
            row["speedups"] = {
                "S_moments": base_metrics["P_moments_total"]["median"] / stats["P_moments_total"]["median"]
                if stats["P_moments_total"]["median"] else math.nan,
                "S_transport": base_metrics["T_transport_total"]["median"] / stats["T_transport_total"]["median"]
                if stats["T_transport_total"]["median"] else math.nan,
                "S_whole": base_metrics["wall_time_s"]["median"] / stats["wall_time_s"]["median"]
                if stats["wall_time_s"]["median"] else math.nan,
            }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replications", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument("--cond-types", nargs="+", choices=("charge", "spin", "orbital"), default=["spin"])
    parser.add_argument("--omp-threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--blas-threads", type=int, nargs="+", default=[1])
    parser.add_argument("--cond-ll", type=int, default=500)
    parser.add_argument("--lld", type=int, default=150)
    parser.add_argument("--channels", type=int, default=2510 - 10)
    parser.add_argument("--rc", type=int, default=20)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--cheb-backends", nargs="+", choices=("legacy", "fast", "fast_dp"), default=["legacy", "fast", "fast_dp"])
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-precisions", nargs="+", choices=("fp32", "fp64"), default=["fp32", "fp64"])
    args = parser.parse_args()

    binary = args.binary.resolve()
    build_dir = args.build_dir.resolve() if args.build_dir else binary.parent.parent
    scratch_root = args.scratch_root.resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []

    for replication in args.replications:
        for cond_type in args.cond_types:
            for backend in args.cheb_backends:
                for omp_threads in args.omp_threads:
                    for blas_threads in args.blas_threads:
                        rows.append(run_one(
                            binary, scratch_root / f"r{replication}_{cond_type}_{backend}_omp{omp_threads}_blas{blas_threads}",
                            replication=replication, cond_type=cond_type, cheb_backend=backend,
                            gpu_plugin=False, gpu_precision="fp64", omp_threads=omp_threads,
                            blas_threads=blas_threads, cond_ll=args.cond_ll, lld=args.lld,
                            channels=args.channels, rc=args.rc, warmups=args.warmups,
                            repetitions=args.repetitions,
                        ))
            if args.gpu:
                for gpu_precision in args.gpu_precisions:
                    for omp_threads in args.omp_threads:
                        rows.append(run_one(
                            binary, scratch_root / f"r{replication}_{cond_type}_cuda_{gpu_precision}_omp{omp_threads}",
                            replication=replication, cond_type=cond_type, cheb_backend="legacy",
                            gpu_plugin=True, gpu_precision=gpu_precision, omp_threads=omp_threads,
                            blas_threads=1, cond_ll=args.cond_ll, lld=args.lld,
                            channels=args.channels, rc=args.rc, warmups=args.warmups,
                            repetitions=args.repetitions,
                        ))

    add_speedups(rows)
    report = {
        "schema": "rslmto.kpm-g1.2.v1",
        "scope": "precision-fair, exclusive-phase KPM/Kubo-Bastin transport benchmark",
        "physics": {
            "material": "real Pt SOC fixture",
            "M": args.cond_ll,
            "NE": args.channels + 10,
            "lld": args.lld,
            "estimator": "per_type",
            "kernel": "Lorentz(alpha=6)",
            "paired_rows_share_input": True,
        },
        "policy": {
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "persistent_process": False,
            "omp_sweep": args.omp_threads,
            "blas_sweep": args.blas_threads,
            "correctness": "performance rows require an attached precision-matched moment/conductivity comparison before publication",
        },
        "environment": capture_environment(ROOT, build_dir, mpi_ranks=1),
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"WROTE {args.output}: {len(rows)} valid rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
