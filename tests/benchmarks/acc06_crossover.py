#!/usr/bin/env python3
"""Run and summarize the ACC-06 reciprocal CPU/CUDA crossover campaign."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any


FIXTURES = (
    ("Si_sp", 1, 1),
    ("bccFe_spd", 1, 2),
    ("two_site_spd", 2, 2),
    ("multisite_4_lmax_2", 4, 2),
)


def median_wall(document: dict[str, Any]) -> float:
    return statistics.median(sample["wall_time_s"] for sample in document["benchmark"]["samples"])


def acc06_metadata(document: dict[str, Any]) -> dict[str, Any]:
    records = document.get("profile_records", [])
    if not records:
        raise ValueError(f"{document['benchmark']['name']}: missing ACC-06 profile record")
    return records[0]["metadata"]


def record_key(document: dict[str, Any]) -> tuple[Any, ...]:
    metadata = acc06_metadata(document)
    return (metadata["fixture"], metadata["matrix_dimension"], metadata["k_points"],
            metadata["tile_size"], metadata["eigenvectors"])


def run_one(args: argparse.Namespace, name: str, binary: Path, build_dir: Path,
            backend: str, strategy: str, sites: int, lmax: int, nk: int,
            tile_size: int, eigenvectors: int, output: Path) -> None:
    command = [str(binary), "--backend", backend, "--strategy", strategy,
               "--sites", str(sites), "--lmax", str(lmax), "--nk", str(nk),
               "--tile-size", str(tile_size), "--eigenvectors", str(eigenvectors)]
    harness = Path(__file__).with_name("benchmark_harness.py")
    subprocess.run([
        sys.executable, str(harness), "run", "--name", name, "--class", "component",
        "--labels", "performance", "component", "reciprocal", "eigensolver", "acc06",
        "--build-dir", str(build_dir), "--warmups", str(args.warmups),
        "--repetitions", str(args.repetitions), "--output", str(output),
        "--command", *command,
    ], check=True)


def summarize(documents: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for document in documents:
        grouped.setdefault(record_key(document), []).append(document)
    rows: list[dict[str, Any]] = []
    for key, entries in sorted(grouped.items(), key=lambda item: tuple(str(value) for value in item[0])):
        fixture, matrix_dimension, nk, tile_size, eigenvectors = key
        cpu = [entry for entry in entries if acc06_metadata(entry).get("backend") == "lapack"]
        gpu = [entry for entry in entries if acc06_metadata(entry).get("backend") == "cuda"]
        if not cpu:
            continue
        best_cpu = min(cpu, key=median_wall)
        cpu_time = median_wall(best_cpu)
        gpu_time = median_wall(gpu[0]) if gpu else None
        rows.append({
            "fixture": fixture, "matrix_dimension": matrix_dimension, "nk": nk,
            "tile_size": tile_size, "eigenvectors": bool(eigenvectors),
            "best_cpu": acc06_metadata(best_cpu).get("strategy", "backend"),
            "cpu_s": cpu_time, "gpu_s": gpu_time,
            "speedup": cpu_time / gpu_time if gpu_time else None,
            "recommended_backend": "cuda" if gpu_time and gpu_time < cpu_time else "lapack",
        })
    return rows


def write_report(path: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "# ACC-06 reciprocal CPU/GPU crossover campaign", "",
        "Wall time is the ACC-00 harness median, including executable setup and backend initialization.",
        "`cpu_s` is the best measured CPU strategy; GPU rows are absent without a CUDA binary.", "",
        "| fixture | matrix | Nk | tile | eigenvectors | best CPU | CPU s | GPU s | speedup | recommended |",
        "|---|---:|---:|---:|:---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        gpu = "n/a" if row["gpu_s"] is None else f"{row['gpu_s']:.6f}"
        speedup = "n/a" if row["speedup"] is None else f"{row['speedup']:.2f}x"
        lines.append(
            f"| {row['fixture']} | {row['matrix_dimension']} | {row['nk']} | {row['tile_size']} | "
            f"{'yes' if row['eigenvectors'] else 'no'} | {row['best_cpu']} | {row['cpu_s']:.6f} | "
            f"{gpu} | {speedup} | {row['recommended_backend']} |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cpu-build-dir", type=Path, required=True)
    parser.add_argument("--cuda-binary", type=Path)
    parser.add_argument("--cuda-build-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--quick", action="store_true", help="only representative sizes and Nk=8")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    nk_values = (8,) if args.quick else (1, 8, 32)
    tile_values = (8,) if args.quick else (1, 8, 16)
    documents: list[dict[str, Any]] = []
    runs = []
    for fixture, sites, lmax in FIXTURES:
        for nk in nk_values:
            for tile_size in tile_values:
                for eigenvectors in (0, 1):
                    runs.append((f"acc06_cpu_{fixture}_{nk}_{tile_size}_{eigenvectors}", args.cpu_binary,
                                 args.cpu_build_dir, "lapack", "backend", sites, lmax, nk, tile_size, eigenvectors))
                    runs.append((f"acc06_cpu_parallel_{fixture}_{nk}_{tile_size}_{eigenvectors}", args.cpu_binary,
                                 args.cpu_build_dir, "lapack", "parallel", sites, lmax, nk, tile_size, eigenvectors))
                    if args.cuda_binary and args.cuda_build_dir:
                        runs.append((f"acc06_cuda_{fixture}_{nk}_{tile_size}_{eigenvectors}", args.cuda_binary,
                                     args.cuda_build_dir, "cuda", "backend", sites, lmax, nk, tile_size, eigenvectors))
    for name, binary, build_dir, backend, strategy, sites, lmax, nk, tile_size, eigenvectors in runs:
        output = args.output_dir / f"{name}.json"
        run_one(args, name, binary, build_dir, backend, strategy, sites, lmax, nk, tile_size, eigenvectors, output)
        documents.append(json.loads(output.read_text(encoding="utf-8")))
    rows = summarize(documents)
    write_report(args.report, rows)
    (args.output_dir / "acc06_summary.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {args.report}: {len(rows)} rows from {len(documents)} runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
