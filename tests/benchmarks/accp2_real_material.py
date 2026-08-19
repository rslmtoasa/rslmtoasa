#!/usr/bin/env python3
"""ACC-P2 vector-first real-material campaign and timing-budget report.

This driver keeps the ACC-P0 persistent-process convention but narrows the
primary case to the reciprocal-SCF workload: eigenvectors requested, Fe
L=1,3,4,5, and five measured repetitions.  FP32 is not selected here unless
the preceding ACC-P1b evidence explicitly accepts it; the current policy is
FP64 Zheevd because the local P1b CUDA gate is unavailable on CPU-only hosts.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
P0 = ROOT / "tests/benchmarks/accp0_real_material.py"
MANDATORY_LENGTHS = (1, 3, 4, 5)


def load_rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))
    document = json.loads(path.read_text(encoding="utf-8"))
    return list(document.get("rows", []))


def number(row: dict[str, Any], key: str) -> float | None:
    value = row.get(key)
    if value in (None, "", "null"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def strategy_matches(row: dict[str, Any], strategy: str) -> bool:
    value = row.get("solver_strategy")
    aliases = {strategy}
    if strategy == "fp64_zheevd":
        aliases.add("zheevd_serial")
    return value in aliases


def matching(rows: list[dict[str, Any]], length: int, backend: str, strategy: str) -> dict[str, Any] | None:
    candidates = [
        row for row in rows
        if str(row.get("L")) == str(length)
        and row.get("backend") == backend
        and strategy_matches(row, strategy)
        and str(row.get("eigenvectors")) in {"1", "True", "true"}
    ]
    return candidates[0] if candidates else None


def budget_row(row: dict[str, Any]) -> dict[str, Any]:
    components = {
        "T_Hk_CPU": number(row, "T_Hk_CPU_s"),
        "T_host_staging": number(row, "T_host_staging_s"),
        "T_H2D": number(row, "H2D_s"),
        "T_solver": number(row, "solver_s"),
        "T_D2H_values": number(row, "T_D2H_values_s"),
        "T_D2H_vectors": number(row, "T_D2H_vectors_s"),
        "T_sync": number(row, "T_sync_s"),
        "T_other_backend": number(row, "T_other_backend_s"),
    }
    total = number(row, "T_total_s") or number(row, "total_steady_s")
    percentages = {
        key: None if value is None or not total else 100.0 * value / total
        for key, value in components.items()
    }
    return {
        "fixture": "bccFe",
        "nmat": int(row.get("nmat", 0)),
        "vectors": True,
        "precision": "FP64",
        "solver": row.get("solver_strategy"),
        "components_seconds": components,
        "percent_of_T_total": percentages,
        "T_total": total,
        "interval_repetitions": int(float(row.get("metric_repetitions", 0) or 0)),
        "H2D_bytes": int(float(row.get("H2D_bytes", 0) or 0)),
        "D2H_values_bytes": int(float(row.get("D2H_values_bytes", 0) or 0)),
        "D2H_vectors_bytes": int(float(row.get("D2H_vectors_bytes", 0) or 0)),
        "resource_counters": {
            key: int(float(row.get(key, 0) or 0))
            for key in (
                "cuda_malloc_count", "cuda_free_count", "workspace_query_count",
                "workspace_reuse_count", "event_create_count", "event_destroy_count",
                "pinned_alloc_count", "pinned_free_count",
            )
        },
        "resource_counters_before_measured_interval": {
            key: int(float(row.get(f"{key}_before", 0) or 0))
            for key in (
                "cuda_malloc_count", "cuda_free_count", "workspace_query_count",
                "workspace_reuse_count", "event_create_count", "event_destroy_count",
                "pinned_alloc_count", "pinned_free_count",
            )
        },
        "support_status": row.get("support_status"),
    }


def before_after(after_rows: list[dict[str, Any]], before_rows: list[dict[str, Any]] | None) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for length in MANDATORY_LENGTHS:
        cpu = matching(after_rows, length, "lapack", "lapack")
        gpu = matching(after_rows, length, "cuda", "fp64_zheevd")
        old = matching(before_rows or [], length, "cuda", "fp64_zheevd")
        after_total = number(gpu or {}, "total_steady_s")
        old_total = number(old or {}, "total_steady_s")
        cpu_total = number(cpu or {}, "total_steady_s")
        table.append({
            "fixture": "Fe primitive" if length == 1 else f"Fe {length}^{3}",
            "nmat": int((gpu or cpu or {}).get("nmat", 0)),
            "vectors": "yes",
            "precision": "FP64",
            "solver": "fp64_zheevd",
            "before_total": old_total,
            "after_total": after_total,
            "improvement": None if old_total is None or not old_total or not after_total else (old_total - after_total) / old_total,
            "CPU_total": cpu_total,
            "final_CPU_GPU_speedup": None if not after_total or not cpu_total else cpu_total / after_total,
            "after_support_status": (gpu or {}).get("support_status", "unsupported"),
        })
    return table


def pinned_comparison(pageable: list[dict[str, Any]], pinned: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for length in (3, 4, 5):
        page = matching(pageable, length, "cuda", "fp64_zheevd")
        pin = matching(pinned, length, "cuda", "fp64_zheevd")
        page_total = number(page or {}, "total_steady_s")
        pin_total = number(pin or {}, "total_steady_s")
        row = {
            "fixture": f"Fe {length}^{3}",
            "nmat": int((pin or page or {}).get("nmat", 0)),
            "pageable_total": page_total,
            "pinned_total": pin_total,
            "total_improvement": None if page_total is None or not page_total or not pin_total else (page_total - pin_total) / page_total,
        }
        for label, key in (("H2D", "H2D_s"), ("D2H_vectors", "T_D2H_vectors_s")):
            page_value = number(page or {}, key)
            pin_value = number(pin or {}, key)
            row[f"{label}_pageable"] = page_value
            row[f"{label}_pinned"] = pin_value
            row[f"{label}_improvement"] = None if page_value is None or not page_value or not pin_value else (page_value - pin_value) / page_value
        result.append(row)
    return result


def run_p0(args: argparse.Namespace, output_dir: Path, pinned: str, skip_cpu: bool = False) -> None:
    command = [
        sys.executable, str(P0),
        "--binary", str(args.binary), "--build-dir", str(args.build_dir),
        "--output-dir", str(output_dir), "--vectors", "--fe-lengths", "1,3,4,5",
        "--meshes", "1", "--tiles", "1", "--warmups", str(args.warmups),
        "--repetitions", str(args.repetitions), "--pinned-host", pinned,
        "--cuda-strategies", "fp64_zheevd", "--skip-cuda-validation",
    ]
    if args.skip_cuda:
        command.append("--skip-cuda")
    if skip_cpu:
        command.append("--skip-cpu")
    subprocess.run(command, cwd=ROOT, check=True, env={**os.environ, "RSLMTO_CUDA_PINNED_HOST": pinned})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--before-results", type=Path, default=None,
                        help="optional prior accp0_results.json/csv for the mandatory before/after table")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--skip-cuda", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="rslmto-accp2-", dir="/tmp") as scratch_name:
        scratch = Path(scratch_name)
        pageable_dir = scratch / "pageable"
        run_p0(args, pageable_dir, "0")
        pinned_dir = None
        if not args.skip_cuda:
            pinned_dir = scratch / "pinned"
            run_p0(args, pinned_dir, "1", skip_cpu=True)

        after_rows = load_rows(pageable_dir / "accp0_results.json")
        pinned_rows = load_rows(pinned_dir / "accp0_results.json") if pinned_dir else []
        before_rows = load_rows(args.before_results) if args.before_results else None
        primary = [budget_row(row) for row in after_rows
                   if row.get("backend") == "cuda" and row.get("solver_strategy") == "fp64_zheevd"
                   and str(row.get("eigenvectors")) == "1"
                   and str(row.get("L")) in {str(value) for value in MANDATORY_LENGTHS}]
        pinned = pinned_comparison(after_rows, pinned_rows)
        result = {
            "schema": "rslmto.accp2-real-material.v1",
            "policy": {
                "primary_workload": "reciprocal SCF, eigenvectors=yes",
                "repetitions": args.repetitions,
                "warmups": args.warmups,
                "production_precision": "FP64",
                "production_solver": "fp64_zheevd",
                "fp32_policy": "excluded pending a supported ACC-P1b scientific gate",
                "pinned_policy": "persistent backend-owned buffers, n>=486 only",
                "hk_gpu_port": False,
            },
            "accp1b_classification": {
                "Si_SCF": "unsupported",
                "metallic_Fe_SCF": "unsupported",
                "basis": "CUDA device unavailable in the completed local ACC-P1b run; no FP32 pass inferred",
            },
            "time_budgets": sorted(primary, key=lambda row: row["nmat"]),
            "pinned_comparison": sorted(pinned, key=lambda row: row["nmat"]),
            "before_after": before_after(after_rows, before_rows),
            "notes": [
                "T_sync is the host wait at the C++ return boundary; device-event transfer intervals are reported separately.",
                "T_D2H_vectors and D2H_vectors_bytes are the ACC-11 handoff evidence.",
                "A missing before-results file leaves before_total null rather than inventing a favorable baseline.",
            ],
        }
    output = args.output_dir / "accp2_results.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
