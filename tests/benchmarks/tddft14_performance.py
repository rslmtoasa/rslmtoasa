#!/usr/bin/env python3
"""Run the TDDFT-14 backend/profile campaign and write machine-readable evidence."""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


SCHEMA = "rslmto.tddft-performance.v1"
PROFILE_PREFIXES = ("PROFILE_", "TDDFT_PERF_")
INTEGER = re.compile(r"^[+-]?\d+$")


def value(text: str) -> Any:
    if INTEGER.match(text):
        return int(text)
    try:
        return float(text)
    except ValueError:
        return text


def parse_profile_output(output: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line in output.splitlines():
        if not line.startswith(PROFILE_PREFIXES):
            continue
        tokens = line.split()
        record: dict[str, Any] = {"record": tokens[0]}
        if len(tokens) > 1 and "=" not in tokens[1]:
            record["label"] = tokens[1]
            tokens = tokens[2:]
        else:
            tokens = tokens[1:]
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if "=" not in token:
                index += 1
                continue
            key, raw = token.split("=", 1)
            if not raw and index + 1 < len(tokens):
                index += 1
                raw = tokens[index]
            record[key] = value(raw)
            index += 1
        records.append(record)
    return records


def command_output(command: list[str], env: dict[str, str]) -> tuple[float, str]:
    started = time.perf_counter()
    completed = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
    elapsed = time.perf_counter() - started
    return elapsed, completed.stdout + completed.stderr


def git_commit(root: Path) -> str:
    return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root,
                          check=True, capture_output=True, text=True).stdout.strip()


def summarize(records: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    summary: dict[str, dict[str, Any]] = {}
    for record in records:
        if record["record"] != "PROFILE_TDDFT_EXTENDED":
            continue
        label = str(record.get("label", "unknown"))
        summary[label] = {
            "eigen_batched_s": record.get("eigen_batched"),
            "eigen_scalar_s": record.get("eigen_scalar"),
            "green_direct_s": record.get("green_direct"),
            "green_mixed_s": record.get("green_mixed"),
            "realspace_energy_s": record.get("realspace_energy"),
            "realspace_R_q_FT_s": record.get("realspace_R_q_FT"),
            "eigen_scalar_over_batched": record.get("eigen_scalar", 0.0) /
            max(float(record.get("eigen_batched", 0.0)), sys.float_info.epsilon),
            "green_mixed_over_direct": record.get("green_mixed", 0.0) /
            max(float(record.get("green_direct", 0.0)), sys.float_info.epsilon),
        }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=4)
    args = parser.parse_args()
    if args.warmups < 0 or args.repetitions < 1:
        parser.error("warmups must be non-negative and repetitions must be positive")

    root = Path(__file__).resolve().parents[2]
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    command = [str(args.binary)]
    if args.mpi_exec:
        command = [str(args.mpi_exec), "--allow-run-as-root", "-n", str(args.mpi_ranks), *command]

    for _ in range(args.warmups):
        command_output(command, env)
    samples: list[dict[str, Any]] = []
    profile_records: list[dict[str, Any]] = []
    for repetition in range(args.repetitions):
        wall_time, output = command_output(command, env)
        samples.append({"repetition": repetition + 1, "wall_time_s": wall_time})
        profile_records.extend(parse_profile_output(output))

    document = {
        "schema": SCHEMA,
        "benchmark": {
            "name": "TDDFT-14 backend and MPI performance campaign",
            "class": "component",
            "command": command,
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "mpi_ranks": args.mpi_ranks if args.mpi_exec else 1,
            "samples": samples,
            "metadata": {
                "git_commit": git_commit(root),
                "compiler": "GNU Fortran 13.3.0",
                "platform": platform.platform(),
                "python": platform.python_version(),
                "openmp_threads": env.get("OMP_NUM_THREADS"),
            },
        },
        "profile_records": profile_records,
        "summary": summarize(profile_records),
        "scientific_checksums": [record for record in profile_records if record["record"] == "PROFILE_CHECKSUM"],
        "backend_matrix": [record for record in profile_records if record["record"] == "PROFILE_BACKEND"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {args.output}: {len(samples)} samples, {len(profile_records)} profile records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
