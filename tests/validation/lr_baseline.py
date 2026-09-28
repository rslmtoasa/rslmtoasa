#!/usr/bin/env python3
"""Run one archived LR baseline case and compare its physics artifact."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


CASES = {
    "rotation": {"output": "rotation_dynamics.dat", "timeout": 300},
    "tddft_lehmann": {"output": "tddft_smoke.dat", "timeout": 180},
    "native_rsgf": {"output": "tddft_native_smoke.dat", "timeout": 300},
}


def parse_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values: list[float] = []
        for token in stripped.split():
            if token == "-":
                # Rotation output uses '-' for the circular-channel field.
                continue
            try:
                value = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                # The final status word in rotation output is intentionally
                # non-physics metadata and ends the row.
                break
            if not math.isfinite(value):
                raise ValueError(f"non-finite physics value in {path}: {token}")
            values.append(value)
        if values:
            rows.append(values)
    if not rows:
        raise ValueError(f"no numeric physics rows found in {path}")
    return rows


def compare_rows(actual: list[list[float]], expected: list[list[float]], abs_tol: float, rel_tol: float) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"row count mismatch: actual={len(actual)} expected={len(expected)}")
    for i, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        if len(actual_row) != len(expected_row):
            raise AssertionError(
                f"column count mismatch at row {i}: actual={len(actual_row)} expected={len(expected_row)}"
            )
        for j, (actual_value, expected_value) in enumerate(zip(actual_row, expected_row)):
            scale = max(abs(actual_value), abs(expected_value), 1.0)
            if abs(actual_value - expected_value) > abs_tol + rel_tol * scale:
                raise AssertionError(
                    f"physics mismatch at row {i}, column {j}: "
                    f"actual={actual_value:.17g} expected={expected_value:.17g}"
                )


def run_case(args: argparse.Namespace) -> None:
    spec = CASES[args.case]
    case_dir = Path(args.cases_root) / args.case
    reference_path = Path(args.references_dir) / f"{args.case}.json"
    scratch_root = Path(args.scratch_root)
    scratch_root.mkdir(parents=True, exist_ok=True)
    scratch_dir = scratch_root / args.case
    if scratch_dir.exists():
        shutil.rmtree(scratch_dir)
    shutil.copytree(case_dir, scratch_dir)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    try:
        completed = subprocess.run(
            [str(args.binary), "input.nml"],
            cwd=scratch_dir,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=args.timeout or spec["timeout"],
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "")[-4000:]
        raise RuntimeError(f"{args.case} timed out; tail of driver output:\n{output}") from exc
    if completed.returncode != 0:
        raise RuntimeError(
            f"{args.case} driver failed with exit code {completed.returncode}; "
            f"tail of driver output:\n{completed.stdout[-4000:]}"
        )

    output_path = scratch_dir / spec["output"]
    if not output_path.is_file():
        raise RuntimeError(f"{args.case} did not produce {spec['output']}")
    actual = parse_rows(output_path)
    reference = json.loads(reference_path.read_text())
    expected = reference["rows"]
    compare_rows(actual, expected, reference.get("abs_tol", 1.0e-9), reference.get("rel_tol", 1.0e-8))
    print(f"{args.case}: {len(actual)} physics row(s) matched archived reference")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", choices=sorted(CASES), required=True)
    parser.add_argument("--cases-root", type=Path, required=True)
    parser.add_argument("--references-dir", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--timeout", type=int)
    args = parser.parse_args()
    try:
        run_case(args)
    except (AssertionError, OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        print(f"lr baseline failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
