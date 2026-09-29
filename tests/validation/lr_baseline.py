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


def parse_named_rows(path: Path, schema: dict[str, object]) -> list[dict[str, object]]:
    """Parse a schema-described output whose rows include text fields."""
    columns = schema.get("columns")
    header = schema.get("header")
    text_columns = set(schema.get("text_columns", []))
    if not isinstance(columns, list) or not all(isinstance(name, str) for name in columns):
        raise ValueError(f"invalid named-column schema in {path}")
    if not isinstance(header, list) or not all(isinstance(name, str) for name in header):
        raise ValueError(f"invalid output-header schema in {path}")

    header_columns: list[str] | None = None
    rows: list[dict[str, object]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            candidate = stripped[1:].strip().split()
            if candidate and candidate[0] == header[0]:
                header_columns = candidate
            continue
        tokens = stripped.split()
        if len(tokens) != len(columns):
            raise ValueError(
                f"column count mismatch in {path}: actual={len(tokens)} expected={len(columns)}"
            )
        row: dict[str, object] = {}
        for name, token in zip(columns, tokens):
            if name in text_columns:
                row[name] = token
            elif token == "-":
                row[name] = None
            else:
                try:
                    value = float(token.replace("D", "E").replace("d", "e"))
                except ValueError as exc:
                    raise ValueError(f"non-numeric value for {name} in {path}: {token}") from exc
                if not math.isfinite(value):
                    raise ValueError(f"non-finite value for {name} in {path}: {token}")
                row[name] = value
        rows.append(row)
    if header_columns != header:
        raise ValueError(f"output schema header mismatch in {path}: {header_columns!r} != {header!r}")
    if not rows:
        raise ValueError(f"no named physics rows found in {path}")
    return rows


def compare_named_rows(
    actual: list[dict[str, object]],
    expected: list[dict[str, object]],
    schema: dict[str, object],
    abs_tol: float,
    rel_tol: float,
) -> None:
    strict_columns = schema.get("strict_columns", [])
    tolerant_columns = schema.get("tolerant_columns", {})
    diagnostic_columns = schema.get("diagnostic_columns", [])
    columns = schema["columns"]
    if len(actual) != len(expected):
        raise AssertionError(f"row count mismatch: actual={len(actual)} expected={len(expected)}")
    if not isinstance(strict_columns, list) or not isinstance(tolerant_columns, dict) or not isinstance(diagnostic_columns, list):
        raise AssertionError("invalid rotation comparison policy")
    covered = set(strict_columns) | set(tolerant_columns) | set(diagnostic_columns)
    if covered != set(columns) or len(covered) != len(strict_columns) + len(tolerant_columns) + len(diagnostic_columns):
        raise AssertionError("rotation comparison policy does not partition the named output schema")

    for i, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        for name in [*strict_columns, *tolerant_columns]:
            actual_value = actual_row.get(name)
            expected_value = expected_row.get(name)
            if isinstance(expected_value, str) or expected_value is None:
                if actual_value != expected_value:
                    raise AssertionError(
                        f"strict rotation mismatch at row {i}, column {name}: "
                        f"actual={actual_value!r} expected={expected_value!r}"
                    )
                continue
            if not isinstance(actual_value, (int, float)):
                raise AssertionError(
                    f"missing numeric rotation value at row {i}, column {name}: {actual_value!r}"
                )
            if name in strict_columns:
                column_abs_tol, column_rel_tol = abs_tol, rel_tol
            else:
                policy = tolerant_columns[name]
                column_abs_tol = float(policy.get("abs_tol", abs_tol))
                column_rel_tol = float(policy.get("rel_tol", rel_tol))
            scale = max(abs(actual_value), abs(float(expected_value)), 1.0)
            if abs(actual_value - float(expected_value)) > column_abs_tol + column_rel_tol * scale:
                raise AssertionError(
                    f"rotation mismatch at row {i}, column {name}: "
                    f"actual={actual_value:.17g} expected={float(expected_value):.17g}"
                )


def parse_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values: list[float] = []
        for token in stripped.split():
            if token == "-":
                # Rotation output uses '-' for the circular-channel field and
                # for optional diagnostics that were not requested.
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
    reference = json.loads(reference_path.read_text())
    abs_tol = reference.get("abs_tol", 1.0e-9)
    rel_tol = reference.get("rel_tol", 1.0e-8)
    if args.case == "rotation":
        schema = reference.get("schema")
        if not isinstance(schema, dict):
            raise ValueError("rotation reference is missing its named comparison schema")
        actual_named = parse_named_rows(output_path, schema)
        expected_named = reference["rows"]
        if not isinstance(expected_named, list) or not all(isinstance(row, dict) for row in expected_named):
            raise ValueError("rotation reference rows must be named objects")
        compare_named_rows(actual_named, expected_named, schema, abs_tol, rel_tol)
        print(f"{args.case}: {len(actual_named)} schema-aware physics row(s) matched archived reference")
    else:
        actual = parse_rows(output_path)
        expected = reference["rows"]
        compare_rows(actual, expected, abs_tol, rel_tol)
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
