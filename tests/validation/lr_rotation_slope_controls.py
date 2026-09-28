#!/usr/bin/env python3
"""Check the canonical q=0 rotation finite-difference step."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


def set_namelist_value(text: str, key: str, value: str) -> str:
    pattern = rf"(?im)^\s*{re.escape(key)}\s*=.*$"
    updated, count = re.subn(pattern, f"{key} = {value}", text, count=1)
    if count == 1:
        return updated
    group_start = updated.lower().find("&linear_response")
    if group_start < 0:
        raise ValueError("rotation case has no &linear_response namelist")
    group_end = updated.find("/", group_start)
    if group_end < 0:
        raise ValueError("&linear_response namelist is unterminated")
    return updated[:group_end] + f"{key} = {value}\n" + updated[group_end:]


def parse_output(path: Path) -> tuple[dict[str, float], list[list[float]]]:
    headers: dict[str, float] = {}
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            match = re.match(r"#\s*([^=]+?)\s*=\s*([^\s]+)", stripped)
            if match:
                try:
                    headers[match.group(1).strip()] = float(match.group(2).replace("D", "E"))
                except ValueError:
                    pass
            continue
        values: list[float] = []
        for token in stripped.split():
            if token == "-":
                continue
            try:
                value = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                break
            if not math.isfinite(value):
                raise ValueError(f"non-finite rotation output value: {token}")
            values.append(value)
        if values:
            rows.append(values)
    if not rows:
        raise ValueError(f"no rotation rows found in {path}")
    required = {"Berry_commutator", "q0_slope_plus", "q0_slope_minus"}
    missing = required.difference(headers)
    if missing:
        raise ValueError(f"missing rotation diagnostics in {path}: {sorted(missing)}")
    return headers, rows


def run_variant(binary: Path, case_dir: Path, scratch_root: Path, name: str, slope_step: str | None) -> Path:
    scratch_dir = scratch_root / name
    if scratch_dir.exists():
        shutil.rmtree(scratch_dir)
    shutil.copytree(case_dir, scratch_dir)
    input_path = scratch_dir / "input.nml"
    input_text = input_path.read_text()
    input_text = set_namelist_value(input_text, "native_turek", ".false.")
    if slope_step is not None:
        input_text = set_namelist_value(input_text, "rotation_probe_omega", "1.0e-5")
        input_text = set_namelist_value(input_text, "rotation_slope_step", slope_step)
    input_path.write_text(input_text)
    completed = subprocess.run(
        [str(binary), "input.nml"],
        cwd=scratch_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{name} failed with exit code {completed.returncode}:\n{completed.stdout[-4000:]}")
    output = scratch_dir / "rotation_dynamics.dat"
    if not output.is_file():
        raise RuntimeError(f"{name} did not produce {output}")
    return output


def compare_values(label: str, left: float, right: float) -> None:
    tolerance = 2.0e-11 + 2.0e-10 * max(abs(left), abs(right), 1.0)
    if abs(left - right) > tolerance:
        raise AssertionError(f"{label} changed: {left:.17g} vs {right:.17g}")


def compare_rows(label: str, left: list[list[float]], right: list[list[float]]) -> None:
    if len(left) != len(right):
        raise AssertionError(f"{label} row count changed: {len(left)} vs {len(right)}")
    for row_index, (left_row, right_row) in enumerate(zip(left, right)):
        if len(left_row) != len(right_row):
            raise AssertionError(f"{label} column count changed at row {row_index}")
        for column_index, (left_value, right_value) in enumerate(zip(left_row, right_row)):
            compare_values(f"{label} row {row_index} column {column_index}", left_value, right_value)


def check_derivative(headers: dict[str, float]) -> None:
    berry = headers["Berry_commutator"]
    plus = headers["q0_slope_plus"]
    minus = headers["q0_slope_minus"]
    scale = max(abs(berry), 1.0e-12)
    if abs(abs(plus) - abs(berry)) > 5.0e-2 * scale or abs(abs(minus) - abs(berry)) > 5.0e-2 * scale:
        raise AssertionError(f"q=0 slope is not Berry-normalized: plus={plus} minus={minus} Berry={berry}")
    if abs(plus + minus) > 5.0e-2 * scale:
        raise AssertionError(f"q=0 circular slopes are not antisymmetric: plus={plus} minus={minus}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        args.scratch_root.mkdir(parents=True, exist_ok=True)
        default_path = run_variant(args.binary.resolve(), args.case.resolve(), args.scratch_root, "default", None)
        equal_path = run_variant(args.binary.resolve(), args.case.resolve(), args.scratch_root, "explicit_equal", "1.0e-5")
        unequal_path = run_variant(args.binary.resolve(), args.case.resolve(), args.scratch_root, "unequal_legacy", "2.0e-5")
        default_headers, default_rows = parse_output(default_path)
        equal_headers, equal_rows = parse_output(equal_path)
        unequal_headers, unequal_rows = parse_output(unequal_path)
        check_derivative(default_headers)
        check_derivative(equal_headers)
        check_derivative(unequal_headers)
        for key in ("Berry_commutator", "q0_slope_plus", "q0_slope_minus"):
            compare_values(f"explicit-equal {key}", default_headers[key], equal_headers[key])
            compare_values(f"unequal-legacy {key}", default_headers[key], unequal_headers[key])
        compare_rows("explicit-equal output", default_rows, equal_rows)
        compare_rows("unequal-legacy output", default_rows, unequal_rows)
        print("rotation slope-control regression: PASS")
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"rotation slope-control regression: FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
