#!/usr/bin/env python3
"""Verify that optional Turek diagnostics do not alter H2 rotation physics."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


def fortran_float(token: str) -> float:
    value = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError(f"non-finite rotation output value: {token}")
    return value


def parse_rotation_output(path: Path) -> tuple[list[list[float]], list[float | None], str]:
    h2_rows: list[list[float]] = []
    turek_values: list[float | None] = []
    lines = path.read_text().splitlines()
    enabled_line = next((line for line in lines if line.startswith("# native_turek_diagnostic_enabled")), "")
    if not enabled_line:
        raise ValueError(f"missing Turek provenance header in {path}")

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < 17:
            raise ValueError(f"malformed rotation row in {path}: {line}")
        # q(3), q magnitude, finite-H, optional Turek, six pole/loss values,
        # circular channel, and four complex-kernel diagnostics plus resolution.
        h2_indices = list(range(5)) + list(range(6, 12)) + list(range(13, 17))
        h2_rows.append([fortran_float(fields[index]) for index in h2_indices])
        turek_values.append(None if fields[5] == "-" else fortran_float(fields[5]))
    if not h2_rows:
        raise ValueError(f"no rotation rows found in {path}")
    return h2_rows, turek_values, enabled_line


def run_variant(binary: Path, case_dir: Path, scratch_root: Path, enabled: bool) -> tuple[Path, str]:
    name = "native_on" if enabled else "native_off"
    scratch_dir = scratch_root / name
    if scratch_dir.exists():
        shutil.rmtree(scratch_dir)
    shutil.copytree(case_dir, scratch_dir)
    input_path = scratch_dir / "input.nml"
    input_text = input_path.read_text()
    replacement = ".true." if enabled else ".false."
    input_text, count = re.subn(
        r"(?im)^\s*native_turek\s*=.*$",
        f"native_turek = {replacement}",
        input_text,
        count=1,
    )
    if count != 1:
        raise ValueError(f"could not set native_turek in {input_path}")
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
        raise RuntimeError(
            f"native_turek={enabled} rotation run failed with exit code {completed.returncode};\n"
            f"{completed.stdout[-4000:]}"
        )
    output = scratch_dir / "rotation_dynamics.dat"
    if not output.is_file():
        raise RuntimeError(f"native_turek={enabled} did not produce {output}")
    return output, completed.stdout


def compare_rows(off_rows: list[list[float]], on_rows: list[list[float]]) -> None:
    if len(off_rows) != len(on_rows):
        raise AssertionError(f"H2 row count changed: off={len(off_rows)} on={len(on_rows)}")
    max_error = 0.0
    for row_index, (off_row, on_row) in enumerate(zip(off_rows, on_rows)):
        if len(off_row) != len(on_row):
            raise AssertionError(f"H2 column count changed at row {row_index}")
        for column_index, (off_value, on_value) in enumerate(zip(off_row, on_row)):
            error = abs(off_value - on_value)
            max_error = max(max_error, error)
            tolerance = 2.0e-11 + 2.0e-10 * max(abs(off_value), abs(on_value), 1.0)
            if error > tolerance:
                raise AssertionError(
                    f"H2 mismatch at row {row_index}, column {column_index}: "
                    f"off={off_value:.17g} on={on_value:.17g} error={error:.3e}"
                )
    print(f"H2 rotation rows identical within roundoff; max absolute difference = {max_error:.3e}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        args.scratch_root.mkdir(parents=True, exist_ok=True)
        off_path, off_stdout = run_variant(args.binary.resolve(), args.case.resolve(), args.scratch_root, False)
        on_path, on_stdout = run_variant(args.binary.resolve(), args.case.resolve(), args.scratch_root, True)
        off_rows, off_turek, off_header = parse_rotation_output(off_path)
        on_rows, on_turek, on_header = parse_rotation_output(on_path)
        if not off_header.endswith("= F"):
            raise AssertionError(f"disabled run did not report native Turek as disabled: {off_header}")
        if not on_header.endswith("= T"):
            raise AssertionError(f"enabled run did not report native Turek as enabled: {on_header}")
        if any(value is not None for value in off_turek):
            raise AssertionError("disabled run serialized Turek data instead of an explicit missing marker")
        if any(value is None for value in on_turek):
            raise AssertionError("enabled run did not serialize Turek diagnostic data")
        if "Fe validation campaign = NOT RUN" not in on_stdout:
            raise AssertionError("native_turek-only run did not report the Fe campaign as disabled")
        for campaign_marker in ("PASS-A", "PASS-B", "Fe material convergence", "Goldstone correction"):
            if campaign_marker in on_stdout:
                raise AssertionError(f"native_turek-only run emitted Fe campaign marker: {campaign_marker}")
        compare_rows(off_rows, on_rows)
        print("rotation Turek-independence regression: PASS")
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"rotation Turek-independence regression: FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
