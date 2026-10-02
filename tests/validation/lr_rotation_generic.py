#!/usr/bin/env python3
"""Check that the generic rotation workflow does not require Fe campaign gates."""

from __future__ import annotations

import argparse
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    scratch = args.scratch_root / "temperature_350_generic"
    try:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(args.case.resolve(), scratch)
        input_path = scratch / "input.nml"
        input_text = input_path.read_text()
        input_text = set_namelist_value(input_text, "native_turek", ".false.")
        input_text = set_namelist_value(input_text, "native_crosscheck", ".false.")
        input_text = set_namelist_value(input_text, "diagnostics", "'none'")
        input_text = set_namelist_value(input_text, "temperature", "350.0")
        input_text = set_namelist_value(input_text, "rotation_n_omega", "3")
        input_path.write_text(input_text)
        completed = subprocess.run(
            [str(args.binary.resolve()), "input.nml"],
            cwd=scratch,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=300,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"generic rotation run failed with exit code {completed.returncode}:\n{completed.stdout[-4000:]}")
        output = scratch / "rotation_dynamics.dat"
        if not output.is_file():
            raise RuntimeError(f"generic rotation run did not produce {output}")
        if "Rotation dynamics implementation = GENERIC" not in completed.stdout:
            raise AssertionError("generic workflow did not report the generic implementation")
        if "Fe validation campaign = NOT RUN" not in completed.stdout:
            raise AssertionError("generic workflow unexpectedly ran the Fe validation campaign")
        if "temperature = 350.0" not in input_path.read_text():
            raise AssertionError("generic fixture temperature was not set")
        grid = (scratch / "rotation_response_grid.dat").read_text()
        expected_grid_header = (
            "# q_fraction_x q_fraction_y q_fraction_z omega_Ry ReK_plus_Ry ImK_plus_Ry "
            "ReK_minus_Ry ImK_minus_Ry minus_Im_Kinv_plus_invRy minus_Im_Kinv_minus_invRy"
        )
        if expected_grid_header not in grid.splitlines():
            raise AssertionError("rotation grid spectral-weight header differs from the public schema")
        rows = 0
        for line in output.read_text().splitlines():
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            try:
                float(fields[0].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            rows += 1
        if rows == 0:
            raise AssertionError("generic rotation output contains no numeric response rows")
        print("generic rotation workflow regression: PASS")
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"generic rotation workflow regression: FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
