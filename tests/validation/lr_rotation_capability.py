#!/usr/bin/env python3
"""Exercise the rotation capability gates and diagnostics-only Fe campaign."""

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


def prepare_case(case: Path, scratch: Path, updates: dict[str, str]) -> Path:
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(case.resolve(), scratch)
    input_path = scratch / "input.nml"
    input_text = input_path.read_text()
    for key, value in updates.items():
        input_text = set_namelist_value(input_text, key, value)
    input_path.write_text(input_text)
    return input_path


def run(binary: Path, scratch: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary.resolve()), "input.nml"],
        cwd=scratch,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        first_order_dir = args.scratch_root / "first_order_rejected"
        prepare_case(args.case, first_order_dir, {"kspace_ham_order": "'first'"})
        first_order = run(args.binary, first_order_dir)
        if first_order.returncode == 0:
            raise AssertionError("first-order rotation workflow unexpectedly succeeded")
        if "rotation capability gate requires kspace_ham_order=second" not in first_order.stdout:
            raise AssertionError(
                "first-order rotation failure did not identify the second-order capability gate:\n"
                f"{first_order.stdout[-4000:]}"
            )

        invariants_dir = args.scratch_root / "invariants_campaign"
        prepare_case(
            args.case,
            invariants_dir,
            {
                "native_turek": ".false.",
                "native_crosscheck": ".false.",
                "diagnostics": "'invariants'",
                "temperature": "300.0",
            },
        )
        invariants = run(args.binary, invariants_dir)
        if invariants.returncode != 0:
            raise RuntimeError(
                f"diagnostics=invariants rotation run failed with exit code {invariants.returncode}:\n"
                f"{invariants.stdout[-4000:]}"
            )
        if "Fe material convergence =" not in invariants.stdout:
            raise AssertionError("diagnostics=invariants did not execute the Fe campaign")
        if "Fe validation campaign = NOT RUN" in invariants.stdout:
            raise AssertionError("diagnostics=invariants reported the Fe campaign as disabled")

        print("rotation capability and diagnostics separation regression: PASS")
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"rotation capability and diagnostics separation regression: FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
