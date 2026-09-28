#!/usr/bin/env python3
"""Smoke-test the shipped bcc-Fe rotation example from a clean workdir."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


def require_all(text: str, expected: list[str], label: str) -> None:
    missing = [item for item in expected if item not in text]
    if missing:
        raise AssertionError(f"{label} is missing required text: {missing}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()

    scratch = args.scratch_root / "bccFe_rotation"
    try:
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(args.case.resolve(), scratch)

        input_text = (scratch / "input.nml").read_text()
        require_all(
            input_text,
            [
                "use_kspace = .true.",
                "nsp = 1",
                "hoh = .true.",
                "reciprocal_mode = 'ham_only'",
                "kspace_ham_order = 'second'",
                "formulation = 'rotation'",
            ],
            "bcc-Fe rotation input",
        )
        if not (scratch / "Fe.nml").is_file():
            raise AssertionError("shipped bcc-Fe rotation example has no Fe.nml")

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OMP_DYNAMIC"] = "FALSE"
        completed = subprocess.run(
            [str(args.binary.resolve()), "input.nml"],
            cwd=scratch,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=args.timeout,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"bcc-Fe rotation example failed with exit code {completed.returncode}:\n"
                f"{completed.stdout[-4000:]}"
            )

        output = scratch / "rotation_dynamics.dat"
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError("bcc-Fe rotation example did not produce rotation_dynamics.dat")
        output_text = output.read_text()
        require_all(
            output_text,
            [
                "# state_provenance = accepted second-order ham_only; SOC off; CCOR off; constraining field zero",
                "# k_mesh = 4 4 4",
            ],
            "rotation_dynamics.dat",
        )
        numeric_rows = 0
        for line in output_text.splitlines():
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            try:
                float(fields[0].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            numeric_rows += 1
        if numeric_rows == 0:
            raise AssertionError("rotation_dynamics.dat contains no numeric response rows")

        print(f"bcc-Fe rotation example smoke: PASS ({numeric_rows} response row(s))")
    except (AssertionError, OSError, RuntimeError, subprocess.TimeoutExpired, ValueError) as exc:
        print(f"bcc-Fe rotation example smoke: FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
