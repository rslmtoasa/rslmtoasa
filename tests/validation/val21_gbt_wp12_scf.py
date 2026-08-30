#!/usr/bin/env python3
"""WP12 controlled q=0 SCF trace and ordinary/GBT closure check.

This is deliberately a small diagnostic campaign.  It runs the same bcc-Fe
two-iteration k-space SCF from one potential for periodic-NC and q=0 GBT,
then checks that the end-to-end trace is complete and the converged state is
the same within the fixed-potential q=0 tolerance.  Finite-q SCF is not
claimed here because the repository's production supercell counterpart is
still unmatched (VAL-16/WP04 evidence).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import f90nml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "gbt_wp01" / "cases" / "bccFe"
RUNNER = ROOT / "tests" / "run_binary.sh"
TRACE = "gbt_scf_diagnostics.dat"
TRACE_NUMERIC_COLUMNS = 51


def patch_case(destination: Path, representation: str) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(FIXTURE, destination)
    input_path = destination / "input.nml"
    patched = destination / "input.patched.nml"
    f90nml.patch(
        str(input_path),
        {
            "calculation": {"pre_processing": "bravais", "post_processing": "none"},
            "self": {"nstep": 2, "use_kspace": True, "conv_thr": 1.0e-10},
            "control": {
                "density_policy": "constrained_spiral",
                "gbt_scf_diagnostics": True,
            },
            "mix": {"beta": 0.2, "mixtype": "linear"},
            "hamiltonian": {
                "magnetic_representation": representation,
                "q_ss": [0.0, 0.0, 0.0],
                "theta_ss": 0.0,
            },
        },
        str(patched),
    )
    patched.replace(input_path)


def run_case(binary: Path, destination: Path, representation: str) -> list[list[str]]:
    patch_case(destination, representation)
    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=destination,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=900,
        check=False,
    )
    log = result.stdout
    (destination / "run.log").write_text(log, encoding="utf-8")
    if result.returncode != 0 or "fatal" in log.lower():
        raise RuntimeError(f"{representation} q=0 SCF failed:\n{log[-6000:]}")

    trace = destination / TRACE
    if not trace.exists():
        raise RuntimeError(f"{representation} q=0 SCF did not emit {TRACE}")
    rows: list[list[str]] = []
    for line in trace.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 4 + TRACE_NUMERIC_COLUMNS:
            raise RuntimeError(f"malformed WP12 trace row ({len(fields)} fields): {line}")
        rows.append(fields)
    if not rows:
        raise RuntimeError(f"empty WP12 trace in {destination}")
    return rows


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=2.0e-8)
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    ordinary = run_case(args.binary.resolve(), args.scratch_root / "ordinary", "periodic_nc")
    gbt = run_case(args.binary.resolve(), args.scratch_root / "gbt", "gbt_single_q")
    if len(ordinary) != len(gbt):
        raise RuntimeError(f"ordinary/GBT trace row count differs: {len(ordinary)} != {len(gbt)}")

    max_abs = 0.0
    for left, right in zip(ordinary, gbt):
        if left[:2] != right[:2]:
            raise RuntimeError(f"ordinary/GBT trace site ordering differs: {left[:2]} != {right[:2]}")
        for left_value, right_value in zip(left[4:], right[4:]):
            max_abs = max(max_abs, abs(number(left_value) - number(right_value)))

    iterations = sorted({int(row[0]) for row in gbt})
    print(f"WP12 q=0 trace rows: {len(gbt)}; iterations: {iterations}")
    print(f"WP12 q=0 ordinary/GBT maximum trace residual: {max_abs:.6e}")
    if max_abs > args.tolerance:
        raise RuntimeError(f"ordinary and GBT q=0 traces differ beyond {args.tolerance:.3e}")
    print("WP12 q=0 SCF diagnostic/closure fixture: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
