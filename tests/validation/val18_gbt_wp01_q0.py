#!/usr/bin/env python3
"""VAL-18/WP01: fixed-potential bcc-Fe q=0 production identity.

The two runs start from the same committed Fe potential.  The ordinary
periodic-NC and GBT single-q builders are then exercised by the production
``kspace_green`` stack at q=0.  The resulting Green-function records are
compared as a compact application-level check; the standalone Fortran WP01
test owns the direct directed-block and H(k) operator metrics.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import f90nml


ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
FIXTURE = ROOT / "tests" / "gbt_wp01" / "cases" / "bccFe"
RUNNER = TESTS / "run_binary.sh"
FLOAT = r"[-+0-9.EeDd]+"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def run_case(binary: Path, workdir: Path, representation: str) -> tuple[Path, float]:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(FIXTURE, workdir)

    input_path = workdir / "input.nml"
    patched_path = workdir / "input.patched.nml"
    f90nml.patch(
        str(input_path),
        {
            "lattice": {"strux_backend": "strux_lib", "strux_want_sdot": False},
            "hamiltonian": {
                "magnetic_representation": representation,
                "q_ss": [0.0, 0.0, 0.0],
                "theta_ss": 0.0,
                "hoh": False,
            },
        },
        str(patched_path),
    )
    patched_path.replace(input_path)

    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"{representation} production run failed:\n{log[-6000:]}")
    if "fatal" in log.lower():
        raise RuntimeError(f"{representation} production run emitted a fatal diagnostic:\n{log[-6000:]}")

    links = []
    bond_path = workdir / "gbt_bonds.out"
    if bond_path.exists():
        for line in bond_path.read_text(errors="replace").splitlines():
            if not line.startswith("link "):
                continue
            fields = line.split()[1:]
            if len(fields) != 8:
                raise RuntimeError(f"malformed GBT link record: {line}")
            links.append([number(field) for field in fields])

    if representation == "gbt_single_q":
        if not links:
            raise RuntimeError("GBT production run did not write any link records")
        identity_error = 0.0
        for link in links:
            # Complex values are written as real/imaginary pairs.
            identity_error = max(
                identity_error,
                abs(link[0] - 1.0), abs(link[1]),
                abs(link[2]), abs(link[3]), abs(link[4]), abs(link[5]),
                abs(link[6] - 1.0), abs(link[7]),
            )
    else:
        identity_error = 0.0
        if links:
            raise RuntimeError("ordinary production run unexpectedly wrote gbt_bonds.out")

    return workdir, identity_error


def numeric_records(path: Path) -> list[list[float]]:
    records: list[list[float]] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        try:
            records.append([number(field) for field in fields])
        except ValueError:
            raise RuntimeError(f"non-numeric data record in {path}: {line}")
    return records


def compare_outputs(ordinary: Path, gbt: Path) -> tuple[float, float, int]:
    files = ("kspace_green_c1.dat", "kspace_green_gf.dat")
    max_abs = 0.0
    max_rel = 0.0
    count = 0
    for filename in files:
        left_path = ordinary / filename
        right_path = gbt / filename
        if not left_path.exists() or not right_path.exists():
            raise RuntimeError(f"missing common production output {filename}")
        left = numeric_records(left_path)
        right = numeric_records(right_path)
        if len(left) != len(right):
            raise RuntimeError(f"{filename} row count differs: {len(left)} != {len(right)}")
        for row_index, (left_row, right_row) in enumerate(zip(left, right), start=1):
            if len(left_row) != len(right_row):
                raise RuntimeError(f"{filename} column count differs on row {row_index}")
            for left_value, right_value in zip(left_row, right_row):
                difference = abs(left_value - right_value)
                max_abs = max(max_abs, difference)
                max_rel = max(max_rel, difference / max(1.0, abs(left_value)))
                count += 1
    return max_abs, max_rel, count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    ordinary, ordinary_link_error = run_case(
        args.binary.resolve(), args.scratch_root.resolve() / "ordinary", "periodic_nc"
    )
    gbt, gbt_link_error = run_case(
        args.binary.resolve(), args.scratch_root.resolve() / "gbt", "gbt_single_q"
    )
    max_abs, max_rel, count = compare_outputs(ordinary, gbt)

    print(f"WP01 bcc-Fe q=0 production records compared: {count}")
    print(f"WP01 bcc-Fe q=0 production residual (abs, rel): {max_abs:.6e} {max_rel:.6e}")
    print(f"WP01 bcc-Fe GBT q=0 link identity residual: {gbt_link_error:.6e}")
    if ordinary_link_error > args.tolerance or gbt_link_error > args.tolerance:
        raise RuntimeError("unexpected link records in the production q=0 fixture")
    if max_abs > args.tolerance and max_rel > args.tolerance:
        raise RuntimeError(
            f"ordinary and GBT q=0 production records differ beyond {args.tolerance:.3e}"
        )
    print("WP01 bcc-Fe q=0 production fixture: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
