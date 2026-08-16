#!/usr/bin/env python3
"""VAL-07: validate the production exchange tensor and two-index routes.

The bcc-Fe fixture has two distinct pair records.  This campaign checks only
the relationships justified by ``docs/dev/EXCHANGE_VALIDATION_MAP.md``:
full-tensor assembly and the documented two-index J/D/A recombinations.

The Gauss-Legendre, auxiliary-GF, native-import, and spin-lattice routes are
not silently compared to this fixture.  VAL-06 records that their inputs,
representations, or physical quantities are not currently matched, or that
the route has no caller; they remain Experimental in the maturity ledger.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path


FLOAT = r"[-+0-9.EeDd]+"
PAIR_COUNT = 2
ALGEBRA_TOL = 5.0e-5  # parts files are printed with six digits after decimal
TENSOR_TOL = 4.0e-6   # full tensor and A/J files are printed with six digits


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def run_fixture(binary: Path, root: Path, scratch_root: Path) -> Path:
    workdir = scratch_root / "Example_exchange_bccFe"
    if workdir.exists():
        shutil.rmtree(workdir)
    fixture = root / "tests/postproc/cases/exchange/bccFe"
    shutil.copytree(fixture, workdir)
    input_path = workdir / "input.nml"
    text = input_path.read_text()
    text = text.replace("&lattice\n", "&lattice\nstrux_backend = 'strux_lib'\n", 1)
    text = re.sub(r"(?im)^\s*nstep\s*=\s*\d+", "nstep = 1", text)
    text = re.sub(r"(?im)^\s*lld\s*=\s*\d+", "lld = 20", text)
    input_path.write_text(text)

    result = subprocess.run(
        ["/bin/bash", str(root / "tests/run_binary.sh"), str(binary)],
        cwd=workdir, capture_output=True, text=True,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"VAL-07 exchange fixture failed:\n{log[-5000:]}")
    if "fatal" in log.lower():
        raise RuntimeError(f"VAL-07 exchange fixture emitted a fatal diagnostic:\n{log[-5000:]}")
    return workdir


def rows(path: Path, columns: int) -> list[list[float]]:
    result: list[list[float]] = []
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != columns:
            raise RuntimeError(f"{path.name}:{line_number}: expected {columns} columns, got {len(fields)}")
        values = [number(field) for field in fields]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"{path.name}:{line_number}: non-finite output")
        result.append(values)
    if len(result) != PAIR_COUNT:
        raise RuntimeError(f"{path.name}: expected {PAIR_COUNT} pair rows, got {len(result)}")
    return result


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("VAL-07 FAIL: " + message)


def close(actual: float, expected: float, tolerance: float) -> bool:
    return abs(actual - expected) <= tolerance


def check_pair_identity(reference: list[float], candidate: list[float], label: str) -> None:
    require(reference[:5] == candidate[:5], f"{label}: pair/site geometry changed")
    require(close(reference[-1], candidate[-1], TENSOR_TOL), f"{label}: bond distance changed")


def parse_tensors(log: str) -> list[list[list[float]]]:
    pattern = re.compile(
        rf"J tensor between pair\s+(\d+)\s+and\s+(\d+)\s+is\s*\n"
        rf"\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\n"
        rf"\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\n"
        rf"\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        re.MULTILINE,
    )
    tensors: list[list[list[float]]] = []
    for match in pattern.finditer(log):
        values = [number(value) for value in match.groups()[2:]]
        tensors.append([
            values[0:3], values[3:6], values[6:9]
        ])
    return tensors


def transpose_flattened(values: list[float]) -> list[list[float]]:
    """Convert Fortran array output order to the displayed 3x3 matrix."""
    return [[values[row + 3 * column] for column in range(3)] for row in range(3)]


def check_outputs(workdir: Path) -> dict[str, float]:
    j = rows(workdir / "jij.out", 7)
    d = rows(workdir / "dij.out", 9)
    a = rows(workdir / "aij.out", 15)
    jso = rows(workdir / "jijso.out", 7)
    jfo = rows(workdir / "jijfo.out", 7)
    dso = rows(workdir / "dijso.out", 9)
    dfo = rows(workdir / "dijfo.out", 9)
    aso = rows(workdir / "aijso.out", 15)
    afo = rows(workdir / "aijfo.out", 15)
    jparts = rows(workdir / "jijparts.out", 10)
    dparts = rows(workdir / "dijparts.out", 12)
    aparts = rows(workdir / "aijparts.out", 24)

    for other in (d, a, jso, jfo, dso, dfo, aso, afo, jparts, dparts, aparts):
        for index, reference in enumerate(j):
            check_pair_identity(reference, other[index], f"row {index + 1}")

    # The two rows are intentionally different pair records.  This catches a
    # site-index/communication error that could otherwise pass algebraically.
    require(j[0][2:5] == [-0.5, -0.5, -0.5], "first pair vector was not preserved")
    require(j[1][2:5] == [0.0, 0.0, -1.0], "second pair vector was not preserved")
    require(abs(j[0][5] - j[1][5]) > 1.0e-3, "distinct pair records collapsed to one J value")

    max_j_error = max_d_error = max_a_error = 0.0
    for index in range(PAIR_COUNT):
        # VAL-06: J_SO = CD - SD + CC - SC; J_FO = CD + SD - CC - SC.
        cd, sd, cc, sc = jparts[index][5:9]
        expected_jso = cd - sd + cc - sc
        expected_jfo = cd + sd - cc - sc
        max_j_error = max(max_j_error, abs(jso[index][5] - expected_jso), abs(jfo[index][5] - expected_jfo))
        require(close(jso[index][5], expected_jso, ALGEBRA_TOL), f"row {index + 1}: J_SO recombination")
        require(close(jfo[index][5], expected_jfo, ALGEBRA_TOL), f"row {index + 1}: J_FO recombination")

        # VAL-06: D_SO = 2*(SC + CC); D_FO = 2*(SC - CC).
        dcc = dparts[index][5:8]
        dsc = dparts[index][8:11]
        for component in range(3):
            expected_dso = 2.0 * (dsc[component] + dcc[component])
            expected_dfo = 2.0 * (dsc[component] - dcc[component])
            max_d_error = max(max_d_error, abs(dso[index][5 + component] - expected_dso), abs(dfo[index][5 + component] - expected_dfo))
            require(close(dso[index][5 + component], expected_dso, ALGEBRA_TOL), f"row {index + 1}: D_SO component {component + 1}")
            require(close(dfo[index][5 + component], expected_dfo, ALGEBRA_TOL), f"row {index + 1}: D_FO component {component + 1}")

        # VAL-06: A_SO = SD + SC; A_FO = -SD + SC.
        asd = aparts[index][5:14]
        asc = aparts[index][14:23]
        for component in range(9):
            expected_aso = asd[component] + asc[component]
            expected_afo = -asd[component] + asc[component]
            max_a_error = max(max_a_error, abs(aso[index][5 + component] - expected_aso), abs(afo[index][5 + component] - expected_afo))
            require(close(aso[index][5 + component], expected_aso, ALGEBRA_TOL), f"row {index + 1}: A_SO component {component + 1}")
            require(close(afo[index][5 + component], expected_afo, ALGEBRA_TOL), f"row {index + 1}: A_FO component {component + 1}")

    # Exact source convention: T = J I + D^skew + A.  Do not infer D from
    # antisym(T), since the raw A tensor is not generally symmetrized.
    tensors = parse_tensors((workdir / "testrun.log").read_text(errors="replace"))
    require(len(tensors) == PAIR_COUNT, "missing one or more printed full tensors")
    for index, row in enumerate(j):
        require(index < len(tensors), f"missing tensor for pair row {index + 1}")
        dvec = d[index][5:8]
        amatrix = transpose_flattened(a[index][5:14])
        expected = [[amatrix[r][c] for c in range(3)] for r in range(3)]
        for axis in range(3):
            expected[axis][axis] += row[5]
        expected[0][1] += dvec[2]
        expected[1][0] -= dvec[2]
        expected[0][2] -= dvec[1]
        expected[2][0] += dvec[1]
        expected[1][2] += dvec[0]
        expected[2][1] -= dvec[0]
        for r in range(3):
            for c in range(3):
                require(close(tensors[index][r][c], expected[r][c], TENSOR_TOL), f"row {index + 1}: full tensor assembly ({r + 1},{c + 1})")

        # This scalar-relativistic bcc-Fe fixture is a DMI-zero limit, not a
        # noncentrosymmetric SOC validation of a nonzero D vector.
        require(max(abs(value) for value in dvec) < 1.0e-5, f"row {index + 1}: unexpected DMI in bcc-Fe limit")

    # A nonzero symmetric-traceless component prevents this from becoming a
    # scalar-only test of the anisotropic output path.
    first_a = transpose_flattened(a[0][5:14])
    trace = sum(first_a[axis][axis] for axis in range(3)) / 3.0
    symmetric_traceless = [
        [0.5 * (first_a[r][c] + first_a[c][r]) - (trace if r == c else 0.0) for c in range(3)]
        for r in range(3)
    ]
    require(max(abs(value) for row in symmetric_traceless for value in row) > 1.0e-3, "anisotropic symmetric-traceless output is absent")
    return {"max_j_error": max_j_error, "max_d_error": max_d_error, "max_a_error": max_a_error}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    workdir = run_fixture(args.binary.resolve(), root, args.scratch_root.resolve())
    metrics = check_outputs(workdir)
    print(
        "VAL-07 PASS: canonical J/D/A tensor and two-index SO/FO identities "
        f"for bcc-Fe two-pair fixture "
        f"(max errors J={metrics['max_j_error']:.3e}, "
        f"D={metrics['max_d_error']:.3e}, A={metrics['max_a_error']:.3e})"
    )
    print("VAL-07 scope: Gauss-Legendre, auxiliary-GF, native-import, and jijk routes remain Experimental pending matched oracles/callers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
