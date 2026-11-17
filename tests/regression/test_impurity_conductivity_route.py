#!/usr/bin/env python3
"""A bulk-impurity conductivity deck must reuse its newclubulk geometry.

Regression for post_processing_conductivity rebuilding a calctype='I' stack as
the surface route. Before the route is forwarded, the second stack calls
build_surf_full on the already generated bulk impurity cluster and fails while
reading it. The case explicitly keeps CCOR disabled and exercises an actual
quadrupole-accumulation output.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/scf/cases/impurity/B2FeCo_nsp1"


def prepare(work: Path) -> None:
    shutil.copytree(FIXTURE, work)
    deck = work / "input.nml"
    text = deck.read_text()
    text = text.replace(
        "pre_processing = 'newclubulk'",
        "pre_processing = 'newclubulk'\npost_processing = 'conductivity'",
        1,
    )
    text = text.replace(
        "calctype = 'I'",
        "calctype = 'I'\ncond_ll = 20\ncond_calctype = 'per_type'\n"
        "linear_in = 'charge'\nlinear_out = 'quadrupole_accumulation'",
        1,
    )
    text = text.replace("hoh = false", "hoh = false\nccor_2c = .false.\npol_alpha = 'Qxy'", 1)
    text = text.replace("hubbard_u_impurity(1,3) = 2.0\n", "")
    text = text.replace("hubbard_j_impurity(1,3) = 0.0\n", "")
    text += "\n&conductivity\nresponse_geometry = 'bulk'\n/\n"
    deck.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    work = args.scratch_root / "bulk_impurity_qxy"
    prepare(work)
    result = subprocess.run(
        [str(args.binary.resolve()), "input.nml"],
        cwd=work,
        text=True,
        capture_output=True,
        env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
    )
    if result.returncode != 0:
        print("FAIL bulk impurity conductivity route")
        print((result.stdout + result.stderr)[-5000:])
        return 1
    output = work / "response_total.out"
    if not output.exists():
        print("FAIL bulk impurity conductivity did not write response_total.out")
        return 1
    data = np.loadtxt(output)
    if data.ndim != 2 or data.shape[1] != 3 or not np.isfinite(data).all():
        print("FAIL bulk impurity response output has an invalid shape or non-finite values")
        return 1
    print(f"PASS bulk impurity conductivity route: {data.shape[0]} finite Qxy response rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
