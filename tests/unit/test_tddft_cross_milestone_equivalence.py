#!/usr/bin/env python3
"""Single CI gate for numerical-equivalence coverage of TDDFT-03 through -10.

Each Fortran test owns an independent numerical oracle.  This runner makes
the required cross-milestone contract explicit and keeps performance changes
from silently omitting a response representation or backend.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

MILESTONES = (
    ("TDDFT-03", "bare chi_KS scalar/GEMM equivalence", "UnitTddftChiKS", "256 eps relative"),
    ("TDDFT-04", "ALSDA/Goldstone diagnostics", "UnitTddftGoldstone", "2e-12 absolute/relative"),
    ("TDDFT-05", "Dyson/Xi/mode solve", "UnitTddftDysonModes", "1e-12 absolute/relative"),
    ("TDDFT-06", "production configuration and dispatch", "UnitTddftConfig", "exact configuration contract"),
    ("TDDFT-07", "transverse validation evidence", None, "fixture acceptance limits"),
    ("TDDFT-08", "longitudinal calibration", "UnitTddftLongitudinal", "1e-11 absolute"),
    ("TDDFT-09", "four-component/SOC reduction", "UnitTddftFourComponent", "2e-12 absolute/relative"),
    ("TDDFT-10", "Green/eigenpair chi_KS backend", "UnitTddftGreenChiKS", "backend fixture tolerance"),
)


def run(command: list[str], label: str) -> None:
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.returncode:
        print(f"{label}: FAIL\n{completed.stdout}")
        raise SystemExit(completed.returncode)
    if "RESULT: PASS" not in completed.stdout:
        print(f"{label}: FAIL (missing PASS marker)\n{completed.stdout}")
        raise SystemExit(1)
    print(f"{label}: PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary-dir", required=True, type=Path)
    args = parser.parse_args()
    for milestone, purpose, executable, tolerance in MILESTONES:
        print(f"{milestone}: {purpose}; tolerance={tolerance}")
        if executable is None:
            run([sys.executable, str(ROOT / "regression" / "tddft_validation" / "test_validation.py")], milestone)
        else:
            run([str(args.binary_dir / executable)], milestone)
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
