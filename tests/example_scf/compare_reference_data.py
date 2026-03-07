#!/usr/bin/env python3

from __future__ import annotations

import argparse
import fnmatch
import os
import sys
from typing import Iterable

from f90nml import read


def read_manifest(manifest_path: str) -> list[dict[str, str]]:
    cases: list[dict[str, str]] = []
    with open(manifest_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [item.strip() for item in line.split(",")]
            if len(parts) != 7:
                raise ValueError(f"Invalid manifest line: {line}")
            cases.append(
                {
                    "test_name": parts[0],
                    "case_rel": parts[1],
                    "nsp": parts[2],
                    "recur": parts[3],
                    "hoh": parts[4],
                    "lld": parts[5],
                    "nstep": parts[6],
                }
            )
    return cases


def select_cases(cases: Iterable[dict[str, str]], patterns: list[str], include_cheb: bool) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for case in cases:
        name = case["test_name"]
        if not include_cheb and case["recur"] == "chebyshev":
            continue
        if patterns and not any(fnmatch.fnmatch(name, pattern) for pattern in patterns):
            continue
        selected.append(case)
    return selected


def value_from_par(nml_path: str, key: str) -> float:
    data = read(nml_path)
    if "par" not in data:
        raise KeyError(f"Section 'par' missing in {nml_path}")
    if key not in data["par"]:
        raise KeyError(f"Key '{key}' missing in par section of {nml_path}")

    value = data["par"][key]
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            raise ValueError(f"Empty sequence for key '{key}' in {nml_path}")
        value = value[-1]

    return float(value)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare example SCF run outputs against stored references.")
    parser.add_argument("--manifest", default=os.path.join(os.path.dirname(__file__), "cases_manifest.csv"))
    parser.add_argument("--references-dir", default=os.path.join(os.path.dirname(__file__), "references"))
    parser.add_argument("--scratch-root", required=True, help="Directory where run_example_scf.sh wrote case outputs")
    parser.add_argument("--keys", default="etot,ws_r,vmad", help="Comma-separated keys in par section")
    parser.add_argument("--abs-tol", type=float, default=1e-6)
    parser.add_argument("--rel-tol", type=float, default=1e-6)
    parser.add_argument("--pattern", action="append", default=[], help="Glob pattern for test names")
    parser.add_argument("--include-cheb", action="store_true")
    args = parser.parse_args()

    keys = [item.strip() for item in args.keys.split(",") if item.strip()]
    if not keys:
        raise ValueError("At least one comparison key must be provided")

    cases = read_manifest(args.manifest)
    selected = select_cases(cases, args.pattern, args.include_cheb)
    if not selected:
        print("No cases selected for comparison")
        return 1

    failures: list[str] = []

    for case in selected:
        name = case["test_name"]
        run_data = os.path.join(args.scratch_root, name, "data.nml")
        ref_data = os.path.join(args.references_dir, name, "ref.nml")

        if not os.path.exists(run_data):
            failures.append(f"{name}: missing run output {run_data}")
            continue
        if not os.path.exists(ref_data):
            failures.append(f"{name}: missing reference {ref_data}")
            continue

        for key in keys:
            run_value = value_from_par(run_data, key)
            ref_value = value_from_par(ref_data, key)
            abs_diff = abs(run_value - ref_value)
            scale = max(abs(ref_value), 1.0)
            rel_diff = abs_diff / scale

            if abs_diff > args.abs_tol and rel_diff > args.rel_tol:
                failures.append(
                    f"{name}: key={key} run={run_value:.12e} ref={ref_value:.12e} "
                    f"abs_diff={abs_diff:.3e} rel_diff={rel_diff:.3e}"
                )

    if failures:
        print("Reference comparison FAILED")
        for item in failures:
            print(f" - {item}")
        return 1

    print(f"Reference comparison PASSED for {len(selected)} case(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
