#!/usr/bin/env python3
"""
Unified runner for RSLMTO example tests.

Reads a test case from a cases.json manifest, sets up a scratch directory,
patches input.nml with the test parameters, runs the binary, checks for
fatal errors, and optionally compares or saves reference data.

Usage (smoke):
  run_test.py --binary BIN --cases-json JSON --case-name NAME --scratch-root DIR

Usage (reference comparison):
  run_test.py ... --compare-ref REF_DIR [--keys k1,k2] [--abs-tol 1e-6] [--rel-tol 1e-6]

Usage (save reference):
  run_test.py ... --gen-ref REF_DIR
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

import f90nml


# ---------------------------------------------------------------------------
# Case loading
# ---------------------------------------------------------------------------

def load_case(cases_json: str, case_name: str) -> dict:
    with open(cases_json) as fh:
        data = json.load(fh)
    for case in data["cases"]:
        if case["name"] == case_name:
            return case
    raise KeyError(f"Case '{case_name}' not found in {cases_json}")


# ---------------------------------------------------------------------------
# Scratch directory setup
# ---------------------------------------------------------------------------

def setup_scratch(case_dir: str, workdir: str) -> None:
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    shutil.copytree(case_dir, workdir)
    for pattern in ("*_out.nml", "data.nml"):
        for path in glob.glob(os.path.join(workdir, pattern)):
            os.remove(path)


# ---------------------------------------------------------------------------
# NML patching
# ---------------------------------------------------------------------------

def patch_input_nml(workdir: str, case: dict) -> None:
    input_path = os.path.join(workdir, "input.nml")
    tmp_path = input_path + ".tmp"
    f90nml.patch(input_path, case["namelists"], tmp_path)
    os.replace(tmp_path, input_path)


# ---------------------------------------------------------------------------
# Binary invocation
# ---------------------------------------------------------------------------

def run_binary(binary: str, workdir: str) -> None:
    run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_binary.sh")
    result = subprocess.run(
        ["/bin/bash", run_script, binary],
        cwd=workdir,
    )
    if result.returncode != 0:
        print(f"ERROR: binary returned non-zero exit code {result.returncode}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Output checks and normalization
# ---------------------------------------------------------------------------

def check_log(workdir: str, case_name: str) -> None:
    log_path = os.path.join(workdir, "testrun.log")
    if not os.path.exists(log_path):
        print(f"ERROR [{case_name}]: testrun.log missing in {workdir}")
        sys.exit(1)
    with open(log_path) as fh:
        content = fh.read()
    if "fatal" in content.lower():
        print(f"ERROR [{case_name}]: fatal message in testrun.log")
        lines = content.splitlines()
        print("\n".join(lines[-50:]))
        sys.exit(1)


def normalize_output(workdir: str, case_name: str) -> str:
    candidates = sorted(glob.glob(os.path.join(workdir, "*_out.nml")))
    if not candidates:
        print(f"ERROR [{case_name}]: no *_out.nml found in {workdir}")
        sys.exit(1)
    if len(candidates) > 1:
        print(f"NOTE [{case_name}]: multiple *_out.nml, using {os.path.basename(candidates[0])}")
    data_path = os.path.join(workdir, "data.nml")
    os.replace(candidates[0], data_path)
    return data_path


# ---------------------------------------------------------------------------
# Reference comparison
# ---------------------------------------------------------------------------

def _read_par_value(nml_path: str, key: str) -> float:
    data = f90nml.read(nml_path)
    if "par" not in data:
        raise KeyError(f"Section 'par' missing in {nml_path}")
    value = data["par"][key]
    if isinstance(value, (list, tuple)):
        value = value[-1]
    return float(value)


def compare_ref(
    data_path: str,
    case_name: str,
    ref_dir: str,
    keys: list[str],
    abs_tol: float,
    rel_tol: float,
) -> None:
    ref_path = os.path.join(ref_dir, case_name, "ref.nml")
    if not os.path.exists(ref_path):
        print(f"FAIL [{case_name}]: missing reference {ref_path}")
        sys.exit(1)

    failures = []
    for key in keys:
        run_val = _read_par_value(data_path, key)
        ref_val = _read_par_value(ref_path, key)
        abs_diff = abs(run_val - ref_val)
        scale = max(abs(ref_val), 1.0)
        rel_diff = abs_diff / scale
        if abs_diff > abs_tol and rel_diff > rel_tol:
            failures.append(
                f"  key={key}  run={run_val:.12e}  ref={ref_val:.12e}"
                f"  |diff|={abs_diff:.3e}  rel={rel_diff:.3e}"
            )

    if failures:
        print(f"FAIL [{case_name}]: reference comparison failed")
        for msg in failures:
            print(msg)
        sys.exit(1)

    print(f"PASS [{case_name}]: reference OK ({len(keys)} keys)")


# ---------------------------------------------------------------------------
# Reference generation
# ---------------------------------------------------------------------------

def save_ref(data_path: str, case_name: str, ref_dir: str, case: dict) -> None:
    dest_dir = os.path.join(ref_dir, case_name)
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy2(data_path, os.path.join(dest_dir, "ref.nml"))
    with open(os.path.join(dest_dir, "meta.json"), "w") as fh:
        json.dump(case, fh, indent=2)
    print(f"REF  [{case_name}]: saved to {dest_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Run one RSLMTO example test case.")
    parser.add_argument("--binary", required=True, help="Path to rslmto.x binary")
    parser.add_argument("--cases-json", required=True, help="Path to cases.json manifest")
    parser.add_argument("--case-name", required=True, help="Name of the case to run")
    parser.add_argument("--scratch-root", required=True, help="Root directory for scratch runs")

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--compare-ref", metavar="REF_DIR", help="Compare output against stored reference")
    mode.add_argument("--gen-ref", metavar="REF_DIR", help="Save output as new reference")

    parser.add_argument("--keys", default="etot,ws_r,vmad", help="Keys in &par to compare")
    parser.add_argument("--abs-tol", type=float, default=1e-6)
    parser.add_argument("--rel-tol", type=float, default=1e-6)
    args = parser.parse_args()

    case = load_case(args.cases_json, args.case_name)
    binary = os.path.abspath(args.binary)
    cases_dir = os.path.join(os.path.dirname(os.path.abspath(args.cases_json)), "cases")
    case_dir = os.path.join(cases_dir, case["case"])
    workdir = os.path.join(args.scratch_root, args.case_name)

    setup_scratch(case_dir, workdir)
    patch_input_nml(workdir, case)
    run_binary(binary, workdir)
    check_log(workdir, args.case_name)
    data_path = normalize_output(workdir, args.case_name)

    if args.compare_ref:
        keys = [k.strip() for k in args.keys.split(",") if k.strip()]
        compare_ref(data_path, args.case_name, args.compare_ref, keys, args.abs_tol, args.rel_tol)
    elif args.gen_ref:
        save_ref(data_path, args.case_name, args.gen_ref, case)
    else:
        print(f"PASS [{args.case_name}]")


if __name__ == "__main__":
    main()
