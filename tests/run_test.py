#!/usr/bin/env python3
"""
Unified runner for RSLMTO example tests.

Reads a test case from a cases.json manifest, sets up a scratch directory,
patches input.nml with the test parameters, runs the binary, checks for
fatal errors, and optionally compares or saves reference data.

Usage (smoke):
  run_test.py --binary BIN --cases-json JSON --case-name NAME --scratch-root DIR

Usage (reference comparison):
  run_test.py ... --compare-ref REF_DIR [--abs-tol 1e-6] [--rel-tol 1e-6]

Usage (save reference):
  run_test.py ... --gen-ref REF_DIR

Reference data is driven by a "checks" dict in cases.json:

  "checks": {
    "nml": [
      {
        "file": "Fe_out.nml",
        "scalars": ["etot", "ws_r", "vmad"],
        "arrays": { "moment": [1, 2], "dos": [1, 5, 10] }
      }
    ],
    "text": [
      { "file": "Fe_dos.out", "rows": [50, 100], "cols": [1, 2] }
    ]
  }

If a case has no "checks" key the test still runs as a smoke test.
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

def run_binary(binary: str, workdir: str, mpi_procs: int = 1) -> None:
    run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_binary.sh")
    cmd = ["/bin/bash", run_script, binary]
    if mpi_procs > 1:
        cmd.append(str(mpi_procs))
    result = subprocess.run(cmd, cwd=workdir)
    if result.returncode != 0:
        print(f"ERROR: binary returned non-zero exit code {result.returncode}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Post-run checks
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



# ---------------------------------------------------------------------------
# Value extraction
# ---------------------------------------------------------------------------

def _nml_flat(filepath: str) -> dict:
    """Read a namelist file and return a flat key->value dict (all sections merged)."""
    nml = f90nml.read(filepath)
    flat: dict = {}
    for section in nml.values():
        if isinstance(section, dict):
            flat.update(section)
    return flat


def extract_nml_values(workdir: str, nml_check: dict) -> dict:
    """Extract scalars and array elements from a namelist output file."""
    filepath = os.path.join(workdir, nml_check["file"])
    flat = _nml_flat(filepath)
    result: dict = {}

    for key in nml_check.get("scalars", []):
        result[key] = float(flat[key])

    for key, indices in nml_check.get("arrays", {}).items():
        arr = flat[key]  # f90nml returns a Python list (0-based)
        result[key] = {str(i): float(arr[i - 1]) for i in indices}

    return result


def extract_text_values(workdir: str, text_check: dict) -> dict:
    """Extract values at specified (row, col) positions from a space-separated file."""
    filepath = os.path.join(workdir, text_check["file"])
    rows = text_check["rows"]
    cols = text_check["cols"]
    row_set = set(rows)

    result: dict = {}
    with open(filepath) as fh:
        for lineno, line in enumerate(fh, start=1):
            if lineno in row_set:
                parts = line.split()
                result[str(lineno)] = {str(c): float(parts[c - 1]) for c in cols}

    return result


def build_ref_data(workdir: str, checks: dict) -> dict:
    """Build a complete reference data dict from the workdir outputs."""
    ref: dict = {}

    if "nml" in checks:
        ref["nml"] = {}
        for nml_check in checks["nml"]:
            ref["nml"][nml_check["file"]] = extract_nml_values(workdir, nml_check)

    if "text" in checks:
        ref["text"] = {}
        for text_check in checks["text"]:
            ref["text"][text_check["file"]] = extract_text_values(workdir, text_check)

    return ref


# ---------------------------------------------------------------------------
# Reference comparison
# ---------------------------------------------------------------------------

def _check_value(
    failures: list,
    label: str,
    run_v: float | None,
    ref_v: float,
    abs_tol: float,
    rel_tol: float,
) -> None:
    if run_v is None:
        failures.append(f"  {label}: missing in run output")
        return
    abs_diff = abs(run_v - ref_v)
    scale = max(abs(ref_v), 1.0)
    rel_diff = abs_diff / scale
    if abs_diff > abs_tol and rel_diff > rel_tol:
        failures.append(
            f"  {label}  run={run_v:.12e}  ref={ref_v:.12e}"
            f"  |diff|={abs_diff:.3e}  rel={rel_diff:.3e}"
        )


def compare_ref(
    workdir: str,
    case_name: str,
    ref_dir: str,
    checks: dict,
    abs_tol: float,
    rel_tol: float,
) -> None:
    if not checks:
        print(f"PASS [{case_name}]: no checks defined (smoke only)")
        return

    ref_path = os.path.join(ref_dir, case_name, "ref.json")
    if not os.path.exists(ref_path):
        print(f"PASS [{case_name}]: no reference found, smoke only (run --gen-ref to create)")
        return

    with open(ref_path) as fh:
        ref_data = json.load(fh)

    run_data = build_ref_data(workdir, checks)
    failures: list[str] = []
    n_checked = 0

    # NML comparisons
    for filename, ref_vals in ref_data.get("nml", {}).items():
        run_vals = run_data.get("nml", {}).get(filename, {})
        for key, ref_val in ref_vals.items():
            if isinstance(ref_val, dict):
                run_val = run_vals.get(key, {})
                for idx, ref_v in ref_val.items():
                    run_v = run_val.get(idx) if isinstance(run_val, dict) else None
                    _check_value(failures, f"{filename}:{key}[{idx}]", run_v, ref_v, abs_tol, rel_tol)
                    n_checked += 1
            else:
                _check_value(failures, f"{filename}:{key}", run_vals.get(key), ref_val, abs_tol, rel_tol)
                n_checked += 1

    # Text comparisons
    for filename, ref_rows in ref_data.get("text", {}).items():
        run_rows = run_data.get("text", {}).get(filename, {})
        for row_str, ref_cols in ref_rows.items():
            run_cols = run_rows.get(row_str, {})
            for col_str, ref_v in ref_cols.items():
                _check_value(
                    failures,
                    f"{filename}:row{row_str}:col{col_str}",
                    run_cols.get(col_str),
                    ref_v,
                    abs_tol,
                    rel_tol,
                )
                n_checked += 1

    if failures:
        print(f"FAIL [{case_name}]: {len(failures)} of {n_checked} value(s) out of tolerance")
        for msg in failures:
            print(msg)
        sys.exit(1)

    print(f"PASS [{case_name}]: reference OK ({n_checked} values checked)")


# ---------------------------------------------------------------------------
# Reference generation
# ---------------------------------------------------------------------------

def save_ref(workdir: str, case_name: str, ref_dir: str, case: dict) -> None:
    checks = case.get("checks", {})
    dest_dir = os.path.join(ref_dir, case_name)
    os.makedirs(dest_dir, exist_ok=True)

    ref_data = build_ref_data(workdir, checks)
    with open(os.path.join(dest_dir, "ref.json"), "w") as fh:
        json.dump(ref_data, fh, indent=2)
    with open(os.path.join(dest_dir, "meta.json"), "w") as fh:
        json.dump(case, fh, indent=2)

    n = sum(
        (len(v) if isinstance(v, dict) else 1)
        for nml_vals in ref_data.get("nml", {}).values()
        for v in nml_vals.values()
    ) + sum(
        len(cols)
        for rows in ref_data.get("text", {}).values()
        for cols in rows.values()
    )
    print(f"REF  [{case_name}]: {n} values saved to {dest_dir}")


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
    mode.add_argument("--gen-ref",     metavar="REF_DIR", help="Save output as new reference")

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
    run_binary(binary, workdir, case.get("mpi_procs", 1))
    check_log(workdir, args.case_name)

    if args.compare_ref:
        compare_ref(
            workdir, args.case_name, args.compare_ref,
            case.get("checks", {}), args.abs_tol, args.rel_tol,
        )
    elif args.gen_ref:
        save_ref(workdir, args.case_name, args.gen_ref, case)
    else:
        print(f"PASS [{args.case_name}]")


if __name__ == "__main__":
    main()
