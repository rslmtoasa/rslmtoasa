#!/usr/bin/env python3
"""
Batch reference generator for RSLMTO example tests.

Runs selected cases from a cases.json manifest using a known-good binary
and saves the outputs as reference data for future comparisons.

Usage:
  generate_references.py --binary BIN --cases-json JSON --references-dir DIR
                         [--scratch-root DIR] [--case NAME [NAME ...]]

Examples:
  # Generate all SCF references
  generate_references.py --binary build/bin/rslmto.x \\
      --cases-json tests/scf/cases.json \\
      --references-dir tests/scf/references

  # Only specific cases
  generate_references.py ... --case Example_bulk_bccFe_nsp2_block_hoh_true
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile


def load_cases(cases_json: str, selected: list[str]) -> list[dict]:
    with open(cases_json) as fh:
        data = json.load(fh)
    cases = data["cases"]
    if selected:
        cases = [c for c in cases if c["name"] in selected]
    return cases


def run_case(
    binary: str,
    cases_json: str,
    case_name: str,
    scratch_root: str,
    ref_dir: str,
) -> bool:
    run_test = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_test.py")
    cmd = [
        sys.executable, run_test,
        "--binary", binary,
        "--cases-json", cases_json,
        "--case-name", case_name,
        "--scratch-root", scratch_root,
        "--gen-ref", ref_dir,
    ]
    result = subprocess.run(cmd)
    return result.returncode == 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate RSLMTO example test references.")
    parser.add_argument("--binary", required=True, help="Path to known-good rslmto.x binary")
    parser.add_argument("--cases-json", required=True, help="Path to cases.json manifest")
    parser.add_argument("--references-dir", required=True, help="Where to write reference data")
    parser.add_argument("--scratch-root", default=None, help="Scratch directory (default: temp dir)")
    parser.add_argument("--case", nargs="+", default=[], metavar="NAME",
                        help="Run only these named cases")
    args = parser.parse_args()

    cases = load_cases(args.cases_json, args.case)
    if not cases:
        print("No cases selected.")
        sys.exit(1)

    use_tmp = args.scratch_root is None
    scratch = args.scratch_root or tempfile.mkdtemp(prefix="rslmto_refgen_")
    os.makedirs(args.references_dir, exist_ok=True)

    failures: list[str] = []
    for case in cases:
        name = case["name"]
        print(f"RUN  [{name}]")
        ok = run_case(args.binary, args.cases_json, name, scratch, args.references_dir)
        if not ok:
            failures.append(name)

    if use_tmp:
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)

    if failures:
        print(f"\nReference generation FAILED for {len(failures)} case(s):")
        for name in failures:
            print(f"  {name}")
        sys.exit(1)

    print(f"\nReference generation complete: {len(cases) - len(failures)} case(s) written to {args.references_dir}")


if __name__ == "__main__":
    main()
