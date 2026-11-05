#!/usr/bin/env python3
"""AppGbtKspaceRsFreeze: real-space GBT spiral and frozen-magnon cases against committed references.

Relation: none (regression freeze). Reference: values committed in tests/gbt/rs_freeze/references,
generated once from the C5 commit, compared with tolerance 0 (equal printed digits).
Catches: any change that reaches the real-space spiral route (hopping phase, frame link, W/C, sign).
Exit 77 (CTest skip) when this platform has no committed reference.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
sys.path.insert(0, TESTS)
import run_test  # noqa: E402

CASES_JSON = os.path.join(TESTS, "scf", "gbt_rs_freeze.json")
REFS = os.path.join(HERE, "rs_freeze", "references")

parser = argparse.ArgumentParser()
parser.add_argument("--binary", required=True)
parser.add_argument("--scratch-root", required=True)
parser.add_argument("--gen-ref", action="store_true")
args = parser.parse_args()

names = [c["name"] for c in run_test.json.load(open(CASES_JSON))["cases"]]
if not args.gen_ref and not all(os.path.exists(run_test._reference_path(REFS, n)) for n in names):
    print("SKIP: no committed reference for this platform")
    sys.exit(77)

failed = []
for name in names:
    mode = ["--gen-ref", REFS] if args.gen_ref else ["--compare-ref", REFS, "--abs-tol", "0", "--rel-tol", "0"]
    cmd = [sys.executable, os.path.join(TESTS, "run_test.py"), "--binary", args.binary, "--cases-json", CASES_JSON,
           "--case-name", name, "--scratch-root", args.scratch_root, "--force-serial"] + mode
    if subprocess.run(cmd).returncode:
        failed.append(name)
if args.gen_ref and not failed and run_test._reference_variant():
    for name in names:  # keep one reference per platform, as tests/scf/references does
        os.replace(os.path.join(REFS, name, "ref.json"), os.path.join(REFS, name, f"ref.{run_test._reference_variant()}.json"))
        os.remove(os.path.join(REFS, name, "meta.json"))
print("FAILED: " + ", ".join(failed) if failed else "PASS [AppGbtKspaceRsFreeze]")
sys.exit(1 if failed else 0)
