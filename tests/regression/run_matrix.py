#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manifest_defaults import apply_manifest_defaults


def load_case(cases_json: Path, case_name: str) -> dict:
    with cases_json.open() as fh:
        data = json.load(fh)
    for case in data["cases"]:
        if case["name"] == case_name:
            return apply_manifest_defaults(data, case)
    raise KeyError(f"case {case_name!r} not found in {cases_json}")


def cmake_option_enabled(binary: Path, option: str) -> bool:
    cache = binary.parent.parent / "CMakeCache.txt"
    if not cache.exists():
        return False
    needle = f"{option}:BOOL="
    for line in cache.read_text().splitlines():
        if line.startswith(needle):
            return line.rsplit("=", 1)[-1].upper() in {"ON", "TRUE", "1"}
    return False


def setup_workdir(base_dir: Path, workdir: Path) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base_dir, workdir)
    for pattern in ("*_out.nml", "data.nml", "ref.nml", "testrun.log"):
        for path in workdir.glob(pattern):
            path.unlink()


def patch_input(workdir: Path, patch: dict) -> None:
    input_path = workdir / "input.nml"
    tmp_path = workdir / "input.nml.tmp"
    f90nml.patch(str(input_path), patch, str(tmp_path))
    tmp_path.replace(input_path)


def run_binary(binary: Path, workdir: Path) -> None:
    result = subprocess.run([str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (workdir / "testrun.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(result.returncode)


def compare_nml(run_file: Path, ref_file: Path, tol: float) -> None:
    run = f90nml.read(str(run_file))
    ref = f90nml.read(str(ref_file))
    for key in ("etot", "ws_r", "vmad"):
        run_value = float(run["par"][key])
        ref_value = float(ref["par"][key])
        if abs(run_value - ref_value) > tol:
            raise AssertionError(f"{key}: run={run_value:.12e} ref={ref_value:.12e}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--cases-json", required=True)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--scratch-root", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--gen-ref", action="store_true")
    parser.add_argument("--abs-tol", type=float, default=1e-6)
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    cases_json = Path(args.cases_json).resolve()
    tests_dir = cases_json.parent
    case = load_case(cases_json, args.case_name)

    required = case.get("requires_cmake_option")
    if required and not cmake_option_enabled(binary, required):
        print(f"SKIP [{case['name']}]: {required}=OFF")
        return 0

    workdir = Path(args.scratch_root).resolve() / case["name"]
    setup_workdir(tests_dir / case["base"], workdir)
    patch_input(workdir, case["namelists"])
    run_binary(binary, workdir)

    output = workdir / "Fe_out.nml"
    if not output.exists():
        print((workdir / "testrun.log").read_text()[-3000:])
        raise SystemExit(f"ERROR [{case['name']}]: Fe_out.nml was not produced")

    ref_dir = Path(args.references).resolve()
    ref_path = ref_dir / f"{case['name']}.nml"
    if args.gen_ref:
        ref_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(output, ref_path)
        print(f"WROTE [{case['name']}]: {ref_path}")
        return 0

    if not ref_path.exists():
        raise SystemExit(f"ERROR [{case['name']}]: missing reference {ref_path}")

    compare_nml(output, ref_path, args.abs_tol)
    print(f"PASS [{case['name']}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
