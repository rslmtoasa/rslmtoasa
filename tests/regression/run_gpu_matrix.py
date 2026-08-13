#!/usr/bin/env python3
"""
GPU-vs-CPU consistency matrix for the CUDA recursion plugin.

Not part of CI (no GPU there) - this is the manual pre-release step described
in tests/README.md, run by hand on a workstation with an NVIDIA GPU and a
build configured with -DENABLE_CUDA_PLUGIN=ON.

For every recur='chebyshev' case in tests/regression/cases.json (gpu_plugin
dispatch bypasses cheb_backend entirely - see recursion_chebyshev.f90's
gpu_plugin_ready() gate - so cheb_backend variants are a CPU-side concern
only), runs the case twice: once with gpu_plugin=false (the case's own
cheb_backend) and once with gpu_plugin=true (CUDA), then diffs
etot/ws_r/vmad between the two runs directly. This does not need or use the
committed CPU references - it is a pairwise CPU/GPU check, not a regression
check against history.

Usage:
  run_gpu_matrix.py --binary build/bin/rslmto.x
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manifest_defaults import apply_manifest_defaults


def load_chebyshev_cases(cases_json: Path) -> list[dict]:
    with cases_json.open() as fh:
        data = json.load(fh)
    cases = [apply_manifest_defaults(data, c) for c in data["cases"]]
    return [c for c in cases if c["namelists"]["control"].get("recur") == "chebyshev"]


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


def run_variant(binary: Path, tests_dir: Path, case: dict, workdir: Path, gpu: bool) -> dict:
    setup_workdir(tests_dir / case["base"], workdir)
    patch = dict(case["namelists"])
    patch["control"] = dict(patch.get("control", {}))
    patch["control"]["gpu_plugin"] = gpu
    tmp_path = workdir / "input.nml.tmp"
    f90nml.patch(str(workdir / "input.nml"), patch, str(tmp_path))
    tmp_path.replace(workdir / "input.nml")

    result = subprocess.run([str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (workdir / "testrun.log").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"binary exited {result.returncode}; see {workdir / 'testrun.log'}")

    output_file = case.get("output_file", "Fe_out.nml")
    out = f90nml.read(str(workdir / output_file))["par"]
    return {"etot": float(out["etot"]), "ws_r": float(out["ws_r"]), "vmad": float(out.get("vmad", 0.0))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, help="Path to a rslmto.x built with -DENABLE_CUDA_PLUGIN=ON")
    parser.add_argument("--cases-json", default=None, help="Defaults to tests/regression/cases.json next to this script")
    parser.add_argument("--scratch-root", default=None, help="Scratch directory (default: temp dir)")
    parser.add_argument("--tol", type=float, default=1e-5, help="Absolute tolerance for CPU-vs-GPU etot/ws_r/vmad")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    cases_json = Path(args.cases_json) if args.cases_json else script_dir / "cases.json"
    binary = Path(args.binary).resolve()

    import tempfile
    use_tmp = args.scratch_root is None
    scratch = Path(args.scratch_root).resolve() if args.scratch_root else Path(tempfile.mkdtemp(prefix="rslmto_gpu_matrix_"))

    if not cmake_option_enabled(binary, "ENABLE_CUDA_PLUGIN"):
        print(f"ERROR: {binary} was not built with -DENABLE_CUDA_PLUGIN=ON (no CMakeCache.txt match).")
        return 1

    cases = load_chebyshev_cases(cases_json)
    rows: list[tuple[str, str, str]] = []
    failures = 0

    for case in cases:
        name = case["name"]
        required = case.get("requires_cmake_option")
        if required and not cmake_option_enabled(binary, required):
            rows.append((name, "SKIP", f"{required}=OFF"))
            continue
        try:
            cpu = run_variant(binary, cases_json.parent, case, scratch / f"{name}_cpu", gpu=False)
            gpu = run_variant(binary, cases_json.parent, case, scratch / f"{name}_gpu", gpu=True)
        except Exception as exc:  # noqa: BLE001 - report and continue the matrix
            rows.append((name, "ERROR", str(exc)))
            failures += 1
            continue

        diffs = {k: abs(cpu[k] - gpu[k]) for k in ("etot", "ws_r", "vmad")}
        worst = max(diffs.values())
        status = "PASS" if worst <= args.tol else "FAIL"
        if status == "FAIL":
            failures += 1
        detail = ", ".join(f"{k}={v:.3e}" for k, v in diffs.items())
        rows.append((name, status, detail))

    print(f"{'Case':<40} {'Status':<7} Detail")
    print("-" * 80)
    for name, status, detail in rows:
        print(f"{name:<40} {status:<7} {detail}")
    print("-" * 80)
    print(f"{len(rows) - failures}/{len(rows)} passed")

    if use_tmp:
        shutil.rmtree(scratch, ignore_errors=True)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
