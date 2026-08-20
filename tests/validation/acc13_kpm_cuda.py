#!/usr/bin/env python3
"""ACC-13 production RS transport CPU/CUDA probe.

This is a focused hardware campaign, separate from the broad VAL-09 physics
campaign.  It runs the same production fcc-Pt SOC fixture once through the
legacy CPU recursion route and once through the CUDA plugin for each of the
charge, spin, and orbital conductivity contracts.  The GPU route is the
existing production moment path; this driver does not reimplement transport.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from val09_kubo_bastin_transport import parse_conductivity, patch_case


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "tests/postproc/cases/conductivity/fccPt"
RUNNER = Path(__file__).resolve().parents[1] / "run_binary.sh"


def run_case(binary: Path, scratch: Path, cond_type: str, gpu: bool,
             cond_ll: int, channels: int, replication: int, rc: int,
             gpu_precision: str, lld: int, cpu_cheb_backend: str) -> dict[str, float | bool]:
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    patch_case(
        BASE,
        scratch,
        cond_type=cond_type,
        va=[1, 0, 0] if cond_type == "charge" else [0, 1, 0],
        vb=[1, 0, 0] if cond_type == "charge" else [1, 0, 0],
        cond_ll=cond_ll,
        channels=channels,
        replication=replication,
        rc=rc,
        gpu_plugin=gpu,
        gpu_backend="csr",
        gpu_precision=gpu_precision,
        cheb_backend=cpu_cheb_backend,
        lld=lld,
    )
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    start = time.perf_counter()
    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=scratch,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - start
    log_path = scratch / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"ACC-13 {cond_type} {'GPU' if gpu else 'CPU'} failed:\n{log[-6000:]}")
    if gpu and "CUDA_DEVICE_MAPPING" not in log:
        raise RuntimeError(f"ACC-13 {cond_type} did not enter the CUDA plugin route")
    observable = parse_conductivity(scratch)
    return {**observable, "elapsed_seconds": elapsed, "gpu": gpu}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path,
                        help="CUDA-plugin build of rslmto.x")
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-3)
    parser.add_argument("--cond-ll", type=int, default=20)
    parser.add_argument("--lld", type=int, default=50)
    parser.add_argument("--channels", type=int, default=120)
    parser.add_argument("--replication", type=int, default=4)
    parser.add_argument("--rc", type=int, default=20)
    parser.add_argument("--gpu-precision", choices=("fp32", "fp64"), default="fp32")
    parser.add_argument("--cpu-cheb-backend", choices=("legacy", "fast", "fast_dp"),
                        default="legacy",
                        help="CPU reference precision route; use fast for FP32 matching")
    parser.add_argument("--cond-types", nargs="+", choices=("charge", "spin", "orbital"),
                        default=("charge", "spin", "orbital"),
                        help="conductivity observables to compare")
    args = parser.parse_args()

    binary = args.binary.resolve()
    scratch = args.scratch_root.resolve()
    results: dict[str, dict[str, dict[str, float | bool]]] = {}
    failures: list[str] = []
    for cond_type in args.cond_types:
        cpu = run_case(binary, scratch / f"{cond_type}_cpu", cond_type, False,
                       args.cond_ll, args.channels, args.replication, args.rc,
                       args.gpu_precision, args.lld, args.cpu_cheb_backend)
        gpu = run_case(binary, scratch / f"{cond_type}_gpu", cond_type, True,
                       args.cond_ll, args.channels, args.replication, args.rc,
                       args.gpu_precision, args.lld, args.cpu_cheb_backend)
        scale = max(abs(float(cpu["real"])), abs(float(gpu["real"])),
                    abs(float(cpu["imag"])), abs(float(gpu["imag"])), 1.0e-12)
        relative_error = (
            (float(cpu["real"]) - float(gpu["real"])) ** 2
            + (float(cpu["imag"]) - float(gpu["imag"])) ** 2
        ) ** 0.5 / scale
        speedup = float(cpu["elapsed_seconds"]) / max(float(gpu["elapsed_seconds"]), 1.0e-12)
        results[cond_type] = {
            "cpu": cpu,
            "gpu": gpu,
            "relative_complex_error": relative_error,
            "cpu_over_gpu_wall_speedup": speedup,
        }
        if relative_error > args.relative_tolerance:
            failures.append(
                f"{cond_type}: relative real-part error {relative_error:.3e} "
                f"> {args.relative_tolerance:.3e}"
            )

    report = {
        "scope": (f"fcc Pt SOC production conductivity fixture; cond_ll={args.cond_ll}, "
                  f"lld={args.lld}, channels={args.channels}, replication={args.replication}, rc={args.rc}"),
        "gpu_precision": f"CUDA {args.gpu_precision} moments, fp64 host outputs",
        "cpu_cheb_backend": args.cpu_cheb_backend,
        "cond_types": args.cond_types,
        "relative_tolerance": args.relative_tolerance,
        "results": results,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if failures:
        print("ACC-13 FAIL")
        print("\n".join(failures))
        return 1
    print("ACC-13 PASS: production charge/spin/orbital CUDA transport contracts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
