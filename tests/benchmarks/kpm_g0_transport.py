#!/usr/bin/env python3
"""KPM-G0 real-material transport baseline recorder.

This driver stages the existing Pt and bcc-Fe conductivity fixtures, runs the
production recursion route, and records the ``KPM_PROFILE`` line emitted by
the executable.  It deliberately does not alter physics inputs beyond the
documented benchmark controls.  Repeated persistent timing belongs to the
later B0 harness; G0's purpose is to establish stage ownership and dimensions.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from benchmark_harness import parse_profile_output

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from val09_kubo_bastin_transport import patch_case


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tests/run_binary.sh"
PT_BASE = ROOT / "tests/postproc/cases/conductivity/fccPt"
FE_BASE = ROOT / "tests/regression/triad_bccFe_conductivity"


def run_case(binary: Path, scratch: Path, *, material: str, base: Path,
             cond_type: str, va: list[int], vb: list[int], replication: int,
             cond_ll: int, lld: int, channels: int, rc: int,
             gpu_plugin: bool, omp_threads: int, cheb_backend: str) -> dict[str, Any]:
    if scratch.exists():
        shutil.rmtree(scratch)
    patch_case(
        base,
        scratch,
        cond_type=cond_type,
        va=va,
        vb=vb,
        cond_ll=cond_ll,
        lld=lld,
        channels=channels,
        replication=replication,
        rc=rc,
        gpu_plugin=gpu_plugin,
        gpu_backend="csr",
        cheb_backend=cheb_backend,
    )
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["OMP_PROC_BIND"] = env.get("OMP_PROC_BIND", "close")
    env["OMP_PLACES"] = env.get("OMP_PLACES", "cores")
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    started = time.perf_counter()
    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=scratch,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    wall_time_s = time.perf_counter() - started
    log_path = scratch / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(
            f"KPM-G0 {material} {'GPU' if gpu_plugin else 'CPU'} failed:\n{log[-8000:]}"
        )
    profiles = [record for record in parse_profile_output(log)
                if record.get("name") == "kpm_transport"]
    if len(profiles) != 1:
        raise RuntimeError(
            f"KPM-G0 {material} did not emit exactly one KPM_PROFILE record; "
            f"found {len(profiles)}"
        )
    return {
        "material": material,
        "backend_requested": "cuda" if gpu_plugin else "cpu",
        "replication": replication,
        "cond_ll": cond_ll,
        "lld": lld,
        "channels_ldos": channels,
        "rc": rc,
        "omp_threads": omp_threads,
        "cheb_backend": cheb_backend,
        "wall_time_s": wall_time_s,
        "profile": profiles[0],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cond-ll", type=int, default=500)
    parser.add_argument("--lld", type=int, default=150)
    parser.add_argument("--channels", type=int, default=2500)
    parser.add_argument("--rc", type=int, default=20)
    parser.add_argument("--replications", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--cheb-backend", choices=("legacy", "fast", "fast_dp"),
                        default="legacy",
                        help="CPU moment route; CUDA bypasses this setting")
    parser.add_argument("--skip-fe", action="store_true")
    parser.add_argument("--skip-gpu", action="store_true",
                        help="record CPU baselines only")
    args = parser.parse_args()

    binary = args.binary.resolve()
    scratch_root = args.scratch_root.resolve()
    rows: list[dict[str, Any]] = []
    routes = [False] if args.skip_gpu else [False, True]
    for gpu_plugin in routes:
        for replication in args.replications:
            rows.append(run_case(
                binary, scratch_root / f"pt_r{replication}_{'gpu' if gpu_plugin else 'cpu'}",
                material="fccPt_SOC", base=PT_BASE, cond_type="spin",
                va=[0, 1, 0], vb=[1, 0, 0], replication=replication,
                cond_ll=args.cond_ll, lld=args.lld, channels=args.channels,
                rc=args.rc, gpu_plugin=gpu_plugin, omp_threads=args.omp_threads,
                cheb_backend=args.cheb_backend,
            ))
        if not args.skip_fe:
            rows.append(run_case(
                binary, scratch_root / f"bccFe_{'gpu' if gpu_plugin else 'cpu'}",
                material="bccFe_magnetic", base=FE_BASE, cond_type="charge",
                va=[1, 0, 0], vb=[1, 0, 0], replication=args.replications[0],
                cond_ll=args.cond_ll, lld=args.lld, channels=args.channels,
                rc=max(args.rc, 80), gpu_plugin=gpu_plugin,
                omp_threads=args.omp_threads,
                cheb_backend=args.cheb_backend,
            ))

    report = {
        "schema": "rslmto.kpm-g0.v1",
        "scope": "real-material KPM/Kubo-Bastin stage profile",
        "policy": {
            "cond_ll": args.cond_ll,
            "lld": args.lld,
            "estimator": "per_type",
            "persistent_process": False,
            "note": "G0 baseline; use B0 for persistent repeated timing",
        },
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
