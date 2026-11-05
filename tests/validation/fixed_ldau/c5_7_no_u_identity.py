#!/usr/bin/env python3
"""Single-thread no-U RS/KS identity gate for C5.7."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PREFLIGHT = ROOT / "docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py"
LAUNCHER = ROOT / "tests/validation/fixed_ldau/c5_7_launch.sh"
BASE = Path("/private/tmp/c52-rd5-no-u")
ARCHIVED_INPUTS = ROOT / "docs/validation/fixed_ldau/rs/C5_2_PART_D_assets/rd5_no_u_identity_checks"
INPUTS = {
    # Use the repository-archived copies of the exact original identity inputs.
    # Their hashes match the C5.7 successful run; the former /private/tmp source
    # paths are disposable and may exist as empty directories after cleanup.
    "rs": ARCHIVED_INPUTS / "rs_rd1",
    "ks": ARCHIVED_INPUTS / "ks_rd1",
}
CAPS = {
    "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "RSLMTO_OMP_THREADS_SERIAL": "1",
    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized(path: Path) -> bytes:
    data = path.read_bytes()
    if path.name == "str.out":
        data = re.sub(rb"(?m)(cpu_s=\s*)[0-9.]+", rb"\1<RUNTIME>", data)
    return data


def run_one(binary: Path, scratch: Path, route: str, baseline_override: Path | None = None) -> dict:
    baseline = baseline_override if baseline_override is not None else BASE / f"{route}_rd1"
    work = scratch / f"{route}_candidate"
    work.mkdir()
    input_source = INPUTS[route]
    for name in ("input.nml", "Fe.nml"):
        shutil.copy2(input_source / name, work / name)
        target = work / name
        target.chmod(target.stat().st_mode | 0o200)
    env = os.environ.copy()
    env.update(CAPS)
    pre = subprocess.run(
        [sys.executable, str(PREFLIGHT), str(work), "--deck", str(work / "input.nml"),
         "--restart", str(work / "Fe.nml"), "--guard", str(LAUNCHER)],
        cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if pre.returncode:
        raise RuntimeError(f"{route} preflight failed: {pre.stdout[-4000:]}")
    (work / "preflight.json").write_text(pre.stdout)
    df_pre = subprocess.check_output(["df", "-k", "/private/tmp"], text=True)
    if shutil.disk_usage("/private/tmp").free < 1_073_741_824:
        raise RuntimeError("1 GiB /private/tmp floor reached before run")
    proc = subprocess.Popen(
        ["/bin/bash", str(LAUNCHER), str(binary)], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True,
    )
    samples = [shutil.disk_usage("/private/tmp").free]
    while proc.poll() is None:
        time.sleep(1)
        samples.append(shutil.disk_usage("/private/tmp").free)
        if samples[-1] < 1_073_741_824:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            raise RuntimeError("1 GiB /private/tmp live guard stopped run")
    if proc.returncode:
        raise RuntimeError(f"{route} no-U candidate failed: {(work/'testrun.log').read_text()[-4000:]}")
    exclude = {"testrun.log", "preflight.json"}
    base_files = {p.name: p for p in baseline.iterdir() if p.is_file() and p.name not in exclude}
    cand_files = {p.name: p for p in work.iterdir() if p.is_file() and p.name not in exclude}
    mismatches = []
    for name in sorted(set(base_files) | set(cand_files)):
        if name not in base_files or name not in cand_files:
            mismatches.append({"name": name, "reason": "name_set"})
        elif normalized(base_files[name]) != normalized(cand_files[name]):
            mismatches.append({"name": name, "reason": "bytes",
                               "baseline_sha256": sha(base_files[name]), "candidate_sha256": sha(cand_files[name])})
    if mismatches:
        raise RuntimeError(f"{route} no-U identity mismatches: {mismatches[:5]}")
    return {
        "status": "PASS", "baseline": str(baseline), "input_source": str(input_source), "candidate": str(work),
        "input_sha256": sha(work / "input.nml"), "restart_sha256": sha(work / "Fe.nml"),
        "artifact_count": len(cand_files), "mismatches": mismatches,
        "preflight_path": str(work / "preflight.json"), "thread_caps": CAPS,
        "df_pre": df_pre, "df_post": subprocess.check_output(["df", "-k", "/private/tmp"], text=True),
        "solver_wrapper_pid": proc.pid, "returncode": proc.returncode,
        "disk_free_bytes_samples": samples, "disk_free_bytes_after": shutil.disk_usage("/private/tmp").free,
        "artifact_sha256": {name: sha(path) for name, path in sorted(cand_files.items())},
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--rs-baseline", type=Path,
                        help="override the historical RS baseline with a versioned artifact directory")
    parser.add_argument("--ks-baseline", type=Path,
                        help="override the historical KS baseline with a versioned artifact directory")
    parser.add_argument("--result", type=Path,
                        default=ROOT / "docs/validation/fixed_ldau/rs/C5_7_NO_U_IDENTITY.json",
                        help="identity result JSON destination")
    args = parser.parse_args()
    binary = args.binary.resolve()
    args.scratch.mkdir(parents=True, exist_ok=False)
    baselines = {"rs": args.rs_baseline, "ks": args.ks_baseline}
    results = {route.upper(): run_one(binary, args.scratch, route, baselines[route]) for route in ("rs", "ks")}
    result = {
        "schema": "c5_7_no_u_identity_v1", "binary": str(binary), "binary_sha256": sha(binary),
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": sha(ROOT / "source/self.f90"),
        "runner_sha256": sha(Path(__file__).resolve()),
        "criterion": "exact artifact names/bytes; exclude testrun.log and preflight record; normalize only cpu_s fields in str.out",
        "results": results,
    }
    destination = args.result
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({route: {"status": row["status"], "artifact_count": row["artifact_count"]}
                      for route, row in results.items()}, indent=2))


if __name__ == "__main__":
    main()
