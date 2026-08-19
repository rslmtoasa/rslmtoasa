#!/usr/bin/env python3
"""Compare validated reciprocal HOH/CCOR physical outputs on CPU and CUDA.

This is an opt-in campaign, not a timing gate.  It uses the existing compact
post-processing fixture and compares downstream DOS data plus canonical
occupation/band observables.  The timing records are retained as evidence,
but no speedup threshold is imposed.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FLOAT = r"[-+0-9.EeDd]+"
CANONICAL_RE = re.compile(
    rf"Canonical k-space occupations: EF=\s*({FLOAT}).*?N=\s*({FLOAT}),"
    rf"\s*dN=\s*({FLOAT}),.*?EBAND=\s*({FLOAT})"
)
TIMING_RE = re.compile(
    rf"ACC05_TIMING backend=(\w+) total_seconds=\s*({FLOAT})"
    rf"\s+host_assembly_seconds=\s*({FLOAT})"
)
CUDA_TIMING_RE = re.compile(
    rf"ACC05_TIMING backend=cuda total_seconds=\s*({FLOAT})"
    rf"\s+host_assembly_seconds=\s*({FLOAT})"
    rf"\s+h2d_seconds=\s*({FLOAT})"
    rf"\s+gpu_solve_seconds=\s*({FLOAT})"
    rf"\s+d2h_seconds=\s*({FLOAT})"
)

ROOT = Path(__file__).resolve().parents[2]
CASE_ROOT = ROOT / "tests/postproc/cases/density_of_states/bccFe"
RUN_BINARY = ROOT / "tests/run_binary.sh"


def as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def variant_patches(backend: str) -> dict[str, dict]:
    common = {
        "self": {"nstep": 1},
        "reciprocal": {
            "nk1": 2,
            "nk2": 2,
            "nk3": 2,
            "n_energy_points": 40,
            "reciprocal_backend": backend,
        },
        "control": {"nsp": 1, "recur": "block", "lld": 12},
    }
    return {
        "hoh": {**common, "hamiltonian": {"hoh": True}},
        "ccor": {
            **common,
            "lattice": {
                "strux_backend": "strux_lib",
                "screening": "fitted",
                "strux_want_sdot": True,
            },
            "hamiltonian": {"ccor_2c": True},
        },
    }


def parse_dos(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("#", 1)[0].split()
        if fields:
            rows.append([float(value.replace("D", "E").replace("d", "e")) for value in fields])
    if not rows:
        raise RuntimeError(f"no numeric DOS rows in {path}")
    return rows


def run_variant(binary: Path, backend: str, variant: str, workdir: Path) -> dict:
    if workdir.exists():
        shutil.rmtree(workdir)
    # This fixture intentionally preserves its converged Fe_out.nml restart,
    # matching the existing post-processing reference cases.
    shutil.copytree(CASE_ROOT, workdir)
    patch_path = workdir / "input.nml.tmp"
    f90nml.patch(
        str(workdir / "input.nml"),
        variant_patches(backend)[variant],
        str(patch_path),
    )
    patch_path.replace(workdir / "input.nml")
    completed = subprocess.run(
        ["bash", str(RUN_BINARY), str(binary.resolve())],
        cwd=workdir,
        env={**os.environ, "RSLMTO_OMP_THREADS_SERIAL": "1"},
        text=True,
        capture_output=True,
        check=False,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else completed.stderr
    if completed.returncode != 0:
        raise RuntimeError(f"{variant}/{backend} failed:\n{log[-4000:]}")
    canonical_matches = list(CANONICAL_RE.finditer(log))
    timing_matches = list(TIMING_RE.finditer(log))
    if not canonical_matches or not timing_matches:
        raise RuntimeError(f"{variant}/{backend} omitted canonical or ACC05 timing diagnostics")
    canonical = canonical_matches[-1]
    timing = timing_matches[-1]
    cuda_timing = CUDA_TIMING_RE.search(log) if backend == "cuda" else None
    timing_values = {
        "total_seconds": as_float(timing.group(2)),
        "host_assembly_seconds": as_float(timing.group(3)),
        "h2d_seconds": 0.0,
        "gpu_solve_seconds": 0.0,
        "d2h_seconds": 0.0,
    }
    if backend == "cuda":
        if cuda_timing is None:
            raise RuntimeError(f"{variant}/{backend} omitted CUDA timing fields")
        timing_values.update(
            h2d_seconds=as_float(cuda_timing.group(3)),
            gpu_solve_seconds=as_float(cuda_timing.group(4)),
            d2h_seconds=as_float(cuda_timing.group(5)),
        )
    return {
        "variant": variant,
        "backend": backend,
        "canonical": {
            "ef": as_float(canonical.group(1)),
            "electron_count": as_float(canonical.group(2)),
            "electron_residual": as_float(canonical.group(3)),
            "band_energy": as_float(canonical.group(4)),
        },
        "timing": {
            "backend": timing.group(1),
            **timing_values,
        },
        "dos": parse_dos(workdir / "dos_kspace.dat"),
    }


def compare(cpu: dict, cuda: dict) -> dict:
    if len(cpu["dos"]) != len(cuda["dos"]):
        raise RuntimeError(f"{cpu['variant']}: CPU/CUDA DOS row counts differ")
    max_dos_error = max(
        abs(left - right)
        for cpu_row, cuda_row in zip(cpu["dos"], cuda["dos"])
        for left, right in zip(cpu_row, cuda_row)
    )
    if max_dos_error > 1.1e-5:
        raise RuntimeError(f"{cpu['variant']}: CPU/CUDA DOS max error {max_dos_error:.3e} exceeds print precision")
    observable_errors = {
        key: abs(cpu["canonical"][key] - cuda["canonical"][key])
        for key in ("ef", "electron_count", "band_energy")
    }
    if max(observable_errors.values()) > 1.0e-6:
        raise RuntimeError(f"{cpu['variant']}: CPU/CUDA canonical observable mismatch {observable_errors}")
    return {"max_dos_abs_error": max_dos_error, "canonical_abs_errors": observable_errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cuda-binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    report = {"variants": {}, "binary": {"cpu": str(args.cpu_binary), "cuda": str(args.cuda_binary)}}
    for variant in ("hoh", "ccor"):
        cpu = run_variant(args.cpu_binary, "lapack", variant, args.scratch_root / f"{variant}-cpu")
        cuda = run_variant(args.cuda_binary, "cuda", variant, args.scratch_root / f"{variant}-cuda")
        report["variants"][variant] = {
            "comparison": compare(cpu, cuda),
            "cpu": {"canonical": cpu["canonical"], "timing": cpu["timing"]},
            "cuda": {"canonical": cuda["canonical"], "timing": cuda["timing"]},
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"ACC-09 reciprocal variants PASS: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
