#!/usr/bin/env python3
"""Run the opt-in ACC-P1b reciprocal-SCF observable campaign.

The driver stages the canonical Si and bcc-Fe fixtures into temporary
directories, runs the benchmark-only physical probe, and records observable
comparisons without changing the source inputs or the production executable.
CUDA rows are retained as unsupported when the host has no usable device.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
STRATEGIES = (
    "fp64_zheevd",
    "fp64_zheevj_batched",
    "fp32_cheevd",
    "fp32_cheevj_batched",
)
OBSERVABLES = (
    "electron_count",
    "fermi_level",
    "band_energy",
    "site_occupation",
    "site_charge_transfer",
    "site_moment",
    "dos_state_count",
    "near_ef_min_abs",
    "near_ef_value_1",
    "near_ef_value_2",
    "near_ef_value_3",
    "near_ef_value_4",
    "scf_residual",
    "scf_iterations",
)
FLOAT = r"[-+0-9.EeDd]+"
TOKEN_RE = re.compile(rf"(?P<key>[A-Za-z0-9_]+)=(?P<value>{FLOAT}|true|false|PASS|UNSUPPORTED)")
ITERATION_RE = re.compile(r"(?:Not converged!|Converged!)")


def parse_strategies(text: str) -> list[str]:
    values = [item.strip() for item in text.split(",") if item.strip()]
    if not values or any(value not in STRATEGIES for value in values):
        raise argparse.ArgumentTypeError(f"strategies must be drawn from {list(STRATEGIES)}")
    return values


def parse_probe_output(text: str, returncode: int) -> dict[str, Any]:
    line = next((line for line in text.splitlines() if line.startswith("ACCP1B_SCF ")), None)
    if line is None:
        return {"status": "ERROR", "returncode": returncode, "error_tail": text[-1200:]}
    values: dict[str, Any] = {match.group("key"): match.group("value") for match in TOKEN_RE.finditer(line)}
    for key in list(values):
        if key in OBSERVABLES:
            values[key] = float(str(values[key]).replace("D", "E").replace("d", "e"))
    values["scf_iterations"] = len(ITERATION_RE.findall(text))
    values["returncode"] = returncode
    values["log_tail"] = text[-1200:] if returncode else None
    return values


def run_probe(
    *, binary: Path, fixture: Path, backend: str, strategy: str, sigma: float, temperature: float, nstep: int,
) -> dict[str, Any]:
    command = [
        str(binary),
        "--input", "input.nml",
        "--backend", backend,
        "--solver-strategy", strategy,
        "--dos-method", "gaussian",
        "--sigma", str(sigma),
        "--temperature", str(temperature),
        "--nstep", str(nstep),
    ]
    completed = subprocess.run(command, cwd=fixture, capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    result = parse_probe_output(text, completed.returncode)
    result.update({
        "backend": backend,
        "solver_strategy": strategy,
        "sigma": sigma,
        "temperature": temperature,
        "nstep_requested": nstep,
    })
    return result


def compare_observables(reference: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    comparison: dict[str, Any] = {
        "reference_backend": reference.get("backend"),
        "candidate_backend": candidate.get("backend"),
        "candidate_strategy": candidate.get("solver_strategy"),
        "status": candidate.get("status"),
        "absolute_differences": {},
    }
    if candidate.get("status") != "PASS" or reference.get("status") != "PASS":
        return comparison
    for key in OBSERVABLES:
        left = reference.get(key)
        right = candidate.get(key)
        comparison["absolute_differences"][key] = None if left is None or right is None else abs(float(right) - float(left))
    if reference.get("converged") != "true" or candidate.get("converged") != "true":
        comparison["status"] = "NOT_CONVERGED"
    return comparison


def classify_scf(rows: list[dict[str, Any]], system: str) -> str:
    candidates = [
        row for row in rows
        if row.get("system") == system
        and row.get("backend") == "cuda"
        and str(row.get("solver_strategy", "")).startswith("fp32_")
    ]
    if not candidates:
        return "not_run"
    supported = [row for row in candidates if row.get("status") == "PASS"]
    if not supported:
        # A missing CUDA device or an explicit unsupported solver route is an
        # unavailable gate, not a failed scientific gate.  Keep FP32 out of
        # production optimization until a supported run exists.
        return "unsupported"
    if any(row.get("converged") != "true" for row in supported):
        return "unacceptable"
    if len(supported) != len(candidates):
        return "partially_evaluated"
    return "acceptable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nstep", type=int, default=12)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--strategies", type=parse_strategies, default=list(STRATEGIES))
    parser.add_argument("--skip-cuda", action="store_true")
    args = parser.parse_args()

    sources = {
        "Si": ROOT / "tests/scf/cases/bulk/diamondSi",
        "bccFe": ROOT / "tests/scf/cases/bulk/bccFe",
    }
    rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    gpu_comparisons: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="rslmto-accp1b-scf-", dir="/tmp") as temporary:
        scratch = Path(temporary)
        for system, source in sources.items():
            sigmas = [0.01] if system == "Si" else [0.005, 0.01, 0.02]
            for sigma in sigmas:
                configurations = [("lapack", "fp64_zheevd")]
                if not args.skip_cuda:
                    configurations.extend(("cuda", strategy) for strategy in args.strategies)
                for backend, strategy in configurations:
                    fixture = scratch / f"{system}_{sigma:g}_{backend}_{strategy}"
                    shutil.copytree(source, fixture)
                    row = run_probe(
                        binary=args.binary,
                        fixture=fixture,
                        backend=backend,
                        strategy=strategy,
                        sigma=sigma,
                        temperature=args.temperature,
                        nstep=args.nstep,
                    )
                    row.update({"system": system, "sigma": sigma})
                    rows.append(row)
                reference = next(
                    row for row in rows
                    if row["system"] == system and row["sigma"] == sigma and row["backend"] == "lapack"
                )
                for row in rows:
                    if row["system"] == system and row["sigma"] == sigma and row["backend"] == "cuda":
                        comparisons.append({"system": system, "sigma": sigma, **compare_observables(reference, row)})
                gpu_reference = next(
                    (
                        row for row in rows
                        if row["system"] == system and row["sigma"] == sigma
                        and row["backend"] == "cuda" and row["solver_strategy"] == "fp64_zheevd"
                    ),
                    None,
                )
                if gpu_reference is not None:
                    for row in rows:
                        if (
                            row["system"] == system and row["sigma"] == sigma
                            and row["backend"] == "cuda" and str(row["solver_strategy"]).startswith("fp32_")
                        ):
                            gpu_comparisons.append({
                                "system": system,
                                "sigma": sigma,
                                **compare_observables(gpu_reference, row),
                            })

    report = {
        "schema": "rslmto.accp1b-physical-scf.v1",
        "policy": {
            "production_executable_unchanged": True,
            "production_precision_unchanged": True,
            "fp32_is_default": False,
            "normal_physics_arrays": "complex(rp)/real(rp)",
            "delicate_small_energy_physics": "not_run",
            "mixed_scf_probe": "not_run",
        },
        "rows": rows,
        "comparisons_to_fp64_cpu": comparisons,
        "comparisons_to_fp64_gpu": gpu_comparisons,
        "classification": {
            "eigensystem_benchmark": "reported_by_accp0_validation",
            "Si_SCF": classify_scf(rows, "Si"),
            "metallic_Fe_SCF": classify_scf(rows, "bccFe"),
            "physical_acceptance_basis": "all supported requested FP32 CUDA rows must report converged=true; unsupported is not a pass",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
