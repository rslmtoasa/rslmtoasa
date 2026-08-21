#!/usr/bin/env python3
"""Fast B0C contract tests; no production benchmark execution."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_harness import (  # noqa: E402
    build_comparison_key,
    capture_environment,
    empty_correctness,
    numeric_mode,
    validate_pairing,
)
from kpm_g12_transport import CSV_COLUMNS, add_speedups, write_outputs  # noqa: E402


def row(*, mode: str = "fp64", material: str = "fccPt_SOC", M: int = 500,
        lld: int = 150, NE: int = 2510, seed: str = "per_type_projectors",
        gpu: bool = False, omp: int = 1, transport: float = 10.0) -> dict:
    moment, reconstruction = {
        "fp64": ("fp64", "fp64"), "fp32": ("fp32", "fp32"), "mixed": ("fp32", "fp64"),
    }[mode]
    result = {
        "row_id": f"{'gpu' if gpu else 'cpu'}-{mode}-{omp}", "material": material,
        "fixture_id": material + "_fixture", "fixture_revision": "fixture-v1",
        "replication": 4, "N": 1152, "nnz": 10000, "nsp": 2, "soc_state": "SOC",
        "cond_type": "spin", "current_operator": {"va": [0, 1, 0], "vb": [1, 0, 0], "spin_orbital_selector": "S_z"},
        "cond_calctype": "random_vec" if seed != "per_type_projectors" else "per_type",
        "Ntrace": 16 if seed != "per_type_projectors" else 1, "projector_count": 1,
        "random_seed": seed, "vector_contract": "uniform_unit_phase_normalized_by_sqrt_sites",
        "M": M, "lld": lld, "NE": NE, "channels": NE - 10, "rc": 20,
        "energy_min": -2.5, "energy_max": 1.2, "fermi": -0.085837,
        "kernel": {"name": "Lorentz", "alpha": 6.0},
        "chebyshev_scaling": {"contract": "(E-b)/a", "window": [-2.5, 1.2]},
        "moment_backend": "cuda" if gpu else "cpu_fast", "moment_precision": moment,
        "reconstruction_backend": "cuda_blas" if gpu else "cpu_cblas", "reconstruction_precision": reconstruction,
        "canonical_output_precision": "fp64",
        "numeric_mode": mode, "gpu_plugin": gpu, "OMP_threads": omp, "BLAS_threads": 1,
        "profile_status": "PASS", "output_mode": "production", "correctness": {
            **empty_correctness("test"), "status": "PASS", "tolerance_set": {"name": mode},
        },
        "statistics": {
            "P_moments_total": {"median": transport / 2}, "T_transport_total": {"median": transport},
            "wall_time_s": {"median": transport + 1},
        },
    }
    result["comparison_key"] = build_comparison_key(result)
    return result


def test_numeric_taxonomy() -> None:
    assert numeric_mode("fp64", "fp64") == "fp64"
    assert numeric_mode("fp32", "fp32") == "fp32"
    assert numeric_mode("fp32", "fp64") == "mixed"


def test_pairing_rejects_precision_physics_seed_correctness_profile_and_no_write() -> None:
    cpu = row(mode="fp64")
    gpu = row(mode="fp64", gpu=True)
    assert validate_pairing(cpu, gpu)["eligible"]
    assert "precision_mismatch" in validate_pairing(cpu, row(mode="fp32", gpu=True))["reasons"]
    assert "physics_mismatch" in validate_pairing(cpu, row(M=501, gpu=True))["reasons"]
    assert "physics_mismatch" in validate_pairing(cpu, row(lld=149, gpu=True))["reasons"]
    assert "physics_mismatch" in validate_pairing(cpu, row(NE=2000, gpu=True))["reasons"]
    assert "physics_mismatch" in validate_pairing(cpu, row(material="bccFe_magnetic", gpu=True))["reasons"]
    assert "seed_mismatch" in validate_pairing(row(seed="11"), row(seed="12", gpu=True))["reasons"]
    failed = row()
    failed["correctness"] = {**empty_correctness("test"), "status": "FAIL"}
    assert "correctness_failed" in validate_pairing(failed, gpu)["reasons"]
    profile_failed = row()
    profile_failed["profile_status"] = "FAIL"
    assert "profile_failed" in validate_pairing(profile_failed, gpu)["reasons"]
    no_write = row()
    no_write["output_mode"] = "benchmark_no_write"
    assert "benchmark_no_write" in validate_pairing(no_write, gpu)["reasons"]
    assert "mixed_precision_not_headline" in validate_pairing(row(mode="mixed"), row(mode="mixed", gpu=True))["reasons"]


def test_best_cpu_uses_only_valid_omp_sweep_rows() -> None:
    rows = [row(omp=1, transport=8.0), row(omp=2, transport=5.0), row(omp=4, transport=3.0), row(omp=8, transport=4.0), row(gpu=True, omp=1, transport=2.0)]
    add_speedups(rows)
    gpu = rows[-1]
    assert gpu["headline_speedup_eligible"]
    assert gpu["best_cpu_row_id"] == "cpu-fp64-4"
    assert gpu["speedups"]["S_transport"] == 1.5


def test_outputs_and_environment_contract() -> None:
    report = {"environment": capture_environment(ROOT), "rows": [row(), row(mode="mixed")], "schema": "rslmto.kpm-b0c.v1"}
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "campaign.json"
        write_outputs(report, output)
        with output.with_suffix(".csv").open(newline="", encoding="utf-8") as handle:
            assert set(CSV_COLUMNS).issubset(csv.DictReader(handle).fieldnames or [])
        markdown = output.with_suffix(".md").read_text(encoding="utf-8")
        for section in ("Environment", "Equal-FP64 pairs", "Equal-FP32 pairs", "Mixed practical rows", "Correctness summary", "Definitions"):
            assert section in markdown
    environment = capture_environment(ROOT)
    for key in ("git_commit", "git_dirty", "compiler", "blas_lapack", "cpu_model", "ram_mib", "cuda_toolkit", "gpu_model", "selected_gpu_index"):
        assert key in environment


def test_material_contract() -> None:
    from kpm_g12_transport import MATERIALS
    assert MATERIALS["pt"].fixture_id == "fccPt_SOC_val09"
    assert MATERIALS["fe"].fixture_id == "bccFe_magnetic_triad"


if __name__ == "__main__":
    for function in (test_numeric_taxonomy, test_pairing_rejects_precision_physics_seed_correctness_profile_and_no_write,
                     test_best_cpu_uses_only_valid_omp_sweep_rows, test_outputs_and_environment_contract, test_material_contract):
        function()
    print("PASS: KPM-B0C harness contract")
