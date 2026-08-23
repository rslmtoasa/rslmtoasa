#!/usr/bin/env python3
"""Focused SCF-B0C contract tests; no real SCF or CUDA run is required."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scf_b0c import (  # noqa: E402
    CSV_FIELDS,
    PROFILE_PHASES,
    SCHEMA,
    build_comparison_key,
    numeric_mode,
    parse_probe_output,
    stage_material_fixture,
    validate_pairing,
    validate_scf_profile,
    write_outputs,
)


def _row(*, material: str = "diamondSi", nk: int = 8, soc: str = "off", mode: str = "fp64",
         start: str = "normal_initial", backend: str = "lapack", omp: int = 1) -> dict:
    precisions = ("fp64", "fp64", "fp64", "fp64", "fp64") if mode == "fp64" else ("fp64", "fp32", "fp32", "fp64", "fp64")
    row = {
        "material": material, "fixture_id": material + "_fixture", "fixture_revision": "v1", "supercell": "1x1x1",
        "Natom": 2, "nmat": 16, "nsp": 1, "SOC": soc, "basis": "sp", "lmax": 1,
        "strux_backend": "strux_lib", "Nk_nominal": "2x2x2", "Nk_unique": nk, "k_mesh": "2x2x2",
        "smearing": {"method": "gaussian", "sigma": 0.01}, "temperature": 300.0, "electron_count": 8.0,
        "fermi_policy": "auto_from_eigenvalues", "mixing_method": "broyden", "mixing_parameters": {"beta": 0.1},
        "convergence_threshold": 5e-10, "starting_state": start, "potential_identity": "fixture",
        "feature_flags": {}, "eigenvectors": "required", "hamiltonian_precision": precisions[0],
        "eigensolver_precision": precisions[1], "eigenvector_precision": precisions[2],
        "density_accumulation_precision": precisions[3], "scf_canonical_precision": precisions[4],
        "numeric_mode": mode, "backend": backend, "fallback_detected": False,
        "profile_status": "PASS", "correctness": {"status": "PASS"}, "starting_state_id": start,
        "solver_strategy": "fp64_zheevd", "OMP_threads": omp, "BLAS_threads": 1,
        "row_id": f"{backend}-{mode}-{omp}", "benchmark_level": "scf_convergence",
        "P_hk_assembly": 1.0, "P_eigensolver": 2.0, "P_occupations_fermi": 1.0,
        "P_density_build": 1.0, "P_potential_update": 1.0, "P_mixing": 1.0, "P_scf_io": 1.0,
        "P_scf_misc": 0.0, "P_scf_iteration_total": 8.0, "steady_iteration_median": 8.0,
        "full_scf_wall": 8.0, "n_scf_iterations": 2,
    }
    row["comparison_key"] = build_comparison_key(row)
    return row


def test_numeric_modes_are_end_to_end() -> None:
    assert numeric_mode("fp64", "fp64", "fp64", "fp64", "fp64") == "fp64"
    assert numeric_mode("fp32", "fp32", "fp32", "fp32", "fp32") == "fp32"
    assert numeric_mode("fp64", "fp32", "fp32", "fp64", "fp64") == "mixed"


def test_profile_closure_and_misc_gate() -> None:
    values = {name: 0.0 for name in PROFILE_PHASES}
    values.update({"P_hamiltonian_prepare": 1.0, "P_eigensolver": 2.0, "P_rs_solver_kernel": 3.0,
                   "P_rs_density_build": 1.0, "P_potential_update": 1.0, "P_mixing": 1.0, "P_scf_io": 1.0})
    good = [{**values, "P_scf_iteration_total": 10.0}]
    assert validate_scf_profile(good)["status"] == "PASS"
    bad = [dict(good[0], P_scf_misc=1.0, P_scf_iteration_total=10.0)]
    assert validate_scf_profile(bad)["status"] == "FAIL"


def test_pairing_contract_rejects_physics_start_precision_and_fallback() -> None:
    cpu = _row()
    gpu = _row(backend="cuda")
    assert validate_pairing(cpu, gpu)["eligible"]
    assert "physics_mismatch" in validate_pairing(cpu, _row(backend="cuda", nk=9))["reasons"]
    assert "physics_mismatch" in validate_pairing(cpu, _row(backend="cuda", soc="on"))["reasons"]
    assert "starting_state_mismatch" in validate_pairing(cpu, _row(backend="cuda", start="restart"))["reasons"]
    assert "numeric_mode_mismatch" in validate_pairing(cpu, _row(backend="cuda", mode="mixed"))["reasons"]
    fallback = _row(backend="cuda")
    fallback["fallback_detected"] = True
    assert "silent_fallback_detected" in validate_pairing(cpu, fallback)["reasons"]


def test_real_space_pairing_and_route_parser() -> None:
    cpu = _row()
    gpu = _row(backend="cuda")
    for row in (cpu, gpu):
        row.update({
            "scf_route": "real_space", "Nk_nominal": None, "Nk_unique": None, "k_mesh": None,
            "rs_solver": "block", "rs_backend": "csr", "recursion_depth": 21, "block_size": 18,
            "terminator": 5, "chebyshev_order": None, "chebyshev_kernel": None,
            "spectral_bounds_policy": None, "rs_kernel_correctness_status": "PASS",
            "P_rs_hamiltonian_prepare": 1.0, "P_rs_solver_kernel": 3.0,
            "P_rs_green_function": 2.0, "P_rs_density_build": 1.0,
            "T_rs_kernel": 3.0,
        })
        row["comparison_key"] = build_comparison_key(row)
    assert validate_pairing(cpu, gpu)["eligible"]
    mismatch = dict(gpu, recursion_depth=22)
    mismatch["comparison_key"] = build_comparison_key(mismatch)
    assert "recursion_depth_mismatch" in validate_pairing(cpu, mismatch)["reasons"]

    parsed = parse_probe_output(
        "SCF_B0C status=UNSUPPORTED reason=real_space_cuda_fallback\n"
        "SCF_B0C_RESULT scf_route=real_space rs_solver=block\n"
    )
    assert parsed["unsupported_reason"] == "real_space_cuda_fallback"


def test_outputs_have_one_canonical_dataset_and_required_fields() -> None:
    campaign = {"schema": SCHEMA, "environment": {}, "rows": [_row(), _row(backend="cuda")], "pairings": [], "skipped_rows": []}
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "campaign.json"
        write_outputs(campaign, output)
        document = json.loads(output.read_text(encoding="utf-8"))
        assert document["schema"] == SCHEMA
        with output.with_suffix(".csv").open(newline="", encoding="utf-8") as handle:
            assert set(CSV_FIELDS).issubset(csv.DictReader(handle).fieldnames or [])
        markdown = output.with_suffix(".md").read_text(encoding="utf-8")
        assert "Profile closure" in markdown
        assert "Headline eligibility" in markdown
        assert list((output.parent / "correctness").glob("*.json"))
        assert list((output.parent / "iteration_history").glob("*.json"))


def test_material_staging_contract() -> None:
    with tempfile.TemporaryDirectory() as directory:
        spec = stage_material_fixture("si", Path(directory) / "si")
        assert spec["fixture_id"] == "diamondSi_reciprocal_scf"
        assert (Path(directory) / "si" / "input.nml").exists()


if __name__ == "__main__":
    for test in (test_numeric_modes_are_end_to_end, test_profile_closure_and_misc_gate,
                 test_pairing_contract_rejects_physics_start_precision_and_fallback,
                 test_outputs_have_one_canonical_dataset_and_required_fields,
                 test_material_staging_contract):
        test()
    print("PASS: SCF-B0C harness contract")
