#!/usr/bin/env python3
"""Small tooling test for the ACC-00 schema/parser/comparison contract."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_harness import SCHEMA, compare_documents, parse_profile_output  # noqa: E402


def main() -> int:
    output = """\
PROFILE_DIMENSIONS bccFe_one_site sites=1 spinor_basis=9 nk=16 mesh=4 x 4 x 1 nw=96
PROFILE_RECIPROCAL bccFe_one_site fourier_assembly=1.0E-02 k_eigensolution=2.0E-02 arbitrary_kq_assembly_eigensolution=3.0E-02 pair_operator_construction=4.0E-02
PROFILE_MEMORY_MIB bccFe_one_site hk=1.0 normal_eigenpairs=2.0 arbitrary_kq_eigenpairs=3.0 pair_operators_and_workspace=4.0 response=5.0 principal_payload=6.0
"""
    records = parse_profile_output(output)
    assert len(records) == 1
    assert records[0]["metadata"]["matrix_dimension"] == 9
    assert records[0]["metrics"]["eigensolver_s"] == 0.02
    assert records[0]["metrics"]["principal_payload_mib"] == 6.0

    acc06 = parse_profile_output(
        "ACC06_DIMENSIONS fixture=Si_sp backend=cuda strategy=backend sites=1 "
        "matrix_dimension=8 nk=16 tile_size=4 eigenvectors=1 lmax=1\n"
        "ACC06_TIMING fixture=Si_sp assembly_s=1.0E-03 solve_s=2.0E-03 total_s=3.0E-03\n"
    )
    assert acc06[0]["metadata"]["backend"] == "cuda"
    assert acc06[0]["metadata"]["tile_size"] == 4
    assert acc06[0]["metrics"]["total_s"] == 3.0e-3

    accp0 = parse_profile_output(
        "ACCP0_DIMENSIONS fixture=diamondSi source=scf workload=crossover backend=cuda "
        "L=1 natom=2 nmat=16 nominal_mesh=8x8x8 actual_unique_nk=512 tile=16 eigenvectors=0\n"
        "ACCP0_TIMING fixture=diamondSi backend=cuda cold_process_wall_s=0.4 "
        "cuda_context_backend_init_s=0.2 first_solve_s=0.1 steady_solve_median_s=0.03 "
        "Hk_CPU_s=0.01 H2D_s=0.002 solver_s=0.02 D2H_s=0.001 total_steady_s=0.04 "
        "memory_estimate_mib=3 memory_free_before_mib=1000 memory_total_mib=16000\n"
    )
    assert accp0[0]["metadata"]["actual_unique_nk"] == 512
    assert accp0[0]["metadata"]["nominal_mesh"] == "8x8x8"
    assert accp0[0]["metrics"]["cold_process_wall_s"] == 0.4
    assert accp0[0]["metrics"]["total_steady_s"] == 0.04

    base = {
        "schema": SCHEMA,
        "benchmark": {
            "name": "cpu",
            "metadata": {"git_commit": "a", "compiler": "gfortran"},
            "samples": [{"wall_time_s": 4.0}, {"wall_time_s": 2.0}],
        },
    }
    candidate = {
        "schema": SCHEMA,
        "benchmark": {
            "name": "gpu",
            "metadata": {"git_commit": "b", "compiler": "gfortran"},
            "samples": [{"wall_time_s": 1.0}, {"wall_time_s": 1.5}],
        },
    }
    report = compare_documents(base, candidate)
    assert report["metrics"][0]["speedup"] == 2.4
    assert report["environment_mismatch_warnings"]

    with tempfile.TemporaryDirectory() as directory:
        output_path = Path(directory) / "record.json"
        command = [
            sys.executable,
            str(ROOT / "tests/benchmarks/benchmark_harness.py"),
            "run",
            "--name",
            "smoke",
            "--warmups",
            "0",
            "--repetitions",
            "2",
            "--output",
            str(output_path),
            "--command",
            sys.executable,
            "-c",
            "print('ok')",
        ]
        subprocess.run(command, check=True, cwd=ROOT, capture_output=True, text=True)
        document = json.loads(output_path.read_text(encoding="utf-8"))
        assert document["schema"] == SCHEMA
        assert len(document["benchmark"]["samples"]) == 2

        persistent_path = Path(directory) / "persistent.json"
        persistent_command = [
            sys.executable,
            str(ROOT / "tests/benchmarks/benchmark_harness.py"),
            "run",
            "--name",
            "persistent-smoke",
            "--persistent",
            "--output",
            str(persistent_path),
            "--command",
            sys.executable,
            "-c",
            "print('ACCP0_DIMENSIONS fixture=smoke backend=lapack L=1 natom=1 nmat=1 nominal_mesh=1x1x1 actual_unique_nk=1 tile=1 eigenvectors=0'); print('ACCP0_TIMING fixture=smoke backend=lapack total_steady_s=0.01')",
        ]
        subprocess.run(persistent_command, check=True, cwd=ROOT, capture_output=True, text=True)
        persistent_document = json.loads(persistent_path.read_text(encoding="utf-8"))
        assert persistent_document["benchmark"]["policy"]["persistent_process"] is True
        assert len(persistent_document["benchmark"]["samples"]) == 1
    print("PASS: ACC-00 benchmark harness contract")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
