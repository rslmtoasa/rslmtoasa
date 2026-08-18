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
    print("PASS: ACC-00 benchmark harness contract")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
