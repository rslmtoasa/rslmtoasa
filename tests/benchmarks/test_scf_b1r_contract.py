from __future__ import annotations

import json
import sys
from pathlib import Path

import scf_b1r
from scf_b1_report import eligible_ratio


def test_release_preflight_uses_last_effective_optimization_flag(tmp_path: Path) -> None:
    build = tmp_path / "build"
    build.mkdir()
    (build / "CMakeCache.txt").write_text(
        "CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran\n",
        encoding="utf-8",
    )
    commands = [
        f"gfortran -O3 -c {scf_b1r.ROOT / 'source' / 'main.f90'} -O0",
        f"gfortran -O3 -c {scf_b1r.ROOT / 'source' / 'green_block.f90'} -O0",
        f"gfortran -O3 -c {scf_b1r.ROOT / 'source' / 'hamiltonian_build.f90'} -O0",
        f"gfortran -O3 -c {scf_b1r.ROOT / 'source' / 'reciprocal.f90'} -O0",
    ]
    (build / "compile_commands.json").write_text(
        json.dumps([{"directory": str(build), "command": command, "file": command.split()[-3]} for command in commands]),
        encoding="utf-8",
    )
    result = scf_b1r.build_preflight(build)
    assert result["status"] == "FAIL"
    assert "trailing -O0" in result["reason"]


def test_ineligible_ratios_are_not_reportable() -> None:
    row = {"S_iteration": 4.0, "equal_precision_eligible": False}
    assert scf_b1r._case_specs()
    assert eligible_ratio(row, "S_iteration") is None


def test_reconfigure_seeds_release_optimization() -> None:
    args = scf_b1r.configure_args_from_cache({})
    assert "-DCMAKE_Fortran_FLAGS_RELEASE=-O3" in args


def test_command_logging_tolerates_non_utf8_output(tmp_path: Path) -> None:
    evidence = tmp_path / "evidence"
    result = scf_b1r.run_logged_command(
        [sys.executable, "-c", "import sys; sys.stdout.buffer.write(bytes([251, 10]))"], evidence
    )
    assert result["returncode"] == 0
    assert "\ufffd" in (evidence / "stdout.log").read_text(encoding="utf-8")
    assert (evidence / "stdout.bin").read_bytes() == bytes([251, 10])


def test_resume_retries_failed_case(tmp_path: Path) -> None:
    state_path = tmp_path / "raw" / "case" / "case.json"
    state_path.parent.mkdir(parents=True)
    state_path.write_text(json.dumps({"completed": True, "aggregate": {"status": "FAIL"}}), encoding="utf-8")
    assert scf_b1r._load_completed_case(tmp_path, "case", False) is None

    valid_but_ineligible = {
        "completed": True,
        "aggregate": {
            "status": "FAIL",
            "repetition_rows": [{"returncode": 0, "final_state": {"nmat": 18}}],
        },
    }
    state_path.write_text(json.dumps(valid_but_ineligible), encoding="utf-8")
    assert scf_b1r._load_completed_case(tmp_path, "case", False) == valid_but_ineligible
