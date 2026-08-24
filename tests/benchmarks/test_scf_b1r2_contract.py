from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import scf_b1r2


def test_case_manifest_contains_fair_chebyshev_ladder_and_reciprocal_gate() -> None:
    specs = scf_b1r2.case_specs()
    cpu = [item for item in specs if item["lane"] == "chebyshev" and item["kind"] == "timing" and item["backend"] == "lapack"]
    assert {item["size"] for item in cpu} == {1, 2, 3, 4}
    assert {item["omp_threads"] for item in cpu} == {1, 2, 4, 8}
    assert all(item["warmups"] == 2 and item["repetitions"] == 5 for item in cpu)
    assert all(item["nstep"] == 1 for item in scf_b1r2.case_specs() if item["kind"] == "timing")
    assert {item["case_id"] for item in specs if item["lane"] == "reciprocal" and item["kind"] == "common_state"} == {
        "recip_fe2_common_cpu", "recip_fe2_common_gpu", "recip_fe3_common_cpu", "recip_fe3_common_gpu",
    }


def test_lean_tier_is_small_and_one_shot() -> None:
    specs = scf_b1r2.case_specs("lean")
    assert len(specs) == 18
    assert all(item["tier"] == "lean" for item in specs)
    assert all(item["warmups"] == 0 and item["repetitions"] == 1 for item in specs)
    assert {item["size"] for item in specs if item["lane"] == "chebyshev"} == {1, 2}
    assert {item["omp_threads"] for item in specs if item["case_id"] == "cheb_si1_cpu_omp1"} == {1}
    assert {item["omp_threads"] for item in specs if item["case_id"].startswith("cheb_si2_cpu")} == {8}


def test_si_supercell_fixture_is_deterministic_and_expands_sites(tmp_path: Path) -> None:
    fixture = scf_b1r2.stage_si_fixture(2, tmp_path / "si2")
    assert fixture["supercell"] == "2x2x2"
    assert fixture["common_state_id"].startswith("si2-")
    input_text = (tmp_path / "si2" / "input.nml").read_text(encoding="utf-8")
    assert "ntype = 16" in input_text
    assert "label(16) = 'Si16'" in input_text
    assert (tmp_path / "si2" / "Si16.nml").is_file()
    lattice = (tmp_path / "si2" / "lattice.nml").read_text(encoding="utf-8")
    assert "ntot = 16" in lattice
    assert "a(:, 1) = 1.0000000000, 1.0000000000, 0.0000000000" in lattice


def test_common_state_dos_gate_accepts_identical_grid() -> None:
    state = {
        "status": "PASS", "Natom": 2,
        "final_state": {"fermi_energy": 0.1, "total_charge": 8.0, "site_charge_transfer": 0.0,
                         "site_moment": 0.0, "final_residual": 1.0e-5, "final_total_energy": -10.0},
        "dos": {"grid": [0.0, 1.0], "dos": [1.0, 2.0]},
    }
    other = json.loads(json.dumps(state))
    result = scf_b1r2.compare_common(state, other, require_dos=True)
    assert result["status"] == "PASS"
    assert result["dos"]["relative_l2"] == 0.0


def test_k1_timer_policy_is_explicit() -> None:
    campaign = scf_b1r2.write_campaign_outputs
    assert callable(campaign)
    assert scf_b1r2.SCHEMA == "rslmto.scf-b1r2.v2"
