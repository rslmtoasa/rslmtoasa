"""Fast, system-independent tests for the WP09 sweep schema and analysis.

The numerical values below are synthetic by design.  They exercise the
analysis contract without turning a material-specific Fe energy into a unit
test golden.  Production runs are collated by the WP09 validation runner and
retain their raw values in the JSON report.
"""

from __future__ import annotations

import importlib.util
import math
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source" / "calculation.f90"
TOOL = ROOT / "tools" / "analyze_gbt_wp09.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("gbt_wp09_analysis", TOOL)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load {TOOL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_rows():
    rows = []
    for mode in ("mft", "mft_constrained"):
        for nk in (8, 12, 16):
            for theta in (2.0, 5.0, 10.0, 15.0, 20.0):
                sin2 = math.sin(math.radians(theta)) ** 2
                e_ref = -1.0 + 2.0e-4 * sin2
                # The finite mesh drift is deliberately q-only and therefore
                # must disappear from the same-q gauge subtraction.
                drift = {8: 2.0e-4, 12: 5.0e-5, 16: 1.0e-5}[nk]
                e_q0 = -1.0 + drift if theta else -1.0
                nonlinear = 1.0 if theta < 20.0 else 1.20
                physical = 4.0e-3 * nonlinear * sin2
                if mode == "mft_constrained":
                    # Corrected MFT has no same-q zero-cone gauge probe in
                    # this schema; its physical observable is the raw
                    # fixed-potential energy difference.
                    e_q0_value = None
                    e_qref0 = None
                    gauge = False
                else:
                    e_q0_value = e_q0
                    e_qref0 = -1.0
                    gauge = True
                rows.append({
                    "source": "synthetic",
                    "line": len(rows) + 1,
                    "mode": mode,
                    "q": (0.0, 0.0, 0.0),
                    "theta_deg": theta,
                    "mesh": (nk, nk, nk),
                    "e_q_theta": e_ref,
                    "e_q0": -1.0 if mode == "mft" else None,
                    "e_qref_theta": e_ref,
                    "e_qref0": e_qref0,
                    "delta_raw": 0.0,
                    "delta_gauge": 0.0 if mode == "mft" else None,
                    "delta_pure": 0.0 if mode == "mft" else None,
                    "sin2theta": sin2,
                    "mtot": 2.0,
                    "omega": 0.0,
                    "fermi_level": -0.04,
                    "electron_count": 8.0,
                    "electron_error": 0.0,
                    "target_electrons": 8.0,
                    "weight_sum": 1.0,
                    "gauge_available": gauge,
                    "q_half_maps_to_mesh": True,
                    "q_half_mapping_source": "synthetic",
                })
                rows.append({
                    "source": "synthetic",
                    "line": len(rows) + 1,
                    "mode": mode,
                    "q": (0.0, 0.0, 0.05),
                    "theta_deg": theta,
                    "mesh": (nk, nk, nk),
                    "e_q_theta": e_q0 + physical if mode == "mft" else e_ref + physical,
                    "e_q0": e_q0_value,
                    "e_qref_theta": e_ref,
                    "e_qref0": e_qref0,
                    "delta_raw": (e_q0 + physical - e_ref) if mode == "mft" else physical,
                    "delta_gauge": physical if mode == "mft" else None,
                    "delta_pure": drift if mode == "mft" else None,
                    "sin2theta": sin2,
                    "mtot": 2.0,
                    "omega": physical / sin2,
                    "fermi_level": -0.04,
                    "electron_count": 8.0,
                    "electron_error": 0.0,
                    "target_electrons": 8.0,
                    "weight_sum": 1.0,
                    "gauge_available": gauge,
                    "q_half_maps_to_mesh": False,
                    "q_half_mapping_source": "synthetic",
                })
    return rows


def find(result, **wanted):
    for item in result:
        if all(item.get(key) == value for key, value in wanted.items()):
            return item
    raise AssertionError(f"no result matching {wanted!r}")


def test_plateau_excludes_nonlinear_angle() -> None:
    tool = load_tool()
    result = tool.analyze_rows(synthetic_rows(), relative_tolerance=0.05)
    plateau = find(result["plateaus"], mode="mft", q=[0.0, 0.0, 0.05], mesh=[12, 12, 12], definition="gauge")
    assert plateau["status"] == "plateau"
    assert plateau["admitted_angles_deg"] == [2.0, 5.0, 10.0, 15.0]
    assert plateau["excluded_angles_deg"] == [20.0]
    assert plateau["spread_rel"] < 0.05

    corrected = find(result["plateaus"], mode="mft_constrained", q=[0.0, 0.0, 0.05], mesh=[12, 12, 12], definition="raw")
    assert corrected["admitted_angles_deg"] == [2.0, 5.0, 10.0, 15.0]
    assert not any(item["definition"] == "gauge" and item["mode"] == "mft_constrained" for item in result["plateaus"])


def test_kgrid_pure_drift_and_electron_consistency_are_reported() -> None:
    tool = load_tool()
    result = tool.analyze_rows(synthetic_rows())
    convergence = find(result["k_grid_convergence"], mode="mft", q=[0.0, 0.0, 0.05], theta_deg=5.0, definition="gauge")
    assert convergence["status"] == "converged"
    assert len(convergence["meshes"]) == 3
    assert convergence["spread_abs"] < 1.0e-12
    raw_convergence = find(result["k_grid_convergence"], mode="mft", q=[0.0, 0.0, 0.05], theta_deg=5.0, definition="raw")
    assert raw_convergence["status"] == "not_converged"

    drift = find(result["pure_gauge_drift"], mode="mft", q=[0.0, 0.0, 0.05], theta_deg=5.0, mesh=[8, 8, 8])
    assert abs(drift["drift"] - 2.0e-4) < 1.0e-14
    consistency = find(result["electron_consistency"], mode="mft", q=[0.0, 0.0, 0.05], theta_deg=5.0)
    assert consistency["max_abs_electron_error"] == 0.0
    assert consistency["fermi_min"] == consistency["fermi_max"] == -0.04


def test_schema_reader_and_production_output_contract() -> None:
    tool = load_tool()
    text = """# schema = gbt_wp09_harmonic_v1
# mode = mft
# columns = q1 q2 q3 theta_deg E_q_theta E_q0 E_qref_theta E_qref0 DeltaE_raw DeltaE_gauge DeltaE_pure sin2theta Mtot omega fermi_level electron_count electron_error target_electrons weight_sum nk1 nk2 nk3 nk_total gauge_available
0 0 0 5 -1.0 -1.0 -1.0 -1.0 0 0 0 0.007596 2 0 -0.04 8 0 8 1 8 8 8 512 1
0 0 .05 5 -0.9999696 -1.0 -1.0 -1.0 0.0000304 0.0000304 0 0.007596 2 0.008  -0.04 8 0 8 1 8 8 8 512 1
"""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "sweep.dat"
        path.write_text(text)
        rows = tool.read_sweep_file(path)
    assert len(rows) == 2
    assert rows[1]["gauge_available"] is True
    assert rows[1]["mesh"] == (8, 8, 8)
    assert abs(rows[1]["delta_gauge"] - 3.04e-5) < 1.0e-12

    source = SOURCE.read_text()
    for needle in (
        "frozen_magnon_harmonic_diagnostics.dat",
        "schema = gbt_wp09_harmonic_v1",
        "fermi_q",
        "electron_count_q",
        "nk_total_q",
        "gauge_available",
    ):
        assert needle in source, f"missing WP09 production-output contract: {needle!r}"


def _small_q_rows():
    rows = []
    sin2 = math.sin(math.radians(5.0)) ** 2
    for q in (0.0, 0.02, -0.02, 0.05, -0.05, 0.10, -0.10, 0.20, -0.20):
        response = sin2 * (0.8 * q**2 + 1.5 * q**4) + 2.0e-12 * q
        rows.append({
            "source": "synthetic_small_q", "line": len(rows) + 1, "mode": "mft",
            "q": (0.0, 0.0, q), "theta_deg": 5.0, "mesh": (12, 12, 12),
            "delta_raw": response, "delta_gauge": response, "sin2theta": sin2,
            "gauge_available": True, "alat_angstrom": 2.8612,
        })
    return rows


def test_small_q_even_odd_fits_and_units() -> None:
    tool = load_tool()
    result = tool.analyze_wp10_rows(_small_q_rows(), direction=(0.0, 0.0, 1.0))
    assert result["schema"] == "gbt_wp10_small_q_v1"
    assert result["q_scale_inv_angstrom"] is not None
    group = result["groups"][0]
    assert group["reversal"]["complete_pair_count"] == 4
    assert group["reversal"]["odd_component"]["status"] == "within_tolerance"
    fit = group["reversal"]["fits"]["sin2"]["quadratic_quartic"]
    assert abs(fit["coefficients"]["A"] - 0.8) < 1.0e-8
    assert abs(fit["coefficients"]["B"] - 1.5) < 1.0e-6
    assert fit["physical_units"]["A_Ry_angstrom2"] > 0.0
    assert len(group["reversal"]["fits"]["sin2"]["quadratic_quartic_windows"]) == 2


def test_small_q_missing_reversal_is_not_fitted() -> None:
    tool = load_tool()
    rows = [row for row in _small_q_rows() if row["q"] != (0.0, 0.0, -0.05)]
    result = tool.analyze_wp10_rows(rows, direction=(0.0, 0.0, 1.0))
    pairs = result["groups"][0]["reversal"]["pairs"]
    missing = next(item for item in pairs if item["q_abs_internal"] == 0.05)
    assert missing["status"] == "missing_negative"
    assert result["groups"][0]["reversal"]["complete_pair_count"] == 3


if __name__ == "__main__":
    test_plateau_excludes_nonlinear_angle()
    test_kgrid_pure_drift_and_electron_consistency_are_reported()
    test_schema_reader_and_production_output_contract()
    test_small_q_even_odd_fits_and_units()
    test_small_q_missing_reversal_is_not_fitted()
    print("GBT WP09 harmonic/k-grid unit checks PASS")
