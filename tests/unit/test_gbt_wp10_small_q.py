"""System-independent WP10 tests using synthetic raw sweep rows."""

from __future__ import annotations

import importlib.util
import math
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "analyze_gbt_wp10.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("gbt_wp10_small_q", TOOL)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load {TOOL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows(odd_coefficient: float = 2.0e-12):
    result = []
    sin2 = math.sin(math.radians(5.0)) ** 2
    for q in (0.0, 0.02, -0.02, 0.05, -0.05, 0.10, -0.10, 0.20, -0.20):
        value = sin2 * (0.8 * q**2 + 1.5 * q**4) + odd_coefficient * q
        result.append({
            "source": "synthetic", "line": len(result) + 1, "mode": "mft",
            "q": (0.0, 0.0, q), "theta_deg": 5.0, "mesh": (12, 12, 12),
            "delta_raw": value, "delta_gauge": value, "sin2theta": sin2,
            "gauge_available": True, "alat_angstrom": 2.8612,
        })
    return result


def test_even_odd_and_quadratic_quartic_fit() -> None:
    tool = load_tool()
    result = tool.analyze_wp10_rows(rows(), direction=(0.0, 0.0, 1.0))
    group = result["groups"][0]
    reversal = group["reversal"]
    assert reversal["complete_pair_count"] == 4
    assert reversal["odd_component"]["status"] == "within_tolerance"
    fit = reversal["fits"]["sin2"]["quadratic_quartic"]
    assert fit["status"] == "fit"
    assert abs(fit["coefficients"]["A"] - 0.8) < 1.0e-8
    assert abs(fit["coefficients"]["B"] - 1.5) < 1.0e-6
    assert fit["physical_units"]["A_Ry_angstrom2"] > 0.0
    assert len(reversal["fits"]["sin2"]["quadratic_quartic_windows"]) == 2


def test_missing_pair_and_odd_red_flag_are_explicit() -> None:
    tool = load_tool()
    missing = [item for item in rows() if item["q"] != (0.0, 0.0, -0.05)]
    result = tool.analyze_wp10_rows(missing, direction=(0.0, 0.0, 1.0))
    reversal = result["groups"][0]["reversal"]
    assert next(item for item in reversal["pairs"] if item["q_abs_internal"] == 0.05)["status"] == "missing_negative"
    assert reversal["complete_pair_count"] == 3

    result = tool.analyze_wp10_rows(rows(1.0e-3), direction=(0.0, 0.0, 1.0))
    assert result["groups"][0]["reversal"]["odd_component"]["status"] == "red_flag"


def test_metadata_alat_is_carried_into_wp10_units() -> None:
    tool = load_tool()
    text = """# schema = gbt_wp09_harmonic_v1
# mode = mft
# alat_angstrom = 2.8612
# columns = q1 q2 q3 theta_deg DeltaE_raw DeltaE_gauge sin2theta nk1 nk2 nk3 gauge_available
0 0 0 5 0 0 0.007596 12 12 12 1
0 0 .02 5 0.0001 0.0001 0.007596 12 12 12 1
0 0 -.02 5 0.0001 0.0001 0.007596 12 12 12 1
"""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "sweep.dat"
        path.write_text(text)
        parsed = tool._ANALYZER.read_sweep_file(path)
    assert parsed[0]["alat_angstrom"] == 2.8612
    result = tool.analyze_wp10_rows(parsed)
    assert result["physical_q_units"] == "1/angstrom"


if __name__ == "__main__":
    test_even_odd_and_quadratic_quartic_fit()
    test_missing_pair_and_odd_red_flag_are_explicit()
    test_metadata_alat_is_carried_into_wp10_units()
    print("GBT WP10 small-q unit checks PASS")
