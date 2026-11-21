"""Infrastructure checks for recursive manifest-level namelist defaults."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manifest_defaults import apply_manifest_defaults, resolve_extends


def test_defaults_merge_recursively_without_mutating_inputs() -> None:
    manifest = {
        "defaults": {
            "namelists": {
                "lattice": {
                    "strux_backend": "strux_lib",
                    "screening": {"mode": "manual", "sigma": 0.7},
                },
                "control": {"nsp": 1, "recur": "block"},
            }
        }
    }
    case = {
        "name": "sentinel",
        "namelists": {
            "lattice": {"strux_backend": "legacy", "screening": {"sigma": 0.9}},
            "control": {"nsp": 2},
            "energy": {"fermi": -0.04},
        },
    }

    merged = apply_manifest_defaults(manifest, case)

    assert merged["namelists"] == {
        "lattice": {
            "strux_backend": "legacy",
            "screening": {"mode": "manual", "sigma": 0.9},
        },
        "control": {"nsp": 2, "recur": "block"},
        "energy": {"fermi": -0.04},
    }
    assert case["namelists"]["lattice"]["screening"] == {"sigma": 0.9}
    assert manifest["defaults"]["namelists"]["lattice"]["screening"]["sigma"] == 0.7


def test_old_manifest_without_defaults_is_unchanged() -> None:
    case = {"name": "old", "namelists": {"control": {"nsp": 2}}}
    merged = apply_manifest_defaults({}, case)
    assert merged == case
    assert merged is not case


def test_extends_merges_child_over_parent() -> None:
    """Bug caught: lists merged instead of replaced, dicts replaced instead of merged, parent mutated."""
    manifest = {"cases": [
        {"name": "base", "labels": ["a", "b"], "tier": "quick", "namelists": {"control": {"nsp": 1, "recur": "block"}}},
        {"name": "mid", "extends": "base", "namelists": {"control": {"nsp": 2}}},
        {"name": "leaf", "extends": "mid", "labels": ["c"], "tier": "standard",
         "namelists": {"reciprocal": {"use_symmetry_reduction": True}}},
    ]}
    cases = {c["name"]: c for c in resolve_extends(manifest)["cases"]}
    assert cases["base"] == manifest["cases"][0]
    assert cases["mid"] == {"name": "mid", "labels": ["a", "b"], "tier": "quick",
                            "namelists": {"control": {"nsp": 2, "recur": "block"}}}
    assert cases["leaf"] == {"name": "leaf", "labels": ["c"], "tier": "standard",
                             "namelists": {"control": {"nsp": 2, "recur": "block"},
                                           "reciprocal": {"use_symmetry_reduction": True}}}
    assert manifest["cases"][0]["labels"] == ["a", "b"] and manifest["cases"][1]["extends"] == "base"


def test_extends_without_extends_is_identity() -> None:
    manifest = {"defaults": {"namelists": {"control": {"nsp": 1}}}, "cases": [{"name": "x", "labels": ["a"]}]}
    assert resolve_extends(manifest) == manifest


def _raises(manifest: dict, text: str) -> None:
    try:
        resolve_extends(manifest)
    except ValueError as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_extends_cycle_and_missing_parent_are_errors() -> None:
    _raises({"cases": [{"name": "a", "extends": "b"}, {"name": "b", "extends": "a"}]}, "cycle")
    _raises({"cases": [{"name": "a", "extends": "a"}]}, "cycle")
    _raises({"cases": [{"name": "a", "extends": "nope"}]}, "missing case 'nope'")


if __name__ == "__main__":
    test_defaults_merge_recursively_without_mutating_inputs()
    test_old_manifest_without_defaults_is_unchanged()
    test_extends_merges_child_over_parent()
    test_extends_without_extends_is_identity()
    test_extends_cycle_and_missing_parent_are_errors()
    print("RESULT: PASS")
