"""Infrastructure checks for recursive manifest-level namelist defaults."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manifest_defaults import apply_manifest_defaults


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


if __name__ == "__main__":
    test_defaults_merge_recursively_without_mutating_inputs()
    test_old_manifest_without_defaults_is_unchanged()
    print("RESULT: PASS")
