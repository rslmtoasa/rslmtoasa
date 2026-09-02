#!/usr/bin/env python3
"""Collect the reproducible TDDFT-11 Fe/Ni probe into one JSON record.

The collector intentionally treats the current material runs as evidence, not
as a release waiver.  It compares the emitted bare chi0 matrix rows from all
three providers and records which higher-level gates were not exercised by the
bounded selected-point decks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROUTES = ("eigenpairs", "kspace_lehmann", "realspace_gf")
MATERIALS = {
    "bccFe": "fe",
    "fccNi": "ni",
}
POINTWISE_TOLERANCE = 1.0e-8
GOLDSTONE_TOLERANCE = 1.0e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("# ") or " = " not in line:
            continue
        key, value = line[2:].split(" = ", 1)
        values[key.strip()] = value.strip()
    return values


def chi0_rows(path: Path) -> dict[tuple[float, int, int], complex]:
    rows: dict[tuple[float, int, int], complex] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 7 or fields[1] != "matrix":
            continue
        omega = float(fields[0])
        left = int(fields[2])
        right = int(fields[3])
        rows[(omega, left, right)] = complex(float(fields[4]), float(fields[5]))
    if not rows:
        raise ValueError(f"no matrix rows found in {path}")
    return rows


def scalar_from_file(path: Path, key: str) -> float | None:
    prefix = key + " "
    for line in path.read_text().splitlines():
        if line.startswith(prefix):
            fields = line.split()
            return float(fields[1])
    return None


def route_prefix(stem: str, route: str) -> str:
    if route == "kspace_lehmann":
        route = "kspace_lehmann"
    return f"tddft11_{stem}_{route}"


def collect_material(root: Path, material: str, stem: str) -> dict[str, Any]:
    directory = root / "tests" / "regression" / "tddft_validation" / "materials" / material
    routes: dict[str, Any] = {}
    matrices: dict[str, list[dict[tuple[float, int, int], complex]]] = {}

    for route in ROUTES:
        prefix = route_prefix(stem, route)
        q_files = sorted(directory.glob(f"{prefix}_q*_chi0.dat"))
        if len(q_files) != 3:
            raise ValueError(f"{material}/{route}: expected 3 q files, found {len(q_files)}")
        manifest = directory / f"{prefix}_manifest.dat"
        if not manifest.is_file():
            raise ValueError(f"missing manifest {manifest}")
        parsed = [chi0_rows(path) for path in q_files]
        first_metadata = metadata(q_files[0])
        routes[route] = {
            "manifest": str(manifest.relative_to(root)),
            "manifest_sha256": sha256(manifest),
            "q_files": [
                {
                    "path": str(path.relative_to(root)),
                    "sha256": sha256(path),
                }
                for path in q_files
            ],
            "q_count": len(q_files),
            "omega_count": int(first_metadata["omega_batch_size"]),
            "backend": first_metadata.get("chi0_backend_canonical"),
            "response_fermi_level_Ry": float(first_metadata["fermi_level_Ry"]),
            "response_k_mesh": [int(value) for value in first_metadata["k_mesh_shape"].split()],
            "real_space_points": int(first_metadata.get("real_space_points", "0")),
            "real_space_cutoff": float(first_metadata.get("real_space_cutoff", "0")),
            "real_space_tail_assessed": first_metadata.get("real_space_tail_assessed", "F") == "T",
        }
        matrices[route] = parsed

    pairwise: dict[str, float] = {}
    for left_index, left_route in enumerate(ROUTES):
        for right_route in ROUTES[left_index + 1 :]:
            maximum = 0.0
            for left_rows, right_rows in zip(matrices[left_route], matrices[right_route]):
                if left_rows.keys() != right_rows.keys():
                    raise ValueError(f"{material}: row grids differ for {left_route}/{right_route}")
                maximum = max(
                    maximum,
                    max(abs(left_rows[key] - right_rows[key]) for key in left_rows),
                )
            pairwise[f"{left_route}_vs_{right_route}"] = maximum

    goldstone: dict[str, Any] = {}
    for route in ("eigenpairs", "kspace_lehmann"):
        path = directory / f"tddft11_{stem}_{route}_goldstone.dat"
        goldstone[route] = {
            "path": str(path.relative_to(root)),
            "sha256": sha256(path),
            "raw_closest_eigenvalue": scalar_from_file(path, "raw_closest_eigenvalue"),
            "raw_residual": scalar_from_file(path, "raw_residual"),
            "raw_ward_residual": scalar_from_file(path, "raw_ward_residual"),
            "raw_dm_residual": scalar_from_file(path, "raw_dm_residual"),
        }

    max_backend_difference = max(pairwise.values())
    goldstone_residual = max(item["raw_residual"] for item in goldstone.values())
    native_route = routes["realspace_gf"]
    return {
        "deck_scope": {
            "bounded_probe": True,
            "q_points": [[0.0, 0.0, 0.0], [0.01, 0.0, 0.0], [0.02, 0.0, 0.0]],
            "omega_Ry": [0.0, 0.001, 0.002],
            "eta_Ry": 0.0005,
            "soc": False,
            "reciprocal_mode": "ham_only",
        },
        "routes": routes,
        "backend_equivalence": {
            "pairwise_max_abs_matrix_difference": pairwise,
            "max_abs_matrix_difference": max_backend_difference,
            "pointwise_tolerance": POINTWISE_TOLERANCE,
            "pass": max_backend_difference <= POINTWISE_TOLERANCE,
        },
        "goldstone": {
            "routes_with_static_diagnostic": goldstone,
            "native_realspace_static_diagnostic": "not available: native provider is bare dynamic chi0 only",
            "maximum_raw_residual": goldstone_residual,
            "gapless_tolerance": GOLDSTONE_TOLERANCE,
            "pass": goldstone_residual <= GOLDSTONE_TOLERANCE,
        },
        "native_realspace_coverage": {
            "real_space_points": native_route["real_space_points"],
            "real_space_cutoff": native_route["real_space_cutoff"],
            "tail_assessed": native_route["real_space_tail_assessed"],
            "status": "finite bounded cluster; tail convergence not assessed",
        },
        "unassessed_gates": {
            "quadratic_small_q_dispersion": "not assessed: bare chi0 selected-point probe has no enhanced magnon mode fit",
            "lkag_frozen_magnon_reference": "not assessed: no independent same-ground-state reference run in this bounded probe",
            "ni_path_backfolding": "not assessed: no dense reciprocal path/backfolding audit in this bounded probe",
            "damping_stoner_convergence": "not assessed: no eta/k-mesh/frequency convergence sweep",
        },
        "release_gate_pass": False,
        "release_gate_reason": (
            "The bounded probe exercises all three bare-chi0 routes, but the native route disagrees "
            "pointwise with the two K-space routes and the K-space raw Ward residual is nonzero; "
            "dense q^2, independent references, Ni path, and convergence sweeps remain unassessed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output JSON (default: results/validation/TDDFT-11_FE_NI/evidence.json)",
    )
    args = parser.parse_args()
    root = args.repo.resolve()
    output = args.output or root / "results" / "validation" / "TDDFT-11_FE_NI" / "evidence.json"
    evidence = {
        "campaign": "TDDFT-11",
        "date": "2026-09-02",
        "source_revision": "working tree before TDDFT-11 commit",
        "materials": {
            material: collect_material(root, material, stem)
            for material, stem in MATERIALS.items()
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
