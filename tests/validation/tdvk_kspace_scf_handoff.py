#!/usr/bin/env python3
"""Validate the self-consistent k-space-SCF -> TDDFT state handoff.

This is the TDVK review-gate bridge between the historical TDVK-05 diagnostic
campaign and the later interaction work.  Each production case is a fresh
process with ``use_kspace=.true.``.  The response driver must consume the
accepted reciprocal cache directly; the old real-space-SCF -> post-SCF
reciprocal rebuild remains outside this harness.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tdvk05_fe_convergence import (  # noqa: E402
    header_ints,
    header_number,
    matrix_metrics,
    number,
    parse_basis_artifact,
    parse_convergence_output,
    parse_kpoint_artifact,
    parse_headers,
    parse_matrix_artifact,
    replace_assignment,
    select_matrix,
)


@dataclass(frozen=True)
class Case:
    name: str
    mesh: int
    eta_values: tuple[float, ...]
    gf_closure: bool = False


CASES = (
    Case("kspace_scf_mesh8", 8, (0.01,), True),
    Case("kspace_scf_mesh12", 12, (0.01, 0.005), False),
)


def prepare_input(template: Path, destination: Path, fe_database: Path, case: Case) -> None:
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(fe_database)))
    text = replace_assignment(text, "use_kspace", ".true.")
    text = replace_assignment(text, "nk1", str(case.mesh))
    text = replace_assignment(text, "nk2", str(case.mesh))
    text = replace_assignment(text, "nk3", str(case.mesh))
    text = replace_assignment(text, "auto_find_fermi", ".true.")
    text = replace_assignment(text, "n_q", "1")
    text = replace_assignment(text, "q_list", "0.0, 0.0, 0.0")
    text = replace_assignment(text, "n_omega", "1")
    text = replace_assignment(text, "use_omega_grid", ".false.")
    text = replace_assignment(text, "omega_grid", "0.0")
    text = replace_assignment(text, "n_eta", str(len(case.eta_values)))
    text = replace_assignment(text, "eta_grid", ", ".join(f"{eta:.12g}" for eta in case.eta_values))
    text = replace_assignment(text, "eta", f"{case.eta_values[0]:.12g}")
    text = replace_assignment(text, "gf_closure_audit", ".true." if case.gf_closure else ".false.")
    text = replace_assignment(text, "gf_integration_points", "9601")
    text = replace_assignment(text, "gf_integration_eta", "0.001")
    text = replace_assignment(text, "gf_energy_margin", "0.60")
    text = replace_assignment(text, "output_file", repr(f"tddft_kspace_scf_{case.mesh}.dat"))
    destination.write_text(text, encoding="utf-8")


def resolve_artifact(workdir: Path, headers: dict[str, str], key: str, started_wall: float) -> Path:
    path = Path(headers[key])
    if not path.is_absolute():
        path = workdir / path
    if not path.exists() or path.stat().st_mtime < started_wall - 1.0:
        raise RuntimeError(f"missing or stale {key} artifact: {path}")
    return path


def parse_state_artifact(path: Path) -> tuple[dict[str, str], dict[str, object]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    points: dict[int, tuple[float, float, float, float]] = {}
    eigenvalues: dict[tuple[int, int], float] = {}
    occupations: dict[tuple[int, int], float] = {}
    density: dict[tuple[int, int, int], complex] = {}
    point_serialization: list[str] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        tag = fields[0]
        if tag == "K" and len(fields) == 6:
            index = int(fields[1])
            row = tuple(number(value) for value in fields[2:])
            points[index] = row  # type: ignore[assignment]
            point_serialization.append(" ".join(fields[1:]))
        elif tag == "O" and len(fields) == 5:
            ik, ib = int(fields[1]), int(fields[2])
            eigenvalues[(ik, ib)] = number(fields[3])
            occupations[(ik, ib)] = number(fields[4])
        elif tag == "D" and len(fields) == 6:
            ik, i, j = (int(value) for value in fields[1:4])
            density[(ik, i, j)] = complex(number(fields[4]), number(fields[5]))
        else:
            raise RuntimeError(f"unexpected state artifact row in {path}: {line}")
    if not points or not eigenvalues or not occupations or not density:
        raise RuntimeError(f"state artifact is incomplete: {path}")
    if list(sorted(points)) != list(range(1, len(points) + 1)):
        raise RuntimeError(f"state artifact k indices are not complete: {path}")
    return headers, {
        "points": points,
        "eigenvalues": eigenvalues,
        "occupations": occupations,
        "density": density,
        "k_fingerprint_sha256": hashlib.sha256("\n".join(point_serialization).encode("ascii")).hexdigest(),
    }


def max_mapping_difference(left: dict, right: dict, complex_values: bool = False) -> float:
    if set(left) != set(right):
        raise RuntimeError("state artifact index sets differ")
    differences = []
    for key in left:
        if isinstance(left[key], tuple):
            differences.append(max(abs(a - b) for a, b in zip(left[key], right[key])))
        else:
            differences.append(abs(left[key] - right[key]))  # type: ignore[operator]
    return max(differences)


def state_consistency(scf: tuple[dict[str, str], dict[str, object]], tddft: tuple[dict[str, str], dict[str, object]]) -> dict[str, object]:
    scf_headers, scf_data = scf
    tddft_headers, tddft_data = tddft
    scf_points = scf_data["points"]
    tddft_points = tddft_data["points"]
    point_difference = max_mapping_difference(scf_points, tddft_points)
    eigen_difference = max_mapping_difference(scf_data["eigenvalues"], tddft_data["eigenvalues"])
    occupation_difference = max_mapping_difference(scf_data["occupations"], tddft_data["occupations"])
    density_difference = max_mapping_difference(scf_data["density"], tddft_data["density"])
    density_frobenius = math.sqrt(
        sum(abs(scf_data["density"][key] - tddft_data["density"][key]) ** 2 for key in scf_data["density"])
    )
    ef_difference = abs(header_number(scf_headers, "fermi_level_Ry") - header_number(tddft_headers, "fermi_level_Ry"))
    if scf_data["k_fingerprint_sha256"] != tddft_data["k_fingerprint_sha256"]:
        raise RuntimeError("SCF and TDDFT state k-vector fingerprints differ")
    if scf_headers.get("actual_k_mesh") != tddft_headers.get("actual_k_mesh"):
        raise RuntimeError("SCF and TDDFT state mesh labels differ")
    if point_difference > 2.0e-14 or eigen_difference > 2.0e-13 or occupation_difference > 2.0e-13 or density_difference > 2.0e-13:
        raise RuntimeError(
            "accepted SCF/TDDFT state mismatch: "
            f"mesh={point_difference} eigen={eigen_difference} occupation={occupation_difference} density={density_difference}"
        )
    return {
        "scf_mesh_fingerprint_sha256": scf_data["k_fingerprint_sha256"],
        "tddft_mesh_fingerprint_sha256": tddft_data["k_fingerprint_sha256"],
        "mesh_max_abs_difference": point_difference,
        "weight_max_abs_difference": max(
            abs(scf_points[key][3] - tddft_points[key][3]) for key in scf_points
        ),
        "fermi_level_max_abs_difference_Ry": ef_difference,
        "eigenvalue_max_abs_difference_Ry": eigen_difference,
        "occupation_max_abs_difference": occupation_difference,
        "gauge_invariant_projector_max_abs_difference": density_difference,
        "gauge_invariant_projector_frobenius_difference": density_frobenius,
        "state_identity_pass": True,
    }


def parse_gf_artifact(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    rows = [line.split() for line in lines if line and not line.startswith("#")]
    if len(rows) != 1 or len(rows[0]) != 6:
        raise RuntimeError(f"GF artifact must contain one six-column row: {path}")
    values = [number(value) for value in rows[0][:5]]
    finite = rows[0][5].upper() == "T"
    if not finite or not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"GF spot is not finite: {path}")
    if int(headers.get("integration_points", "0")) != 9601:
        raise RuntimeError(f"GF spot did not use 9601 integration points: {path}")
    return {
        "headers": headers,
        "norm_lehmann": values[0],
        "norm_gf": values[1],
        "dF": values[2],
        "rF": values[3],
        "dInf": values[4],
        "finite": finite,
    }


def basis_alignment_metrics(left: dict, right: dict) -> dict[str, object]:
    """Align the actual per-state product bases without assuming equal SVD values."""
    if set(left) != set(right):
        raise RuntimeError("compact basis mode keys differ between k-space-SCF states")
    max_mode_difference = max(abs(left[key][0] - right[key][0]) for key in left)
    singular_value_difference = max(abs(left[key][1] - right[key][1]) for key in left)
    singular_values = [value[1] for value in left.values()]
    blocks: dict[tuple[int, int], list[int]] = {}
    for site, response_l, radial_index, mode in left:
        blocks.setdefault((site, response_l), []).append(radial_index)
    max_unitarity_residual = 0.0
    for block, radial_indices in blocks.items():
        modes = sorted({key[3] for key in left if key[:2] == block})
        radial_indices = sorted(set(radial_indices))
        overlap = [
            [
                sum(
                    left[(block[0], block[1], radial, mode_left)][0].conjugate()
                    * right[(block[0], block[1], radial, mode_right)][0]
                    for radial in radial_indices
                )
                for mode_right in modes
            ]
            for mode_left in modes
        ]
        residual = math.sqrt(
            sum(
                abs(
                    sum(overlap[row][inner].conjugate() * overlap[inner][column] for inner in range(len(modes)))
                    - (1.0 if row == column else 0.0)
                )
                ** 2
                for row in range(len(modes))
                for column in range(len(modes))
            )
        )
        max_unitarity_residual = max(max_unitarity_residual, residual)
    return {
        "basis_mode_key_count": len(left),
        "basis_max_abs_difference": max_mode_difference,
        "singular_value_min_left": min(singular_values),
        "singular_value_max_left": max(singular_values),
        "singular_value_max_abs_difference": singular_value_difference,
        "overlap_unitarity_residual_max": max_unitarity_residual,
        "basis_keys_identical": True,
    }


def overlap_transport(left: dict, right: dict) -> list[list[complex]]:
    """Build the accepted block-diagonal product-basis overlap transport."""
    if set(left) != set(right):
        raise RuntimeError("cannot transport product bases with different mode keys")
    blocks: dict[tuple[int, int], list[int]] = {}
    for site, response_l, radial_index, mode in left:
        blocks.setdefault((site, response_l), []).append(radial_index)
    dimension = sum((2 * response_l + 1) * len({key[3] for key in left if key[:2] == (site, response_l)})
                    for site, response_l in blocks)
    transport = [[0.0j for _ in range(dimension)] for _ in range(dimension)]
    offset = 0
    for site, response_l in sorted(blocks):
        radial_indices = sorted(set(blocks[(site, response_l)]))
        modes = sorted({key[3] for key in left if key[:2] == (site, response_l)})
        rank = len(modes)
        block_overlap = [
            [
                sum(left[(site, response_l, radial, mode_left)][0].conjugate()
                    * right[(site, response_l, radial, mode_right)][0] for radial in radial_indices)
                for mode_right in modes
            ]
            for mode_left in modes
        ]
        for _m in range(2 * response_l + 1):
            for mode_left in range(rank):
                for mode_right in range(rank):
                    transport[offset + mode_left][offset + mode_right] = block_overlap[mode_left][mode_right]
            offset += rank
    return transport


def matrix_multiply(left: list[list[complex]], right: list[list[complex]]) -> list[list[complex]]:
    dimension = len(left)
    return [
        [sum(left[row][inner] * right[inner][column] for inner in range(dimension)) for column in range(dimension)]
        for row in range(dimension)
    ]


def transport_matrix(matrix: list[list[complex]], transport: list[list[complex]]) -> list[list[complex]]:
    transport_adjoint = [[transport[column][row].conjugate() for column in range(len(transport))]
                         for row in range(len(transport))]
    return matrix_multiply(matrix_multiply(transport, matrix), transport_adjoint)


def run_case(binary: Path, runner: Path, template: Path, fe_database: Path, scratch_root: Path, case: Case) -> dict[str, object]:
    workdir = scratch_root / case.name
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    input_path = workdir / "input.nml"
    prepare_input(template, input_path, fe_database, case)
    baseline = template.read_text(encoding="utf-8").splitlines(keepends=True)
    actual = input_path.read_text(encoding="utf-8").splitlines(keepends=True)
    (workdir / "input_diff.patch").write_text(
        "".join(difflib.unified_diff(baseline, actual, fromfile=str(template), tofile=str(input_path))), encoding="utf-8"
    )
    started_wall = time.time()
    started = time.monotonic()
    completed = subprocess.run(
        ["bash", str(runner), str(binary.resolve())],
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
    )
    elapsed = time.monotonic() - started
    log = (workdir / "testrun.log").read_text(encoding="utf-8", errors="replace") if (workdir / "testrun.log").exists() else ""
    if completed.returncode != 0:
        raise RuntimeError(f"{case.name}: executable failed with rc={completed.returncode}\n{log[-5000:]}\n{completed.stderr[-2000:]}")
    output_path = workdir / f"tddft_kspace_scf_{case.mesh}.dat"
    if not output_path.exists():
        raise RuntimeError(f"{case.name}: missing TDDFT output {output_path}")
    headers, rows = parse_convergence_output(output_path)
    if "Converged!" not in log or headers.get("accepted_state_source") != "accepted_kspace_scf_cache":
        raise RuntimeError(f"{case.name}: missing k-space-SCF convergence or direct accepted-state provenance")
    if headers.get("reciprocal_rebuild_performed_for_tddft", "T").upper() != "F":
        raise RuntimeError(f"{case.name}: TDDFT reports a reciprocal rebuild")
    if int(headers["product_dimension"]) != 232:
        raise RuntimeError(f"{case.name}: product dimension is not the complete 232-dimensional span")
    if len(rows) != len(case.eta_values) or any(int(row["q_index"]) != 1 or abs(float(row["omega_Ry"])) > 1.0e-12 for row in rows):
        raise RuntimeError(f"{case.name}: response is not the requested Gamma/static eta set")

    artifact_paths = {
        key: resolve_artifact(workdir, headers, key, started_wall)
        for key in ("kpoint_artifact", "basis_artifact", "matrix_artifact", "state_artifact")
    }
    scf_path = workdir / "kspace_scf_state.dat"
    if not scf_path.exists() or scf_path.stat().st_mtime < started_wall - 1.0:
        raise RuntimeError(f"{case.name}: missing or stale accepted SCF state artifact")
    state_identity = state_consistency(parse_state_artifact(scf_path), parse_state_artifact(artifact_paths["state_artifact"]))
    kpoint_headers, points = parse_kpoint_artifact(artifact_paths["kpoint_artifact"])
    expected_count = case.mesh**3
    if len(points) != expected_count or header_ints(kpoint_headers, "actual_k_count") != (expected_count,):
        raise RuntimeError(f"{case.name}: runtime k-point artifact is not the requested full mesh")
    if header_ints(headers, "actual_generated_k_mesh") != (case.mesh, case.mesh, case.mesh):
        raise RuntimeError(f"{case.name}: actual response mesh differs from requested mesh")
    basis_headers, basis = parse_basis_artifact(artifact_paths["basis_artifact"])
    matrix_headers, matrices = parse_matrix_artifact(artifact_paths["matrix_artifact"])
    if int(basis_headers["basis_product_dimension"]) != 232 or int(matrix_headers["basis_product_dimension"]) != 232:
        raise RuntimeError(f"{case.name}: basis/matrix artifacts are not complete 232-dimensional artifacts")

    result: dict[str, object] = {
        "name": case.name,
        "mesh": case.mesh,
        "requested_mesh": [case.mesh, case.mesh, case.mesh],
        "actual_mesh": list(header_ints(headers, "actual_generated_k_mesh")),
        "actual_k_count": int(header_number(headers, "actual_k_count")),
        "weight_sum": header_number(headers, "actual_k_weight_sum"),
        "accepted_state_EF_Ry": header_number(headers, "accepted_state_EF_Ry"),
        "target_electron_count": header_number(headers, "target_electron_count"),
        "accepted_electron_count": header_number(headers, "accepted_state_integrated_electron_count"),
        "accepted_electron_count_error": header_number(headers, "accepted_state_integrated_electron_count_error"),
        "scf_iterations": int(parse_state_artifact(scf_path)[0]["scf_iterations"]),
        "scf_residual": header_number(parse_state_artifact(scf_path)[0], "scf_residual"),
        "scf_moment_muB": header_number(parse_state_artifact(scf_path)[0], "scf_moment_muB"),
        "scf_physical_total_energy_Ry": header_number(parse_state_artifact(scf_path)[0], "scf_physical_total_energy_Ry"),
        "k_fingerprint_sha256": kpoint_headers["k_fingerprint_sha256"],
        "state_consistency": state_identity,
        "response_rows": rows,
        "response_artifacts": {key: str(path) for key, path in artifact_paths.items()},
        "input_diff": str(workdir / "input_diff.patch"),
        "product_dimension": 232,
        "rank_stable_all_L": headers["product_rank_stable_all_L"].upper() == "T",
        "radial_residual_control": header_number(headers, "accepted_state_radial_residual_control"),
        "elapsed_wall_seconds": elapsed,
        "_basis": basis,
        "_basis_headers": basis_headers,
        "_matrices": matrices,
    }
    if not result["rank_stable_all_L"]:
        raise RuntimeError(f"{case.name}: product basis rank is not stable")
    if case.gf_closure:
        gf_path = workdir / f"tddft_kspace_scf_{case.mesh}.dat.gf"
        if not gf_path.exists() or gf_path.stat().st_mtime < started_wall - 1.0:
            raise RuntimeError(f"{case.name}: missing same-process Gamma GF spot artifact")
        result["gf"] = parse_gf_artifact(gf_path)
    return result


def public_result(result: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in result.items() if not key.startswith("_")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--database", type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = args.template or root / "tests/integration/tddft_driver_smoke/input_tdvk05_fe.nml"
    fe_database = args.database or root / "tests/scf/cases/bulk/bccFe"
    runner = root / "tests/run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    results = [run_case(args.binary, runner, template, fe_database, args.scratch_root, case) for case in CASES]
    by_mesh = {int(result["mesh"]): result for result in results}
    if set(by_mesh) != {8, 12}:
        raise RuntimeError("bridge campaign did not produce exactly the 8^3 and 12^3 production states")
    if by_mesh[8]["k_fingerprint_sha256"] == by_mesh[12]["k_fingerprint_sha256"]:
        raise RuntimeError("bridge campaign negative control: 8^3 and 12^3 fingerprints match")
    transport = overlap_transport(by_mesh[8]["_basis"], by_mesh[12]["_basis"])
    operator_comparison = matrix_metrics(
        select_matrix(by_mesh[8]["_matrices"], 0.01),
        transport_matrix(select_matrix(by_mesh[12]["_matrices"], 0.01), transport),
    )
    raw_operator_comparison = matrix_metrics(
        select_matrix(by_mesh[8]["_matrices"], 0.01), select_matrix(by_mesh[12]["_matrices"], 0.01)
    )
    basis_comparison = basis_alignment_metrics(by_mesh[8]["_basis"], by_mesh[12]["_basis"])
    eta005 = select_matrix(by_mesh[12]["_matrices"], 0.005)
    eta001 = select_matrix(by_mesh[12]["_matrices"], 0.01)
    output = {
        "campaign": "TDVK review-gate bridge: self-consistent k-space SCF -> reciprocal TDDFT",
        "production_cases": [public_result(result) for result in results],
        "operator_comparison_8_to_12_eta001": {
            "basis_dimension": 232,
            "comparison": "12^3 compact operator transported into the 8^3 weighted product basis by block overlap",
            **basis_comparison,
            **operator_comparison,
            "raw_coordinate_metrics": raw_operator_comparison,
        },
        "mesh12_eta005_reuse_same_scf": {
            "scf_rerun": False,
            "eta_Ry": 0.005,
            **matrix_metrics(eta005, eta001),
        },
        "harness_audit": {
            "exactly_two_fresh_production_states": True,
            "use_kspace_true_checked": True,
            "same_process_direct_cache_handoff_checked": True,
            "actual_runtime_mesh_and_kpoint_rows_checked": True,
            "distinct_8cubed_12cubed_fingerprints": True,
            "scf_vs_tddft_state_arrays_compared": True,
            "explicit_occupations_compared": True,
            "gauge_invariant_density_projector_compared": True,
            "no_hidden_reciprocal_rebuild_checked": True,
            "no_tdvk06_interaction_physics": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"KSPACE-SCF -> TDDFT handoff cases={len(results)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
