#!/usr/bin/env python3
"""Run the TDVK-05 compact Fe bare-response convergence campaign.

This is an evidence campaign, not a default quick test.  Each case starts
from the same physical Fe input and launches its own isolated executable
process.  The response driver accepts the resulting real-space state, then
regenerates the requested reciprocal mesh and evaluates the compact Lehmann
service.  Only the selected 8^3 reference case performs the prescribed
finite-q reciprocal-GF spot checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


FLOAT = r"[-+0-9.EeDd]+"
HEADER_RE = re.compile(rf"^#\s*([^=]+?)\s*=\s*(.*?)\s*$")
METRIC_COLUMNS = (
    "eta_index",
    "eta_Ry",
    "q_index",
    "qx",
    "qy",
    "qz",
    "omega_Ry",
    "frobenius_norm",
    "max_abs_element",
    "trace_real",
    "trace_imag",
    "transitions_evaluated",
    "occupation_skips",
    "runtime_cpu_seconds",
    "finite",
)
GF_COLUMNS = (
    "q_index",
    "qx",
    "qy",
    "qz",
    "integration_points",
    "integration_eta",
    "h_over_integration_eta",
    "norm_lehmann",
    "norm_gf",
    "dF",
    "rF",
    "dInf",
    "wall_seconds",
    "finite",
)


@dataclass(frozen=True)
class Case:
    name: str
    mesh: int
    response_lmax: int
    eta_values: tuple[float, ...]
    backend: str = "product_convergence"


CASES = (
    Case("mesh4_full", 4, -1, (0.01,)),
    Case("mesh8_full_eta_ladder", 8, -1, (0.02, 0.01, 0.005)),
    Case("mesh12_full", 12, -1, (0.01,)),
    Case("mesh12_full_eta005_corner", 12, -1, (0.005,)),
    Case("mesh8_reduced_lmax2", 8, 2, (0.01,)),
)


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(name)}\s*=\s*).*$")
    text, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"could not patch namelist assignment {name}")
    return text


def prepare_input(template: Path, destination: Path, fe_database: Path, case: Case, gf: bool) -> None:
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(fe_database)))
    text = replace_assignment(text, "nk1", str(case.mesh))
    text = replace_assignment(text, "nk2", str(case.mesh))
    text = replace_assignment(text, "nk3", str(case.mesh))
    text = replace_assignment(text, "response_lmax", str(case.response_lmax))
    text = replace_assignment(text, "backend", repr(case.backend))
    text = replace_assignment(text, "n_eta", str(len(case.eta_values)))
    text = replace_assignment(text, "eta_grid", ", ".join(f"{eta:.12g}" for eta in case.eta_values))
    text = replace_assignment(text, "eta", f"{case.eta_values[0]:.12g}")
    if gf:
        text = replace_assignment(text, "n_q", "4")
        text = replace_assignment(text, "q_list", "0.0, 0.0, 0.0, 0.125, 0.0, 0.0, -0.125, 0.0, 0.0, 0.23, 0.07, -0.11")
        text = replace_assignment(text, "n_omega", "1")
        text = replace_assignment(text, "eta", "0.01")
        text = replace_assignment(text, "n_eta", "1")
        text = replace_assignment(text, "eta_grid", "0.01")
        text = replace_assignment(text, "gf_integration_points", "9601")
        text = replace_assignment(text, "backend", repr("product_finite_q"))
        text = replace_assignment(text, "output_file", repr("tddft_tdvk05_gf_spot.dat"))
    else:
        text = replace_assignment(text, "n_q", "2")
        text = replace_assignment(text, "q_list", "0.0, 0.0, 0.0, 0.125, 0.0, 0.0")
        text = replace_assignment(text, "n_omega", "2")
        text = replace_assignment(text, "output_file", repr("tddft_tdvk05_fe.dat"))
    destination.write_text(text, encoding="utf-8")


def parse_headers(lines: list[str]) -> dict[str, str]:
    headers: dict[str, str] = {}
    for line in lines:
        match = HEADER_RE.match(line)
        if match:
            headers[match.group(1).strip()] = match.group(2).strip()
    return headers


def parse_convergence_output(path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    rows: list[dict[str, object]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != len(METRIC_COLUMNS):
            raise RuntimeError(f"unexpected TDVK-05 metric row in {path}: {line}")
        row: dict[str, object] = {}
        for index, column in enumerate(METRIC_COLUMNS):
            if column in {"eta_index", "q_index", "transitions_evaluated", "occupation_skips"}:
                row[column] = int(fields[index])
            elif column == "finite":
                row[column] = fields[index].upper() == "T"
            else:
                row[column] = number(fields[index])
        rows.append(row)
    if not rows:
        raise RuntimeError(f"TDVK-05 output contains no metric rows: {path}")
    if not all(bool(row["finite"]) for row in rows):
        raise RuntimeError(f"TDVK-05 output contains a non-finite response: {path}")
    return headers, rows


def header_ints(headers: dict[str, str], key: str) -> tuple[int, ...]:
    try:
        return tuple(int(value) for value in headers[key].split())
    except (KeyError, ValueError) as exc:
        raise RuntimeError(f"missing or malformed integer header {key}") from exc


def header_number(headers: dict[str, str], key: str) -> float:
    try:
        return number(headers[key].split()[0])
    except (KeyError, IndexError, ValueError) as exc:
        raise RuntimeError(f"missing or malformed numeric header {key}") from exc


def parse_kpoint_artifact(path: Path) -> tuple[dict[str, str], list[tuple[int, float, float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = parse_headers(lines)
    points: list[tuple[int, float, float, float, float]] = []
    data_lines: list[str] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 5:
            raise RuntimeError(f"unexpected k-point artifact row in {path}: {line}")
        row = (int(fields[0]), *(number(value) for value in fields[1:]))
        points.append(row)
        data_lines.append(" ".join(fields))
    if not points:
        raise RuntimeError(f"k-point artifact contains no points: {path}")
    if [row[0] for row in points] != list(range(1, len(points) + 1)):
        raise RuntimeError(f"k-point artifact indices are not a complete ordered list: {path}")
    headers["k_fingerprint_sha256"] = hashlib.sha256("\n".join(data_lines).encode("ascii")).hexdigest()
    return headers, points


def parse_basis_artifact(path: Path) -> tuple[dict[str, str], dict[tuple[int, int, int, int], tuple[complex, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = parse_headers(lines)
    modes: dict[tuple[int, int, int, int], tuple[complex, float]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 7:
            raise RuntimeError(f"unexpected product-basis artifact row in {path}: {line}")
        site, response_l, radial_index, mode = (int(value) for value in fields[:4])
        modes[(site, response_l, radial_index, mode)] = (complex(number(fields[4]), number(fields[5])), number(fields[6]))
    if not modes:
        raise RuntimeError(f"product-basis artifact contains no modes: {path}")
    return headers, modes


def parse_matrix_artifact(path: Path) -> tuple[dict[str, str], dict[tuple[float, float], list[list[complex]]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = parse_headers(lines)
    grouped: dict[tuple[float, float], dict[tuple[int, int], complex]] = {}
    dimension = int(headers.get("basis_product_dimension", "0"))
    if dimension < 1:
        raise RuntimeError(f"matrix artifact has no valid basis dimension: {path}")
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 7:
            raise RuntimeError(f"unexpected compact matrix artifact row in {path}: {line}")
        eta_index = int(fields[0])
        eta = number(fields[1])
        omega = number(fields[2])
        row = int(fields[3])
        column = int(fields[4])
        if eta_index < 1 or not (1 <= row <= dimension and 1 <= column <= dimension):
            raise RuntimeError(f"invalid compact matrix index in {path}: {line}")
        grouped.setdefault((eta, omega), {})[(row, column)] = complex(number(fields[5]), number(fields[6]))
    matrices: dict[tuple[float, float], list[list[complex]]] = {}
    for key, values in grouped.items():
        expected = dimension * dimension
        if len(values) != expected:
            raise RuntimeError(f"incomplete compact matrix {key} in {path}: {len(values)}/{expected} entries")
        matrices[key] = [[values[(row, column)] for column in range(1, dimension + 1)] for row in range(1, dimension + 1)]
    if not matrices:
        raise RuntimeError(f"compact matrix artifact contains no matrices: {path}")
    return headers, matrices


def select_matrix(matrices: dict[tuple[float, float], list[list[complex]]], eta: float, omega: float = 0.0) -> list[list[complex]]:
    matches = [(key, matrix) for key, matrix in matrices.items() if abs(key[0] - eta) <= 1.0e-10 and abs(key[1] - omega) <= 1.0e-10]
    if len(matches) != 1:
        raise RuntimeError(f"expected one matrix at eta={eta}, omega={omega}, found {len(matches)}")
    return matches[0][1]


def matrix_metrics(left: list[list[complex]], right: list[list[complex]]) -> dict[str, float]:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise RuntimeError("cannot compare compact matrices with different dimensions")
    n = len(left)
    norm_left = sum(abs(value) ** 2 for row in left for value in row) ** 0.5
    norm_right = sum(abs(value) ** 2 for row in right for value in row) ** 0.5
    delta = [[left[row][column] - right[row][column] for column in range(n)] for row in range(n)]
    d_f = sum(abs(value) ** 2 for row in delta for value in row) ** 0.5
    trace_left = sum(left[index][index] for index in range(n))
    trace_right = sum(right[index][index] for index in range(n))
    trace_delta = trace_left - trace_right
    return {
        "norm_left": norm_left,
        "norm_right": norm_right,
        "dF": d_f,
        "relative_frobenius": d_f / max(norm_left, norm_right, 1.0e-300),
        "dInf": max(abs(value) for row in delta for value in row),
        "trace_left_real": trace_left.real,
        "trace_left_imag": trace_left.imag,
        "trace_right_real": trace_right.real,
        "trace_right_imag": trace_right.imag,
        "trace_difference_real": trace_delta.real,
        "trace_difference_imag": trace_delta.imag,
        "trace_difference_abs": abs(trace_delta),
    }


def basis_metrics(left: dict[tuple[int, int, int, int], tuple[complex, float]], right: dict[tuple[int, int, int, int], tuple[complex, float]]) -> dict[str, float | int | bool]:
    if set(left) != set(right):
        raise RuntimeError("compact basis mode keys differ between meshes")
    max_difference = max(abs(left[key][0] - right[key][0]) for key in left)
    singular_values = [value[1] for value in left.values()]
    if any(abs(value[1] - right[key][1]) > 1.0e-12 for key, value in left.items()):
        raise RuntimeError("compact basis singular values differ beyond serialization precision")
    blocks: dict[tuple[int, int], list[int]] = {}
    for site, response_l, radial_index, mode in left:
        blocks.setdefault((site, response_l), []).append(radial_index)
    max_unitarity_residual = 0.0
    for block, radial_indices in blocks.items():
        modes = sorted({key[3] for key in left if key[:2] == block})
        radial_indices = sorted(set(radial_indices))
        overlap = [[sum(left[(block[0], block[1], radial, mode_left)][0].conjugate() * right[(block[0], block[1], radial, mode_right)][0] for radial in radial_indices)
                    for mode_right in modes] for mode_left in modes]
        residual = sum(abs(sum(overlap[row][inner].conjugate() * overlap[inner][column] for inner in range(len(modes))) - (1.0 if row == column else 0.0)) ** 2
                       for row in range(len(modes)) for column in range(len(modes))) ** 0.5
        max_unitarity_residual = max(max_unitarity_residual, residual)
    return {
        "basis_mode_key_count": len(left),
        "basis_max_abs_difference": max_difference,
        "singular_value_min": min(singular_values),
        "singular_value_max": max(singular_values),
        "overlap_unitarity_residual_max": max_unitarity_residual,
        "basis_keys_identical": True,
    }


def parse_gf_output(path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    rows: list[dict[str, object]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) == len(GF_COLUMNS):
            row: dict[str, object] = {}
            for index, column in enumerate(GF_COLUMNS):
                if column in {"q_index", "integration_points"}:
                    row[column] = int(fields[index])
                elif column == "finite":
                    row[column] = fields[index].upper() == "T"
                else:
                    row[column] = number(fields[index])
            rows.append(row)
    if not rows:
        raise RuntimeError(f"TDVK-05 GF spot output contains no GF metric rows: {path}")
    if not all(bool(row["finite"]) for row in rows):
        raise RuntimeError(f"TDVK-05 GF spot output contains a non-finite response: {path}")
    return headers, rows


def run_case(binary: Path, runner: Path, template: Path, fe_database: Path, scratch_root: Path, case: Case, gf: bool) -> dict[str, object]:
    workdir = scratch_root / case.name
    if gf:
        workdir = scratch_root / "gf_spot_reference"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    input_path = workdir / "input.nml"
    prepare_input(template, input_path, fe_database, case, gf)
    started = time.monotonic()
    started_wall = time.time()
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
        raise RuntimeError(f"{case.name}: executable failed with rc={completed.returncode}\n{log[-4000:]}")
    output_name = "tddft_tdvk05_gf_spot.dat" if gf else "tddft_tdvk05_fe.dat"
    output_path = workdir / output_name
    if not output_path.exists():
        raise RuntimeError(f"{case.name}: missing {output_name}")
    if gf:
        headers, rows = parse_gf_output(output_path)
    else:
        headers, rows = parse_convergence_output(output_path)
    if "Converged!" not in log:
        raise RuntimeError(f"{case.name}: SCF did not report convergence")
    result: dict[str, object] = {
        "name": case.name,
        "mesh": case.mesh,
        "response_lmax_requested": case.response_lmax,
        "backend": "product_finite_q" if gf else case.backend,
        "elapsed_wall_seconds": elapsed,
        "headers": headers,
        "rows": rows,
    }
    if gf:
        actual_mesh = header_ints(headers, "accepted_state_nk")
        expected_count = case.mesh ** 3
        if actual_mesh != (expected_count,):
            raise RuntimeError(f"{case.name}: GF state has nk={actual_mesh}, expected {expected_count}")
        if not rows or any(int(row["integration_points"]) != 9601 for row in rows):
            raise RuntimeError(f"{case.name}: GF quadrature did not use the controlled 9601-point grid")
        result["actual_generated_k_mesh"] = [case.mesh, case.mesh, case.mesh]
        result["actual_k_count"] = expected_count
        result["quadrature_controlled"] = all(float(row["h_over_integration_eta"]) <= 0.5 for row in rows)
        if not bool(result["quadrature_controlled"]):
            raise RuntimeError(f"{case.name}: GF quadrature h/integration_eta exceeds 0.5")
        return result

    requested_mesh = (case.mesh, case.mesh, case.mesh)
    actual_mesh = header_ints(headers, "actual_generated_k_mesh")
    if header_ints(headers, "requested_k_mesh") != requested_mesh or actual_mesh != requested_mesh:
        raise RuntimeError(f"{case.name}: requested and actual mesh provenance do not match: {requested_mesh}/{actual_mesh}")
    actual_count = int(header_number(headers, "actual_k_count"))
    if actual_count != case.mesh ** 3 or int(headers.get("accepted_state_nk", "0")) != actual_count:
        raise RuntimeError(f"{case.name}: actual k-point count is not the requested full mesh")
    weight_sum = header_number(headers, "actual_k_weight_sum")
    if not (abs(weight_sum - 1.0) < 1.0e-10):
        raise RuntimeError(f"{case.name}: k-point weights do not normalize to one: {weight_sum}")

    artifact_paths: dict[str, Path] = {}
    artifact_data: dict[str, object] = {}
    for key in ("kpoint_artifact", "basis_artifact", "matrix_artifact"):
        artifact_path = Path(headers[key])
        if not artifact_path.is_absolute():
            artifact_path = workdir / artifact_path
        if not artifact_path.exists() or artifact_path.stat().st_mtime < started_wall - 1.0:
            raise RuntimeError(f"{case.name}: missing or stale runtime artifact {artifact_path}")
        artifact_paths[key] = artifact_path
    k_headers, points = parse_kpoint_artifact(artifact_paths["kpoint_artifact"])
    if len(points) != actual_count or header_ints(k_headers, "actual_k_count") != (actual_count,):
        raise RuntimeError(f"{case.name}: k-point artifact is not the actual complete mesh")
    basis_headers, basis = parse_basis_artifact(artifact_paths["basis_artifact"])
    matrix_headers, matrices = parse_matrix_artifact(artifact_paths["matrix_artifact"])
    if int(basis_headers["basis_product_dimension"]) != int(headers["product_dimension"]):
        raise RuntimeError(f"{case.name}: basis artifact dimension disagrees with response output")
    if int(matrix_headers["basis_product_dimension"]) != int(headers["product_dimension"]):
        raise RuntimeError(f"{case.name}: matrix artifact dimension disagrees with response output")
    result.update(
        {
            "requested_k_mesh": list(requested_mesh),
            "actual_generated_k_mesh": list(actual_mesh),
            "actual_k_count": actual_count,
            "actual_k_weight_sum": weight_sum,
            "k_fingerprint_sha256": k_headers["k_fingerprint_sha256"],
            "representative_k_vectors": {
                "first": [number(value) for value in headers["representative_k_first"].split()],
                "middle": [number(value) for value in headers["representative_k_middle"].split()],
                "last": [number(value) for value in headers["representative_k_last"].split()],
            },
            "eigenvalue_summary_Ry": {
                "min": header_number(headers, "eigenvalue_min_Ry"),
                "max": header_number(headers, "eigenvalue_max_Ry"),
                "mean": header_number(headers, "eigenvalue_mean_Ry"),
            },
            "fixed_EF_electron_count": header_number(headers, "fixed_EF_electron_count"),
            "target_electron_count": header_number(headers, "target_electron_count"),
            "fixed_EF_electron_count_error": header_number(headers, "fixed_EF_electron_count_error"),
            "diagnostic_mesh_EF_Ry": header_number(headers, "diagnostic_mesh_EF_Ry"),
            "diagnostic_mesh_EF_shift_Ry": header_number(headers, "diagnostic_mesh_EF_shift_Ry"),
            "product_dimension": int(headers["product_dimension"]),
            "product_unpruned_dimension": int(headers["product_unpruned_dimension"]),
            "rank_stable_all_L": headers["product_rank_stable_all_L"].upper() == "T",
            "artifacts": {key: str(path) for key, path in artifact_paths.items()},
            "_basis": basis,
            "_basis_headers": basis_headers,
            "_matrices": matrices,
        }
    )
    if not result["rank_stable_all_L"]:
        raise RuntimeError(f"{case.name}: product rank is not stable across the prescribed SVD tests")
    return result


def compare_operator_cases(left: dict[str, object], right: dict[str, object], eta: float = 0.01) -> dict[str, object]:
    left_headers = left["_basis_headers"]
    right_headers = right["_basis_headers"]
    left_dimension = int(left_headers["basis_product_dimension"])
    right_dimension = int(right_headers["basis_product_dimension"])
    if left_dimension != right_dimension:
        raise RuntimeError(f"operator comparison dimensions differ: {left_dimension}/{right_dimension}")
    basis = basis_metrics(left["_basis"], right["_basis"])
    matrices_left = left["_matrices"]
    matrices_right = right["_matrices"]
    metrics = matrix_metrics(select_matrix(matrices_left, eta), select_matrix(matrices_right, eta))
    return {
        "left_case": left["name"],
        "right_case": right["name"],
        "q": [0.0, 0.0, 0.0],
        "omega_Ry": 0.0,
        "eta_Ry": eta,
        "basis_dimension_left": left_dimension,
        "basis_dimension_right": right_dimension,
        "rank_provenance_left": left_headers["basis_rank_stable_all_L"],
        "rank_provenance_right": right_headers["basis_rank_stable_all_L"],
        **basis,
        **metrics,
    }


def compare_cross_corner(reference: dict[str, object], other: dict[str, object], reference_eta: float, other_eta: float) -> dict[str, object]:
    metrics = matrix_metrics(
        select_matrix(reference["_matrices"], reference_eta),
        select_matrix(other["_matrices"], other_eta),
    )
    return {
        "reference_case": reference["name"],
        "reference_eta_Ry": reference_eta,
        "other_case": other["name"],
        "other_eta_Ry": other_eta,
        "q": [0.0, 0.0, 0.0],
        "omega_Ry": 0.0,
        **metrics,
    }


def public_result(result: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in result.items() if not key.startswith("_")}


def audit_results(results: list[dict[str, object]]) -> dict[str, object]:
    convergence = [result for result in results if result["backend"] == "product_convergence"]
    by_mesh = {int(result["mesh"]): result for result in convergence if int(result["response_lmax_requested"]) == -1}
    required_meshes = (4, 8, 12)
    if set(by_mesh) != set(required_meshes):
        raise RuntimeError(f"missing complete-product mesh evidence: {sorted(by_mesh)}")
    fingerprints = {mesh: by_mesh[mesh]["k_fingerprint_sha256"] for mesh in required_meshes}
    if len(set(fingerprints.values())) != len(required_meshes):
        raise RuntimeError(f"negative control: distinct meshes have identical k fingerprints: {fingerprints}")
    for result in convergence:
        if result["actual_generated_k_mesh"] != [result["mesh"]] * 3:
            raise RuntimeError(f"negative control: requested mesh was not reached for {result['name']}")
    workdirs = [str(Path(result["artifacts"]["kpoint_artifact"]).parent) for result in convergence]
    if len(workdirs) != len(set(workdirs)):
        raise RuntimeError("negative control: convergence cases share a work directory")
    artifact_paths = [path for result in convergence for path in result["artifacts"].values()]
    if len(artifact_paths) != len(set(artifact_paths)):
        raise RuntimeError("negative control: runtime artifacts are reused across cases")
    return {
        "complete_meshes_compared": list(required_meshes),
        "distinct_k_fingerprints": fingerprints,
        "same_mesh_controls": {
            "mesh8_reduced_lmax2_matches_mesh8_k_fingerprint": next(
                result["k_fingerprint_sha256"] == fingerprints[8]
                for result in convergence
                if result["name"] == "mesh8_reduced_lmax2"
            ),
            "mesh12_eta005_matches_mesh12_k_fingerprint": next(
                result["k_fingerprint_sha256"] == fingerprints[12]
                for result in convergence
                if result["name"] == "mesh12_full_eta005_corner"
            ),
        },
        "requested_mesh_checked_against_actual_runtime": True,
        "actual_kpoint_artifact_rows_checked": True,
        "actual_eigenvalue_and_occupation_headers_checked": True,
        "fresh_artifacts_in_unique_workdirs": True,
        "no_cross_mesh_result_reuse": True,
        "no_baseline_scalar_substitution": True,
        "no_expected_fingerprint_generation": True,
        "cache_policy": "one fresh subprocess and fresh artifact set per case; no response cache is accepted without matching runtime provenance",
    }


def validate_gf_against_reference(gf: dict[str, object], reference: dict[str, object]) -> None:
    for key in ("accepted_state_EF_Ry", "accepted_state_temperature_K", "accepted_state_moment_muB"):
        if abs(header_number(gf["headers"], key) - header_number(reference["headers"], key)) > 1.0e-10:
            raise RuntimeError(f"GF state provenance differs from the 8^3 reference for {key}")
    for row in gf["rows"]:
        if abs(float(row["integration_eta"]) - 0.001) > 1.0e-12 or not math.isfinite(float(row["rF"])):
            raise RuntimeError("GF control has an unexpected integration or response parameter")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-gf", action="store_true", help="skip the single reference reciprocal-GF spot run")
    parser.add_argument("--gf-only", action="store_true", help="run only the single reference reciprocal-GF spot case")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = root / "tests/integration/tddft_driver_smoke/input_tdvk05_fe.nml"
    fe_database = root / "tests/scf/cases/bulk/bccFe"
    runner = root / "tests/run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    gf_case = Case("gf_spot_reference", 8, -1, (0.01,))
    if args.gf_only:
        gf_result = run_case(args.binary, runner, template, fe_database, args.scratch_root, gf_case, True)
        results: list[dict[str, object]] = [gf_result]
        output: dict[str, object] = {"cases": [public_result(gf_result)]}
    else:
        results = [run_case(args.binary, runner, template, fe_database, args.scratch_root, case, False) for case in CASES]
        convergence_audit = audit_results(results)
        complete = {result["name"]: result for result in results if result["response_lmax_requested"] == -1}
        operator_comparison = compare_operator_cases(complete["mesh8_full_eta_ladder"], complete["mesh12_full"])
        cross_corner = [
            compare_cross_corner(complete["mesh12_full_eta005_corner"], complete["mesh12_full"], 0.005, 0.01),
            compare_cross_corner(complete["mesh12_full_eta005_corner"], complete["mesh8_full_eta_ladder"], 0.005, 0.005),
            compare_cross_corner(complete["mesh12_full_eta005_corner"], complete["mesh8_full_eta_ladder"], 0.005, 0.01),
        ]
        output = {
            "cases": [public_result(result) for result in results],
            "operator_comparison": operator_comparison,
            "cross_corner": cross_corner,
            "harness_audit": convergence_audit,
        }
        if not args.no_gf:
            gf_result = run_case(args.binary, runner, template, fe_database, args.scratch_root, gf_case, True)
            validate_gf_against_reference(gf_result, complete["mesh8_full_eta_ladder"])
            output["gf"] = public_result(gf_result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"TDVK-05 cases={len(results)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
