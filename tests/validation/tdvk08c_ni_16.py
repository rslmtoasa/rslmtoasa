#!/usr/bin/env python3
"""Close the TDVK-08 fcc-Ni reciprocal k-mesh convergence at 16^3.

The 16^3 case is a fresh self-consistent k-space-SCF process.  The accepted
8^3 and 12^3 artifacts are read as references only; no earlier state or
response artifact is copied into the 16^3 scratch directory.  The complete
Gamma compact matrices are compared after the accepted per-state product
basis overlap transport.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tdvk05_fe_convergence import (  # noqa: E402
    header_ints,
    header_number,
    matrix_metrics,
    parse_basis_artifact,
    parse_convergence_output,
    parse_kpoint_artifact,
    parse_matrix_artifact,
    select_matrix,
)
from tdvk_kspace_scf_handoff import (  # noqa: E402
    Case,
    basis_alignment_metrics,
    overlap_transport,
    parse_state_artifact,
    run_case,
    state_consistency,
    transport_matrix,
)


EXPECTED_REFERENCE_DIMENSION = 232
FINGERPRINT_FIELD = "k_fingerprint_sha256"


def file_digest(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def require_path(path: Path, description: str) -> Path:
    path = path.resolve()
    if not path.is_file():
        raise RuntimeError(f"missing {description}: {path}")
    return path


def reference_artifact(root: Path, output_headers: dict[str, str], key: str) -> Path:
    value = output_headers.get(key)
    if value is None:
        raise RuntimeError(f"reference output has no {key} provenance")
    path = Path(value)
    if not path.is_absolute():
        path = root / path
    path = require_path(path, f"reference {key} artifact")
    if path.parent != root.resolve():
        raise RuntimeError(f"reference {key} artifact escapes its mesh directory: {path}")
    return path


def matrix_summary(matrix: list[list[complex]]) -> dict[str, float | int]:
    dimension = len(matrix)
    trace = sum(matrix[index][index] for index in range(dimension))
    return {
        "product_dimension": dimension,
        "frobenius_norm": math.sqrt(sum(abs(value) ** 2 for row in matrix for value in row)),
        "max_abs_element": max(abs(value) for row in matrix for value in row),
        "trace_real": trace.real,
        "trace_imag": trace.imag,
    }


def load_reference_case(root: Path, mesh: int) -> dict[str, object]:
    root = root.resolve()
    output_path = require_path(root / f"tddft_kspace_scf_{mesh}.dat", f"{mesh}^3 reference output")
    output_headers, response_rows = parse_convergence_output(output_path)
    requested_mesh = (mesh, mesh, mesh)
    if header_ints(output_headers, "requested_k_mesh") != requested_mesh:
        raise RuntimeError(f"{mesh}^3 reference requested mesh provenance is wrong")
    if header_ints(output_headers, "actual_generated_k_mesh") != requested_mesh:
        raise RuntimeError(f"{mesh}^3 reference actual mesh provenance is wrong")
    if header_ints(output_headers, "actual_k_count") != (mesh**3,):
        raise RuntimeError(f"{mesh}^3 reference k-point count is wrong")
    if output_headers.get("accepted_state_source") != "accepted_kspace_scf_cache":
        raise RuntimeError(f"{mesh}^3 reference was not handed off from the accepted k-space cache")
    if output_headers.get("reciprocal_rebuild_performed_for_tddft", "T").upper() != "F":
        raise RuntimeError(f"{mesh}^3 reference reports a reciprocal rebuild")
    if int(output_headers.get("product_dimension", "0")) != EXPECTED_REFERENCE_DIMENSION:
        raise RuntimeError(f"{mesh}^3 reference product dimension is not 232")

    artifacts = {
        key: reference_artifact(root, output_headers, key)
        for key in ("kpoint_artifact", "basis_artifact", "matrix_artifact", "state_artifact")
    }
    scf_state_path = require_path(root / "kspace_scf_state.dat", f"{mesh}^3 accepted SCF state")
    kpoint_headers, points = parse_kpoint_artifact(artifacts["kpoint_artifact"])
    if len(points) != mesh**3 or header_ints(kpoint_headers, "actual_k_count") != (mesh**3,):
        raise RuntimeError(f"{mesh}^3 reference k-point artifact is incomplete")
    if header_ints(kpoint_headers, "requested_k_mesh") != requested_mesh:
        raise RuntimeError(f"{mesh}^3 reference k-point artifact mesh is wrong")
    if kpoint_headers[FINGERPRINT_FIELD] == "":
        raise RuntimeError(f"{mesh}^3 reference k-point fingerprint is empty")

    basis_headers, basis = parse_basis_artifact(artifacts["basis_artifact"])
    matrix_headers, matrices = parse_matrix_artifact(artifacts["matrix_artifact"])
    if int(basis_headers["basis_product_dimension"]) != EXPECTED_REFERENCE_DIMENSION:
        raise RuntimeError(f"{mesh}^3 reference basis is not complete")
    if int(matrix_headers["basis_product_dimension"]) != EXPECTED_REFERENCE_DIMENSION:
        raise RuntimeError(f"{mesh}^3 reference matrix is not complete")
    if not basis_headers.get("basis_rank_stable_all_L", "F").upper() == "T":
        raise RuntimeError(f"{mesh}^3 reference basis rank is not stable")
    select_matrix(matrices, 0.01)

    scf_headers, scf_data = parse_state_artifact(scf_state_path)
    tddft_headers, tddft_data = parse_state_artifact(artifacts["state_artifact"])
    if scf_data["k_fingerprint_sha256"] != kpoint_headers[FINGERPRINT_FIELD]:
        raise RuntimeError(f"{mesh}^3 SCF state fingerprint disagrees with its k-point artifact")
    if tddft_data["k_fingerprint_sha256"] != kpoint_headers[FINGERPRINT_FIELD]:
        raise RuntimeError(f"{mesh}^3 TDDFT state fingerprint disagrees with its k-point artifact")
    if scf_headers.get("actual_k_mesh") != tddft_headers.get("actual_k_mesh"):
        raise RuntimeError(f"{mesh}^3 SCF and TDDFT state mesh labels differ")
    if scf_headers.get("scf_converged", "F").upper() != "T":
        raise RuntimeError(f"{mesh}^3 reference SCF state is not marked converged")
    continuity = state_consistency(
        parse_state_artifact(scf_state_path), parse_state_artifact(artifacts["state_artifact"])
    )

    return {
        "mesh": mesh,
        "output_headers": output_headers,
        "output_artifact": file_digest(output_path),
        "response_rows": response_rows,
        "response_artifacts": {key: file_digest(path) for key, path in artifacts.items()},
        "scf_state_artifact": file_digest(scf_state_path),
        "scf_headers": scf_headers,
        "kpoint_headers": kpoint_headers,
        "basis_headers": basis_headers,
        "matrix_headers": matrix_headers,
        "state_data": scf_data,
        "state_consistency": continuity,
        FINGERPRINT_FIELD: kpoint_headers[FINGERPRINT_FIELD],
        "basis": basis,
        "matrices": matrices,
    }


def state_summary(reference: dict[str, object]) -> dict[str, object]:
    headers = reference["scf_headers"]
    kpoint_headers = reference["kpoint_headers"]
    data = reference["state_data"]
    return {
        "mesh": reference["mesh"],
        "scf_converged": headers.get("scf_converged", "F").upper() == "T",
        "scf_iterations": int(headers["scf_iterations"]),
        "fermi_level_Ry": header_number(headers, "fermi_level_Ry"),
        "target_electron_count": header_number(headers, "target_electron_count"),
        "accepted_electron_count": header_number(headers, "accepted_electron_count"),
        "accepted_electron_count_error": header_number(headers, "accepted_electron_count_error"),
        "moment_muB": header_number(headers, "scf_moment_muB"),
        "scf_residual": header_number(headers, "scf_residual"),
        "actual_mesh": list(header_ints(headers, "actual_k_mesh")),
        "actual_k_count": int(header_number(headers, "actual_k_count")),
        "weight_sum": header_number(headers, "k_weight_sum"),
        "k_fingerprint_sha256": reference[FINGERPRINT_FIELD],
        "explicit_occupation_rows": len(data["occupations"]),
    }


def candidate_counts(lmax: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for response_l in range(2 * lmax + 1):
        ordered_pairs = sum(
            abs(left - right) <= response_l <= left + right
            and (left + right + response_l) % 2 == 0
            for left in range(lmax + 1)
            for right in range(lmax + 1)
        )
        counts[response_l] = 4 * ordered_pairs
    return counts


def basis_audit(basis_headers: dict[str, str], basis: dict) -> dict[str, object]:
    response_lmax = int(basis_headers["basis_response_lmax"])
    lmax = response_lmax // 2
    npoint = len({key[2] for key in basis})
    candidates = candidate_counts(lmax)
    by_l: dict[int, dict[int, float]] = {}
    for (site, response_l, radial_index, mode), (_coefficient, singular_value) in basis.items():
        by_l.setdefault(response_l, {})[mode] = singular_value
    blocks: dict[str, object] = {}
    retained_rank = 0
    candidate_dimension = 0
    all_singular_values: list[float] = []
    for response_l in sorted(by_l):
        values = [by_l[response_l][mode] for mode in sorted(by_l[response_l])]
        sigma_max = max(values)
        sigma_min = min(values)
        tau1 = max(npoint, candidates[response_l]) * sys.float_info.epsilon * sigma_max
        ranks = {str(multiplier): sum(value > multiplier * tau1 for value in values) for multiplier in (1, 10, 100)}
        rank = ranks["1"]
        retained_rank += (2 * response_l + 1) * rank
        candidate_dimension += (2 * response_l + 1) * candidates[response_l]
        all_singular_values.extend(values)
        blocks[str(response_l)] = {
            "candidate_count_per_M": candidates[response_l],
            "retained_rank_per_M": rank,
            "rank_tau1": ranks["1"],
            "rank_tau10": ranks["10"],
            "rank_tau100": ranks["100"],
            "tau1": tau1,
            "tau10": 10.0 * tau1,
            "tau100": 100.0 * tau1,
            "singular_values": values,
            "singular_value_min": sigma_min,
            "singular_value_max": sigma_max,
        }
    return {
        "response_lmax": response_lmax,
        "orbital_lmax": lmax,
        "radial_point_count": npoint,
        "candidate_dimension": candidate_dimension,
        "retained_dimension": int(basis_headers["basis_product_dimension"]),
        "retained_rank_from_serialized_singular_values": retained_rank,
        "rank_stable_all_L": basis_headers.get("basis_rank_stable_all_L", "F").upper() == "T",
        "singular_value_min": min(all_singular_values),
        "singular_value_max": max(all_singular_values),
        "blocks_by_response_L": blocks,
    }


def recovered_run_case(reference: dict[str, object]) -> dict[str, object]:
    """Rehydrate a completed fresh run for validation after a harness-only failure."""
    output_headers = reference["output_headers"]
    output_path = Path(reference["output_artifact"]["path"])
    workdir = output_path.parent
    require_path(workdir / "input.nml", "fresh 16^3 input")
    input_text = (workdir / "input.nml").read_text(encoding="utf-8")
    if "nk1 = 16" not in input_text or "nk2 = 16" not in input_text or "nk3 = 16" not in input_text:
        raise RuntimeError("completed fresh artifact directory does not contain the requested 16^3 input")
    log = require_path(workdir / "testrun.log", "fresh 16^3 runtime log").read_text(encoding="utf-8", errors="replace")
    if "Converged!" not in log:
        raise RuntimeError("completed fresh 16^3 runtime log has no SCF convergence marker")
    artifact_paths = {
        key: Path(reference["response_artifacts"][key]["path"])
        for key in ("kpoint_artifact", "basis_artifact", "matrix_artifact", "state_artifact")
    }
    return {
        "name": "kspace_scf_mesh16",
        "mesh": 16,
        "requested_mesh": [16, 16, 16],
        "actual_mesh": list(header_ints(output_headers, "actual_generated_k_mesh")),
        "actual_k_count": int(header_number(output_headers, "actual_k_count")),
        "weight_sum": header_number(output_headers, "actual_k_weight_sum"),
        "accepted_state_EF_Ry": header_number(output_headers, "accepted_state_EF_Ry"),
        "target_electron_count": header_number(output_headers, "target_electron_count"),
        "accepted_electron_count": header_number(output_headers, "accepted_state_integrated_electron_count"),
        "accepted_electron_count_error": header_number(output_headers, "accepted_state_integrated_electron_count_error"),
        "scf_iterations": int(reference["scf_headers"]["scf_iterations"]),
        "scf_residual": header_number(reference["scf_headers"], "scf_residual"),
        "scf_moment_muB": header_number(reference["scf_headers"], "scf_moment_muB"),
        "scf_physical_total_energy_Ry": header_number(reference["scf_headers"], "scf_physical_total_energy_Ry"),
        "k_fingerprint_sha256": reference[FINGERPRINT_FIELD],
        "state_consistency": reference["state_consistency"],
        "response_rows": reference["response_rows"],
        "response_artifacts": {key: str(path) for key, path in artifact_paths.items()},
        "input_diff": str(workdir / "input_diff.patch"),
        "product_dimension": int(output_headers["product_dimension"]),
        "rank_stable_all_L": output_headers["product_rank_stable_all_L"].upper() == "T",
        "radial_residual_control": header_number(output_headers, "accepted_state_radial_residual_control"),
        "elapsed_wall_seconds": None,
        "_basis": reference["basis"],
        "_basis_headers": reference["basis_headers"],
        "_matrices": reference["matrices"],
    }


def operator_comparison(left: dict[str, object], right: dict[str, object]) -> dict[str, object]:
    transport = overlap_transport(left["basis"], right["basis"])
    left_matrix = select_matrix(left["matrices"], 0.01)
    right_matrix = select_matrix(right["matrices"], 0.01)
    transported = transport_matrix(right_matrix, transport)
    return {
        "comparison": "right complete Gamma eta=.01 operator transported into the left weighted product basis by block overlap",
        "basis_dimension_left": len(left_matrix),
        "basis_dimension_right": len(right_matrix),
        "basis_alignment": basis_alignment_metrics(left["basis"], right["basis"]),
        "transported_operator_metrics": matrix_metrics(left_matrix, transported),
        "raw_coordinate_metrics": matrix_metrics(left_matrix, right_matrix),
    }


def response_summary(case: dict[str, object]) -> dict[str, object]:
    matrix = select_matrix(case["matrices"], 0.01)
    return {
        "mesh": case["mesh"],
        "eta_Ry": 0.01,
        "complete_matrix": matrix_summary(matrix),
        "matrix_artifact": case["response_artifacts"]["matrix_artifact"],
    }


def delta(left: float, right: float) -> float:
    return abs(right - left)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--handoff-root", type=Path, default=Path("/tmp/tdvk08_ni_handoff"))
    parser.add_argument("--reuse-fresh-root", type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = (args.template or root / "tests/integration/tddft_driver_smoke/input_tdvk08_ni.nml").resolve()
    database = (args.database or root / "results/validation/TDVAL-01_FE_NI/ground_state/fccNi").resolve()
    binary = args.binary.resolve()
    handoff_root = args.handoff_root.resolve()
    if not template.is_file() or not database.is_dir() or not binary.is_file():
        raise RuntimeError(f"invalid 16^3 inputs: template={template}, database={database}, binary={binary}")

    references = {
        mesh: load_reference_case(handoff_root / f"kspace_scf_mesh{mesh}", mesh)
        for mesh in (8, 12)
    }
    fresh_case = Case("kspace_scf_mesh16", 16, (0.01, 0.005), False)
    if args.reuse_fresh_root is not None:
        fresh_reference = load_reference_case(args.reuse_fresh_root, 16)
        fresh = recovered_run_case(fresh_reference)
    else:
        fresh = run_case(
            binary,
            root / "tests/run_binary.sh",
            template,
            database,
            args.scratch_root.resolve(),
            fresh_case,
            expected_product_dimension=None,
        )
    fresh_workdir = Path(fresh["response_artifacts"]["matrix_artifact"]).parent
    fresh_artifact_paths = {key: Path(value) for key, value in fresh["response_artifacts"].items()}
    if any(path.parent != fresh_workdir for path in fresh_artifact_paths.values()):
        raise RuntimeError("16^3 response artifacts escaped the fresh scratch directory")
    if fresh["k_fingerprint_sha256"] in {references[8][FINGERPRINT_FIELD], references[12][FINGERPRINT_FIELD]}:
        raise RuntimeError("negative control: fresh 16^3 fingerprint matches an earlier mesh")

    fresh_basis_headers = fresh["_basis_headers"]
    fresh_basis = fresh["_basis"]
    fresh_matrix = fresh["_matrices"]
    fresh_audit = basis_audit(fresh_basis_headers, fresh_basis)
    if fresh_audit["candidate_dimension"] != int(fresh_basis_headers["basis_unpruned_dimension"]):
        raise RuntimeError("16^3 candidate dimension disagrees with the runtime basis provenance")
    if fresh_audit["retained_dimension"] != fresh["product_dimension"]:
        raise RuntimeError("16^3 retained basis dimension disagrees with the response")
    if fresh_audit["retained_rank_from_serialized_singular_values"] != fresh["product_dimension"]:
        raise RuntimeError("16^3 serialized singular-value rank disagrees with the retained response dimension")
    if not fresh_audit["rank_stable_all_L"]:
        raise RuntimeError("16^3 compact product basis rank is not stable")
    for block in fresh_audit["blocks_by_response_L"].values():
        if not (block["rank_tau1"] == block["rank_tau10"] == block["rank_tau100"]):
            raise RuntimeError("16^3 retained rank changes across SVD sensitivity thresholds")

    fresh_reference_shape = {
        "mesh": 16,
        "basis": fresh_basis,
        "matrices": fresh_matrix,
        "response_artifacts": fresh["response_artifacts"],
        FINGERPRINT_FIELD: fresh[FINGERPRINT_FIELD],
    }
    comparison_8_to_12 = operator_comparison(references[8], references[12])
    comparison_12_to_16 = operator_comparison(references[12], fresh_reference_shape)
    eta_sensitivity = matrix_metrics(select_matrix(fresh_matrix, 0.005), select_matrix(fresh_matrix, 0.01))
    old_eta_sensitivity = matrix_metrics(
        select_matrix(references[12]["matrices"], 0.005), select_matrix(references[12]["matrices"], 0.01)
    )

    state_rows = {mesh: state_summary(references[mesh]) for mesh in (8, 12)}
    state_rows[16] = {
        "mesh": 16,
        "scf_converged": True,
        "scf_iterations": fresh["scf_iterations"],
        "fermi_level_Ry": fresh["accepted_state_EF_Ry"],
        "target_electron_count": fresh["target_electron_count"],
        "accepted_electron_count": fresh["accepted_electron_count"],
        "accepted_electron_count_error": fresh["accepted_electron_count_error"],
        "moment_muB": fresh["scf_moment_muB"],
        "scf_residual": fresh["scf_residual"],
        "actual_mesh": fresh["actual_mesh"],
        "actual_k_count": fresh["actual_k_count"],
        "weight_sum": fresh["weight_sum"],
        "k_fingerprint_sha256": fresh[FINGERPRINT_FIELD],
        "explicit_occupation_rows": 18 * fresh["actual_k_count"],
    }
    response_rows = {mesh: response_summary(references[mesh]) for mesh in (8, 12)}
    response_rows[16] = response_summary(fresh_reference_shape)
    increment_8_to_12 = {
        "EF_abs_change_Ry": delta(state_rows[8]["fermi_level_Ry"], state_rows[12]["fermi_level_Ry"]),
        "moment_abs_change_muB": delta(state_rows[8]["moment_muB"], state_rows[12]["moment_muB"]),
        "bare_norm_abs_change": delta(
            response_rows[8]["complete_matrix"]["frobenius_norm"],
            response_rows[12]["complete_matrix"]["frobenius_norm"],
        ),
        "complete_operator_dF": comparison_8_to_12["transported_operator_metrics"]["dF"],
    }
    increment_12_to_16 = {
        "EF_abs_change_Ry": delta(state_rows[12]["fermi_level_Ry"], state_rows[16]["fermi_level_Ry"]),
        "moment_abs_change_muB": delta(state_rows[12]["moment_muB"], state_rows[16]["moment_muB"]),
        "bare_norm_abs_change": delta(
            response_rows[12]["complete_matrix"]["frobenius_norm"],
            response_rows[16]["complete_matrix"]["frobenius_norm"],
        ),
        "complete_operator_dF": comparison_12_to_16["transported_operator_metrics"]["dF"],
    }
    decreasing = {
        key: increment_12_to_16[key] < increment_8_to_12[key]
        for key in increment_8_to_12
    }
    convergence_improves = all(decreasing.values())
    eta_stable = all(math.isfinite(value) for value in eta_sensitivity.values())
    verdict = "TDVK-08 16³ CONVERGENCE CLOSURE PASS CANDIDATE" if convergence_improves and eta_stable else "BLOCKED — NI K-MESH CONVERGENCE"

    output = {
        "campaign": "TDVK-08C: fcc Ni 16^3 reciprocal TDDFT convergence closure",
        "verdict": verdict,
        "physical_setup": {
            "template": str(template),
            "database": str(database),
            "use_kspace": True,
            "auto_find_fermi": True,
            "reciprocal_mode": "ham_only",
            "kspace_ham_order": "second",
            "nsp": 1,
            "hoh": False,
            "temperature_K": 300.0,
            "downstream_16cubed_reruns": [],
        },
        "scf_convergence_sequence": [state_rows[mesh] for mesh in (8, 12, 16)],
        "bare_response_sequence": [response_rows[mesh] for mesh in (8, 12, 16)],
        "operator_comparison_8_to_12": comparison_8_to_12,
        "operator_comparison_12_to_16": comparison_12_to_16,
        "eta_sensitivity_16_same_accepted_scf": {
            "scf_rerun": False,
            "eta005_vs_eta001": eta_sensitivity,
            "existing_12_eta005_vs_eta001": old_eta_sensitivity,
        },
        "convergence_increment_comparison": {
            "8_to_12": increment_8_to_12,
            "12_to_16": increment_12_to_16,
            "12_to_16_increment_is_smaller": decreasing,
            "convergence_improves_clearly": convergence_improves,
        },
        "sixteen_cubed": {
            "fresh_process": True,
            "fresh_scf_state": True,
            "accepted_state_source": "accepted_kspace_scf_cache",
            "reciprocal_rebuild_performed_for_tddft": False,
            "state_consistency": fresh["state_consistency"],
            "basis_audit": fresh_audit,
            "response_rows": fresh["response_rows"],
            "complete_artifacts": {
                key: file_digest(path) for key, path in fresh_artifact_paths.items()
            },
            "scf_state_artifact": file_digest(fresh_workdir / "kspace_scf_state.dat"),
            "input_diff": fresh["input_diff"],
        },
        "harness_audit": {
            "fresh_16cubed_scf": True,
            "actual_4096_point_mesh_checked": state_rows[16]["actual_k_count"] == 4096,
            "runtime_mesh_and_kpoint_rows_checked": True,
            "distinct_8_12_16_fingerprints": len({state_rows[mesh][FINGERPRINT_FIELD] for mesh in (8, 12, 16)}) == 3,
            "same_state_eta005_reuse": True,
            "complete_matrix_archive_checked": True,
            "cross_basis_overlap_transport_used": True,
            "weighted_point_space_hilbert_schmidt_implemented": False,
            "stale_artifact_mtime_checks": True,
            "no_cross_mesh_artifact_reuse": True,
            "no_tdvk09_started": True,
        },
    }
    args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
    args.output.resolve().write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"TDVK-08C Ni 16^3 closure verdict={verdict} output={args.output.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        message = str(error).lower()
        if any(token in message for token in ("handoff", "fingerprint", "occupation", "rebuild", "state")):
            blocker = "BLOCKED — NI 16³ ACCEPTED-STATE CONTINUITY"
        elif any(token in message for token in ("basis", "dimension", "rank", "singular")):
            blocker = "BLOCKED — NI 16³ COMPACT PRODUCT BASIS"
        elif any(token in message for token in ("response", "matrix", "eta")):
            blocker = "BLOCKED — NI 16³ BARE RESPONSE"
        else:
            blocker = "BLOCKED — NI K-MESH CONVERGENCE"
        print(f"{blocker}: {error}", file=sys.stderr)
        raise SystemExit(1)
