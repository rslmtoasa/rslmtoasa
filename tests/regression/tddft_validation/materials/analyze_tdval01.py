#!/usr/bin/env python3
"""Check the physical Fe/Ni TDVAL-01 fixtures.

This checker deliberately makes a material claim harder than an executable
smoke test.  It checks the physical periodic setup, q=0 Goldstone residual,
resolved finite-q collective signals, their energy window, and an origin-
constrained omega versus |q_cartesian|^2 fit.  A Stoner feature or an
unresolved Lorentzian is reported as diagnostic evidence, not as a magnon.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any


def parse_float(token: str) -> float:
    normalized = token.replace("D", "E").replace("d", "e")
    try:
        return float(normalized)
    except ValueError:
        match = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)", normalized)
        if not match:
            raise
        return float(f"{match.group(1)}e{match.group(2)}")


def metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        record = line.strip()
        if not record.startswith("#") or "=" not in record:
            continue
        key, value = record[1:].split("=", 1)
        values[key.strip()] = value.strip()
    return values


def vector(values: str) -> tuple[float, float, float]:
    fields = values.split()
    if len(fields) < 3:
        raise ValueError(f"expected three vector components, got {values!r}")
    return tuple(parse_float(value) for value in fields[:3])  # type: ignore[return-value]


def norm(values: tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def read_modes(path: Path) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "candidate" and len(fields) >= 7:
            index = int(fields[1])
            records.setdefault(index, {})["candidate"] = {
                "grid_omega": parse_float(fields[2]),
                "mode": int(fields[3]),
                "unity_distance": parse_float(fields[4]),
                "branch_overlap": parse_float(fields[5]),
            }
        elif fields[0] == "crossing" and len(fields) >= 10:
            index = int(fields[1])
            records.setdefault(index, {})["crossing"] = {
                "present": fields[2].upper().startswith("T"),
                "omega": parse_float(fields[3]),
                "imaginary_part": parse_float(fields[4]),
                "branch_overlap": parse_float(fields[5]),
                "eigenvalue_step": parse_float(fields[6]),
                "projected_weight": parse_float(fields[7]),
                "condition_number": parse_float(fields[8]),
                "exceptional_warning": fields[9].upper().startswith("T"),
            }
        elif fields[0] == "fit" and len(fields) >= 8:
            index = int(fields[1])
            records.setdefault(index, {})["fit"] = {
                "accepted": fields[2].upper().startswith("T"),
                "center": parse_float(fields[3]),
                "fwhm": parse_float(fields[4]),
                "hwhm": parse_float(fields[5]),
                "relative_residual": parse_float(fields[6]),
                "reason": " ".join(fields[7:]),
            }
    return records


def source_contract(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    errors: list[str] = []
    for key in ("pbc", "b1", "b2", "b3"):
        if not re.search(rf"(?im)^\s*{key}\s*=\s*\.?true\.?\s*$", text):
            errors.append(f"{key} is not explicitly true")
    if not re.search(r"(?im)^\s*strux_backend\s*=\s*['\"]strux_lib['\"]\s*$", text):
        errors.append("strux_backend is not explicitly 'strux_lib'")
    return errors


def q_points(run_dir: Path) -> list[dict[str, Any]]:
    q_file = run_dir / "q_points.dat"
    if not q_file.is_file():
        raise ValueError(f"{run_dir}: q_points.dat is missing")
    q_lines = [line.split() for line in q_file.read_text().splitlines() if line.strip()]
    if not q_lines or len(q_lines[0]) != 1:
        raise ValueError(f"{q_file}: malformed q-point count")
    count = int(q_lines[0][0])
    if len(q_lines) != count + 1 or any(len(line) != 3 for line in q_lines[1:]):
        raise ValueError(f"{q_file}: q-point count and rows disagree")
    points: list[dict[str, Any]] = []
    for index in range(1, count + 1):
        candidates = sorted(
            path
            for path in run_dir.glob(f"*_q{index:06d}_chi0.dat")
            if "minus_plus" not in path.name
        )
        if len(candidates) != 1:
            raise ValueError(f"{run_dir}: expected one primary chi0 file for q index {index}")
        path = candidates[0]
        data = metadata(path)
        direct = vector(data["q_direct"])
        requested = tuple(parse_float(value) for value in q_lines[index])
        if max(abs(left - right) for left, right in zip(direct, requested)) > 1.0e-12:
            raise ValueError(f"{path}: emitted q_direct does not match q_points.dat")
        cartesian = vector(data["q_cartesian"])
        points.append(
            {
                "path": str(path),
                "q_direct": direct,
                "q_cartesian": cartesian,
                "q_norm": norm(cartesian),
                "eta_Ry": parse_float(data["eta_Ry"]),
                "omega_min_Ry": parse_float(data["omega_min_Ry"]),
                "omega_max_Ry": parse_float(data["omega_max_Ry"]),
                "omega_count": int(data["omega_batch_size"]),
                "k_mesh": data.get("k_mesh_shape"),
                "backend": data.get("chi0_backend_canonical"),
            }
        )
    return points


def goldstone(run_dir: Path) -> dict[str, Any]:
    files = sorted(run_dir.glob("*_goldstone.dat"))
    files = [path for path in files if "minus_plus" not in path.name]
    if not files:
        raise ValueError(f"{run_dir}: q=0 Goldstone output is missing")
    data = metadata(files[0])
    return {
        "path": str(files[0]),
        "raw_r_Xi": parse_float(data["raw_r_Xi"]),
        "closest_eigenvalue": parse_float(data["raw_closest_eigenvalue"].split()[0]),
        "signed_moment": parse_float(data["raw_magnetization_norm"]),
    }


def fit_origin(points: list[tuple[float, float]]) -> dict[str, float]:
    denominator = sum(q2 * q2 for q2, _ in points)
    stiffness = sum(q2 * omega for q2, omega in points) / denominator
    scale = max(max(abs(omega) for _, omega in points), 1.0e-30)
    residual = math.sqrt(sum((omega - stiffness * q2) ** 2 for q2, omega in points) / len(points)) / scale
    return {"D_Ry_A2": stiffness, "relative_residual": residual, "points": float(len(points))}


def analyse_run(record: dict[str, Any], root: Path, goldstone_limit: float, dispersion_limit: float,
                minimum_points: int) -> dict[str, Any]:
    run_dir = Path(record["directory"])
    if not run_dir.is_absolute():
        run_dir = root / run_dir
    result: dict[str, Any] = {
        "material": record["material"],
        "mesh": record["mesh"],
        "q_set": record["q_set"],
        "directory": str(run_dir),
        "status": "PASS",
        "reasons": [],
    }
    if record.get("status") != "PASS":
        result["status"] = "FAIL"
        result["reasons"].append(f"executable returned {record.get('returncode')}")
        return result

    input_errors = source_contract(run_dir / "source_input.nml")
    if input_errors:
        result["status"] = "FAIL"
        result["reasons"].extend(input_errors)
    try:
        points = q_points(run_dir)
        modes_files = sorted(path for path in run_dir.glob("*_pair_modes.dat") if "minus_plus" not in path.name)
        if len(modes_files) != 1:
            raise ValueError(f"expected one primary pair-mode file, found {len(modes_files)}")
        modes = read_modes(modes_files[0])
        zero = goldstone(run_dir)
    except (KeyError, OSError, ValueError) as error:
        result["status"] = "FAIL"
        result["reasons"].append(str(error))
        return result

    result["goldstone"] = zero
    zero_pass = zero["raw_r_Xi"] <= goldstone_limit and abs(zero["closest_eigenvalue"] - 1.0) <= 10.0 * goldstone_limit
    if not zero_pass:
        result["status"] = "FAIL"
        result["reasons"].append("q=0 pair-Xi Goldstone residual/eigenvalue is outside tolerance")
    backends = {point["backend"] for point in points}
    if backends != {"eigenpairs"}:
        result["status"] = "FAIL"
        result["reasons"].append(f"unexpected response backend metadata: {sorted(backends)}")

    signals: list[dict[str, Any]] = []
    resolved: list[tuple[float, float]] = []
    point_reports: list[dict[str, Any]] = []
    for index, point in enumerate(points, start=1):
        mode = modes.get(index, {})
        crossing = mode.get("crossing", {})
        fit = mode.get("fit", {})
        omega_step = (point["omega_max_Ry"] - point["omega_min_Ry"]) / max(point["omega_count"] - 1, 1)
        has_signal = (
            crossing.get("present", False)
            and not crossing.get("exceptional_warning", True)
            and abs(crossing.get("imaginary_part", math.inf)) <= 0.05
            and crossing.get("projected_weight", 0.0) > 0.0
            and crossing.get("omega", -1.0) > max(5.0 * point["eta_Ry"], 3.0 * omega_step)
            and crossing.get("omega", math.inf) < point["omega_max_Ry"] - 3.0 * omega_step
        )
        report = {
            "index": index,
            "q_direct": point["q_direct"],
            "q_cartesian_A-1": point["q_cartesian"],
            "q_norm_A-1": point["q_norm"],
            "crossing_omega_Ry": crossing.get("omega"),
            "fit_center_Ry": fit.get("center"),
            "fit_accepted": fit.get("accepted", False),
            "fit_reason": fit.get("reason"),
            "signal": has_signal,
        }
        point_reports.append(report)
        if has_signal:
            signals.append(report)
        if has_signal and fit.get("accepted", False) and point["q_norm"] > 0.0:
            resolved.append((point["q_norm"] ** 2, fit["center"]))

    result["points"] = point_reports
    result["finite_q_signal_count"] = len(signals)
    result["resolved_dispersion_points"] = len(resolved)
    if not signals:
        result["status"] = "FAIL"
        result["reasons"].append("no resolved positive-frequency collective signal")

    if len(resolved) < minimum_points:
        if result["status"] == "PASS":
            result["status"] = "INCONCLUSIVE"
        result["reasons"].append(
            f"only {len(resolved)} resolved q points; need {minimum_points} for an independent q^2 fit"
        )
    else:
        stiffness = fit_origin(resolved)
        result["stiffness_fit"] = stiffness
        if stiffness["D_Ry_A2"] <= 0.0 or stiffness["relative_residual"] > dispersion_limit:
            result["status"] = "FAIL"
            result["reasons"].append("resolved low-q energies do not support a stable positive Dq^2 fit")

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument("--results", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--goldstone-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--dispersion-relative-residual", type=float, default=0.20)
    parser.add_argument("--minimum-points", type=int, default=3)
    parser.add_argument("--allow-inconclusive", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    results_root = args.results or repo / "results" / "validation" / "TDVAL-01_FE_NI" / "runs"
    if not results_root.is_absolute():
        results_root = repo / results_root
    results_root = results_root.resolve()
    summary_path = results_root / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    analyses = [
        analyse_run(record, results_root, args.goldstone_tolerance, args.dispersion_relative_residual, args.minimum_points)
        for record in summary.get("runs", [])
    ]
    report = {
        "summary": str(summary_path),
        "goldstone_tolerance": args.goldstone_tolerance,
        "dispersion_relative_residual_limit": args.dispersion_relative_residual,
        "minimum_dispersion_points": args.minimum_points,
        "runs": analyses,
    }
    report_path = args.report or results_root / "physics_analysis.json"
    if not report_path.is_absolute():
        report_path = repo / report_path
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for item in analyses:
        label = f"{item['material']} {item['q_set']} N={item['mesh'][0]}"
        suffix = "; ".join(item["reasons"])
        print(f"{item['status']:12s} {label}: {suffix or 'Goldstone, signal, and q^2 checks passed'}")
        if "stiffness_fit" in item:
            fit = item["stiffness_fit"]
            print(f"  D={fit['D_Ry_A2']:.8g} Ry A^2, relative residual={fit['relative_residual']:.4g}")
    bad = [item for item in analyses if item["status"] == "FAIL"]
    inconclusive = [item for item in analyses if item["status"] == "INCONCLUSIVE"]
    if bad or (inconclusive and not args.allow_inconclusive):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
