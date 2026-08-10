#!/usr/bin/env python3
"""Evidence checker for the transverse LR-TDDFT validation campaign.

This utility consumes the deliberately plain-text files written by the
production susceptibility route.  It makes no physics calculation and never
changes a GBT result: GBT and Jij stiffnesses are read only as independently
produced comparison values.

The checker is intentionally standard-library-only so its synthetic fixture
can run in every CI configuration.  A material benchmark is accepted only
when its manifest supplies raw Goldstone data, KS and enhanced spectra, mode
fits, convergence records, and (when available) GBT/Jij comparison values.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


class ValidationError(RuntimeError):
    """A missing or scientifically inconsistent validation record."""


@dataclass(frozen=True)
class StiffnessFit:
    value: float
    relative_residual: float
    points: int


@dataclass(frozen=True)
class EtaExtrapolation:
    zero_eta_fwhm: float
    slope: float
    relative_residual: float


def _records(path: Path) -> list[list[str]]:
    if not path.is_file():
        raise ValidationError(f"required output is missing: {path}")
    return [line.split() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]


def _metadata(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line.startswith("#") or "=" not in line:
            continue
        key, value = line[1:].split("=", 1)
        data[key.strip()] = value.strip()
    return data


def _scalar_record(path: Path, label: str) -> float:
    for columns in _records(path):
        if columns[0] == label and len(columns) >= 2:
            return float(columns[-1])
    raise ValidationError(f"{path}: missing {label!r} record")


def read_goldstone(path: Path) -> dict[str, float | bool]:
    metadata = _metadata(path)
    result: dict[str, float | bool] = {
        "sum_rule_requested": metadata.get("sum_rule_requested", "F").upper().startswith("T"),
        "sum_rule_applied": metadata.get("sum_rule_applied", "F").upper().startswith("T"),
    }
    for columns in _records(path):
        if len(columns) < 2:
            continue
        if columns[0] in {"raw_residual", "raw_magnetization_overlap", "sum_rule_corrected_residual"}:
            result[columns[0]] = float(columns[-1])
    for required in ("raw_residual", "raw_magnetization_overlap"):
        if required not in result:
            raise ValidationError(f"{path}: raw Goldstone diagnostics are incomplete ({required} missing)")
    if result["sum_rule_applied"] and "sum_rule_corrected_residual" not in result:
        raise ValidationError(f"{path}: sum-rule correction overwrote rather than supplemented raw diagnostics")
    return result


def read_mode_fits(path: Path) -> dict[int, dict[str, float | bool | str]]:
    fits: dict[int, dict[str, float | bool | str]] = {}
    candidates: set[int] = set()
    crossings: dict[int, dict[str, float | bool]] = {}
    for columns in _records(path):
        if columns[0] == "candidate" and len(columns) >= 7:
            candidates.add(int(columns[1]))
        elif columns[0] == "fit" and len(columns) >= 7:
            fits[int(columns[1])] = {
                "accepted": columns[2].upper().startswith("T"),
                "center": float(columns[3]),
                "fwhm": float(columns[4]),
                "hwhm": float(columns[5]),
                "relative_residual": float(columns[6]),
                "reason": " ".join(columns[7:]),
            }
        elif columns[0] == "crossing" and len(columns) >= 10:
            crossings[int(columns[1])] = {
                "present": columns[2].upper().startswith("T"),
                "omega": float(columns[3]),
                "imaginary_part": float(columns[4]),
                "branch_overlap": float(columns[5]),
                "eigenvalue_step": float(columns[6]),
                "projected_weight": float(columns[7]),
                "condition_number": float(columns[8]),
                "exceptional_warning": columns[9].upper().startswith("T"),
            }
    if not fits or candidates != set(fits) or candidates != set(crossings):
        raise ValidationError(f"{path}: candidate/crossing/mode-fit records are incomplete")
    for index, fit in fits.items():
        fit["crossing"] = crossings[index]
    return fits


def read_trace_spectrum(path: Path, label: str) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for columns in _records(path):
        # chi_KS records begin with omega, whereas Dyson records begin with
        # their kind.  Supporting both is intentional: both writers are
        # production interfaces and the campaign must inspect the bare and
        # enhanced spectrum separately.
        if columns[0] == label and len(columns) >= 3:
            rows.append((float(columns[1]), float(columns[-1])))
        elif len(columns) >= 3 and columns[1] == label:
            rows.append((float(columns[0]), float(columns[-1])))
    if not rows:
        raise ValidationError(f"{path}: no {label} spectrum records")
    return rows


def fit_stiffness(q_norms: list[float], omega: list[float]) -> StiffnessFit:
    if len(q_norms) != len(omega) or len(q_norms) < 2:
        raise ValidationError("stiffness fit requires at least two matched q/omega records")
    x = [q * q for q in q_norms]
    active = [(q2, w) for q2, w in zip(x, omega) if q2 > 0.0]
    if len(active) < 2:
        raise ValidationError("stiffness fit requires at least two nonzero q points")
    denominator = sum(q2 * q2 for q2, _ in active)
    if denominator == 0.0:
        raise ValidationError("degenerate q values in stiffness fit")
    value = sum(q2 * w for q2, w in active) / denominator
    scale = max(max(abs(w) for _, w in active), 1.0e-30)
    residual = math.sqrt(sum((w - value * q2) ** 2 for q2, w in active) / len(active)) / scale
    return StiffnessFit(value, residual, len(active))


def extrapolate_zero_eta(eta: list[float], fwhm: list[float]) -> EtaExtrapolation:
    if len(eta) != len(fwhm) or len(eta) < 3 or any(x <= 0.0 for x in eta):
        raise ValidationError("eta extrapolation needs three or more positive eta/observed-FWHM pairs")
    n = float(len(eta))
    mean_x, mean_y = sum(eta) / n, sum(fwhm) / n
    denom = sum((x - mean_x) ** 2 for x in eta)
    if denom == 0.0:
        raise ValidationError("eta sweep contains only one eta value")
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(eta, fwhm)) / denom
    intercept = mean_y - slope * mean_x
    scale = max(max(abs(y) for y in fwhm), 1.0e-30)
    residual = math.sqrt(sum((y - (intercept + slope * x)) ** 2 for x, y in zip(eta, fwhm)) / len(eta)) / scale
    return EtaExtrapolation(intercept, slope, residual)


def _q_norm(q: list[float]) -> float:
    if len(q) != 3:
        raise ValidationError("each q_direct entry must contain exactly three fractional reciprocal coordinates")
    return math.sqrt(sum(float(component) ** 2 for component in q))


def analyse_campaign(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text())
    root = manifest_path.parent
    goldstone = read_goldstone(root / manifest["goldstone_file"])
    modes = read_mode_fits(root / manifest["modes_file"])
    q_points = manifest["q_points"]
    q_norms, centers = [], []
    for index, point in enumerate(q_points, start=1):
        fit = modes.get(index)
        if fit is None or not fit["accepted"]:
            raise ValidationError(f"q index {index}: no accepted coherent-mode fit; do not use it for stiffness")
        crossing = fit["crossing"]
        assert isinstance(crossing, dict)
        if not crossing["present"] or crossing["exceptional_warning"]:
            raise ValidationError(f"q index {index}: fit lacks a well-conditioned Xi unity crossing")
        if float(crossing["projected_weight"]) <= 0.0:
            raise ValidationError(f"q index {index}: Xi crossing has no positive mode-projected enhanced weight")
        q_norms.append(_q_norm(point["q_direct"]))
        centers.append(float(fit["center"]))
    stiffness = fit_stiffness(q_norms, centers)

    chi0 = read_trace_spectrum(root / manifest["chi0_file"], "trace")
    enhanced = read_trace_spectrum(root / manifest["dyson_file"], "trace_loss")
    if min(weight for _, weight in chi0) < -1.0e-12:
        raise ValidationError("KS/Stoner trace has negative reported spectral weight")
    if max(weight for _, weight in enhanced) <= max(weight for _, weight in chi0):
        raise ValidationError("enhanced spectrum has no spectral feature above the bare KS/Stoner trace")

    eta_data = manifest["eta_sweep"]
    eta_fit = extrapolate_zero_eta([float(row["eta_Ry"]) for row in eta_data], [float(row["observed_fwhm_Ry"]) for row in eta_data])
    if eta_fit.zero_eta_fwhm < -1.0e-12:
        raise ValidationError("zero-eta linewidth extrapolates to an unphysical negative width")

    convergence = manifest["convergence"]
    required_axes = {"k_mesh", "band_window", "response_projection", "electronic_smearing", "eta", "frequency_grid"}
    found_axes = {row["axis"] for row in convergence}
    missing_axes = required_axes - found_axes
    if missing_axes:
        raise ValidationError(f"missing documented convergence axes: {', '.join(sorted(missing_axes))}")

    limits = manifest["acceptance"]
    if float(goldstone["raw_residual"]) > float(limits["max_raw_goldstone_residual"]):
        raise ValidationError("raw Goldstone residual exceeds campaign limit")
    if float(goldstone["raw_magnetization_overlap"]) < float(limits["min_magnetization_overlap"]):
        raise ValidationError("Goldstone unity-mode overlap is unstable")
    if stiffness.relative_residual > float(limits["max_stiffness_relative_residual"]):
        raise ValidationError("small-q dispersion is not quadratic within the campaign limit")
    if eta_fit.relative_residual > float(limits["max_eta_fit_relative_residual"]):
        raise ValidationError("linewidth cannot be cleanly extrapolated over the selected eta sweep")

    routes: dict[str, dict[str, float]] = {"LR-TDDFT": {"D_Ry_per_q2": stiffness.value}}
    for route in ("GBT", "Jij"):
        if route in manifest.get("independent_routes", {}):
            reference = float(manifest["independent_routes"][route]["D_Ry_per_q2"])
            delta = abs(stiffness.value - reference) / max(abs(reference), 1.0e-30)
            routes[route] = {"D_Ry_per_q2": reference, "relative_delta_to_LR": delta}

    return {
        "material": manifest["material"],
        "goldstone": goldstone,
        "stiffness": {"D_Ry_per_q2": stiffness.value, "relative_residual": stiffness.relative_residual, "points": stiffness.points},
        "linewidth": {"zero_eta_fwhm_Ry": eta_fit.zero_eta_fwhm, "slope_per_Ry": eta_fit.slope,
                      "relative_residual": eta_fit.relative_residual,
                      "interpretation": "observed widths are numerical-eta dependent; only the zero-eta intercept is reported"},
        "routes": routes,
        "convergence_axes": sorted(found_axes),
        "status": "evidence fixture passed; material credibility still requires a real benchmark manifest",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="JSON manifest for one material or deterministic fixture")
    parser.add_argument("--report", type=Path, help="write the machine-readable evidence summary here")
    args = parser.parse_args()
    try:
        summary = analyse_campaign(args.manifest.resolve())
    except (OSError, ValueError, KeyError, ValidationError) as error:
        print(f"FAIL: {error}")
        return 1
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
