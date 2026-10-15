#!/usr/bin/env python3
"""Build a backend-matched M0.2 floor from accepted no-U captures.

The central reference is the signed discrepancy D = sum(q0*c) - ATOMSC SUMEV
from the baseline accepted U=0 state.  Its absolute value is the systematic
ASA inversion/relaxation floor.  Uncertainty is a separate maximum absolute
deviation from that central D across deterministic repeats and resolution
refinements.  Do not use max(abs(D_i)) as the floor and then add deviations
again: that would count the observed variation twice.

The resulting `floor+uncertainty` is an engineering allowance for comparing
the independent scalar-v control.  It is not an ETOT correction and does not
establish M0.3 or production energy closure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import local_moment_capture_oracle as m0


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolve_input_hash(manifest: dict[str, Any], suffix: str) -> str:
    matches = [value for name, value in manifest.get("input_sha256", {}).items()
               if Path(name).name == suffix]
    if len(matches) != 1:
        raise ValueError(f"manifest must identify exactly one {suffix} hash")
    return str(matches[0])


def _resolution(manifest: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    backend = str(manifest["backend"]).upper()
    if backend == "KS":
        return {"kind": "kpoint_count", "value": int(manifest["kpoint_count"])}
    energies: dict[int, float] = {}
    for row in rows:
        if row["kind"] == "green_function" and row["stage"] == "RS_GREEN":
            k = m0._i(row, "k")
            e = m0._f(row, "energy")
            if k in energies and energies[k] != e:
                raise ValueError(f"RS energy index {k} has multiple energies")
            energies[k] = e
    if not energies:
        raise ValueError("RS capture has no Green-function energy grid")
    grid = [energies[k] for k in sorted(energies)]
    return {"kind": "rs_energy_grid", "value": len(grid), "emin_ry": min(grid), "emax_ry": max(grid)}


def _accepted_log(directory: Path) -> float:
    log_path = directory / "testrun.log"
    if not log_path.is_file():
        raise ValueError(f"missing convergence log {log_path}")
    log = log_path.read_text(encoding="utf-8", errors="replace")
    matches = list(re.finditer(r"Converged!\s*([0-9.eEdD+-]+)", log))
    if not matches:
        raise ValueError(f"{log_path} has no final `Converged! <residual>` marker")
    last = matches[-1]
    if re.search(r"Not converged!", log[last.end():]):
        raise ValueError(f"{log_path} has a nonconverged update after its converged marker")
    residual = float(last.group(1).replace("D", "E").replace("d", "e"))
    if not (residual >= 0.0):
        raise ValueError(f"invalid convergence residual in {log_path}: {residual}")
    return residual


def _pre_input(directory: Path) -> str:
    path = directory / "input.nml.pre"
    if not path.is_file():
        raise ValueError(f"missing immutable pre-run input {path}")
    return path.read_text(encoding="utf-8")


def _without_resolution(text: str, backend: str) -> str:
    """Remove only the preregistered resolution field(s), retaining all else."""
    names = (r"nk[123]" if backend == "KS" else r"channels_ldos")
    return re.sub(rf"(?im)^\s*(?:{names})\s*=\s*[^!\n]*(?:!.*)?$", "<resolution-field>", text)


def _input_number(text: str, key: str) -> float:
    match = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*([-+0-9.eEdD]+)", text)
    if not match:
        raise ValueError(f"input.nml.pre is missing {key}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def _load_run(directory: Path) -> dict[str, Any]:
    manifest, rows = m0.load_capture(directory)
    if int(manifest.get("nsp", -1)) != 1 or int(manifest.get("lmax", -1)) != 2:
        raise ValueError(f"{directory}: floor protocol requires nsp=1, lmax=2")
    if manifest.get("metric") != "identity":
        raise ValueError(f"{directory}: unsupported/nontrivial metric cannot define this floor")
    result = m0.run(directory, mode="m02-no-u")
    d = float(result["band_m1_sum_q0c_ry"] - result["production_atomsc_sumev_ry"])
    output_seed = directory / "Fe_out.nml"
    if not output_seed.is_file():
        raise ValueError(f"missing accepted output restart seed {output_seed}")
    pre_input = _pre_input(directory)
    threshold = _input_number(pre_input, "conv_thr")
    convergence_residual = _accepted_log(directory)
    if threshold <= 0.0 or convergence_residual > threshold:
        raise ValueError(f"{directory}: final convergence residual {convergence_residual:.6g} exceeds conv_thr {threshold:.6g}")
    run_manifest_path = directory / "run_manifest.json"
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8")) if run_manifest_path.exists() else {}
    if run_manifest and (run_manifest.get("hubbard_fields_present") is not False
                         or run_manifest.get("hubbard_U_J_verified_zero_in_restart_seed", True) is not True):
        raise ValueError(f"{directory}: run manifest does not verify a no-U restart/setup")
    return {
        "directory": str(directory.resolve()),
        "backend": str(manifest["backend"]).upper(),
        "run_id": manifest["run_id"],
        "iteration": int(manifest["iteration"]),
        "csv_sha256": manifest["csv_sha256"],
        "executable_sha256": manifest["executable_sha256"],
        "source_sha": manifest.get("source_sha"),
        "input_sha256": manifest.get("input_sha256", {}),
        "fe_seed_sha256": _resolve_input_hash(manifest, "Fe.nml.pre"),
        "fe_output_seed_sha256": _hash(output_seed),
        "input_deck_sha256": _resolve_input_hash(manifest, "input.nml.pre"),
        "input_without_resolution": _without_resolution(pre_input, str(manifest["backend"]).upper()),
        "declared_fermi_ry": _input_number(pre_input, "fermi"),
        "resolution": _resolution(manifest, rows),
        "electron_count": float(manifest["electron_count"]),
        "fermi_ry": float(manifest["fermi_ry"]),
        "temperature_ry": float(manifest["temperature_ry"]),
        "convergence_threshold": threshold,
        "convergence_residual": convergence_residual,
        "D_m1_minus_sumev_ry": d,
        "abs_D_ry": abs(d),
        "m02_metrics": result,
        "log_sha256": _hash(directory / "testrun.log"),
    }


def build_floor(baseline_dir: Path, repeats: list[Path], refinements: list[Path]) -> dict[str, Any]:
    base = _load_run(baseline_dir)
    if not repeats:
        raise ValueError("at least one deterministic repeat is required to estimate uncertainty")
    if not refinements:
        raise ValueError("at least one numerical refinement is required to estimate uncertainty")
    repeat_runs = [_load_run(p) for p in repeats]
    refinement_runs = [_load_run(p) for p in refinements]
    all_runs = [base, *repeat_runs, *refinement_runs]
    backend = base["backend"]
    for run in all_runs:
        for key in ("backend", "executable_sha256", "source_sha", "electron_count",
                    "temperature_ry", "declared_fermi_ry"):
            if run[key] != base[key]:
                raise ValueError(f"run {run['run_id']} differs from baseline in {key}")
        if run["backend"] != backend:
            raise ValueError("mixed-backend M0.2 floor is forbidden")
    for run in repeat_runs:
        if run["fe_seed_sha256"] != base["fe_output_seed_sha256"]:
            raise ValueError(f"repeat {run['run_id']} is not seeded from the accepted baseline Fe_out")
        if run["input_deck_sha256"] != base["input_deck_sha256"] or run["resolution"] != base["resolution"]:
            raise ValueError(f"repeat {run['run_id']} is not the same-input/same-resolution control")
    for run in refinement_runs:
        if run["fe_seed_sha256"] != base["fe_output_seed_sha256"]:
            raise ValueError(f"refinement {run['run_id']} is not seeded from the accepted baseline Fe_out")
        if run["resolution"] == base["resolution"]:
            raise ValueError(f"refinement {run['run_id']} does not increase numerical resolution")
        if run["resolution"]["kind"] != base["resolution"]["kind"]:
            raise ValueError("refinement changes numerical resolution kind")
        if run["resolution"]["value"] <= base["resolution"]["value"]:
            raise ValueError(f"refinement {run['run_id']} is not finer than baseline")
        if run["input_without_resolution"] != base["input_without_resolution"]:
            raise ValueError(f"refinement {run['run_id']} changes physical input beyond the declared resolution field")

    d0 = base["D_m1_minus_sumev_ry"]
    repeat_deviations = [abs(run["D_m1_minus_sumev_ry"] - d0) for run in repeat_runs]
    refinement_deviations = [abs(run["D_m1_minus_sumev_ry"] - d0) for run in refinement_runs]
    uncertainty = max([0.0, *repeat_deviations, *refinement_deviations])
    floor = abs(d0)
    return {
        "schema": "rslmto-m02-u0-floor-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "CALIBRATION_ENVELOPE_ONLY",
        "backend": backend,
        "definition": {
            "D_ry": "sum of same-state Re Tr(M1) over all recorded shells/sites minus production ATOMSC SUMEV",
            "floor_ry": "abs(D_baseline), the central systematic inversion/relaxation offset",
            "uncertainty_ry": "max abs(D_i-D_baseline) over deterministic repeats and numerical refinements",
            "repeat_deviation_ry": "max abs(D_repeat-D_baseline); input deck, resolution, seed, source, executable must match",
            "refinement_deviation_ry": "max abs(D_refined-D_baseline); same seed, source, executable and physical setup, finer spectral resolution",
            "allowance_ry": "floor+uncertainty; max(abs(D_i)) is reported separately and is not added again",
        },
        "max_sum_ev_floor_ry": floor,
        "uncertainty_ry": uncertainty,
        "allowance_ry": floor + uncertainty,
        "max_abs_D_diagnostic_ry": max(run["abs_D_ry"] for run in all_runs),
        "max_repeat_deviation_ry": max(repeat_deviations, default=0.0),
        "max_refinement_deviation_ry": max(refinement_deviations, default=0.0),
        "distinct_residuals_not_included_in_allowance": {
            "max_m1_independent_residual_ry": max(run["m02_metrics"]["max_production_independent_m1_error_ry"] for run in all_runs),
            "max_trace_q0c_residual_ry": max(run["m02_metrics"]["max_m1_trace_q0c_error_ry"] for run in all_runs),
            "max_replay_sumev_residual_ry": max(abs(run["m02_metrics"]["newrho_replay_difference_from_atomsc_ry"]) for run in all_runs),
        },
        "runs": all_runs,
        "limitations": [
            "A floor is not available until the baseline, at least one deterministic repeat, and a refinement run all pass their independent M0.2 checks.",
            "This scalar U=0 envelope does not prove M0.3 carrier validity, total-energy closure, or an ETOT correction.",
            "Scalar-v controls and U-active results must be compared against this backend-matched allowance without retuning it.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--repeat", required=True, action="append", type=Path)
    parser.add_argument("--refinement", action="append", default=[], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        report = build_floor(args.baseline, args.repeat, args.refinement)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"m02_floor_oracle: FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({k: report[k] for k in (
        "schema", "status", "backend", "max_sum_ev_floor_ry", "uncertainty_ry",
        "allowance_ry", "max_abs_D_diagnostic_ry", "max_repeat_deviation_ry",
        "max_refinement_deviation_ry")}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
