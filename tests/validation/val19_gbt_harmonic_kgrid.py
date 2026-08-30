#!/usr/bin/env python3
"""WP09 runner for raw cone-angle and reciprocal k-grid convergence.

This is an intentionally data-producing validation.  It runs the small bcc-Fe
fixture for bare MFT at several cone angles and meshes, and can also run the
corrected constrained-MFT reference fixture.  The latter uses the stable WP08
theta=90 degree, q=(0,0,+/-0.025) convention; it is included for comparison
and is not silently presented as a small-angle harmonic plateau.  The runner
collates production diagnostic files without changing their energies and
delegates all analysis to ``tools/analyze_gbt_wp09.py``.  No material-specific
energy golden is embedded here.  Use ``--require-plateau`` when a campaign
should be promoted from an evidence run to a gate.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
MFT_DECK = ROOT / "tests" / "regression" / "wp9_validation" / "gammaH_sweep" / "base_kspace"
CMFT_DECK = ROOT / "example" / "bulk" / "bccFe"
FIXTURE = ROOT / "tests" / "regression" / "wp9_validation" / "harmonic_kgrid"
ANALYZER = ROOT / "tools" / "analyze_gbt_wp09.py"
CORRECTED_REFERENCE = FIXTURE / "q_points_constrained_reference.dat"

NK_RE = re.compile(r"^(\s*nk([123])\s*=\s*)\d+", re.MULTILINE)
THETA_RE = re.compile(r"^(\s*theta_ss\s*=\s*)[0-9.eEdD+-]+", re.MULTILINE)
MODE_RE = re.compile(r"^(\s*mode\s*=\s*)'[^']*'", re.MULTILINE)
QFILE_RE = re.compile(r"^(\s*q_file\s*=\s*)'[^']*'", re.MULTILINE)

CSV_COLUMNS = [
    "mode", "theta_deg", "q1", "q2", "q3", "e_q_theta", "e_q0", "e_qref_theta", "e_qref0",
    "delta_raw", "delta_gauge", "delta_pure", "sin2theta", "mtot", "omega", "fermi_level",
    "electron_count", "electron_error", "target_electrons", "weight_sum", "nk1", "nk2", "nk3",
    "nk_total", "gauge_available", "q_half_maps_to_mesh",
]


def _load_analyzer():
    spec = importlib.util.spec_from_file_location("gbt_wp09_analyzer", ANALYZER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {ANALYZER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _patch_input(text: str, nk: int, theta: float, mode: str) -> str:
    if mode == "mft_constrained":
        text, count = THETA_RE.subn(rf"\g<1>{theta}", text)
        if count != 1:
            raise RuntimeError("expected one theta_ss in the constrained fixture")
        text, count = MODE_RE.subn(rf"\g<1>'{mode}'", text)
        if count != 1:
            raise RuntimeError("expected one frozen-magnon mode in the constrained fixture")
        return text
    text, count = re.subn(r"^(\s*nk[123]\s*=\s*)\d+", rf"\g<1>{nk}", text, flags=re.MULTILINE)
    if count != 3:
        raise RuntimeError(f"expected nk1/nk2/nk3 in the fixture, found {count}")
    text, count = THETA_RE.subn(rf"\g<1>{theta}", text)
    if count != 1:
        raise RuntimeError("expected one theta_ss in the fixture")
    text, count = MODE_RE.subn(rf"\g<1>'{mode}'", text)
    if count != 1:
        raise RuntimeError("expected one frozen-magnon mode in the fixture")
    text, count = QFILE_RE.subn(r"\g<1>'q_points.dat'", text)
    if count != 1:
        raise RuntimeError("expected one q_file in the fixture")
    if mode == "mft_constrained":
        text += "\n&constraints\nconstraints_enable = .true.\nconstraints_i_cons = 3\nconstraints_tolerance = 1.0d-6\nconstraints_diagnostics = .false.\n/\n"
    return text


def _run(binary: Path, workdir: Path, nk: int, theta: float, mode: str, q_file: Path, timeout: int) -> list[dict[str, Any]]:
    if workdir.exists():
        shutil.rmtree(workdir)
    source_deck = CMFT_DECK if mode == "mft_constrained" else MFT_DECK
    shutil.copytree(source_deck, workdir)
    shutil.copy2(q_file, workdir / "q_points.dat")
    if mode == "mft_constrained":
        input_path = workdir / "input_frozen_magnon_constrained.nml"
        input_text = input_path.read_text()
        input_path.rename(workdir / "input.nml")
        input_path = workdir / "input.nml"
        shutil.copy2(q_file, workdir / "q_wp08.dat")
    else:
        input_path = workdir / "input.nml"
        input_text = input_path.read_text()
    input_path.write_text(_patch_input(input_text, nk, theta, mode))
    result = subprocess.run([str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=timeout)
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{workdir.name} exited {result.returncode}:\n{result.stdout[-4000:]}")
    analyzer = _load_analyzer()
    output = workdir / "frozen_magnon_harmonic_diagnostics.dat"
    if output.exists():
        return analyzer.read_sweep_file(output)
    # The corrected WP08 reference fixture uses the real-space route.  Keep
    # its raw eband rows in the same schema, with mesh/Fermi/gauge fields
    # unavailable rather than treating a real-space run as a k-grid sample.
    legacy = workdir / "frozen_magnon.dat"
    if legacy.exists() and mode == "mft_constrained":
        return analyzer.read_sweep_file(legacy, mode=mode, theta_deg=theta)
    raise RuntimeError(f"{workdir}: production WP09 diagnostic output was not written")


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            values: dict[str, Any] = {
                "mode": row["mode"], "theta_deg": row.get("theta_deg"),
                "q1": row["q"][0], "q2": row["q"][1], "q3": row["q"][2],
                "gauge_available": int(bool(row["gauge_available"])) if row.get("gauge_available") is not None else None,
                "q_half_maps_to_mesh": int(bool(row["q_half_maps_to_mesh"])) if row.get("q_half_maps_to_mesh") is not None else None,
            }
            for key in CSV_COLUMNS:
                if key in values:
                    continue
                if key in {"nk1", "nk2", "nk3"}:
                    mesh = row.get("mesh")
                    values[key] = mesh[int(key[-1]) - 1] if mesh else None
                elif key == "nk_total":
                    values[key] = row.get("nk_total")
                else:
                    values[key] = row.get(key)
            writer.writerow(values)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--angles", type=float, nargs="+", default=[2.0, 5.0, 10.0, 15.0, 20.0])
    parser.add_argument("--meshes", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--modes", default="mft,mft_constrained", help="comma-separated frozen-magnon modes")
    parser.add_argument("--corrected-theta", type=float, default=90.0,
                        help="reference angle for the corrected constrained-MFT campaign")
    parser.add_argument("--require-plateau", action="store_true", help="fail if a mode has no admitted angle plateau")
    args = parser.parse_args(argv)
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    modes = [item.strip() for item in args.modes.split(",") if item.strip()]
    analyzer = _load_analyzer()
    rows: list[dict[str, Any]] = []

    for mode in modes:
        if mode == "mft":
            for nk in args.meshes:
                for theta in args.angles:
                    rows.extend(_run(args.binary, args.scratch_root / f"{mode}_nk{nk}_theta{theta:g}", nk, theta,
                                     mode, FIXTURE / "q_points_cone.dat", args.timeout))
                # One finite-q control outside the cone sweep makes the q/2
                # mapping and k-grid energy convergence visible at each mesh.
                rows.extend(_run(args.binary, args.scratch_root / f"{mode}_nk{nk}_grid", nk, args.angles[0], mode,
                                 FIXTURE / "q_points_grid.dat", args.timeout))
        elif mode == "mft_constrained":
            # The current corrected-MFT implementation is demonstrably stable
            # for the original WP08 reference convention.  Keep that evidence
            # explicit; a failed attempt to infer a corrected small-angle
            # plateau must not be hidden by this WP09 runner.
            nk = args.meshes[-1]
            rows.extend(_run(args.binary, args.scratch_root / f"{mode}_reference_theta{args.corrected_theta:g}", nk,
                             args.corrected_theta, mode, CORRECTED_REFERENCE, args.timeout))
        else:
            raise ValueError(f"unsupported WP09 mode {mode!r}; use mft or mft_constrained")

    output_csv = args.scratch_root / "wp09_sweep.csv"
    output_json = args.scratch_root / "wp09_analysis.json"
    _write_csv(output_csv, rows)
    result = analyzer.analyze_rows(rows)
    output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(analyzer.format_report(result))
    print(f"\nraw sweep: {output_csv}\nanalysis JSON: {output_json}")

    if args.require_plateau:
        required = {mode for mode in modes}
        found = {item["mode"] for item in result["plateaus"] if item["status"] == "plateau" and item["admitted_angles_deg"]}
        missing = sorted(required - found)
        if missing:
            print(f"FAIL: no admitted harmonic plateau for mode(s): {', '.join(missing)}")
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
