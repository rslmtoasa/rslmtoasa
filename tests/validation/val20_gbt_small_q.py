#!/usr/bin/env python3
"""Run and analyse the WP10 symmetric small-q GBT sweep.

The runner reuses the committed WP09 bcc-Fe force-theorem fixture and its
raw-output collation helper.  It changes only the mesh, cone angle, and
symmetric q list.  The resulting CSV is retained before the WP10 analyzer is
called, so the fit never replaces the underlying energies.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
WP09_RUNNER = ROOT / "tests" / "validation" / "val19_gbt_harmonic_kgrid.py"
WP10_TOOL = ROOT / "tools" / "analyze_gbt_wp10.py"
Q_FILE = ROOT / "tests" / "regression" / "wp10_validation" / "small_q" / "q_points_symmetric.dat"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--angles", type=float, nargs="+", default=[2.0, 5.0, 10.0, 15.0])
    parser.add_argument("--meshes", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--modes", default="mft", help="comma-separated modes; mft is the primary WP10 route")
    parser.add_argument("--alat", type=float, default=None, help="lattice parameter in Angstrom; infer from output when available")
    parser.add_argument("--require-odd-within-tolerance", action="store_true")
    args = parser.parse_args(argv)
    args.binary = args.binary.resolve()
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    modes = [item.strip() for item in args.modes.split(",") if item.strip()]
    wp09 = _load(WP09_RUNNER, "gbt_wp09_runner")
    wp10 = _load(WP10_TOOL, "gbt_wp10_tool")
    rows: list[dict[str, Any]] = []

    for mode in modes:
        if mode not in {"mft", "mft_constrained"}:
            raise ValueError(f"unsupported mode {mode!r}; WP10 runner supports mft or mft_constrained")
        meshes = args.meshes if mode == "mft" else [args.meshes[-1]]
        for nk in meshes:
            for theta in args.angles:
                workdir = args.scratch_root / f"{mode}_nk{nk}_theta{theta:g}"
                rows.extend(wp09._run(args.binary, workdir, nk, theta, mode, Q_FILE, args.timeout))

    output_csv = args.scratch_root / "wp10_sweep.csv"
    wp09._write_csv(output_csv, rows)
    alat_values = {row.get("alat_angstrom") for row in rows if row.get("alat_angstrom") is not None}
    if args.alat is None and len(alat_values) == 1:
        args.alat = float(alat_values.pop())
    if args.alat is not None:
        output_csv.write_text(f"# alat_angstrom = {args.alat:.16g}\n" + output_csv.read_text())
    result = wp10.analyze_wp10_rows(rows, alat_angstrom=args.alat, direction=(0.0, 0.0, 1.0))
    output_json = args.scratch_root / "wp10_analysis.json"
    output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(wp10.format_wp10_report(result))
    print(f"\nraw sweep: {output_csv}\nanalysis JSON: {output_json}")

    if args.require_odd_within_tolerance:
        failures = [item for item in result["odd_component_status"] if item["status"] != "within_tolerance"]
        if failures:
            print(f"FAIL: {len(failures)} group(s) have an odd q component outside tolerance")
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
