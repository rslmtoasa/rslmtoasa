#!/usr/bin/env python3
"""CLI for the WP10 small-q GBT analysis.

The implementation is shared with ``tools/analyze_gbt_wp09.py`` so the WP09
machine-readable input contract has one parser.  This entry point makes the
WP10 deliverable explicit and provides the symmetric q-reversal/quadratic-fit
workflow requested by the Luna prompt.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
WP09 = ROOT / "tools" / "analyze_gbt_wp09.py"


def _load_wp09():
    spec = importlib.util.spec_from_file_location("gbt_wp09_analysis", WP09)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {WP09}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_ANALYZER = _load_wp09()
analyze_wp10_rows = _ANALYZER.analyze_wp10_rows
format_wp10_report = _ANALYZER.format_wp10_report


def _expand_inputs(inputs: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for value in inputs:
        path = Path(value)
        if path.is_dir():
            candidates = sorted(path.rglob("frozen_magnon_harmonic_diagnostics.dat"))
            candidates += sorted(path.rglob("*.csv"))
            if not candidates:
                raise ValueError(f"{path}: directory contains no WP09-compatible sweep files")
            paths.extend(candidates)
        else:
            paths.append(path)
    return paths


def analyze(paths: Sequence[str | Path], **kwargs: Any) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.extend(_ANALYZER.read_sweep_file(path))
    return analyze_wp10_rows(rows, **kwargs)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="WP09 sweep files, CSV files, or directories")
    parser.add_argument("--alat", type=float, dest="alat_angstrom", help="lattice parameter in Angstrom")
    parser.add_argument("--direction", type=float, nargs=3, metavar=("QX", "QY", "QZ"),
                        help="Cartesian q direction; default is the dominant axis of the first nonzero q")
    parser.add_argument("--energy-definition", choices=("auto", "raw", "gauge"), default="auto")
    parser.add_argument("--odd-absolute-tolerance", type=float, default=1.0e-10,
                        help="absolute odd-energy tolerance in Ry")
    parser.add_argument("--odd-relative-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--json", action="store_true", help="print complete JSON")
    parser.add_argument("--json-out", type=Path, help="write complete JSON to this file")
    args = parser.parse_args(argv)
    try:
        result = analyze(_expand_inputs(args.inputs), alat_angstrom=args.alat_angstrom,
                         direction=args.direction, energy_definition=args.energy_definition,
                         odd_absolute_tolerance=args.odd_absolute_tolerance,
                         odd_relative_tolerance=args.odd_relative_tolerance)
    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(payload + "\n")
    print(payload if args.json else format_wp10_report(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
