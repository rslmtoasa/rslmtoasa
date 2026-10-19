#!/usr/bin/env python3
"""Extract one direct d-shell P_OUT and evaluate fixed-P production E_Hub."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def export_p(trace: Path, iteration: int, site: int, out: Path) -> dict:
    found: dict[tuple[int, int, int], tuple[str, str]] = {}
    with trace.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if (row.get("kind") != "test_density_matrix" or row.get("stage") != "P_OUT" or
                    int(row.get("iteration", "-1")) != iteration or int(row.get("site", "-1")) != site or
                    int(row.get("l", "-1")) != 2):
                continue
            spin, i, j = int(row["spin"]), int(row["row"]), int(row["col"])
            if spin not in (0, 1) or i not in range(5) or j not in range(5):
                raise ValueError(f"invalid d-shell matrix index: {row}")
            re_value, im_value = row["re"].strip(), row["im"].strip()
            if not math.isfinite(float(re_value)) or not math.isfinite(float(im_value)):
                raise ValueError("non-finite P_OUT entry")
            key = (spin, i, j)
            if key in found:
                raise ValueError(f"duplicate P_OUT matrix entry: {key}")
            found[key] = (re_value, im_value)
    if len(found) != 50:
        raise ValueError(f"expected exactly 50 d-shell P_OUT entries, got {len(found)}")
    with out.open("w") as stream:
        for (spin, i, j), (real, imag) in sorted(found.items()):
            stream.write(f"{spin + 1} {i + 1} {j + 1} {real} {imag}\n")
    return {"matrix_entry_count": len(found), "P_matrix_path": str(out), "P_matrix_sha256": sha(out)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("iteration", type=int)
    parser.add_argument("u_eV", type=float)
    parser.add_argument("j_eV", type=float)
    parser.add_argument("--site", type=int, default=0)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    matrix = args.outdir / "P_OUT_d_shell.txt"
    meta = export_p(args.trace, args.iteration, args.site, matrix)
    proc = subprocess.run([str(args.driver), str(matrix), f"{args.u_eV:.16g}", f"{args.j_eV:.16g}"],
                          check=False, capture_output=True, text=True)
    stdout = args.outdir / "fixed_P_evaluator.stdout"
    stderr = args.outdir / "fixed_P_evaluator.stderr"
    stdout.write_text(proc.stdout)
    stderr.write_text(proc.stderr)
    result = {
        "schema": "c3-rc2-fixed-p-energy-v1",
        "trace_path": str(args.trace), "trace_sha256": sha(args.trace),
        "iteration": args.iteration, "site_zero_based": args.site,
        "u_center_eV": args.u_eV, "j_center_eV": args.j_eV,
        **meta,
        "driver_path": str(args.driver), "driver_sha256": sha(args.driver),
        "command": [str(args.driver), str(matrix), f"{args.u_eV:.16g}", f"{args.j_eV:.16g}"],
        "returncode": proc.returncode,
        "stdout_path": str(stdout), "stdout_sha256": sha(stdout),
        "stderr_path": str(stderr), "stderr_sha256": sha(stderr),
    }
    (args.outdir / "FIXED_P_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))
    if proc.returncode:
        raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
