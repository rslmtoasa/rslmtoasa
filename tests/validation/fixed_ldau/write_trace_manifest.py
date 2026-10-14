#!/usr/bin/env python3
"""Create a provenance manifest for an opt-in fixed-U trace capture.

CSV rows are scanned incrementally; large RS captures need not be loaded into
memory. Input conditions that cannot be inferred safely (N, kT, EF, lambda)
are required on the command line.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

SCHEMA = "rslmto-fixed-u-trace-v1"
FIELDS = ("kind", "run_id", "lambda", "iteration", "stage", "site", "type", "l", "spin",
          "k", "band", "row", "col", "order", "weight", "energy", "scalar", "re", "im")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path, help="directory containing records.csv")
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--backend", choices=("RS", "KS"), required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--input", type=Path, action="append", required=True,
                        help="input deck(s) to hash; may be repeated")
    parser.add_argument("--lambda", dest="lam", type=float, required=True)
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--nsp", type=int, required=True)
    parser.add_argument("--temperature-ry", type=float, required=True)
    parser.add_argument("--electron-count", type=float, required=True)
    parser.add_argument("--fermi-ry", type=float, required=True)
    parser.add_argument("--metric", choices=("identity", "overlap"), default="identity")
    parser.add_argument("--lmax", type=int, help="local basis lmax (needed if the trace predates basis_map rows)")
    parser.add_argument("--norb-per-spin", type=int, help="local orbital count, e.g. 9 for spd")
    parser.add_argument("--source-id", help="honest source-build identifier when the working tree changed after the binary")
    parser.add_argument("--executable-sha256", help="record a previously measured executable hash when its file is unavailable")
    parser.add_argument("--output", type=Path, help="manifest path (default capture/manifest.json)")
    args = parser.parse_args()

    csv_path = args.capture / "records.csv"
    git = lambda *cmd: subprocess.check_output(["git", "-C", str(args.repo), *cmd], text=True).strip()
    revision = git("rev-parse", "HEAD")
    diff = subprocess.check_output(["git", "-C", str(args.repo), "diff", "--binary", "HEAD"])
    diff_sha = hashlib.sha256(diff).hexdigest()
    run_id = None
    max_index = -1
    k_weights: set[int] = set()
    min_iteration = max_iteration = None
    first_lambda = None
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise SystemExit(f"unexpected CSV header: {reader.fieldnames!r}")
        for row in reader:
            run_id = run_id or row["run_id"]
            if row["run_id"] != run_id:
                raise SystemExit("CSV contains multiple run_id values")
            lam = float(row["lambda"])
            first_lambda = lam if first_lambda is None else first_lambda
            if abs(lam - args.lam) > 1e-12:
                raise SystemExit(f"CSV lambda={lam} disagrees with requested {args.lam}")
            iteration = int(row["iteration"])
            min_iteration = iteration if min_iteration is None else min(min_iteration, iteration)
            max_iteration = iteration if max_iteration is None else max(max_iteration, iteration)
            if row["kind"] == "basis_map":
                max_index = max(max_index, int(row["row"]))
            elif row["kind"] in ("h_total", "s_overlap", "green_function", "eigenvector"):
                max_index = max(max_index, int(row["row"]), int(row["col"]))
            elif row["kind"] == "eigenvalue":
                max_index = max(max_index, int(row["band"]))
            elif row["kind"] == "k_weight":
                k_weights.add(int(row["k"]))
    if run_id is None or min_iteration is None:
        raise SystemExit("CSV contains no records")
    if not min_iteration <= args.iteration <= max_iteration:
        raise SystemExit(f"requested iteration {args.iteration} outside CSV range {min_iteration}..{max_iteration}")
    input_hashes = {str(path): sha256(path) for path in args.input}
    manifest = {
        "schema": SCHEMA,
        "run_id": run_id,
        "source_sha": args.source_id or f"{revision}+working-diff:{diff_sha}",
        "source_revision": revision,
        "working_tree_diff_sha256": None if args.source_id else diff_sha,
        "executable_sha256": args.executable_sha256 or sha256(args.executable),
        "executable_hash_source": "reported_measurement" if args.executable_sha256 else "file",
        "input_sha256": input_hashes,
        "backend": args.backend,
        "lambda": args.lam,
        "iteration": args.iteration,
        "nsp": args.nsp,
        "temperature_ry": args.temperature_ry,
        "electron_count": args.electron_count,
        "fermi_ry": args.fermi_ry,
        "kpoint_count": len(k_weights) if args.backend == "KS" else 1,
        "basis_size": max_index + 1,
        "metric": args.metric,
        "lmax": args.lmax,
        "norb_per_spin": args.norb_per_spin,
        "csv_sha256": sha256(csv_path),
    }
    output = args.output or args.capture / "manifest.json"
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
