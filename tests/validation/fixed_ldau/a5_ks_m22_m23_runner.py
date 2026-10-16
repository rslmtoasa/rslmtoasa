#!/usr/bin/env python3
"""Run the preregistered A5 KS U/J grid without interpreting energy outcomes."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
U_GRID = (1.6, 1.8, 1.9, 1.95, 2.0, 2.05, 2.1, 2.2, 2.4, 4.0)
TIGHT_U = (1.95, 2.0, 2.05)
J_GRID = (0.4, 0.45, 0.5, 0.55, 0.6, 0.8, 0.85, 0.9, 0.95, 1.0)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def edit_group(deck: str, group: str, assignments: dict[str, str]) -> str:
    lines = deck.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip().lower() == f"&{group}"), None)
    if start is None:
        lines.extend([f"&{group}", "/"])
        start = len(lines) - 2
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == "/"), None)
    if end is None:
        raise ValueError(f"unterminated &{group} namelist")
    remaining = dict(assignments)
    kept: list[str] = []
    for line in lines[start + 1:end]:
        stripped = line.lstrip()
        key = stripped.split("=", 1)[0].strip().lower() if "=" in stripped and not stripped.startswith("!") else ""
        if key in remaining:
            kept.append(f"{key} = {remaining.pop(key)}")
        else:
            kept.append(line)
    kept.extend(f"{key} = {value}" for key, value in remaining.items())
    lines[start + 1:end] = kept
    return "\n".join(lines) + "\n"


def cases() -> list[dict]:
    out = []
    for mode in ("spectral", "local"):
        for u in U_GRID:
            out.append({"name": f"u{u:g}_{mode}_standard", "u": u, "j": 0.0,
                        "mode": mode, "tight": False})
        for j in J_GRID:
            out.append({"name": f"u2_j{j:g}_{mode}_standard", "u": 2.0, "j": j,
                        "mode": mode, "tight": False})
        for u in TIGHT_U:
            out.append({"name": f"u{u:g}_{mode}_tight", "u": u, "j": 0.0,
                        "mode": mode, "tight": True})
    for mode in ("spectral", "local"):
        out.append({"name": f"u0_{mode}_standard", "u": 0.0, "j": 0.0,
                    "mode": mode, "tight": False})
    return out


def run_case(cfg: dict, args: argparse.Namespace, assets_hashes: dict[str, str],
             source_hashes: dict[str, str], binary_hash: str) -> dict:
    case_dir = args.output / cfg["name"]
    case_dir.mkdir(parents=False, exist_ok=False)
    deck = args.input.read_text()
    deck = edit_group(deck, "hamiltonian", {
        "hubbard_u_general(1,3)": f"{cfg['u']:.16g}",
        "hubbard_j_general(1,3)": f"{cfg['j']:.16g}",
        "hubbard_u_potential_form": "'liechtenstein'",
        "hubbard_radial_moments": f"'{cfg['mode']}'",
        "hubbard_test_operator": "'w1'",
        "hubbard_test_lambda": "0.0",
        "hubbard_test_spin": "0",
        "fixed_u_test_tighten_r05": ".true." if cfg["tight"] else ".false.",
    })
    if cfg["tight"]:
        deck = edit_group(deck, "self", {"conv_thr": "1.0d-10"})
    input_path = case_dir / "input.nml"
    input_path.write_text(deck)
    seed_path = case_dir / "Fe_out.nml"
    shutil.copyfile(args.seed, seed_path)
    for asset in ASSETS:
        shutil.copyfile(args.assets / asset, case_dir / asset)
    input_sha = sha(input_path)
    seed_sha = sha(seed_path)
    log_path = case_dir / "stdout.log"
    trace_path = case_dir / "records.csv"
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1", "RSLMTO_FIXED_U_TRACE": str(trace_path),
        "RSLMTO_FIXED_U_TRACE_RUN_ID": cfg["name"], "RSLMTO_FIXED_U_LAMBDA": "0",
        "RSLMTO_FIXED_U_TRACE_ITERATION": "0", "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
        "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1",
    })
    with log_path.open("wb") as log:
        result = subprocess.run([str(args.binary)], cwd=case_dir, env=env,
                                stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                check=False)
    log_text = log_path.read_text(errors="replace")
    outputs = {p.name: sha(p) for p in case_dir.iterdir() if p.is_file() and
               p.name not in {"input.nml", "Fe_out.nml"} and
               (p.name in {"stdout.log", "records.csv", "Fe_scf.nml", "Fe.nml"} or p.name.endswith("_out.nml"))}
    row = {
        **cfg, "returncode": result.returncode,
        "converged_marker_count": log_text.count("Converged!"),
        "input_sha256": input_sha, "seed_input_sha256": seed_sha,
        "runtime_Fe_out_sha256": sha(seed_path) if seed_path.exists() else None,
        "output_sha256": outputs,
        "binary_sha256": binary_hash, "source_sha256": source_hashes,
        "asset_input_sha256": assets_hashes,
    }
    (case_dir / "receipt.json").write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True, help="accepted common-state input.nml")
    parser.add_argument("--seed", type=Path, required=True, help="immutable accepted Fe_out.nml")
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--case-names", nargs="*", help="optional exact subset of preregistered case names")
    args = parser.parse_args()
    for path in (args.binary, args.input, args.seed, *(args.assets / name for name in ASSETS)):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=False)
    source_paths = [
        "source/hamiltonian.f90", "source/hamiltonian_build.f90",
        "source/hamiltonian_hubbard.f90", "source/include_codes/namelists/hamiltonian.f90",
        "source/self.f90",
    ]
    source_hashes = {path: sha(Path(path)) for path in source_paths}
    assets_hashes = {name: sha(args.assets / name) for name in ASSETS}
    binary_hash = sha(args.binary)
    work = cases()
    if args.case_names is not None:
        selected = set(args.case_names)
        if len(selected) != len(args.case_names):
            parser.error("--case-names must not contain duplicates")
        known = {case["name"] for case in work}
        unknown = selected - known
        if unknown:
            parser.error(f"unknown case name(s): {', '.join(sorted(unknown))}")
        by_name = {case["name"]: case for case in work}
        work = [by_name[name] for name in args.case_names]
        if not work:
            parser.error("--case-names selected no cases")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_case, cfg, args, assets_hashes, source_hashes, binary_hash) for cfg in work]
        rows = [future.result() for future in futures]
    summary = {"binary_sha256": binary_hash, "source_sha256": source_hashes,
               "case_count": len(rows), "cases": rows}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"case_count": len(rows), "converged": sum(r["converged_marker_count"] > 0 for r in rows),
                      "nonzero_exit": sum(r["returncode"] != 0 for r in rows),
                      "summary": str(args.output / "summary.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
