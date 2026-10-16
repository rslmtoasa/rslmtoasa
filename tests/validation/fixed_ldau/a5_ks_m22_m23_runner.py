#!/usr/bin/env python3
"""Run the preregistered A5 KS U/J grid without interpreting energy outcomes."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
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


def has_converged_marker(log: str) -> bool:
    """Match the positive marker without accepting `Not converged!`."""
    return re.search(r"(?<!Not )Converged!", log, re.IGNORECASE) is not None


def fortran_true(value: str) -> bool:
    """Parse the Fortran logical spellings emitted in fixed-U residual rows."""
    return value.strip().lower() in {"true", ".true."}


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


def configure_case_deck(base_text: str, cfg: dict) -> str:
    """Build a campaign deck that explicitly selects native reciprocal KS."""
    deck = edit_group(base_text, "self", {"use_kspace": ".true."})
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
    return deck


def run_case(cfg: dict, args: argparse.Namespace, assets_hashes: dict[str, str],
             source_hashes: dict[str, str], binary_hash: str) -> dict:
    case_dir = args.output / cfg["name"]
    case_dir.mkdir(parents=False, exist_ok=False)
    deck = configure_case_deck(args.input.read_text(), cfg)
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
    converged_markers = re.findall(r"(?<!Not )Converged!", log_text, re.IGNORECASE)
    converged = has_converged_marker(log_text)
    ks_branch = "run_recursion: use_kspace=.true., skipping recursion solver stage" in log_text
    thermo_pairs = re.findall(r"entropy=\s*([^\s]+).*?TS=\s*([^\s]+)", log_text)
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    finite_thermo = bool(thermo_pairs) and all(re.fullmatch(number, x) for x in thermo_pairs[-1])
    r05_required = abs(float(cfg["u"])) > 1.0e-10 or abs(float(cfg["j"])) > 1.0e-10
    r05_lines = re.findall(r"FIXED_U_RESIDUAL ([^\r\n]*)", log_text)
    r05_values = dict(re.findall(r"([A-Za-z0-9_]+)=([^\s]+)", r05_lines[-1])) if r05_lines else {}
    r05_streak = int(r05_values["streak"]) if "streak" in r05_values else None
    r05_required_streak = int(r05_values["required_streak"]) if "required_streak" in r05_values else None
    r05_physical = fortran_true(r05_values.get("physical", ""))
    r05_ok = (not r05_required) or (r05_streak is not None and r05_required_streak is not None
                                   and r05_streak >= r05_required_streak and r05_physical)
    row = {
        **cfg, "returncode": result.returncode,
        "converged_marker_count": len(converged_markers),
        "converged": converged,
        "r05_required": r05_required, "r05_status": "PASS" if r05_ok and r05_required else
            ("NOT_APPLICABLE" if not r05_required else "MISSING_OR_NOT_SATISFIED"),
        "r05_latest_streak": r05_streak, "r05_required_streak": r05_required_streak,
        "r05_latest_physical": r05_physical if r05_required else None,
        "ks_branch": ks_branch, "finite_entropy_TS": finite_thermo,
        "case_status": "CONVERGED" if result.returncode == 0 and converged and r05_ok and ks_branch and finite_thermo else "NONCONVERGED",
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
    parser.add_argument("--fail-fast", action="store_true",
                        help="with --workers 1, stop before launching cases after the first rejected receipt")
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
    if args.fail_fast and args.workers != 1:
        parser.error("--fail-fast requires --workers 1")
    unstarted = []
    if args.fail_fast:
        rows = []
        for index, cfg in enumerate(work):
            row = run_case(cfg, args, assets_hashes, source_hashes, binary_hash)
            rows.append(row)
            if row["case_status"] != "CONVERGED":
                unstarted = [item["name"] for item in work[index + 1:]]
                break
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_case, cfg, args, assets_hashes, source_hashes, binary_hash) for cfg in work]
            rows = [future.result() for future in futures]
    complete = all(row["case_status"] == "CONVERGED" for row in rows)
    status = ("COMPLETE" if complete and not unstarted else
              "STOPPED_AFTER_REJECTED_CASE" if unstarted else "INCOMPLETE_NONCONVERGED_CASES")
    summary = {"status": status, "unstarted_cases": unstarted,
               "binary_sha256": binary_hash, "source_sha256": source_hashes,
               "case_count": len(rows), "cases": rows}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": summary["status"], "case_count": len(rows),
                      "converged": sum(r["case_status"] == "CONVERGED" for r in rows),
                      "nonzero_exit": sum(r["returncode"] != 0 for r in rows),
                      "summary": str(args.output / "summary.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
