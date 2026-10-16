#!/usr/bin/env python3
"""One-step native-KS check for the A5 moments-only fixed-U capture path."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
REQUIRED = {"test_density_matrix", "test_q0_trace", "test_trace_wp", "ql", "radial_moment"}
FORBIDDEN = {
    "h_total", "eigenvalue", "occupation", "eigenvector", "s_overlap", "basis_map",
    "green_function", "fermi_step", "rs_hopping", "rs_green", "hamiltonian",
}


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
    kept = []
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.binary, args.input, args.seed, *(args.assets / a for a in ASSETS)):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=False)
    deck = edit_group(args.input.read_text(), "self", {"nstep": "1", "use_kspace": ".true."})
    deck = edit_group(deck, "hamiltonian", {
        "hubbard_u_general(1,3)": "0.0", "hubbard_j_general(1,3)": "0.0",
        "hubbard_u_potential_form": "'liechtenstein'", "hubbard_radial_moments": "'local'",
        "hubbard_test_operator": "'w1'", "hubbard_test_lambda": "0.0",
        "hubbard_test_spin": "0", "hoh": ".false.",
    })
    input_path = args.output / "input.nml"
    seed_path = args.output / "Fe.nml"
    input_path.write_text(deck)
    shutil.copyfile(args.seed, seed_path)
    input_hash = sha(input_path)
    seed_hash = sha(seed_path)
    for asset in ASSETS:
        shutil.copyfile(args.assets / asset, args.output / asset)
    trace_path = args.output / "moments_trace.csv"
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1", "RSLMTO_FIXED_U_TRACE": str(trace_path),
        "RSLMTO_FIXED_U_TRACE_ITERATION": "0", "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
        "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1", "RSLMTO_FIXED_U_LAMBDA": "0",
        "RSLMTO_FIXED_U_RUN_ID": "a5_moments_only_smoke",
    })
    result = subprocess.run([str(args.binary)], cwd=args.output, env=env,
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    stdout = args.output / "stdout.log"
    stdout.write_bytes(result.stdout)
    if result.returncode != 0 or not trace_path.is_file():
        raise RuntimeError(f"native KS smoke failed, returncode={result.returncode}; see {stdout}")
    log = result.stdout.decode(errors="replace")
    if "run_recursion: use_kspace=.true., skipping recursion solver stage" not in log:
        raise RuntimeError("native KS routing marker absent; this smoke cannot certify the KS path")
    thermo_pairs = re.findall(r"entropy=\s*([^\s]+).*?TS=\s*([^\s]+)", log)
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    if not thermo_pairs or not all(re.fullmatch(number, value) for value in thermo_pairs[-1]):
        raise RuntimeError("finite entropy/TS marker absent; this smoke cannot certify the KS path")
    with trace_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    kinds = {row["kind"] for row in rows}
    missing = sorted(REQUIRED - kinds)
    present_forbidden = sorted((FORBIDDEN & kinds) | {kind for kind in kinds if kind.startswith("radial_replay_")})
    p_rows = [row for row in rows if row["kind"] == "test_density_matrix"]
    q_rows = [row for row in rows if row["kind"] == "test_q0_trace"]
    wp_rows = [row for row in rows if row["kind"] == "test_trace_wp"]
    sites = {row["site"] for row in p_rows}
    expected_sites_spins = {(site, spin) for site in sites for spin in ("0", "1")}
    p_counts = {}
    for row in p_rows:
        key = (row["site"], row["spin"])
        p_counts[key] = p_counts.get(key, 0) + 1
    q_pairs = {(row["site"], row["spin"]) for row in q_rows}
    wp_pairs = {(row["site"], row["spin"]) for row in wp_rows}
    p_out_ok = (bool(sites) and all(count == 25 for count in p_counts.values()) and
                set(p_counts) == expected_sites_spins and
                {row["stage"] for row in p_rows} == {"P_OUT"} and
                {row["spin"] for row in p_rows} == {"0", "1"} and
                {row["stage"] for row in q_rows + wp_rows} == {"P_OUT"} and
                q_pairs == expected_sites_spins and wp_pairs == expected_sites_spins)
    if missing or present_forbidden or not p_out_ok:
        raise RuntimeError(f"trace contract failed; missing={missing}; forbidden={present_forbidden}; "
                           f"P_OUT_spin_capture={p_out_ok}; P/q0/WP_rows={len(p_rows)}/{len(q_rows)}/{len(wp_rows)}")
    receipt = {
        "status": "PASS", "returncode": result.returncode,
        "ks_branch": True, "finite_entropy_TS": True,
        "binary_sha256": sha(args.binary), "input_sha256": input_hash,
        "seed_input_sha256": seed_hash, "stdout_sha256": sha(stdout),
        "trace_output_sha256": sha(trace_path), "trace_bytes": trace_path.stat().st_size,
        "required_kinds": sorted(REQUIRED), "observed_kinds": sorted(kinds),
        "forbidden_kinds_absent": sorted(FORBIDDEN) + ["radial_replay_*"],
        "runtime_Fe_nml_sha256": sha(seed_path),
        "runtime_Fe_out_nml_sha256": sha(args.output / "Fe_out.nml") if (args.output / "Fe_out.nml").is_file() else None,
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
