#!/usr/bin/env python3
"""Exercise the A5 test-only HOH bypass at the real namelist boundary.

The positive case uses one SCF step and must pass Hamiltonian construction.
Negative cases must stop in input validation before any SCF iteration. This is
an execution receipt for the frozen A5 binary, not a unit test of copied logic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
CASES = {
    "positive_w1_local_u0_hoh": {
        "hoh": True, "radial": "local", "operator": "w1", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.0, "bypass": True, "expect": "A5 test-only HOH guard bypass active",
    },
    "negative_bypass_off": {
        "hoh": True, "radial": "local", "operator": "w1", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.0, "bypass": False,
        "expect": "A4 test operator requires scalar-collinear nsp=1, lmax=2, and hoh=.false..",
    },
    "negative_nonzero_u": {
        "hoh": True, "radial": "local", "operator": "w1", "lambda": 0.0,
        "u_ev": 0.1, "j_ev": 0.0, "bypass": True,
        "expect": "A5 HOH bypass requires exactly U=J=0 on every atom and shell.",
    },
    "negative_nonzero_j": {
        "hoh": True, "radial": "local", "operator": "w1", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.1, "bypass": True,
        "expect": "A5 HOH bypass requires exactly U=J=0 on every atom and shell.",
    },
    "negative_w2": {
        "hoh": True, "radial": "local", "operator": "w2", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.0, "bypass": True,
        "expect": "A5 HOH bypass is only valid for explicit general U=J=0 local W1 test runs.",
    },
    "negative_no_w_path": {
        "hoh": True, "radial": "local", "operator": "none", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.0, "bypass": True,
        "expect": "A5 HOH bypass is only valid for explicit general U=J=0 local W1 test runs.",
    },
    "negative_spectral": {
        "hoh": True, "radial": "spectral", "operator": "w1", "lambda": 0.0,
        "u_ev": 0.0, "j_ev": 0.0, "bypass": True,
        "expect": "A5 HOH bypass is only valid for explicit general U=J=0 local W1 test runs.",
    },
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


def prepare_case(name: str, cfg: dict, base_input: Path, base_seed: Path,
                 assets_dir: Path, output_root: Path, binary: Path) -> dict:
    case_dir = output_root / name
    case_dir.mkdir(parents=True, exist_ok=False)
    deck = base_input.read_text()
    deck = edit_group(deck, "self", {"nstep": "1", "use_kspace": ".true."})
    ham = {
        "hoh": ".true." if cfg["hoh"] else ".false.",
        "hubbard_u_general(1,3)": f"{cfg['u_ev']:.16g}",
        "hubbard_j_general(1,3)": f"{cfg['j_ev']:.16g}",
        "hubbard_u_potential_form": "'liechtenstein'",
        "hubbard_radial_moments": f"'{cfg['radial']}'",
        "hubbard_test_operator": f"'{cfg['operator']}'",
        "hubbard_test_lambda": f"{cfg['lambda']:.16g}",
    }
    if cfg["bypass"]:
        ham["fixed_u_test_allow_hoh_w1_u0"] = ".true."
    deck = edit_group(deck, "hamiltonian", ham)
    (case_dir / "input.nml").write_text(deck)
    input_sha256 = sha(case_dir / "input.nml")
    shutil.copyfile(base_seed, case_dir / "Fe.nml")
    seed_input_sha256 = sha(case_dir / "Fe.nml")
    for asset in ASSETS:
        shutil.copyfile(assets_dir / asset, case_dir / asset)
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"})
    result = subprocess.run([str(binary)], cwd=case_dir, env=env,
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    (case_dir / "stdout.log").write_bytes(result.stdout)
    log_path = case_dir / "testrun.log"
    log = result.stdout.decode(errors="replace")
    if log_path.is_file():
        log += "\n" + log_path.read_text(errors="replace")
    passed = cfg["expect"] in log
    ks_branch = "run_recursion: use_kspace=.true., skipping recursion solver stage" in log
    thermo_pairs = re.findall(r"entropy=\s*([^\s]+).*?TS=\s*([^\s]+)", log)
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    finite_thermo = bool(thermo_pairs) and all(
        re.fullmatch(number, value) for value in thermo_pairs[-1])
    iteration_markers = ("SCF_B0C_ITER", "HUBBARD_MOMENT_ENERGY_TOTAL")
    if name.startswith("positive_"):
        passed = passed and result.returncode == 0 and "fatal" not in log.lower() and ks_branch and finite_thermo
    else:
        passed = passed and result.returncode != 0 and not any(marker in log for marker in iteration_markers)
    row = {
        "case": name, "expected_log_text": cfg["expect"], "passed": passed,
        "returncode": result.returncode,
        "ks_branch": ks_branch, "finite_entropy_TS": finite_thermo,
        "binary_sha256": sha(binary), "input_sha256": input_sha256,
        "seed_input_sha256": seed_input_sha256,
        "stdout_sha256": sha(case_dir / "stdout.log"),
    }
    if (case_dir / "Fe.nml").is_file():
        row["Fe_nml_output_sha256"] = sha(case_dir / "Fe.nml")
    if (case_dir / "Fe_out.nml").is_file():
        row["Fe_out_nml_output_sha256"] = sha(case_dir / "Fe_out.nml")
    if log_path.is_file():
        row["testrun_log_sha256"] = sha(log_path)
    (case_dir / "receipt.json").write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    if not passed:
        raise RuntimeError(f"{name}: expected boundary behavior not observed; see {case_dir}")
    return row


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
    rows = [prepare_case(n, cfg, args.input, args.seed, args.assets, args.output, args.binary)
            for n, cfg in CASES.items()]
    summary = {"status": "PASS" if all(row["passed"] for row in rows) else "FAIL",
               "binary_sha256": sha(args.binary), "cases": rows}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": summary["status"], "case_count": len(rows),
                      "summary": str(args.output / "summary.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
