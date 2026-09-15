#!/usr/bin/env python3
"""Run and validate the TDVK-08 fcc-Ni compact Dyson/loss campaign.

The parser and numerical archive checks are shared with the accepted TDVK-07
workflow; this wrapper changes only the material template and database.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from tdvk07_fe_dyson_loss import headers, replace_assignment, validate_output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--database", type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = (args.template or root / "tests/integration/tddft_driver_smoke/input_tdvk08_ni_dyson.nml").resolve()
    database = (args.database or root / "results/validation/TDVAL-01_FE_NI/ground_state/fccNi").resolve()
    workdir = args.scratch_root / "run"
    workdir.mkdir(parents=True, exist_ok=True)
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(database)))
    text = replace_assignment(text, "output_file", repr("tdvk08_ni.dat"))
    input_path = workdir / "input.nml"
    input_path.write_text(text, encoding="utf-8")
    completed = subprocess.run(
        [str(args.binary.resolve()), input_path.name],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=7200,
        check=False,
    )
    (workdir / "tdvk08_ni.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"TDVK-08 Ni executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
    result_path = workdir / "tdvk08_ni.dat"
    state_path = workdir / "tdvk08_ni.dat.state"
    if not result_path.exists() or not state_path.exists():
        raise RuntimeError("TDVK-08 Ni output or state artifact is missing")
    result = validate_output(result_path)
    state_meta = headers(state_path.read_text(encoding="utf-8", errors="replace").splitlines())
    if state_meta.get("direct_accepted_state_handoff") != "T":
        raise RuntimeError("TDVK-08 Ni state artifact does not certify direct handoff")
    if state_meta.get("reciprocal_mode") != "ham_only":
        raise RuntimeError("TDVK-08 Ni response is not ham_only")
    if state_meta.get("kspace_hamiltonian_order") != "second":
        raise RuntimeError("TDVK-08 Ni response is not second-order")
    result["material"] = "fcc Ni"
    result["input_template"] = str(template)
    result["database"] = str(database)
    result["state_artifact"] = str(state_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(result_path), "raw_matrix_rows": result["raw_matrix_rows"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
