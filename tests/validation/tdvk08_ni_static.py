#!/usr/bin/env python3
"""Run the TDVK-08 fcc-Ni direct ALSDA static diagnostic."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from tdvk06_fe_static_interactions import number, parse_headers, replace_assignment


def validate_ni_static_output(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    metadata = parse_headers(path)
    direct: list[list[str]] = []
    blocks: list[list[str]] = []
    radial: list[list[str]] = []
    gsr: list[list[str]] = []
    statuses: list[str] = []
    for line in lines:
        if not line or line.startswith("#"):
            if line.startswith("# gsr_status"):
                statuses.append(line.split("=", 1)[1].strip())
            continue
        fields = line.split()
        if len(fields) == 9 and fields[-1] == "EXECUTED":
            direct.append(fields)
        elif len(fields) == 4:
            blocks.append(fields)
        elif len(fields) == 5:
            radial.append(fields)
        elif len(fields) == 24:
            gsr.append(fields)
        else:
            raise RuntimeError(f"unexpected Ni static data row: {line}")
    if metadata.get("backend") != "static_interactions" or metadata.get("state_source") != "accepted_kspace_scf_cache":
        raise RuntimeError("Ni static backend did not certify accepted-state execution")
    if metadata.get("accepted_state_cache_reused") != "T" or metadata.get("channel") != "chi_plus":
        raise RuntimeError("Ni static accepted-state metadata is incomplete")
    if tuple(int(value) for value in metadata["actual_k_mesh"].split()) != (12, 12, 12):
        raise RuntimeError("Ni static actual mesh is not 12^3")
    if int(metadata.get("actual_k_count", "0")) != 1728 or int(metadata.get("product_dimension", "0")) != 232:
        raise RuntimeError("Ni static state/product dimensions are incomplete")
    if len(direct) != 2 or len(blocks) != 10 or len(radial) != 4940 or len(gsr) != 2 or len(statuses) != 2:
        raise RuntimeError(f"unexpected Ni static row counts: direct={len(direct)} block={len(blocks)} radial={len(radial)} gsr={len(gsr)}")
    for row in direct:
        values = [number(value) for value in row[:7]]
        if not all(value == value and abs(value) != float("inf") for value in values) or row[7].upper() != "T":
            raise RuntimeError("Ni direct ALSDA residual diagnostic is non-finite")
    if not all(status.startswith("BLOCKED:") for status in statuses):
        raise RuntimeError(f"unexpected Ni auxiliary GSR status: {statuses}")
    if "no regularization" not in metadata.get("gsr_solve_policy", ""):
        raise RuntimeError("Ni auxiliary GSR policy is not explicitly unregularized")
    if "# TDVK-06 PASS CANDIDATE" not in lines:
        raise RuntimeError("Ni static pass marker is missing")
    return {"headers": metadata, "direct_rows": direct, "gsr_rows": gsr, "gsr_status": statuses}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--database", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    template = (args.template or root / "tests/integration/tddft_driver_smoke/input_tdvk08_ni_static.nml").resolve()
    database = (args.database or root / "results/validation/TDVAL-01_FE_NI/ground_state/fccNi").resolve()
    workdir = args.scratch_root / "run"
    workdir.mkdir(parents=True, exist_ok=True)
    text = replace_assignment(template.read_text(encoding="utf-8"), "database", repr(str(database)))
    text = replace_assignment(text, "output_file", repr("tdvk08_ni_static.dat"))
    (workdir / "input.nml").write_text(text, encoding="utf-8")
    completed = subprocess.run(
        [str(args.binary.resolve()), "input.nml"], cwd=workdir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=7200, check=False,
    )
    (workdir / "tdvk08_ni_static.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"TDVK-08 Ni static executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
    result_path = workdir / "tdvk08_ni_static.dat"
    state_path = workdir / "tdvk08_ni_static.dat.state"
    if not result_path.exists() or not state_path.exists():
        raise RuntimeError("Ni static output or state artifact is missing")
    result = validate_ni_static_output(result_path)
    state_meta = parse_headers(state_path)
    if state_meta.get("direct_accepted_state_handoff") != "T":
        raise RuntimeError("Ni static state artifact does not certify direct handoff")
    result.update(material="fcc Ni", output_artifact=str(result_path), state_artifact=str(state_path), gsr_used_only_as_unregularized_auxiliary=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(result_path), "gsr_status": result["gsr_status"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
