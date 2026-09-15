#!/usr/bin/env python3
"""Run and validate the TDVK-06 Fe static interaction diagnostics.

This is an evidence campaign rather than a default quick test.  It launches
one fresh self-consistent 12^3 k-space SCF process, consumes that accepted
cache in the static-only TDDFT backend, and accepts either PASS or an honest
GSR BLOCKED status.  A blocked GSR solve is evidence; it is never repaired by
the harness.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from pathlib import Path


FLOAT = r"[-+0-9.EeDd]+"
HEADER_RE = re.compile(r"^#\s*([^=]+?)\s*=\s*(.*?)\s*$")


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(name)}\s*=\s*).*$")
    text, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"could not patch namelist assignment {name}")
    return text


def parse_headers(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = HEADER_RE.match(line)
        if match:
            headers[match.group(1).strip()] = match.group(2).strip()
    return headers


def validate_output(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(path)
    direct_rows: list[list[str]] = []
    block_rows: list[list[str]] = []
    gsr_rows: list[list[str]] = []
    statuses: list[str] = []
    for line in lines:
        if not line or line.startswith("#"):
            if line.startswith("# gsr_status"):
                statuses.append(line.split("=", 1)[1].strip())
            continue
        fields = line.split()
        if len(fields) == 9 and fields[-1] == "EXECUTED":
            direct_rows.append(fields)
        elif len(fields) == 4:
            block_rows.append(fields)
        elif len(fields) == 5:
            continue
        elif len(fields) == 24 and fields[-1] in {"T", "F"}:
            gsr_rows.append(fields)
        else:
            raise RuntimeError(f"unexpected TDVK-06 data row: {line}")

    if headers.get("backend") != "static_interactions":
        raise RuntimeError("static interaction backend metadata is missing")
    if headers.get("state_source") != "accepted_kspace_scf_cache":
        raise RuntimeError("TDVK-06 did not consume the accepted k-space SCF cache")
    if headers.get("accepted_state_cache_reused") != "T":
        raise RuntimeError("accepted-state cache reuse was not certified")
    if tuple(int(value) for value in headers["actual_k_mesh"].split()) != (12, 12, 12):
        raise RuntimeError("TDVK-06 actual mesh is not 12^3")
    if int(headers.get("actual_k_count", "0")) != 1728:
        raise RuntimeError("TDVK-06 actual k-point count is not 1728")
    if int(headers.get("product_dimension", "0")) != 232:
        raise RuntimeError("TDVK-06 compact product dimension is not 232")
    if headers.get("channel") != "chi_plus":
        raise RuntimeError("TDVK-06 channel is not chi_plus")
    if headers.get("pauli_magnetization_label") != "pauli_projected":
        raise RuntimeError("TDVK-06 Pauli magnetization label is not pauli_projected")
    for key in (
        "pauli_magnetization_weighted_norm_m00",
        "pauli_magnetization_projection_residual_norm_m00",
        "pauli_magnetization_projection_relative_residual_m00",
        "pauli_valence_projection_relative_residual_m00",
        "pauli_core_projection_relative_residual_m00",
    ):
        if not math.isfinite(number(headers.get(key, "nan"))):
            raise RuntimeError(f"missing or non-finite magnetization projection diagnostic: {key}")
    if not lines or "# TDVK-06 PASS CANDIDATE" not in lines:
        raise RuntimeError("TDVK-06 PASS CANDIDATE marker is missing")
    if not any(line.startswith("# TDVK-06 GSR CLOSURE PASS CANDIDATE") for line in lines):
        raise RuntimeError("TDVK-06 GSR closure candidate marker is missing")
    radial_rows = [line.split() for line in lines if line and not line.startswith("#") and len(line.split()) == 5]
    if len(direct_rows) != 2 or len(gsr_rows) != 2 or len(block_rows) != 10 or len(radial_rows) != 4940 or len(statuses) != 2:
        raise RuntimeError(
            f"unexpected diagnostic row counts: direct={len(direct_rows)} block={len(block_rows)} "
            f"radial={len(radial_rows)} gsr={len(gsr_rows)} status={len(statuses)}"
        )
    direct_values = [[number(value) for value in row[:7]] for row in direct_rows]
    if not all(math.isfinite(value) for row in direct_values for value in row):
        raise RuntimeError("direct ALSDA residual diagnostics are not finite")
    for row in gsr_rows:
        if int(row[1]) != 232:
            raise RuntimeError("GSR compact dimension is inconsistent with the production product space")
        values = [number(value) for index, value in enumerate(row) if index not in {5, 23}]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("GSR diagnostic contains a non-finite value")
        if number(row[17]) > 1.0e-8:
            raise RuntimeError("GSR solve and reconstructed residuals are not mutually consistent")
    if not all(status.startswith(("PASS:", "BLOCKED:")) for status in statuses):
        raise RuntimeError(f"unexpected GSR status: {statuses}")
    return {
        "headers": headers,
        "direct_rows": direct_rows,
        "gsr_rows": gsr_rows,
        "gsr_status": statuses,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = root / "tests/integration/tddft_driver_smoke/input_tdvk06_fe.nml"
    fe_database = root / "tests/scf/cases/bulk/bccFe"
    workdir = args.scratch_root / "run"
    workdir.mkdir(parents=True, exist_ok=True)
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(fe_database)))
    text = replace_assignment(text, "output_file", repr("tdvk06_static.dat"))
    input_path = workdir / "input.nml"
    input_path.write_text(text, encoding="utf-8")

    completed = subprocess.run(
        [str(args.binary), input_path.name],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=7200,
        check=False,
    )
    (workdir / "tdvk06.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"TDVK-06 executable failed with {completed.returncode}:\n{completed.stdout[-8000:]}")
    result_path = workdir / "tdvk06_static.dat"
    if not result_path.exists():
        raise RuntimeError(f"TDVK-06 output was not written: {result_path}")
    result = validate_output(result_path)
    state_path = workdir / "tdvk06_static.dat.state"
    if not state_path.exists():
        raise RuntimeError("TDVK-06 state artifact is missing")
    state_headers = parse_headers(state_path)
    if state_headers.get("direct_accepted_state_handoff") != "T":
        raise RuntimeError("TDVK-06 state artifact does not certify direct handoff")
    result["state_artifact"] = str(state_path)
    result["output_artifact"] = str(result_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(result_path), "gsr_status": result["gsr_status"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
