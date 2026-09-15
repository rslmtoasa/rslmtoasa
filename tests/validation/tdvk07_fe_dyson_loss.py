#!/usr/bin/env python3
"""Run and validate the TDVK-07 compact Fe Dyson/loss evidence campaign."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from pathlib import Path


HEADER_RE = re.compile(r"^#\s*([^=]+?)\s*=\s*(.*?)\s*$")


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(name)}\s*=\s*).*$")
    text, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"could not patch namelist assignment {name}")
    return text


def headers(lines: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in lines:
        match = HEADER_RE.match(line)
        if match:
            values[match.group(1).strip()] = match.group(2).strip()
    return values


def validate_output(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    meta = headers(lines)
    data = [line.split() for line in lines if line and not line.startswith("#")]
    dynamic = [row for row in data if len(row) == 15]
    static = [row for row in data if len(row) == 8]
    covariance = [row for row in data if len(row) == 7]
    gf = [row for row in data if len(row) == 11]
    raw_matrix = [row for row in data if len(row) == 10]

    if meta.get("backend") != "compact_dyson":
        raise RuntimeError("TDVK-07 compact backend metadata is missing")
    if meta.get("accepted_state_cache_reused") != "T" or meta.get("state_source") != "accepted_kspace_scf_cache":
        raise RuntimeError("TDVK-07 did not certify accepted-state reuse")
    if tuple(int(value) for value in meta["accepted_k_mesh"].split()) != (12, 12, 12):
        raise RuntimeError("TDVK-07 accepted mesh is not 12^3")
    if int(meta.get("accepted_k_count", "0")) != 1728:
        raise RuntimeError("TDVK-07 accepted k-point count is not 1728")
    if meta.get("interaction_route") != "direct_alsda":
        raise RuntimeError("TDVK-07 did not use direct ALSDA")
    if meta.get("BES_GCR_goldstone_correction") != "OFF":
        raise RuntimeError("TDVK-07 Goldstone correction policy changed")
    if int(meta.get("frequency_count", "0")) != 4 or int(meta.get("q_count", "0")) != 3:
        raise RuntimeError("TDVK-07 q/frequency grid metadata is incomplete")
    if len(static) != 2 or len(dynamic) != 12 or len(covariance) != 4 or len(gf) != 1:
        raise RuntimeError(f"unexpected TDVK-07 row counts: static={len(static)} dynamic={len(dynamic)} covariance={len(covariance)} gf={len(gf)}")
    expected_q = {1: (0.0, 0.0, 0.0), 2: (1.0 / 12.0, 0.0, 0.0), 3: (-1.0 / 12.0, 0.0, 0.0)}
    expected_omega = (0.0, 0.01, 0.02, 0.05)
    for row in dynamic:
        q_index = int(row[0])
        q = tuple(number(value) for value in row[1:4])
        omega = number(row[4])
        if q_index not in expected_q or max(abs(a - b) for a, b in zip(q, expected_q[q_index])) > 1e-12:
            raise RuntimeError(f"unexpected literal q in dynamic row: {row}")
        if min(abs(omega - target) for target in expected_omega) > 1e-12:
            raise RuntimeError(f"unexpected omega in dynamic row: {row}")
        values = [number(value) for value in row[5:]]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("non-finite dynamic diagnostic")
        if values[-3] > 1e-8 or values[-1] > 1e-8:
            raise RuntimeError("uncontrolled compact Dyson residual")
    for row in static:
        values = [number(value) for value in row]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("non-finite static denominator diagnostic")
    for row in covariance:
        values = [number(value) for value in row[2:]]
        if not all(math.isfinite(value) for value in values) or values[1] > 5e-8 or values[2] > 5e-8:
            raise RuntimeError("interacting q/-q covariance failed")
    for row in gf:
        values = [number(value) for value in row[1:]]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("non-finite GF spot-check diagnostic")
        if int(values[4]) != 6401 or values[6] >= 1.0:
            raise RuntimeError("GF spot-check quadrature is not TDVK-03 controlled")
    if len(raw_matrix) != 3 * 4 * 232 * 232:
        raise RuntimeError(f"complete raw matrix archive is incomplete: {len(raw_matrix)} rows")
    if "# TDVK-07 PASS CANDIDATE" not in lines:
        raise RuntimeError("TDVK-07 PASS CANDIDATE marker is missing")
    return {
        "output": str(path),
        "static_rows": static,
        "dynamic_rows": dynamic,
        "covariance_rows": covariance,
        "gf_rows": gf,
        "raw_matrix_rows": len(raw_matrix),
        "headers": meta,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = root / "tests/integration/tddft_driver_smoke/input_tdvk07_fe.nml"
    fe_database = root / "tests/scf/cases/bulk/bccFe"
    workdir = args.scratch_root / "run"
    workdir.mkdir(parents=True, exist_ok=True)
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(fe_database)))
    text = replace_assignment(text, "output_file", repr("tdvk07_fe.dat"))
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
    (workdir / "tdvk07.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"TDVK-07 executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
    result_path = workdir / "tdvk07_fe.dat"
    state_path = workdir / "tdvk07_fe.dat.state"
    if not result_path.exists() or not state_path.exists():
        raise RuntimeError("TDVK-07 output or state artifact is missing")
    result = validate_output(result_path)
    state_meta = headers(state_path.read_text(encoding="utf-8", errors="replace").splitlines())
    if state_meta.get("direct_accepted_state_handoff") != "T":
        raise RuntimeError("TDVK-07 state artifact does not certify direct handoff")
    result["state_artifact"] = str(state_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(result_path), "raw_matrix_rows": result["raw_matrix_rows"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
