#!/usr/bin/env python3
"""Run and validate the TDVK-08 fcc-Ni finite-q bare-response seam."""

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


def parse_headers(lines: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in lines:
        match = HEADER_RE.match(line)
        if match:
            result[match.group(1).strip()] = match.group(2).strip()
    return result


def validate(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    meta = parse_headers(lines)
    rows = [line.split() for line in lines if line and not line.startswith("#")]
    lehmann = [row for row in rows if len(row) == 10]
    endpoint = [row for row in rows if len(row) == 18]
    covariance = [row for row in rows if len(row) == 6]
    gf = [row for row in rows if len(row) == 14]
    if meta.get("backend") != "product_finite_q":
        raise RuntimeError("Ni finite-q backend metadata is missing")
    if int(meta.get("product_dimension", "0")) != 232 or int(meta.get("accepted_state_nk", "0")) != 1728:
        raise RuntimeError("Ni finite-q response does not use complete 12^3 compact state")
    if len(lehmann) != 4 or len(endpoint) != 4 or len(covariance) != 1 or len(gf) != 2:
        raise RuntimeError(f"unexpected finite-q row counts: lehmann={len(lehmann)} endpoint={len(endpoint)} covariance={len(covariance)} gf={len(gf)}")
    expected_q = [(0.0, 0.0, 0.0), (1.0 / 12.0, 0.0, 0.0), (-1.0 / 12.0, 0.0, 0.0), (0.23, 0.07, -0.11)]
    for row, expected in zip(lehmann, expected_q):
        q = tuple(number(value) for value in row[1:4])
        if max(abs(a - b) for a, b in zip(q, expected)) > 1.0e-12:
            raise RuntimeError(f"unexpected literal q: {q}")
        if not all(math.isfinite(number(value)) for value in row[4:9]):
            raise RuntimeError("non-finite finite-q Lehmann diagnostic")
    for row in endpoint:
        values = [number(value) for value in row[1:-1]]
        if not all(math.isfinite(value) for value in values) or abs(values[6]) > 2.0e-11:
            raise RuntimeError("finite-q endpoint folding provenance failed")
    covariance_values = [number(value) for value in covariance[0][2:4]]
    if covariance[0][0:2] != ["2", "3"] or not all(math.isfinite(value) for value in covariance_values):
        raise RuntimeError("finite-q covariance pair is missing or non-finite")
    if covariance_values[1] > 5.0e-8:
        raise RuntimeError("finite-q bare q/-q covariance failed")
    for row in gf:
        values = [number(value) for value in row[4:13]]
        if not all(math.isfinite(value) for value in values) or int(row[4]) % 2 == 0 or row[13].upper() != "T":
            raise RuntimeError("finite-q auxiliary GF row is not finite")
    return {
        "material": "fcc Ni",
        "output": str(path),
        "headers": meta,
        "q_coordinates": [list(q) for q in expected_q],
        "lehmann_rows": lehmann,
        "endpoint_rows": endpoint,
        "covariance_rows": covariance,
        "gf_rows": gf,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--database", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    template = (args.template or root / "tests/integration/tddft_driver_smoke/input_tdvk08_ni_finite_q.nml").resolve()
    database = (args.database or root / "results/validation/TDVAL-01_FE_NI/ground_state/fccNi").resolve()
    workdir = args.scratch_root / "run"
    workdir.mkdir(parents=True, exist_ok=True)
    text = replace_assignment(template.read_text(encoding="utf-8"), "database", repr(str(database)))
    text = replace_assignment(text, "output_file", repr("tdvk08_ni_finite_q.dat"))
    (workdir / "input.nml").write_text(text, encoding="utf-8")
    completed = subprocess.run(
        [str(args.binary.resolve()), "input.nml"], cwd=workdir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=7200, check=False,
    )
    (workdir / "tdvk08_ni_finite_q.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"TDVK-08 Ni finite-q executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
    result_path = workdir / "tdvk08_ni_finite_q.dat"
    state_path = workdir / "tdvk08_ni_finite_q.dat.state"
    if not result_path.exists() or not state_path.exists():
        raise RuntimeError("Ni finite-q output or state artifact is missing")
    result = validate(result_path)
    state_meta = parse_headers(state_path.read_text(encoding="utf-8", errors="replace").splitlines())
    if state_meta.get("direct_accepted_state_handoff") != "T" or state_meta.get("reciprocal_mode") != "ham_only":
        raise RuntimeError("Ni finite-q state artifact does not certify direct ham_only handoff")
    result["state_artifact"] = str(state_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(result_path), "covariance_residual": result["covariance_rows"][0][3]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
