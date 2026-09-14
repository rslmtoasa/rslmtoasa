#!/usr/bin/env python3
"""Run the TDVK-05 compact Fe bare-response convergence campaign.

This is an evidence campaign, not a default quick test.  Each case starts
from the same physical Fe input and performs its own SCF run.  The response
driver then evaluates the compact Lehmann service on the requested q/omega
grid.  Only the selected 8^3 reference case performs the prescribed finite-q
reciprocal-GF spot checks.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


FLOAT = r"[-+0-9.EeDd]+"
HEADER_RE = re.compile(rf"^#\s*([^=]+?)\s*=\s*(.*?)\s*$")
METRIC_COLUMNS = (
    "eta_index",
    "eta_Ry",
    "q_index",
    "qx",
    "qy",
    "qz",
    "omega_Ry",
    "frobenius_norm",
    "max_abs_element",
    "trace_real",
    "trace_imag",
    "transitions_evaluated",
    "occupation_skips",
    "runtime_cpu_seconds",
    "finite",
)
GF_COLUMNS = (
    "q_index",
    "qx",
    "qy",
    "qz",
    "integration_points",
    "integration_eta",
    "h_over_integration_eta",
    "norm_lehmann",
    "norm_gf",
    "dF",
    "rF",
    "dInf",
    "wall_seconds",
    "finite",
)


@dataclass(frozen=True)
class Case:
    name: str
    mesh: int
    response_lmax: int
    eta_values: tuple[float, ...]
    backend: str = "product_convergence"


CASES = (
    Case("mesh4_full", 4, -1, (0.01,)),
    Case("mesh8_full_eta_ladder", 8, -1, (0.02, 0.01, 0.005)),
    Case("mesh12_full", 12, -1, (0.01,)),
    Case("mesh8_reduced_lmax2", 8, 2, (0.01,)),
)


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(name)}\s*=\s*).*$")
    text, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise RuntimeError(f"could not patch namelist assignment {name}")
    return text


def prepare_input(template: Path, destination: Path, fe_database: Path, case: Case, gf: bool) -> None:
    text = template.read_text(encoding="utf-8")
    text = replace_assignment(text, "database", repr(str(fe_database)))
    text = replace_assignment(text, "nk1", str(case.mesh))
    text = replace_assignment(text, "nk2", str(case.mesh))
    text = replace_assignment(text, "nk3", str(case.mesh))
    text = replace_assignment(text, "response_lmax", str(case.response_lmax))
    text = replace_assignment(text, "backend", repr(case.backend))
    text = replace_assignment(text, "n_eta", str(len(case.eta_values)))
    text = replace_assignment(text, "eta_grid", ", ".join(f"{eta:.12g}" for eta in case.eta_values))
    text = replace_assignment(text, "eta", f"{case.eta_values[0]:.12g}")
    if gf:
        text = replace_assignment(text, "n_q", "4")
        text = replace_assignment(text, "q_list", "0.0, 0.0, 0.0, 0.125, 0.0, 0.0, -0.125, 0.0, 0.0, 0.23, 0.07, -0.11")
        text = replace_assignment(text, "n_omega", "1")
        text = replace_assignment(text, "eta", "0.01")
        text = replace_assignment(text, "n_eta", "1")
        text = replace_assignment(text, "eta_grid", "0.01")
        text = replace_assignment(text, "backend", repr("product_finite_q"))
        text = replace_assignment(text, "output_file", repr("tddft_tdvk05_gf_spot.dat"))
    else:
        text = replace_assignment(text, "n_q", "2")
        text = replace_assignment(text, "q_list", "0.0, 0.0, 0.0, 0.125, 0.0, 0.0")
        text = replace_assignment(text, "n_omega", "2")
        text = replace_assignment(text, "output_file", repr("tddft_tdvk05_fe.dat"))
    destination.write_text(text, encoding="utf-8")


def parse_headers(lines: list[str]) -> dict[str, str]:
    headers: dict[str, str] = {}
    for line in lines:
        match = HEADER_RE.match(line)
        if match:
            headers[match.group(1).strip()] = match.group(2).strip()
    return headers


def parse_convergence_output(path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    rows: list[dict[str, object]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != len(METRIC_COLUMNS):
            raise RuntimeError(f"unexpected TDVK-05 metric row in {path}: {line}")
        row: dict[str, object] = {}
        for index, column in enumerate(METRIC_COLUMNS):
            if column in {"eta_index", "q_index", "transitions_evaluated", "occupation_skips"}:
                row[column] = int(fields[index])
            elif column == "finite":
                row[column] = fields[index].upper() == "T"
            else:
                row[column] = number(fields[index])
        rows.append(row)
    if not rows:
        raise RuntimeError(f"TDVK-05 output contains no metric rows: {path}")
    if not all(bool(row["finite"]) for row in rows):
        raise RuntimeError(f"TDVK-05 output contains a non-finite response: {path}")
    return headers, rows


def parse_gf_output(path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    headers = parse_headers(lines)
    rows: list[dict[str, object]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) == len(GF_COLUMNS):
            row: dict[str, object] = {}
            for index, column in enumerate(GF_COLUMNS):
                if column in {"q_index", "integration_points"}:
                    row[column] = int(fields[index])
                elif column == "finite":
                    row[column] = fields[index].upper() == "T"
                else:
                    row[column] = number(fields[index])
            rows.append(row)
    if not rows:
        raise RuntimeError(f"TDVK-05 GF spot output contains no GF metric rows: {path}")
    if not all(bool(row["finite"]) for row in rows):
        raise RuntimeError(f"TDVK-05 GF spot output contains a non-finite response: {path}")
    return headers, rows


def run_case(binary: Path, runner: Path, template: Path, fe_database: Path, scratch_root: Path, case: Case, gf: bool) -> dict[str, object]:
    workdir = scratch_root / case.name
    if gf:
        workdir = scratch_root / "gf_spot_reference"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    input_path = workdir / "input.nml"
    prepare_input(template, input_path, fe_database, case, gf)
    started = time.monotonic()
    completed = subprocess.run(
        ["bash", str(runner), str(binary.resolve())],
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
    )
    elapsed = time.monotonic() - started
    log = (workdir / "testrun.log").read_text(encoding="utf-8", errors="replace") if (workdir / "testrun.log").exists() else ""
    if completed.returncode != 0:
        raise RuntimeError(f"{case.name}: executable failed with rc={completed.returncode}\n{log[-4000:]}")
    output_name = "tddft_tdvk05_gf_spot.dat" if gf else "tddft_tdvk05_fe.dat"
    output_path = workdir / output_name
    if not output_path.exists():
        raise RuntimeError(f"{case.name}: missing {output_name}")
    if gf:
        headers, rows = parse_gf_output(output_path)
    else:
        headers, rows = parse_convergence_output(output_path)
    if "Converged!" not in log:
        raise RuntimeError(f"{case.name}: SCF did not report convergence")
    result: dict[str, object] = {
        "name": case.name,
        "mesh": case.mesh,
        "response_lmax_requested": case.response_lmax,
        "backend": "product_finite_q" if gf else case.backend,
        "elapsed_wall_seconds": elapsed,
        "headers": headers,
        "rows": rows,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-gf", action="store_true", help="skip the single reference reciprocal-GF spot run")
    parser.add_argument("--gf-only", action="store_true", help="run only the single reference reciprocal-GF spot case")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = root / "tests/integration/tddft_driver_smoke/input_tdvk05_fe.nml"
    fe_database = root / "tests/scf/cases/bulk/bccFe"
    runner = root / "tests/run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    if args.gf_only:
        results = [run_case(args.binary, runner, template, fe_database, args.scratch_root, Case("gf_spot_reference", 8, -1, (0.01,)), True)]
    else:
        results = [run_case(args.binary, runner, template, fe_database, args.scratch_root, case, False) for case in CASES]
    if not args.no_gf and not args.gf_only:
        results.append(run_case(args.binary, runner, template, fe_database, args.scratch_root, Case("gf_spot_reference", 8, -1, (0.01,)), True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"TDVK-05 cases={len(results)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
