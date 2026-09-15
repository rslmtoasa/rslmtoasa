#!/usr/bin/env python3
"""Black-box acceptance tests for the generic compact-Dyson production path.

The positive cases are fresh executions of the real bcc-Fe k-space-SCF ->
TDDFT workflow.  They deliberately use small omega grids and an 8^3 accepted
mesh so the test exercises both the production q-list contract and the
non-12^3 runtime path without becoming a material-validation campaign.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
from pathlib import Path

from tdvk07_fe_dyson_loss import headers, number, replace_assignment


CASES = {
    "positive_q_only": {
        "q": ((0.0, 0.0, 0.0), (0.03, 0.0, 0.0), (0.06, 0.0, 0.0), (0.09, 0.0, 0.0)),
        "static": ".false.",
        "covariance": ".false.",
    },
    "no_gamma": {
        "q": ((0.03, 0.0, 0.0), (0.06, 0.0, 0.0)),
        "static": ".false.",
        "covariance": ".false.",
    },
    "covariance_pair": {
        # The covariance harness keeps the established commensurate finite-q
        # point on the 8^3 mesh; production-only cases above intentionally use
        # arbitrary positive q values.
        "q": ((0.125, 0.0, 0.0), (-0.125, 0.0, 0.0)),
        "static": ".false.",
        "covariance": ".true.",
    },
}


def q_literal(points: tuple[tuple[float, float, float], ...]) -> str:
    return ", ".join(f"{value:.12g}" for point in points for value in point)


def run(binary: Path, workdir: Path, input_text: str, output_name: str) -> subprocess.CompletedProcess[str]:
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "input.nml"
    input_path.write_text(input_text, encoding="utf-8")
    completed = subprocess.run(
        [str(binary.resolve()), input_path.name],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=7200,
        check=False,
    )
    (workdir / f"{output_name}.log").write_text(completed.stdout, encoding="utf-8")
    return completed


def validate_positive(workdir: Path, points: tuple[tuple[float, float, float], ...], expect_covariance: bool) -> dict[str, object]:
    output = workdir / "compact.dat"
    if not output.exists():
        raise RuntimeError(f"compact-Dyson output is missing: {output}")
    lines = output.read_text(encoding="utf-8", errors="replace").splitlines()
    meta = headers(lines)
    if lines[0] != "# TDDFT compact Dyson response":
        raise RuntimeError("generic compact-Dyson provenance header is missing")
    if any("TDVK" in line or "PASS CANDIDATE" in line for line in lines):
        raise RuntimeError("campaign-specific production-output label leaked into compact-Dyson output")
    if meta.get("backend") != "compact_dyson" or meta.get("accepted_state_cache_reused") != "T":
        raise RuntimeError("accepted compact-Dyson handoff metadata is incomplete")
    if meta.get("accepted_k_mesh") != "8 8 8" or meta.get("accepted_k_count") != "512":
        raise RuntimeError("the non-12^3 accepted-state contract was not exercised")
    if meta.get("output_storage_mode") != "streaming_summary" or meta.get("maximum_retained_q_batches") != "1 (q-local bare, Dyson, and loss matrices)":
        raise RuntimeError("bounded streaming storage provenance is missing")
    if meta.get("dyson_static_audit") != "F" or meta.get("gf_closure_audit") != "F":
        raise RuntimeError("optional compact audits were not disabled")
    if expect_covariance:
        if meta.get("validate_interacting_covariance") != "T" or not any("# interacting_covariance =" in line for line in lines):
            raise RuntimeError("requested covariance audit is missing")
    else:
        if meta.get("validate_interacting_covariance") != "F" or any("# interacting_covariance =" in line for line in lines):
            raise RuntimeError("covariance was performed or serialized when disabled")
    if any("# static_denominator =" in line or "# gf_spots columns:" in line for line in lines):
        raise RuntimeError("disabled static or GF audit was serialized")

    data = [line.split() for line in lines if line and not line.startswith("#")]
    dynamic = [row for row in data if len(row) == 15]
    covariance = [row for row in data if len(row) == 7]
    expected_omega = (0.0, 0.01, 0.015)
    if len(dynamic) != len(points) * len(expected_omega):
        raise RuntimeError(f"unexpected dynamic row count: {len(dynamic)}")
    for row in dynamic:
        q_index = int(row[0])
        q = tuple(number(value) for value in row[1:4])
        if q_index < 1 or q_index > len(points) or max(abs(a - b) for a, b in zip(q, points[q_index - 1])) > 1.0e-12:
            raise RuntimeError(f"unexpected q in compact-Dyson row: {row}")
        if min(abs(number(row[4]) - target) for target in expected_omega) > 1.0e-12:
            raise RuntimeError(f"unexpected omega in compact-Dyson row: {row}")
        values = [number(value) for value in row[5:]]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("non-finite compact-Dyson diagnostic")
        if values[-3] > 1.0e-8 or values[-1] > 1.0e-8:
            raise RuntimeError("compact-Dyson solve residual exceeds the production threshold")
    if expect_covariance and len(covariance) != len(expected_omega):
        raise RuntimeError(f"unexpected covariance row count: {len(covariance)}")
    if not expect_covariance and covariance:
        raise RuntimeError("unexpected covariance data in production-only output")
    return {
        "q_count": len(points),
        "omega_count": len(expected_omega),
        "dynamic_rows": len(dynamic),
        "covariance_rows": len(covariance),
        "mesh": [8, 8, 8],
        "product_dimension": int(meta["product_dimension"]),
        "maximum_retained_q_batches": 1,
        "output": str(output),
    }


def validate_multisite(workdir: Path) -> dict[str, object]:
    output = workdir / "compact_multisite_fe.dat"
    if not output.exists():
        raise RuntimeError(f"multi-site compact-Dyson output is missing: {output}")
    lines = output.read_text(encoding="utf-8", errors="replace").splitlines()
    meta = headers(lines)
    if lines[0] != "# TDDFT compact Dyson response" or meta.get("accepted_k_mesh") != "4 4 4":
        raise RuntimeError("multi-site compact-Dyson provenance is incomplete")
    if meta.get("accepted_k_count") != "64" or int(meta.get("product_dimension", "0")) == 232:
        raise RuntimeError("multi-site runtime product dimension did not differ from 232")
    data = [line.split() for line in lines if line and not line.startswith("#")]
    dynamic = [row for row in data if len(row) == 15]
    if len(dynamic) != 1 or not all(math.isfinite(number(value)) for value in dynamic[0][5:]):
        raise RuntimeError("multi-site compact-Dyson response is incomplete or non-finite")
    if meta.get("dyson_static_audit") != "F" or meta.get("validate_interacting_covariance") != "F" or meta.get("gf_closure_audit") != "F":
        raise RuntimeError("multi-site production case unexpectedly enabled an audit")
    return {
        "mesh": [4, 4, 4],
        "q_count": 1,
        "omega_count": 1,
        "dynamic_rows": 1,
        "product_dimension": int(meta["product_dimension"]),
        "output": str(output),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    template = root / "tests/integration/tddft_driver_smoke/input_tdvk07_fe.nml"
    database = root / "tests/scf/cases/bulk/bccFe"
    template_text = template.read_text(encoding="utf-8")
    results: dict[str, object] = {}

    for name, case in CASES.items():
        workdir = args.scratch_root / name
        text = template_text
        text = replace_assignment(text, "database", repr(str(database)))
        text = replace_assignment(text, "nk1", "8")
        text = replace_assignment(text, "nk2", "8")
        text = replace_assignment(text, "nk3", "8")
        text = replace_assignment(text, "n_q", str(len(case["q"])))
        text = replace_assignment(text, "q_list", q_literal(case["q"]))
        text = replace_assignment(text, "n_omega", "3")
        text = replace_assignment(text, "omega_grid", "0.0, 0.01, 0.015")
        text = replace_assignment(text, "omega_min", "0.0")
        text = replace_assignment(text, "omega_max", "0.015")
        text = replace_assignment(text, "gf_closure_audit", ".false.")
        text = replace_assignment(text, "dyson_static_audit", case["static"])
        text = replace_assignment(text, "validate_interacting_covariance", case["covariance"])
        text = replace_assignment(text, "write_full_matrix", ".false.")
        text = replace_assignment(text, "output_file", repr("compact.dat"))
        completed = run(args.binary, workdir, text, "compact")
        if completed.returncode != 0:
            raise RuntimeError(f"{name} executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
        results[name] = validate_positive(workdir, case["q"], case["covariance"] == ".true.")

    for name, flag, q_values, expected in (
        ("static_audit_without_gamma", "dyson_static_audit", "0.03, 0.0, 0.0", "dyson_static_audit requires Gamma"),
        ("covariance_without_pair", "validate_interacting_covariance", "0.03, 0.0, 0.0", "requires an exact +q/-q pair"),
    ):
        workdir = args.scratch_root / name
        text = template_text
        text = replace_assignment(text, "database", repr(str(database)))
        text = replace_assignment(text, "n_q", "1")
        text = replace_assignment(text, "q_list", q_values)
        text = replace_assignment(text, "dyson_static_audit", ".false.")
        text = replace_assignment(text, "validate_interacting_covariance", ".false.")
        text = replace_assignment(text, flag, ".true.")
        text = replace_assignment(text, "output_file", repr("must_not_exist.dat"))
        completed = run(args.binary, workdir, text, "preflight")
        if completed.returncode == 0 or expected not in completed.stdout:
            raise RuntimeError(f"{name} did not fail in input preflight as expected:\n{completed.stdout[-12000:]}")
        if (workdir / "must_not_exist.dat").exists() or (workdir / "must_not_exist.dat.state").exists():
            raise RuntimeError(f"{name} performed response work before rejecting the input")
        results[name] = {"returncode": completed.returncode, "preflight_error": expected}

    multisite_workdir = args.scratch_root / "non232_multisite"
    multisite_source = root / "tests/regression/wp9_validation/commensurate_supercell/gbt_supercell/q050/super_scf"
    multisite_text = (root / "tests/integration/tddft_driver_smoke/input_tddft_multisite_fe.nml").read_text(encoding="utf-8")
    multisite_text = replace_assignment(multisite_text, "database", repr(str(multisite_source)))
    multisite_workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(multisite_source / "lattice.nml", multisite_workdir / "lattice.nml")
    completed = run(args.binary, multisite_workdir, multisite_text, "compact_multisite_fe")
    if completed.returncode != 0:
        raise RuntimeError(f"non232_multisite executable failed with {completed.returncode}:\n{completed.stdout[-12000:]}")
    results["non232_multisite"] = validate_multisite(multisite_workdir)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
