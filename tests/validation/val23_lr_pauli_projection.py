#!/usr/bin/env python3
"""Live LR-02N numerical closure for converged scalar-relativistic bcc Fe.

This is a consistency oracle, not a physics pass/fail threshold.  It checks
that the emitted Pauli density is built from the accepted LR-01 snapshot,
that the declared pointwise and integrated quantities agree, and that the
occupied coefficient sum is separately recorded.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
from pathlib import Path

import f90nml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "scf" / "cases" / "bulk" / "bccFe"
RUNNER = ROOT / "tests" / "run_binary.sh"
SR_FILE = "radial_ground_state_Fe_1.dat"
P_FILE = "lr02n_pauli_projection_Fe_1.dat"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def close_to(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise RuntimeError(f"{label}: {actual:.16e} != {expected:.16e} (tol {tolerance:.3e})")


def patch_case(destination: Path, nstep: int, nk: int) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(FIXTURE, destination)
    input_path = destination / "input.nml"
    patched = destination / "input.patched.nml"
    f90nml.patch(
        str(input_path),
        {
            "calculation": {"pre_processing": "bravais", "post_processing": "pauli_projection"},
            "self": {"nstep": nstep, "conv_thr": 2.0e-6, "use_kspace": False},
            "control": {"nsp": 1},
            "reciprocal": {
                "nk1": nk,
                "nk2": nk,
                "nk3": nk,
                "use_symmetry_reduction": False,
                "use_time_reversal": True,
                "temperature": 300.0,
                "auto_find_fermi": False,
                "reciprocal_mode": "ham_only",
                "kspace_ham_order": "first",
            },
        },
        str(patched),
    )
    patched.replace(input_path)


def parse_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if " = " in line and not line.startswith("#"):
            key, value = line.split(" = ", 1)
            values[key.strip()] = value.strip()
    return values


def parse_radial_file(path: Path, expected_columns: int) -> tuple[dict[str, str], list[list[float]]]:
    metadata: dict[str, str] = {}
    rows: list[list[float]] = []
    data_started = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# columns:"):
            data_started = True
            continue
        if not data_started:
            if line.startswith("# ") and " = " in line:
                key, value = line[2:].split(" = ", 1)
                metadata[key] = value.strip()
            continue
        if line.strip():
            fields = line.split()
            require(len(fields) == expected_columns, f"{path.name}: expected {expected_columns} columns, got {len(fields)}")
            rows.append([number(value) for value in fields])
    require(len(rows) >= 3, f"{path.name}: too few radial rows")
    require(all(math.isfinite(value) for row in rows for value in row), f"{path.name}: non-finite radial value")
    return metadata, rows


def simpson_integral(values: list[float], radii: list[float], a: float, b: float) -> float:
    total = 0.0
    nr = len(values)
    for index, value in enumerate(values, start=1):
        weight = 2.0 * ((index + 1) % 2 + 1) / 3.0
        if index == 1 or index == nr:
            weight = 1.0 / 3.0
        total += weight * a * (radii[index - 1] + b) * value
    return total


def volume_integral(values: list[float], radii: list[float], a: float, b: float) -> float:
    return simpson_integral(
        [4.0 * math.pi * radius * radius * value for radius, value in zip(radii, values)], radii, a, b
    )


def run_case(binary: Path, destination: Path, nstep: int, nk: int) -> tuple[dict[str, str], dict[str, str], list[list[float]], str]:
    patch_case(destination, nstep, nk)
    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=destination,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=1800,
        check=False,
    )
    log_path = destination / "testrun.log"
    log = log_path.read_text(encoding="utf-8") if log_path.exists() else result.stdout
    (destination / "run.log").write_text(log, encoding="utf-8")
    require(result.returncode == 0 and "fatal" not in log.lower(), f"LR-02N SCF failed:\n{log[-8000:]}")
    require("Converged!" in log, "LR-02N case did not reach outer SCF convergence")

    sr_metadata, _ = parse_radial_file(destination / SR_FILE, 12)
    p_metadata, p_rows = parse_radial_file(destination / P_FILE, 13)
    summary = parse_key_value_file(destination / (P_FILE + ".summary"))
    return sr_metadata, {**p_metadata, **summary}, p_rows, log


def check_case(destination: Path, sr_metadata: dict[str, str], metadata: dict[str, str], rows: list[list[float]]) -> None:
    require(metadata.get("reference") == "exact accepted LR-01 radial ground-state snapshot", "wrong SR reference")
    require("large component only" in metadata.get("pauli_state", ""), "Pauli route is not declared large-component-only")
    require("no GFAC" in metadata.get("pauli_state", ""), "Pauli route advertises GFAC contamination")
    require("lower component" in metadata.get("pauli_state", ""), "Pauli route does not document lower-component exclusion")
    require(metadata.get("reciprocal_mode") == "ham_only", "unexpected reciprocal eigenproblem mode")
    require(metadata.get("occupation_rule") == "Fermi-Dirac at reciprocal temperature", "occupation provenance missing")

    a, b = (number(value) for value in sr_metadata["mesh_a mesh_b"].split())
    radii = [row[0] for row in rows]
    for index, radius in enumerate(radii):
        expected = 0.0 if index == 0 else b * (math.exp(a * index) - 1.0)
        close_to(radius, expected, 2.0e-12, f"logarithmic mesh row {index + 1}")

    sr_up = [row[1] for row in rows]
    sr_down = [row[2] for row in rows]
    p_up = [row[3] for row in rows]
    p_down = [row[4] for row in rows]
    d_up = [row[5] for row in rows]
    d_down = [row[6] for row in rows]
    m_sr = [row[7] for row in rows]
    m_p = [row[8] for row in rows]
    d_m = [row[9] for row in rows]
    charge_sr = [row[10] for row in rows]
    charge_p = [row[11] for row in rows]
    d_charge = [row[12] for row in rows]

    for index in range(len(rows)):
        close_to(d_up[index], sr_up[index] - p_up[index], 3.0e-10, f"delta n_up row {index + 1}")
        close_to(d_down[index], sr_down[index] - p_down[index], 3.0e-10, f"delta n_down row {index + 1}")
        close_to(m_sr[index], sr_up[index] - sr_down[index], 3.0e-10, f"m_SR row {index + 1}")
        close_to(m_p[index], p_up[index] - p_down[index], 3.0e-10, f"m_P row {index + 1}")
        close_to(d_m[index], m_sr[index] - m_p[index], 3.0e-10, f"delta m row {index + 1}")
        close_to(charge_sr[index], sr_up[index] + sr_down[index], 3.0e-10, f"charge_SR row {index + 1}")
        close_to(charge_p[index], p_up[index] + p_down[index], 3.0e-10, f"charge_P row {index + 1}")
        close_to(d_charge[index], charge_sr[index] - charge_p[index], 3.0e-10, f"delta charge row {index + 1}")

    sr_integrals = [volume_integral(sr_up, radii, a, b), volume_integral(sr_down, radii, a, b)]
    p_integrals = [volume_integral(p_up, radii, a, b), volume_integral(p_down, radii, a, b)]
    delta_integrals = [sr_integrals[0] - p_integrals[0], sr_integrals[1] - p_integrals[1]]
    for label, actual, key in (
        ("N_up_SR", sr_integrals[0], "N_up_SR"),
        ("N_down_SR", sr_integrals[1], "N_down_SR"),
        ("N_up_P", p_integrals[0], "N_up_P"),
        ("N_down_P", p_integrals[1], "N_down_P"),
        ("Delta_N_up", delta_integrals[0], "Delta_N_up"),
        ("Delta_N_down", delta_integrals[1], "Delta_N_down"),
    ):
        close_to(actual, number(metadata[key]), 3.0e-8, label)

    close_to(sr_integrals[0], number(sr_metadata["integrated_n_up"]), 3.0e-8, "SR up reproduces LR-01 charge")
    close_to(sr_integrals[1], number(sr_metadata["integrated_n_down"]), 3.0e-8, "SR down reproduces LR-01 charge")
    close_to(delta_integrals[0] - delta_integrals[1], number(metadata["Delta_M"]), 3.0e-8, "Delta_M")

    norm_m_sr = math.sqrt(volume_integral([value * value for value in m_sr], radii, a, b))
    norm_dm = math.sqrt(volume_integral([value * value for value in d_m], radii, a, b))
    norm_q_sr = math.sqrt(volume_integral([value * value for value in charge_sr], radii, a, b))
    norm_dq = math.sqrt(volume_integral([value * value for value in d_charge], radii, a, b))
    close_to(norm_dm / norm_m_sr, number(metadata["relative_L2_magnetization"]), 3.0e-8, "relative magnetization L2")
    close_to(norm_dq / norm_q_sr, number(metadata["relative_L2_total_charge"]), 3.0e-8, "relative charge L2")

    direct_total = number(metadata["direct_coefficient_weight_up"]) + number(metadata["direct_coefficient_weight_down"])
    close_to(direct_total, number(metadata["occupied_valence_weight"]), 3.0e-8, "occupied coefficient sum")
    require(int(metadata["magnetic_region_last_point"]) >= 2, "objective magnetic region is empty")
    require(number(metadata["magnetic_region_fraction"]) == 0.9, "magnetic region fraction is not declared")

    for index in (1, len(rows) // 2, len(rows) - 1):
        print(
            f"LR-02N sample row {index + 1}: r={radii[index]:.8e} "
            f"delta_n(up/down)={d_up[index]:.8e}/{d_down[index]:.8e} delta_m={d_m[index]:.8e}"
        )

    print(
        f"LR-02N Fe N_SR(up/down)={sr_integrals[0]:.10f}/{sr_integrals[1]:.10f}; "
        f"N_P(up/down)={p_integrals[0]:.10f}/{p_integrals[1]:.10f}"
    )
    print(
        f"LR-02N Fe Delta_N(up/down)={delta_integrals[0]:.10e}/{delta_integrals[1]:.10e}; "
        f"Delta_M={number(metadata['Delta_M']):.10e}"
    )
    print(
        f"LR-02N Fe relative L2(magnetization/charge)="
        f"{number(metadata['relative_L2_magnetization']):.10e}/"
        f"{number(metadata['relative_L2_total_charge']):.10e}"
    )
    print(
        f"LR-02N Fe magnetic region: {metadata['magnetic_region_fraction']} of |m_SR|, "
        f"rmax={number(metadata['magnetic_region_rmax_bohr']):.8e} bohr"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--nstep", type=int, default=16)
    parser.add_argument("--nk", type=int, default=8)
    parser.add_argument("--repeat", type=int, default=1, help="repeat in separate directories for determinism checks")
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    fingerprints: list[tuple[str, str]] = []
    for repeat in range(args.repeat):
        destination = args.scratch_root / f"bccFe_{repeat + 1}"
        sr_metadata, metadata, rows, _ = run_case(args.binary.resolve(), destination, args.nstep, args.nk)
        check_case(destination, sr_metadata, metadata, rows)
        fingerprints.append(
            ((destination / P_FILE).read_text(encoding="utf-8"), (destination / (P_FILE + ".summary")).read_text(encoding="utf-8"))
        )

    if len(fingerprints) > 1:
        require(fingerprints[0] == fingerprints[1], "repeated Fe LR-02N execution is not deterministic")
        print("LR-02N repeated Fe execution: byte-identical radial data and summary")
    print("Val23LrPauliProjection: PASS (reproducible numerical closure; no physical discrepancy threshold imposed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
