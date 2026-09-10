#!/usr/bin/env python3
"""LR-01 live radial ground-state oracle for the magnetic bcc-Fe fixture.

The executable emits the snapshot only after ATOMSC has accepted its final
inner iteration.  This test independently integrates the emitted weighted
density, checks the exact XC/Pauli decomposition, and compares the radial
spin number with the production magnetic-moment report.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "scf" / "cases" / "bulk" / "bccFe"
RUNNER = ROOT / "tests" / "run_binary.sh"
SNAPSHOT = "radial_ground_state_Fe_1.dat"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def patch_case(destination: Path, nstep: int, txc: int) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(FIXTURE, destination)
    input_path = destination / "input.nml"
    patched = destination / "input.patched.nml"
    f90nml.patch(
        str(input_path),
        {
            "calculation": {"pre_processing": "bravais", "post_processing": "none"},
            "self": {"nstep": nstep, "conv_thr": 1.0e-6},
            "control": {"nsp": 2, "txc": txc},
            "mix": {"beta": 0.01, "mixtype": "broyden"},
        },
        str(patched),
    )
    patched.replace(input_path)


def run_case(binary: Path, destination: Path, nstep: int, txc: int) -> tuple[dict[str, str], list[list[str]], str]:
    patch_case(destination, nstep, txc)
    result = subprocess.run(
        ["/bin/bash", str(RUNNER), str(binary)],
        cwd=destination,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=900,
        check=False,
    )
    log_path = destination / "testrun.log"
    log = log_path.read_text(encoding="utf-8") if log_path.exists() else result.stdout
    (destination / "run.log").write_text(log, encoding="utf-8")
    if result.returncode != 0 or "fatal" in log.lower():
        raise RuntimeError(f"bcc-Fe LR-01 SCF failed:\n{log[-6000:]}")

    snapshot = destination / SNAPSHOT
    if not snapshot.exists():
        raise RuntimeError(f"SCF did not emit {SNAPSHOT}")
    metadata: dict[str, str] = {}
    rows: list[list[str]] = []
    data_started = False
    for line in snapshot.read_text(encoding="utf-8").splitlines():
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
            if len(fields) != 12:
                raise RuntimeError(f"malformed snapshot row ({len(fields)} fields): {line}")
            rows.append(fields)
    if len(rows) < 3:
        raise RuntimeError("snapshot has fewer than three radial points")
    return metadata, rows, log


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def close_to(actual: float, expected: float, tolerance: float, label: str) -> None:
    if abs(actual - expected) > tolerance:
        raise RuntimeError(f"{label}: {actual:.16e} != {expected:.16e} (tol {tolerance:.3e})")


def simpson_integral(values: list[float], radii: list[float], a: float, b: float) -> float:
    total = 0.0
    nr = len(values)
    for index, value in enumerate(values, start=1):
        weight = 2.0 * ((index + 1) % 2 + 1) / 3.0
        if index == 1 or index == nr:
            weight = 1.0 / 3.0
        total += weight * a * (radii[index - 1] + b) * value
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--nstep", type=int, default=8)
    parser.add_argument("--txc", type=int, default=1,
                        help="XC selector; use 101 to run the same oracle with libXC enabled")
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    metadata, fields, log = run_case(args.binary.resolve(), args.scratch_root / "bccFe", args.nstep, args.txc)
    require(metadata.get("valid") == "true", "snapshot is not marked valid")
    require(metadata.get("accepted") == "true", "snapshot is not marked accepted")
    require(int(metadata["accepted_inner_iteration"]) > 0, "accepted inner iteration was not retained")
    require(number(metadata["accepted_residual_control"]) <= 1.0e-6, "accepted radial residual exceeds tolerance")
    require(
        metadata.get("channel_convention") == "up=RHO(:,1)=Vxc(:,1); down=RHO(:,2)=Vxc(:,2)",
        "unexpected radial channel convention",
    )
    expected_libxc = args.txc >= 100
    require(metadata.get("xc_use_libxc") == ("T" if expected_libxc else "F"),
            "snapshot XC backend provenance does not match the selected route")
    require(metadata.get("xc_backend") == ("libXC" if expected_libxc else "legacy RS-LMTO"),
            "snapshot XC backend name is wrong")
    if expected_libxc:
        require(int(metadata.get("xc_component_count", "0")) > 0, "libXC component provenance is missing")
    else:
        require("Barth-Hedin" in metadata.get("xc_functional", ""), "XC functional provenance is missing")

    a, b = (number(value) for value in metadata["mesh_a mesh_b"].split())
    columns = [[number(value) for value in row[1:]] for row in fields]
    radii = [row[0] for row in columns]
    rho_up = [row[1] for row in columns]
    rho_down = [row[2] for row in columns]
    n_up = [row[3] for row in columns]
    n_down = [row[4] for row in columns]
    vxc_up = [row[5] for row in columns]
    vxc_down = [row[6] for row in columns]
    delta = [row[7] for row in columns]
    scalar = [row[8] for row in columns]
    bxc = [row[9] for row in columns]
    total_split = [row[10] for row in columns]
    b_fsm = number(metadata["constraining_field_ry"])

    for index, radius in enumerate(radii):
        expected_radius = 0.0 if index == 0 else b * (math.exp(a * index) - 1.0)
        close_to(radius, expected_radius, 2.0e-12, f"logarithmic mesh row {index + 1}")
    require(all(math.isfinite(value) for row in columns for value in row), "snapshot contains a non-finite value")

    independent_up = simpson_integral(rho_up, radii, a, b)
    independent_down = simpson_integral(rho_down, radii, a, b)
    expected_origin_up = (
        rho_up[1] / radii[1] ** 2 * radii[2] - rho_up[2] / radii[2] ** 2 * radii[1]
    ) / (4.0 * math.pi * (radii[2] - radii[1]))
    expected_origin_down = (
        rho_down[1] / radii[1] ** 2 * radii[2] - rho_down[2] / radii[2] ** 2 * radii[1]
    ) / (4.0 * math.pi * (radii[2] - radii[1]))
    close_to(n_up[0], expected_origin_up, 1.0e-8, "up origin extrapolation")
    close_to(n_down[0], expected_origin_down, 1.0e-8, "down origin extrapolation")
    physical_up = [4.0 * math.pi * radius * radius * density for radius, density in zip(radii, n_up)]
    physical_down = [4.0 * math.pi * radius * radius * density for radius, density in zip(radii, n_down)]
    close_to(independent_up, simpson_integral(physical_up, radii, a, b), 2.0e-10, "up charge quadrature")
    close_to(independent_down, simpson_integral(physical_down, radii, a, b), 2.0e-10, "down charge quadrature")
    close_to(independent_up, number(metadata["integrated_n_up"]), 2.0e-10, "stored up charge")
    close_to(independent_down, number(metadata["integrated_n_down"]), 2.0e-10, "stored down charge")
    close_to(independent_up + independent_down, 26.0, 5.0e-4, "independent Fe charge")
    require(independent_up - independent_down > 0.5, "magnetic oracle is not spin polarized")
    close_to(independent_up - independent_down, number(metadata["integrated_spin_number"]), 2.0e-10,
             "stored spin number")

    for index in range(len(fields)):
        close_to(delta[index], vxc_up[index] - vxc_down[index], 2.0e-11, f"delta Vxc row {index + 1}")
        close_to(scalar[index], 0.5 * (vxc_up[index] + vxc_down[index]), 2.0e-11,
                 f"scalar Vxc row {index + 1}")
        close_to(bxc[index], 0.5 * delta[index], 2.0e-11, f"Pauli Bxc row {index + 1}")
        close_to(total_split[index], delta[index] + 2.0 * b_fsm, 2.0e-11,
                 f"total spin splitting row {index + 1}")
    require(max(abs(value) for value in delta) > 1.0e-5, "live Vxc spin splitting is identically zero")

    report = (args.scratch_root / "bccFe" / "report.out").read_text(encoding="utf-8")
    match = re.search(r"Spin moment projections of atom\s+1:\s*([+-]?\S+)\s+([+-]?\S+)\s+([+-]?\S+)", report)
    require(match is not None, "native magnetic-moment report is missing")
    reported = [number(value) for value in match.groups()]
    reported_metadata = [number(value) for value in metadata["reported_moment_xyz_muB"].split()]
    for index in range(3):
        close_to(reported[index], reported_metadata[index], 2.0e-6, f"reported moment component {index}")
    close_to(reported[2], independent_up - independent_down, 2.0e-4, "reported z moment")
    close_to(reported[2], number(metadata["integrated_moment_muB"]), 2.0e-4, "reported/integrated moment")

    print(f"LR-01 bcc-Fe accepted inner iteration: {metadata['accepted_inner_iteration']}")
    print(
        f"LR-01 independent N(up/down/total) = {independent_up:.10f} "
        f"{independent_down:.10f} {independent_up + independent_down:.10f}"
    )
    print(
        f"LR-01 independent spin number / reported z moment = "
        f"{independent_up - independent_down:.10f} / {reported[2]:.10f} mu_B"
    )
    print(f"LR-01 max |Vxc(up)-Vxc(down)| = {max(abs(value) for value in delta):.6e} Ry")
    print(f"LR-01 provenance = {metadata['xc_backend']}; {metadata['xc_functional']}; txc={metadata['xc_txc']}")
    for index in (0, 1, len(fields) // 2, len(fields) - 1):
        print(
            f"LR-01 sample row {index + 1}: r={radii[index]:.8e} "
            f"n(up/down)={n_up[index]:.8e}/{n_down[index]:.8e} "
            f"deltaVxc={delta[index]:.8e} Ry"
        )
    require("Converged!" in log, "real oracle did not reach outer SCF convergence")
    print("Val22LrRadialGroundState: PASS (live converged magnetic oracle)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
