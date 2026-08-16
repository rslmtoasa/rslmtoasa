#!/usr/bin/env python3
"""VAL-13: deterministic end-to-end LMTO -> abspinlib spin dynamics.

The oracle is deliberately small and production-facing: after the electronic
refresh, a one-site bcc-Fe scalar-relativistic collinear state has a
longitudinal LMTO field and zero transverse torque.  At zero temperature and
zero damping, the Depondt update must preserve the unit direction and produce
no spurious transverse motion.  The same short trajectory is repeated at
smaller timesteps, and the final electronic energy is compared with a zero-step
run as a bounded feedback diagnostic rather than as an independent energy
functional oracle.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FLOAT = r"[-+0-9.EeDd]+"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def run_case(binary: Path, root: Path, base: Path, workdir: Path,
             *, dt: float, steps: int) -> dict[str, object]:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)

    input_path = workdir / "input.nml"
    patched_input = workdir / "input.patched.nml"
    f90nml.patch(
        str(input_path),
        {
            "calculation": {"pre_processing": "bravais", "processing": "sd"},
            "control": {"nsp": 2, "recur": "block", "lld": 20},
            "self": {"nstep": 1, "soc_scale": 0.0},
            "hamiltonian": {"hoh": False},
            "sd": {"asd_step": steps, "dt": dt, "alpha": 0.0, "sd_temp": 0.0},
        },
        str(patched_input),
    )
    patched_input.replace(input_path)

    atom_path = workdir / "Fe.nml"
    patched_atom = workdir / "Fe.patched.nml"
    f90nml.patch(str(atom_path), {"par": {"mom": [0.0, 0.0, 1.0]}}, str(patched_atom))
    patched_atom.replace(atom_path)

    result = subprocess.run(
        ["/bin/bash", str(root / "tests/run_binary.sh"), str(binary)],
        cwd=workdir, capture_output=True, text=True,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"VAL-13 run failed (dt={dt}, steps={steps}):\n{log[-6000:]}")

    trajectory: list[list[float]] = []
    trajectory_path = workdir / "output.lammpstrj"
    lines = trajectory_path.read_text().splitlines() if trajectory_path.exists() else []
    for index, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            fields = lines[index + 1].split()
            trajectory.append([number(value) for value in fields[4:7]])
    if steps > 0 and not trajectory:
        raise RuntimeError("VAL-13: production trajectory is empty")

    field_matches = re.findall(
        r"Magnetic field on atom\s+1 is\s+(%s)\s+(%s)\s+(%s)" % (FLOAT, FLOAT, FLOAT), log
    )
    energy_match = re.findall(r"Total energy of system:\s*(%s)" % FLOAT,
                              (workdir / "report.out").read_text(errors="replace"))
    if not field_matches or not energy_match:
        raise RuntimeError("VAL-13: missing production field or energy output")

    fields = [[number(value) for value in match] for match in field_matches]
    return {
        "dt_ps": dt,
        "steps": steps,
        "trajectory": trajectory,
        "field_trace_T": fields,
        "final_field_T": fields[-1],
        "energy_Ry": number(energy_match[-1]),
    }


def check_trajectory(case: dict[str, object], label: str) -> dict[str, float]:
    trajectory = case["trajectory"]
    assert isinstance(trajectory, list)
    norms = [norm(vector) for vector in trajectory]
    norm_error = max(abs(value - 1.0) for value in norms)
    continuity = min(dot(left, right) for left, right in zip(trajectory, trajectory[1:])) if len(trajectory) > 1 else 1.0
    final = trajectory[-1]
    require(norm_error < 2.0e-11, f"{label}: spin norm error {norm_error:.3e}")
    require(continuity > 0.999999999, f"{label}: trajectory discontinuity ({continuity:.12f})")
    return {"norm_error": norm_error, "minimum_step_dot": continuity,
            "final_x": final[0], "final_y": final[1], "final_z": final[2]}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("VAL-13 FAIL: " + message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    base = root / "tests/scf/cases/bulk/bccFe"
    scratch = args.scratch_root.resolve()
    scratch.mkdir(parents=True, exist_ok=True)

    # Fixed short interval: 0.01 ps.  The finest run is also the primary
    # trajectory; the other two runs establish short-interval timestep trends.
    cases = [
        run_case(args.binary.resolve(), root, base, scratch / "dt_0p01", dt=0.01, steps=1),
        run_case(args.binary.resolve(), root, base, scratch / "dt_0p005", dt=0.005, steps=2),
        run_case(args.binary.resolve(), root, base, scratch / "dt_0p0025", dt=0.0025, steps=4),
        run_case(args.binary.resolve(), root, base, scratch / "zero_step", dt=0.01, steps=0),
    ]
    coarse, reference, fine, zero_step = cases

    evidence = {
        "coarse": check_trajectory(coarse, "dt=0.01 ps"),
        "reference": check_trajectory(reference, "dt=0.005 ps"),
        "fine": check_trajectory(fine, "dt=0.0025 ps"),
    }

    refreshed_field = reference["final_field_T"]
    reference_final = reference["trajectory"][-1]
    require(abs(refreshed_field[2]) > 1.0e3 and norm(refreshed_field[:2]) < 1.0e-2,
            f"refreshed LMTO equilibrium field is not longitudinal: B={refreshed_field}")
    require(reference_final[2] > 0.999999999 and norm(reference_final[:2]) < 1.0e-8,
            f"zero-torque direction was not preserved: m={reference_final}")

    coarse_final = coarse["trajectory"][-1]
    fine_final = fine["trajectory"][-1]
    coarse_fine_distance = norm([a - b for a, b in zip(coarse_final, fine_final)])
    reference_fine_distance = norm([a - b for a, b in zip(reference_final, fine_final)])
    timestep_envelope = max(coarse_fine_distance, reference_fine_distance)
    require(timestep_envelope < 1.0e-8,
            f"timestep trajectories disagree beyond the zero-torque envelope: {timestep_envelope:.3e}")

    energy_drift = abs(reference["energy_Ry"] - zero_step["energy_Ry"])
    require(energy_drift < 2.0e-4,
            f"zero-damping electronic energy feedback drift is too large: {energy_drift:.3e} Ry")

    evidence.update({
        "refreshed_field_T": refreshed_field,
        "coarse_fine_distance": coarse_fine_distance,
        "reference_fine_distance": reference_fine_distance,
        "timestep_envelope": timestep_envelope,
        "zero_damping_energy_drift_Ry": energy_drift,
        "oracle": "LMTO longitudinal-field/zero-torque limit, Depondt norm/continuity, and short-interval timestep refinement",
    })
    print(json.dumps(evidence, indent=2, sort_keys=True))
    print("VAL-13 PASS: deterministic ab-initio spin-dynamics loop validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
