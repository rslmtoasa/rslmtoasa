#!/usr/bin/env python3
"""VAL-12: validate the production LMTO field/torque feedback seam.

The campaign deliberately does not reproduce the LMTO field expression.  It
uses the production ``report.out`` and ``testrun.log`` records, checking only
invariants and the published field/torque relation.  The constrained-field
finite-difference oracle remains the VAL-10 ``UnitConstrainingField`` test.
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


def patch_case(base: Path, workdir: Path, *, nsp: int, soc_scale: float,
               moment: list[float]) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)

    input_path = workdir / "input.nml"
    input_patch: dict[str, dict[str, object]] = {
        "calculation": {"pre_processing": "bravais", "processing": "none"},
        "control": {"nsp": nsp, "recur": "block", "lld": 20},
        "self": {"nstep": 1, "soc_scale": soc_scale},
        "hamiltonian": {"hoh": False},
    }
    patched_input = workdir / "input.patched.nml"
    f90nml.patch(str(input_path), input_patch, str(patched_input))
    patched_input.replace(input_path)

    atom_path = workdir / "Fe.nml"
    patched_atom = workdir / "Fe.patched.nml"
    f90nml.patch(str(atom_path), {"par": {"mom": moment}}, str(patched_atom))
    patched_atom.replace(atom_path)


def run_case(binary: Path, root: Path, base: Path, workdir: Path, **kwargs: object) -> dict:
    patch_case(base, workdir, **kwargs)
    result = subprocess.run(
        ["/bin/bash", str(root / "tests/run_binary.sh"), str(binary)],
        cwd=workdir, capture_output=True, text=True,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"VAL-12 run failed ({kwargs}):\n{log[-6000:]}")
    return {
        "parameters": kwargs,
        "report": parse_report(workdir / "report.out", workdir / "Fe_out.nml"),
        "log": parse_log(log),
    }


def parse_report(path: Path, atom_path: Path) -> dict[str, list[float] | float]:
    text = path.read_text(errors="replace")
    energy_match = re.search(r"Total energy of system:\s*(%s)" % FLOAT, text)
    moment_match = re.search(
        r"Spin moment projections of atom\s*1\s*:\s*(%s)\s+(%s)\s+(%s)" % (FLOAT, FLOAT, FLOAT), text
    )
    field_match = re.search(
        r"Magnetic force on atom\s*1\s*:\s*(%s)\s+(%s)\s+(%s)" % (FLOAT, FLOAT, FLOAT), text
    )
    if not (energy_match and moment_match and field_match):
        raise RuntimeError(f"VAL-12: incomplete production report {path}")
    atom = f90nml.read(str(atom_path))
    direction = [float(value) for value in atom["par"]["mom"]]
    return {
        "energy": number(energy_match.group(1)),
        "moment": [number(value) for value in moment_match.groups()],
        "direction": direction,
        "field": [number(value) for value in field_match.groups()],
    }


def parse_log(text: str) -> dict[str, list[float]]:
    torque_matches = re.findall(
        r"Magnetic torque on atom\s*1 is\s*(%s)\s+(%s)\s+(%s)" % (FLOAT, FLOAT, FLOAT), text
    )
    if not torque_matches:
        raise RuntimeError("VAL-12: production torque diagnostic was not emitted")
    return {
        "torque": [number(value) for value in torque_matches[-1]],
    }


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def cross(left: list[float], right: list[float]) -> list[float]:
    return [left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("VAL-12 FAIL: " + message)


def check_production_route(case: dict, label: str) -> dict[str, float]:
    report = case["report"]
    log = case["log"]
    moment = report["direction"]
    field = report["field"]
    torque = log["torque"]
    require(all(math.isfinite(value) for value in moment + field + torque),
            f"{label}: non-finite field/torque observable")
    require(norm(moment) > 0.5, f"{label}: magnetic moment was not established")
    expected_torque = cross(moment, field)
    relation_error = norm([a - b for a, b in zip(torque, expected_torque)])
    require(relation_error < 5.0e-3,
            f"{label}: production torque is inconsistent with m x B ({relation_error:.3e})")
    return {"field_norm_T": norm(field), "torque_norm_T": norm(torque),
            "m_cross_B_error_T": relation_error}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    base = root / "tests/scf/cases/bulk/bccFe"
    scratch = args.scratch_root.resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    z = [0.0, 0.0, 1.0]
    x = [1.0, 0.0, 0.0]

    cases: dict[str, dict] = {}
    cases["collinear_soc_off"] = run_case(
        args.binary.resolve(), root, base, scratch / "collinear_soc_off",
        nsp=2, soc_scale=0.0, moment=z)
    cases["global_rotation_z"] = run_case(
        args.binary.resolve(), root, base, scratch / "global_rotation_z",
        nsp=3, soc_scale=0.0, moment=z)
    cases["global_rotation_x"] = run_case(
        args.binary.resolve(), root, base, scratch / "global_rotation_x",
        nsp=3, soc_scale=0.0, moment=x)
    cases["soc_global_rotation_z"] = run_case(
        args.binary.resolve(), root, base, scratch / "soc_global_rotation_z",
        nsp=4, soc_scale=1.0, moment=z)
    cases["soc_global_rotation_x"] = run_case(
        args.binary.resolve(), root, base, scratch / "soc_global_rotation_x",
        nsp=4, soc_scale=1.0, moment=x)
    evidence = {
        "collinear_soc_off": check_production_route(cases["collinear_soc_off"], "collinear SOC-off"),
        "global_rotation_z": check_production_route(cases["global_rotation_z"], "global rotation z"),
        "global_rotation_x": check_production_route(cases["global_rotation_x"], "global rotation x"),
        "soc_global_rotation_z": check_production_route(cases["soc_global_rotation_z"], "SOC rotation z"),
        "soc_global_rotation_x": check_production_route(cases["soc_global_rotation_x"], "SOC rotation x"),
    }

    collinear = cases["collinear_soc_off"]["log"]["torque"]
    require(norm(collinear[:2]) < 5.0e-3,
            "collinear SOC-off transverse torque is not zero within the printed precision")

    rotation_energy_error = abs(cases["global_rotation_z"]["report"]["energy"] -
                                cases["global_rotation_x"]["report"]["energy"])
    require(rotation_energy_error < 2.0e-5,
            f"SOC-free global rotation changed energy by {rotation_energy_error:.3e} Ry")
    rotation_torque_error = norm([
        a - b for a, b in zip(cases["global_rotation_z"]["log"]["torque"],
                              cases["global_rotation_x"]["log"]["torque"])
    ])
    require(rotation_torque_error < 5.0e-3,
            f"SOC-free global rotation changed the torque invariant by {rotation_torque_error:.3e} T")

    soc_torque = cases["soc_global_rotation_x"]["log"]["torque"]
    require(norm(soc_torque) > 1.0e-4,
            "SOC case was incorrectly treated as a zero-torque state")

    result = {
        "scope": "bcc-Fe one-site production SCF/report cases; no trajectory integration",
        "route": "bands%calculate_magnetic_torques -> mag_for (Ry*ry2tesla) and m x I",
        "evidence": evidence,
        "soc_free_global_rotation_energy_error_Ry": rotation_energy_error,
        "soc_free_global_rotation_torque_difference_T": rotation_torque_error,
        "controlled_canting_oracle": "UnitConstrainingField (VAL-10): canted restoring field and penalty finite difference",
        "interpretation": (
            "LMTO field/torque feedback and SOC separation pass in the scoped fixture. "
            "VAL-10 supplies the constrained-field finite-difference penalty oracle; "
            "no independent fixed-potential electronic-energy angle sweep is claimed."
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    print("VAL-12 PASS: LMTO field/torque invariants and SOC separation checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
