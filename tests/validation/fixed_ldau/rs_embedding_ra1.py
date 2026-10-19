#!/usr/bin/env python3
"""Check one-step RS impurity energy diagnostics against production evaluator."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import f90nml


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "tests/scf/cases/impurity/B2FeCo_nsp1"
LAUNCHER = ROOT / "tests/run_binary.sh"
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--evaluator", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="rs_embedding_ra1_"))
    scratch.mkdir(parents=True, exist_ok=True)
    work = scratch / "impurity_one_step"
    shutil.copytree(FIXTURE, work)
    input_path = work / "input.nml"
    text = input_path.read_text()
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", text)
    input_path.write_text(text)
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    proc = subprocess.run(
        ["/bin/bash", str(LAUNCHER), str(args.binary.resolve())], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    log_path = work / "testrun.log"
    log = log_path.read_text() if log_path.exists() else proc.stdout
    if proc.returncode != 0:
        raise RuntimeError(f"one-step impurity run failed:\n{log[-5000:]}")

    rows = re.findall(
        rf"HUBBARD_MOMENT_ENERGY_SITE site=\s*1 type=\s*(\d+).*?E_int=\s*({NUMBER})"
        rf".*?E_DC=\s*({NUMBER}).*?E_Hub=\s*({NUMBER})",
        log,
    )
    if not rows:
        raise RuntimeError("impurity fixed-U energy diagnostic row is missing")
    type_id, *reported_text = rows[-1]
    # This B2FeCo fixture has two bulk symbolic types and one recursive atom.
    # The diagnostic must use nbulk + ia (=3), not lattice%ib(1), which is
    # the geometric representative index (16 for this cluster).
    if int(type_id) != 3:
        raise RuntimeError(f"impurity diagnostic used type={type_id}; expected symbolic type nbulk+ia=3")
    # The single impurity is represented by the Fe symbolic type, so its
    # serialized recursion potential is Fe_out.nml (the input label retains
    # the `-imp` suffix and does not produce a separate *_out file).
    out_path = work / "Fe_out.nml"
    if not out_path.exists():
        raise RuntimeError("recursion atom potential output Fe_out.nml is missing")
    output = f90nml.read(out_path)["par"]
    ldm = output["ldm_flatten"]
    matrix_path = work / "impurity_p_matrix.txt"
    with matrix_path.open("w") as stream:
        for spin in range(2):
            for i in range(5):
                for j in range(5):
                    stream.write(f"{spin + 1} {i + 1} {j + 1} {float(ldm[(i * 5 + j)][spin][2]):.17e} 0.0\n")
    evaluator = subprocess.run(
        [str(args.evaluator.resolve()), str(matrix_path), "2.0", "0.0"],
        cwd=work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if evaluator.returncode != 0:
        raise RuntimeError(f"production energy evaluator failed:\n{evaluator.stdout}")
    names = ("E_int_Ry", "E_DC_Ry", "E_Hub_Ry")
    calculated = []
    for name in names:
        match = re.search(rf"{name}=\s*({NUMBER})", evaluator.stdout)
        if not match:
            raise RuntimeError(f"evaluator omitted {name}: {evaluator.stdout}")
        calculated.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    reported = [float(value.replace("D", "E").replace("d", "e")) for value in reported_text]
    differences = [abs(a - b) for a, b in zip(reported, calculated)]
    if max(differences) > 1.0e-12:
        raise RuntimeError(
            "impurity diagnostic disagrees with evaluate_hubbard_energy: "
            f"type={type_id}, reported={reported}, independent={calculated}, delta={differences}"
        )
    print(f"RA1 PASS: impurity type={type_id}, max |E diagnostic - evaluate_hubbard_energy|={max(differences):.3e} Ry")


if __name__ == "__main__":
    main()
