#!/usr/bin/env python3
"""VAL-04: physical/scope checks for ordinary onsite Liechtenstein U/J.

The campaign deliberately uses the compact bcc-Fe reciprocal-SCF fixture for
the material response.  It runs fresh calculations from the same atomic
starting state, so the U sweep is a controlled finite-SCF response envelope,
not a comparison of unrelated saved fixed points.  The RS/KS common-state
identity is protected separately by UnitLdaUHamiltonian.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import f90nml


RY_PER_EV = 1.0 / 13.605703976
ANSI = re.compile(r"\x1b\[[0-9;]*m")
SPLIT_RE = re.compile(r"HUBBARD_VDIAG.*l=\s*2.*split\(up-dn\)=\s*([-+0-9.EeDd]+)")
RESIDUAL_RE = re.compile(r"(?:Not converged! Diff=|Converged!\s*)\s*([-+0-9.EeDd]+)")
LDM_CHECK_RE = re.compile(
    r"HUBBARD_LDM_CHECK atom=\s*1 max_antiherm=\s*([-+0-9.EeDd]+)"
)


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def run_case(binary: Path, fixture: Path, scratch_root: Path, name: str, u_ev: float | None, nstep: int) -> dict:
    workdir = scratch_root / name
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(fixture, workdir)
    input_path = workdir / "input.nml"
    text = input_path.read_text()
    text = text.replace("nstep = 2", f"nstep = {nstep}")
    text = text.replace("nsp = 1", "nsp = 2")
    if u_ev is not None:
        text += (
            "\n&hamiltonian\n"
            f"hubbard_u_general(1,3) = {u_ev:.12g}\n"
            "hubbard_j_general(1,3) = 0.0\n"
            "hubbard_u_potential_form = 'liechtenstein'\n/\n"
        )
    input_path.write_text(text)

    launcher = Path(__file__).resolve().parents[1] / "run_binary.sh"
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    result = subprocess.run(
        ["/bin/bash", str(launcher), str(binary)], cwd=workdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if result.returncode != 0:
        log = workdir / "testrun.log"
        detail = log.read_text()[-4000:] if log.exists() else result.stdout
        raise RuntimeError(f"{name} failed:\n{detail}")

    log_text = ANSI.sub("", (workdir / "testrun.log").read_text())
    if "fatal" in log_text.lower():
        raise RuntimeError(f"{name} contains a fatal diagnostic")

    output = f90nml.read(workdir / "Fe_out.nml")["par"]
    ql = output["ql"]
    d_up = float(ql[0][2][0])
    d_down = float(ql[1][2][0])
    ldm = output["ldm_flatten"]
    # f90nml reverses the Fortran index order: ldm_flatten(k, spin, l)
    # is represented as ldm_flatten[k][spin][l].
    d_matrices = [[ldm[k][spin][2] for k in range(25)] for spin in range(2)]
    herm_residual = max(
        abs(matrix[i * 5 + j] - matrix[j * 5 + i])
        for matrix in d_matrices for i in range(5) for j in range(5)
    )
    splits = [number(v) for v in SPLIT_RE.findall(log_text)]
    residuals = [number(v) for v in RESIDUAL_RE.findall(log_text)]
    ldm_checks = [number(v) for v in LDM_CHECK_RE.findall(log_text)]

    return {
        "u_ev": 0.0 if u_ev is None else u_ev,
        "u_internal_ry": float(output["hubbard_u"][2]),
        "d_trace_up": sum(d_matrices[0][i * 5 + i] for i in range(5)),
        "d_trace_down": sum(d_matrices[1][i * 5 + i] for i in range(5)),
        "ldm_hermiticity_residual": herm_residual,
        "ldm_checks": ldm_checks,
        "d_moment": d_up - d_down,
        "d_charge": d_up + d_down,
        "etot": float(output["etot"]),
        "sumev": float(output["sumev"]),
        "hubbard_split_ry": splits[-1] if splits else 0.0,
        "residuals": residuals,
        "workdir": str(workdir),
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("VAL-04 FAIL: " + message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path, default=None)
    args = parser.parse_args()

    args.binary = args.binary.resolve()
    root = Path(__file__).resolve().parents[2]
    fixture = root / "tests/scf/cases/k_space_scf/bccFe"
    scratch_root = args.scratch_root or Path(tempfile.mkdtemp(prefix="val04_lda_u_"))
    scratch_root.mkdir(parents=True, exist_ok=True)

    baseline = run_case(args.binary, fixture, scratch_root, "baseline", None, 4)
    u0 = run_case(args.binary, fixture, scratch_root, "u0", 0.0, 4)
    u2 = run_case(args.binary, fixture, scratch_root, "u2", 2.0, 4)
    u4 = run_case(args.binary, fixture, scratch_root, "u4", 4.0, 4)
    u2_conv = run_case(args.binary, fixture, scratch_root, "u2_convergence", 2.0, 8)

    # U=0 recovery is against a calculation with no Hubbard namelist at all.
    require(abs(u0["etot"] - baseline["etot"]) < 1.0e-9, "U=0 total-energy recovery")
    require(abs(u0["d_charge"] - baseline["d_charge"]) < 1.0e-9, "U=0 correlated charge recovery")
    require(abs(u0["u_internal_ry"]) < 1.0e-12, "U=0 internal parameter is zero")

    # Input conversion and the direct physical response of the d channel.
    require(abs(u2["u_internal_ry"] - 2.0 * RY_PER_EV) < 1.0e-10, "2 eV is converted to Ry")
    require(abs(u4["u_internal_ry"] - 4.0 * RY_PER_EV) < 1.0e-10, "4 eV is converted to Ry")
    require(abs(u4["hubbard_split_ry"]) > abs(u2["hubbard_split_ry"]) > 1.0e-3,
            "d-level spin splitting grows with U")
    require(abs(u4["d_moment"] - u0["d_moment"]) > 1.0e-2, "d-shell moment responds to U")
    require(abs(u4["d_charge"] - u0["d_charge"]) > 1.0e-2, "d-shell charge responds to U")

    for result in (u2, u4, u2_conv):
        require(0.0 <= result["d_trace_up"] <= 5.0 and 0.0 <= result["d_trace_down"] <= 5.0,
                "correlated-shell traces remain physical")
        require(result["ldm_hermiticity_residual"] < 1.0e-7, "stored occupation matrix is Hermitian")
        require(result["ldm_checks"], "runtime LDM Hermiticity diagnostic is present")
        require(max(result["ldm_checks"]) < 1.0e-6,
                "Green-function occupation matrix is Hermitian to integration tolerance")
        require(math.isfinite(result["etot"]) and math.isfinite(result["sumev"]),
                "total and band-energy observables are finite")
        require("HUBBARD +V" not in Path(result["workdir"], "testrun.log").read_text(),
                "onsite campaign did not activate Hubbard-V")

    require(u2_conv["residuals"], "SCF residual record is present")
    require(min(u2_conv["residuals"]) < u2_conv["residuals"][0],
            "additional SCF iterations reduce the residual envelope")
    require(abs(u4["etot"] - u0["etot"]) > 1.0e-4, "total energy responds under the active convention")
    require(abs(u4["sumev"] - u0["sumev"]) > 1.0e-4, "band energy responds under the active convention")

    print("VAL-04 PASS: bcc-Fe onsite Liechtenstein U/J response")
    for result in (u0, u2, u4):
        print(
            f"U={result['u_ev']:.1f} eV  d_trace=({result['d_trace_up']:.6f},"
            f" {result['d_trace_down']:.6f})  split={result['hubbard_split_ry']:.6f} Ry"
            f"  d_moment={result['d_moment']:.6f}  d_charge={result['d_charge']:.6f}"
            f"  etot={result['etot']:.9f}  sumev={result['sumev']:.9f}"
        )
    print(
        f"convergence U=2 eV: first={u2_conv['residuals'][0]:.6e}"
        f" min(8 steps)={min(u2_conv['residuals']):.6e}"
    )


if __name__ == "__main__":
    main()
