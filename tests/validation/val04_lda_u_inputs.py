#!/usr/bin/env python3
"""Exercise the fixed-U namelist boundary on supported scalar-collinear inputs."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/scf/cases/k_space_scf/bccFe"
LAUNCHER = ROOT / "tests/run_binary.sh"
RY_PER_EV = 1.0 / 13.605703976


def run(binary: Path, scratch: Path, name: str, assignments: str, *, succeeds: bool,
        nsp: int = 1, lmax: int = 2, hoh: bool = False, expected_message: str | tuple[str, ...] = "") -> str:
    work = scratch / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(FIXTURE, work)
    inp = work / "input.nml"
    text = inp.read_text()
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", text)
    text = re.sub(r"(?m)^nsp\s*=\s*\d+", f"nsp = {nsp}", text)
    if re.search(r"(?m)^\s*lmax\s*=", text):
        text = re.sub(r"(?m)^\s*lmax\s*=\s*\d+", f"lmax = {lmax}", text, count=1)
    else:
        text = re.sub(r"(?m)^(&control\s*$)", rf"\1\nlmax = {lmax}", text, count=1)
    if not re.search(r"(?m)^\s*lmax\s*=\s*\d+", text):
        raise RuntimeError("failed to set control%lmax in generated input")
    text += f"\n&hamiltonian\nhoh = {'.true.' if hoh else '.false.'}\n{assignments}\n/\n"
    inp.write_text(text)
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    proc = subprocess.run(
        ["/bin/bash", str(LAUNCHER), str(binary)], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    log = work / "testrun.log"
    output = log.read_text() if log.exists() else proc.stdout
    output = re.sub(r"\x1b\[[0-9;]*m", "", output)
    if succeeds:
        if proc.returncode != 0 or "fatal" in output.lower():
            raise RuntimeError(f"{name} should succeed, exit={proc.returncode}:\n{output[-4000:]}")
    else:
        expected = (expected_message,) if isinstance(expected_message, str) else expected_message
        if proc.returncode == 0 or "fatal" not in output.lower() or not any(message.lower() in output.lower() for message in expected):
            raise RuntimeError(f"{name} should be rejected with a fatal diagnostic:\n{output[-4000:]}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch_obj = tempfile.TemporaryDirectory(prefix="val04_lda_u_inputs_") if args.scratch_root is None else None
    scratch = args.scratch_root or Path(scratch_obj.name)
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()

    no_u = run(binary, scratch, "no_u", "", succeeds=True)
    # SOC mode remains available when fixed-U is inactive. This is a regression
    # boundary: the fixed-U guard must not reject unrelated no-U SOC runs.
    no_u_nsp2 = run(binary, scratch, "no_u_nsp2", "", succeeds=True, nsp=2)
    zero_uj = run(binary, scratch, "zero_uj", "hubbard_u_general(1,3) = 0.0\nhubbard_j_general(1,3) = 0.0", succeeds=True)
    local_zero_uj = run(binary, scratch, "local_zero_uj",
                        "hubbard_radial_moments = 'local'\nhubbard_u_general(1,3) = 0.0\nhubbard_j_general(1,3) = 0.0",
                        succeeds=True)
    run(binary, scratch, "local_without_explicit_uj", "hubbard_radial_moments = 'local'", succeeds=False,
        expected_message="requires an explicit fixed Liechtenstein U/J input")
    no_u_spdf = run(binary, scratch, "no_u_spdf", "", succeeds=True, lmax=3)
    zero_uj_spdf = run(binary, scratch, "zero_uj_spdf",
                       "hubbard_u_general(1,3) = 0.0\nhubbard_j_general(1,3) = 0.0",
                       succeeds=True, lmax=3)
    if any("HUBBARD summary:" in output or "HUBBARD fixed U/J" in output
           for output in (no_u, no_u_nsp2, zero_uj, local_zero_uj, no_u_spdf, zero_uj_spdf)):
        raise RuntimeError("no-U and explicit U=J=0 inputs must leave Hubbard inactive")

    u_only = run(binary, scratch, "u_only", "hubbard_u_general(1,3) = 2.0", succeeds=True)
    if "hubbard_j_general" not in u_only.lower() and "HUBBARD fixed U/J" not in u_only:
        raise RuntimeError("U-only input did not reach the fixed-U parser summary")
    # The parser prints internal Ry values. A missing J must become exact zero;
    # U=2 eV must undergo one eV-to-Ry conversion.
    u_line = next((line for line in u_only.splitlines() if "HUBBARD fixed U/J" in line), "")
    u_text = re.search(r"\bU=(.*?)\s+J=", u_line)
    j_text = re.search(r"\bJ=(.*)$", u_line)
    values = lambda part: [float(v.replace("D", "E").replace("d", "e")) for v in re.findall(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?", part)]
    u_values = values(u_text.group(1)) if u_text else []
    j_values = values(j_text.group(1)) if j_text else []
    rounded_u = float(f"{2.0 * RY_PER_EV:.6f}")
    if not u_values or abs(u_values[2] - rounded_u) > 1e-12 or not j_values or any(abs(v) > 1e-12 for v in j_values):
        raise RuntimeError(f"U-only parsing did not produce U=2 eV in Ry and J=0: {u_line}")

    uj = run(binary, scratch, "u_and_j", "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5", succeeds=True)
    if "hubbard_j_general" not in uj.lower() and "HUBBARD fixed U/J" not in uj:
        raise RuntimeError("U+J input did not reach the fixed-U parser summary")
    uj_line = next((line for line in uj.splitlines() if "HUBBARD fixed U/J" in line), "")
    u_text = re.search(r"\bU=(.*?)\s+J=", uj_line)
    j_text = re.search(r"\bJ=(.*)$", uj_line)
    u_values = values(u_text.group(1)) if u_text else []
    j_values = values(j_text.group(1)) if j_text else []
    rounded_u = float(f"{2.0 * RY_PER_EV:.6f}")
    rounded_j = float(f"{0.5 * RY_PER_EV:.6f}")
    if not u_values or abs(u_values[2] - rounded_u) > 1e-12 or not j_values or abs(j_values[2] - rounded_j) > 1e-12:
        raise RuntimeError(f"U+J parsing did not convert eV to Ry once: {uj_line}")

    run(binary, scratch, "j_without_u", "hubbard_j_general(1,3) = 0.5", succeeds=False,
        expected_message="requires a corresponding hubbard_u_general")
    run(binary, scratch, "unsupported_p", "hubbard_u_general(1,2) = 2.0", succeeds=False,
        expected_message="validated only for the d channel")
    run(binary, scratch, "invalid_type_index", "hubbard_u_general(999,3) = 2.0", succeeds=False,
        expected_message="invalid &hamiltonian namelist")
    run(binary, scratch, "unsupported_sp_basis", "hubbard_u_general(1,3) = 2.0", succeeds=False, lmax=1,
        expected_message=("requires at least an spd basis", "channel index exceeds the available potential channels"))
    # The fixed tensor is in the complex spherical m basis. The n=16 spdf
    # hcpx path currently has no Cartesian-to-spherical transformation, so
    # an active d-shell U there would contract cubic-harmonic LDM entries with
    # the wrong tensor basis. No-U spdf behavior is deliberately unaffected.
    run(binary, scratch, "unsupported_spdf_basis", "hubbard_u_general(1,3) = 2.0", succeeds=False, lmax=3,
        expected_message="the spdf d-shell basis transform is unsupported")
    run(binary, scratch, "unsupported_spdf_basis_u_j",
        "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5", succeeds=False, lmax=3,
        expected_message="the spdf d-shell basis transform is unsupported")
    run(binary, scratch, "unsupported_soc", "hubbard_u_general(1,3) = 2.0", succeeds=False, nsp=2,
        expected_message="SOC and noncollinear spinors are unsupported")
    run(binary, scratch, "unsupported_noncollinear", "hubbard_u_general(1,3) = 2.0", succeeds=False, nsp=3,
        expected_message="SOC and noncollinear spinors are unsupported")
    run(binary, scratch, "unsupported_noncollinear_soc", "hubbard_u_general(1,3) = 2.0", succeeds=False, nsp=4,
        expected_message="SOC and noncollinear spinors are unsupported")
    run(binary, scratch, "unsupported_hoh", "hubbard_u_general(1,3) = 2.0", succeeds=False, hoh=True,
        expected_message="unsupported with hoh=.true.")
    print("VAL-04 PASS: fixed-U inputs/units/guards; local mode requires explicit U/J and accepts zero-U no-op")


if __name__ == "__main__":
    main()
