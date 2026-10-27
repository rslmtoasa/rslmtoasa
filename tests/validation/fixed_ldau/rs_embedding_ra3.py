#!/usr/bin/env python3
"""Exercise scalar U fixtures and retained nsp=2 guard decks."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
LAUNCHER = ROOT / "tests/run_binary.sh"


def run(binary: Path, fixture: Path, work: Path, *, active_impurity: bool = False,
        radial_mode: str | None = None, should_pass: bool, expected_error: str | None = None) -> str:
    shutil.copytree(fixture, work)
    path = work / "input.nml"
    text = path.read_text()
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", text)
    additions = []
    if active_impurity:
        additions.extend(("hubbard_u_impurity(1,3) = 2.0", "hubbard_j_impurity(1,3) = 0.0"))
    if radial_mode:
        additions.append(f"hubbard_radial_moments = '{radial_mode}'")
    if additions:
        match = re.search(r"(?ims)^\s*&hamiltonian\b(.*?)^\s*/\s*$", text)
        if match:
            content = match.group(1).rstrip() + "\n" + "\n".join(additions) + "\n"
            text = text[:match.start()] + "&hamiltonian\n" + content + "/\n" + text[match.end():]
        else:
            text += "\n&hamiltonian\n" + "\n".join(additions) + "\n/\n"
    path.write_text(text)
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    proc = subprocess.run(
        ["/bin/bash", str(LAUNCHER), str(binary)], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    log_path = work / "testrun.log"
    log = log_path.read_text() if log_path.exists() else proc.stdout
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if should_pass and (proc.returncode != 0 or "fatal" in log.lower()):
        raise RuntimeError(f"{fixture.name} should pass:\n{log[-4000:]}")
    if not should_pass and (proc.returncode == 0 or (expected_error and expected_error.lower() not in log.lower())):
        raise RuntimeError(f"{fixture.name} should hit a documented fixed-U input guard:\n{log[-4000:]}")
    return log


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="rs_embedding_ra3_"))
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    cases = (
        (ROOT / "tests/scf/cases/bulk/bccFe_hubbard_nsp1", "nsp1_bulk", False, None, True, None),
        (ROOT / "tests/scf/cases/impurity/B2FeCo_nsp1", "nsp1_impurity", False, None, True, None),
        (ROOT / "tests/scf/cases/bulk/bccFe_hubbard", "nsp2_bulk", False, None, True, None),
        (ROOT / "tests/scf/cases/bulk/bccFe_hubbard", "nsp2_bulk_local", False, "local", True, None),
        (ROOT / "tests/scf/cases/impurity/B2FeCo", "nsp2_impurity_soc_guard", True, None, False,
         "SOC+U embedded runs (calctype I or S) are not yet supported"),
        (ROOT / "tests/scf/cases/impurity/B2FeCo", "nsp2_impurity_mode_guard", True, "local", False,
         "SOC+U embedded runs (calctype I or S) are not yet supported"),
    )
    for fixture, name, active_impurity, radial_mode, should_pass, expected_error in cases:
        run(binary, fixture, scratch / name, active_impurity=active_impurity, radial_mode=radial_mode,
            should_pass=should_pass, expected_error=expected_error)
    print("RA3 PASS: nsp=1 and nsp=2 bulk and nsp=1 impurity U fixtures run; nsp=2 impurity SOC and local-mode guards are exercised")


if __name__ == "__main__":
    main()
