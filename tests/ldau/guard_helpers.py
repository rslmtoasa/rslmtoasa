"""Deck helpers for the LDA+U guard scenarios that need a host run (host-nsp and host-axis mismatch)."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HOST_FIXTURE = ROOT / "tests/scf/cases/k_space_scf/bccFe"
IMPURITY_FIXTURE = ROOT / "tests/scf/cases/impurity/B2FeCo_nsp1"
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}


def set_line(path: Path, pattern: str, replacement: str) -> None:
    text = path.read_text()
    new, n = re.subn(pattern, replacement, text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError(f"{path}: pattern {pattern!r} not found")
    path.write_text(new)


def fresh(fixture: Path, root: Path, name: str) -> Path:
    work = root / name
    shutil.copytree(fixture, work)
    return work


def configure_host(work: Path) -> None:
    """One-step KS bulk Fe host, local radial moments, U=2 J=0.5 on type 3."""
    path = work / "input.nml"
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", path.read_text())
    if re.search(r"(?m)^\s*use_kspace\s*=", text):
        text = re.sub(r"(?m)^\s*use_kspace\s*=.*$", "use_kspace = .true.", text)
    else:
        text = re.sub(r"(?m)^&self\s*$", "&self\nuse_kspace = .true.", text, count=1)
    hamiltonian = ("&hamiltonian\nhubbard_radial_moments = 'local'\n"
                   "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n/\n")
    path.write_text(text + "\n" + hamiltonian)


def configure_impurity(work: Path, *, impurity_u: bool) -> None:
    path = work / "input.nml"
    text = path.read_text()
    text = re.sub(r"(?ims)(^\s*&hamiltonian\b.*?)(^\s*/\s*$)",
                  lambda m: re.sub(r"(?im)^\s*hubbard_radial_moments\s*=.*\n", "", m.group(1))
                  + "hubbard_radial_moments = 'local'\n" + m.group(2), text, count=1)
    if not re.search(r"(?im)^\s*hubbard_u_general\(1,3\)", text):
        text = re.sub(r"(?im)^\s*&hamiltonian\s*$",
                      "&hamiltonian\nhubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n", text, count=1)
    path.write_text(text)
    set_line(path, r"^nsp\s*=\s*\d+.*$", "nsp = 2")
    if impurity_u:
        text = path.read_text()
        path.write_text(re.sub(r"(?im)^\s*&hamiltonian\s*$",
                               "&hamiltonian\nhubbard_u_impurity(1,3) = 2.0\nhubbard_j_impurity(1,3) = 0.5\n",
                               text, count=1))


def copy_host_potential(host: Path, destination: Path, name: str = "Fe.nml") -> None:
    """Copy the host potential and, when it names one, its dm sidecar."""
    shutil.copy2(host / "Fe_out.nml", destination / name)
    match = re.search(r"(?m)^\s*dm_file\s*=\s*'([^']*)'", (host / "Fe_out.nml").read_text())
    if match:
        shutil.copy2(host / match.group(1), destination / match.group(1))


def run_binary(binary: Path, work: Path, timeout: int = 600) -> tuple[int, str]:
    proc = subprocess.run([str(binary), "input.nml"], cwd=work, env={**os.environ, **SINGLE_THREAD},
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    return proc.returncode, proc.stdout
