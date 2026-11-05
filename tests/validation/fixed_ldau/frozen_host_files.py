#!/usr/bin/env python3
"""Embedded Ni(001) U=3 run must not touch its frozen host files, and a rerun must pass provenance."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DECK = ROOT / "tests/validation/fixed_ldau/decks/ni001_u3"
LAUNCHER = ROOT / "tests/run_binary.sh"
INPUTS = ["input.nml", "ES.nml", "Ni-S.nml", "Ni-S-1.nml", "Ni.nml", "Ni_dm.nml"]
HOST_INPUTS = ["Ni.nml", "Ni_dm.nml"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary: Path, work: Path) -> str:
    proc = subprocess.run(["/bin/bash", str(LAUNCHER), str(binary.resolve())], cwd=work,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_path = work / "testrun.log"
    log = log_path.read_text() if log_path.exists() else proc.stdout
    if proc.returncode != 0:
        raise RuntimeError(f"run failed:\n{log[-5000:]}")
    return log


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="frozen_host_files_"))
    work = scratch / "ni001_u3_two_steps"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for name in INPUTS:
        shutil.copy(DECK / name, work / name)
    # Two steps; reduced recursion and the legacy strux backend (strux_lib exhausts memory on this deck).
    text = (work / "input.nml").read_text()
    text = re.sub(r"(?m)^nstep\s*=.*$", "nstep = 2", text)
    text = re.sub(r"(?m)^lld\s*=.*$", "lld = 50", text)
    text = re.sub(r"(?m)^channels_ldos\s*=.*$", "channels_ldos = 500", text)
    text = re.sub(r"(?m)^rc\s*=.*$", "rc = 100\nstrux_backend = 'legacy'", text)
    (work / "input.nml").write_text(text)

    before = {name: sha256(work / name) for name in INPUTS}
    run(args.binary, work)
    changed = [name for name in INPUTS if sha256(work / name) != before[name]]
    if changed:
        raise RuntimeError(f"input files changed by the run: {changed}")
    updated = [name for name in ("Ni-S_out.nml", "Ni-S-1_dm.nml", "ES_scf.nml") if not (work / name).exists()]
    if updated:
        raise RuntimeError(f"SCF-updated atom files missing: {updated}")

    log = run(args.binary, work)
    if not re.search(r"provenance=PASS", log):
        raise RuntimeError("rerun in the same directory did not report provenance=PASS")
    changed = [name for name in HOST_INPUTS if sha256(work / name) != before[name]]
    if changed:
        raise RuntimeError(f"host input files changed by the rerun: {changed}")
    print("frozen host files unchanged after run and rerun; rerun provenance=PASS")


if __name__ == "__main__":
    main()
