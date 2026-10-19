#!/usr/bin/env python3
"""Exercise fixed-U potential provenance, mismatch and legacy behavior."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "tests/scf/cases/k_space_scf/bccFe"
LAUNCHER = ROOT / "tests/run_binary.sh"


def configure(work: Path, *, backend: str, mode: str | None, uj: bool = True) -> None:
    path = work / "input.nml"
    text = path.read_text()
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", text)
    if re.search(r"(?m)^\s*use_kspace\s*=", text):
        text = re.sub(r"(?m)^\s*use_kspace\s*=.*$", f"use_kspace = {'.true.' if backend == 'ks' else '.false.'}", text)
    else:
        text = re.sub(r"(?m)^&self\s*$", "&self\nuse_kspace = " + (".true." if backend == "ks" else ".false."), text, count=1)
    radial = f"hubbard_radial_moments = '{mode}'\n" if mode else ""
    hamiltonian = "&hamiltonian\n" + radial
    if uj:
        hamiltonian += "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n"
    hamiltonian += "/\n"
    path.write_text(text + "\n" + hamiltonian)


def run(binary: Path, work: Path) -> tuple[int, str]:
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    proc = subprocess.run(
        ["/bin/bash", str(LAUNCHER), str(binary)], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    log = work / "testrun.log"
    output = log.read_text() if log.exists() else proc.stdout
    output = re.sub(r"\x1b\[[0-9;]*m", "", output)
    return proc.returncode, output


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh_work(root: Path, name: str) -> Path:
    work = root / name
    shutil.copytree(FIXTURE, work)
    return work


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="rs_embedding_ra2_"))
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()

    # A fresh U run writes its resolved mode and backend. This fixture is a
    # reciprocal KS bulk run using an ordinary no-U atomic seed.
    seed = fresh_work(scratch, "ks_local_seed")
    configure(seed, backend="ks", mode="local")
    status, log = run(binary, seed)
    if status != 0:
        raise RuntimeError(f"KS local seed failed:\n{log[-4000:]}")
    potential = (seed / "Fe_out.nml").read_text().lower()
    if "hubbard_radial_moments" not in potential or "'local'" not in potential:
        raise RuntimeError("U-carrying potential output omitted its resolved local mode")
    if "hubbard_backend" not in potential or "'ks'" not in potential:
        raise RuntimeError("U-carrying potential output omitted its KS backend")

    # A matching local-mode RS restart is accepted; the backend is informational.
    matching = fresh_work(scratch, "rs_matching_local")
    shutil.copy2(seed / "Fe_out.nml", matching / "Fe_out.nml")
    configure(matching, backend="rs", mode="local")
    status, log = run(binary, matching)
    if status != 0:
        raise RuntimeError(f"matching KS-host/RS local-mode run failed:\n{log[-4000:]}")

    # Backend provenance is informational: an RS consumer must accept the
    # same local-mode host even when the backend field is absent.
    missing_backend = fresh_work(scratch, "rs_missing_backend_local")
    backendless_text = re.sub(r"(?im)^\s*hubbard_backend\s*=.*\n", "", (seed / "Fe_out.nml").read_text())
    (missing_backend / "Fe_out.nml").write_text(backendless_text)
    configure(missing_backend, backend="rs", mode="local")
    status, log = run(binary, missing_backend)
    if status != 0:
        raise RuntimeError(f"local-mode potential without backend provenance was rejected:\n{log[-4000:]}")

    # Auto resolves to spectral for RS bulk. A KS-local host must stop here.
    mismatch = fresh_work(scratch, "rs_mode_mismatch")
    shutil.copy2(seed / "Fe_out.nml", mismatch / "Fe_out.nml")
    configure(mismatch, backend="rs", mode=None)
    status, log = run(binary, mismatch)
    if status == 0 or "fixed-u potential moment-mode mismatch" not in log.lower():
        raise RuntimeError(f"local-host/spectral-RS mismatch was not rejected:\n{log[-4000:]}")

    # A legacy U-carrying file without the new fields is classified as spectral.
    legacy = fresh_work(scratch, "rs_legacy_spectral")
    legacy_text = (seed / "Fe_out.nml").read_text()
    legacy_text = re.sub(r"(?im)^\s*hubbard_radial_moments\s*=.*\n", "", legacy_text)
    legacy_text = re.sub(r"(?im)^\s*hubbard_backend\s*=.*\n", "", legacy_text)
    (legacy / "Fe_out.nml").write_text(legacy_text)
    configure(legacy, backend="rs", mode="spectral")
    status, log = run(binary, legacy)
    if status != 0 or "has no hubbard_radial_moments field; treating it as spectral" not in log.lower():
        raise RuntimeError(f"legacy potential did not warn and run in spectral mode:\n{log[-4000:]}")

    # No-U potential output must not acquire provenance fields.
    no_u = fresh_work(scratch, "no_u")
    configure(no_u, backend="ks", mode=None, uj=False)
    status, log = run(binary, no_u)
    if status != 0:
        raise RuntimeError(f"no-U run failed:\n{log[-4000:]}")
    no_u_text = (no_u / "Fe_out.nml").read_text().lower()
    if "hubbard_radial_moments" in no_u_text or "hubbard_backend" in no_u_text:
        raise RuntimeError("no-U output gained fixed-U provenance fields")

    scenarios = {}
    for name in ("ks_local_seed", "rs_matching_local", "rs_missing_backend_local", "rs_mode_mismatch",
                 "rs_legacy_spectral", "no_u"):
        work = scratch / name
        record = {"input": str(work / "input.nml"), "input_sha256": sha256(work / "input.nml"),
                  "log": str(work / "testrun.log"), "log_sha256": sha256(work / "testrun.log")}
        potential_path = work / "Fe_out.nml"
        if potential_path.exists():
            record["potential_output"] = str(potential_path)
            record["potential_output_sha256"] = sha256(potential_path)
        scenarios[name] = record
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unavailable"
    manifest = {
        "stage": "RA2",
        "source_commit": commit,
        "executable": str(binary),
        "executable_sha256": sha256(binary),
        "command": ["python3", str(Path(__file__).resolve()), "--binary", str(binary), "--scratch-root", str(scratch)],
        "raw_output_location": str(scratch),
        "scenarios": scenarios,
    }
    (scratch / "RA2_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("RA2 PASS: stored mode/backend, cross-backend restart, missing backend accepted, mode mismatch fatal, "
          "legacy warning, no-U omission")


if __name__ == "__main__":
    main()
