#!/usr/bin/env python3
"""Regression for hcp Co impurity cluster cutoff-array sizing."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()

    fixture = (
        args.repo_root
        / "docs/validation/fixed_ldau/rs/RA4_assets/co_impurity_templates/u0_from_u0_seed"
    )
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    run_dir = args.scratch_root / "run"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir()

    input_text = (fixture / "input.template.nml").read_text()
    input_text, nstep_count = re.subn(r"(?m)^nstep\s*=\s*\d+", "nstep = 0", input_text)
    if nstep_count != 1:
        raise RuntimeError(f"expected one nstep setting, found {nstep_count}")
    (run_dir / "input.nml").write_text(input_text)

    lattice_text = (fixture / "lattice.nml").read_text()
    lattice_text, ntype_count = re.subn(r"(?mi)^\s*ntype\s*=\s*\d+\s*\n", "", lattice_text)
    if ntype_count != 1:
        raise RuntimeError(f"expected one stale host ntype setting, found {ntype_count}")
    (run_dir / "lattice.nml").write_text(lattice_text)

    for name in ("Co-imp.nml", "Co1-bulk.nml", "Co2-bulk.nml"):
        shutil.copy2(fixture / name, run_dir / name)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(
        [str(args.binary.resolve())],
        cwd=run_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (args.scratch_root / "console.log").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"hcp Co impurity setup exited {result.returncode}; see "
            f"{args.scratch_root / 'console.log'}\n{result.stdout[-4000:]}"
        )

    outnewclu = run_dir / "outnewclu"
    if not outnewclu.is_file():
        raise RuntimeError("outnewclu was not written; impurity cluster setup did not complete")
    text = outnewclu.read_text(errors="replace")
    header = re.search(r"NTYPE\s*=\s*(\d+).*?NBAS\s*=\s*(\d+).*?NREC\s*=\s*(\d+)", text)
    if not header:
        raise RuntimeError(f"could not parse outnewclu control header:\n{text[-2000:]}")
    ntype, nbas, nrec = map(int, header.groups())
    if ntype != 3 or nbas < 1 or nrec != 1:
        raise RuntimeError(f"unexpected hcp Co cluster sizes: ntype={ntype}, nbas={nbas}, nrec={nrec}")
    if "50000" in text:
        raise RuntimeError("ndim leaked into the emitted cluster type metadata")

    print(f"hcp Co impurity cluster setup passed: ntype={ntype}, nbas={nbas}, nrec={nrec}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
