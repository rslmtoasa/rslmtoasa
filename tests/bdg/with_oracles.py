#!/usr/bin/env python3
"""Run a command with a freshly generated BdG oracle directory.

usage: with_oracles.py OUTDIR -- COMMAND... ({oracles} in COMMAND is replaced by OUTDIR)

OUTDIR receives a copy of tests/bdg/oracles and the reference files written there by make_references.py.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ORACLES = Path(__file__).resolve().parent / "oracles"


def main() -> int:
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        print(__doc__)
        return 2
    out = Path(sys.argv[1])
    command = sys.argv[3:]
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(ORACLES, out)
    done = subprocess.run([sys.executable, str(ORACLES / "make_references.py"), str(out)],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if done.returncode != 0:
        print(f"oracle generator failed:\n{done.stdout}")
        return done.returncode
    if "FAIL" in done.stdout:
        print(f"oracle generator reported a failed self-check:\n{done.stdout}")
        return 1
    return subprocess.run([str(out) if a == "{oracles}" else a for a in command]).returncode


if __name__ == "__main__":
    sys.exit(main())
