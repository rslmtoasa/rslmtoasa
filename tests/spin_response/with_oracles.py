#!/usr/bin/env python3
"""Run a command with a freshly generated spin-response oracle directory.

usage: with_oracles.py OUTDIR -- COMMAND... ({oracles} in COMMAND is replaced by OUTDIR)

OUTDIR receives a copy of tests/spin_response/oracles (tolerances, LSWT tables, generators) and the
reference .nml files written by the three generators there (numpy only, about 3 s in total).
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ORACLES = Path(__file__).resolve().parent / "oracles"
GENERATORS = ("make_references.py", "make_references_3d.py", "make_references_nc.py")


def main() -> int:
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        print(__doc__)
        return 2
    out = Path(sys.argv[1])
    command = sys.argv[3:]
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(ORACLES, out)
    for generator in GENERATORS:
        done = subprocess.run([sys.executable, str(ORACLES / generator), str(out)],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if done.returncode != 0:
            print(f"oracle generator {generator} failed:\n{done.stdout}")
            return done.returncode
        if "FAIL" in done.stdout:
            print(f"oracle generator {generator} reported a failed self-check:\n{done.stdout}")
            return 1
    return subprocess.run([str(out) if a == "{oracles}" else a for a in command]).returncode


if __name__ == "__main__":
    sys.exit(main())
