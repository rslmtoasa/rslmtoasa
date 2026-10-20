#!/usr/bin/env python3
"""AppGbtKspaceHamOrderInvalid: an unknown kspace_ham_order must stop the run with a message naming the value and the accepted set.
Deck: tests/spin_response/fe_bcc with kspace_ham_order = 'bogus'. Bug caught: the old warning plus silent fallback to 'auto' (exit 0).
usage: test_kspace_ham_order_invalid.py <rslmto.x> <scratch directory>"""
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

from gbt_app_common import SR

binary, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(SR / "fe_bcc", scratch)
nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
nml["reciprocal"].update(kspace_ham_order="bogus")
nml.write(str(scratch / "input.nml"), force=True)
run = subprocess.run([str(binary), "input.nml"], cwd=scratch, capture_output=True, text=True)
out = run.stdout + run.stderr
ok = run.returncode != 0 and "kspace_ham_order = 'bogus' is not one of 'auto', 'first', 'second', 'orthogonal'" in out
print(f"exit code {run.returncode}; message found: {'bogus' in out}")
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
