#!/usr/bin/env python3
"""AppGbtKspaceHamOrderInvalid: an unknown kspace_ham_order must stop the run with a message naming the value and the accepted set, and
reciprocal_mode = 'generalized_overlap_kanpur' must stop it with "not implemented".
Deck: tests/spin_response/fe_bcc. Bug caught: the old warning plus silent fallback (exit 0) in either case.
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
ok = True
for key, value, message in (("kspace_ham_order", "bogus", "kspace_ham_order = 'bogus' is not one of 'auto', 'first', 'second', 'orthogonal'"),
                            ("reciprocal_mode", "generalized_overlap_kanpur", "reciprocal_mode='generalized_overlap_kanpur' is not implemented")):
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["reciprocal"].update({key: value})
    nml.write(str(scratch / "input.nml"), force=True)
    run = subprocess.run([str(binary), "input.nml"], cwd=scratch, capture_output=True, text=True)
    out = run.stdout + run.stderr
    print(f"{key} = {value}: exit code {run.returncode}; message found: {message in out}")
    ok = ok and run.returncode != 0 and message in out
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
