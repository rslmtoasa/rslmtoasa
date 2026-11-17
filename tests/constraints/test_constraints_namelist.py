#!/usr/bin/env python3
"""AppConstraintsNamelist{Mom,Malformed}: constraints_mom_ref read from the &constraints namelist.
Deck: example/bulk/bccFe Fe.nml, one SCF step, constraints_i_cons = 3, constraints_diagnostics = .true.
mom: constraints_mom_ref(1:3,1) = 0.6, 0, 0.8 in the documented per-component syntax. Passes if the run exits 0 and
  constraint_diagnostics.dat (iteration 1, site 1) has target_x,y,z = 0.6, 0, 0.8.
  Bug caught: the namelist read fails, constraints stay off and the file is not written (or the run stops at the zero-length target).
malformed: constraints_mom_ref(4,1) = 1.0 (row index outside 1:3). Passes if the run exits non-zero with
  'Error while reading constraints namelist: iostat = ...' in the output.
  Bug caught: the read error is only logged and the run continues unconstrained.
usage: test_constraints_namelist.py <mom|malformed> <rslmto.x> <scratch directory>"""
import shutil
import subprocess
import sys
from pathlib import Path

case, binary, scratch = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3])
root = Path(__file__).resolve().parents[2]
mom_ref = {"mom": "constraints_mom_ref(1,1) = 0.6\nconstraints_mom_ref(2,1) = 0.0\nconstraints_mom_ref(3,1) = 0.8\n",
           "malformed": "constraints_mom_ref(4,1) = 1.0\n"}[case]
deck = f"""&calculation
processing = 'sd'
verbose = T
/
&lattice
rc = 080
alat = 2.86120
crystal_sym = 'bcc'
wav = 1.40880
ntype = 1
ct(1) = 3.0d0
r2 = 9.00d0
strux_backend = 'strux_lib'
/
&atoms
database = './'
label(1) = 'Fe'
/
&self
nstep = 1
/
&energy
fermi = -0.042265
energy_min = -2.0
energy_max = 0.8
channels_ldos = 500
/
&control
calctype = 'B'
nsp = 4
lld = 20
recur = 'block'
/
&constraints
constraints_enable = .true.
constraints_i_cons = 3
constraints_diagnostics = .true.
{mom_ref}/
"""
shutil.rmtree(scratch, ignore_errors=True)
scratch.mkdir(parents=True)
shutil.copy(root / "example/bulk/bccFe/Fe.nml", scratch / "Fe.nml")
(scratch / "input.nml").write_text(deck)
run = subprocess.run([str(binary), "input.nml"], cwd=scratch, capture_output=True, text=True)
out = run.stdout + run.stderr
if case == "malformed":
    ok = run.returncode != 0 and "Error while reading constraints namelist: iostat = " in out
    print(f"exit code {run.returncode}; message found: {'Error while reading constraints namelist: iostat = ' in out}")
else:
    diag = scratch / "constraint_diagnostics.dat"
    target = None
    if diag.exists():
        rows = [ln.split() for ln in diag.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        if rows:
            target = [float(x) for x in rows[0][2:5]]
    ok = run.returncode == 0 and target is not None and all(abs(a - b) < 1e-12 for a, b in zip(target, (0.6, 0.0, 0.8)))
    print(f"exit code {run.returncode}; target(iteration 1, site 1) = {target}; expected [0.6, 0.0, 0.8]")
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
