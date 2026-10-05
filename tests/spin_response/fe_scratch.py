"""Copy the bcc-Fe deck to a scratch directory and run post_processing='spin_response' there.

Usage: fe_scratch.py <rslmto.x> <deck directory> <scratch directory>
The deck directory is only read.
"""
import pathlib
import shutil
import subprocess
import sys

# q = 0 and two mesh-commensurate q along Gamma-H (xi = 1/6, 1/3 on the 12x12x12 mesh); eta = 5e-4 Ry.
SPIN_RESPONSE = """
&spin_response
n_q = 3
q_list(:, 1) = 0.0, 0.0, 0.0
q_list(:, 2) = -0.083333333333333333, 0.083333333333333333, 0.083333333333333333
q_list(:, 3) = -0.16666666666666667, 0.16666666666666667, 0.16666666666666667
omega_min = -0.004
omega_max = 0.02
n_omega = 193
eta = 5.0e-4
/
"""

binary, deck, scratch = sys.argv[1:4]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
old = "pre_processing = 'bravais'"
text = pathlib.Path(deck, "input.nml").read_text()
if old not in text:
    sys.exit("fe_scratch.py: " + old + " not found in " + deck + "/input.nml")
nml = pathlib.Path(scratch, "input_driver.nml")
nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + SPIN_RESPONSE)
subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, "driver.log"), "w"),
               stderr=subprocess.STDOUT)
