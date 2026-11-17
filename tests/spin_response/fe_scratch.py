"""Copy the bcc-Fe deck to a scratch directory and run post_processing='spin_response' there.

Usage: fe_scratch.py <rslmto.x> <deck directory> <scratch directory>
The deck directory is only read. The response mesh is set to 12^3 (the deck's nk is its SCF mesh).
"""
import pathlib
import re
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

# The same two q (q = 0 and xi = 1/6) as the static value: one omega = 0, eta = 0, files static_q001.dat, static_q002.dat.
STATIC = """
&spin_response
output_prefix = 'static'
n_q = 2
q_list(:, 1) = 0.0, 0.0, 0.0
q_list(:, 2) = -0.083333333333333333, 0.083333333333333333, 0.083333333333333333
omega_min = 0.0
n_omega = 1
eta = 0.0
/
"""

binary, deck, scratch = sys.argv[1:4]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
old = "pre_processing = 'bravais'"
text = pathlib.Path(deck, "input.nml").read_text()
if old not in text:
    sys.exit("fe_scratch.py: " + old + " not found in " + deck + "/input.nml")
for key in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{key} = \d+$", f"{key} = 12", text, flags=re.M)
    if n != 1:
        sys.exit("fe_scratch.py: " + key + " not found in " + deck + "/input.nml")
pathlib.Path(scratch, "input.nml").write_text(text)  # the Fortran test programs read ./input.nml
# the same two runs with channels = 'spd' (output prefixes spd and spd_static)
SPD = SPIN_RESPONSE.replace("&spin_response\n", "&spin_response\nchannels = 'spd'\noutput_prefix = 'spd'\n")
SPD_STATIC = STATIC.replace("output_prefix = 'static'", "channels = 'spd'\noutput_prefix = 'spd_static'")
for name, block in (("driver", SPIN_RESPONSE), ("static", STATIC), ("spd", SPD), ("spd_static", SPD_STATIC)):
    nml = pathlib.Path(scratch, f"input_{name}.nml")
    nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + block)
    subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, f"{name}.log"), "w"),
                   stderr=subprocess.STDOUT)
