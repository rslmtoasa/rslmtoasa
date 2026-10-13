"""With nsp = 2 and xi_p = xi_d = 0 the 12^3 Fe response must equal the nsp = 1 response byte for byte; with the deck's xi it must not.

Usage: soc_zero_identity.py <rslmto.x> <deck directory> <scratch directory>
Catches a SOC branch that changes the xi = 0 result, and an nsp = 2 run in which the L.S term is silently absent.
"""
import filecmp
import pathlib
import re
import shutil
import subprocess
import sys

FILES = ["spin_response_q001.dat", "spin_response_q002.dat", "spin_response_q003.dat",
         "spin_response_dispersion.dat", "spin_response_adiabatic.dat"]
BLOCK = ("\n&spin_response\nmethod = 'juelich'\nn_q = 3\nq_list(:, 1) = 0.0, 0.0, 0.0\n"
         "q_list(:, 2) = 0.25, -0.25, 0.25\nq_list(:, 3) = 0.3333333333333333, -0.3333333333333333, 0.3333333333333333\n"
         "omega_min = 0.0\nomega_max = 0.05\nn_omega = 41\neta = 2.0e-3\n/\n")


def run(binary, deck, scratch, nsp, zero_xi):
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.copytree(deck, scratch)
    text = pathlib.Path(deck, "input.nml").read_text()
    old = "pre_processing = 'bravais'"
    if old not in text:
        sys.exit("soc_zero_identity.py: " + old + " not found in " + deck + "/input.nml")
    for key in ("nk1", "nk2", "nk3"):
        text, n = re.subn(rf"^{key} = \d+$", f"{key} = 12", text, flags=re.M)
        if n != 1:
            sys.exit("soc_zero_identity.py: " + key + " not found in " + deck + "/input.nml")
    text, n = re.subn(r"^nsp = 1", f"nsp = {nsp}", text, flags=re.M)
    if n != 1:
        sys.exit("soc_zero_identity.py: nsp = 1 not found in " + deck + "/input.nml")
    if zero_xi:
        atom = pathlib.Path(scratch, "Fe.nml")
        body, n = re.subn(r"^( xi_[pd] =).*$", r"\1 0.0, 0.0", atom.read_text(), flags=re.M)
        if n != 2:
            sys.exit("soc_zero_identity.py: xi_p, xi_d not found in " + str(atom))
        atom.write_text(body)
    nml = pathlib.Path(scratch, "input_driver.nml")
    nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + BLOCK)
    subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def state_lines(scratch):
    """_state.dat without the SOC-only lines."""
    return [ln for ln in pathlib.Path(scratch, "spin_response_state.dat").read_text().splitlines()
            if not ln.startswith(("u_m0_nosoc", "# SOC"))]


binary, deck, scratch = sys.argv[1:]
run(binary, deck, scratch + "_nsp1", 1, False)
run(binary, deck, scratch + "_nsp2_xi0", 2, True)
run(binary, deck, scratch + "_nsp2", 2, False)
same = lambda a, b, f: filecmp.cmp(pathlib.Path(scratch + a, f), pathlib.Path(scratch + b, f), shallow=False)
bad = [f for f in FILES if not same("_nsp1", "_nsp2_xi0", f)]
if state_lines(scratch + "_nsp1") != state_lines(scratch + "_nsp2_xi0"):
    bad.append("spin_response_state.dat")
print("nsp = 1 against nsp = 2 with xi = 0, files differing:", bad if bad else "none")
moved = [f for f in FILES if not same("_nsp1", "_nsp2", f)]
print("nsp = 1 against nsp = 2 with the deck's xi, files differing:", moved if moved else "none")
sys.exit(1 if bad or len(moved) != len(FILES) else 0)
