"""Run the 12^3 Fe driver case with OMP_NUM_THREADS = 1 and 4; the q files and the state file must be byte-identical.

Usage: omp_identity.py <rslmto.x> <deck directory> <scratch directory>
Catches a race or a thread-dependent summation order in accumulate_chi0.
"""
import filecmp
import os
import pathlib
import shutil
import subprocess
import sys

Q_T = [3 / 12, 4 / 12]
FILES = ["spin_response_q001.dat", "spin_response_q002.dat", "spin_response_state.dat"]


def run(binary, deck, scratch, threads):
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.copytree(deck, scratch)
    text = pathlib.Path(deck, "input.nml").read_text()
    old = "pre_processing = 'bravais'"
    if old not in text:
        sys.exit("omp_identity.py: " + old + " not found in " + deck + "/input.nml")
    q = "".join(f"q_list(:, {i + 1}) = {t!r}, {-t!r}, {t!r}\n" for i, t in enumerate(Q_T))
    block = f"\n&spin_response\nmethod = 'juelich'\nn_q = 2\n{q}omega_min = 0.0\nomega_max = 0.05\nn_omega = 41\neta = 2.0e-3\n/\n"
    nml = pathlib.Path(scratch, "input_driver.nml")
    nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'") + block)
    # the threaded MKL eigensolver is not bit-reproducible across thread counts and is not under test
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS="1")
    subprocess.run([binary, nml.name], cwd=scratch, check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


binary, deck, scratch = sys.argv[1:]
run(binary, deck, scratch + "_1", 1)
run(binary, deck, scratch + "_4", 4)
bad = [f for f in FILES if not filecmp.cmp(pathlib.Path(scratch + "_1", f), pathlib.Path(scratch + "_4", f), shallow=False)]
print("files differing between 1 and 4 threads:", bad if bad else "none")
sys.exit(1 if bad else 0)
