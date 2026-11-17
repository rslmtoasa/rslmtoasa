"""Copy the bcc-Fe deck to a scratch directory as a two-atom cell (A1 = a1, A2 = a2, A3 = 2 a3, the same Fe.nml on both sites)
and run post_processing = 'spin_response' there, which writes the structure files the Fortran S2' test reads.

Usage: af_scratch.py <rslmto.x> <deck directory> <scratch directory>
The deck directory is only read. The response mesh is 6 x 6 x 6.
"""
import pathlib
import re
import shutil
import subprocess
import sys

TAU_CART = [(0, 0, 0), (0.5, 0.5, -0.5)]   # units of alat

binary, deck, scratch = sys.argv[1:4]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
text = pathlib.Path(deck, "input.nml").read_text()
for old, new in (("pre_processing = 'bravais'", "pre_processing = 'none'\npost_processing = 'spin_response'"),
                 ("crystal_sym = 'bcc'", "crystal_sym = 'file'"), ("ntype = 1", "ntype = 2"), ("ct(:) = 4.0d0", "ct(:) = 2*4.0d0"),
                 ("label(1) = 'Fe'", "label(1) = 'Fe'\nlabel(2) = 'Fe'")):
    if old not in text:
        sys.exit("af_scratch.py: " + old + " not found in " + deck + "/input.nml")
    text = text.replace(old, new)
for key in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{key} = \d+$", f"{key} = 6", text, flags=re.M)
    if n != 1:
        sys.exit("af_scratch.py: " + key + " not found in " + deck + "/input.nml")
lat = "&lattice\n nbulk_bulk = 2\n ntot = 2\n nbas = 2\n nrec = 2\n a(:, 1) = -0.5, 0.5, 0.5\n a(:, 2) = 0.5, -0.5, 0.5\n a(:, 3) = 1.0, 1.0, -1.0\n"
lat += "".join(f" crd(:, {i + 1}) = {t[0]}, {t[1]}, {t[2]}\n" for i, t in enumerate(TAU_CART))
lat += "".join(f" {k}({i}) = {i}\n" for k in ("izp", "no", "iu", "ib", "irec") for i in (1, 2)) + "/\n"
pathlib.Path(scratch, "lattice.nml").write_text(lat)
block = "\n&spin_response\nmethod = 'juelich'\nn_q = 1\nq_list(:, 1) = 0.0, 0.0, 0.0\nn_omega = 1\neta = 0.0\n/\n"
pathlib.Path(scratch, "input.nml").write_text(text + block)
subprocess.run([binary, "input.nml"], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, "driver.log"), "w"),
               stderr=subprocess.STDOUT)
