"""Scratch copy of the Ni(001) U=3 surface deck at nsp=2 with two inputs for the surface spinor-U operator test.

Usage: soc_u_surface_scratch.py <deck directory> <scratch directory>
input_u.nml is the deck at nsp=2 (U=3 eV, J=0.7 eV on types 1, 3, 4); input_base.nml is the same without U, from base/ (copies of the potential files with the stored host U removed).
Reduced lattice and recursion settings: only the operator build is exercised, no SCF.
"""
import pathlib
import re
import shutil
import sys

deck, scratch = sys.argv[1:3]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
text = pathlib.Path(deck, "input.nml").read_text()
text = re.sub(r"(?m)^nsp\s*=.*$", "nsp = 2", text)
text = re.sub(r"(?m)^rc\s*=.*$", "rc = 100\nstrux_backend = 'legacy'", text)
if "nsp = 2" not in text or "legacy" not in text:
    sys.exit("soc_u_surface_scratch.py: input.nml has no nsp or rc line")
pathlib.Path(scratch, "input_u.nml").write_text(text)
base = re.sub(r"(?m)^hubbard_(?:[uj]_general|radial_moments).*\n", "", text)
base = base.replace("database = './'", "database = './base/'")
pathlib.Path(scratch, "input_base.nml").write_text(base)
# The frozen host file carries U=3 eV; the base build needs the same host potential without U.
shutil.copytree(scratch, pathlib.Path(scratch, "base"), ignore=shutil.ignore_patterns("base", "input*.nml"))
host = pathlib.Path(scratch, "base", "Ni.nml")
host.write_text(re.sub(r"(?m)^ (?:hubbard_\w+|dm_file|dm_sha256) = .*\n", "", host.read_text()))
