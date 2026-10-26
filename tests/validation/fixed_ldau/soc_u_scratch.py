"""Scratch copy of the bcc-Fe nsp=2 deck with two inputs for the spinor-U insertion test.

Usage: soc_u_scratch.py <deck directory> <scratch directory>
input_base.nml is the deck unchanged; input_u.nml adds fixed U=2 eV, J=0.5 eV and
hubbard_soc_experimental to the existing &hamiltonian group.
"""
import pathlib
import shutil
import sys

deck, scratch = sys.argv[1:3]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
text = pathlib.Path(deck, "input.nml").read_text()
group = "&hamiltonian\n"
if text.count(group) != 1:
    sys.exit("soc_u_scratch.py: expected one &hamiltonian group in " + deck + "/input.nml")
extra = "hubbard_soc_experimental = .true.\nhubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n"
pathlib.Path(scratch, "input_base.nml").write_text(text)
pathlib.Path(scratch, "input_u.nml").write_text(text.replace(group, group + extra))
