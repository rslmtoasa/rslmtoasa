"""Copy the bcc-Fe deck to a scratch directory and run post_processing='spin_response' there.

Usage: fe_scratch.py <rslmto.x> <deck directory> <scratch directory>
The deck directory is only read.
"""
import pathlib
import shutil
import subprocess
import sys

binary, deck, scratch = sys.argv[1:4]
shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
old = "pre_processing = 'bravais'"
text = pathlib.Path(deck, "input.nml").read_text()
if old not in text:
    sys.exit("fe_scratch.py: " + old + " not found in " + deck + "/input.nml")
nml = pathlib.Path(scratch, "input_driver.nml")
nml.write_text(text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'"))
subprocess.run([binary, nml.name], cwd=scratch, check=True, stdout=open(pathlib.Path(scratch, "driver.log"), "w"),
               stderr=subprocess.STDOUT)
