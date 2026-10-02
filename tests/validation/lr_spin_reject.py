#!/usr/bin/env python3
"""Verify that unsupported reference spin structures fail before SCF setup."""
from pathlib import Path
import subprocess
import sys
import tempfile

binary, formulation, spin = sys.argv[1:]
nsp = {'soc': 2, 'noncollinear': 3}[spin]
route = "formulation='rotation', n_q_points=2, q_list=0,0,0,0.05,0,0"
if formulation == 'tddft':
    route = "formulation='tddft', representation='product_compact', interaction='alsda', n_q=1"
elif formulation == 'projected':
    route = "formulation='projected', interaction='mills_1u', projection='d', n_q=1"
expected = 'SOC transverse response' if spin == 'soc' else 'noncollinear reference response'
with tempfile.TemporaryDirectory(prefix='lr-spin-reject-') as directory:
    case = Path(directory)
    (case / 'input.nml').write_text(f"&calculation pre_processing='bravais', post_processing='linear_response' /\n"
                                  f"&control calctype='B', recur='lanczos', nsp={nsp} /\n&linear_response {route} /\n")
    result = subprocess.run([str(Path(binary).resolve()), 'input.nml'], cwd=case,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
    if result.returncode == 0 or expected not in result.stdout:
        raise AssertionError(f'expected early {expected} rejection; exit={result.returncode}\n{result.stdout}')
    if (case / 'kspace_scf_state.dat').exists():
        raise AssertionError('unsupported spin request reached SCF')
print(f'PASS: {formulation}: early {spin} rejection')
