#!/usr/bin/env python3
"""Run a pinned-main transport_operator_test on the fixed pytest input fixtures.
Usage: compare_open_main.py BASELINE_EXECUTABLE PYTEST_BASEDIR OUTPUT_CSV
Requires test_compiled_actions to have completed in PYTEST_BASEDIR.
"""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np

exe,fixtures,output=map(Path,sys.argv[1:])
rows=[]
for fixture in sorted(fixtures.glob('test_compiled_actions*')):
    if fixture.is_symlink() or not (fixture/'result.bin').exists(): continue
    with tempfile.TemporaryDirectory() as directory:
        run=Path(directory)
        shutil.copy2(fixture/'fixture.bin',run/'fixture.bin')
        subprocess.run([str(exe.resolve())],cwd=run,check=True,capture_output=True,
                       env={**os.environ,'OMP_NUM_THREADS':'2'})
        a=np.fromfile(fixture/'result.bin',np.complex128)
        b=np.fromfile(run/'result.bin',np.complex128)
        np.testing.assert_allclose(a,b,atol=5e-13,rtol=5e-13)
        rows.append([fixture.name,a.size,float(np.max(np.abs(a-b))),
                     float(np.linalg.norm(a-b)/np.linalg.norm(b))])
assert len(rows)==4, 'Expected four nonperiodic fixtures (two overlaps, two thread counts)'
with output.open('w') as f:
    writer=csv.writer(f,lineterminator="\n")
    writer.writerow(['fixture','complex_elements','max_abs_error','relative_l2_error'])
    writer.writerows(rows)
print(rows)
