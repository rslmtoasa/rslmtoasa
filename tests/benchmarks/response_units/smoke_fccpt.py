#!/usr/bin/env python3
"""Frozen-potential output checks. Usage: smoke_fccpt.py EXECUTABLE NEW_DIRECTORY.
Requires numpy and f90nml. No SCF; layer case is only a geometry/units smoke test.
"""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import f90nml
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
exe=Path(sys.argv[1]).resolve(); dest=Path(sys.argv[2]).resolve()
dest.mkdir(parents=True,exist_ok=False)
rows=[]
for name,operator,geometry,n3 in [('charge_xx','charge','bulk',12),
                                ('spin_hall','spin','bulk',12),
                                ('orbital_accumulation','orbital_accumulation','bulk',12),
                                ('charge_sheet','charge','sheet',1)]:
    run=dest/name; run.mkdir()
    source=ROOT/'example/conductivity/fccPt/sigma_xy_z-spin'
    for p in source.glob('*.nml'): shutil.copy2(p,run/p.name)
    patch={'control':{'cond_ll':50,'linear_out':operator},
           'lattice':{'n1':12,'n2':12,'n3':n3,'b1':True,'b2':True,'b3':n3>1},
           'hamiltonian':{'v_alpha':[1,0,0] if operator=='charge' else [0,1,0]}}
    f90nml.patch(run/'input.nml',patch,run/'patched.nml'); (run/'patched.nml').replace(run/'input.nml')
    with (run/'input.nml').open('a') as f: f.write(f"\n&conductivity\nresponse_geometry='{geometry}'\n/\n")
    start=time.monotonic()
    with (run/'run.log').open('w') as log:
        p=subprocess.run([str(exe)],cwd=run,stdout=log,stderr=subprocess.STDOUT,
            env={**os.environ,'OMP_NUM_THREADS':'4','OMP_DYNAMIC':'FALSE','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OMP_STACKSIZE':'64M'})
    assert p.returncode==0,(name,p.returncode)
    elapsed=time.monotonic()-start
    raw=np.loadtxt(run/'cond_total.out'); physical=np.loadtxt(run/'response_total.out')
    assert np.isfinite(raw).all() and np.isfinite(physical).all()
    text=(run/'response_units.txt').read_text()
    weight=float(next(line.split('=')[1] for line in text.splitlines() if line.startswith('trace weights')).split()[0])
    conversion=float(next(line.split('=')[1] for line in text.splitlines() if line.startswith('output conversion')))
    expected=raw[:,1:]*weight*conversion
    error=np.max(np.abs(expected-physical[:,1:]))
    scale=np.max(np.abs(physical[:,1:]))
    assert error<2e-6*max(scale,1e-100),(name,error,scale)
    orbitals=np.loadtxt(run/'response_orbitals_real.out')
    np.testing.assert_allclose(orbitals[:,1:].sum(axis=1),physical[:,1],rtol=2e-13,atol=1e-12)
    typed=np.loadtxt(run/'response_type_1.out')
    np.testing.assert_allclose(typed[:,0],physical[:,0],rtol=2e-14)
    # Different summation orders cancel near-zero components differently.
    assert np.max(np.abs(typed[:,1:]-physical[:,1:])) < 1e-12*max(scale,1e-100)
    unit=next(line.split('=',1)[1].strip() for line in text.splitlines() if line.startswith('units ='))
    rows.append([name,geometry,unit,len(raw),elapsed,error,error/max(scale,1e-100)])
    print(rows[-1],flush=True)
with (dest/'smoke_results.csv').open('w') as f:
    writer=csv.writer(f,lineterminator='\n')
    writer.writerow(['case','geometry','unit','energy_points','seconds','max_conversion_error','scaled_conversion_error'])
    writer.writerows(rows)
