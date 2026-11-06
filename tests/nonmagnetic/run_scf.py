#!/usr/bin/env python3
"""Fe-Co B2: constrain Co, retain magnetic Fe; check absent == all-false."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import f90nml

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--workdir', type=Path, required=True)
p.add_argument('--nsp', type=int, default=4)
p.add_argument('--steps', type=int, default=4)
p.add_argument('--mpi-exec')
p.add_argument('--ranks', type=int, default=2)
p.add_argument('--baseline', type=Path)
p.add_argument('--mixtype', choices=('linear', 'broyden'), default='linear')
p.add_argument('--cold', action='store_true')
a = p.parse_args()
source = Path(__file__).resolve().parents[1] / 'scf/cases/impurity/B2FeCo'
results = {}


def flatten(x):
    if isinstance(x, (list, tuple)):
        return [v for row in x for v in flatten(row)]
    return [x]


def run(name, flags, binary):
    # Fresh paths prevent any generated checkpoint from affecting repeat runs.
    a.workdir.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix=name + '-', dir=a.workdir.resolve()))
    for label in ('Fe', 'Co'):
        shutil.copy2(source / (label + '.nml'), directory)
        if a.nsp == 4:
            # Preserve the zero-based orbital indices/colon slices in the fixture.
            original = directory / (label + '.nml')
            patched = directory / (label + '.tmp')
            f90nml.patch(original, {'par': {'mom': [0.6, 0.0, 0.8]}}, patched)
            patched.replace(original)
    n = f90nml.read(source / 'input.nml')
    n['calculation']['pre_processing'] = 'bravais'
    n['lattice']['ntype'] = 2
    n['lattice']['rc'] = 12
    n['lattice']['ct'] = [3.0, 3.0]
    for key in ('nclu', 'inclu'):
        del n['lattice'][key]
    n['atoms']['label'] = ['Fe', 'Co']
    n['self']['nstep'] = a.steps
    n['self']['conv_thr'] = 1e-12
    if flags is not None:
        n['self']['force_nonmagnetic'] = flags
    n['control']['calctype'] = 'B'
    n['control']['nsp'] = a.nsp
    n['control']['lld'] = 10
    n['energy']['channels_ldos'] = 800
    n['mix']['beta'] = 0.05
    n['mix']['mixtype'] = a.mixtype
    n['self']['cold'] = a.cold
    n['hamiltonian']['hoh'] = True
    n.write(directory / 'input.nml', force=True)
    command = [str(binary.resolve())]
    if a.mpi_exec:
        command = [a.mpi_exec, '-np', str(a.ranks)] + command
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    with (directory / 'run.log').open('w') as out:
        subprocess.run(command, cwd=directory, env=env, stdout=out,
                       stderr=subprocess.STDOUT, check=True, timeout=240)
    potentials = {label: dict(f90nml.read(directory / (label + '_out.nml'))['par'])
                  for label in ('Fe', 'Co')}
    for par in potentials.values():
        assert all(math.isfinite(v) for v in flatten(list(par.values())) if isinstance(v, (int, float)))
    assert 'Error while reading namelist' not in (directory / 'run.log').read_text()
    report = (directory / 'report.out').read_text()
    moments = [float(v) for v in re.findall(r'Spin moment of atom\s+\d+:\s*([\d.Ee+-]+)', report)]
    assert len(moments) == 2
    assert moments[0] > 0.1, f'Fe lost magnetism: {moments}'
    if flags and any(flags):
        par = potentials['Co']
        for key in ('ql', 'pl', 'center_band', 'width_band', 'shifted_band', 'obar'):
            # f90nml stores Fortran arrays with the last index outermost.
            assert flatten(par[key][0]) == flatten(par[key][1]), (key, par[key])
        line = re.search(r'Nonmagnetic check atom\s+2:(.*)', report)
        assert line, 'Missing runtime build_pot diagnostics'
        checks = [float(v) for v in line.group(1).split()]
        assert len(checks) == 6 and max(map(abs, checks)) < 1e-14, checks
        assert moments[1] == 0.0
        assert any(x != y for x, y in zip(potentials['Fe']['center_band'][0], potentials['Fe']['center_band'][1]))
        for filename in ('minfo.out', 'linfo.out', 'angles_magmom.out', 'angles_lmom.out'):
            text = (directory / filename).read_text().lower()
            assert not re.search(r'\b(nan|inf)\b', text), filename
        results[name] = dict(zip(('Q_up-Q_down', 'cx1', 'wx1', 'cex1', 'obx1', 'spin_moment'), checks))
        results[name]['Fe_spin_moment'] = moments[0]
    else:
        assert moments[1] > 0.1, f'Co should be magnetic without constraint: {moments}'
        results[name] = {'Fe_spin_moment': moments[0], 'Co_spin_moment': moments[1]}
    results[name]['directory'] = directory.name
    return potentials


absent = run('absent', None, a.binary)
all_false = run('all_false', [False, False], a.binary)
assert absent == all_false, 'All-false changed default behavior'
if a.baseline:
    baseline = run('baseline', None, a.baseline)
    max_difference = 0.0
    for label in absent:
        for key in absent[label]:
            for x, y in zip(flatten(absent[label][key]), flatten(baseline[label][key])):
                if isinstance(x, (int, float)):
                    max_difference = max(max_difference, abs(x-y))
                    assert math.isclose(x, y, abs_tol=1e-7, rel_tol=1e-7), (label, key, x, y)
                else:
                    assert x == y, (label, key)
    results['baseline_max_abs_difference'] = max_difference
run('constrained', [False, True], a.binary)
(a.workdir / 'checks.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
