#!/usr/bin/env python3
"""LR-METHOD-05R accepted Fe evidence and projected bare selector regression.

A: raw material Ward action and SR/Pauli discrepancy; no fitted tolerance.
C: executed d/spd/both selector semantics. Artifacts retain finite-eta labels.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / 'tests/integration/tddft_driver_smoke'
parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True, type=Path)
parser.add_argument('--scratch', required=True, type=Path)
parser.add_argument('--selectors', action='store_true')
args = parser.parse_args()
env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
text = (CASE / 'input_compact_alsda.nml').read_text()
text = re.sub(r"database\s*=\s*'([^']+)'", lambda m: "database = '" + str((CASE / m[1]).resolve()) + "/'", text)

def run(name, response):
    folder = args.scratch.resolve() / name
    folder.mkdir(parents=True, exist_ok=True)
    # Keep the accepted gate's lattice, SCF, reciprocal mesh and temperature.
    deck = text[:text.index('&linear_response')] + response
    (folder / 'input.nml').write_text(deck)
    with (folder / 'run.log').open('w') as log:
        subprocess.run([str(args.binary.resolve()), 'input.nml'], cwd=folder, env=env,
                       stdout=log, stderr=subprocess.STDOUT, timeout=720, check=True)
    log = (folder / 'run.log').read_text()
    if 'Converged!' not in log:
        raise RuntimeError('accepted Fe SCF did not converge')
    return folder

folder = run('material', """&linear_response
formulation='tddft', representation='product_compact', bare_response='lehmann'
interaction='alsda', diagnostics='invariants', channel='chi_plus'
n_q=3, q_list=0.0,0.0,0.0, 0.25,0.0,0.0, -0.25,0.0,0.0
n_omega=1, use_omega_grid=.false., omega_min=0.0, omega_max=0.0
eta=0.02, n_eta=1, response_lmax=-1, write_full_matrix=.false.
output_file='alsda_05r.dat'
/
""")
rows = (folder / 'alsda_05r.dat').read_text().splitlines()
metadata = {}
ward = []
for line in rows:
    if line.startswith('# ') and ' = ' in line:
        key, value = line[2:].split(' = ', 1)
        if key == 'raw_Ward_eta_abs_rel_inf_overlap':
            ward.append([float(x) for x in value.split()])
        else:
            metadata[key] = value
static_dyson = [[float(x) for x in line.split()] for line in rows
                if line.strip() and not line.startswith('#') and len(line.split()) == 8]
scf_metadata = {}
for line in (folder / 'kspace_scf_state.dat').read_text().splitlines():
    if line.startswith('# ') and ' = ' in line:
        key, value = line[2:].split(' = ', 1)
        scf_metadata[key] = value.strip()
report = {'branch': subprocess.check_output(['git','branch','--show-current'], cwd=ROOT, text=True).strip(),
          'commit': subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
          'working_tree_changes': True, 'evidence_category': 'A',
          'metadata': metadata, 'accepted_scf_metadata': scf_metadata,
          'finite_eta_Ward_eta_abs_rel_inf_overlap': ward,
          'static_Dyson_eta_minSV_maxSV_cond_minEig_abs_rel_inf': static_dyson}
if len(ward) != 5:
    raise RuntimeError('missing raw independent Ward eta ladder')
(args.scratch / 'material_record.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
if args.selectors:
    for selector, expected in [('d', {'d'}), ('spd', {'spd'}), ('both', {'d','spd'})]:
        folder = run('selector_'+selector, f"""&linear_response
formulation='projected', representation='radial_points', bare_response='lehmann'
interaction='none', projection='{selector}', diagnostics='none', channel='chi_plus'
n_q=1, q_list=0.0,0.0,0.0, n_omega=1, use_omega_grid=.false.
omega_min=0.0, omega_max=0.0, eta=0.02, n_eta=1, response_lmax=4
output_file='selector.dat'
/
""")
        selected = {line.split(' = ',1)[1].strip() for line in (folder/'selector.dat').read_text().splitlines()
                    if line.startswith('# projection = ')}
        if selected != expected:
            raise RuntimeError(f'{selector}: got {selected}, expected {expected}')
        print(f'projected bare {selector}: PASS ({sorted(selected)})')
