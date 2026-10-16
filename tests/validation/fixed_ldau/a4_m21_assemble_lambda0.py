#!/usr/bin/env python3
"""Hash-first assembly of raw W1 lambda-zero captures for the frozen analyzer.

The inventory phase binds source execution manifests and raw CSV/stdout hashes
without reading numerical contents. The assemble phase refuses unless the
inventory hash is explicitly supplied, then preserves exact CSV/log tokens.
It never launches the production executable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from decimal import Decimal, localcontext, ROUND_CEILING
from collections import defaultdict
from pathlib import Path

ROOT = Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
EXEC_MANIFEST = Path('/private/tmp/a4-m21-ks-lambda0-execution/execution_manifest.json')
EXPECTED_EXEC_SHA = 'cfea7226f0e2bcbfae778d7a45a92205e6ff6dc70a01d66c2935a8f5f5bbb461'
EXPECTED_BINARY_SHA = '536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4'
EXPECTED_SOURCE_SHA = '2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010'
CALIBRATION_PREREG = ROOT / 'docs/validation/fixed_ldau/evidence/a4_m21_ks_lambda0_calibration_preregistration.json'
EXPECTED_CALIBRATION_PREREG_SHA = '50f24c646af2e7c1eedc050d6b0b3c14b7cb55b5d6134e02cfa081913e1ba80f'
BASE_RUN_MANIFEST = Path('/private/tmp/rslmto-localmom-m02-release/ks/run_manifest.json')
EXPECTED_BASE_RUN_SHA = '8be28890f4926cf8cab889b7995bd28ba3b10a85523db2a978ea1ca23ce7502a'
BASE_STR_OUT = Path('/private/tmp/rslmto-localmom-m02-release/ks/str.out')
EXPECTED_BASE_STR_SHA = '36cdc02b7182a27f94568858df248c21a5ef076045edb81d0967262b97030fc3'
ANALYZER = ROOT / 'tests/validation/fixed_ldau/a4_m21_w1_richardson.py'
KIND_ROWS = ('eigenvalue', 'occupation', 'eigenvector', 'k_weight', 'basis_map')
FIELDS = ('kind','run_id','lambda','iteration','stage','site','type','l','spin','k','band',
          'row','col','order','weight','energy','scalar','re','im')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def token(row: dict[str, str], field: str) -> str:
    value = row[field].strip()
    if not re.fullmatch(r'[+-]?(?:\d+\.?\d*|\.\d+)[EeDd][+-]?\d+', value):
        raise ValueError(f'{field} is not an ES24.16 numeric token: {value!r}')
    return value.replace('D', 'E').replace('d', 'E')


def read_inventory_source() -> tuple[dict, str]:
    if sha(EXEC_MANIFEST) != EXPECTED_EXEC_SHA:
        raise ValueError('lambda-zero execution manifest hash mismatch')
    if sha(ROOT/'source/self.f90') != EXPECTED_SOURCE_SHA or sha(Path('/private/tmp/a4-m2-serial/bin/rslmto.x')) != EXPECTED_BINARY_SHA:
        raise ValueError('frozen production source or binary hash mismatch')
    if sha(CALIBRATION_PREREG) != EXPECTED_CALIBRATION_PREREG_SHA:
        raise ValueError('lambda-zero calibration preregistration hash mismatch')
    batch = json.loads(EXEC_MANIFEST.read_text())
    if batch.get('status') != 'ALL_PROCESSES_EXITED' or len(batch.get('cases', {})) != 8:
        raise ValueError('expected the reviewed eight-case lambda-zero operational batch')
    return batch, sha(EXEC_MANIFEST)


def reciprocal_mesh_from_pre(path: Path, expected_sha: str) -> tuple[int, list[int]]:
    """Derive the exact full-BZ k domain from the hash-bound PRE reciprocal group."""
    if not path.is_file() or sha(path) != expected_sha:
        raise ValueError('PRE input deck absent or hash differs from execution binding')
    text=path.read_text(errors='strict')
    match=re.search(r'(?ims)^\s*&reciprocal\s*(.*?)^\s*/',text)
    if not match: raise ValueError('PRE deck has no reciprocal namelist')
    group=match.group(1)
    def assignment(name: str) -> str:
        found=re.findall(rf'(?im)^\s*{name}\s*=\s*([^!\n,]+)',group)
        if len(found)!=1: raise ValueError(f'PRE reciprocal group must define {name} exactly once')
        return found[0].strip().lower()
    dims=[int(assignment(name)) for name in ('nk1','nk2','nk3')]
    reduction=assignment('use_symmetry_reduction').split()[0]
    if reduction not in ('.false.','false'):
        raise ValueError('expected full unreduced Monkhorst-Pack mesh for the locked calibration')
    nk=math.prod(dims)
    if nk < 1: raise ValueError('invalid PRE mesh dimensions')
    return nk, list(range(nk))


def site_count_from_frozen_base(prereg: dict) -> tuple[int, dict]:
    """Bind site cardinality to the accepted base run's structure output, not CSV rows."""
    if sha(BASE_RUN_MANIFEST) != EXPECTED_BASE_RUN_SHA or sha(BASE_STR_OUT) != EXPECTED_BASE_STR_SHA:
        raise ValueError('accepted base run manifest/structure output hash mismatch')
    base=json.loads(BASE_RUN_MANIFEST.read_text())
    common=prereg['common_state']
    if (base.get('run_id')!='m02_ks_u0' or base.get('status')!='COMPLETE_CONVERGED' or
        base.get('Fe_out_sha256')!=common['base_seed_sha256'] or base.get('input_sha256')!=common['base_deck_sha256']):
        raise ValueError('accepted base run identity does not match lambda-zero common seed/deck')
    lines=BASE_STR_OUT.read_text(errors='strict').splitlines()
    try: start=next(i for i,line in enumerate(lines) if line.strip()=='LATTICE COORDINATES')+1
    except StopIteration as exc: raise ValueError('base str.out lacks LATTICE COORDINATES') from exc
    site_lines=[]
    for line in lines[start:]:
        if re.fullmatch(r'\s*\d+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*',line): site_lines.append(line)
        elif site_lines: break
    ids=[int(line.split()[0]) for line in site_lines]
    if not ids or ids!=list(range(1,len(ids)+1)):
        raise ValueError('base structure site list is empty or noncontiguous')
    return len(ids), {'base_run_manifest_path':str(BASE_RUN_MANIFEST),'base_run_manifest_sha256':EXPECTED_BASE_RUN_SHA,
      'base_structure_path':str(BASE_STR_OUT),'base_structure_sha256':EXPECTED_BASE_STR_SHA,
      'base_seed_sha256':common['base_seed_sha256'],'base_deck_sha256':common['base_deck_sha256'],
      'site_count_source':'one-based contiguous LATTICE COORDINATES rows from hash-bound accepted M0.2 base str.out'}


def inventory() -> dict:
    batch, batch_sha = read_inventory_source()
    prereg=json.loads(CALIBRATION_PREREG.read_text())
    expected_site_count, site_provenance=site_count_from_frozen_base(prereg)
    prereg_cases={case['run_id']:case for case in prereg['cases']}
    cases = {}
    for run_id in sorted(k for k in batch['cases'] if k.startswith('w1_')):
        entry = batch['cases'][run_id]
        manifest_path = Path(entry['run_manifest'])
        run = json.loads(manifest_path.read_text())
        if entry.get('status') != 'PROCESS_EXITED' or entry.get('returncode') != 0 or run.get('returncode') != 0:
            raise ValueError(f'{run_id}: process status is not successful')
        if run.get('run_id') != run_id or run.get('status') != 'PROCESS_EXITED':
            raise ValueError(f'{run_id}: run manifest identity/status mismatch')
        bindings = run.get('bindings', {})
        if bindings.get('binary_sha256') != EXPECTED_BINARY_SHA:
            raise ValueError(f'{run_id}: binary binding mismatch')
        if run_id not in prereg_cases: raise ValueError(f'{run_id}: absent from frozen calibration preregistration')
        registered=prereg_cases[run_id]
        pre_path=Path(registered['input_nml_pre_path'])
        if bindings.get('input_nml_pre_sha256') != registered['input_nml_pre_sha256']:
            raise ValueError(f'{run_id}: execution PRE hash differs from frozen preregistration')
        if bindings.get('Fe_nml_pre_sha256') != registered['Fe_nml_pre_sha256']:
            raise ValueError(f'{run_id}: execution seed hash differs from frozen preregistration')
        expected_nk, expected_k_indices=reciprocal_mesh_from_pre(pre_path, registered['input_nml_pre_sha256'])
        post = run.get('post_files', {})
        for name in ('records.csv', 'stdout.log'):
            path = manifest_path.parent/name
            expected = post.get(name, {}).get('sha256')
            if not expected or not path.is_file() or sha(path) != expected:
                raise ValueError(f'{run_id}: {name} absent or hash differs from run manifest')
        cases[run_id] = {
            'run_id':run_id,
            'mode': 'local' if '_local_' in run_id else 'spectral',
            'threshold': 'tight' if run_id.endswith('_tight') else 'standard',
            'run_manifest_path': str(manifest_path),
            'run_manifest_sha256': sha(manifest_path),
            'records_csv_path': str(manifest_path.parent/'records.csv'),
            'records_csv_sha256': post['records.csv']['sha256'],
            'records_csv_size': post['records.csv']['size_bytes'],
            'stdout_path': str(manifest_path.parent/'stdout.log'),
            'stdout_sha256': post['stdout.log']['sha256'],
            'stdout_size': post['stdout.log']['size_bytes'],
            'input_nml_pre_sha256': bindings['input_nml_pre_sha256'],
            'input_nml_pre_path': str(pre_path),
            'expected_n_kpoints': expected_nk,
            'expected_k_indices': expected_k_indices,
            'expected_site_count': expected_site_count,
            'site_provenance': site_provenance,
            'basis_domain_rationale': f"Accepted base structure independently records {expected_site_count} site(s); source basis map is site-major with two spin blocks and 9 spd rows/site/spin.",
            'seed_sha256': bindings['Fe_nml_pre_sha256'],
            'source_sha256': EXPECTED_SOURCE_SHA,
            'binary_sha256': EXPECTED_BINARY_SHA,
        }
    if set(cases) != {'w1_spectral_standard','w1_spectral_tight','w1_local_standard','w1_local_tight'}:
        raise ValueError('W1 anchor case set is incomplete or unexpected')
    prereg=json.loads(CALIBRATION_PREREG.read_text())
    site_count, site_provenance=site_count_from_frozen_base(prereg)
    return {'schema':'rslmto-a4-m21-lambda0-raw-inventory-v1',
            'status':'HASH_INVENTORY_ONLY_NO_NUMERIC_CONTENT_READ',
            'execution_manifest_path':str(EXEC_MANIFEST), 'execution_manifest_sha256':batch_sha,
            'calibration_preregistration_path':str(CALIBRATION_PREREG),
            'calibration_preregistration_sha256':EXPECTED_CALIBRATION_PREREG_SHA,
            'expected_site_count':site_count,'site_provenance':site_provenance,
            'source_sha256':EXPECTED_SOURCE_SHA, 'binary_sha256':EXPECTED_BINARY_SHA,
            'cases':cases}


def parse_total_log(text: str) -> dict[str, str]:
    lines=text.splitlines()
    converged=[i for i,line in enumerate(lines) if 'Converged!' in line]
    diagnostic_indices=[i for i,line in enumerate(lines) if 'HUBBARD_MOMENT_ENERGY_TOTAL ' in line]
    if len(converged)!=1 or not diagnostic_indices or diagnostic_indices[-1]>=converged[0]:
        raise ValueError('logger must contain a unique final convergence marker after the last total-energy diagnostic')
    final_line=lines[diagnostic_indices[-1]]
    fields=dict(re.findall(r'([A-Za-z0-9_]+)=\s*([^\s]+)',final_line.split('HUBBARD_MOMENT_ENERGY_TOTAL ',1)[1]))
    required = ('physical_total_energy','TS','candidate_status')
    if any(k not in fields for k in required):
        raise ValueError('total diagnostic line misses required physical energy / TS / status fields')
    if 'candidate_status=PROVISIONAL_ONLY;' not in final_line:
        raise ValueError('unexpected candidate status or logger line')
    fields['_diagnostic_count']=str(len(diagnostic_indices))
    fields['_selected_line_index']=str(diagnostic_indices[-1])
    return fields


def assemble_case(meta: dict) -> dict:
    csv_path, log_path = Path(meta['records_csv_path']), Path(meta['stdout_path'])
    if sha(csv_path) != meta['records_csv_sha256'] or sha(log_path) != meta['stdout_sha256']:
        raise ValueError('raw capture changed after hash inventory')
    manifest_path=Path(meta['run_manifest_path'])
    if sha(manifest_path)!=meta['run_manifest_sha256']:
        raise ValueError('run manifest changed after hash inventory')
    run = json.loads(manifest_path.read_text())
    if run.get('run_id') != meta['run_id']:
        raise ValueError('run manifest identity differs from inventory case')
    if run.get('bindings',{}).get('binary_sha256') != meta['binary_sha256'] or run.get('bindings',{}).get('input_nml_pre_sha256') != meta['input_nml_pre_sha256']:
        raise ValueError('run manifest source/input binding differs from inventory')
    if 'input_nml_pre_path' in meta:
        nk, kdomain=reciprocal_mesh_from_pre(Path(meta['input_nml_pre_path']), meta['input_nml_pre_sha256'])
        if nk != meta['expected_n_kpoints'] or kdomain != meta['expected_k_indices']:
            raise ValueError('inventory k domain differs from hash-bound PRE mesh')
    env = run['environment']
    run_id = run['run_id']
    if env.get('RSLMTO_FIXED_U_TRACE_ITERATION') != '0' or env.get('RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY') != '1':
        raise ValueError(f'{run_id}: trace is not a capture-only final-state trace')

    basis_rows_by_iteration=defaultdict(list); scalar_rows_by_iteration=defaultdict(lambda:defaultdict(list))
    iterations=set(); row_count=0
    with csv_path.open(newline='', encoding='utf-8') as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != FIELDS: raise ValueError('trace CSV column schema mismatch')
        for r in reader:
            row_count += 1
            if r['run_id'] != run_id or float(r['lambda']) != float(env['RSLMTO_FIXED_U_LAMBDA']):
                raise ValueError(f'{run_id}: trace run/lambda identity mismatch')
            iteration=int(r['iteration']); iterations.add(iteration)
            if r['stage']=='EIGENPAIRS' and r['kind']=='basis_map': basis_rows_by_iteration[iteration].append(r)
            if r['stage']=='EIGENPAIRS' and r['kind'] in ('fermi_level','temperature_k','temperature_ry'):
                scalar_rows_by_iteration[iteration][r['kind']].append(r)
    if row_count == 0: raise ValueError(f'{run_id}: empty raw CSV')
    final_iteration=max(iterations)
    if iterations not in ({final_iteration},set(range(1,final_iteration+1))):
        raise ValueError(f'{run_id}: CSV iteration history has gaps or does not reach a complete final iteration')
    basis_rows=basis_rows_by_iteration[final_iteration]
    basis = {}
    for r in basis_rows:
        idx = int(r['row']); key = (int(r['site']),int(r['spin']), int(r['l']), int(r['col']))
        if idx in basis or idx < 0:
            raise ValueError(f'{run_id}: duplicate/out-of-range basis row')
        basis[idx] = key
    sites=sorted({key[0] for key in basis.values()})
    if sites != list(range(meta['expected_site_count'])):
        raise ValueError(f'{run_id}: site map differs from independently preregistered site domain')
    basis_size=18*len(sites)
    if set(basis) != set(range(basis_size)):
        raise ValueError(f'{run_id}: basis map does not cover the complete scalar-collinear spd basis')
    drows = {}
    for idx, (site,spin,l,m) in basis.items():
        if l == 2:
            key=(site,spin,m)
            if key in drows: raise ValueError(f'{run_id}: duplicate site/spin/d-m basis mapping {key}')
            drows[key] = idx
    for idx,(site,spin,l,m) in basis.items():
        expected=site*18+spin*9+l*l+(m+l)
        if spin not in (0,1) or l not in (0,1,2) or not (-l<=m<=l) or idx!=expected:
            raise ValueError(f'{run_id}: basis_map row does not match source site/spin/l/m ordering')
    if set(drows) != {(site,s,m) for site in sites for s in (0,1) for m in range(-2,3)}:
        raise ValueError(f'{run_id}: d-shell spin/orbital map incomplete')

    weights: dict[int, dict[str,str]] = {}
    eig: dict[tuple[int,int], str] = {}
    occ: dict[tuple[int,int], str] = {}
    vectors: dict[tuple[int,int,int], tuple[str,str]] = {}
    vector_keys_all:set[tuple[int,int,int]]=set()
    drow_set=set(drows.values())
    with csv_path.open(newline='', encoding='utf-8') as f:
        reader=csv.DictReader(f)
        for r in reader:
            if r['stage']!='EIGENPAIRS' or int(r['iteration'])!=final_iteration: continue
            kind=r['kind']
            if kind=='k_weight':
                k=int(r['k'])
                if k in weights: raise ValueError(f'{run_id}: duplicate k weight')
                weights[k]={'weight_text':token(r,'weight')}
            elif kind in ('eigenvalue','occupation'):
                k,b=int(r['k']),int(r['band']); key=(k,b)
                target=eig if kind=='eigenvalue' else occ
                if key in target: raise ValueError(f'{run_id}: duplicate {kind} row {key}')
                target[key]=token(r,'energy' if kind=='eigenvalue' else 'scalar')
            elif kind=='eigenvector':
                k,b,i=int(r['k']),int(r['band']),int(r['row'])
                full_key=(k,b,i)
                if full_key in vector_keys_all: raise ValueError(f'{run_id}: duplicate eigenvector row {full_key}')
                vector_keys_all.add(full_key)
                if i in drow_set:
                    key=(k,b,i)
                    vectors[key]=(token(r,'re'),token(r,'im'))
    n_k = meta['expected_n_kpoints']
    n_band = 18*meta['expected_site_count']
    if set(b for _,b in eig) != set(range(n_band)):
        raise ValueError(f'{run_id}: observed global bands differ from source/PRE-derived domain')
    if set(weights) != set(meta['expected_k_indices']) or set(eig) != {(k,b) for k in meta['expected_k_indices'] for b in range(n_band)} or set(occ) != set(eig):
        raise ValueError(f'{run_id}: incomplete zero-based k/band eigenvalue and occupation domains')
    if n_band != basis_size:
        raise ValueError(f'{run_id}: global spin-inclusive band count {n_band} != basis dimension {basis_size}')
    expected_vec_keys = {(k,b,i) for k in range(n_k) for b in range(n_band) for i in drows.values()}
    if set(vectors) != expected_vec_keys:
        raise ValueError(f'{run_id}: missing/extra d-shell eigenvector components')
    expected_all_vector_keys={(k,b,i) for k in meta['expected_k_indices'] for b in range(n_band) for i in range(basis_size)}
    if vector_keys_all != expected_all_vector_keys:
        raise ValueError(f'{run_id}: complete eigenvector row domain differs from PRE/source basis dimensions')

    # Build separate spin P matrices using the source projector P_ab=psi_a psi*_b.
    pspin = {s:[[0j for _ in range(5)] for _ in range(5)] for s in (0,1)}
    bspin = {s:[[0.0 for _ in range(5)] for _ in range(5)] for s in (0,1)}
    for s in (0,1):
        for i,m in enumerate(range(-2,3)):
            for j,n in enumerate(range(-2,3)):
                site_values=[]; site_bounds=[]
                for site in sites:
                    ia=drows[(site,s,m)]; ib=drows[(site,s,n)]
                    terms=[]; input_bound=0.0
                    for k in range(n_k):
                        wt_text=weights[k]['weight_text']; wt=float(wt_text); qw=token_radius(wt_text)
                        for b in range(n_band):
                            ft=occ[(k,b)]; f=float(ft); qf=token_radius(ft)
                            ar,ai=vectors[(k,b,ia)]; br,bi=vectors[(k,b,ib)]
                            a=complex(float(ar),float(ai)); bb=complex(float(br),float(bi))
                            qa=math.nextafter(math.hypot(token_radius(ar),token_radius(ai)),math.inf)
                            qb=math.nextafter(math.hypot(token_radius(br),token_radius(bi)),math.inf)
                            term=wt*f*a*bb.conjugate(); terms.append(term)
                            term_bound=complex_product_token_error_bound(wt,f,a,bb,qw,qf,qa,qb)
                            input_bound=up_add(input_bound,term_bound)
                    val=sum(terms,0j)
                    abs_scale=math.nextafter(math.fsum(abs(z) for z in terms),math.inf)
                    prod_bound=math.nextafter(gamma_bound(16)*abs_scale,math.inf)
                    sum_bound=math.nextafter(gamma_bound(2*len(terms))*abs_scale,math.inf)
                    site_values.append(val); site_bounds.append(up_add(up_add(input_bound,prod_bound),sum_bound))
                pspin[s][i][j]=sum(site_values,0j)
                site_sum_bound=math.nextafter(math.fsum(site_bounds),math.inf)
                site_scale=math.nextafter(math.fsum(abs(z) for z in site_values),math.inf)
                site_sum_fp=math.nextafter(gamma_bound(2*len(site_values))*site_scale,math.inf)
                bspin[s][i][j]=up_add(site_sum_bound,site_sum_fp)
    ptotal=[[pspin[0][i][j]+pspin[1][i][j] for j in range(5)] for i in range(5)]
    btotal=[[up_add(up_add(bspin[0][i][j],bspin[1][i][j]),
             up_mul(gamma_bound(2),up_add(abs(pspin[0][i][j]),abs(pspin[1][i][j])))) for j in range(5)] for i in range(5)]
    states=[]
    for k in range(n_k):
        wt=weights[k]['weight_text']
        for b in range(n_band):
            states.append({'k_index':k,'band_index':b,'eigenvalue_text':eig[(k,b)],
                           'occupation_text':occ[(k,b)],'weight_text':wt})
    scalar_rows=scalar_rows_by_iteration[final_iteration]['fermi_level']
    tempk_rows=scalar_rows_by_iteration[final_iteration]['temperature_k']
    tempry_rows=scalar_rows_by_iteration[final_iteration]['temperature_ry']
    for name,items in (('fermi_level',scalar_rows),('temperature_k',tempk_rows),('temperature_ry',tempry_rows)):
        if len(items)!=1: raise ValueError(f'{run_id}: expected exactly one {name} scalar trace row')
    log_fields=parse_total_log(log_path.read_text(errors='strict'))
    if int(log_fields['_diagnostic_count']) != final_iteration:
        raise ValueError(f'{run_id}: final CSV iteration does not match SCF diagnostics before convergence')
    if meta['mode']=='spectral' and log_fields.get('candidate_a')=='UNAVAILABLE':
        raise ValueError(f'{run_id}: spectral mode lacks candidate_a')
    if meta['mode']=='local' and log_fields.get('candidate_b')=='UNAVAILABLE':
        raise ValueError(f'{run_id}: local mode lacks candidate_b')
    return {'mode':meta['mode'],'threshold':meta['threshold'],'lambda_ry':0.0,
       'lambda_runtime_text':env['RSLMTO_FIXED_U_LAMBDA'],
       'run_manifest_path':meta['run_manifest_path'],'run_manifest_sha256':meta['run_manifest_sha256'],
       'records_csv_path':str(csv_path),'records_csv_sha256':meta['records_csv_sha256'],
       'stdout_path':str(log_path),'stdout_sha256':meta['stdout_sha256'],
       'source_sha256':meta['source_sha256'],'binary_sha256':meta['binary_sha256'],
       'input_nml_pre_path':meta['input_nml_pre_path'],'input_nml_pre_sha256':meta['input_nml_pre_sha256'],
       'expected_n_kpoints':n_k,'expected_k_indices':meta['expected_k_indices'],
       'expected_site_count':meta['expected_site_count'],'site_provenance':meta.get('site_provenance'),
       'expected_band_indices':list(range(n_band)),
       'final_iteration':final_iteration,'captured_iteration_domains':sorted(iterations),
       'logger_selection_rule':'last HUBBARD_MOMENT_ENERGY_TOTAL before unique final Converged!; source emits one diagnostic per SCF loop iteration',
       'band_index_semantics':'zero-based global spin-inclusive row of eigenvalues(:,ik)',
       'P_d':[[[ptotal[i][j].real,ptotal[i][j].imag] for j in range(5)] for i in range(5)],
       'P_d_component_error_bound':btotal,
       'P_d_by_spin':{str(s):[[[pspin[s][i][j].real,pspin[s][i][j].imag] for j in range(5)] for i in range(5)] for s in (0,1)},
       'entropy_inputs':{'expected_n_kpoints':n_k,'expected_band_indices':list(range(n_band)),
        'band_index_semantics':'zero-based global spin-inclusive row of eigenvalues(:,ik)',
        'temperature_k_text':token(tempk_rows[0],'scalar'),'temperature_ry_text':token(tempry_rows[0],'scalar'),
        'fermi_level_text':token(scalar_rows[0],'scalar'),'states':states},
       'ETOT_spectral_text':log_fields['physical_total_energy'],
       'ETOT_local_text':log_fields['physical_total_energy'],'TS_text':log_fields['TS'],
       'source_log_fields':{k:log_fields[k] for k in ('physical_total_energy','TS','entropy','candidate_a','candidate_b','test_F_local','test_F_spectral') if k in log_fields}}


def half_quantum(text: str) -> float:
    x=text.strip().upper().replace('D','E'); m,e=x.split('E'); ndec=len(m.partition('.')[2])
    with localcontext() as ctx:
        ctx.prec=80; ctx.rounding=ROUND_CEILING
        exact=Decimal(5).scaleb(int(e)-ndec-1)
        value=float(exact)
    return math.nextafter(value,math.inf) if value else math.nextafter(0.0,math.inf)


def token_radius(text: str) -> float:
    """Decimal serialization half-quantum plus exact decimal-to-binary parse displacement."""
    normalized=text.strip().upper().replace('D','E')
    binary=float(normalized)
    with localcontext() as ctx:
        ctx.prec=100
        parse_displacement=abs(Decimal(normalized)-Decimal.from_float(binary))
        displacement=float(parse_displacement)
    displacement=math.nextafter(displacement,math.inf) if displacement else 0.0
    return up_add(half_quantum(normalized),displacement)


def up_add(a: float, b: float) -> float:
    return math.nextafter(math.fsum((a,b)),math.inf)


def up_mul(a: float, b: float) -> float:
    value=a*b
    return math.nextafter(value,math.inf) if value else math.nextafter(0.0,math.inf)


def complex_product_token_error_bound(w: float, f: float, a: complex, b: complex,
                                      qw: float, qf: float, qa: float, qb: float) -> float:
    """Disk/telescoping input error for w*f*a*conj(b), including phase changes."""
    A=math.nextafter(abs(a),math.inf); B=math.nextafter(abs(b),math.inf)
    F=math.nextafter(abs(f),math.inf); W=math.nextafter(abs(w),math.inf)
    Au=up_add(A,qa); Bu=up_add(B,qb); Fu=up_add(F,qf)
    first=up_mul(up_mul(up_mul(qw,Fu),Au),Bu)
    second=up_mul(up_mul(up_mul(W,qf),Au),Bu)
    third=up_mul(up_mul(W,F),up_add(up_mul(qa,Bu),up_mul(A,qb)))
    return up_add(up_add(first,second),third)


def gamma_bound(n: int) -> float:
    eps=2.220446049250313e-16; ne=n*eps
    if ne>=1: raise ValueError('gamma operation count too large')
    return math.nextafter(ne/(1.0-ne),math.inf) if ne else 0.0


def main() -> int:
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='command',required=True)
    inv=sub.add_parser('inventory'); inv.add_argument('--output',type=Path,required=True)
    ext=sub.add_parser('assemble'); ext.add_argument('inventory',type=Path); ext.add_argument('--inventory-sha256',required=True); ext.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists(): raise SystemExit(f'refusing to overwrite {args.output}')
    if args.command=='inventory':
        result=inventory(); args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'status':result['status'],'path':str(args.output),'sha256':sha(args.output)},sort_keys=True)); return 0
    if sha(args.inventory)!=args.inventory_sha256: raise SystemExit('inventory hash mismatch')
    inv=json.loads(args.inventory.read_text())
    if inv.get('status')!='HASH_INVENTORY_ONLY_NO_NUMERIC_CONTENT_READ': raise SystemExit('unexpected inventory status')
    captures=[assemble_case(inv['cases'][key]) for key in sorted(inv['cases'])]
    payload={'schema':'rslmto-a4-m21-lambda0-assembled-v1','status':'RAW_TOKENS_ASSEMBLED; NO_DERIVATIVE_RESPONSE_ANALYSIS',
             'inventory_path':str(args.inventory),'inventory_sha256':args.inventory_sha256,
             'execution_manifest_sha256':inv['execution_manifest_sha256'],
             'source_sha256':EXPECTED_SOURCE_SHA,'binary_sha256':EXPECTED_BINARY_SHA,'captures':captures}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':payload['status'],'path':str(args.output),'sha256':sha(args.output),'capture_count':len(captures)},sort_keys=True)); return 0


if __name__=='__main__': sys.exit(main())
