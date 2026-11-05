#!/usr/bin/env python3
"""Hash-bound accepted-U2 frozen-state M0.3 replay adapter.

This adapter never runs an SCF calculation. It exports exact rows from one
accepted CAPTURE_ONLY CSV to the private Fortran driver, which invokes the
production calculate_pl/NEWRHO and POISS0/VXC0SP fixture seams.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path

import numpy as np

PL_TOL = 1.0e-10
RHO_TOL = 1.0e-10
ENERGY_TOL = 1.0e-8
FIT_PATH = Path('/private/tmp/a3-rs-fit-compat.json')
FIT_SHA = '48330ec758b680979eb08f1404eabecb9a685be3fa4412a63e14a4936365fb3d'
PRE_INPUT_SHA = 'a05761fe70accf9cf31755c568ac7bda9eff4a79469c6872d297ce7cc8ff040b'
PRE_SEED_SHA = 'f4482119ecb33e1bd443ef9d9481f87c0b99c18a3fea4ee7acde61d49d572e3a'
REPLAY_SOURCE_SHA = 'd91a468a4493da21710fc7f60e745b5f1bef8f20d5d667ae939c31d5ba3061ef'
SHIFTS = (0.0, 0.01)
_BY_KIND_STAGE = {}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def one(rows, kind, stage, *, spin=None, ell=None, order=None, row=None, col=None, band=None):
    bucket = _BY_KIND_STAGE.get((kind, stage))
    if bucket is None:
        bucket = [r for r in rows if r['kind'] == kind and r['stage'] == stage]
    got = [r for r in bucket
           if (spin is None or int(r['spin']) == spin)
           and (ell is None or int(r['l']) == ell)
           and (order is None or int(r['order']) == order)
           and (row is None or int(r['row']) == row)
           and (col is None or int(r['col']) == col)
           and (band is None or int(r['band']) == band)]
    if len(got) != 1:
        raise ValueError(f'expected exactly one {kind}/{stage} row, got {len(got)}')
    return got[0]


def matrix(rows, kind, stage, spin, ell, dim):
    out = np.zeros((dim, dim), dtype=np.complex128)
    seen = set()
    bucket = _BY_KIND_STAGE.get((kind, stage))
    if bucket is None: bucket = [r for r in rows if r['kind'] == kind and r['stage'] == stage]
    for r in bucket:
        if int(r['spin']) == spin and int(r['l']) == ell:
            ij = (int(r['row']), int(r['col']))
            if ij in seen or not all(0 <= x < dim for x in ij):
                raise ValueError(f'duplicate/out-of-range {kind} matrix key {ij}')
            seen.add(ij)
            out[ij] = complex(float(r['re']), float(r['im']))
    if len(seen) != dim * dim:
        raise ValueError(f'incomplete {kind}/{stage} matrix spin={spin} l={ell}')
    return out


def hub_matrix(rows, stage):
    out=np.zeros((18,18),dtype=np.complex128); seen=set()
    bucket = _BY_KIND_STAGE.get(('v_hub', stage))
    if bucket is None: bucket = [r for r in rows if r['kind'] == 'v_hub' and r['stage'] == stage]
    for r in bucket:
        if int(r['type'])==1:
            ij=(int(r['row']),int(r['col']))
            if ij in seen or not all(0<=x<18 for x in ij): raise ValueError(f'duplicate hub matrix key {ij}')
            seen.add(ij); out[ij]=complex(float(r['re']),float(r['im']))
    if len(seen)!=324: raise ValueError(f'incomplete {stage} 18x18 V matrix')
    return out


def ordered(rows, kind, stage, field, n, spin=None):
    found = {}
    bucket = _BY_KIND_STAGE.get((kind, stage))
    if bucket is None: bucket = [r for r in rows if r['kind'] == kind and r['stage'] == stage]
    for r in bucket:
        if spin is None or int(r['spin']) == spin:
            i = int(r['row'])
            if i in found:
                raise ValueError(f'duplicate {kind} index {i}')
            found[i] = r[field]
    if set(found) != set(range(n)):
        raise ValueError(f'incomplete {kind}/{stage}: got {len(found)} expected {n}')
    return [found[i] for i in range(n)]


def export_rows(capture: Path):
    csvp, mp = capture / 'records.csv', capture / 'manifest.final.json'
    manifest = json.loads(mp.read_text())
    if manifest.get('status') != 'COMPLETE_CONVERGED' or manifest.get('backend') != 'RS':
        raise ValueError('capture is not a complete converged KS state')
    if manifest.get('output', {}).get('csv_sha256') != sha(csvp):
        raise ValueError('manifest does not bind records.csv')
    if (manifest.get('input_pre_sha256') != PRE_INPUT_SHA or manifest.get('seed_pre_sha256') != PRE_SEED_SHA or
        manifest.get('source_self_sha256') != REPLAY_SOURCE_SHA):
        raise ValueError('capture does not bind the original U2 PRE deck/seed and reviewed trace source')
    pre=Path(manifest['input_pre_path']); seed=Path(manifest['seed_pre_path'])
    if sha(pre)!=PRE_INPUT_SHA or sha(seed)!=PRE_SEED_SHA:
        raise ValueError('original U2 PRE deck or Fe seed changed since authorization')
    conv=manifest.get('convergence',{})
    tail=conv.get('fixed_u_residual_tail',[])
    if (len(tail)<3 or not all(x.get('physical') is True for x in tail[-3:]) or
        int(tail[-1].get('streak',0))<3):
        raise ValueError('accepted capture lacks the existing converged R05 P/V residual streak')
    if manifest.get('status')!='COMPLETE_CONVERGED':
        raise ValueError('selected U2 run did not record actual convergence')
    rows = list(csv.DictReader(csvp.open(newline='')))
    global _BY_KIND_STAGE
    _BY_KIND_STAGE={}
    for r in rows: _BY_KIND_STAGE.setdefault((r['kind'],r['stage']),[]).append(r)
    for r in rows:
        for key in ('weight','energy','scalar','re','im'):
            if not math.isfinite(float(r[key])):
                raise ValueError(f'non-finite capture value {r["kind"]}/{r["stage"]}/{key}')
    if any(r['kind'].startswith('m0_replay_') or r['stage'].startswith(('A1_SHIFT_','A2_SHIFT_')) for r in rows):
        raise ValueError('capture-only source emitted replay/scaling rows; refuse to inspect it as eligible baseline')
    iters = {int(r['iteration']) for r in rows}
    if len(iters) != 1:
        raise ValueError(f'capture must be one selected iteration, got {iters}')
    it = next(iter(iters))
    selected=manifest.get('output', {}).get('captured_iterations')
    if selected is not None and int(selected[0]) != it:
        raise ValueError('manifest iteration mismatch')
    if int(manifest.get('output', {}).get('row_count', len(rows))) != len(rows):
        raise ValueError('manifest CSV row-count mismatch')
    get = lambda k: one(rows, k, 'RADIAL_REPLAY')['scalar']
    nr, lmax = int(float(get('radial_NR'))), int(float(get('radial_lmax')))
    ncore = len([r for r in rows if r['kind'] == 'radial_replay_ec_seed' and r['stage'] == 'RADIAL_REPLAY'])
    nval = len([r for r in rows if r['kind'] == 'radial_replay_eigen_seed' and r['stage'] == 'RADIAL_REPLAY'])
    if (nr,lmax,ncore,nval) != (501,2,10,6):
        raise ValueError(f'unsupported exact dimensions {nr,lmax,ncore,nval}')
    return rows, manifest, it, nr, ncore, nval, csvp, mp


def prepare(rows, nr, ncore, nval, mode):
    # Exact live operands and all shell inputs are from the one selected trace.
    rofi = ordered(rows, 'radial_replay_mesh', 'RADIAL_REPLAY', 'energy', nr, spin=0)
    # mesh rows repeat mesh per spin; use spin 0 only.
    v = [[one(rows,'radial_replay_v','RADIAL_REPLAY',spin=s, row=i)['scalar'] for s in (0,1)] for i in range(nr)]
    rho_in = [[one(rows,'radial_replay_rho_in','RADIAL_REPLAY',spin=s, row=i)['scalar'] for s in (0,1)] for i in range(nr)]
    ec_seed = ordered(rows,'radial_replay_ec_seed','RADIAL_REPLAY','scalar',ncore)
    evrows = sorted([r for r in rows if r['kind']=='radial_replay_eigen_seed' and r['stage']=='RADIAL_REPLAY'], key=lambda r:int(r['band']))
    if len(evrows)!=nval or [int(r['band']) for r in evrows] != list(range(nval)):
        raise ValueError('incomplete/duplicate EV seeds')
    ev_seed=[r['scalar'] for r in evrows]
    gc=np.zeros((3,2)); c=np.zeros((3,2)); qp=np.zeros((3,2)); sd=np.zeros((3,2)); pl=np.zeros((3,2)); pl_live=np.zeros((3,2)); ql_live=np.zeros((3,3,2))
    for s in (0,1):
        for ell in range(3):
            gc[ell,s]=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=ell,order=0)['scalar'])
            c[ell,s]=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=ell,order=1)['scalar'])
            qp[ell,s]=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=ell,order=2)['scalar'])
            sd[ell,s]=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=ell,order=3)['scalar'])
            pl[ell,s]=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=ell,order=4)['scalar'])
            pl_live[ell,s]=float(one(rows,'radial_replay_pl','RADIAL_REPLAY',spin=s,ell=ell)['scalar'])
            # LIVE_STATE bypasses calculate_pl and uses this exact accepted
            # return; calculated cases start from the captured oldPL operand.
            for m in range(3): ql_live[m,ell,s]=float(one(rows,'radial_replay_ql','RADIAL_REPLAY',spin=s,ell=ell,order=m)['scalar'])
    # The first call is strictly LIVE_STATE only; CLI baseline mode refuses all cases.
    if mode == 'baseline':
        return {'rofi':rofi,'v':v,'rho_in':rho_in,'ec_seed':ec_seed,'ev_seed':ev_seed,'gc':gc,'c':c,'qp':qp,'sd':sd,
                'pl_seed':pl,'pl_live':pl_live,'cases':[('LIVE_STATE',np.zeros(2),False,ql_live)],
                'params':rows}
    # Bind fitted response and allowance to the independently locked A2.2 artifact.
    fit=json.loads(FIT_PATH.read_text())
    if sha(FIT_PATH)!=FIT_SHA: raise ValueError('A2.2 artifact hash mismatch')
    # Derive target/control moments from actual accepted P_OUT and Vin/P_OUT matrices.
    vin=[]; vout=[]; p=[]; m1=[]; q0=[]; center=[]; q2=[]
    for s in (1,2):
        vin.append(hub_matrix(rows,'Pin_H')[4+9*(s-1):9+9*(s-1),4+9*(s-1):9+9*(s-1)])
        vout.append(hub_matrix(rows,'P_OUT')[4+9*(s-1):9+9*(s-1),4+9*(s-1):9+9*(s-1)])
        p.append(matrix(rows,'ldm','P_OUT',s,2,5)); m1.append(matrix(rows,'ldm_m1','M1_OUT',s,2,5))
        q0.append(float(one(rows,'radial_moment','P_OUT',spin=s,ell=2,order=0)['scalar']))
        center.append(float(one(rows,'radial_moment','P_OUT',spin=s,ell=2,order=1)['energy']))
        q2.append(float(one(rows,'radial_moment','P_OUT',spin=s,ell=2,order=2)['energy'])-center[-1]**2)
    # The frozen production replay preserves every captured QL operand and
    # changes only d-shell q2 for corrected/control cases. P_OUT spectral
    # q0/c/q2 are independent target operands used for A2.1/allowances below;
    # their raw-vs-live gaps are reported, never silently copied into QL.
    ql_raw=ql_live.copy()
    shift_local=[]; shift_vout=[]; shift_pt=[]; ql_local=ql_raw.copy(); ql_vout=ql_raw.copy(); ql_pt=ql_raw.copy()
    for idx in range(2):
        def correction(vmat, pmat):
            t=float(np.trace(vmat@pmat).real); cc=center[idx]-t/q0[idx]
            q2loc=q2[idx]-2*float(np.trace(vmat@(m1[idx]-center[idx]*pmat)).real)+float(np.trace(vmat@vmat@pmat).real)-t*t/q0[idx]
            return t,cc,q2loc
        t,cc,q2loc=correction(vin[idx],p[idx]); shift_local.append(t/q0[idx]); ql_local[2,2,idx]=q2loc
        tvout,ccv,q2v=correction(vout[idx],p[idx]); shift_vout.append(tvout/q0[idx]); ql_vout[2,2,idx]=q2v
        pt=p[idx].T; tpt,cpt,q2pt=correction(vin[idx],pt); shift_pt.append(tpt/q0[idx]); ql_pt[2,2,idx]=q2pt
    if min(ql_local[2,2,:]) < 0.0 or min(ql_vout[2,2,:]) < 0.0 or min(ql_pt[2,2,:]) < 0.0:
        raise ValueError('negative q2_local/control moment; refuse to launch M0.3 replay')
    if max(abs(x) for group in (shift_local,shift_vout,shift_pt) for x in group) > 0.12:
        raise ValueError('a requested control shift lies outside the registered A2 scan range')
    d0=float(fit['m03_allowance_at_accepted_u2_targets']['D0_ry'])
    terms={}; a21=0.0; true_r=0.0
    for idx,s in enumerate((1,2)):
        key=f'site=0,type=0,spin={s}'
        coeff=fit['per_site_spin'][key]['r_cubic_coefficients_ascending']
        a0=float(coeff[0]); a21 += abs(a0*float(np.trace(vin[idx]@p[idx]).real))
        vv=float(np.trace(vin[idx]@p[idx]).real)/q0[idx]
        rv=sum(float(coeff[j])*vv**j for j in range(len(coeff)))
        true_r += abs(rv)
        terms[key]={'a0_ry':a0,'q0_ry':q0[idx],'t_ry':float(np.trace(vin[idx]@p[idx]).real),
                    'target_shift_ry':vv,'r_fit_at_target_ry':rv,'q2_local_ry2':float(ql_local[2,2,idx])}
    scalar_r=sum(abs(sum(float(fit['per_site_spin'][f'site=0,type=0,spin={s}']['r_cubic_coefficients_ascending'][j])*.01**j
                           for j in range(len(fit['per_site_spin'][f'site=0,type=0,spin={s}']['r_cubic_coefficients_ascending']))))
                  for s in (1,2))
    if a21>d0: raise ValueError(f'A2.1 failed for this accepted U2 state: {a21:.12e} > D0 {d0:.12e}')
    # Distinct response cases: RAW, true LOCAL, scalar +.01, and wrong V/P controls.
    return {'rofi':rofi,'v':v,'rho_in':rho_in,'ec_seed':ec_seed,'ev_seed':ev_seed,'gc':gc,'c':c,'qp':qp,'sd':sd,
       'pl_seed':pl,'pl_live':pl_live,'fit':fit,'cases':[
       ('RAW_PL',np.zeros(2),True,ql_raw),('LOCAL_PL',np.array(shift_local),True,ql_local),
       ('SCALAR_V_PLUS_001',np.array([.01,.01]),True,ql_raw),('VOUT_CONTROL',np.array(shift_vout),True,ql_vout),
       ('P_TRANSPOSE_CONTROL',np.array(shift_pt),True,ql_pt)], 'preregistration':{'a21_sum_abs_a0t_ry':a21,'D0_ry':d0,
       'a21_pass':a21<=d0,'m03_allowance_at_targets_ry':true_r+d0,'scalar_v_ry':.01,
       'scalar_allowance_ry':scalar_r+d0,'terms':terms}, 'target':{'q0':q0,'center':center,'q2':q2,
       'vin':vin,'vout':vout,'p':p,'m1':m1,'shift_local':shift_local,'shift_vout':shift_vout,'shift_pt':shift_pt}}


def text_input(d, rows, nr, ncore, nval):
    required=('radial_A','radial_B','radial_Z','radial_tolerance','radial_f_core','radial_txc',
              'radial_VHRMAX','radial_B_fsm','radial_hyperfine')
    # Seed arrays and radial samples may also use row=-1; only the explicit
    # scalar-control vocabulary belongs in the unique-control map.
    param_rows=[r for r in rows if r['stage']=='RADIAL_REPLAY' and r['kind'] in required and r['row']=='-1']
    params={}
    for r in param_rows:
        if r['kind'] in params: raise ValueError(f'duplicate captured control {r["kind"]}')
        params[r['kind']]=r['scalar']
    if set(required)-set(params): raise ValueError(f'missing captured runtime operands: {sorted(set(required)-set(params))}')
    a=float(params['radial_A']); b=float(params['radial_B']); z=float(params['radial_Z']); tol=float(params['radial_tolerance']); fcore=int(float(params['radial_f_core']))
    txc=int(float(params['radial_txc'])); vhrmax=float(params['radial_VHRMAX']); bfsm=float(params['radial_B_fsm'])
    if float(params['radial_hyperfine']) != 0.0: raise ValueError('captured hyperfine=true is unsupported by this replay seam')
    lines=[f"2 {nr} {ncore} {nval} {len(d['cases'])} {fcore} {txc}",f"{z:.17g} {a:.17g} {b:.17g} {tol:.17g} {vhrmax:.17g} {bfsm:.17g}",
       ' '.join(d['rofi'])]
    lines += [' '.join([*d['v'][i],*d['rho_in'][i]]) for i in range(nr)]
    lines += [' '.join(d['ec_seed']),' '.join(d['ev_seed'])]
    for s in range(2):
        for ell in range(3): lines.append(' '.join(f'{x:.17g}' for x in (d['gc'][ell,s],d['c'][ell,s],d['qp'][ell,s],d['sd'][ell,s],d['pl_seed'][ell,s],d['pl_live'][ell,s])))
    for _,sh,calc,ql in d['cases']:
        lines.append(f"{sh[0]:.17g} {sh[1]:.17g} {'.true.' if calc else '.false.'}")
        for s in range(2):
            for ell in range(3):
                for m in range(3): lines.append(f'{ql[m,ell,s]:.17g}')
    return '\n'.join(lines)+'\n'


def run(driver: Path, payload: str):
    with tempfile.TemporaryDirectory(prefix='accepted-m03-') as td:
        inp=Path(td)/'input.dat'; out=Path(td)/'output.csv'; inp.write_text(payload)
        p=subprocess.run([str(driver),str(inp),str(out)],capture_output=True,text=True)
        if p.returncode: raise RuntimeError(f'M0.3 driver failed: {p.stderr[-2000:]}')
        return out.read_bytes(), sha(inp)

def run_persistent(driver: Path, payload: str, outdir: Path):
    outdir.mkdir(parents=True,exist_ok=True)
    inp=outdir/'driver.input.txt'; out=outdir/'driver.output.csv'; inp.write_text(payload)
    p=subprocess.run([str(driver),str(inp),str(out)],capture_output=True,text=True)
    (outdir/'driver.stdout.log').write_text(p.stdout)
    (outdir/'driver.stderr.log').write_text(p.stderr)
    if p.returncode: raise RuntimeError(f'M0.3 driver failed: {p.stderr[-2000:]}')
    return out.read_bytes(), sha(inp), sha(out)


def output_cases(payload: bytes):
    result={}
    for row in csv.reader(payload.decode().splitlines()):
        kind,icase=row[0],int(row[1]); c=result.setdefault(icase,{'pl':{},'sumec':{},'sumev':{},'rvh_in':{},'rvh_out':{},'rho':{},'d_rho':{},'components':{},'ev':{}})
        if kind=='PL': c['pl'][(int(row[2]),int(row[3]))]=float(row[4])
        elif kind in ('SUMEC','SUMEV','RVH_IN','RVH_OUT'): c[kind.lower()][int(row[2])]=float(row[3])
        elif kind=='RHO': c['rho'][(int(row[2]),int(row[3]))]=float(row[4]); c['d_rho'][(int(row[2]),int(row[3]))]=float(row[5])
        elif kind=='COMPONENT': c['components'][int(row[2])]=[float(row[3]),float(row[4])]
        elif kind=='EV': c['ev'][int(row[2])-1]=float(row[3])
        else: raise ValueError(f'unexpected M0.3 driver row {kind}')
    return result


def radial_integral(values, rofi, a, b):
    if len(values)!=len(rofi) or len(values)<3 or len(values)%2==0: raise ValueError('invalid Simpson radial extent')
    total=0.0; n=len(values)
    for i,x in enumerate(values):
        ir=i+1; w=1/3 if ir==1 or ir==n else (4/3 if ir%2==0 else 2/3)
        total += w*a*(rofi[i]+b)*x
    return total


def baseline_check(rows, case):
    checks={}
    for name,kind in (('PL','radial_replay_pl'),('rho','radial_replay_rho_out'),
                      ('SUMEC','radial_replay_sumec'),('SUMEV','radial_replay_sumev'),
                      ('EV','radial_replay_eigenvalue')):
        refs=[r for r in rows if r['kind']==kind and r['stage']=='RADIAL_REPLAY']
        if name=='PL':
            ref={(int(r['l']),int(r['spin'])+1):float(r['scalar']) for r in refs}
            got=case['pl']
        elif name=='rho': ref={(int(r['row']),int(r['spin'])+1):float(r['scalar']) for r in refs}; got=case['rho']
        elif name in ('SUMEC','SUMEV'): ref={int(r['spin'])+1:float(r['scalar']) for r in refs}; got=case[name.lower()]
        else: ref={int(r['band']):float(r['scalar']) for r in refs}; got=case['ev']
        if set(ref)!=set(got): raise ValueError(f'LIVE baseline {name} keys differ: {len(ref)} vs {len(got)}')
        err=max((abs(ref[k]-got[k]) for k in ref),default=0.0)
        tol=PL_TOL if name=='PL' else RHO_TOL if name=='rho' else ENERGY_TOL
        checks[name]={'max_abs':err,'tolerance':tol,'pass':err<=tol}
    comp_order=('sumev','sumec','rhovh','rhomu','vnucl','zvnucl','ekin','utot','rhoeps','etot')
    comp_errors={}
    for i,key in enumerate(comp_order,1):
        ref=float(one(rows,'energy_component_'+key,'ATOMSC_ETOT')['scalar'])
        err=abs(case['components'][i][0]-ref)
        comp_errors[key]={'max_abs_ry':err,'tolerance_ry':ENERGY_TOL,'pass':err<=ENERGY_TOL}
    comp_errors['pass']=all(x['pass'] for x in comp_errors.values())
    checks['captured_ATOMSC_input_density_components']=comp_errors
    checks['pass']=all(x['pass'] for x in checks.values() if isinstance(x,dict) and 'pass' in x)
    params={k:float(one(rows,k,'RADIAL_REPLAY')['scalar']) for k in ('radial_A','radial_B')}
    rofi=[float(r['energy']) for r in sorted([r for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0],key=lambda r:int(r['row']))]
    full={s:radial_integral([case['rho'][(i,s)] for i in range(len(rofi))],rofi,params['radial_A'],params['radial_B']) for s in (1,2)}
    dcharge={s:radial_integral([case['d_rho'][(i,s)] for i in range(len(rofi))],rofi,params['radial_A'],params['radial_B']) for s in (1,2)}
    checks['LIVE_STATE_integrals']={'sphere_charge_by_spin':full,'sphere_spin_moment':full[1]-full[2],
       'actual_l2_rho_contribution_charge_by_spin':dcharge,'actual_l2_rho_contribution_total':sum(dcharge.values()),
       'method':'integrate NEWRHO-observed rho-after-minus-before l=2 with ATOMSC logarithmic Simpson metric'}
    ql_gaps=[]; pl_gaps=[]
    for s in (0,1):
        for m in (0,2):
            live=float(one(rows,'radial_replay_ql','RADIAL_REPLAY',spin=s,ell=2,order=m)['scalar'])
            spectral=(float(one(rows,'radial_moment','P_OUT',spin=s+1,ell=2,order=0)['scalar']) if m==0 else
                      float(one(rows,'radial_moment','P_OUT',spin=s+1,ell=2,order=2)['energy'])-
                      float(one(rows,'radial_moment','P_OUT',spin=s+1,ell=2,order=1)['energy'])**2)
            ql_gaps.append({'spin':s+1,'moment':m,'live_ql':live,'p_out_spectral':spectral,'difference':spectral-live})
        pl_live=float(one(rows,'radial_replay_pl','RADIAL_REPLAY',spin=s,ell=2)['scalar'])
        pl_candidate=float(one(rows,'pl_map_input','CALCULATE_PL',spin=s+1,ell=2,order=5)['scalar'])
        pl_gaps.append({'spin':s+1,'live_pl':pl_live,'p_out_calculate_pl':pl_candidate,'difference':pl_candidate-pl_live})
    checks['RAW_PL_vs_LIVE_STATE_operands_diagnostic_only']={'d_shell_ql_q0_q2':ql_gaps,
       'd_shell_pl':pl_gaps,'use':'Report mixing/operand gap; RAW_PL replay uses P_OUT moments while LIVE_STATE baseline uses captured actual PL/QL.'}
    return checks


def cases_summary(rows, d, cases):
    names=[x[0] for x in d['cases']]
    index={name:i+1 for i,name in enumerate(names)}
    required=('RAW_PL','LOCAL_PL','SCALAR_V_PLUS_001','VOUT_CONTROL','P_TRANSPOSE_CONTROL')
    if set(index)!=set(required) or set(cases)!=set(range(1,1+len(required))):
        raise ValueError('M0.3 case table incomplete or contains unexpected cases')
    params={k:float(one(rows,k,'RADIAL_REPLAY')['scalar']) for k in ('radial_A','radial_B')}
    rofi=[float(r['energy']) for r in sorted([r for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0],key=lambda r:int(r['row']))]
    def integrate(case, field, spin):
        return radial_integral([case[field][(i,spin)] for i in range(len(rofi))],rofi,params['radial_A'],params['radial_B'])
    record={}
    component_names={1:'SUMEV',2:'SUMEC',3:'RHOVH',4:'RHOMU',5:'VNUCL',6:'ZVNUCL',
                     7:'EKIN',8:'UTOT',9:'RHOEPS',10:'ETOT'}
    for name in required:
        case=cases[index[name]]
        charge={s:integrate(case,'rho',s) for s in (1,2)}
        dcharge={s:integrate(case,'d_rho',s) for s in (1,2)}
        record[name]={'sumec_by_spin_ry':case['sumec'],'sumev_by_spin_ry':case['sumev'],
          'delta_sumec_from_raw_ry':sum(case['sumec'].values())-sum(cases[index['RAW_PL']]['sumec'].values()),
          'delta_sumev_from_raw_ry':sum(case['sumev'].values())-sum(cases[index['RAW_PL']]['sumev'].values()),
          'sphere_charge_by_spin':charge,'sphere_spin_moment':charge[1]-charge[2],
          'actual_l2_radial_charge_by_spin':dcharge,'actual_l2_radial_charge_total':sum(dcharge.values()),
          'frozen_input_rho_casewise_components':{component_names[i]:v[0] for i,v in case['components'].items()},
          'output_rho_diagnostic_components':{component_names[i]:v[1] for i,v in case['components'].items()},
          'rvh_input_by_spin':case['rvh_in'],'rvh_output_by_spin':case['rvh_out'],
          'ev_by_state_ry':case['ev']}
    raw=sum(cases[index['RAW_PL']]['sumev'].values()); local=sum(cases[index['LOCAL_PL']]['sumev'].values())
    delta=raw-local; ttotal=sum(float(np.trace(v@p).real) for v,p in zip(d['target']['vin'],d['target']['p']))
    allowance=float(d['preregistration']['m03_allowance_at_targets_ry'])
    carrier_error=abs(delta-ttotal)
    scalar=sum(cases[index['RAW_PL']]['sumev'].values())-sum(cases[index['SCALAR_V_PLUS_001']]['sumev'].values())
    scalar_expected=.01*sum(d['target']['q0']); scalar_error=abs(scalar-scalar_expected)
    scalar_allowance=float(d['preregistration']['scalar_allowance_ry'])
    wrong={}
    for name in ('VOUT_CONTROL','P_TRANSPOSE_CONTROL'):
        dsum=raw-sum(cases[index[name]]['sumev'].values())
        dev=abs(dsum-ttotal)
        wrong[name]={'delta_sumev_ry':dsum,'true_target_t_ry':ttotal,'absolute_deviation_from_true_carrier_ry':dev,
                     'allowance_ry':allowance,'discriminating':dev>allowance}
    q2min=min(float(d['target']['q2'][s]-2*np.trace(d['target']['vin'][s]@(d['target']['m1'][s]-d['target']['center'][s]*d['target']['p'][s])).real
                    +np.trace(d['target']['vin'][s]@d['target']['vin'][s]@d['target']['p'][s]).real
                    -(np.trace(d['target']['vin'][s]@d['target']['p'][s]).real**2)/d['target']['q0'][s]) for s in range(2))
    gates={'q2_local_nonnegative':q2min>=0.0,
      'carrier_comparison_with_A2_allowance':carrier_error<=allowance,
      'scalar_shift_response_with_registered_allowance':scalar_error<=scalar_allowance,
      'vout_control_discriminating':wrong['VOUT_CONTROL']['discriminating'],
      'p_transpose_control_discriminating':wrong['P_TRANSPOSE_CONTROL']['discriminating']}
    gates['m03_all_gates_pass']=all(gates.values())
    return {'cases':record,'carrier':{'delta_sumev_raw_minus_local_ry':delta,'true_tr_vin_pout_ry':ttotal,
       'signed_ratio':delta/ttotal if ttotal else None,'absolute_error_ry':carrier_error,'allowance_ry':allowance},
       'scalar_control':{'delta_sumev_ry':scalar,'expected_v_sum_q0_ry':scalar_expected,
       'absolute_error_ry':scalar_error,'allowance_ry':scalar_allowance,'q2_unchanged_by_construction':True},
       'wrong_controls':wrong,'q2_local_min_ry2':q2min,'gates':gates}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--capture',type=Path,required=True); ap.add_argument('--driver',type=Path,required=True)
    ap.add_argument('--mode',choices=('baseline','preregister','cases'),required=True)
    ap.add_argument('--baseline-report',type=Path)
    ap.add_argument('--prereg-report',type=Path)
    ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
    rows,manifest,it,nr,ncore,nval,csvp,mf=export_rows(args.capture)
    source_path=Path(__file__).resolve().parents[3]/'source'/'self.f90'
    if sha(source_path)!=REPLAY_SOURCE_SHA:
        raise ValueError('working source hash differs from the reviewed frozen d91 build')
    if args.mode in ('preregister','cases'):
        if args.baseline_report is None: ap.error('--mode preregister/cases requires --baseline-report')
        gate=json.loads(args.baseline_report.read_text())
        if (gate.get('status')!='PASS_LIVE_BASELINE' or gate.get('capture_csv_sha256')!=sha(csvp) or
            gate.get('manifest_sha256')!=sha(mf) or gate.get('driver_sha256')!=sha(args.driver)):
            raise ValueError('M0.3 case replay requires a matching passing hash-bound LIVE baseline report')
    if args.mode=='cases':
        if args.prereg_report is None: ap.error('--mode cases requires --prereg-report')
        preg=json.loads(args.prereg_report.read_text())
        if (preg.get('status')!='PASS_PREREGISTERED_M03' or preg.get('capture_csv_sha256')!=sha(csvp) or
            preg.get('manifest_sha256')!=sha(mf) or preg.get('driver_sha256')!=sha(args.driver) or
            preg.get('baseline_report_sha256')!=sha(args.baseline_report) or
            preg.get('fit_sha256')!=FIT_SHA):
            raise ValueError('M0.3 cases require the same accepted capture and a hash-bound preregistration PASS')
    d=prepare(rows,nr,ncore,nval,'baseline' if args.mode=='baseline' else 'cases')
    if args.mode=='preregister':
        report={'schema':'rslmto-accepted-u2-private-m03-preregistration-v1','status':'PASS_PREREGISTERED_M03',
          'iteration':it,'capture_csv_sha256':sha(csvp),'manifest_sha256':sha(mf),'driver_sha256':sha(args.driver),
          'source_sha256':sha(Path(__file__).resolve().parents[3]/'source'/'self.f90'),
          'baseline_report_sha256':sha(args.baseline_report),'fit_path':str(FIT_PATH),'fit_sha256':FIT_SHA,
          'convergence_and_R05':manifest.get('convergence'),
          'pre_case_gate':d['preregistration'],'target_operands':{'q0':d['target']['q0'],'center':d['target']['center'],
          'q2':d['target']['q2'],'shift_local':d['target']['shift_local'],'shift_vout':d['target']['shift_vout'],
          'shift_pt':d['target']['shift_pt']}}
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'mode':args.mode,'status':report['status'],'iteration':it,'output':str(args.output),'output_sha256':sha(args.output)},sort_keys=True))
        return
    raw,input_hash,driver_output_hash=run_persistent(args.driver.resolve(),text_input(d,rows,nr,ncore,nval),args.output.parent)
    case_data=output_cases(raw)
    if args.mode=='baseline':
        if set(case_data)!={1}: raise ValueError('baseline phase emitted unexpected extra replay cases')
        checks=baseline_check(rows,case_data[1])
        csvout=args.output.with_suffix('.csv'); csvout.write_bytes(raw)
        report={'schema':'rslmto-accepted-u2-private-m03-live-baseline-v1','status':'PASS_LIVE_BASELINE' if checks['pass'] else 'FAIL_LIVE_BASELINE',
          'iteration':it,'capture_csv_sha256':sha(csvp),'manifest_sha256':sha(mf),'driver_sha256':sha(args.driver),
          'source_sha256':sha(Path(__file__).resolve().parents[3]/'source'/'self.f90'),'input_sha256':input_hash,
          'PRE_input_sha256':manifest.get('input_pre_sha256'),'PRE_seed_sha256':manifest.get('seed_pre_sha256'),
          'convergence_and_R05':manifest.get('convergence'),
          'runtime_radial_controls':{k:next(r['scalar'] for r in rows if r['kind']==k and r['stage']=='RADIAL_REPLAY')
             for k in ('radial_VHRMAX','radial_B_fsm','radial_txc','radial_hyperfine')},
          'baseline_output_path':str(csvout.resolve()),'baseline_output_sha256':sha(csvout),'checks':checks}
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    else:
        if sha(FIT_PATH)!=FIT_SHA: raise ValueError('A2.2 model artifact hash mismatch')
        csvout=args.output.with_suffix('.csv'); csvout.write_bytes(raw)
        summary=cases_summary(rows,d,case_data)
        report={'schema':'rslmto-accepted-u2-private-m03-cases-v1','status':'MEASURED_REPLAY_DIAGNOSTIC',
          'iteration':it,'capture_csv_sha256':sha(csvp),'manifest_sha256':sha(mf),'driver_sha256':sha(args.driver),
          'source_sha256':sha(Path(__file__).resolve().parents[3]/'source'/'self.f90'),'input_sha256':input_hash,
          'baseline_report_sha256':sha(args.baseline_report),'prereg_report_sha256':sha(args.prereg_report),
          'convergence_and_R05':manifest.get('convergence'),
          'fit_path':str(FIT_PATH),'fit_sha256':FIT_SHA,
          'pre_case_gate':d['preregistration'],
          'measurements':summary,
          'case_names':[x[0] for x in d['cases']],'case_output_path':str(csvout.resolve()),'case_output_sha256':sha(csvout)}
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'mode':args.mode,'status':report['status'],'iteration':it,'output':str(args.output),
      'output_sha256':sha(args.output),'capture_csv_sha256':sha(csvp),'driver_sha256':sha(args.driver)},sort_keys=True))


# RS-specific A3.5 controls. The KS runner above remains untouched in its source artifact.
def prepare_rs(rows, nr, ncore, nval):
    d=prepare(rows,nr,ncore,nval,'baseline')
    ql_raw=d['cases'][0][3].copy()
    vin=[]; p=[]; m1=[]; q0=[]; center=[]; q2=[]; t=[]; q2local=[]
    for sp in (1,2):
        vm=hub_matrix(rows,'Pin_H'); off=4+9*(sp-1)
        vin.append(vm[off:off+5,off:off+5])
        p.append(matrix(rows,'ldm','P_OUT',sp,2,5))
        m1.append(matrix(rows,'ldm_m1','M1_OUT',sp,2,5))
        row=one(rows,'m0_carrier_input','FROZEN_REPLAY',spin=sp,ell=2)
        q0.append(float(row['scalar']));center.append(float(row['energy']));q2.append(float(row['re']))
        tv=float(np.trace(vin[-1]@p[-1]).real);t.append(tv)
        ql=q2[-1]-2*float(np.trace(vin[-1]@(m1[-1]-center[-1]*p[-1])).real)+float(np.trace(vin[-1]@vin[-1]@p[-1]).real)-tv*tv/q0[-1]
        q2local.append(ql)
        observed=float(one(rows,'m0_carrier_local','FROZEN_REPLAY',spin=sp,ell=2)['re'])
        if abs(observed-ql)>5e-12: raise ValueError(f'RS local q2 reconstruction mismatch spin {sp}: {observed-ql}')
        ql_raw[2,2,sp-1]=ql
    if min(q2local)<0: raise ValueError('negative RS local q2')
    d0=5.056500278044496e-6
    fit=json.loads(FIT_PATH.read_text())
    shifts=[t[i]/q0[i] for i in range(2)]
    ql_wrong=ql_raw.copy(); wrong_t=[]; wrong_q2=[]; wrong_shift=[]
    for i,sp in enumerate((1,2)):
        ov=vin[1-i];tw=float(np.trace(ov@p[i]).real);dd=tw/q0[i]
        qw=q2[i]-2*float(np.trace(ov@(m1[i]-center[i]*p[i])).real)+float(np.trace(ov@ov@p[i]).real)-tw*tw/q0[i]
        if qw<0: raise ValueError(f'negative opposite-spin q2 spin {sp}')
        ql_wrong[2,2,i]=qw;wrong_t.append(tw);wrong_q2.append(qw);wrong_shift.append(dd)
    ql_local=ql_raw.copy(); ql_local[2,2,:]=q2local
    ql_sign=ql_local.copy()
    cases=[('RAW_PL',np.zeros(2),True,ql_raw),('LOCAL_PL',np.array(shifts),True,ql_local),
      ('SIGN_WRONG',-np.array(shifts),True,ql_sign),('OPPOSITE_SPIN',np.array(wrong_shift),True,ql_wrong),
      ('SCALAR_V_PLUS_001',np.array([.01,.01]),True,ql_raw.copy())]
    # M0.3 allowance uses the RS fit at each true x plus D0; scalar allowance is residual-only plus D0.
    Atrue=Ascalar=Asign=Aspin=0.0; A21=0.0; terms={}
    for i,sp in enumerate((1,2)):
        coeff=fit['per_site_spin'][f'site=0,type=0,spin={sp}']['r_cubic_coefficients_ascending']
        a0=fit['per_site_spin'][f'site=0,type=0,spin={sp}']['a0_ry']
        r=lambda x:sum(float(coeff[j])*x**j for j in range(4))
        x=shifts[i];Atrue+=abs(r(x));Ascalar+=abs(r(.01));Asign+=abs(r(-x));Aspin+=abs(r(wrong_shift[i]));A21+=abs(float(a0)*t[i])
        terms[str(sp)]={'t_ry':t[i],'q0':q0[i],'x':x,'rfit_x_ry':r(x),'predicted_true_delta_ry':t[i]+r(x),
          'A1_old8_intercept_ry':a0,'opposite_spin_t_ry':wrong_t[i],'opposite_spin_shift_ry':wrong_shift[i],
          'opposite_spin_q2_ry2':wrong_q2[i],'opposite_spin_rfit_ry':r(wrong_shift[i]),
          'opposite_spin_predicted_delta_ry':wrong_t[i]+r(wrong_shift[i]),
          'sign_wrong_shift_ry':-x,'sign_wrong_rfit_ry':r(-x),'sign_wrong_predicted_delta_ry':-t[i]+r(-x),
          'scalar_plus_001_predicted_delta_ry':.01*q0[i]+r(.01),'q2_local_ry2':q2local[i]}
    if A21>d0: raise ValueError(f'RS A2.1 old-eight criterion failed: {A21} > {d0}')
    return {**d,'fit':fit,'cases':cases,'preregistration':{'A2_1_old8_product_ry':A21,'D0_RS_ry':d0,'A_true_ry':Atrue+d0,'A_scalar_plus_001_ry':Ascalar+d0,'A_sign_ry':Asign+d0,'A_opposite_spin_ry':Aspin+d0,'terms':terms},
       'target':{'vin':vin,'p':p,'m1':m1,'t':t,'q0':q0,'center':center,'q2':q2,'q2_local':q2local,'wrong_t':wrong_t,'wrong_q2':wrong_q2,'wrong_shift':wrong_shift,'local_shift':shifts}}

def cases_summary_rs(rows,d,cases):
    names=[x[0] for x in d['cases']]; index={n:i+1 for i,n in enumerate(names)}
    if set(cases)!=set(range(1,len(names)+1)): raise ValueError('incomplete RS A3 response cases')
    params={k:float(one(rows,k,'RADIAL_REPLAY')['scalar']) for k in ('radial_A','radial_B')}
    rofi=[float(r['energy']) for r in sorted([r for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0],key=lambda r:int(r['row']))]
    def integ(case,field,sp): return radial_integral([case[field][(i,sp)] for i in range(len(rofi))],rofi,params['radial_A'],params['radial_B'])
    names_comp={1:'SUMEV',2:'SUMEC',3:'RHOVH',4:'RHOMU',5:'VNUCL',6:'ZVNUCL',7:'EKIN',8:'UTOT',9:'RHOEPS',10:'ETOT'}
    rec={}
    for name in names:
        c=cases[index[name]];rho={sp:integ(c,'rho',sp) for sp in (1,2)};drho={sp:integ(c,'d_rho',sp) for sp in (1,2)}
        rec[name]={'sumec_by_spin_ry':c['sumec'],'sumev_by_spin_ry':c['sumev'],
          'sphere_charge_by_spin':rho,'sphere_spin_moment':rho[1]-rho[2],
          'actual_l2_radial_charge_by_spin':drho,'actual_l2_radial_charge_total':sum(drho.values()),
          'frozen_input_rho_casewise_components':{names_comp[i]:v[0] for i,v in c['components'].items()},
          'output_rho_diagnostic_components':{names_comp[i]:v[1] for i,v in c['components'].items()},'ev_by_state_ry':c['ev']}
    raw=cases[index['RAW_PL']]['sumev']; local=cases[index['LOCAL_PL']]['sumev']; true_t=d['target']['t']
    deltas={sp:raw[sp]-local[sp] for sp in (1,2)}; errs={sp:deltas[sp]-true_t[sp-1] for sp in (1,2)}
    scalar={sp:raw[sp]-cases[index['SCALAR_V_PLUS_001']]['sumev'][sp] for sp in (1,2)}
    expected_scalar={sp:d['preregistration']['terms'][str(sp)]['scalar_plus_001_predicted_delta_ry'] for sp in (1,2)}
    wrong={}
    for name in ('SIGN_WRONG','OPPOSITE_SPIN'):
        got={sp:raw[sp]-cases[index[name]]['sumev'][sp] for sp in (1,2)}
        exp={sp:d['preregistration']['terms'][str(sp)]['sign_wrong_predicted_delta_ry' if name=='SIGN_WRONG' else 'opposite_spin_predicted_delta_ry'] for sp in (1,2)}
        wrong[name]={'delta_sumev_by_spin_ry':got,'predicted_delta_by_spin_ry':exp,'errors_by_spin_ry':{sp:got[sp]-exp[sp] for sp in (1,2)}}
    carrier_abs_error=abs(sum(deltas.values())-sum(true_t))
    ideal_scalar={sp:.01*d['target']['q0'][sp-1] for sp in (1,2)}
    scalar_ideal_errors={sp:scalar[sp]-ideal_scalar[sp] for sp in (1,2)}
    scalar_ideal_total_error=abs(sum(scalar.values())-sum(ideal_scalar.values()))
    wrong={}; true_total=sum(true_t)
    for name in ('SIGN_WRONG','OPPOSITE_SPIN'):
        got={sp:raw[sp]-cases[index[name]]['sumev'][sp] for sp in (1,2)}
        exp={sp:d['preregistration']['terms'][str(sp)]['sign_wrong_predicted_delta_ry' if name=='SIGN_WRONG' else 'opposite_spin_predicted_delta_ry'] for sp in (1,2)}
        fit_error=abs(sum(got.values())-sum(exp.values()))
        detectability=abs(sum(got.values())-true_total)
        allowance=d['preregistration']['A_sign_ry' if name=='SIGN_WRONG' else 'A_opposite_spin_ry']
        wrong[name]={'delta_sumev_by_spin_ry':got,'predicted_delta_by_spin_ry':exp,
          'errors_by_spin_ry':{sp:got[sp]-exp[sp] for sp in (1,2)},'residual_inclusive_fit_error_ry':fit_error,
          'relevant_allowance_ry':allowance,'fit_within_allowance':fit_error<=allowance,
          'true_carrier_total_ry':true_total,'absolute_separation_from_true_carrier_ry':detectability,
          'detectability_allowance_A_true_ry':d['preregistration']['A_true_ry'],
          'detectable':detectability>d['preregistration']['A_true_ry']}
    carrier_pass=carrier_abs_error<=d['preregistration']['A_true_ry']
    scalar_pass=scalar_ideal_total_error<=d['preregistration']['A_scalar_plus_001_ry']
    q2min=min(d['target']['q2_local'])
    gates={'carrier_literal_error_within_A_true':carrier_pass,'scalar_deviation_from_vq0_within_A_scalar':scalar_pass,
      'q2_local_nonnegative':q2min>=0.0,'sign_wrong_control_detectable':wrong['SIGN_WRONG']['detectable'],
      'opposite_spin_control_detectable':wrong['OPPOSITE_SPIN']['detectable'],
      'sign_wrong_response_matches_registered_prediction':wrong['SIGN_WRONG']['fit_within_allowance'],
      'opposite_spin_response_matches_registered_prediction':wrong['OPPOSITE_SPIN']['fit_within_allowance']}
    gates['all_required_gates_pass']=all(gates.values())
    return {'cases':rec,'carrier':{'delta_sumev_by_spin_ry':deltas,'true_trace_by_spin_ry':{1:true_t[0],2:true_t[1]},'errors_by_spin_ry':errs,
       'signed_ratio_by_spin':{sp:deltas[sp]/true_t[sp-1] for sp in (1,2)},'absolute_error_ry':carrier_abs_error,
       'allowance_ry':d['preregistration']['A_true_ry'],'pass':carrier_pass},
      'scalar_plus_001':{'delta_sumev_by_spin_ry':scalar,'ideal_vq0_by_spin_ry':ideal_scalar,
       'deviation_from_vq0_by_spin_ry':scalar_ideal_errors,'absolute_deviation_from_v_sum_q0_ry':scalar_ideal_total_error,
       'allowance_ry':d['preregistration']['A_scalar_plus_001_ry'],'pass':scalar_pass,
       'residual_inclusive_prediction_by_spin_ry':expected_scalar,
       'deviation_from_residual_inclusive_prediction_by_spin_ry':{sp:scalar[sp]-expected_scalar[sp] for sp in (1,2)},
       'q2_input_policy':'captured RAW q2 invariant under V=vI'},
      'wrong_controls':wrong,'q2_local_min_ry2':q2min,'gates':gates}



def _main_rs():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--capture',type=Path,required=True)
    ap.add_argument('--driver',type=Path,required=True)
    ap.add_argument('--registration',type=Path,required=True)
    ap.add_argument('--mode',choices=('baseline','cases'),required=True)
    ap.add_argument('--baseline-report',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    reg=json.loads(args.registration.read_text())
    def sh(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if reg.get('schema')!='a3-rs-m03-preregistration-v1' or reg.get('status')!='LOCKED_BEFORE_RESPONSES':
        raise ValueError('missing finalized RS M0.3 preregistration')
    if reg.get('runner_sha256')!=sh(__file__) or reg.get('driver_sha256')!=sh(args.driver):
        raise ValueError('runner/driver hash differs from frozen preregistration')
    if sh(reg['candidate_path'])!=reg.get('candidate_sha256'):
        raise ValueError('prediction candidate hash differs from frozen preregistration')
    if sh(FIT_PATH)!=FIT_SHA or reg.get('fit_adapter_sha256')!=FIT_SHA:
        raise ValueError('RS fit adapter changed')
    rows,manifest,it,nr,ncore,nval,csvp,mf=export_rows(args.capture)
    src=Path(__file__).resolve().parents[3]/'source'/'self.f90'
    if manifest.get('source_self_sha256')!=REPLAY_SOURCE_SHA: raise ValueError('capture manifest source hash mismatch')
    if manifest.get('executable_sha256')!='b7ca883ec736211552ca1b8df7e2d519c7806459df9a6238dd597c58f3c49f7f': raise ValueError('capture manifest binary hash mismatch')
    if it!=6 or manifest.get('run_id')!='a3_rs_u2_m03_d91_state_3e8': raise ValueError('capture selected iteration/run ID mismatch')
    if manifest.get('convergence',{}).get('accepted') is not True: raise ValueError('capture lacks accepted convergence')
    if args.mode=='baseline':
        d=prepare(rows,nr,ncore,nval,'baseline')
        raw,input_hash,driver_output_hash=run_persistent(args.driver.resolve(),text_input(d,rows,nr,ncore,nval),args.output.parent)
        cdata=output_cases(raw)
        if set(cdata)!={1}: raise ValueError('baseline phase returned extra cases')
        checks=baseline_check(rows,cdata[1]); csvout=args.output.with_suffix('.csv');csvout.write_bytes(raw)
        report={'schema':'a3-rs-m03-live-anchor-v1','status':'PASS_LIVE_BASELINE' if checks['pass'] else 'FAIL_LIVE_BASELINE',
          'iteration':it,'capture_csv_sha256':sh(csvp),'manifest_sha256':sh(mf),'driver_sha256':sh(args.driver),
          'runner_sha256':sh(__file__),'source_sha256':manifest.get('source_self_sha256'),'current_source_file_sha256':sh(src),'input_sha256':input_hash,
          'capture_convergence':manifest.get('convergence'),'checks':checks,
          'raw_output_path':str(csvout),'raw_output_sha256':sh(csvout),'driver_output_sha256':driver_output_hash}
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    else:
        if args.baseline_report is None: ap.error('--mode cases requires --baseline-report')
        base=json.loads(args.baseline_report.read_text())
        if (base.get('status')!='PASS_LIVE_BASELINE' or base.get('capture_csv_sha256')!=sh(csvp) or
            base.get('manifest_sha256')!=sh(mf) or base.get('driver_sha256')!=sh(args.driver)):
            raise ValueError('RS A3 cases require the exact passing LIVE anchor')
        d=prepare_rs(rows,nr,ncore,nval)
        raw,input_hash,driver_output_hash=run_persistent(args.driver.resolve(),text_input(d,rows,nr,ncore,nval),args.output.parent); cdata=output_cases(raw)
        result=cases_summary_rs(rows,d,cdata)
        csvout=args.output.with_suffix('.csv');csvout.write_bytes(raw)
        report={'schema':'a3-rs-m03-cases-v1','status':'MEASURED_RS_A3_M03','iteration':it,
          'capture_csv_sha256':sh(csvp),'manifest_sha256':sh(mf),'driver_sha256':sh(args.driver),
          'runner_sha256':sh(__file__),'baseline_report_sha256':sh(args.baseline_report),
          'source_sha256':manifest.get('source_self_sha256'),'current_source_file_sha256':sh(src),'fit_adapter_sha256':sh(FIT_PATH),'input_sha256':input_hash,
          'precase':d['preregistration'],'measurements':result,
          'case_names':[x[0] for x in d['cases']],'raw_output_path':str(csvout),'raw_output_sha256':sh(csvout),'driver_output_sha256':driver_output_hash}
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'mode':args.mode,'status':report['status'],'output':str(args.output),
      'output_sha256':sh(args.output),'capture_csv_sha256':sh(csvp),'driver_sha256':sh(args.driver)},sort_keys=True))

if __name__=='__main__':
    _main_rs()
