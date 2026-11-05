#!/usr/bin/env python3
"""RD2 Ni frozen radial replay adapter, bound to actual LOCAL_HALL Vin.

Modes are ordered baseline -> preregister -> cases. Baseline and cases invoke
the frozen radial helper; every mode fails closed without a separate exact-hash
authorization file. This adapter never launches an SCF solver.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from collections import Counter
from pathlib import Path

REPO = Path('/Users/andersb/Jobb/rslmto_fable_v4b')
CAPTURE = Path('/private/tmp/c52-rd2-ni-u3j07-capture/run')
BINDING = Path('/private/tmp/c52-rd2-private/vin_binding_nr501.json')
DRIVER = Path('/private/tmp/c52-rd2-private/nr501_adapter/build_proven/accepted_m03_replay_driver_nr501')
PLAN = Path('/private/tmp/c52-rd2-private/rd2_replay_plan.json')
AUTH = Path('/private/tmp/c52-rd2-private/execution_authorization.json')
FIT = REPO / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_fit_results.json'
FIT_PREREG = REPO / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_preregistration_v3.json'
FIT_SHA = '0f0b543f81f53bc1c7fdfc46d3b404adef1ca455b8366a90e87154afd2b0987c'
FIT_PREREG_SHA = '46ca7042678dc5b74ad15a7c2bc5f66786ca9191f70d2395405e78c69b5b1525'
LOCAL_ORACLE = REPO / 'tests/validation/fixed_ldau/local_moment_capture_oracle.py'
LOCAL_ORACLE_SHA = '98d14387be47c3cea4a4803ce331a1c92cb8c0d9d5d3ca412b5491634b68645d'
RUN_ID = 'rd2_ni_u3j07_endpoint23_fullgreen_i1'
ITERATION = 1
TARGET_CLUSTER_SITE = 1
TARGET_TYPE = 2
BUILD_ID = 1


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_json(path: Path):
    return json.loads(path.read_text())


def require_authorization(args, mode: str) -> dict:
    """Fail closed unless a separately reviewed, hash-bound authorization exists."""
    auth_path = args.authorization
    need(auth_path.is_file(), f'no explicit authorization file for mode={mode}')
    auth = load_json(auth_path)
    need(auth.get('status') == 'ROOT_AND_AUDITOR_CLEARED', 'authorization is not cleared')
    need(auth.get('mode') == mode, f'authorization mode mismatch: expected {mode}')
    need(Path(auth.get('plan_path', '')).resolve() == PLAN.resolve(), 'authorization plan path mismatch')
    need(sha(PLAN) == auth.get('plan_sha256'), 'authorization does not bind current plan')
    need(Path(auth.get('adapter_path', '')).resolve() == Path(__file__).resolve(), 'authorization adapter path mismatch')
    need(sha(Path(__file__)) == auth.get('adapter_sha256'), 'authorization does not bind current adapter')
    need(Path(auth.get('binding_path', '')).resolve() == args.binding.resolve(), 'authorization binding path mismatch')
    need(sha(args.binding) == auth.get('binding_sha256'), 'authorization does not bind current binding')
    need(Path(auth.get('driver_path', '')).resolve() == args.driver.resolve(), 'authorization driver path mismatch')
    need(sha(args.driver) == auth.get('driver_sha256'), 'authorization does not bind current driver')
    if mode in ('preregister', 'cases'):
        need(args.baseline_report is not None and args.baseline_report.is_file(), 'matching baseline report required')
        need(auth.get('baseline_report_sha256') == sha(args.baseline_report), 'authorization baseline report mismatch')
    if mode == 'cases':
        need(args.prereg_report is not None and args.prereg_report.is_file(), 'matching preregistration required')
        need(auth.get('prereg_report_sha256') == sha(args.prereg_report), 'authorization preregistration mismatch')
    return auth


def module_at(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def modules():
    runner = module_at('rd2_b1_radial_runner_nr501', Path('/private/tmp/c52-rd2-private/nr501_adapter/rs_a3_m03_runner_nr501.py'))
    c2 = module_at('rd2_c2_audit_parser', REPO / 'tests/validation/fixed_ldau/c2_hubbard_operator_audit.py')
    return runner, c2


def local_moment_oracle():
    need(sha(LOCAL_ORACLE) == LOCAL_ORACLE_SHA, 'reviewed local-moment formula helper changed')
    return module_at('rd2_reviewed_local_moment_oracle', LOCAL_ORACLE)


def load_bound_rows(capture_dir: Path, binding_path: Path):
    runner, c2 = modules()
    binding = load_json(binding_path)
    need(binding.get('status') == 'STRUCTURAL_BINDING_ONLY_NO_RD2_GATE_OR_REPLAY_EXECUTED',
         'binding must be the reviewed structural-only artifact')
    csv_path = capture_dir / 'records.csv'
    audit_path = capture_dir / 'c2_hubbard_audit.dat'
    need(sha(csv_path) == binding['iteration_and_capture']['csv_sha256'], 'capture CSV differs from binding')
    need(sha(audit_path) == binding['source']['sha256'], 'LOCAL_HALL audit differs from binding')
    rows = list(csv.DictReader(csv_path.open(newline='', encoding='utf-8')))
    need(len(rows) == 552196 and {r['run_id'] for r in rows} == {RUN_ID}, 'wrong capture row count/run id')
    need({int(r['iteration']) for r in rows} == {ITERATION}, 'wrong capture iteration')
    for n, row in enumerate(rows, 1):
        for key in ('weight', 'energy', 'scalar', 're', 'im'):
            need(math.isfinite(float(row[key])), f'nonfinite capture field row={n} key={key}')
    runner._BY_KIND_STAGE = {}
    for row in rows:
        runner._BY_KIND_STAGE.setdefault((row['kind'], row['stage']), []).append(row)
    nr = int(float(runner.one(rows, 'radial_NR', 'RADIAL_REPLAY')['scalar']))
    lmax = int(float(runner.one(rows, 'radial_lmax', 'RADIAL_REPLAY')['scalar']))
    ncore = sum(r['kind'] == 'radial_replay_ec_seed' and r['stage'] == 'RADIAL_REPLAY' for r in rows)
    nval = sum(r['kind'] == 'radial_replay_eigen_seed' and r['stage'] == 'RADIAL_REPLAY' for r in rows)
    need((nr, lmax, ncore, nval) == (501, 2, 10, 6), 'unsupported frozen radial replay dimensions')
    rec = c2.choose_record(c2.read_records(audit_path), 'LOCAL_HALL', TARGET_CLUSTER_SITE, BUILD_ID)
    need(int(rec['ptype']) == TARGET_TYPE and rec['symbol'] == 'Ni-imp' and int(rec['count']) == 1,
         'LOCAL_HALL target mapping mismatch')
    need(len(rec['block']) == 50, 'LOCAL_HALL exact d correction block incomplete')
    vin = []
    for spin in (1, 2):
        mat = [[0j] * 5 for _ in range(5)]
        for (s, i, j), entry in rec['block'].items():
            if s == spin:
                mat[i-1][j-1] = entry['correction']
        need(sum(1 for (s, _, _), entry in rec['block'].items() if s == spin and entry['count'] == 1) == 25,
             f'incomplete LOCAL_HALL target for spin={spin}')
        vin.append(runner.np.asarray(mat, dtype=complex))
    return runner, rows, nr, ncore, nval, vin, rec, binding


def baseline_check_dynamic(runner, rows, case, nr):
    """Same accepted baseline checks, using the capture's actual radial extent."""
    checks = {}
    for name, kind in (('PL','radial_replay_pl'),('rho','radial_replay_rho_out'),
                       ('SUMEC','radial_replay_sumec'),('SUMEV','radial_replay_sumev'),
                       ('EV','radial_replay_eigenvalue')):
        refs=[r for r in rows if r['kind']==kind and r['stage']=='RADIAL_REPLAY']
        if name=='PL':
            ref={(int(r['l']),int(r['spin'])+1):float(r['scalar']) for r in refs}; got=case['pl']
        elif name=='rho':
            ref={(int(r['row']),int(r['spin'])+1):float(r['scalar']) for r in refs}; got=case['rho']
        elif name in ('SUMEC','SUMEV'):
            ref={int(r['spin'])+1:float(r['scalar']) for r in refs}; got=case[name.lower()]
        else:
            ref={int(r['band']):float(r['scalar']) for r in refs}; got=case['ev']
        if set(ref)!=set(got): raise ValueError(f'baseline {name} key sets differ')
        err=max((abs(ref[k]-got[k]) for k in ref),default=0.0)
        tol=runner.PL_TOL if name=='PL' else runner.RHO_TOL if name=='rho' else runner.ENERGY_TOL
        checks[name]={'max_abs':err,'tolerance':tol,'pass':err<=tol}
    order=('sumev','sumec','rhovh','rhomu','vnucl','zvnucl','ekin','utot','rhoeps','etot')
    comp_errors={}
    for i,key in enumerate(order,1):
        ref=float(runner.one(rows,'energy_component_'+key,'ATOMSC_ETOT')['scalar'])
        err=abs(case['components'][i][0]-ref)
        comp_errors[key]={'max_abs_ry':err,'tolerance_ry':runner.ENERGY_TOL,'pass':err<=runner.ENERGY_TOL}
    comp_errors['pass']=all(x['pass'] for x in comp_errors.values() if isinstance(x,dict))
    checks['captured_ATOMSC_input_density_components']=comp_errors
    checks['pass']=all(v['pass'] for v in checks.values() if isinstance(v,dict) and 'pass' in v)
    params={k:float(runner.one(rows,k,'RADIAL_REPLAY')['scalar']) for k in ('radial_A','radial_B')}
    rofi=[float(r['energy']) for r in sorted([r for r in rows if r['kind']=='radial_replay_mesh' and int(r['spin'])==0],key=lambda r:int(r['row']))]
    if len(rofi)!=nr: raise ValueError(f'captured radial mesh size {len(rofi)} != radial_NR {nr}')
    full={s:runner.radial_integral([case['rho'][(i,s)] for i in range(nr)],rofi,params['radial_A'],params['radial_B']) for s in (1,2)}
    dcharge={s:runner.radial_integral([case['d_rho'][(i,s)] for i in range(nr)],rofi,params['radial_A'],params['radial_B']) for s in (1,2)}
    checks['LIVE_STATE_integrals']={'sphere_charge_by_spin':full,'sphere_spin_moment':full[1]-full[2],
        'actual_l2_rho_contribution_charge_by_spin':dcharge,'actual_l2_rho_contribution_total':sum(dcharge.values()),
        'method':'same captured logarithmic Simpson metric, with the exact 501-point mesh extent'}
    return checks


def baseline(args):
    require_authorization(args, 'baseline')
    runner, rows, nr, ncore, nval, _, _, binding = load_bound_rows(args.capture_dir, args.binding)
    data = runner.prepare(rows, nr, ncore, nval, 'baseline')
    payload = runner.text_input(data, rows, nr, ncore, nval)
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=False)
    raw, input_sha, output_sha = runner.run_persistent(args.driver, payload, outdir)
    cases = runner.output_cases(raw)
    checks = baseline_check_dynamic(runner, rows, cases[1], nr)
    report = {'schema': 'c5-rd2-frozen-replay-baseline-v1',
              'status': 'PASS' if checks['pass'] else 'FAIL',
              'purpose': 'B1-linked frozen radial helper reproduction of captured raw state; not an RD2 response case',
              'binding_sha256': sha(args.binding), 'capture_csv_sha256': binding['iteration_and_capture']['csv_sha256'],
              'driver_sha256': sha(args.driver), 'driver_input_sha256': input_sha,
              'driver_output_sha256': output_sha, 'authorization_sha256': sha(args.authorization), 'checks': checks}
    (outdir / 'baseline.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'report': str(outdir/'baseline.json'),
                      'report_sha256': sha(outdir/'baseline.json')}, sort_keys=True))


def preregister(args):
    require_authorization(args, 'preregister')
    runner, rows, nr, ncore, nval, vin, local_rec, binding = load_bound_rows(args.capture_dir, args.binding)
    oracle = local_moment_oracle()
    baseline_report = load_json(args.baseline_report)
    need(baseline_report.get('status') == 'PASS', 'RD2 cases require a passing frozen raw baseline')
    need(baseline_report.get('binding_sha256') == sha(args.binding), 'baseline used different capture binding')
    need(baseline_report.get('capture_csv_sha256') == binding['iteration_and_capture']['csv_sha256'],
         'baseline used different capture CSV')
    fit = load_json(FIT); fit_pre = load_json(FIT_PREREG)
    need(sha(FIT) == FIT_SHA and sha(FIT_PREREG) == FIT_PREREG_SHA, 'locked RB3/A3.2 fit hashes changed')
    need(fit.get('prereg_sha256') == sha(FIT_PREREG), 'fit does not bind its preregistration')
    d0 = float(fit_pre['accepted_no_u_state']['D0_RS_ry'])
    shaped = {'schema': 'rslmto-rb3-fit-shape-adapter-v1', 'source_fit_sha256': sha(FIT),
              'source_prereg_sha256': sha(FIT_PREREG), 'D0_RS_ry': d0, 'per_site_spin': {}}
    for spin in ('1', '2'):
        model = fit['per_spin'][spin]['fit']
        shaped['per_site_spin'][f'site=0,type=0,spin={spin}'] = {
            'r_cubic_coefficients_ascending': model['r_cubic_coefficients_ascending'],
            'a0_ry': model['e_even_fit_coefficients_ascending'][0]}
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=False)
    fit_shape = outdir / 'rb3_fit_shape_adapter.json'
    fit_shape.parent.mkdir(parents=True, exist_ok=True)
    fit_shape.write_text(json.dumps(shaped, indent=2, sort_keys=True) + '\n')
    runner.FIT_PATH = fit_shape; runner.FIT_SHA = sha(fit_shape)

    # The shared runner supplies the exact captured radial/PL inputs. Its
    # generic type-1 selector is deliberately not called for Vin here.
    data = runner.prepare(rows, nr, ncore, nval, 'baseline')
    ql_raw = data['cases'][0][3].copy()
    p_mats, m1_mats, q0, center, q2 = [], [], [], [], []
    for spin in (1, 2):
        p = oracle.matrix(rows, kind='ldm', stage='P_OUT', iteration=ITERATION,
                          site=0, ell=2, spin=spin, dim=5)
        m1 = oracle.matrix(rows, kind='ldm_m1', stage='M1_OUT', iteration=ITERATION,
                           site=0, ell=2, spin=spin, dim=5)
        qrow = runner.one(rows, 'm0_carrier_input', 'FROZEN_REPLAY', spin=spin, ell=2)
        need(int(qrow['site']) == 0 and int(qrow['type']) == TARGET_TYPE-1,
             f'captured q0/c/q2 not bound to recursion site/type spin={spin}')
        q0.append(float(qrow['scalar'])); center.append(float(qrow['energy'])); q2.append(float(qrow['re']))
        need(abs(float(runner.np.trace(p).real)-q0[-1]) <= 1.0e-4, f'P trace/q0 mismatch spin={spin}')
        p_mats.append(p); m1_mats.append(m1)
    carrier, shifts, q2prime = [], [], []
    ql_local = ql_raw.copy()
    for idx in range(2):
        corr = oracle.local_correction_no_controls(p_mats[idx], m1_mats[idx], vin[idx], q0[idx], center[idx], q2[idx])
        t = float(runner.np.trace(vin[idx] @ p_mats[idx]).real)
        need(math.isfinite(t), f'nonfinite target carrier operand spin={idx+1}')
        carrier.append(t); shifts.append(t/q0[idx]); q2prime.append(float(corr['q2_local_ry2']))
        ql_local[2, 2, idx] = q2prime[idx]
    # Preserve the accepted RB3 prepare_rs construction: ql_raw carries the
    # exact-V q2prime, and RAW/LOCAL/SCALAR cases share it; only center shift
    # differs for LOCAL and SCALAR.
    ql_raw = ql_local
    ql_local = ql_raw.copy()
    wrong_shift, wrong_q2, wrong_t = [], [], []
    for idx in range(2):
        other = vin[1-idx]
        t = float(runner.np.trace(other @ p_mats[idx]).real)
        corr = oracle.local_correction_no_controls(p_mats[idx], m1_mats[idx], other, q0[idx], center[idx], q2[idx])
        wrong_t.append(t); wrong_shift.append(t/q0[idx]); wrong_q2.append(float(corr['q2_local_ry2']))
    ql_sign=ql_local.copy(); ql_swap=ql_raw.copy()
    for idx in range(2): ql_swap[2,2,idx]=wrong_q2[idx]
    fit_terms={}; a_true=a_scalar=a_sign=a_swap=0.0
    for idx, spin in enumerate(('1','2')):
        coeff=shaped['per_site_spin'][f'site=0,type=0,spin={spin}']['r_cubic_coefficients_ascending']
        residual=lambda x: sum(float(c)*(x**j) for j,c in enumerate(coeff))
        a_true += abs(residual(shifts[idx])); a_scalar += abs(residual(.01))
        a_sign += abs(residual(-shifts[idx])); a_swap += abs(residual(wrong_shift[idx]))
        fit_terms[spin]={'x_true_ry':shifts[idx], 'r_true_ry':residual(shifts[idx]),
                         'x_spin_swap_ry':wrong_shift[idx], 'r_spin_swap_ry':residual(wrong_shift[idx]),
                         'r_sign_wrong_ry':residual(-shifts[idx]), 'r_scalar_plus_001_ry':residual(.01),
                         't_true_ry':carrier[idx], 't_spin_swap_ry':wrong_t[idx],
                         'sign_wrong_predicted_delta_ry':-carrier[idx]+residual(-shifts[idx]),
                         'opposite_spin_predicted_delta_ry':wrong_t[idx]+residual(wrong_shift[idx]),
                         'scalar_plus_001_predicted_delta_ry':.01*q0[idx]+residual(.01)}
    prereg = {'schema':'c5-rd2-response-preregistration-v1',
        'status':'PREPARED_AFTER_RAW_BASELINE; AWAIT_AUDITOR_REVIEW_BEFORE_CASE_REPLAY',
        'baseline_report_path':str(args.baseline_report), 'baseline_report_sha256':sha(args.baseline_report),
        'binding_path':str(args.binding), 'binding_sha256':sha(args.binding),
        'capture_csv_sha256':binding['iteration_and_capture']['csv_sha256'],
        'Vin_source':'LOCAL_HALL additive BLOCK operand at runtime nlim=1,type2; not v_hub type1',
        'Vin_exact_replay_operand_path':str(Path('/private/tmp/c52-rd2-private/vin_binding.json')),
        'Vin_exact_replay_operand_sha256':sha(Path('/private/tmp/c52-rd2-private/vin_binding.json')),
        'q0_center_q2_source':'captured recursion-site m0_carrier_input; validated type=1 (zero based maps to type2) and site0; P_OUT P trace matches q0',
        'P_M1_source':'P_OUT and M1_OUT matrices at recursion-site site0, per spin',
        'case_definitions':[
            {'name':'RAW_PL','shift_by_spin_ry':[0.0,0.0],'q2':'exact-V corrected q2prime, matching accepted prepare_rs control construction'},
            {'name':'LOCAL_PL','shift_by_spin_ry':shifts,'q2':'correct_d_shell_moments formula with exact LOCAL_HALL Vin'},
            {'name':'SIGN_WRONG','shift_by_spin_ry':[-x for x in shifts],'q2':'same corrected q2; wrong sign control'},
            {'name':'OPPOSITE_SPIN','shift_by_spin_ry':wrong_shift,'q2':'same formula using opposite-spin Vin; spin-swap control'},
            {'name':'SCALAR_V_PLUS_001','shift_by_spin_ry':[.01,.01],'q2':'same exact-V corrected q2prime as RAW_PL/LOCAL_PL; scalar changes center shift only'}],
        'target_operands':{'q0':q0,'center_ry':center,'q2_raw_ry2':q2,'carrier_trace_vin_p_ry':carrier,
                           'local_shift_by_spin_ry':shifts,'q2_prime_by_spin_ry2':q2prime,
                           'spin_swap_carrier_ry':wrong_t,'spin_swap_shift_ry':wrong_shift,
                           'spin_swap_q2_prime_ry2':wrong_q2},
        'registered_fit_control':{'fit_sha256':sha(FIT),'prereg_sha256':sha(FIT_PREREG),'D0_RS_ry':d0,
            'per_spin_fit_responses':fit_terms,'A_true_ry':a_true+d0,'A_scalar_plus_001_ry':a_scalar+d0,
            'A_sign_ry':a_sign+d0,'A_opposite_spin_ry':a_swap+d0,
            'mapping':'per_spin fit coefficients go only to generic RB3 spin slots; old site/type labels do not select Vin'},
        'summary_seed':{'case_names':['RAW_PL','LOCAL_PL','SIGN_WRONG','OPPOSITE_SPIN','SCALAR_V_PLUS_001'],
            'target':{'t':carrier,'q0':q0,'q2_local':q2prime},
            'preregistration':{'A_true_ry':a_true+d0,'A_scalar_plus_001_ry':a_scalar+d0,
                'A_sign_ry':a_sign+d0,'A_opposite_spin_ry':a_swap+d0,'terms':fit_terms}},
        'case_driver_input_path':str(outdir/'driver.input.txt'),
        'authorization_sha256':sha(args.authorization),
        'no_gate_claim':'These are preregistered operands/allowances only; no driver case output has been measured.'}
    data['fit']=shaped; data['preregistration']={'A_true_ry':a_true+d0,'A_scalar_plus_001_ry':a_scalar+d0,
        'A_sign_ry':a_sign+d0,'A_opposite_spin_ry':a_swap+d0,'terms':fit_terms}
    data['target']={'vin':vin,'p':p_mats,'m1':m1_mats,'t':carrier,'q0':q0,'center':center,'q2':q2,
                    'q2_local':q2prime,'wrong_t':wrong_t,'wrong_q2':wrong_q2,'wrong_shift':wrong_shift,
                    'local_shift':shifts}
    data['cases']=[('RAW_PL',runner.np.zeros(2),True,ql_raw.copy()),
                   ('LOCAL_PL',runner.np.asarray(shifts),True,ql_local),
                   ('SIGN_WRONG',-runner.np.asarray(shifts),True,ql_sign),
                   ('OPPOSITE_SPIN',runner.np.asarray(wrong_shift),True,ql_swap),
                   ('SCALAR_V_PLUS_001',runner.np.array([.01,.01]),True,ql_raw.copy())]
    payload=runner.text_input(data,rows,nr,ncore,nval)
    (outdir/'driver.input.txt').write_text(payload)
    prereg['case_driver_input_sha256'] = sha(outdir/'driver.input.txt')
    (outdir/'preregistration.json').write_text(json.dumps(prereg,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':prereg['status'],'preregistration':str(outdir/'preregistration.json'),
                      'preregistration_sha256':sha(outdir/'preregistration.json'),
                      'driver_input':str(outdir/'driver.input.txt'),'driver_input_sha256':sha(outdir/'driver.input.txt')},sort_keys=True))


def cases(args):
    require_authorization(args, 'cases')
    runner, rows, nr, ncore, nval, vin, _, binding = load_bound_rows(args.capture_dir, args.binding)
    prereg=load_json(args.prereg_report)
    need(prereg.get('status') == 'PREPARED_AFTER_RAW_BASELINE; AWAIT_AUDITOR_REVIEW_BEFORE_CASE_REPLAY',
         'case execution requires the audited preregistration status')
    need(prereg.get('binding_sha256') == sha(args.binding), 'preregistration binding mismatch')
    input_path=Path(prereg['case_driver_input_path'])
    need(sha(input_path)==prereg['case_driver_input_sha256'], 'case input payload hash mismatch')
    outdir=args.output_dir; outdir.mkdir(parents=True,exist_ok=False)
    p=__import__('subprocess').run([str(args.driver),str(input_path),str(outdir/'driver.output.csv')],
                                  capture_output=True,text=True)
    (outdir/'driver.stdout.log').write_text(p.stdout); (outdir/'driver.stderr.log').write_text(p.stderr)
    need(p.returncode==0, f'frozen radial replay helper exited {p.returncode}')
    raw=(outdir/'driver.output.csv').read_bytes()
    case_data=runner.output_cases(raw)
    names=prereg['summary_seed']['case_names']
    summary_data={'cases':[(name,None,None,None) for name in names],
                  'target':prereg['summary_seed']['target'],
                  'preregistration':prereg['summary_seed']['preregistration']}
    measurements=runner.cases_summary_rs(rows,summary_data,case_data)
    ratios=measurements['carrier']['signed_ratio_by_spin']
    carrier_pass=all(.99 <= float(ratios[s]) <= 1.01 for s in (1,2))
    q2_pass=float(measurements['q2_local_min_ry2']) >= 0.0
    scalar=measurements['scalar_plus_001']
    scalar_pass=float(scalar['absolute_deviation_from_v_sum_q0_ry']) <= float(scalar['allowance_ry'])
    carrier_error=float(measurements['carrier']['absolute_error_ry'])
    separation={name:float(measurements['wrong_controls'][name]['absolute_separation_from_true_carrier_ry'])
                for name in ('SIGN_WRONG','OPPOSITE_SPIN')}
    controls_pass=all(value > 10.0*carrier_error for value in separation.values())
    gates={'carrier_ratio_each_spin_within_1_percent':carrier_pass,
           'q2_prime_nonnegative':q2_pass,
           'scalar_control_within_registered_fit_allowance':scalar_pass,
           'sign_and_spin_swap_each_exceed_10x_carrier_error':controls_pass,
           'wrong_control_separation_ry':separation,
           'measured_carrier_error_ry':carrier_error}
    gates['all_rd2_gate_conditions_pass']=all(gates[k] for k in (
        'carrier_ratio_each_spin_within_1_percent','q2_prime_nonnegative',
        'scalar_control_within_registered_fit_allowance','sign_and_spin_swap_each_exceed_10x_carrier_error'))
    report={'schema':'c5-rd2-measured-replay-cases-v1',
            'status':'GATE_MEASURED_PASS' if gates['all_rd2_gate_conditions_pass'] else 'GATE_MEASURED_FAIL',
            'binding_sha256':sha(args.binding),'prereg_report_sha256':sha(args.prereg_report),
            'baseline_report_sha256':sha(args.baseline_report),'driver_sha256':sha(args.driver),
            'authorization_sha256':sha(args.authorization),'measurements':measurements,'formal_gates':gates,
            'driver_output_path':str(outdir/'driver.output.csv'),'driver_output_sha256':sha(outdir/'driver.output.csv'),
            'stdout_sha256':sha(outdir/'driver.stdout.log'),'stderr_sha256':sha(outdir/'driver.stderr.log')}
    (outdir/'rd2_cases_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'report':str(outdir/'rd2_cases_report.json'),
                      'report_sha256':sha(outdir/'rd2_cases_report.json'),
                      'driver_output_sha256':report['driver_output_sha256']},sort_keys=True))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=('baseline','preregister','cases'),required=True)
    p.add_argument('--capture-dir',type=Path,default=CAPTURE)
    p.add_argument('--binding',type=Path,default=BINDING)
    p.add_argument('--driver',type=Path,default=DRIVER)
    p.add_argument('--baseline-report',type=Path)
    p.add_argument('--prereg-report',type=Path)
    p.add_argument('--authorization',type=Path,default=AUTH)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    need(sha(args.driver)==load_json(args.binding)['driver']['sha256'], 'driver does not match binding')
    if args.mode=='baseline': baseline(args)
    elif args.mode=='preregister':
        need(args.baseline_report is not None,'--baseline-report required')
        preregister(args)
    else:
        need(args.prereg_report is not None,'--prereg-report required')
        cases(args)


if __name__=='__main__':
    main()
