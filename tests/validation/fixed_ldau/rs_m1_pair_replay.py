#!/usr/bin/env python3
"""Hash-bound RS local-M1 LIVE anchor and P_OUT spectral/local pair.

The adapter only serializes captured operands to the unchanged accepted-M03
Fortran fixture driver.  It does not run SCF or alter production source.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CAPTURE = Path('/private/tmp/a4-rs-m1-local-u2')
TRACE = CAPTURE / 'manifest.trace.final.json'
SCF = CAPTURE / 'scf_result.json'
FIT = ROOT / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_fit_results.json'
FIT_CSV = ROOT / 'docs/validation/fixed_ldau/evidence/a3_rs_a22_six_shift_responses.csv'
DRIVER = Path('/private/tmp/rslmto-localmom-a2-u2-d91-bin/accepted_m03_replay_driver')
OUTROOT = Path('/private/tmp/a4-rs-m1-pair')
PRE = ROOT / 'docs/validation/fixed_ldau/evidence/a4_rs_m1_pair_preregistration.json'
PL_TOL, RHO_TOL, ENERGY_TOL = 1.0e-10, 1.0e-10, 1.0e-8
D0_RS = 5.056500278044496e-6
E_NUM = 2.0e-8
ARITH = 2.0e-12

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load_runner():
    p=ROOT/'tests/validation/fixed_ldau/rs_a3_m03_runner.py'
    spec=importlib.util.spec_from_file_location('rs_a3_runner',p)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def load_rows():
    rows=list(csv.DictReader((CAPTURE/'records.csv').open(newline='')))
    m=json.loads(TRACE.read_text()); s=json.loads(SCF.read_text())
    if m['backend']!='RS' or m['csv_sha256']!=sha(CAPTURE/'records.csv'):
        raise ValueError('capture manifest/CSV mismatch or non-RS backend')
    if m['iteration']!=103 or s['status']!='SCF_CAPTURE_PASS; PAIRED_DUAL_REPLAY_PENDING':
        raise ValueError('unexpected local RS SCF capture')
    if s['records_csv_sha256']!=sha(CAPTURE/'records.csv') or not s['capture_only']:
        raise ValueError('capture result does not bind capture-only records')
    if sha(DRIVER)!='0c1c3c12e8501f4642af87168c611d261fb5a34cc4747aaef9c8420c2c6449c4':
        raise ValueError('production replay driver hash mismatch')
    if sha(FIT)!='0f0b543f81f53bc1c7fdfc46d3b404adef1ca455b8366a90e87154afd2b0987c':
        raise ValueError('RS scalar fit report hash mismatch')
    if sha(FIT_CSV)!='293db40ac192edb5f4428e7a87fff85ad99061c6c51d602c1291a67aac46cdb0':
        raise ValueError('RS scalar response CSV hash mismatch')
    return rows,m,s

def build_payloads():
    r=load_runner(); rows,m,s=load_rows(); nr,ncore,nval=495,10,6
    r._BY_KIND_STAGE={}
    for row in rows: r._BY_KIND_STAGE.setdefault((row['kind'],row['stage']),[]).append(row)
    common=r.prepare(rows,nr,ncore,nval,'baseline')
    # Stage-specific LIVE_STATE anchor: exact captured RADIAL_REPLAY QL and PL.
    anchor=common.copy(); anchor['cases']=[('LIVE_STATE',np.zeros(2),False,common['cases'][0][3])]
    anchor_txt=r.text_input(anchor,rows,nr,ncore,nval)

    # P_OUT operands are selected independently from the same iteration.  ql
    # is (q0,0,q2); source mapping leaves ql(2)=0 and stores corrected q2 in
    # radial_parameter order=3.  Bind every shell, then vary only d center/q2.
    ql=np.zeros((3,3,2)); q0=[]; cl=[]; q2l=[]; cs=[]; q2s=[]; t=[]; shifts=[]
    vins=[]; ps=[]; m1s=[]; shell_checks=[]
    for sp in (1,2):
        for ell in range(3):
            q0v=float(r.one(rows,'radial_moment','P_OUT',spin=sp,ell=ell,order=0)['scalar'])
            ql1=float(r.one(rows,'ql','P_OUT',spin=sp,ell=ell,order=0)['scalar'])
            q2v=float(r.one(rows,'radial_parameter','P_OUT',spin=sp,ell=ell,order=3)['scalar'])
            ql[0,ell,sp-1]=ql1; ql[1,ell,sp-1]=0.; ql[2,ell,sp-1]=q2v
            live0=float(r.one(rows,'radial_replay_ql','RADIAL_REPLAY',spin=sp-1,ell=ell,order=0)['scalar'])
            live2=float(r.one(rows,'radial_replay_ql','RADIAL_REPLAY',spin=sp-1,ell=ell,order=2)['scalar'])
            if abs(q0v-ql1)>2e-14: raise ValueError('P_OUT q0 disagreement between QL and radial_moment')
            shell_checks.append({'spin':sp,'l':ell,'pout_q0':ql1,'live_q0':live0,'q0_pout_minus_live':ql1-live0,
                                 'pout_q2':q2v,'live_q2':live2,'q2_pout_minus_live':q2v-live2,
                                 'pout_ql_middle':0.0})
        q0v=float(r.one(rows,'radial_moment','P_OUT',spin=sp,ell=2,order=0)['scalar'])
        center=float(r.one(rows,'radial_moment','P_OUT',spin=sp,ell=2,order=1)['energy'])
        q2loc=float(r.one(rows,'radial_parameter','P_OUT',spin=sp,ell=2,order=3)['scalar'])
        gc=float(r.one(rows,'pl_map_input','CALCULATE_PL',spin=sp,ell=2,order=0)['scalar'])
        if abs(gc-center)>2e-14: raise ValueError('P_OUT gravity center and radial center differ')
        vm=r.hub_matrix(rows,'Pin_H'); off=4+9*(sp-1); vin=vm[off:off+5,off:off+5]
        p=r.matrix(rows,'ldm','P_OUT',sp,2,5); m1=r.matrix(rows,'ldm_m1','M1_OUT',sp,2,5)
        tv=float(np.trace(vin@p).real); c_spec=center+tv/q0v
        q2spec=q2loc+2*float(np.trace(vin@(m1-c_spec*p)).real)-float(np.trace(vin@vin@p).real)+tv*tv/q0v
        trc=float(np.trace(m1).real/q0v)
        if abs(trc-c_spec)>1e-12 or q2loc<0 or q2spec<0: raise ValueError('P_OUT inverse moment identity failed')
        q0.append(q0v); cl.append(center); q2l.append(q2loc); cs.append(c_spec); q2s.append(q2spec); t.append(tv)
        shifts.append(-tv/q0v); vins.append(vin); ps.append(p); m1s.append(m1)
        # The P_OUT center is also exactly the center passed to calculate_pl.
    ql_local=ql.copy(); ql_spectral=ql.copy()
    for i in range(2): ql_spectral[2,2,i]=q2s[i]
    pair=common.copy(); pair['cases']=[('P_OUT_SPECTRAL',np.array(shifts),True,ql_spectral),
                                       ('P_OUT_LOCAL',np.zeros(2),True,ql_local)]
    pair_txt=r.text_input(pair,rows,nr,ncore,nval)
    fit=json.loads(FIT.read_text()); per={}; rr=[]
    for i,sp in enumerate((1,2)):
        coeff=fit['per_spin'][str(sp)]['fit']['r_cubic_coefficients_ascending']
        x=t[i]/q0[i]; residual=sum(float(a)*x**j for j,a in enumerate(coeff)); rr.append(residual)
        per[str(sp)]={'q0':q0[i],'c_local_ry':cl[i],'q2_local_ry2':q2l[i],'t_ry':t[i],
                      'c_spectral_ry':cs[i],'q2_spectral_ry2':q2s[i],
                      'tr_M1_over_q0_ry':float(np.trace(m1s[i]).real/q0[i]),
                      'driver_center_shift_ry':shifts[i], 'driver_applied_center_change_ry':-shifts[i],
                      'x_t_over_q0_ry':x,'r_fit_ry':residual,'predicted_delta_spectral_minus_local_ry':t[i]+residual,
                      'trace_vin_p_ry':t[i]}
    predictions={'per_spin':per,'sum_t_ry':sum(t),'sum_r_fit_ry':sum(rr),
        'predicted_delta_sumev_spectral_minus_local_ry':sum(t)+sum(rr),
        'literal_carrier_error_allowance_A_RS_ry':sum(abs(x) for x in rr)+D0_RS,
        'D0_RS_once_ry':D0_RS,'numerical_guard_two_case_energy_tol_ry':E_NUM,
        'arithmetic_crosscheck_guard_ry':ARITH,
        'criteria':{'literal_carrier_abs_delta_minus_sum_t_le_A_RS_plus_E_num':True,
                    'model_prediction_abs_delta_minus_sum_t_plus_r_le_E_num_plus_arithmetic':True,
                    'q2_spectral_and_local_nonnegative':True}}
    provenance={'capture_manifest_sha256':sha(TRACE),'capture_csv_sha256':sha(CAPTURE/'records.csv'),
       'scf_result_sha256':sha(SCF),'capture_pre_manifest_sha256':sha(CAPTURE/'manifest.pre.json'),
       'seed_pre_sha256':'f4482119ecb33e1bd443ef9d9481f87c0b99c18a3fea4ee7acde61d49d572e3a',
       'input_pre_sha256':'651263443ff8e121b2831f7d5a4d13a8f8f1c5a2f70e18191ff679f2df8b4862',
       'driver_sha256':sha(DRIVER),'fit_report_sha256':sha(FIT),'fit_response_csv_sha256':sha(FIT_CSV),
       'capture_source_revision':m.get('source_sha'),'capture_executable_sha256':m.get('executable_sha256'),
       'capture_input_sha256':m.get('input_sha256')}
    return rows,m,s,anchor_txt,pair_txt,provenance,predictions,shell_checks

def execute(payload, outdir):
    outdir.mkdir(parents=True,exist_ok=True); inp=outdir/'driver.input.txt'; out=outdir/'driver.output.csv'
    inp.write_text(payload)
    import subprocess
    cp=subprocess.run([str(DRIVER),str(inp),str(out)],capture_output=True,text=True)
    (outdir/'stdout.txt').write_text(cp.stdout); (outdir/'stderr.txt').write_text(cp.stderr)
    if cp.returncode: raise RuntimeError(cp.stderr[-2000:])
    return out.read_bytes()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=('lock','anchor','pair')); args=ap.parse_args()
    rows,m,s,anchor,pair,prov,pred,shell=build_payloads()
    if args.mode=='lock':
        lock={'schema':'a4-rs-m1-pout-pair-prereg-v1','status':'LOCKED_BEFORE_RESPONSE',
              'scope':'RS Gate M1 paired frozen P_OUT replay; separate actual LIVE_STATE anchor first',
              'provenance':prov,'payloads':{'live_anchor_sha256':hashlib.sha256(anchor.encode()).hexdigest(),
                 'p_out_pair_sha256':hashlib.sha256(pair.encode()).hexdigest()},
              'shell_q0_q2_timing_diagnostics':shell,'p_out_inverse_reconstruction':pred,
              'method':{'center_convention':'c_local=c_spectral−t/q0; inverse c_spectral=c_local+t/q0. The driver applies c_input−center_shift, so spectral center_shift=−t/q0 and local center_shift=0.',
                 'q2_inverse':'q2_spectral=q2_local+2 ReTr[V_in(M1−c_spectral P)]−ReTr[V_in²P]+t²/q0.',
                 'V_in':'Captured Pin_H 5×5 d block; P and M1 are captured P_OUT and M1_OUT from the same iteration.',
                 'pair_changes':'Only d-shell center and q2 differ. All P_OUT q0 and QL middle components, all other shells, and common RADIAL_REPLAY mesh/V/rho/EC/EV/PL operands are byte-identical in the generated pair payload.',
                 'time_levels':'LIVE anchor uses only RADIAL_REPLAY arrays. P_OUT pair is a distinct frozen P_OUT reconstruction; P_OUT−LIVE moment gaps are diagnostic and are not compared to LIVE anchor tolerances.'},
              'guards':{'pl':PL_TOL,'rho':RHO_TOL,'energy':ENERGY_TOL,'two_case_energy_numerical_guard':E_NUM,'arithmetic_guard':ARITH},
              'fit_model':'RS own locked A2.2 six-shift/cubic residual model plus D0_RS once; A_true=sum(abs(r_fit(x_s)))+D0_RS.',
              'adapter_sha256':sha(Path(__file__)),'driver_sha256':prov['driver_sha256'],
              'source_revision':'f122a001abd30d8213867b199b702d2d73ba915d',
              'response_status':'NO_RESPONSES_RUN'}
        PRE.write_text(json.dumps(lock,indent=2,sort_keys=True)+'\n'); print(json.dumps({'lock':str(PRE),'sha256':sha(PRE)})); return
    lock=json.loads(PRE.read_text())
    if lock['status']!='LOCKED_BEFORE_RESPONSE' or lock['provenance']!=prov: raise ValueError('lock provenance mismatch')
    if lock.get('adapter_sha256')!=sha(Path(__file__)): raise ValueError('adapter changed after preregistration')
    if hashlib.sha256(anchor.encode()).hexdigest()!=lock['payloads']['live_anchor_sha256'] or hashlib.sha256(pair.encode()).hexdigest()!=lock['payloads']['p_out_pair_sha256']:
        raise ValueError('generated payload differs from preregistered lock')
    r=load_runner()
    if args.mode=='anchor':
        raw=execute(anchor,OUTROOT/'live_anchor'); case=r.output_cases(raw)[1]
        checks=r.baseline_check(rows,case)
        report={'schema':'a4-rs-m1-live-anchor-results-v1','status':'PASS' if checks['pass'] else 'FAIL',
          'pre_registration_sha256':sha(PRE),'capture_manifest_sha256':prov['capture_manifest_sha256'],
          'driver_sha256':prov['driver_sha256'],'driver_input_sha256':sha(OUTROOT/'live_anchor/driver.input.txt'),
          'driver_output_sha256':sha(OUTROOT/'live_anchor/driver.output.csv'),'checks':checks}
        dest=OUTROOT/'live_anchor/report.json'; dest.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n'); print(json.dumps({'report':str(dest),'sha256':sha(dest),'status':report['status']}))
        if not checks['pass']: raise SystemExit('LIVE anchor failed; P_OUT pair remains unrun')
    else:
        anchor_report=OUTROOT/'live_anchor/report.json'
        if not anchor_report.exists() or json.loads(anchor_report.read_text()).get('status')!='PASS': raise ValueError('LIVE anchor has not passed')
        raw=execute(pair,OUTROOT/'pout_pair'); cases=r.output_cases(raw); a,b=cases[1],cases[2]
        diff=sum(a['sumev'].values())-sum(b['sumev'].values()); sumt=lock['p_out_inverse_reconstruction']['sum_t_ry']; pred=lock['p_out_inverse_reconstruction']['predicted_delta_sumev_spectral_minus_local_ry']; A=lock['p_out_inverse_reconstruction']['literal_carrier_error_allowance_A_RS_ry']
        q2min=min(lock['p_out_inverse_reconstruction']['per_spin'][s]['q2_spectral_ry2'] for s in ('1','2'))
        q2localmin=min(lock['p_out_inverse_reconstruction']['per_spin'][s]['q2_local_ry2'] for s in ('1','2'))
        report={'schema':'a4-rs-m1-pout-pair-results-v1','status':'MEASURED',
          'pre_registration_sha256':sha(PRE),'live_anchor_report_sha256':sha(anchor_report),
          'capture_manifest_sha256':prov['capture_manifest_sha256'],'driver_sha256':prov['driver_sha256'],
          'driver_input_sha256':sha(OUTROOT/'pout_pair/driver.input.txt'),'driver_output_sha256':sha(OUTROOT/'pout_pair/driver.output.csv'),
          'delta_sumev_spectral_minus_local_ry':diff,'sum_t_ry':sumt,'predicted_delta_sumev_ry':pred,
          'literal_carrier_error_ry':abs(diff-sumt),'A_RS_ry':A,'numerical_guard_ry':E_NUM,
          'model_prediction_error_ry':abs(diff-pred),'model_prediction_guard_ry':E_NUM+ARITH,
          'q2_spectral_min_ry2':q2min,'q2_local_min_ry2':q2localmin,'criteria':{'literal_carrier':abs(diff-sumt)<=A+E_NUM,
           'fit_prediction':abs(diff-pred)<=E_NUM+ARITH,'q2_nonnegative':q2min>=0 and q2localmin>=0},
          'cases':{name:{'sumec_by_spin_ry':c['sumec'],'sumev_by_spin_ry':c['sumev'],
           'pl':{f'{sp}:{ell}':v for (ell,sp),v in c['pl'].items()},
           'sphere_charge_by_spin':{sp:r.radial_integral([c['rho'][(i,sp)] for i in range(495)],
            [float(x['energy']) for x in sorted([x for x in rows if x['kind']=='radial_replay_mesh' and int(x['spin'])==0],key=lambda x:int(x['row']))],
            float(r.one(rows,'radial_A','RADIAL_REPLAY')['scalar']),float(r.one(rows,'radial_B','RADIAL_REPLAY')['scalar'])) for sp in (1,2)},
           'actual_l2_rho_contribution_charge_by_spin':{sp:r.radial_integral([c['d_rho'][(i,sp)] for i in range(495)],
            [float(x['energy']) for x in sorted([x for x in rows if x['kind']=='radial_replay_mesh' and int(x['spin'])==0],key=lambda x:int(x['row']))],
            float(r.one(rows,'radial_A','RADIAL_REPLAY')['scalar']),float(r.one(rows,'radial_B','RADIAL_REPLAY')['scalar'])) for sp in (1,2)},
           'sphere_spin_moment':None,'frozen_input_rho_components':{
             {1:'SUMEV',2:'SUMEC',3:'RHOVH',4:'RHOMU',5:'VNUCL',6:'ZVNUCL',7:'EKIN',8:'UTOT',9:'RHOEPS',10:'ETOT'}[k]:v[0] for k,v in c['components'].items()},
           'output_rho_components_diagnostic':{
             {1:'SUMEV',2:'SUMEC',3:'RHOVH',4:'RHOMU',5:'VNUCL',6:'ZVNUCL',7:'EKIN',8:'UTOT',9:'RHOEPS',10:'ETOT'}[k]:v[1] for k,v in c['components'].items()},'ev_by_state_ry':c['ev']}
          for name,c in [('P_OUT_SPECTRAL',a),('P_OUT_LOCAL',b)]}}
        for c in report['cases'].values(): c['sphere_spin_moment']=c['sphere_charge_by_spin'][1]-c['sphere_charge_by_spin'][2]
        dest=OUTROOT/'pout_pair/report.json'; dest.parent.mkdir(parents=True,exist_ok=True); dest.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n'); print(json.dumps({'report':str(dest),'sha256':sha(dest),'criteria':report['criteria']}))

if __name__=='__main__': main()
