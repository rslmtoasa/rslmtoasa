#!/usr/bin/env python3
"""Run the C5.4 eight-point scalar scan on the frozen NR501 RD2 state."""
import csv, gzip, importlib.util, json, os, shutil, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
ASSET=ROOT/'docs/validation/fixed_ldau/rs/C5_2_PART_D_assets/ni_u3j07_rd2_replay/replay'
PRIVATE=Path('/private/tmp/c52-rd2-private')
BASE=ASSET/'prereg_005/driver.input.txt'
PRE=ASSET/'prereg_005/preregistration.json'
BIND=PRIVATE/'vin_binding_nr501.json'; DRIVER=PRIVATE/'nr501_adapter/build_proven/accepted_m03_replay_driver_nr501'
CAP=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run')
CAP_EXTRACT=Path('/private/tmp/c54-rd2-capture-t1')
OUT=Path('/private/tmp/c54-rd2-scalar-scan-010')
PREF=ROOT/'docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py'
REPORT=ROOT/'docs/validation/fixed_ldau/rs/C5_4_T1_RD2.json'
V=(-.02,-.01,-.005,-.0025,.0025,.005,.01,.02)
THREADS={'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','RSLMTO_OMP_THREADS_SERIAL':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def sha(p):
 import hashlib
 return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 sys.path.insert(0,str(ROOT/'tests/validation/fixed_ldau'))
 adapter=load('rd2_nr501_adapter',ASSET/'nr501_adapter/rd2_replay_adapter_nr501.py')
 capture=CAP
 if not (CAP/'records.csv').exists():
  if not CAP_EXTRACT.exists():
   CAP_EXTRACT.mkdir()
   with gzip.open(CAP/'records.csv.gz','rb') as src, (CAP_EXTRACT/'records.csv').open('wb') as dst: shutil.copyfileobj(src,dst)
   shutil.copy2(CAP/'c2_hubbard_audit.dat',CAP_EXTRACT/'c2_hubbard_audit.dat')
  capture=CAP_EXTRACT
 runner,rows,nr,ncore,nval,vin,_,binding=adapter.load_bound_rows(capture,BIND)
 pre=json.loads(PRE.read_text()); old=BASE.read_text().splitlines(); h=old[0].split()
 if sha(BASE)!=pre['case_driver_input_sha256'] or h[4]!='5': raise SystemExit('frozen base payload mismatch')
 start=3+int(h[1])+2+6; block=19; ql=old[start+4*block+1:start+5*block]
 if len(ql)!=18: raise SystemExit(f'expected 18 frozen ql values, got {len(ql)}')
 newh=h.copy(); newh[4]='8'; lines=[' '.join(newh)]+old[1:start]
 for v in V: lines += [f'{v:.17g} {v:.17g} .true.']+ql
 payload=Path('/private/tmp/c54-rd2-scalar-scan-010.input.txt')
 if payload.exists() or OUT.exists(): raise SystemExit('refusing to overwrite a prior T1 run')
 payload.write_text('\n'.join(lines)+'\n')
 env=os.environ.copy(); env.update(THREADS)
 subprocess.run([sys.executable,str(PREF),str(OUT),'--offline',str(payload),str(BASE),str(DRIVER),str(BIND)],check=True,env=env)
 result=subprocess.run([str(DRIVER),str(OUT/'driver.input.txt'),str(OUT/'driver.output.csv')],cwd=OUT,env=env,capture_output=True,text=True)
 (OUT/'stdout.log').write_text(result.stdout); (OUT/'stderr.log').write_text(result.stderr)
 if result.returncode: raise SystemExit(f'frozen radial driver failed: {result.stderr[-1000:]}')
 cases=runner.output_cases((OUT/'driver.output.csv').read_bytes()); q0=pre['target_operands']['q0']
 prior=json.loads((ASSET/'cases_005/rd2_cases_report.json').read_text())
 identical=json.loads((ASSET/'diagnostic_006/replay_006/replay_result.json').read_text())
 if not identical.get('output_byte_identical') or identical.get('replay_output_sha256')!=prior['driver_output_sha256']:
  raise SystemExit('sign-control reuse lacks byte-identical same-state replay evidence')
 raw=prior['measurements']['cases']['RAW_PL']['sumev_by_spin_ry']
 fit=load('a2_a11',ROOT/'tests/validation/fixed_ldau/m0_a11_scaling_oracle.py')
 data={}; coeff={}; residual={}
 for si,sp in enumerate((1,2)):
  data[str(sp)]={}
  for i,v in enumerate(V,1):
   delta=float(raw[str(sp)])-cases[i]['sumev'][sp]; r=delta-v*float(q0[si])
   data[str(sp)][str(v)]={'delta_sumev_ry':delta,'r_ry':r}; residual.setdefault(sp,{})[v]=r
  coeff[str(sp)]=[float(x) for x in fit._fit(fit.np.asarray(V),fit.np.asarray([residual[sp][v] for v in V]),3)]
 scalar=.01; d0=float(json.loads((ROOT/'docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_preregistration_v3.json').read_text())['accepted_no_u_state']['D0_RS_ry'])
 rfit={sp:sum(c*scalar**j for j,c in enumerate(coeff[str(sp)])) for sp in (1,2)}
 # The +0.01 response is taken from the existing exact same-state five-case replay.
 scalar_dev=float(prior['measurements']['scalar_plus_001']['absolute_deviation_from_v_sum_q0_ry'])
 allowance=sum(abs(x) for x in rfit.values())+d0
 car=prior['measurements']['carrier']; carrier_error=float(car['absolute_error_ry']); threshold=10*carrier_error
 target=pre['target_operands']['carrier_trace_vin_p_ry']; swap=pre['target_operands']['spin_swap_carrier_ry']
 predicted={'SIGN_WRONG':abs(sum(-x for x in target)-sum(target)),
            'OPPOSITE_SPIN':abs(sum(swap)-sum(target))}
 measured_controls={k:float(v['absolute_separation_from_true_carrier_ry']) for k,v in prior['measurements']['wrong_controls'].items()}
 controls={k:{'predicted_separation_ry':v,'measured_separation_ry':measured_controls[k],
   'threshold_10x_carrier_error_ry':threshold,'classification':'NULL' if v<threshold else ('ADMISSIBLE_PASS' if measured_controls[k]>threshold else 'ADMISSIBLE_FAIL')}
   for k,v in predicted.items()}
 ratios=prior['measurements']['carrier']['signed_ratio_by_spin']
 ratios_pass=all(.99<=float(ratios[str(sp)])<=1.01 for sp in (1,2))
 q2pass=min(pre['target_operands']['q2_prime_by_spin_ry2'])>=0
 report={'schema':'c5-4-rd2-ni-scalar-scan-v1','status':'PASS' if ratios_pass and q2pass and scalar_dev<=allowance and all(x['classification']!='ADMISSIBLE_FAIL' for x in controls.values()) else 'FAIL',
  'offline_only':True,'scf_solver_invocations':0,'preflight_path':str(PREF),'preflight_sha256':sha(PREF),
  'base_payload_sha256':sha(BASE),'scan_payload_sha256':sha(OUT/'driver.input.txt'),
  'driver_sha256':sha(DRIVER),'binding_sha256':sha(BIND),'scan_driver_output_sha256':sha(OUT/'driver.output.csv'),
  'carrier_ratio_per_spin':prior['measurements']['carrier']['signed_ratio_by_spin'],'carrier_error_ry':carrier_error,
  'carrier_ratio_pass':ratios_pass,'q2prime_nonnegative_pass':q2pass,
  'q2prime_per_spin_ry2':pre['target_operands']['q2_prime_by_spin_ry2'],
  'scalar_fit':{'method':'A2.2 cubic least-squares r(v)=Σa_j v^j on these eight signed points, per Ni state/spin','v_ry':list(V),
    'data':data,'coefficients_ascending_by_power':coeff,'scalar_control_v_ry':scalar,
    'fitted_residual_by_spin_ry':rfit,'measured_deviation_from_existing_same_state_plus001_ry':scalar_dev,
    'D0_ry':d0,'allowance_ry':allowance,'pass':scalar_dev<=allowance},
  'controls':controls,'sign_control_reuse':{'case':'SIGN_WRONG in frozen five-case replay','case_report_sha256':sha(ASSET/'cases_005/rd2_cases_report.json'),
    'identical_replay_report_sha256':sha(ASSET/'diagnostic_006/replay_006/rd2_cases_report.json'),
    'identical_replay_status':identical['status'],'byte_identical':identical['output_byte_identical']},
  'historical_RD2_FAIL_superseded':True}
 REPORT.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':report['status'],'report':str(REPORT),'scf_solver_invocations':0,
   'scalar_allowance_ry':allowance,'scalar_deviation_ry':scalar_dev,'controls':controls},indent=2,sort_keys=True))

if __name__=='__main__': main()
