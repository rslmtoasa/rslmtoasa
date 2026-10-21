import hashlib,json,os,platform
from pathlib import Path
B=Path('/private/tmp/c52-rd2-private'); D=B/'diagnostic_006'; CAP=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run')
paths={
 'capture_csv':CAP/'records.csv','local_hall':CAP/'c2_hubbard_audit.dat',
 'local_hall_parser':Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/c2_hubbard_operator_audit.py'),
 'production_moment_source':Path('/Users/andersb/Jobb/rslmto_fable_v4b/source/math.f90'),
 'binding':B/'vin_binding_nr501.json','plan':B/'rd2_replay_plan.json',
 'adapter':B/'nr501_adapter/rd2_replay_adapter_nr501.py',
 'derived_runner':B/'nr501_adapter/rs_a3_m03_runner_nr501.py',
 'driver':B/'nr501_adapter/build_proven/accepted_m03_replay_driver_nr501',
 'baseline_report':B/'baseline_003/baseline.json','prereg_report':B/'prereg_005/preregistration.json',
 'case_input':B/'prereg_005/driver.input.txt','first_output':B/'cases_005/driver.output.csv',
 'first_report':B/'cases_005/rd2_cases_report.json','first_authorization':B/'cases_005.authorization.json',
 'payload_diagnostic_script':D/'audit_payload.py','payload_diagnostic_report':D/'payload_phase_spin_audit.json',
 'replay_script':D/'run_identical_replay.py'}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
comps={k:{'path':str(v),'sha256':sha(v),'bytes':v.stat().st_size} for k,v in paths.items()}
pre={
 'schema':'c5-rd2-identical-debug-replay-preflight-v1',
 'status':'PREPARED_FOR_STATIC_AUDIT; no replay executed',
 'hypothesis':'A remaining payload/phase/spin binding issue, rather than the measured Ni response, may explain the two formal gate failures; audit the frozen case payload and repeat it exactly once.',
 'attempt':{'run_id':'rd2_identical_debug_replay_006','one_replay_only':True,'new_output_dir':str(D/'replay_006'),'same_case_input_path':str(B/'prereg_005/driver.input.txt'),'same_case_input_sha256':comps['case_input']['sha256'],'same_frozen_driver_path':str(paths['driver']),'same_frozen_driver_sha256':comps['driver']['sha256'],'no_scf_solver':True},
 'diagnostic':{'script_path':str(paths['payload_diagnostic_script']),'script_sha256':comps['payload_diagnostic_script']['sha256'],'report_path':str(paths['payload_diagnostic_report']),'report_sha256':comps['payload_diagnostic_report']['sha256'],'checks':'five decoded case records and 501-point input payload; P_OUT and M1_OUT matrices versus all 1,567 RS Green blocks for P/M1/M2; exact LOCAL_HALL Vin versus v_hub type2; true/opposite-spin Tr(VP); P_OUT, m0, and saved radial QL q0; same-phase PL gc and q2 payload.'},
 'diagnostic_summary':{'Vin_up_minus_down_frobenius':9.932164440462364e-09,'Vin_up_minus_down_max_entry':5.330290642113009e-09,'relative_Vin_difference':5.86373438259428e-08,'true_total_TrVP_Ry':-0.7178940285707095,'swapped_total_TrVP_Ry':-0.7178940285707092,'total_swap_delta_Ry':3.3306690738754696e-16,'q0_Pout_minus_saved_QL_abs_by_spin':[2.412443578236889e-10,4.233676875742276e-08]},
 'formal_first_result':{'report_path':str(paths['first_report']),'report_sha256':comps['first_report']['sha256'],'raw_csv_sha256':comps['first_output']['sha256'],'status':'GATE_MEASURED_FAIL','formal_causes':['scalar residual exceeds registered A_scalar','OPPOSITE_SPIN separation below 10x carrier error'],'legacy_diagnostics_not_gates':True},
 'commands':{
  'entrypoint':[str(Path('/opt/homebrew/opt/python@3.11/bin/python3.11')),
                str(paths['replay_script'])],
  'actual_adapter_argv':[str(Path('/opt/homebrew/opt/python@3.11/bin/python3.11')),str(paths['adapter']),'--mode','cases','--capture-dir',str(CAP),'--binding',str(paths['binding']),'--driver',str(paths['driver']),'--authorization',str(D/'replay_authorization.json'),'--baseline-report',str(paths['baseline_report']),'--prereg-report',str(paths['prereg_report']),'--output-dir',str(D/'replay_006')],
  'cwd':'/Users/andersb/Jobb/rslmto_fable_v4b',
  'effective_environment':{'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','RSLMTO_OMP_THREADS_SERIAL':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'},
  'python_version':'3.11.17', 'process_model':'Adapter invokes frozen helper exactly once with argv [driver, immutable prereg_005/driver.input.txt, new replay_006/driver.output.csv].'
 },
 'authorization':{'required_path':str(D/'replay_authorization.json'),'required_mode':'cases','required_scope':'one exact five-case reproducibility replay; no SCF solver','must_bind_preflight_sha256':True,'not_created_by_worker':True},
 'components':comps}
(D/'replay_preflight.json').write_text(json.dumps(pre,indent=2,sort_keys=True)+'\n')
print(json.dumps({'path':str(D/'replay_preflight.json'),'sha256':sha(D/'replay_preflight.json'),'components':len(comps)},indent=2))
