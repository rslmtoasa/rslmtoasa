from pathlib import Path
import hashlib,json
base=Path('/private/tmp/c52-rd2-private')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
paths={
'adapter':base/'nr501_adapter/rd2_replay_adapter_nr501.py',
'binding':base/'vin_binding_nr501.json',
'binding_script':base/'rd2_binding_preflight.py',
'capture_csv':Path('/private/tmp/c52-rd2-ni-u3j07-capture/run/records.csv'),
'local_hall':Path('/private/tmp/c52-rd2-ni-u3j07-capture/run/c2_hubbard_audit.dat'),
'driver':base/'nr501_adapter/build_proven/accepted_m03_replay_driver_nr501',
'driver_source':base/'nr501_adapter/build_proven/accepted_m03_replay_driver_nr501.f90'
,'original_driver_source':Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/accepted_m03_replay_driver.f90')
,'derived_runner':base/'nr501_adapter/rs_a3_m03_runner_nr501.py'
,'variant_manifest':base/'nr501_adapter/variant_manifest.json'
,'driver_build_log':base/'nr501_adapter/build_proven/build.stderr.txt'
,'driver_build_stdout':base/'nr501_adapter/build_proven/build.stdout.txt'
,'driver_build_invocation':base/'nr501_adapter/build_proven/build_invocation.json'
,'builder':base/'write_plan_phasefix.py',
'driver_manifest':Path('/private/tmp/c52-c-private/helper-build/source_build_manifest.json'),
'driver_build_record':Path('/private/tmp/c52-c-private/helper-build/helper_build_record.json'),
'runner':base/'nr501_adapter/rs_a3_m03_runner_nr501.py'
,'original_runner':Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/rs_a3_m03_runner.py'),
'c2_parser':Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/c2_hubbard_operator_audit.py'),
'local_formula_helper':Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/local_moment_capture_oracle.py'),
'no_driver_smoke':base/'oracle_helper_smoke.py',
'attempt_001':base/'baseline_attempt_001.json',
'attempt_001_stderr':base/'baseline_attempt_001.stderr.txt',
'attempt_001_authorization':base/'baseline_attempt_001.authorization.json',
'attempt_002':base/'baseline_002.exit_code.txt',
'attempt_002_stderr':base/'baseline_002.stderr.txt',
'attempt_002_stdout':base/'baseline_002.stdout.txt',
'attempt_002_authorization':base/'baseline_attempt_002.authorization.json',
'fit':Path('/Users/andersb/Jobb/rslmto_fable_v4b/docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_fit_results.json'),
'fit_prereg':Path('/Users/andersb/Jobb/rslmto_fable_v4b/docs/validation/fixed_ldau/evidence/a3_rs_a22_scalar_preregistration_v3.json'), 'baseline_report':base/'baseline_003/baseline.json', 'prereg_attempt_003_stderr':base/'prereg_003.stderr.txt', 'prereg_attempt_003_exit':base/'prereg_003.exit_code.txt', 'prereg_attempt_003_auth':base/'prereg_003.authorization.json', 'phase_proof':base/'green_phase_proof_002.json'}
for k,p in list(paths.items()):
 if not p.is_file(): raise SystemExit(f'missing {k}: {p}')
 paths[k]={'path':str(p),'sha256':sha(p)}
script='/opt/homebrew/opt/python@3.11/bin/python3.11'
common=['--capture-dir','/private/tmp/c52-rd2-ni-u3j07-capture/run','--binding',str(paths['binding']['path']),'--driver',str(paths['driver']['path']),'--authorization','/private/tmp/c52-rd2-private/execution_authorization.json']
cmds={
'baseline':[script,str(paths['adapter']['path']),'--mode','baseline',*common,'--output-dir','/private/tmp/c52-rd2-private/baseline_003'],
'preregister':[script,str(paths['adapter']['path']),'--mode','preregister',*common,'--baseline-report','/private/tmp/c52-rd2-private/baseline_003/baseline.json','--output-dir','/private/tmp/c52-rd2-private/prereg_005'],
'cases':[script,str(paths['adapter']['path']),'--mode','cases',*common,'--baseline-report','/private/tmp/c52-rd2-private/baseline_003/baseline.json','--prereg-report','/private/tmp/c52-rd2-private/prereg_005/preregistration.json','--output-dir','/private/tmp/c52-rd2-private/cases_005']}
plan={
 'schema':'c5-rd2-ni-frozen-replay-plan-v1',
 'status':'BASELINE_003_PASS; PREREG_003_REJECTED_WRITER_TYPE_PHASE_BINDING_BEFORE_CASE_DRIVER; PREREG_004_PREPARED_WITH_RAW_GREEN_MOMENTS_PENDING_AUDIT',
 'repository':{'path':'/Users/andersb/Jobb/rslmto_fable_v4b','branch':'fix/ldau-rs-impurity','head':'683eff937241a463d2dacb395485d6b7c0b43264'},
 'capture_binding':{
  'capture_status':'RD2_CAPTURE_INTEGRITY_PASS; no RD2 numerical replay/gate executed',
  'run_id':'rd2_ni_u3j07_endpoint23_fullgreen_i1','iteration':1,
  'accepted_impurity_checkpoint_iteration':23,
  'accepted_state_sha256':'7ea1b7b8bd6e06be39816251a1de218a9009c39f2ad73d8fd25c52655476e60f',
  'accepted_trace_sha256':'247fe7012d829eefa8f9765688dc27e5abb972c9591029a9a43d697f679133aa',
  'capture_csv':paths['capture_csv'],'local_hall':paths['local_hall'],
  'radial_replay_dimensions':{'NR':501,'lmax':2,'core_seed_count':10,'valence_seed_count':6,'evidence':'unique RADIAL_REPLAY rows in records.csv; prior adapter expectation NR=495 was a preflight-only defect'},
  'binding':paths['binding'],
  'runtime_map':'recursion ordinal 1 -> irec(1)=1 -> cluster nlim/site 1 -> iz(1)=type 2 (Ni-imp); sites 2..19 are frozen host type 1',
  'Vin':'Exact additive LOCAL_HALL BLOCK correction (after-before), site=1/type=2, 50 complex entries; independently equals captured v_hub type=2 Pin_H d-block exactly. Do not use type-1 host selector.',
  'capture_erratum':{'path':'/private/tmp/c52-rd2-ni-u3j07-capture/analysis/runtime_mapping_erratum.json','sha256':'2357bbfa41cc429063904191a1f01aa8e1c9df438b6e7d1c56211cc2cf8b702c'},
  'original_preflight_preserved':{'path':'/private/tmp/c52-rd2-ni-u3j07-capture/provenance/preflight.json','sha256':'75b944fa34fb896ff365c0c984b7d328ec0ca6498463d5af950a2b5425aaaf1f','erratum':'Original expected nlim2 is superseded by append-only runtime mapping erratum; no raw record rewritten.'}},
 'components':paths,
 'nr501_harness_variant':paths['variant_manifest'],
 'attempt_history':{'attempt_001':'adapter rejected NR=501 before driver invocation; zero helper/SCF runs; stderr and auth preserved', 'attempt_002':'baseline rejected before helper invocation due stale driver binding; preserved separately', 'baseline_003':'frozen helper baseline PASS; report SHA dd6eb6919ae1673459afe8f636289b9fc010b02bdc11dbdfc34745f6ded0153f', 'prereg_003':'prereg calculations rejected before response driver because m0 type is writer metadata 19, not physical type2; no case helper invoked'},
 'formula':{
  'carrier_t_s_ry':'Re Tr(Vin_s P_s)',
  'carrier_shift_x_s_ry':'t_s / q0_s',
  'corrected_center_cprime_s_ry':'c_s - x_s',
  'corrected_q2_q2prime_s_ry2':'q2_s - 2 Re Tr(Vin_s (M1_s - c_s P_s)) + Re Tr(Vin_s^2 P_s) - t_s^2/q0_s',
  'Vin_source':'Captured exact additive LOCAL_HALL insertion at recursive site1/type2. Historical generic RB3 spin-fit keys receive per-spin coefficients only; they do not select Vin/site/type.', 'moment_phase':'P_OUT/m0 is post-bands local correction; prereg raw c comes from Tr(M1)/q0 cross-checked with full RS_GREEN p2/p1, and raw q2 from p3-2*c*p2+c^2*p1. Apply Vin correction exactly once; RAW replay center is reconstructed c_raw, LOCAL shifts c_raw by t/q0, and RAW/LOCAL/SCALAR share corrected q2prime.'},
 'locked_fit_control':{
  'fit':paths['fit'],'fit_preregistration':paths['fit_prereg'],
  'D0_RS_ry':5.056500278044496e-06,
  'fit_response':'per-spin cubic residual r_s(x) from locked A3 RS A2.2 fit; registered D0 is added to absolute fit-residual sums for allowances',
  'fit_mapping':'coefficients copied by spin into generic compatibility keys site=0,type=0,spin=1/2 only; no type/site inference'},
 'registered_gate_and_cases':{
  'carrier_ratio':'per-spin measured local replay carrier / exact captured true t_s, inclusive [0.99,1.01]',
  'q2prime':'both spin q2prime >= 0',
  'scalar_control':'absolute scalar +0.01 replay response deviation from registered ideal V*q0 response <= locked scalar fit allowance (RB3 formula)',
  'wrong_controls':'SIGN_WRONG and OPPOSITE_SPIN absolute separation from measured true carrier each > 10 × measured carrier error',
  'cases':[
   {'name':'RAW_PL','shift_ry':'[0,0]','q2':'exact-V corrected q2prime, matching accepted runner.prepare_rs construction'},
   {'name':'LOCAL_PL','shift_ry':'[x_up,x_down] from exact Vin','q2':'exact corrected formula above'},
   {'name':'SCALAR_V_PLUS_001','shift_ry':'[0.01,0.01]','q2':'same exact-V corrected q2prime as RAW_PL/LOCAL_PL'},
   {'name':'SIGN_WRONG','shift_ry':'[-x_up,-x_down]','q2':'same corrected q2 as LOCAL_PL'},
   {'name':'OPPOSITE_SPIN','shift_ry':'[t_down/q0_up,t_up/q0_down]','q2':'corrected formula using opposite-spin Vin'}]},
 'execution_environment':{'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','RSLMTO_OMP_THREADS_SERIAL':'1',
  'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
  'purpose':'single-thread deterministic frozen-helper execution; no SCF solver process'},
 'command_sequence':{
  'baseline':{'argv':cmds['baseline'],'expected_report':'/private/tmp/c52-rd2-private/baseline_003/baseline.json','purpose':'frozen radial helper baseline against captured raw-state outputs; this is a numerical replay and requires separately hash-bound authorization'},
  'preregister':{'argv':cmds['preregister'],'expected_outputs':['/private/tmp/c52-rd2-private/prereg_005/preregistration.json','/private/tmp/c52-rd2-private/prereg_005/driver.input.txt','/private/tmp/c52-rd2-private/prereg_005/rb3_fit_shape_adapter.json'],'purpose':'derive exact operands/case payloads and fit allowances only after baseline PASS and separate authorization'},
  'cases':{'argv':cmds['cases'],'expected_outputs':['/private/tmp/c52-rd2-private/cases_005/driver.output.csv','/private/tmp/c52-rd2-private/cases_005/driver.stdout.log','/private/tmp/c52-rd2-private/cases_005/driver.stderr.log'],'purpose':'five frozen radial helper cases only; not an SCF solver run'},
  'authorization_contract':'Each mode requires exact-hash authorization. The current execution_authorization.json is present and scoped to preregister only, bound to the passing baseline. No cases-mode authorization exists; response-driver execution is not authorized.'},
 'case_phase_contract':'Following accepted runner.prepare_rs, all cases share exact-V corrected q2prime. Reconstruct raw c and q2 from RS_GREEN/M1 before correction. RAW_PL uses the raw center with vmad removed via same-phase CALCULATE_PL gc plus (raw c−captured local m0 c); LOCAL_PL shifts by x=t/q0 once and returns the captured local gc. SCALAR_V_PLUS_001 shifts by +0.01, SIGN_WRONG applies -x, and OPPOSITE_SPIN applies the opposite-spin shift/q2. Captured P_OUT/m0 are post-correction cross-checks, not raw operands. m0 type=19 is lattice-ib writer metadata; physical Vin type2 comes from runtime map.',
 'execution_status':{'baseline_attempts':3,'baseline_helper_runs':1,'baseline_attempt_001':'PREFLIGHT_REJECTED_BEFORE_DRIVER because adapter expected NR=495, capture NR=501; exact traceback and stale authorization preserved separately','baseline_attempt_002':'PREFLIGHT_REJECTED_BEFORE_DRIVER due driver SHA mismatch against original binding; no helper invocation, exact evidence preserved','baseline_attempt_003':'PASS; helper baseline report is frozen and separately audited', 'prereg_attempt_003':'FAILED_CLOSED_BEFORE_RESPONSE_DRIVER: obsolete m0 type assertion; no case helper invocation','preregistration_calculations':2,'prereg_004':'successful preparation preserved under prior plan/auth; output retained as historical, no response driver invoked','prereg_005':'pending execution under current plan-bound prereg-only authorization', 'case_replays':0,'SCF_solver_launches':0,'active_authorization_file_exists':(base/'execution_authorization.json').exists(),'active_authorization_mode':'preregister','cases_authorization_present':False},
 'limitations':['Capture is a one-iteration full-Green integrity capture, not itself a response measurement.','No baseline/helper replay, response preregistration arithmetic, case replay, or RD2 gate result has been executed in this preparation.','The private adapter fails closed behind a separate exact-hash authorization file; this plan grants no execution authorization.']}
out=base/'rd2_replay_plan.json'
out.write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
print(str(out),sha(out))
print('adapter',sha(paths['adapter']['path']))
