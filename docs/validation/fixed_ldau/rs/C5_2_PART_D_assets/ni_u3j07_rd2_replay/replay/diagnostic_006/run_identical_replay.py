#!/usr/bin/env python3
"""One hash-locked identical RD2 response replay; never launches an SCF solver."""
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path

ROOT=Path('/private/tmp/c52-rd2-private')
PREFLIGHT=ROOT/'diagnostic_006/replay_preflight.json'
AUTH=ROOT/'diagnostic_006/replay_authorization.json'
OUT=ROOT/'diagnostic_006/replay_006'
REPO=Path('/Users/andersb/Jobb/rslmto_fable_v4b')
PYTHON=Path('/opt/homebrew/opt/python@3.11/bin/python3.11')
ADAPTER=ROOT/'nr501_adapter/rd2_replay_adapter_nr501.py'
CAPTURE=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run')
BINDING=ROOT/'vin_binding_nr501.json'
DRIVER=ROOT/'nr501_adapter/build_proven/accepted_m03_replay_driver_nr501'
BASELINE=ROOT/'baseline_003/baseline.json'
PREREG=ROOT/'prereg_005/preregistration.json'
FIRST=ROOT/'cases_005'
INPUT=ROOT/'prereg_005/driver.input.txt'


def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def readj(p): return json.loads(Path(p).read_text())
def need(ok,msg):
 if not ok: raise RuntimeError(msg)

def validate(preflight_only=False):
 p=readj(PREFLIGHT)
 need(p['schema']=='c5-rd2-identical-debug-replay-preflight-v1','preflight schema mismatch')
 for name,item in p['components'].items():
  path=Path(item['path'])
  need(path.is_file(),f'missing component {name}: {path}')
  need(sha(path)==item['sha256'],f'component hash mismatch {name}')
 need(not OUT.exists(),'replay output path already exists; refusing repeat')
 if not preflight_only:
  a=readj(AUTH)
  need(a.get('status')=='ROOT_AND_AUDITOR_CLEARED','debug replay authorization not cleared')
  need(a.get('mode')=='cases','adapter requires exact cases-mode authorization')
  need(a.get('scope')=='one exact five-case reproducibility replay; no SCF solver','authorization scope mismatch')
  need(a.get('preflight_sha256')==sha(PREFLIGHT),'authorization does not bind preflight')
  for key,path in [('plan_path',ROOT/'rd2_replay_plan.json'),('adapter_path',ADAPTER),
                   ('binding_path',BINDING),('driver_path',DRIVER),
                   ('baseline_report_path',BASELINE),('prereg_report_path',PREREG),
                   ('driver_input_path',INPUT)]:
   need(Path(a[key]).resolve()==path.resolve(),f'authorization path mismatch {key}')
  for key,path in [('plan_sha256',ROOT/'rd2_replay_plan.json'),('adapter_sha256',ADAPTER),
                   ('binding_sha256',BINDING),('driver_sha256',DRIVER),
                   ('baseline_report_sha256',BASELINE),('prereg_report_sha256',PREREG),
                   ('driver_input_sha256',INPUT)]:
   need(a[key]==sha(path),f'authorization hash mismatch {key}')
  need(Path(a['output_dir']).resolve()==OUT.resolve(),'authorization output directory mismatch')
  need(a['first_output_sha256']==sha(FIRST/'driver.output.csv'),'authorization first-output hash mismatch')
 return p

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--preflight-only',action='store_true')
 args=ap.parse_args()
 pref=validate(args.preflight_only)
 if args.preflight_only:
  print(json.dumps({'status':'PREFLIGHT_COMPONENTS_VALID; NO DRIVER INVOKED',
                    'preflight_sha256':sha(PREFLIGHT),'runner_sha256':sha(Path(__file__)),
                    'output_path_absent':not OUT.exists()},sort_keys=True))
  return
 env=os.environ.copy()
 env.update({'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','RSLMTO_OMP_THREADS_SERIAL':'1',
             'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1',
             'BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'})
 command=[str(PYTHON),str(ADAPTER),'--mode','cases','--capture-dir',str(CAPTURE),
  '--binding',str(BINDING),'--driver',str(DRIVER),'--authorization',str(AUTH),
  '--baseline-report',str(BASELINE),'--prereg-report',str(PREREG),'--output-dir',str(OUT)]
 p=subprocess.run(command,cwd=REPO,env=env,capture_output=True,text=True)
 OUT.mkdir(parents=False,exist_ok=True) if not OUT.exists() else None
 (OUT/'adapter.stdout.log').write_text(p.stdout)
 (OUT/'adapter.stderr.log').write_text(p.stderr)
 (OUT/'adapter.exit_code.txt').write_text(f'{p.returncode}\n')
 # Retain an exact copy proving the input used by the adapter remains immutable.
 (OUT/'driver.input.copy.txt').write_bytes(INPUT.read_bytes())
 need(sha(OUT/'driver.input.copy.txt')==pref['components']['case_input']['sha256'],'copied payload hash mismatch')
 need((OUT/'rd2_cases_report.json').is_file(),'adapter failed before report creation; preserve files and stop')
 replay=readj(OUT/'rd2_cases_report.json'); first=readj(FIRST/'rd2_cases_report.json')
 out_sha=sha(OUT/'driver.output.csv'); first_sha=sha(FIRST/'driver.output.csv')
 same_gates=(replay.get('status')==first.get('status') and
             replay.get('formal_gates')==first.get('formal_gates') and
             replay.get('measurements')==first.get('measurements'))
 summary={'schema':'c5-rd2-identical-debug-replay-result-v1',
  'status':'EXACT_OUTPUT_AND_VERDICT_REPRODUCED' if out_sha==first_sha and same_gates else 'REPLAY_OUTPUT_OR_VERDICT_DIFFERS',
  'one_replay_only':True,'solver_launched':False,'preflight_sha256':sha(PREFLIGHT),
  'authorization_sha256':sha(AUTH),'adapter_exit_code':p.returncode,
  'effective_environment':{k:env[k] for k in ('OMP_NUM_THREADS','OMP_DYNAMIC','RSLMTO_OMP_THREADS_SERIAL','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
  'first_output_sha256':first_sha,'replay_output_sha256':out_sha,'output_byte_identical':out_sha==first_sha,
  'first_formal_status':first.get('status'),'replay_formal_status':replay.get('status'),
  'formal_gates_identical':replay.get('formal_gates')==first.get('formal_gates'),
  'measurements_identical':replay.get('measurements')==first.get('measurements'),
  'first_report_sha256':sha(FIRST/'rd2_cases_report.json'),
  'replay_report_sha256':sha(OUT/'rd2_cases_report.json'),
  'first_payload_sha256':sha(INPUT),'replay_payload_copy_sha256':sha(OUT/'driver.input.copy.txt'),
  'adapter_stdout_sha256':sha(OUT/'adapter.stdout.log'),'adapter_stderr_sha256':sha(OUT/'adapter.stderr.log'),
  'adapter_exit_status_path':str(OUT/'adapter.exit_code.txt')}
 (OUT/'replay_result.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,sort_keys=True))
 if p.returncode: raise SystemExit(p.returncode)

if __name__=='__main__': main()
