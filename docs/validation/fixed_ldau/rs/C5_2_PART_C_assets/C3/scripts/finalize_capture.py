#!/usr/bin/env python3
"""Validate and hash a future one-step C3 full capture; never launches a solver."""
import argparse,csv,hashlib,json,math,os,sys
from pathlib import Path
sys.path.insert(0,'/private/tmp/c52-c-private')
import c3_c2_adapter as c3

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--capture-dir',type=Path,required=True);a=ap.parse_args()
 d=a.capture_dir.resolve(); binding=json.loads((d/'capture_binding.json').read_text())
 rc=int((d/'exit_code.txt').read_text().strip())
 if rc!=0: raise SystemExit(f'solver exit code is {rc}, capture not accepted as structurally complete')
 pid=int((d/'solver.pid').read_text().strip())
 try: os.kill(pid,0)
 except ProcessLookupError: pass
 except PermissionError: raise SystemExit('cannot verify solver PID is closed')
 else: raise SystemExit(f'solver PID {pid} still exists')
 csvp=d/'records.csv'
 with csvp.open(newline='') as f:
  reader=csv.DictReader(f); header=tuple(reader.fieldnames or ()); rows=list(reader)
 if header!=c3.rb3.CSV_FIELDS: raise SystemExit('unexpected C3 trace CSV header')
 if not rows: raise SystemExit('empty capture CSV')
 runs={r['run_id'] for r in rows}; iters={int(r['iteration']) for r in rows}
 wanted='c3_fe_u2_cheb_fullcapture_i1'
 if runs!={wanted}: raise SystemExit(f'observed run ID {runs} differs from requested {wanted}; reject capture')
 if iters!={1}: raise SystemExit(f'capture must contain selected iteration 1 only, observed {iters}')
 for i,row in enumerate(rows,1):
  for key in ('lambda','weight','energy','scalar','re','im'):
   if not math.isfinite(float(row[key])): raise SystemExit(f'nonfinite {key} at CSV row {i}')
 kinds={r['kind'] for r in rows}
 required=set(c3.REQUIRED_KINDS)|{'radial_replay_rho_out','radial_replay_sumec','radial_replay_sumev',
  'radial_replay_rho_core','radial_replay_wave','radial_replay_nre','radial_replay_eigenvalue',
  'radial_replay_pnu','radial_replay_enu','radial_replay_gravity_center','radial_replay_band_center',
  'radial_replay_c','radial_replay_qpar','radial_replay_srdel','m0_p_availability',
  'm0_carrier_identity','m0_trace_consistency','m1_trace_consistency'}
 missing=sorted(required-kinds)
 if missing: raise SystemExit(f'full writer contract incomplete; missing kinds {missing}')
 shape=c3.validate_full_green_rows(rows)
 lock_binding=binding
 manifest={'schema':c3.TRACE_SCHEMA,'status':'SELECTED_ITERATION_CAPTURED','backend':'RS',
  'run_id':wanted,'iteration':1,'lambda':1.0,'nsp':1,'lmax':2,'norb_per_spin':9,
  'source_sha':binding['source_sha256'],'source_revision':binding.get('source_revision'),
  'executable_sha256':binding['binary_sha256'],'csv_sha256':sha(csvp),
  'input_sha256':{'accepted_checkpoint':binding['accepted_checkpoint_sha256'],
                  'capture_input':binding['capture_template_deck_sha256']},
  'capture_execution':{'run_id':wanted,'iteration':1,'input_sha256':binding['capture_template_deck_sha256'],
    'accepted_endpoint_lock_sha256':binding['lock_sha256'],'accepted_checkpoint_sha256':binding['accepted_checkpoint_sha256'],
    'exit_code':rc,'solver_pid':pid,'wrapper_pid':int((d/'wrapper.pid').read_text().strip()),
    'capture_only':True,'selected_phase':'P_OUT at outer iteration 1; radial replay rows from terminal LAST radial pass'},
  'capture_only':True,'ql_only':False,'moments_only':False,'full_green_shape':shape,
  'observed_run_id_verified':True,'required_writer_row_kinds':sorted(required),
  'notice':'Capture integrity only; no C3 scientific gate or Co U>0 admission asserted.'}
 path=d/'manifest.json';path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':manifest['status'],'manifest_path':str(path),'manifest_sha256':sha(path),
  'csv_sha256':sha(csvp),'row_count':len(rows),'full_green_shape':shape},indent=2))
if __name__=='__main__':main()
