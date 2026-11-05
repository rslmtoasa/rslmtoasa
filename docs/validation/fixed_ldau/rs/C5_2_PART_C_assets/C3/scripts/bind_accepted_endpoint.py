#!/usr/bin/env python3
"""Bind a future C3 capture directory to a source-schema-checked C2 endpoint lock."""
from __future__ import annotations
import argparse, hashlib, json, shutil, sys, os
from pathlib import Path
sys.path.insert(0,'/private/tmp/c52-c-private')
import c3_c2_adapter as c3
from fortran_namelist_schema import source_schema,validate_input
import f90nml
ROOT=Path('/private/tmp/c52-c2-stage/c3_debug_capture_template')
REPO=Path('/Users/andersb/Jobb/rslmto_fable_v4b')
BINARY=Path('/private/tmp/c52-c2-stage/crash_triage_preflight/frozen/debug_bounds/b8b00eff0f1a4c3b0536006f7593c53e229cac7f0c77cbfed29126125af26b69/rslmto.x')
SOURCE=REPO/'source/self.f90'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--lock',type=Path,required=True,help='audited C2 accepted Chebyshev U2 local endpoint lock')
 ap.add_argument('--capture-dir',type=Path,required=True,help='fresh dedicated output directory; no default path')
 a=ap.parse_args(); lockp=a.lock.resolve(); out=a.capture_dir.resolve()
 if not lockp.is_file(): raise SystemExit('missing accepted-endpoint lock')
 if out.exists() and any(out.iterdir()): raise SystemExit('capture directory must be fresh and empty')
 out.mkdir(parents=True,exist_ok=True)
 lock=c3.load_json(lockp); binding=c3.validate_c2_lock(lock,lockp)
 if binding['executable']['binary_sha256']!=sha(BINARY): raise SystemExit('accepted executable differs from frozen C3 capture executable')
 if binding['executable']['source_sha256']!=sha(SOURCE): raise SystemExit('accepted source differs from current frozen source')
 endpoint=binding['endpoint']; checkpoint=Path(endpoint['checkpoint_path']).resolve()
 if endpoint['checkpoint_sha256']!=sha(checkpoint): raise SystemExit('accepted checkpoint hash changed')
 schema=source_schema(REPO); validate_input(ROOT/'input.nml',schema)
 deck=f90nml.read(str(ROOT/'input.nml'))
 if deck['self'].get('use_kspace') is not False: raise SystemExit('capture deck must set self.use_kspace=.false.')
 if int(deck['self']['nstep'])!=1: raise SystemExit('capture deck must be one SCF step')
 if deck['control']['recur'].lower() not in {'cheb','chebyshev'}: raise SystemExit('capture must use Chebyshev recurrence')
 immutable=out/'immutable_inputs'; immutable.mkdir(parents=True,exist_ok=True)
 immutable_seed=immutable/'Fe_out.nml'; shutil.copyfile(checkpoint,immutable_seed); os.chmod(immutable_seed,0o444)
 immutable_input=immutable/'input.nml'; shutil.copyfile(ROOT/'input.nml',immutable_input); os.chmod(immutable_input,0o444)
 seed=out/'Fe_out.nml'; shutil.copyfile(immutable_seed,seed); os.chmod(seed,0o644)
 working_input=out/'input.nml'; shutil.copyfile(immutable_input,working_input); os.chmod(working_input,0o644)
 if sha(seed)!=sha(checkpoint) or sha(immutable_seed)!=sha(checkpoint): raise SystemExit('staged endpoint copy changed bytes')
 if not os.access(seed,os.W_OK) or not os.access(working_input,os.W_OK): raise SystemExit('working inputs must be writable')
 result={'schema':'c52-c3-capture-endpoint-binding-v1','status':'BOUND_NOT_LAUNCHED',
  'lock_path':str(lockp),'lock_sha256':sha(lockp),'run_id':endpoint['run_id'],
  'accepted_iteration':endpoint['accepted_iteration'],'accepted_checkpoint_path':str(checkpoint),
  'accepted_checkpoint_sha256':sha(checkpoint),'staged_Fe_out_path':str(seed),'staged_Fe_out_sha256':sha(seed),'immutable_Fe_out_path':str(immutable_seed),'immutable_Fe_out_sha256':sha(immutable_seed),'staged_input_path':str(working_input),'staged_input_sha256':sha(working_input),'immutable_input_path':str(immutable_input),'immutable_input_sha256':sha(immutable_input),'working_inputs_writable':True,
  'accepted_deck_path':endpoint['deck_path'],'accepted_deck_sha256':endpoint['deck_sha256'],
  'block_chebyshev_comparison_status':binding['acceptance']['block_chebyshev_comparison_status'],
  'capture_template_deck_path':str(ROOT/'input.nml'),'capture_template_deck_sha256':sha(ROOT/'input.nml'),
  'capture_deck_schema_valid':True,'binary_path':str(BINARY),'binary_sha256':sha(BINARY),
  'source_path':str(SOURCE),'source_sha256':sha(SOURCE),'driver_path':'/private/tmp/c52-c-private/helper-build/accepted_m03_replay_driver',
  'driver_sha256':sha('/private/tmp/c52-c-private/helper-build/accepted_m03_replay_driver'),
  'note':'This operation binds a future capture to a C2 endpoint; it does not launch the solver or establish a C3 result.'}
 (out/'capture_binding.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
