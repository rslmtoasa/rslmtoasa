#!/usr/bin/env python3
"""Hash-first assembly of the 24 W1 signed response captures.

This tool consumes only the exact reviewed registration, execution receipts,
raw CSV/stdout, and frozen lambda-zero anchors. It never launches SCF or
classifies the response. It reuses the certified lambda-zero assembler's
identical token/radius formulas, with a memoized float/radius layer whose
synthetic equivalence test is separate.
"""
from __future__ import annotations
import argparse
import builtins
import csv
import functools
import hashlib
import importlib.util
import json
import math
import re
import tempfile
from pathlib import Path

ROOT=Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
ASSEMBLER=ROOT/'tests/validation/fixed_ldau/a4_m21_assemble_lambda0.py'
ASSEMBLER_SHA='d63d10b2c64b93a4f05830f08dd2be905f24c6783c44cb53b8446cc144de6cc0'
SOURCE=ROOT/'source/self.f90'
SOURCE_SHA='2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010'
STRUX_WRITER=ROOT/'source/lattice_strux.f90'
STRUX_WRITER_SHA='50a3c1270ad6b54025e503b99d89f5906a62f1e72b02b830ea83a971b390345d'
STRUX_REPLACE=ROOT/'source/include_codes/strux_lib/strux_structb_compare.f90'
STRUX_REPLACE_SHA='dd80cb8f9ad2448af4fe9403121d74281672b2da2f5d1982bd9abdf6d121b450'
BINARY=Path('/private/tmp/a4-m2-serial/bin/rslmto.x')
BINARY_SHA='536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4'
CALIBRATION_PREREG=ROOT/'docs/validation/fixed_ldau/evidence/a4_m21_ks_lambda0_calibration_preregistration.json'
CALIBRATION_PREREG_SHA='50f24c646af2e7c1eedc050d6b0b3c14b7cb55b5d6134e02cfa081913e1ba80f'
BASE_PREP=Path('/private/tmp/a4-m21-ks-lambda0-prep')
EXPECTED_RUNS={f'w1_{mode}_{sign}{h:.4f}_{threshold}'
  for mode in ('local','spectral') for h in (.0025,.005,.01)
  for sign in ('p','m') for threshold in ('standard','tight')}


def sha(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load_assembler():
    if sha(ASSEMBLER)!=ASSEMBLER_SHA: raise ValueError('certified assembler hash mismatch')
    spec=importlib.util.spec_from_file_location('a4_certified_assembler',ASSEMBLER)
    if not spec or not spec.loader: raise ValueError('cannot load certified assembler')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod._uncached_float=builtins.float
    mod._uncached_token_radius=mod.token_radius
    return mod

def install_exact_caches(mod):
    """Memoize exact existing conversions; no arithmetic/radius formula changes."""
    cached_float=functools.lru_cache(maxsize=None)(mod._uncached_float)
    cached_radius=functools.lru_cache(maxsize=None)(mod._uncached_token_radius)
    mod.float=cached_float
    mod.token_radius=cached_radius
    def restore():
        cached_float.cache_clear(); cached_radius.cache_clear()
        mod.float=mod._uncached_float
        mod.token_radius=mod._uncached_token_radius
    return {'float_cache':cached_float,'radius_cache':cached_radius,'restore':restore}

def validate_static_assets(case, receipt, base_prep:Path=BASE_PREP):
    """Permit only the source-proven runtime replacement of str.out."""
    registered=case.get('static_assets',{})
    if receipt.get('static_assets')!=registered:
        raise ValueError(f"{case['run_id']}: prelaunch static-assets receipt differs from registration")
    work=Path(case['pre_dir']); base=base_prep/f"w1_{case['mode']}_{case['threshold']}"/'str.out'
    if not base.is_file() or sha(base)!=registered.get('str.out'):
        raise ValueError(f"{case['run_id']}: original registered str.out preimage is not preserved")
    for name,digest in registered.items():
        path=work/name
        if not path.is_file(): raise ValueError(f"{case['run_id']}: registered asset absent after run: {name}")
        if name!='str.out' and sha(path)!=digest:
            raise ValueError(f"{case['run_id']}: non-str.out static asset changed: {name}")
    current=work/'str.out'; post_sha=sha(current)
    return {'path':str(current),'post_run_sha256':post_sha,'post_run_size_bytes':current.stat().st_size,
      'registered_pre_sha256':registered['str.out'],'preserved_preimage_path':str(base),
      'preserved_preimage_sha256':sha(base),'changed_after_launch':post_sha!=registered['str.out']}

def verify_capture(case, receipt_path, receipt_sha, reg_sha, assembler):
    if not receipt_path.is_file() or sha(receipt_path)!=receipt_sha:
        raise ValueError(f"{case['run_id']}: receipt hash mismatch")
    run=json.loads(receipt_path.read_text())
    if run.get('run_id')!=case['run_id'] or run.get('status')!='PROCESS_EXITED' or run.get('returncode')!=0:
        raise ValueError(f"{case['run_id']}: unsuccessful or misidentified receipt")
    if run.get('registration_sha256')!=reg_sha or run.get('binary_sha256')!=BINARY_SHA:
        raise ValueError(f"{case['run_id']}: registration/binary binding mismatch")
    if run.get('input_nml_pre_sha256')!=case['input_nml_pre_sha256'] or run.get('seed_sha256')!=case['seed_sha256']:
        raise ValueError(f"{case['run_id']}: PRE input/seed binding mismatch")
    if run.get('pre_dir')!=str(Path(case['pre_dir']).resolve()):
        raise ValueError(f"{case['run_id']}: receipt PRE directory mismatch")
    out=case['outputs']; outputs=run.get('outputs',{}); post=run.get('post_files',{})
    required=('stdout','stderr','trace','final_state')
    for name in required:
        path=Path(out[name]).resolve()
        if outputs.get(name)!=str(path) or name not in post or not path.is_file() or sha(path)!=post[name].get('sha256'):
            raise ValueError(f"{case['run_id']}: missing/hash-mismatched output {name}")
    if Path(out['full_scf_log']).resolve()!=Path(out['stdout']).resolve():
        raise ValueError(f"{case['run_id']}: full log must alias stdout")
    env=run.get('environment',{})
    for name,value in case['environment'].items():
        if env.get(name)!=value: raise ValueError(f"{case['run_id']}: environment mismatch for {name}")
    pre=Path(case['pre_dir'])/'input.nml.pre'
    seed=Path(case['pre_dir'])/'Fe.nml.pre'
    if sha(pre)!=case['input_nml_pre_sha256'] or sha(seed)!=case['seed_sha256']:
        raise ValueError(f"{case['run_id']}: frozen PRE files changed")
    validate_static_assets(case,run)
    if run.get('environment',{}).get('RSLMTO_FIXED_U_TRACE_ITERATION')!='0' or env.get('RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY')!='1':
        raise ValueError(f"{case['run_id']}: not an all-iteration capture-only trace")
    nk,kdomain=assembler.reciprocal_mesh_from_pre(pre,case['input_nml_pre_sha256'])
    if nk!=1728 or kdomain!=list(range(1728)):
        raise ValueError(f"{case['run_id']}: response mesh differs from frozen full-BZ 12^3 model")
    return run,post

def check_locked_inputs(registration:Path, registration_sha:str, execution:Path, execution_sha:str):
    if sha(registration)!=registration_sha or sha(execution)!=execution_sha:
        raise ValueError('registration/execution hash mismatch')
    if (sha(SOURCE)!=SOURCE_SHA or sha(BINARY)!=BINARY_SHA or
        sha(STRUX_WRITER)!=STRUX_WRITER_SHA or sha(STRUX_REPLACE)!=STRUX_REPLACE_SHA):
        raise ValueError('frozen production source/binary/writer provenance changed')
    if sha(CALIBRATION_PREREG)!=CALIBRATION_PREREG_SHA: raise ValueError('calibration prereg changed')
    reg=json.loads(registration.read_text()); batch=json.loads(execution.read_text())
    if reg.get('schema')!='rslmto-a4-m21-w1-response-prereg-v1' or reg.get('status')!='PRE_RESPONSE_LOCKED_NO_NONZERO_RUNS':
        raise ValueError('unexpected response registration schema/status')
    cases=reg.get('response_cases',[])
    if len(cases)!=24 or {c.get('run_id') for c in cases}!=EXPECTED_RUNS: raise ValueError('response registration is not the exact 24-case W1 set')
    if batch.get('status')!='ALL_PROCESSES_EXITED' or batch.get('registration_sha256')!=registration_sha or batch.get('case_count')!=24:
        raise ValueError('execution batch is incomplete or not bound to registration')
    if set(batch.get('case_manifests',{}))!=EXPECTED_RUNS: raise ValueError('execution receipt set differs from locked response set')
    if set(batch.get('cases',{}))!=EXPECTED_RUNS or any(batch['cases'][rid].get('status')!='PROCESS_EXITED' or batch['cases'][rid].get('returncode')!=0 for rid in EXPECTED_RUNS):
        raise ValueError('one or more W1 response processes failed')
    return reg,batch

def make_inventory(registration:Path, registration_sha:str, execution:Path, execution_sha:str):
    reg,batch=check_locked_inputs(registration,registration_sha,execution,execution_sha)
    cases={c['run_id']:c for c in reg['response_cases']}; receipts={}
    for runid in sorted(EXPECTED_RUNS):
        case=cases[runid]; receipt_path=Path(batch['case_manifests'][runid])
        if not receipt_path.is_file(): raise ValueError(f'{runid}: missing case receipt')
        receipt_sha=sha(receipt_path); run=json.loads(receipt_path.read_text())
        if run.get('run_id')!=runid or run.get('registration_sha256')!=registration_sha or run.get('status')!='PROCESS_EXITED' or run.get('returncode')!=0:
            raise ValueError(f'{runid}: case receipt is not successful or bound')
        runtime_asset=validate_static_assets(case,run)
        outputs={}
        for name in ('stdout','stderr','trace','final_state'):
            path=Path(case['outputs'][name])
            post=run.get('post_files',{}).get(name,{})
            if not path.is_file() or not post.get('sha256') or sha(path)!=post['sha256']:
                raise ValueError(f'{runid}: {name} missing or hash differs from receipt')
            outputs[name]={'path':str(path),'sha256':post['sha256'],'size_bytes':post.get('size_bytes')}
        if Path(case['outputs']['full_scf_log']).resolve()!=Path(case['outputs']['stdout']).resolve():
            raise ValueError(f'{runid}: full_scf_log does not alias stdout')
        receipts[runid]={'case_manifest_path':str(receipt_path),'case_manifest_sha256':receipt_sha,
          'input_nml_pre_sha256':case['input_nml_pre_sha256'],'seed_sha256':case['seed_sha256'],'outputs':outputs,
          'runtime_mutated_assets':{'str.out':runtime_asset}}
    return {'schema':'rslmto-a4-m21-w1-response-hash-inventory-v1',
      'status':'HASH_INVENTORY_ONLY_NO_NUMERIC_CONTENT_READ',
      'registration_path':str(registration),'registration_sha256':registration_sha,
      'execution_manifest_path':str(execution),'execution_manifest_sha256':execution_sha,
      'case_count':24,'cases':receipts}

def validate_inventory(inventory, registration_sha:str, execution_sha:str):
    if inventory.get('schema')!='rslmto-a4-m21-w1-response-hash-inventory-v1' or inventory.get('status')!='HASH_INVENTORY_ONLY_NO_NUMERIC_CONTENT_READ':
        raise ValueError('unexpected response inventory schema/status')
    if inventory.get('registration_sha256')!=registration_sha or inventory.get('execution_manifest_sha256')!=execution_sha or inventory.get('case_count')!=24:
        raise ValueError('response hash inventory is not bound to exact registration/execution')
    cases=inventory.get('cases',{})
    if set(cases)!=EXPECTED_RUNS: raise ValueError('response hash inventory case set mismatch')
    for runid,row in cases.items():
        if not row.get('case_manifest_path') or not re.fullmatch(r'[0-9a-f]{64}',row.get('case_manifest_sha256','')):
            raise ValueError(f'{runid}: missing/invalid pinned receipt path/hash')
        if not re.fullmatch(r'[0-9a-f]{64}',row.get('input_nml_pre_sha256','')) or not re.fullmatch(r'[0-9a-f]{64}',row.get('seed_sha256','')):
            raise ValueError(f'{runid}: missing/invalid PRE/seed hashes in inventory')
        outs=row.get('outputs',{})
        if set(outs)!={'stdout','stderr','trace','final_state'}:
            raise ValueError(f'{runid}: inventory output set mismatch')
        for name,item in outs.items():
            if not item.get('path') or not re.fullmatch(r'[0-9a-f]{64}',item.get('sha256','')):
                raise ValueError(f'{runid}: missing/invalid {name} path/hash')
        runtime=row.get('runtime_mutated_assets',{}).get('str.out',{})
        if not runtime.get('path') or not runtime.get('preserved_preimage_path') or not isinstance(runtime.get('changed_after_launch'),bool):
            raise ValueError(f'{runid}: missing str.out pre/post provenance')
        for field in ('post_run_sha256','registered_pre_sha256','preserved_preimage_sha256'):
            if not re.fullmatch(r'[0-9a-f]{64}',runtime.get(field,'')):
                raise ValueError(f'{runid}: missing/invalid str.out {field}')
        if not isinstance(runtime.get('post_run_size_bytes'),int) or runtime['post_run_size_bytes']<0:
            raise ValueError(f'{runid}: invalid post-run str.out size')
    return cases

def assemble(registration:Path, registration_sha:str, execution:Path, execution_sha:str,
             inventory_path:Path, inventory_sha:str, zero_payload:Path, zero_sha:str):
    reg,batch=check_locked_inputs(registration,registration_sha,execution,execution_sha)
    if sha(inventory_path)!=inventory_sha: raise ValueError('response hash inventory digest mismatch')
    inventory=json.loads(inventory_path.read_text())
    inventory_cases=validate_inventory(inventory,registration_sha,execution_sha)
    if sha(zero_payload)!=zero_sha: raise ValueError('zero anchor payload hash mismatch')
    zero=json.loads(zero_payload.read_text())
    if zero.get('status')!='RAW_TOKENS_ASSEMBLED; NO_DERIVATIVE_RESPONSE_ANALYSIS' or len(zero.get('captures',[]))!=4:
        raise ValueError('zero anchor payload is not the reviewed W1 four-capture payload')
    if zero_sha!='ce78b15e3984f007b76d61fff5d032c48d64ae13bc60efcd9149b1607de7bb7e':
        raise ValueError('unexpected frozen W1 zero-anchor payload')
    prereg=json.loads(CALIBRATION_PREREG.read_text())
    # The certified assembler derives site count from the accepted base manifest/structure.
    assembler=load_assembler()
    nsite,siteprov=assembler.site_count_from_frozen_base(prereg)
    assembled=[]
    cases=reg['response_cases']
    case_by_id={c['run_id']:c for c in cases}
    with tempfile.TemporaryDirectory(prefix='a4_w1_response_normalized_manifests_') as td:
      td=Path(td)
      for run_id in sorted(EXPECTED_RUNS):
        case=case_by_id[run_id]
        receipt_path=Path(batch['case_manifests'][run_id])
        invrow=inventory_cases[run_id]
        if str(receipt_path)!=invrow['case_manifest_path']:
            raise ValueError(f'{run_id}: receipt path differs from hash inventory')
        receipt_sha=invrow['case_manifest_sha256']
        if invrow['input_nml_pre_sha256']!=case['input_nml_pre_sha256'] or invrow['seed_sha256']!=case['seed_sha256']:
            raise ValueError(f'{run_id}: inventory PRE/seed hashes differ from response registration')
        runtime=invrow['runtime_mutated_assets']['str.out']
        expected_pre_strout=case['static_assets']['str.out']
        if runtime['registered_pre_sha256']!=expected_pre_strout or runtime['preserved_preimage_sha256']!=expected_pre_strout:
            raise ValueError(f'{run_id}: registered/original str.out preimage mismatch')
        if Path(runtime['preserved_preimage_path'])!=BASE_PREP/f"w1_{case['mode']}_{case['threshold']}"/'str.out':
            raise ValueError(f'{run_id}: unexpected original str.out preimage path')
        if Path(runtime['path'])!=Path(case['pre_dir'])/'str.out':
            raise ValueError(f'{run_id}: unexpected runtime-mutated str.out path')
        if Path(runtime['path'])!=Path(case['pre_dir'])/'str.out':
            raise ValueError(f'{run_id}: unexpected runtime-mutated str.out path')
        if sha(Path(runtime['preserved_preimage_path']))!=expected_pre_strout or sha(Path(runtime['path']))!=runtime['post_run_sha256']:
            raise ValueError(f'{run_id}: str.out pre/post hashes changed after inventory')
        for output_name in ('stdout','stderr','trace','final_state'):
            frozen=invrow.get('outputs',{}).get(output_name,{})
            expected_path=str(Path(case['outputs'][output_name]).resolve())
            if frozen.get('path')!=expected_path or not frozen.get('sha256') or sha(Path(expected_path))!=frozen['sha256']:
                raise ValueError(f'{run_id}: {output_name} differs from frozen response inventory')
        run,post=verify_capture(case,receipt_path,receipt_sha,registration_sha,assembler)
        env=run['environment']
        nk,kdomain=assembler.reciprocal_mesh_from_pre(Path(case['pre_dir'])/'input.nml.pre',case['input_nml_pre_sha256'])
        mode=case['mode']; threshold=case['threshold']; lam=float(case['lambda_ry'])
        runid=case['run_id']
        # Adapt the executor receipt to the certified assembler's manifest schema
        # in a private temporary file; original receipt remains hash-bound below.
        normalized={
          'run_id':runid,'status':'PROCESS_EXITED','returncode':0,
          'bindings':{'binary_sha256':run['binary_sha256'],
                      'input_nml_pre_sha256':run['input_nml_pre_sha256']},
          'environment':env,
          'post_files':{'records.csv':post['trace'],'stdout.log':post['stdout'],
                        'Fe_out.nml':post['final_state']}}
        normalized_path=td/f'{runid}.json'
        normalized_path.write_text(json.dumps(normalized,sort_keys=True)+'\n')
        meta={
          'run_id':runid,'mode':mode,'threshold':threshold,
          'run_manifest_path':str(normalized_path),'run_manifest_sha256':sha(normalized_path),
          'records_csv_path':case['outputs']['trace'],'records_csv_sha256':post['trace']['sha256'],
          'stdout_path':case['outputs']['stdout'],'stdout_sha256':post['stdout']['sha256'],
          'input_nml_pre_path':str(Path(case['pre_dir'])/'input.nml.pre'),
          'input_nml_pre_sha256':case['input_nml_pre_sha256'],
          'expected_n_kpoints':nk,'expected_k_indices':kdomain,'expected_site_count':nsite,
          'site_provenance':siteprov,'source_sha256':SOURCE_SHA,'binary_sha256':BINARY_SHA}
        caches=install_exact_caches(assembler)
        try:
            row=assembler.assemble_case(meta)
        finally:
            caches['restore']()
        envlam=env.get('RSLMTO_FIXED_U_LAMBDA')
        if envlam is None or float(envlam)!=lam: raise ValueError(f'{runid}: lambda token/runtime mismatch')
        if row['lambda_runtime_text']!=envlam: raise ValueError(f'{runid}: trace lambda mismatch')
        row['lambda_ry']=lam
        row['operator']='w1'
        row['run_manifest_path']=str(receipt_path)
        row['run_manifest_sha256']=receipt_sha
        row['run_receipt_sha256']=receipt_sha
        row['execution_registration_sha256']=registration_sha
        row['response_case_input_sha256']=case['input_nml_pre_sha256']
        assembled.append(row)
    zero_rows=[r for r in zero['captures'] if r['mode'] in ('local','spectral') and r['lambda_ry']==0.0]
    if len(zero_rows)!=4: raise ValueError('expected four W1 zero anchor rows')
    all_rows=zero_rows+assembled
    seen={(r['mode'],float(r['lambda_ry']),r['threshold']) for r in all_rows}
    expected={(m,l,t) for m in ('local','spectral') for l in (0.0,-.0025,.0025,-.005,.005,-.01,.01) for t in ('standard','tight')}
    if seen!=expected or len(all_rows)!=28: raise ValueError('combined 28-case domain has gaps or duplicates')
    return {'schema':'rslmto-a4-m21-w1-responses-assembled-v1',
      'status':'RAW_RESPONSE_TOKENS_ASSEMBLED_NO_GATE_EVALUATION',
      'registration_path':str(registration),'registration_sha256':registration_sha,
      'execution_manifest_path':str(execution),'execution_manifest_sha256':execution_sha,
      'response_hash_inventory_path':str(inventory_path),'response_hash_inventory_sha256':inventory_sha,
      'zero_anchor_path':str(zero_payload),'zero_anchor_sha256':zero_sha,
      'assembler_path':str(ASSEMBLER),'assembler_sha256':ASSEMBLER_SHA,
      'str_out_writer_sha256':STRUX_WRITER_SHA,'str_out_replace_open_sha256':STRUX_REPLACE_SHA,
      'cache_policy':'Memoize exact built-in float(token) and certified token_radius(token) results; formulas unchanged. Synthetic test compares complete assembled object byte-for-byte with uncached implementation.',
      'response_capture_count':24,'zero_anchor_capture_count':4,'capture_count':28,
      'captures':all_rows}

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='command',required=True)
    pi=sub.add_parser('inventory'); pi.add_argument('registration',type=Path); pi.add_argument('--registration-sha256',required=True); pi.add_argument('execution_manifest',type=Path); pi.add_argument('--execution-sha256',required=True); pi.add_argument('--output',type=Path,required=True)
    pa=sub.add_parser('assemble'); pa.add_argument('registration',type=Path); pa.add_argument('--registration-sha256',required=True); pa.add_argument('execution_manifest',type=Path); pa.add_argument('--execution-sha256',required=True); pa.add_argument('inventory',type=Path); pa.add_argument('--inventory-sha256',required=True); pa.add_argument('zero_anchor',type=Path); pa.add_argument('--zero-anchor-sha256',required=True); pa.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise SystemExit(f'refusing overwrite: {a.output}')
    if a.command=='inventory': result=make_inventory(a.registration,a.registration_sha256,a.execution_manifest,a.execution_sha256)
    else: result=assemble(a.registration,a.registration_sha256,a.execution_manifest,a.execution_sha256,a.inventory,a.inventory_sha256,a.zero_anchor,a.zero_anchor_sha256)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'case_count':result.get('case_count'),'capture_count':result.get('capture_count'),'path':str(a.output),'sha256':sha(a.output)},sort_keys=True))
if __name__=='__main__': main()
