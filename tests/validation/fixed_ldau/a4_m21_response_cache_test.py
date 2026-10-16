#!/usr/bin/env python3
"""Prove cached W1 response assembly is byte-for-byte scalar-equivalent."""
from __future__ import annotations
import csv, hashlib, importlib.util, json, tempfile
from pathlib import Path

ROOT=Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
TEST_FIXTURE=ROOT/'tests/validation/fixed_ldau/a4_m21_assembler_test.py'
RESPONSE_ASSEMBLER=Path('/private/tmp/a4_h2_stage/a4_m21_assemble_w1_responses.py')

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path); assert spec and spec.loader
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module

def main():
    fixture=load('certified_assembler_fixture',TEST_FIXTURE)
    response=load('response_collector_under_test',RESPONSE_ASSEMBLER)
    fake_cases={rid:{'case_manifest_path':f'/receipts/{rid}.json','case_manifest_sha256':'a'*64,
      'input_nml_pre_sha256':'b'*64,'seed_sha256':'c'*64,
      'outputs':{name:{'path':f'/outputs/{rid}/{name}','sha256':'d'*64,'size_bytes':1}
                 for name in ('stdout','stderr','trace','final_state')},
      'runtime_mutated_assets':{'str.out':{'path':f'/outputs/{rid}/str.out','post_run_sha256':'1'*64,
       'post_run_size_bytes':1,'registered_pre_sha256':'2'*64,'preserved_preimage_path':f'/base/{rid}/str.out',
       'preserved_preimage_sha256':'2'*64,'changed_after_launch':True}}} for rid in response.EXPECTED_RUNS}
    fake_inventory={'schema':'rslmto-a4-m21-w1-response-hash-inventory-v1',
      'status':'HASH_INVENTORY_ONLY_NO_NUMERIC_CONTENT_READ','registration_sha256':'e'*64,
      'execution_manifest_sha256':'f'*64,'case_count':24,'cases':fake_cases}
    assert set(response.validate_inventory(fake_inventory,'e'*64,'f'*64))==response.EXPECTED_RUNS
    missing_hash=json.loads(json.dumps(fake_inventory))
    missing_hash['cases'][sorted(response.EXPECTED_RUNS)[0]].pop('case_manifest_sha256')
    try: response.validate_inventory(missing_hash,'e'*64,'f'*64)
    except ValueError as exc: assert 'receipt path/hash' in str(exc)
    else: raise AssertionError('hash inventory accepted missing receipt digest')
    with tempfile.TemporaryDirectory(prefix='a4_strout_mutation_guard_') as td:
        d=Path(td); work=d/'work'; orig=d/'base'/'w1_spectral_standard'; work.mkdir(); orig.mkdir(parents=True)
        names=('clust','map','sbar','view.sbar','str.out','ves.out','mad.mat')
        original={name:(b'original-strout' if name=='str.out' else f'input-{name}'.encode()) for name in names}
        for name,content in original.items(): (orig/name if name=='str.out' else work/name).write_bytes(content)
        (work/'str.out').write_bytes(b'program-replaced-strout')
        static={name:hashlib.sha256(original[name]).hexdigest() for name in names}
        case={'run_id':'synthetic','mode':'spectral','threshold':'standard','pre_dir':str(work),'static_assets':static}
        receipt={'static_assets':dict(static)}
        mutation=response.validate_static_assets(case,receipt,base_prep=d/'base')
        assert mutation['changed_after_launch'] is True and mutation['preserved_preimage_sha256']==static['str.out']
        (work/'ves.out').write_bytes(b'unexpected-edit')
        try: response.validate_static_assets(case,receipt,base_prep=d/'base')
        except ValueError as exc: assert 'non-str.out static asset changed' in str(exc)
        else: raise AssertionError('unexpected non-str.out runtime mutation accepted')
        (work/'ves.out').write_bytes(original['ves.out']); del receipt['static_assets']['str.out']
        try: response.validate_static_assets(case,receipt,base_prep=d/'base')
        except ValueError as exc: assert 'prelaunch static-assets receipt' in str(exc)
        else: raise AssertionError('missing prelaunch str.out receipt entry accepted')
    with tempfile.TemporaryDirectory(prefix='a4_w1_cache_equivalence_') as td:
        d=Path(td); csvp=d/'records.csv'; logp=d/'stdout.log'; mp=d/'run_manifest.json'; pre=d/'input.nml.pre'
        fixture.fixture_csv(csvp)
        pre.write_text('&reciprocal\nnk1=1\nnk2=1\nnk3=2\nuse_symmetry_reduction=.false.\n/\n&symmetry\nuse_symmetry_reduction=.true.\n/\n')
        logp.write_text('''HUBBARD_MOMENT_ENERGY_TOTAL physical_total_energy= -1.0000000000000000E+00 TS= 0.1000000000000000E+00 entropy= 0.0 candidate_a= -1.0 candidate_b=UNAVAILABLE test_F_local=UNAVAILABLE test_F_spectral= -1.1 candidate_status=PROVISIONAL_ONLY;\nHUBBARD_MOMENT_ENERGY_TOTAL physical_total_energy= -2.0000000000000000E+00 TS= 0.2000000000000000E+00 entropy= 0.0 candidate_a= -2.0 candidate_b=UNAVAILABLE test_F_local=UNAVAILABLE test_F_spectral= -2.2 candidate_status=PROVISIONAL_ONLY;\nConverged!\n''')
        pre_sha=hashlib.sha256(pre.read_bytes()).hexdigest(); csv_sha=hashlib.sha256(csvp.read_bytes()).hexdigest(); log_sha=hashlib.sha256(logp.read_bytes()).hexdigest()
        mp.write_text(json.dumps({'run_id':'w1_spectral_standard','bindings':{'binary_sha256':'binary','input_nml_pre_sha256':pre_sha},'environment':{'RSLMTO_FIXED_U_TRACE_ITERATION':'0','RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY':'1','RSLMTO_FIXED_U_LAMBDA':'0.0'}}))
        msha=hashlib.sha256(mp.read_bytes()).hexdigest()
        meta={'records_csv_path':str(csvp),'records_csv_sha256':csv_sha,'stdout_path':str(logp),'stdout_sha256':log_sha,
          'run_manifest_path':str(mp),'run_manifest_sha256':msha,'run_id':'w1_spectral_standard','binary_sha256':'binary',
          'input_nml_pre_sha256':pre_sha,'input_nml_pre_path':str(pre),'expected_n_kpoints':2,'expected_k_indices':[0,1],
          'expected_site_count':2,'mode':'spectral','threshold':'standard','source_sha256':'source'}
        uncached=fixture.asm.assemble_case(meta)
        cached_mod=response.load_assembler()
        for _case in range(2):
            caches=response.install_exact_caches(cached_mod)
            assert caches['float_cache'].cache_info().currsize==0 and caches['radius_cache'].cache_info().currsize==0
            cached=cached_mod.assemble_case(meta)
            assert cached==uncached, 'cached token conversion changed any assembled field'
            assert caches['float_cache'].cache_info().hits>0 and caches['radius_cache'].cache_info().hits>0
            caches['restore']()
            assert cached_mod.float is cached_mod._uncached_float
            assert cached_mod.token_radius is cached_mod._uncached_token_radius
        # The phase-only perturbation included by the certified fixture remains enclosed.
        token='1.2345678901234567E+00'
        assert cached_mod.token_radius(token)==fixture.asm.token_radius(token)
    print('PASS response inventory receipts; source-proven str.out overwrite only; cached assembly equals certified scalar path and clears per case')

if __name__=='__main__': main()
