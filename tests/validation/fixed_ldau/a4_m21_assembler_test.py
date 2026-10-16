#!/usr/bin/env python3
"""Deterministic multi-iteration raw-token assembler fixture; no production data."""
from __future__ import annotations
import csv, hashlib, importlib.util, json, tempfile
from pathlib import Path

ROOT=Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
SCRIPT=ROOT/'tests/validation/fixed_ldau/a4_m21_assemble_lambda0.py'
spec=importlib.util.spec_from_file_location('a4_assembler',SCRIPT); assert spec and spec.loader
asm=importlib.util.module_from_spec(spec); spec.loader.exec_module(asm)
FIELDS=asm.FIELDS

def row(kind: str, iteration: int, **kw) -> dict[str,str]:
    base={k:'-1' for k in FIELDS}
    base.update({'kind':kind,'run_id':'w1_spectral_standard','lambda':'0.0000000000000000E+00',
      'iteration':str(iteration),'stage':'EIGENPAIRS','weight':'0.0000000000000000E+00',
      'energy':'0.0000000000000000E+00','scalar':'0.0000000000000000E+00',
      're':'0.0000000000000000E+00','im':'0.0000000000000000E+00'})
    base.update({k:str(v) for k,v in kw.items()}); return base

def fixture_csv(path: Path, *, missing_last_k=False, duplicate_d=False) -> None:
    rows=[]; nsite=2
    for it in (1,2):
        for site in range(nsite):
            for spin in (0,1):
                base=site*18+spin*9
                for l in range(3):
                    for m in range(-l,l+1):
                        idx=base+l*l+(m+l)
                        col=m
                        if duplicate_d and site==1 and spin==1 and l==2 and m==2: col-=1
                        rows.append(row('basis_map',it,row=idx,spin=spin,l=l,col=col,site=site,type=1))
        for k in (0,1):
            if missing_last_k and k==1: continue
            rows.append(row('k_weight',it,k=k,weight='5.0000000000000000E-01'))
            for b in range(36):
                eval_text='0.0000000000000000E+00' if it==2 else '1.2500000000000000E-01'
                occ_text='5.0000000000000000E-01' if it==2 else '2.5000000000000000E-01'
                rows.append(row('eigenvalue',it,k=k,band=b,energy=eval_text))
                rows.append(row('occupation',it,k=k,band=b,scalar=occ_text))
                for i in range(18*nsite):
                    rows.append(row('eigenvector',it,k=k,band=b,row=i,
                        re='1.0000000000000000E+00' if i==b else '0.0000000000000000E+00'))
        rows.extend([row('fermi_level',it,scalar='0.0000000000000000E+00'),
            row('temperature_k',it,scalar='3.0000000000000000E+02'),
            row('temperature_ry',it,scalar='1.8990000000000000E-03')])
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); w.writerows(rows)

def main() -> None:
  with tempfile.TemporaryDirectory(prefix='a4_assembler_synth_') as td:
    d=Path(td); csvp=d/'records.csv'; logp=d/'stdout.log'; mp=d/'run_manifest.json'; pre=d/'input.nml.pre'
    pre.write_text('&reciprocal\nnk1=1\nnk2=1\nnk3=2\nuse_symmetry_reduction=.false.\n/\n&symmetry\nuse_symmetry_reduction=.true.\n/\n')
    fixture_csv(csvp)
    logp.write_text('''HUBBARD_MOMENT_ENERGY_TOTAL physical_total_energy= -1.0000000000000000E+00 TS= 0.1000000000000000E+00 entropy= 0.0 candidate_a= -1.0 candidate_b=UNAVAILABLE test_F_local=UNAVAILABLE test_F_spectral= -1.1 candidate_status=PROVISIONAL_ONLY;\nHUBBARD_MOMENT_ENERGY_TOTAL physical_total_energy= -2.0000000000000000E+00 TS= 0.2000000000000000E+00 entropy= 0.0 candidate_a= -2.0 candidate_b=UNAVAILABLE test_F_local=UNAVAILABLE test_F_spectral= -2.2 candidate_status=PROVISIONAL_ONLY;\nConverged!\n''')
    mp.write_text(json.dumps({'run_id':'w1_spectral_standard','bindings':{'binary_sha256':'binary','input_nml_pre_sha256':hashlib.sha256(pre.read_bytes()).hexdigest()},'environment':{
      'RSLMTO_FIXED_U_TRACE_ITERATION':'0','RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY':'1','RSLMTO_FIXED_U_LAMBDA':'0.0'}}))
    meta={'records_csv_path':str(csvp),'records_csv_sha256':hashlib.sha256(csvp.read_bytes()).hexdigest(),
      'stdout_path':str(logp),'stdout_sha256':hashlib.sha256(logp.read_bytes()).hexdigest(),
      'run_manifest_path':str(mp),'run_manifest_sha256':hashlib.sha256(mp.read_bytes()).hexdigest(),
      'run_id':'w1_spectral_standard','binary_sha256':'binary','input_nml_pre_sha256':hashlib.sha256(pre.read_bytes()).hexdigest(),
      'input_nml_pre_path':str(pre),'input_nml_pre_sha256':hashlib.sha256(pre.read_bytes()).hexdigest(),
      'expected_n_kpoints':2,'expected_k_indices':[0,1],'expected_site_count':2,
      'mode':'spectral','threshold':'standard','source_sha256':'source','binary_sha256':'binary'}
    result=asm.assemble_case(meta)
    assert result['expected_n_kpoints']==2 and result['expected_k_indices']==[0,1]
    assert result['expected_band_indices']==list(range(36)) and result['final_iteration']==2
    assert result['ETOT_spectral_text']=='-2.0000000000000000E+00' and result['TS_text']=='0.2000000000000000E+00'
    assert result['P_d'][0][0][0]==2.0 and result['P_d'][4][4][0]==2.0
    assert result['P_d_by_spin']['0'][0][0][0]==1.0 and result['P_d_by_spin']['1'][4][4][0]==1.0
    assert result['entropy_inputs']['states'][0]['eigenvalue_text']=='0.0000000000000000E+00'
    # Phase-only component uncertainty must be bounded even when |a| changes only quadratically.
    q=asm.half_quantum('1.0000000000000000E+00')
    w,f,a,b=0.5,0.5,complex(1.0,0.0),complex(1.0,0.0)
    qterm=asm.complex_product_token_error_bound(w,f,a,b,q,q,2**0.5*q,2**0.5*q)
    hidden_a=complex(1.0,q)
    phase_delta=abs(w*f*hidden_a*b.conjugate()-w*f*a*b.conjugate())
    assert 0.0 < phase_delta <= qterm
    token_text='1.2345678901234567E+00'
    assert asm.token_radius(token_text) >= asm.half_quantum(token_text)
    # Missing trailing global k must fail against PRE-derived exact range.
    miss=d/'missing.csv'; fixture_csv(miss,missing_last_k=True)
    bad=dict(meta,records_csv_path=str(miss),records_csv_sha256=hashlib.sha256(miss.read_bytes()).hexdigest())
    try: asm.assemble_case(bad)
    except ValueError as exc: assert 'incomplete zero-based k/band' in str(exc)
    else: raise AssertionError('missing final k block accepted')
    # Duplicate spherical d label must fail instead of overwriting a mapping.
    dup=d/'duplicate.csv'; fixture_csv(dup,duplicate_d=True)
    bad=dict(meta,records_csv_path=str(dup),records_csv_sha256=hashlib.sha256(dup.read_bytes()).hexdigest())
    try: asm.assemble_case(bad)
    except ValueError as exc: assert 'duplicate site/spin/d-m' in str(exc)
    else: raise AssertionError('duplicate d label accepted')
    shifted=d/'shifted.pre'; shifted.write_text('&reciprocal\nnk1=1\nnk2=1\nnk3=2\nuse_symmetry_reduction=.true.\n/\n')
    try: asm.reciprocal_mesh_from_pre(shifted,hashlib.sha256(shifted.read_bytes()).hexdigest())
    except ValueError as exc: assert 'unreduced' in str(exc)
    else: raise AssertionError('symmetry-reduced mesh accepted as full grid')
  print('PASS synthetic multi-iteration assembly, PRE-derived exact k domain, 2-site/2-spin band map, d-label uniqueness, last-log-line selection, and exact-token retention')

if __name__=='__main__': main()
