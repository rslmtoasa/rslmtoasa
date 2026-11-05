import importlib.util, json
from pathlib import Path
P=Path('/private/tmp/c52-rd2-private/nr501_adapter/rd2_replay_adapter_nr501.py')
spec=importlib.util.spec_from_file_location('rd2adapter',P); a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
cap=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run'); bind=Path('/private/tmp/c52-rd2-private/vin_binding_nr501.json')
runner,rows,nr,ncore,nval,vin,rec,binding=a.load_bound_rows(cap,bind); a.runner_np=runner.np
oracle=a.local_moment_oracle(); baseline=runner.prepare(rows,nr,ncore,nval,'baseline')
vmad=float(runner.one(rows,'radial_replay_vmad','RADIAL_REPLAY')['scalar']); out=[]
for spin in (1,2):
 idx=spin-1
 p=oracle.matrix(rows,kind='ldm',stage='P_OUT',iteration=1,site=0,ell=2,spin=spin,dim=5)
 m1=oracle.matrix(rows,kind='ldm_m1',stage='M1_OUT',iteration=1,site=0,ell=2,spin=spin,dim=5)
 gm=a.raw_green_moments(rows,site=0,ell=2,spin=spin,expected_energy_points=binding['capture_completeness']['energy_points'])
 qrow=runner.one(rows,'m0_carrier_input','FROZEN_REPLAY',spin=spin,ell=2)
 q0=float(runner.np.trace(p).real); local_gc=float(baseline['gc'][2,idx]); c_local=float(qrow['energy'])
 correction=oracle.local_correction_no_controls(p,m1,vin[idx],q0,gm['center_raw_ry'],gm['q2_raw_ry2'])
 x=float(runner.np.trace(vin[idx]@p).real)/q0
 vmad_pre=c_local-local_gc
 raw_gc=local_gc+(gm['center_raw_ry']-c_local)
 out.append({'spin':spin,'q0_Ptrace':q0,'green_p1_p2_p3':[gm['p1'],gm['p2'],gm['p3']],
 'raw_center_green':gm['center_raw_ry'],'raw_center_M1_trace':float(runner.np.trace(m1).real)/q0,
 'raw_q2_green':gm['q2_raw_ry2'],'vin_shift_x':x,'local_center_m0_including_vmad':c_local,
 'captured_local_calculate_pl_gc':local_gc,'same_phase_vmad_pre':vmad_pre,
 'raw_replay_gc_excluding_vmad':raw_gc,'local_replay_gc_after_shift':raw_gc-x,
 'later_atomsc_vmad_diagnostic_only':vmad,'later_vmad_minus_pre':vmad-vmad_pre,
 'corrected_center_formula':correction['c_local_ry'],'corrected_q2_formula':correction['q2_local_ry2'],
 'captured_local_q2':float(qrow['re'])})
proof={'schema':'rd2-green-raw-phase-proof-v2','capture_csv_sha256':a.sha(cap/'records.csv'),'binding_sha256':a.sha(bind),
 'vin_site_type':'LOCAL_HALL site/nlim1 type2 (Ni-imp); records site0; m0 type19 is lattice-ib writer metadata',
 'source_phase':'P_OUT/m0 are post-local-correction; raw moments from Green/M1. Gravity-center slot excludes vmad; construct from same-phase CALCULATE_PL local gc plus (raw c - captured m0 local c), not later ATOMSC vmad.',
 'records':out}
p=Path('/private/tmp/c52-rd2-private/green_phase_proof_002.json'); p.write_text(json.dumps(proof,indent=2,sort_keys=True)+'\n')
print(json.dumps({'path':str(p),'sha256':a.sha(p),'records':out},indent=2))
