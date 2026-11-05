import importlib.util, json, math
from pathlib import Path
P=Path('/private/tmp/c52-rd2-private/nr501_adapter/rd2_replay_adapter_nr501.py')
spec=importlib.util.spec_from_file_location('rd2adapter',P); a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
runner,rows,nr,ncore,nval,vin,rec,binding=a.load_bound_rows(Path('/private/tmp/c52-rd2-ni-u3j07-capture/run'),Path('/private/tmp/c52-rd2-private/vin_binding_nr501.json'))
a.runner_np=runner.np
oracle=a.local_moment_oracle()
records=[]
for spin in (1,2):
 p=oracle.matrix(rows,kind='ldm',stage='P_OUT',iteration=1,site=0,ell=2,spin=spin,dim=5)
 m1=oracle.matrix(rows,kind='ldm_m1',stage='M1_OUT',iteration=1,site=0,ell=2,spin=spin,dim=5)
 gm=a.raw_green_moments(rows,site=0,ell=2,spin=spin,expected_energy_points=binding['capture_completeness']['energy_points'])
 q0=float(runner.np.trace(p).real); center_m1=float(runner.np.trace(m1).real)/q0
 qrow=runner.one(rows,'m0_carrier_input','FROZEN_REPLAY',spin=spin,ell=2)
 result=oracle.local_correction_no_controls(p,m1,vin[spin-1],q0,gm['center_raw_ry'],gm['q2_raw_ry2'])
 records.append({'spin':spin,**gm,'p_trace':q0,'center_from_m1':center_m1,'m0_writer_type':int(qrow['type']),
   'm0_site':int(qrow['site']),'m0_post_local_center':float(qrow['energy']),'m0_post_local_q2':float(qrow['re']),
   'vin_trace_p':float(runner.np.trace(vin[spin-1]@p).real),
   'center_local_reconstructed':result['c_local_ry'],'q2_local_reconstructed':result['q2_local_ry2'],
   'p1_diff_green_vs_pout':gm['p1']-q0,'raw_center_diff_m1_vs_green':center_m1-gm['center_raw_ry'],
   'post_local_center_diff':result['c_local_ry']-float(qrow['energy']),
   'post_local_q2_diff':result['q2_local_ry2']-float(qrow['re'])})
out={'schema':'rd2-green-raw-phase-proof-v1','capture_csv_sha256':a.sha(Path('/private/tmp/c52-rd2-ni-u3j07-capture/run/records.csv')),
 'binding_sha256':a.sha(Path('/private/tmp/c52-rd2-private/vin_binding_nr501.json')),
 'Vin_map':'LOCAL_HALL runtime nlim=1/type2; records are recursion site0; m0 type=19 is lattice-ib writer metadata',
 'phase_order':'bands correction precedes P_OUT/m0 export; compare reconstructed one-Vin correction with captured post-correction center/q2',
 'records':records}
path=Path('/private/tmp/c52-rd2-private/green_phase_proof_001.json'); path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'output':str(path),'sha256':a.sha(path),'records':records},indent=2))
