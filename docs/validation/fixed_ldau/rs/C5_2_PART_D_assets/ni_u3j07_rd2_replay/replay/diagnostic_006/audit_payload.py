#!/usr/bin/env python3
import csv, json, hashlib, math, sys
from pathlib import Path
from collections import defaultdict
import numpy as np

CAP=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run')
BASE=Path('/private/tmp/c52-rd2-private')
BIND=BASE/'vin_binding_nr501.json'
PREG=BASE/'prereg_005/preregistration.json'
INPUT=BASE/'prereg_005/driver.input.txt'
FIRST=BASE/'cases_005/driver.output.csv'
AUD=CAP/'c2_hubbard_audit.dat'
CSV=CAP/'records.csv'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def cplx(r): return complex(float(r['re']),float(r['im']))
def read_json(p):return json.loads(p.read_text())

binding=read_json(BIND); preg=read_json(PREG)
PARSER=Path('/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau/c2_hubbard_operator_audit.py')
MATH=Path('/Users/andersb/Jobb/rslmto_fable_v4b/source/math.f90')
if sha(CSV)!=binding['iteration_and_capture']['csv_sha256']:raise ValueError('capture CSV hash differs from immutable binding')
if sha(AUD)!=binding['source']['sha256']:raise ValueError('LOCAL_HALL audit hash differs from immutable binding')
if sha(PARSER)!=binding['source']['parser_sha256']:raise ValueError('LOCAL_HALL parser source hash differs from binding')
if sha(MATH)!=binding['math_binding']['production_sha256']:raise ValueError('production moment source hash differs from binding')
if sha(INPUT)!=preg['case_driver_input_sha256']:raise ValueError('actual case payload differs from preregistration hash')
# Stream the immutable capture once into only the selected exact matrices and Green blocks.
selected=[]; gblocks=defaultdict(list); bmap={}
with CSV.open(newline='') as f:
 for r in csv.DictReader(f):
  k=r['kind']; st=r['stage']; site=int(r['site']); ell=int(r['l']); sp=int(r['spin'])
  if k=='basis_map' and st=='RS_GREEN' and site==0 and ell==2:
   bmap.setdefault(sp,{})[int(r['col'])]=int(r['row'])
  if (k in ('ldm','ldm_m1') and st in ('P_OUT','M1_OUT') and site==0 and ell==2 and sp in (1,2)):
   selected.append(r)
  if k=='v_hub' and st=='Pin_H' and int(r['type'])==2:
   selected.append(r)
  if k=='pl_map_input' and st=='CALCULATE_PL' and site==0 and ell in (0,1,2) and sp in (1,2):
   selected.append(r)
  if k=='radial_replay_pl' and st=='RADIAL_REPLAY' and ell in (0,1,2) and sp in (0,1):
   selected.append(r)
  if k=='radial_moment' and st=='P_OUT' and site==0 and ell==2 and sp in (1,2): selected.append(r)
  if k=='radial_replay_ql' and st=='RADIAL_REPLAY' and site==0 and ell in (0,1,2) and sp in (0,1): selected.append(r)
  if k=='m0_carrier_input' and st=='FROZEN_REPLAY' and site==0 and ell==2 and sp in (1,2): selected.append(r)
  if k=='green_function' and st=='RS_GREEN' and site==0: gblocks[int(r['k'])].append(r)
# exact P/M1 from CSV
P={s:np.zeros((5,5),complex) for s in (1,2)}; M1={s:np.zeros((5,5),complex) for s in (1,2)}
VH=np.zeros((18,18),complex); vh_seen=set()
seen={('P',s):set() for s in (1,2)}|{('M1',s):set() for s in (1,2)}
rad={s:{} for s in (1,2)}; liveq={s:{} for s in (1,2)}; m0={s:{} for s in (1,2)}
plmap={s:{ell:{} for ell in (0,1,2)} for s in (1,2)}; pllive={s:{} for s in (1,2)}
for r in selected:
 k=r['kind']; st=r['stage']; s=int(r['spin'])
 if k=='ldm':
  ij=(int(r['row']),int(r['col'])); P[s][ij]=cplx(r); seen[('P',s)].add(ij)
 elif k=='ldm_m1':
  ij=(int(r['row']),int(r['col'])); M1[s][ij]=cplx(r); seen[('M1',s)].add(ij)
 elif k=='v_hub':
  ij=(int(r['row']),int(r['col'])); VH[ij]=cplx(r); vh_seen.add(ij)
 elif k=='pl_map_input':
  plmap[s][int(r['l'])][int(r['order'])]=float(r['scalar'])
 elif k=='radial_replay_pl':
  pllive[s+1][int(r['l'])]=float(r['scalar'])
 elif k=='radial_moment':
  rad[s][int(r['order'])]={'energy':float(r['energy']),'scalar':float(r['scalar'])}
 elif k=='radial_replay_ql': liveq[s+1][(int(r['l']),int(r['order']))]=float(r['scalar'])
 elif k=='m0_carrier_input': m0[s][0]={'scalar':float(r['scalar']),'re':float(r['re']),'energy':float(r['energy']),'type':int(r['type'])}
for s in (1,2):
 if len(seen['P',s])!=25 or len(seen['M1',s])!=25:raise ValueError(f'incomplete captured P/M1 {s}')
 if not np.allclose(P[s],P[s].conj().T,atol=1e-10,rtol=0):raise ValueError('P nonhermitian')
 if not np.allclose(M1[s],M1[s].conj().T,atol=1e-10,rtol=0):raise ValueError('M1 nonhermitian')
if len(vh_seen)!=324:raise ValueError('incomplete type2 v_hub Pin_H 18x18 matrix')
# d-shell Green moment matrices P, M1, M2 independently from raw G(E).
expected_m=list(range(-2,3)); maps={}
for s in (0,1):
 z=bmap.get(s,{})
 if sorted(z)!=expected_m:raise ValueError(f'bad d basis map spin={s}: {z}')
 maps[s]=[z[m] for m in expected_m]
Gmom={s:[np.zeros((5,5),complex) for _ in range(3)] for s in (1,2)}
energy_weights={}
for ix,block in gblocks.items():
 G=np.zeros((18,18),complex); seenij=set(); es=set(); ws=set()
 for r in block:
  ij=(int(r['row']),int(r['col']))
  if ij in seenij:raise ValueError('duplicate G entry')
  seenij.add(ij);G[ij]=cplx(r);es.add(float(r['energy']));ws.add(float(r['weight']))
 if len(seenij)!=324 or len(es)!=1 or len(ws)!=1:raise ValueError(f'incomplete Green {ix}')
 E=next(iter(es)); w=next(iter(ws));energy_weights[ix]=(E,w)
 S=1j*(G-G.conj().T)/(2*math.pi)
 for s in (1,2):
  sub=S[np.ix_(maps[s-1],maps[s-1])]
  for n in range(3):Gmom[s][n]+=w*(E**n)*sub
if len(gblocks)!=1567:raise ValueError('Green block count mismatch')
# exact operator capture, via audited parser
sys.path.insert(0,'/Users/andersb/Jobb/rslmto_fable_v4b/tests/validation/fixed_ldau')
import c2_hubbard_operator_audit as c2
rec=c2.choose_record(c2.read_records(AUD),'LOCAL_HALL',1,1)
V={s:np.zeros((5,5),complex) for s in (1,2)}
for (s,i,j),e in rec['block'].items():V[s][i-1,j-1]=e['correction']
if rec['symbol']!='Ni-imp' or int(rec['ptype'])!=2 or len(rec['block'])!=50:raise ValueError('LOCAL_HALL binding mismatch')
vh_v={1:VH[4:9,4:9].copy(),2:VH[13:18,13:18].copy()}
vhub_matches={str(s):float(np.max(np.abs(V[s]-vh_v[s]))) for s in (1,2)}
# q0 comparison captured P trace vs radial P_OUT vs replay live QL vs input payload radial q0
spincheck={}
green_vin_checks={}
for s in (1,2):
 q0p=float(np.trace(P[s]).real);q0g=float(np.trace(Gmom[s][0]).real)
 q0rad=rad[s][0]['scalar'];q0live=liveq[s][(2,0)];q0m0=m0[s][0]['scalar']
 pG=Gmom[s][0];m1G=Gmom[s][1];m2G=Gmom[s][2]
 c_raw=float(np.trace(m1G).real/q0g)
 q2raw=float(np.trace(m2G).real-2*c_raw*np.trace(m1G).real+c_raw*c_raw*q0g)
 tv=float(np.trace(V[s]@pG).real)
 q2local=float(q2raw-2*np.trace(V[s]@(m1G-c_raw*pG)).real+
               np.trace(V[s]@V[s]@pG).real-tv*tv/q0g)
 other=V[3-s]
 tvswap=float(np.trace(other@pG).real)
 q2swap=float(q2raw-2*np.trace(other@(m1G-c_raw*pG)).real+
               np.trace(other@other@pG).real-tvswap*tvswap/q0g)
 green_vin_checks[str(s)]={'q0':q0g,'center_raw_from_M1':c_raw,'M2_trace':float(np.trace(m2G).real),
   'q2_raw_from_Green_M2':q2raw,'Tr_VP_true':tv,'center_after_true_V':c_raw-tv/q0g,
   'q2_after_true_V':q2local,'Tr_Vopposite_Psame':tvswap,
   'center_after_swapped_V':c_raw-tvswap/q0g,'q2_after_swapped_V':q2swap,
   'P_OUT_center_crosscheck':rad[s][1]['energy'],'P_OUT_q2_crosscheck':rad[s][2]['energy']-rad[s][1]['energy']**2,
   'payload_true_local_q2':preg['target_operands']['q2_prime_by_spin_ry2'][s-1],
   'payload_opposite_spin_q2':preg['target_operands']['spin_swap_q2_prime_ry2'][s-1]}
 spincheck[str(s)]={
  'q0_P_OUT_matrix_trace':q0p,'q0_RS_GREEN_matrix_trace':q0g,'q0_radial_moment_P_OUT':q0rad,
  'q0_radial_replay_QL':q0live,'q0_m0_carrier_input':q0m0,
  'P_vs_Green_max_abs':float(np.max(np.abs(P[s]-Gmom[s][0]))),
  'M1_vs_Green_max_abs':float(np.max(np.abs(M1[s]-Gmom[s][1]))),
  'q0_P_to_radial_abs':abs(q0p-q0rad),'q0_P_to_live_QL_abs':abs(q0p-q0live),
  'q0_P_to_m0_abs':abs(q0p-q0m0),
  'raw_M2_trace':float(np.trace(Gmom[s][2]).real),
  'P_OUT_q2':rad[s][2]['energy']-rad[s][1]['energy']**2,
  'Green_q2_raw':float(np.trace(Gmom[s][2]).real - 2*(np.trace(Gmom[s][1]).real/q0g)*np.trace(Gmom[s][1]).real + (np.trace(Gmom[s][1]).real/q0g)**2*q0g),
  'P_OUT_after_local_center':rad[s][1]['energy'],'P_OUT_after_local_q2':rad[s][2]['energy']-rad[s][1]['energy']**2,
  'live_q2':liveq[s][(2,2)], 'm0_type_writer_metadata':m0[s][0]['type']}
# operator and spin-swap exact traces
traces={}
for s in (1,2):
 traces[str(s)]={'Tr_VsPs':float(np.trace(V[s]@P[s]).real),
                 'Tr_Vopposite_Psame':float(np.trace(V[3-s]@P[s]).real),
                 'V_frobenius':float(np.linalg.norm(V[s],'fro')),
                 'V_max_abs':float(np.max(np.abs(V[s]))),
                 'P_trace':float(np.trace(P[s]).real),
                 'Vin_antihermitian_max':float(np.max(np.abs(V[s]-V[s].conj().T)))}
Vdiff=V[1]-V[2]
# q1/traces if true/swap cross products; real and complex
traces['Vspin_difference']={'frobenius_norm':float(np.linalg.norm(Vdiff,'fro')),
                            'max_abs_entry':float(np.max(np.abs(Vdiff))),
                            'relative_to_mean_V_frobenius':float(np.linalg.norm(Vdiff,'fro')/(0.5*(np.linalg.norm(V[1],'fro')+np.linalg.norm(V[2],'fro'))))}
traces['total_true_Tr_VP']=sum(traces[str(s)]['Tr_VsPs'] for s in (1,2))
traces['total_swapped_Tr_VP']=sum(traces[str(s)]['Tr_Vopposite_Psame'] for s in (1,2))
traces['total_swap_delta']=traces['total_swapped_Tr_VP']-traces['total_true_Tr_VP']
traces['per_spin_true_minus_swapped_TrVP']={str(s):traces[str(s)]['Tr_VsPs']-traces[str(s)]['Tr_Vopposite_Psame'] for s in (1,2)}
# Parse exact Fortran driver input payload: 501 radial lines then seeds/PL/case fields.
lines=INPUT.read_text().splitlines()
head=list(map(int,lines[0].split())); nr,ncore,nval,ncase=head[1:5]
if head != [2,501,10,6,5,0,1] or (nr,ncore,nval,ncase)!=(501,10,6,5):raise ValueError(f'unexpected payload header {head}')
idx=3+nr+2
pl_payload={s:{ell:[] for ell in (0,1,2)} for s in (1,2)}
for s in (1,2):
 for ell in (0,1,2):
  vals=list(map(float,lines[idx].split()));idx+=1
  pl_payload[s][ell]=vals
  if len(vals)!=6:raise ValueError(f'bad PL input width spin={s} ell={ell}')
  expected=[plmap[s][ell][o] for o in range(4)]+[plmap[s][ell][4],pllive[s][ell]]
  if ell==2:
   raw_center=float(np.trace(Gmom[s][1]).real/np.trace(Gmom[s][0]).real)
   expected[0]=plmap[s][ell][0]+raw_center-m0[s][0]['energy']
  if max(abs(a-b) for a,b in zip(vals,expected))>1e-14:
   raise ValueError(f'PL input/capture mismatch spin={s} ell={ell}: {vals} vs {expected}')
parsed_cases=[]
case_names=['RAW_PL','LOCAL_PL','SIGN_WRONG','OPPOSITE_SPIN','SCALAR_V_PLUS_001']
for cname in case_names:
 sh=list(map(float,lines[idx].split()[:2])); calc=lines[idx].split()[2]; idx+=1
 ql=np.zeros((3,3,2))
 for s in range(2):
  for ell in range(3):
   for m in range(3): ql[m,ell,s]=float(lines[idx]);idx+=1
 parsed_cases.append({'name':cname,'shift_by_spin_ry':sh,'calc_flag':calc,
  'payload_d_q0':[float(lines[0].split()[0])]*0,
  'payload_d_q2_by_spin':[ql[2,2,s] for s in range(2)],
  'payload_d_q0_by_spin':[ql[0,2,s] for s in range(2)],
  'payload_d_q1_by_spin':[ql[1,2,s] for s in range(2)]})
if idx!=len(lines):raise ValueError(f'payload decoder leftover lines {len(lines)-idx}')
# Compare actual input payload q2 and shifts with prereg report definitions.
expected_shifts={x['name']:x['shift_by_spin_ry'] for x in preg['case_definitions']}
for c in parsed_cases:
 if c['shift_by_spin_ry']!=expected_shifts[c['name']]:raise ValueError(f"payload shift mismatch {c['name']}")
 c['registered_shift_by_spin_ry']=expected_shifts[c['name']]
# q0 matching output P_OUT/m0 and input QL values; all three source rows are separately reported.
for c in parsed_cases:
 c['q0_d_payload_vs_P_OUT_abs_by_spin']=[abs(a-b) for a,b in zip(c['payload_d_q0_by_spin'],[rad[s][0]['scalar'] for s in (1,2)])]
 c['q0_d_payload_vs_saved_QL_abs_by_spin']=[abs(a-liveq[s][(2,0)]) for a,s in zip(c['payload_d_q0_by_spin'],(1,2))]
 c['q1_d_payload_vs_saved_QL_abs_by_spin']=[abs(a-liveq[s][(2,1)]) for a,s in zip(c['payload_d_q1_by_spin'],(1,2))]
# Independently bind payload ql arrays to stored radial QL, applying only preregistered q2 values.
payload_ql_checks=[]
for ci,c in enumerate(parsed_cases):
  shift=c['shift_by_spin_ry']
  expected_shift=expected_shifts[c['name']]
  if shift!=expected_shift or c['calc_flag']!='.true.':raise ValueError('case header mismatch')
  # q0/q1 and non-d q2 are unchanged from stored QL; d q2 equals the locked
  # local q2 except the opposite-spin case, whose q2 comes from swapped Vin.
  for s in (1,2):
   for ell in (0,1,2):
    # ql payload was already decoded; recover the same slot by the known order.
    start=512+ci*19+1+(s-1)*9+ell*3
    qvals=[float(lines[start+j]) for j in range(3)]
    if abs(qvals[0]-liveq[s][(ell,0)])>1e-14 or abs(qvals[1]-liveq[s][(ell,1)])>1e-14:
     raise ValueError(f'payload changed saved q0/q1 unexpectedly: {c["name"]} s{s} l{ell}')
    if ell<2 and abs(qvals[2]-liveq[s][(ell,2)])>1e-14:
     raise ValueError(f'payload changed non-d q2 unexpectedly: {c["name"]} s{s} l{ell}')
   target_q2=(preg['target_operands']['spin_swap_q2_prime_ry2'][s-1]
              if c['name']=='OPPOSITE_SPIN' else preg['target_operands']['q2_prime_by_spin_ry2'][s-1])
   actual_q2=c['payload_d_q2_by_spin'][s-1]
   if abs(actual_q2-target_q2)>1e-14:raise ValueError(f'payload d q2 mismatch case={c["name"]} spin={s}')
  payload_ql_checks.append({'case':c['name'],'d_q2_payload':c['payload_d_q2_by_spin'],'registered_expected_d_q2':[preg['target_operands']['spin_swap_q2_prime_ry2'][s-1] if c['name']=='OPPOSITE_SPIN' else preg['target_operands']['q2_prime_by_spin_ry2'][s-1] for s in (1,2)],'q0_q1_and_non_d_q2_match_saved_QL':True})
# Case output energies/outputs can reconstruct observed carrier deltas in report; bind case-specific values.
report=read_json(BASE/'cases_005/rd2_cases_report.json')
measure=report['measurements']
# payload SHA and first run artifacts
out={
 'schema':'c5-rd2-debugging-attempt-payload-phase-spin-audit-v1',
 'status':'PREPARED_READ_ONLY_DIAGNOSTIC; no replay of helper yet',
 'hashes':{str(p):sha(p) for p in [CSV,AUD,PARSER,MATH,BIND,PREG,INPUT,FIRST,BASE/'cases_005/rd2_cases_report.json']},
 'binding_runtime_map':binding['runtime_map'],'capture_dimensions':binding['capture_completeness'],
 'spin_matrix_checks':spincheck,
 'independent_Green_P_M1_M2_and_Vin_moment_reconstruction':green_vin_checks,
 'local_hall_meta':{'symbol':rec['symbol'],'physical_type':int(rec['ptype']),'site':1,'entry_count':len(rec['block']), 'correction_source':'exact after-before additive LOCAL_HALL BLOCK'},
 'LOCAL_HALL_vs_v_hub_type2_P_insertion_max_abs_by_spin':vhub_matches,
 'Vin_spin_traces_and_identifiability':traces,
 'payload_decode':{'path':str(INPUT),'sha256':sha(INPUT),'header':head,'line_count':len(lines),'pl_input_by_spin_ell':{str(s):{str(ell):vals for ell,vals in es.items()} for s,es in pl_payload.items()},'case_ql_input_checks':payload_ql_checks,'cases':parsed_cases},
 'formal_case_summary':{'status':report['status'],'carrier':measure['carrier'],'scalar_plus_001':measure['scalar_plus_001'],'wrong_controls':measure['wrong_controls'],'formal_gates':report['formal_gates']},
 'notes':['P_OUT P/M1, Green P/M1/M2, radial QL, and payload QL are checked separately; the latter can be a different phase.','Opposite-spin binding uses V from LOCAL_HALL blocks keyed by spin and P/M1/Green matrices keyed by captured spin; no renaming/reordering occurs.','Reproduction command is to be run only after independent preflight audit and root authorization.']}
path=Path('/private/tmp/c52-rd2-private/diagnostic_006/payload_phase_spin_audit.json')
path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'output':str(path),'sha256':sha(path),'green_blocks':len(gblocks),'payload_cases':len(parsed_cases),'Vin_norm_diff':traces['Vspin_difference']['frobenius_norm'],'true_TrVP':traces['total_true_Tr_VP'],'swapped_TrVP':traces['total_swapped_Tr_VP'],'q0_checks':{s:{k:v for k,v in z.items() if k.startswith('q0_')} for s,z in spincheck.items()}},indent=2))
