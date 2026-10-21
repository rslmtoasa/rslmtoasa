#!/usr/bin/env python3
"""Read-only Ni host acceptance scan using the registered C1 limits only."""
import argparse,csv,hashlib,json,pathlib,re
from decimal import Decimal
p=argparse.ArgumentParser(); p.add_argument('--manifest',required=True); p.add_argument('--log',required=True); p.add_argument('--output',required=True); a=p.parse_args()
manifest=pathlib.Path(a.manifest); logp=pathlib.Path(a.log); outp=pathlib.Path(a.output)
entries=[json.loads(x) for x in manifest.read_text().splitlines() if x.strip()]
selected={}
for e in entries:
 if e['run_id']!='rd_ni_i_u3j07_local_from_rs_bulk_u3_rc80': raise SystemExit('wrong run id in collector manifest')
 it=int(e['trace_iteration'])
 if it not in selected: selected[it]=e
log=logp.read_text(errors='replace'); marks=list(re.finditer(r'Perform recursion at step\s+(\d+)',log)); diffs={}
for j,m in enumerate(marks):
 n=int(m.group(1)); end=marks[j+1].start() if j+1<len(marks) else len(log); seg=log[m.end():end]
 found=list(re.finditer(r'(?:Not converged!\s*Diff=\s*|Converged!\s*)([0-9.Ee+-]+)',seg))
 diffs[n]=Decimal(found[0].group(1)) if len(found)==1 else None
values={}
for it,e in selected.items():
 for path_key,hash_key in [('trace_snapshot','trace_sha256'),('state_snapshot','state_sha256')]:
  fp=pathlib.Path(e[path_key]); raw=fp.read_bytes()
  if hashlib.sha256(raw).hexdigest()!=e[hash_key]: raise SystemExit('snapshot hash mismatch '+str(fp))
 raw=pathlib.Path(e['trace_snapshot']).read_bytes(); rows=list(csv.DictReader(raw.decode().splitlines()))
 if {int(r['iteration']) for r in rows}!={it}: raise SystemExit('trace iteration mismatch')
 post=[r for r in rows if r['stage']=='POST_SCF_RADIAL']; ql=[r for r in rows if r['stage']=='P_OUT' and r['kind']=='ql']
 if len(post)!=6 or post[-1]['kind']!='etot_total' or len(ql)!=6: raise SystemExit(f'incomplete per-iteration trace families at {it}')
 q={(int(r['l']),int(r['spin'])):Decimal(r['scalar']) for r in ql}
 if set(q)!={(l,s) for l in range(3) for s in (1,2)}: raise SystemExit(f'P_OUT ql map mismatch {it}')
 values[it]={'entry':e,'q0':q,'moment':abs(sum(q[(l,1)]-q[(l,2)] for l in range(3)))}
def span(xs):return max(xs)-min(xs)
results=[]
for endpoint in sorted(selected):
 ids=list(range(endpoint-19,endpoint+1))
 if ids[0]<1 or any(i not in selected or diffs.get(i) is None for i in ids): continue
 w=[values[i] for i in ids]; d1=span([x['q0'][(2,1)] for x in w]); d2=span([x['q0'][(2,2)] for x in w]); m=span([x['moment'] for x in w]); diff=diffs[endpoint]
 ok=diff<=Decimal('3e-8') and d1<=Decimal('5e-7') and d2<=Decimal('5e-7') and m<=Decimal('1e-6')
 results.append({'endpoint_iteration':endpoint,'same_loop_diff':str(diff),'window_iterations':ids,'d_q0_spin1_spread':str(d1),'d_q0_spin2_spread':str(d2),'physical_q0_moment_spread':str(m),'passes_registered_C1_limits':ok,'trace_hashes':[values[i]['entry']['trace_sha256'] for i in ids],'endpoint_state_sha256':values[endpoint]['entry']['state_sha256']})
qualified=[r for r in results if r['passes_registered_C1_limits']]
out={'status':'READ_ONLY_SCAN','run_id':'rd_ni_i_u3j07_local_from_rs_bulk_u3_rc80','capture_rule':'first complete stable collector manifest entry per trace iteration','trace_to_diff_rule':'trace ID i binds only to exactly one Diff inside its Perform recursion at step i segment; missing/multiple rejected','occupancy_phase':'P_OUT fresh/unmixed ql(1,l,spin), l=0,1,2','physical_moment':'abs(sum_l(q0_up-q0_down)), computed as Decimal from emitted P_OUT tokens','registered_limits':{'outer_diff_max':'3e-8','each_d_q0_spin_spread_max':'5e-7','q0_physical_moment_spread_max':'1e-6'},'no_denominator_or_extra_tolerance':True,'captured_iterations':sorted(selected),'diff_complete_ids':sorted(i for i,v in diffs.items() if v is not None),'windows':results,'earliest_qualifying':qualified[0] if qualified else None,'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(logp.read_bytes()).hexdigest()}
outp.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({'status':out['status'],'captured':len(selected),'diff_complete':len(out['diff_complete_ids']),'earliest_qualifying':out['earliest_qualifying'],'output':str(outp)},indent=2))
