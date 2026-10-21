#!/usr/bin/env python3
"""Read-only iteration collector for the C5.3 Ni U0 bulk host run."""
import argparse, csv, datetime, hashlib, json, os, pathlib, signal, sys, time
p=argparse.ArgumentParser()
p.add_argument('--run-id',required=True)
g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--solver-pid',type=int)
g.add_argument('--solver-pid-file')
p.add_argument('--run-dir',required=True)
p.add_argument('--state-file',required=True)
p.add_argument('--output-dir',required=True)
p.add_argument('--ready-file',required=True)
p.add_argument('--poll-seconds',type=float,default=.05)
a=p.parse_args()
run=pathlib.Path(a.run_dir).resolve(); out=pathlib.Path(a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
trace=run/'fixed_u_trace.csv'; state=run/a.state_file; manifest=out/'trace_state_manifest.jsonl'; ready=pathlib.Path(a.ready_file).resolve(); stop=False
pidfile=pathlib.Path(a.solver_pid_file).resolve() if a.solver_pid_file else None
solver_pid=a.solver_pid
def sha(b): return hashlib.sha256(b).hexdigest()
def on_signal(sig,frame):
 global stop
 stop=True
signal.signal(signal.SIGTERM,on_signal); signal.signal(signal.SIGINT,on_signal)
# Initialize all collector state and record startup before readiness is signaled.
manifest.touch(exist_ok=True)
started={'record_type':'read_only_trace_watcher_start','run_id':a.run_id,'solver_pid':solver_pid,'solver_pid_file':str(pidfile) if pidfile else None,'watcher_process_pid':os.getpid(),'run_dir':str(run),'trace_source':str(trace),'state_source':str(state),'state_file':a.state_file,'output_dir':str(out),'manifest_path':str(manifest),'ready_file':str(ready),'actual_argv':sys.argv,'started_local':datetime.datetime.now().astimezone().isoformat()}
(out/'watcher_started.json').write_text(json.dumps(started,indent=2,sort_keys=True)+'\n')
with ready.open('w') as f: f.write('ready\n'); f.flush(); os.fsync(f.fileno())
while not stop and solver_pid is None:
 try:
  solver_pid=int(pidfile.read_text().strip())
  (out/'solver_attached.json').write_text(json.dumps({'run_id':a.run_id,'solver_pid':solver_pid,'attached_local':datetime.datetime.now().astimezone().isoformat()},indent=2,sort_keys=True)+'\n')
 except (OSError,ValueError): time.sleep(a.poll_seconds)
seen=set()
while not stop and solver_pid is not None:
 try: os.kill(solver_pid,0); alive=True
 except ProcessLookupError: alive=False
 except PermissionError: alive=True
 if not alive: break
 try:
  t1=trace.read_bytes(); time.sleep(.01); t2=trace.read_bytes(); s1=state.read_bytes(); time.sleep(.01); s2=state.read_bytes()
  if not t1 or t1!=t2 or s1!=s2: time.sleep(a.poll_seconds); continue
  lines=t2.decode().splitlines()
  if not lines or not lines[0].startswith('kind,run_id,'): time.sleep(a.poll_seconds); continue
  rows=list(csv.DictReader(lines)); ids={r['iteration'].strip() for r in rows}
  if len(ids)!=1 or any(r['run_id'].strip()!=a.run_id for r in rows): time.sleep(a.poll_seconds); continue
  post=[r for r in rows if r['stage'].strip()=='POST_SCF_RADIAL']; pql=[r for r in rows if r['stage'].strip()=='P_OUT' and r['kind'].strip()=='ql']
  if len(post)!=6 or post[-1]['kind'].strip()!='etot_total' or len(pql)!=6: time.sleep(a.poll_seconds); continue
  ls={(int(r['l']),int(r['spin'])) for r in pql}
  if ls!={(l,s) for l in range(3) for s in (1,2)}: time.sleep(a.poll_seconds); continue
  it=next(iter(ids)); th=sha(t2); sh=sha(s2); key=(it,th,sh)
  if key in seen: time.sleep(a.poll_seconds); continue
  stem=f'{a.run_id}_traceiter{int(float(it)):04d}'
  td=out/(stem+'_'+th[:12]+'.csv'); sd=out/(stem+'_'+state.stem+'_'+sh[:12]+'.nml')
  if not td.exists(): td.write_bytes(t2)
  if not sd.exists(): sd.write_bytes(s2)
  rec={'run_id':a.run_id,'solver_pid':solver_pid,'trace_iteration':it,'trace_sha256':th,'trace_bytes':len(t2),'trace_snapshot':str(td),'trace_source':str(trace),'state_file':a.state_file,'state_sha256':sh,'state_bytes':len(s2),'state_snapshot':str(sd),'state_source':str(state),'state_source_mtime_ns':state.stat().st_mtime_ns,'captured_local':datetime.datetime.now().astimezone().isoformat(),'stable_double_reads':True,'post_scf_radial_count':len(post),'post_scf_radial_last_kind':post[-1]['kind'].strip(),'P_OUT_ql_count':len(pql),'P_OUT_ql_keys':sorted([f"l{r['l'].strip()}_spin{r['spin'].strip()}" for r in pql]),'paired_state_ql_phase':'P_MIX (following mixed state; not P_OUT)'}
  with manifest.open('a') as f: f.write(json.dumps(rec,sort_keys=True)+'\n'); f.flush(); os.fsync(f.fileno())
  seen.add(key); print(json.dumps(rec,sort_keys=True),flush=True)
 except (OSError,UnicodeDecodeError,ValueError,KeyError) as e:
  print(json.dumps({'watcher_read_error':type(e).__name__,'message':str(e),'time':datetime.datetime.now().astimezone().isoformat()}),flush=True)
 time.sleep(a.poll_seconds)
ended={'record_type':'read_only_trace_watcher_stop','run_id':a.run_id,'solver_pid':solver_pid,'stop_requested':stop,'ended_local':datetime.datetime.now().astimezone().isoformat(),'manifest_path':str(manifest),'captured_entries':sum(1 for _ in manifest.open()) if manifest.exists() else 0}
(out/'watcher_finished.json').write_text(json.dumps(ended,indent=2,sort_keys=True)+'\n'); print(json.dumps(ended,sort_keys=True),flush=True)
