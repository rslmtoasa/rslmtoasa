from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.util, json
adapter_path=Path('/private/tmp/c52-rd2-private/nr501_adapter/rd2_replay_adapter_nr501.py')
spec=importlib.util.spec_from_file_location('rd2_adapter_preflight',adapter_path)
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
args=SimpleNamespace(mode='baseline',capture_dir=Path('/private/tmp/c52-rd2-ni-u3j07-capture/run'),binding=Path('/private/tmp/c52-rd2-private/vin_binding_nr501.json'),driver=Path('/private/tmp/c52-rd2-private/nr501_adapter/build_proven/accepted_m03_replay_driver_nr501'),authorization=Path('/private/tmp/c52-rd2-private/execution_authorization.json'),output_dir=Path('/private/tmp/c52-rd2-private/baseline_preflight_smoke_output'))
assert not args.output_dir.exists()
assert sha(args.driver)==a.load_json(args.binding)['driver']['sha256']
called=[]
class StopBeforeDriver(RuntimeError): pass
def recorder(driver,payload,outdir):
    called.append({'driver':str(driver),'driver_sha256':sha(driver),'output_dir':str(outdir),'payload_sha256':hashlib.sha256(payload.encode()).hexdigest(),'payload_bytes':len(payload.encode())})
    raise StopBeforeDriver('intentional smoke barrier: no executable launch')
runner,c2=a.modules(); runner.run_persistent=recorder
a.modules=lambda:(runner,c2)
try:
    a.baseline(args)
except StopBeforeDriver:
    pass
else:
    raise SystemExit('smoke did not reach driver barrier')
assert len(called)==1
r=called[0]
report={'schema':'c5-rd2-baseline-nonnumerical-preflight-smoke-v1','status':'PASS_SETUP_REACHED_DRIVER_BARRIER_NO_EXECUTION','plan_sha256':sha(Path('/private/tmp/c52-rd2-private/rd2_replay_plan.json')),'adapter_sha256':sha(adapter_path),'binding_sha256':sha(args.binding),'driver_sha256':sha(args.driver),'authorization_sha256':sha(args.authorization),'capture_csv_sha256':sha(args.capture_dir/'records.csv'),'NR':501,'lmax':2,'ncore':10,'nval':6,'baseline_prepare':'completed using exact derived runner','text_input':'generated and intercepted before process launch','driver_invocations':0,'intercepted_call':r,'smoke_output_dir_created':args.output_dir.is_dir(),'scope':'authorization, capture/binding hashes, full capture parsing, baseline prepare, text_input generation; run_persistent replaced by recorder raising before subprocess; no radial driver/helper, no preregistration, no cases, no SCF solver'}
out=Path('/private/tmp/c52-rd2-private/nr501_adapter/baseline_preflight_smoke.json');out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,sort_keys=True))
