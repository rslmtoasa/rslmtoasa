#!/usr/bin/env python3
"""C5.4 run staging/preflight; copies only explicitly whitelisted inputs."""
import argparse, csv, hashlib, json, os, shutil, sys, f90nml
from pathlib import Path
REPO = Path(__file__).resolve().parents[4]
SCHEMA_MODULE_DIR = REPO / "docs/validation/fixed_ldau/rs/C5_2_PART_C_assets/C3/scripts"
sys.path.insert(0, str(SCHEMA_MODULE_DIR))
from fortran_namelist_schema import source_schema
THREADS = {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
 "RSLMTO_OMP_THREADS_SERIAL": "1", "OPENBLAS_NUM_THREADS": "1",
 "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
 "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
SCHEMA_ROOT=Path(os.environ.get("C5_4_PREFLIGHT_SOURCE_ROOT", REPO)); SCHEMA = source_schema(SCHEMA_ROOT)
def validate(nml, path):
    ignored=[]
    for group, values in nml.items():
        group=str(group).lower()
        if group not in SCHEMA: ignored.append(group); continue  # solver namelist reads leave undeclared groups unconsumed
        bad=set(map(str.lower, values))-SCHEMA[group]
        if bad: raise SystemExit(f"{path}: unknown {group} members {sorted(bad)}")
    return sorted(ignored)
def stage(run, sources):
    if run.exists(): raise SystemExit(f"fresh run directory already exists: {run}")
    run.mkdir(parents=True)
    for src in sources:
        src = Path(src).resolve()
        if not src.is_file(): raise SystemExit(f"missing whitelist input: {src}")
        shutil.copy2(src, run / src.name)
    for p in run.iterdir(): p.chmod(p.stat().st_mode | 0o200)
def check_host_fermi(nml, trace_path, tolerance=1e-6):
    """Compare input &energy fermi to the unique fermi_level row in host trace."""
    energy = nml.get("energy", {})
    if "fermi" not in energy: raise SystemExit("input &energy.fermi missing for host-EF check")
    host_ef = []
    with trace_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("kind", "").strip() == "fermi_level":
                host_ef.append(float(row["scalar"]))
    if len(host_ef) != 1: raise SystemExit(f"{trace_path}: expected one fermi_level row, found {len(host_ef)}")
    delta = float(energy["fermi"]) - host_ef[0]
    if abs(delta) > tolerance:
        raise SystemExit(f"input fermi {energy['fermi']} differs from host EF {host_ef[0]} by {delta} Ry (limit {tolerance})")
    return {"status": "PASS", "host_ef_trace": str(trace_path), "host_ef": host_ef[0],
      "input_fermi": float(energy["fermi"]), "delta_ry": delta, "tolerance_ry": tolerance}
def check(run, deck, restart, guard, host_ef_trace=None):
    if not deck.is_file() or not os.access(deck, os.W_OK): raise SystemExit("working deck missing/not writable")
    nml = f90nml.read(str(deck))
    ignored=validate(nml, deck)
    host_ef_check = check_host_fermi(nml, host_ef_trace) if host_ef_trace else None
    ctl = nml.get("control", {}); selfn = nml.get("self", {})
    if ctl.get("recur") not in ("block", "cheby"): raise SystemExit("recur missing/unsupported")
    use_kspace=selfn.get("use_kspace", False); default_source="explicit"
    if "use_kspace" not in selfn:
        lines=(REPO/"source/self.f90").read_text().splitlines()
        line=next((i for i,s in enumerate(lines,1) if "this%use_kspace = .false." in s),0)
        if not line: raise SystemExit("cannot resolve self.use_kspace source default")
        default_source=f"source/self.f90:{line}"
    if restart and (not restart.is_file() or not os.access(restart, os.W_OK)):
        raise SystemExit("restart source missing/not writable working copy")
    if not run.is_dir() or not os.access(run, os.W_OK): raise SystemExit("run directory missing/not writable")
    if guard:
        if not guard.is_file(): raise SystemExit("launcher missing")
        text = guard.read_text()
        if not text.startswith("#!/bin/bash"):
            raise SystemExit("launcher must be guarded by a bash script")
        if any(k not in text for k in THREADS): raise SystemExit("launcher missing thread cap")
    if restart: ignored+=validate(f90nml.read(str(restart)), restart)
    elem = f90nml.read(str(restart)) if restart else {}
    if restart and (elem.get("element", {}).get("num_quant_d", 0) <= 0 or elem.get("par", {}).get("lmax", -1) < 2):
        raise SystemExit("Ni d orbital channel absent from actual restart source")
    print(json.dumps({"status":"PASS", "deck":str(deck), "recur":ctl["recur"],
      "use_kspace":use_kspace, "use_kspace_source":default_source,
      "unconsumed_groups":sorted(set(ignored)), "schema_source":str(SCHEMA_ROOT), "restart_source":str(restart) if restart else None,
      "host_ef_check":host_ef_check,
      "available_orbital_channels": ["d"] if restart else [],
      "fresh_run_dir":str(run), "thread_caps":THREADS}, indent=2))
def offline(run, payload, base, driver, binding):
    if run.exists(): raise SystemExit("offline output directory must be fresh")
    b=json.loads(binding.read_text())
    if hashlib.sha256(driver.read_bytes()).hexdigest()!=b["driver"]["sha256"]:
        raise SystemExit("offline driver hash differs from frozen binding")
    src=base.read_text().splitlines(); dst=payload.read_text().splitlines()
    h=list(map(int,src[0].split())); start=3+h[1]+2+6
    sh,dh=src[0].split(),dst[0].split()
    if sh[:4]+sh[5:]!=dh[:4]+dh[5:] or dst[1:start]!=src[1:start] or int(dh[4])!=8 or len(dst)!=start+8*19:
        raise SystemExit("offline scalar scan changed frozen payload prefix/dimensions")
    if any(os.environ.get(k)!=v for k,v in THREADS.items()): raise SystemExit("offline thread caps missing")
    run.mkdir(parents=True); shutil.copy2(payload,run/"driver.input.txt")
    print(json.dumps({"status":"PASS_OFFLINE","driver_sha256":b["driver"]["sha256"],
      "fresh_run_dir":str(run),"thread_caps":THREADS}))
if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("run_dir", type=Path); ap.add_argument("--stage", nargs="+", type=Path)
    ap.add_argument("--deck", type=Path); ap.add_argument("--restart", type=Path); ap.add_argument("--guard", type=Path)
    ap.add_argument("--host-ef-trace", type=Path, help="host fixed_u_trace.csv with one fermi_level row")
    ap.add_argument("--offline", nargs=4, metavar=("PAYLOAD","BASE","DRIVER","BINDING"), type=Path)
    a=ap.parse_args()
    if a.stage: stage(a.run_dir, a.stage)
    if a.deck: check(a.run_dir, a.deck, a.restart, a.guard, a.host_ef_trace)
    if a.offline: offline(a.run_dir, *a.offline)
    elif not a.stage and not a.deck: ap.error("provide --stage whitelist files and/or --deck")
