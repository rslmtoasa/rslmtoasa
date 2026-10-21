#!/usr/bin/env python3
"""C5.4 run staging/preflight; copies only explicitly whitelisted inputs."""
import argparse, json, os, shutil, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, "/private/tmp/c52-c-private")
from fortran_namelist_schema import source_schema, validate_input
import f90nml

THREADS = {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
 "RSLMTO_OMP_THREADS_SERIAL": "1", "OPENBLAS_NUM_THREADS": "1",
 "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
 "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
SCHEMA = source_schema(REPO)

def stage(run, sources):
    if run.exists(): raise SystemExit(f"fresh run directory already exists: {run}")
    run.mkdir(parents=True)
    for src in sources:
        src = Path(src).resolve()
        if not src.is_file(): raise SystemExit(f"missing whitelist input: {src}")
        shutil.copy2(src, run / src.name)
    for p in run.iterdir(): p.chmod(p.stat().st_mode | 0o200)

def check(run, deck, restart, guard):
    if not deck.is_file() or not os.access(deck, os.W_OK): raise SystemExit("working deck missing/not writable")
    validate_input(deck, SCHEMA)
    nml = f90nml.read(str(deck))
    ctl = nml.get("control", {}); selfn = nml.get("self", {})
    if ctl.get("recur") not in ("block", "cheby"): raise SystemExit("recur missing/unsupported")
    if "use_kspace" not in selfn: raise SystemExit("self.use_kspace missing")
    if restart and (not restart.is_file() or not os.access(restart, os.W_OK)):
        raise SystemExit("restart source missing/not writable working copy")
    if not run.is_dir() or not os.access(run, os.W_OK): raise SystemExit("run directory missing/not writable")
    if guard:
        if not guard.is_file(): raise SystemExit("launcher missing")
        text = guard.read_text()
        if not text.startswith("#!/bin/bash"):
            raise SystemExit("launcher must be guarded by a bash script")
        if any(k not in text for k in THREADS): raise SystemExit("launcher missing thread cap")
    if restart: validate_input(restart, SCHEMA)
    elem = f90nml.read(str(restart)) if restart else {}
    if restart and (elem.get("element", {}).get("num_quant_d", 0) <= 0 or
                    elem.get("par", {}).get("lmax", -1) < 2):
        raise SystemExit("Ni d orbital channel absent from actual restart source")
    print(json.dumps({"status":"PASS", "deck":str(deck), "recur":ctl["recur"],
      "use_kspace":selfn["use_kspace"], "restart_source":str(restart) if restart else None,
      "available_orbital_channels": ["d"] if restart else [],
      "fresh_run_dir":str(run), "thread_caps":THREADS}, indent=2))

if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("run_dir", type=Path); ap.add_argument("--stage", nargs="+", type=Path)
    ap.add_argument("--deck", type=Path); ap.add_argument("--restart", type=Path); ap.add_argument("--guard", type=Path)
    a=ap.parse_args()
    if a.stage: stage(a.run_dir, a.stage)
    if a.deck: check(a.run_dir, a.deck, a.restart, a.guard)
    elif not a.stage: ap.error("provide --stage whitelist files and/or --deck")
