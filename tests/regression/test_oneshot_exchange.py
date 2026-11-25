#!/usr/bin/env python3
"""One-shot inter-site recursion vs the four-recursions-per-pair route on bcc Fe (J_ij).

Runs intersite_mode = 'pairs' (serial), 'one_shot' (serial) and 'one_shot' (MPI) on the
triad_bccFe_exchange fixture and requires J_ij to agree to rtol 1e-8 (absolute floor 1e-10 mRy).
"""
import argparse
import sys
import time
from pathlib import Path
import subprocess

import f90nml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_triad import extract_jij, setup_workdir

RTOL, ATOL = 1e-8, 1e-10
FIXTURE = Path(__file__).resolve().parent / "triad_bccFe_exchange"


def run(binary, workdir, mode, mpi_exec=None, ranks=1):
    setup_workdir(FIXTURE, workdir)
    tmp = workdir / "input.nml.tmp"
    f90nml.patch(str(workdir / "input.nml"), {"control": {"intersite_mode": mode}}, str(tmp))
    tmp.replace(workdir / "input.nml")
    command = [str(binary)] if mpi_exec is None else [mpi_exec, "-np", str(ranks), str(binary)]
    t0 = time.time()
    result = subprocess.run(command, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)
    elapsed = time.time() - t0
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(f"ERROR: {mode} ranks={ranks} exited {result.returncode}")
    return extract_jij(workdir, result.stdout), elapsed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--mpi-exec", required=True)
    p.add_argument("--ranks", type=int, default=2)
    a = p.parse_args()
    ref, _ = run(a.binary.resolve(), a.workdir.resolve() / "pairs", "pairs")
    cases = {"one_shot serial": run(a.binary.resolve(), a.workdir.resolve() / "one_shot", "one_shot")[0],
             f"one_shot np={a.ranks}": run(a.binary.resolve(), a.workdir.resolve() / "one_shot_mpi", "one_shot",
                                           a.mpi_exec, a.ranks)[0]}
    failed = False
    for name, jij in cases.items():
        for pair, ref_val in ref.items():
            diff = abs(jij[pair] - ref_val)
            ok = diff <= ATOL + RTOL * abs(ref_val)
            failed |= not ok
            print(f"{name} pair {pair}: J={jij[pair]:.12e} pairs={ref_val:.12e} |diff|={diff:.3e} {'ok' if ok else 'FAIL'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
