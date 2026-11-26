#!/usr/bin/env python3
"""One-shot inter-site recursion on the complex contour (PAOFLOW bccFe exchange_p2rs fixture).

(i) recur = chebyshev, lld = 150: legacy ('pairs') vs one_shot J_ij and D_ij, rtol 1e-8 (abs floor 1e-10).
(ii) recur = block, one_shot: serial vs mpirun -np N, same tolerance.
(iii) recur = block, lld = 21, sym_term = F, self pair (1,1): legacy ('pairs') vs one_shot J_11 and D_11, same tolerance.
Haydock one_shot vs legacy for inter-site pairs is reported by hand-run measurements, not checked here.
--gpu (needs a CUDA device, else exit 77): same fixture, lld = 21, sym_term = F, pairs plus the self pair (1,1):
(iv) one_shot block with gpu_plugin = T at gpu_precision = fp64 vs one_shot block on the CPU, all pairs, same tolerance;
     J_11 and D_11 also vs legacy block on the CPU, same tolerance.
(v) the same at gpu_precision = fp32: differences are printed only.
"""
import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import f90nml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_triad import setup_workdir

RTOL, ATOL = 1e-8, 1e-10
CASES = Path(__file__).resolve().parents[1] / "postproc" / "cases" / "paoflow"


def run(binary, workdir, patch, mpi_exec=None, ranks=1, self_pair=False):
    setup_workdir(CASES / "bccFe_exchange", workdir)
    shutil.copy2(CASES / "paoham.dat", workdir / "paoham.dat")
    patch = {"lattice": {"strux_backend": "strux_lib"}, **patch}
    patch["control"] = {"cheb_backend": "legacy", **patch["control"]}
    if self_pair:
        text = (workdir / "input.nml").read_text()
        text = text.replace("njij = 2", "njij = 3").replace("ijpair(2, :) = 1, 2635\n", "ijpair(2, :) = 1, 2635\nijpair(3, :) = 1, 1\n")
        (workdir / "input.nml").write_text(text)
    tmp = workdir / "input.nml.tmp"
    f90nml.patch(str(workdir / "input.nml"), patch, str(tmp))
    tmp.replace(workdir / "input.nml")
    command = [str(binary)] if mpi_exec is None else [mpi_exec, "-np", str(ranks), str(binary)]
    t0 = time.time()
    result = subprocess.run(command, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)
    elapsed = time.time() - t0
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(f"ERROR: {patch['control']} ranks={ranks} exited {result.returncode}")
    values = {}
    for m in re.finditer(r"([JD])ij between pair\s+(\d+)\s+and\s+(\d+)\s+is\s+(.*)", result.stdout):
        for c, x in enumerate(m.group(4).split()):
            values[f"{m.group(1)}ij {m.group(2)}-{m.group(3)}[{c + 1}]"] = float(x)
    n_expected = 12 if self_pair else 8
    if len(values) != n_expected:
        raise SystemExit(f"ERROR: expected {n_expected} J/D values in the log, found {len(values)}")
    return values, elapsed


def compare(label, got, ref):
    failed = False
    for key, ref_val in ref.items():
        diff = abs(got[key] - ref_val)
        ok = diff <= ATOL + RTOL * abs(ref_val)
        failed |= not ok
        print(f"{label} {key}: {got[key]:.12e} ref={ref_val:.12e} |diff|={diff:.3e} {'ok' if ok else 'FAIL'}")
    return failed


def has_cuda_device():
    try:
        out = subprocess.run(["nvidia-smi", "-L"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    except OSError:
        return False
    return out.returncode == 0 and "GPU" in out.stdout


def main_gpu(a):
    if not has_cuda_device():
        print("SKIP: no CUDA device is available")
        return 77
    binary, work = a.binary.resolve(), a.workdir.resolve()
    scf = {"nstep": 1}
    ctl = {"recur": "block", "lld": 21, "nsp": 2, "sym_term": False}
    cpu, _ = run(binary, work / "one_shot_cpu", {"control": {**ctl, "intersite_mode": "one_shot"}, "self": scf}, self_pair=True)
    leg, _ = run(binary, work / "pairs_cpu", {"control": {**ctl, "intersite_mode": "pairs"}, "self": scf}, self_pair=True)
    gpu = {}
    for prec in ("fp64", "fp32"):
        gpu[prec], _ = run(binary, work / f"one_shot_gpu_{prec}",
                           {"control": {**ctl, "intersite_mode": "one_shot", "gpu_plugin": True, "gpu_precision": prec}, "self": scf},
                           self_pair=True)
    is11 = lambda d: {k: v for k, v in d.items() if re.match(r"[JD]ij 1-1\[", k)}
    failed = compare("one_shot GPU fp64 vs one_shot CPU", gpu["fp64"], cpu)
    failed |= compare("one_shot GPU fp64 vs legacy block CPU self pair", is11(gpu["fp64"]), is11(leg))
    for key, ref_val in cpu.items():
        print(f"one_shot GPU fp32 vs one_shot CPU {key}: {gpu['fp32'][key]:.12e} ref={ref_val:.12e} |diff|={abs(gpu['fp32'][key] - ref_val):.3e} (report only)")
    return 1 if failed else 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--mpi-exec")
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--ranks", type=int, default=2)
    a = p.parse_args()
    if a.gpu:
        return main_gpu(a)
    if not a.mpi_exec:
        p.error("--mpi-exec is required without --gpu")
    binary, work = a.binary.resolve(), a.workdir.resolve()
    cheb = {"recur": "chebyshev", "lld": 150, "nsp": 2}
    block = {"recur": "block", "lld": 21, "nsp": 2, "intersite_mode": "one_shot"}
    scf = {"nstep": 1}
    self_ref = {"recur": "block", "lld": 21, "nsp": 2, "sym_term": False}
    ref, _ = run(binary, work / "cheb_pairs", {"control": {**cheb, "intersite_mode": "pairs"}, "self": scf})
    one, _ = run(binary, work / "cheb_one_shot", {"control": {**cheb, "intersite_mode": "one_shot"}, "self": scf})
    serial, _ = run(binary, work / "block_serial", {"control": block, "self": scf})
    mpi, _ = run(binary, work / "block_mpi", {"control": block, "self": scf}, a.mpi_exec, a.ranks)
    leg11, _ = run(binary, work / "self_pairs", {"control": {**self_ref, "intersite_mode": "pairs"}, "self": scf}, self_pair=True)
    one11, _ = run(binary, work / "self_one_shot", {"control": {**self_ref, "intersite_mode": "one_shot"}, "self": scf}, self_pair=True)
    failed = compare("chebyshev one_shot vs pairs", one, ref)
    failed |= compare(f"block one_shot np={a.ranks} vs serial", mpi, serial)
    failed |= compare("block one_shot vs pairs self pair", {k: one11[k] for k in one11 if "1-1[" in k}, {k: leg11[k] for k in leg11 if "1-1[" in k})
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
