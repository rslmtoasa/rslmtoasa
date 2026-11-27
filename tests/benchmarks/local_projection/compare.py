"""Compare frozen Pt/FeCo production spectra and collect per-process resources."""
import argparse
import csv
import os
from pathlib import Path
import sys
import time

import f90nml
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "compact_conductivity"))
from benchmark import prepare, outputs

ROOT = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser()
parser.add_argument("baseline", type=Path)
parser.add_argument("feature", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--orders", type=int, nargs="+", default=[20, 50, 250])
parser.add_argument("--size", type=int, default=4)
parser.add_argument("--repeat", type=int, default=2)
parser.add_argument("--strux", choices=["legacy", "strux_lib"], default="strux_lib")
parser.add_argument("--materials", nargs="+", choices=["Pt", "FeCo"], default=["Pt","FeCo"])
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
fixture = ROOT / "tests/postproc/cases/conductivity/fccPt"
observables = dict(charge="z", spin="z", orbital="z", spin_torque="x",
                   spin_soc_torque="x", spin_accumulation="x", orbital_accumulation="z",
                   quadrupole_accumulation="Qxy", octupole_accumulation="Mx3y2r2")
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
summary = args.output/"summary.csv"
done = set()
if summary.exists():
    with summary.open() as previous:
        done = {(r["material"], int(r["order"]), r["observable"], r["projection"], r["binary"], int(r["repeat"]))
                for r in csv.DictReader(previous)}
with summary.open("a") as f:
    writer = csv.writer(f)
    if not done:
        writer.writerow(["material", "order", "observable", "projection", "binary", "repeat",
                     "seconds", "peak_RSS_KiB", "max_abs_diff", "max_scaled_diff"])
    for material in args.materials:
        for order in args.orders:
            for observable, pol in {**observables, "spin_torque_onsite": "x"}.items():
                cost = observable in ("spin_torque", "spin_torque_onsite")
                if not cost and order != args.orders[0]:
                    continue
                cases = [("baseline", ""), ("feature", ""), ("feature", "seed")]
                if observable == "spin_torque_onsite":
                    cases = [("feature", "seed"), ("feature", "local")]
                elif observable.endswith("_accumulation"):
                    cases.append(("feature", "local"))
                for rep in range(args.repeat if cost else 1):
                    reference = None
                    for binary, projection in cases:
                        run = args.output/f"{material}_{order}_{observable}_{binary}_{projection}_{rep}"
                        if (material, order, observable, projection, binary, rep) in done:
                            if reference is None:
                                reference = run
                            continue
                        if run.exists():
                            raise FileExistsError(f"Move the unfinished run before retrying: {run}")
                        prepare(fixture, run, args.size, order, 160, "multitype" if material == "FeCo" else "spin")
                        inp = run/"input.nml"
                        nml = f90nml.read(inp)
                        nml["lattice"]["strux_backend"] = args.strux
                        nml["control"]["linear_out"] = observable
                        nml["hamiltonian"]["pol_alpha"] = pol
                        if projection:
                            nml["control"]["response_projection"] = projection
                        nml.write(inp, force=True)
                        start = time.monotonic()
                        pid = os.fork()
                        if pid == 0:
                            os.chdir(run)
                            with open("run.log", "w") as log:
                                os.dup2(log.fileno(), 1)
                                os.dup2(log.fileno(), 2)
                                os.execve(getattr(args, binary), [str(getattr(args,binary))], env)
                        _, status, usage = os.wait4(pid, 0)
                        elapsed = time.monotonic()-start
                        if status or not (run/"response_total.out").exists():
                            raise RuntimeError(f"Run failed in {run}")
                        with (run/"clust").open() as cluster:
                            sites = int(cluster.readline().split("=")[1])
                        assert sites == args.size**3 * (2 if material == "FeCo" else 1)
                        if reference is None:
                            reference = run
                        filenames = ["response_total.out"] if projection == "local" else outputs(reference)
                        absolute = scaled = 0.
                        for filename in filenames:
                            orig = np.loadtxt(reference/filename)
                            val = np.loadtxt(run/filename)
                            error = np.max(abs(val-orig), axis=0)
                            absolute = max(absolute, np.max(error))
                            scaled = max(scaled, np.max(error/np.maximum(1, np.max(abs(orig), axis=0))))
                        writer.writerow([material, order, observable, projection, binary, rep, elapsed,
                                         usage.ru_maxrss, absolute, scaled])
                        f.flush()
