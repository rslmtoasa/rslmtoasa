#!/usr/bin/env python3
"""Run the bounded ordinary-q=0 magnetic SCF closure campaign.

The campaign deliberately keeps the production input decks intact.  Each run
is made in a private scratch directory, uses ``fresh_start=.true.``, and
records only scalar observables plus the production diagnostics.  The
``magnetic_seed_*`` controls are a temporary B_fsm perturbation; they are not
fixed-spin constraints.

Example::

    python3 tests/magnetic_scf/run_magnetic_scf_regression.py \
      --binary build-xcr/bin/rslmto.x \
      --output-dir tests/magnetic_scf/results
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FCC_ALAT = 3.683
FCC_WS_ANGSTROM = 1.439302849868756
ANGSTROM_TO_BOHR = 1.8897261246257702
FCC_WS_BOHR = FCC_WS_ANGSTROM * ANGSTROM_TO_BOHR
FCC_SEEDS = (0.0, 0.1, 0.5, 1.0, 2.0, 3.0)
FUNCTIONALS = (1, 101, 5, 105)
FUNCTIONAL_NAMES = {
    1: "legacy Barth-Hedin",
    101: "libXC Slater + von Barth-Hedin",
    5: "legacy PW/PBE-LDA branch",
    105: "libXC Slater + Perdew-Wang",
}


def patch_namelist(path: Path, values: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    f90nml.patch(str(path), values, str(temporary))
    temporary.replace(path)


def prepare_case(base: Path, workdir: Path, kind: str, txc: int, seed: float, steps: int,
                 fresh_start: bool = True) -> float:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)
    for output in workdir.glob("*_out.nml"):
        output.unlink()

    field = 0.0 if seed == 0.0 else 0.005 * seed
    lattice = {}
    if kind == "fcc_fe":
        lattice = {
            "crystal_sym": "fcc",
            "alat": FCC_ALAT,
            "wav": FCC_WS_ANGSTROM,
            "strux_backend": "strux_lib",
        }
        patch_namelist(workdir / "Fe.nml", {"par": {"ws_r": FCC_WS_BOHR}})
    elif kind == "bcc_fe":
        lattice = {"strux_backend": "strux_lib"}
    elif kind == "fcc_cu":
        lattice = {"strux_backend": "strux_lib"}
    else:
        raise ValueError(kind)

    control = {
        "txc": txc,
        "nsp": 2,
        "recur": "block",
        "lld": 21,
        "fresh_start": fresh_start,
    }
    self_values = {
        "nstep": steps,
        "conv_thr": 1.0e-6,
        "magnetic_scf_diagnostics": True,
        "magnetic_seed_enable": seed != 0.0,
        "magnetic_seed_steps": 4 if seed != 0.0 else 0,
        "magnetic_seed_field": field,
    }
    if kind == "fcc_cu":
        control.update({"nsp": 2, "recur": "chebyshev", "lld": 80})
    patch_namelist(workdir / "input.nml", {"lattice": lattice, "control": control, "self": self_values})
    return field


def run_binary(binary: Path, workdir: Path) -> None:
    wrapper = Path(__file__).resolve().parents[1] / "run_binary.sh"
    completed = subprocess.run(
        ["/bin/bash", str(wrapper), str(binary)],
        cwd=workdir,
        env={**os.environ, "RSLMTO_OMP_THREADS_SERIAL": "1"},
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"{binary} failed in {workdir}; see {workdir / 'testrun.log'}")


def diagnostic_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append(line.split())
    return rows


def residual_summary(path: Path) -> dict[str, float]:
    rows = diagnostic_rows(path)
    if not rows:
        raise RuntimeError(f"no atomic residual audit in {path.parent}")
    return {
        "residual_unweighted_first": float(rows[0][3]),
        "residual_integrated_first": float(rows[0][4]),
        "residual_unweighted_last": float(rows[-1][3]),
        "residual_integrated_last": float(rows[-1][4]),
    }


def summarize_run(workdir: Path, kind: str, txc: int, seed: float, field: float, steps: int,
                 fresh_start: bool) -> dict:
    trace = diagnostic_rows(workdir / "magnetic_scf_diagnostics.dat")
    if not trace:
        raise RuntimeError(f"no magnetic trace in {workdir}")
    last_iteration = max(int(row[0]) for row in trace)
    last = next(row for row in trace if int(row[0]) == last_iteration)
    d_channel = next(row for row in trace if int(row[0]) == last_iteration and int(row[2]) == 2)
    log = (workdir / "testrun.log").read_text(errors="replace")
    output = f90nml.read(str(workdir / ("Fe_out.nml" if kind != "fcc_cu" else "Cu_out.nml")))
    par = output.get("par", {})
    provenance = re.findall(r"SCF input provenance: ([^\n]+)", log)
    if not provenance:
        raise RuntimeError(f"missing source provenance in {workdir / 'testrun.log'}")
    moment = float(last[17])
    release_iteration = min(4, last_iteration) if field != 0.0 else 1
    release_row = next(row for row in trace if int(row[0]) == release_iteration)
    unconstrained_iteration = min(release_iteration + 1, last_iteration)
    unconstrained_row = next(row for row in trace if int(row[0]) == unconstrained_iteration)
    result = {
        "system": {"fcc_fe": "fcc Fe", "bcc_fe": "bcc Fe", "fcc_cu": "fcc Cu"}[kind],
        "xc": f"TXC={txc}",
        "functional": FUNCTIONAL_NAMES.get(txc, "unknown"),
        "backend": "libXC" if txc >= 100 else "legacy",
        "seed_mub_nominal": seed,
        "seed_bfsm_ry": field,
        "moment_at_seed_release_mub": float(release_row[17]),
        "first_unconstrained_M_mub": float(unconstrained_row[17]),
        "final_M_mub": moment,
        "delta_Cd": float(d_channel[13]),
        "energy_ry": float(par.get("etot", "nan")),
        "band_energy_ry": float(par.get("sumev", "nan")),
        "iterations": last_iteration,
        "status": "CONVERGED" if "Converged!" in log else "BOUND_NOT_CONVERGED",
        "provenance": provenance[0].strip(),
        "fresh_start_requested": fresh_start,
        "diagnostic_rows": len(trace),
    }
    result.update(residual_summary(workdir / "atomic_scf_residuals.dat"))
    return result


def run_one(base: Path, scratch: Path, binary: Path, kind: str, txc: int, seed: float, steps: int,
            output_dir: Path, fresh_start: bool = True, keep: bool = False) -> dict:
    name = f"{kind}_txc{txc}_seed{seed:g}_{'fresh' if fresh_start else 'restart'}"
    workdir = scratch / name
    field = prepare_case(base, workdir, kind, txc, seed, steps, fresh_start=fresh_start)
    run_binary(binary, workdir)
    result = summarize_run(workdir, kind, txc, seed, field, steps, fresh_start)
    result["case"] = name
    if keep:
        trace_dir = output_dir / "traces"
        trace_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(workdir / "magnetic_scf_diagnostics.dat", trace_dir / f"{name}.dat")
        shutil.copy2(workdir / "atomic_scf_residuals.dat", trace_dir / f"{name}.residuals.dat")
        shutil.copy2(workdir / "testrun.log", trace_dir / f"{name}.log")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--case-root", type=Path, default=Path("tests/scf/cases"))
    parser.add_argument("--scratch-root", type=Path, default=Path("/tmp/rslmto_xcr04_magnetic"))
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument("--keep-traces", action="store_true")
    args = parser.parse_args()

    args.binary = args.binary.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.scratch_root.exists():
        shutil.rmtree(args.scratch_root)
    args.scratch_root.mkdir(parents=True)
    rows: list[dict] = []
    bases = {
        "fcc_fe": args.case_root / "bulk" / "bccFe",
        "bcc_fe": args.case_root / "bulk" / "bccFe",
        "fcc_cu": args.case_root / "bulk" / "fccCu",
    }

    for txc in FUNCTIONALS:
        for seed in FCC_SEEDS:
            rows.append(run_one(bases["fcc_fe"], args.scratch_root, args.binary, "fcc_fe", txc, seed,
                                args.max_steps, args.output_dir, keep=args.keep_traces))
    for seed in (0.0, 2.0):
        rows.append(run_one(bases["bcc_fe"], args.scratch_root, args.binary, "bcc_fe", 1, seed,
                            args.max_steps, args.output_dir, keep=args.keep_traces))
    rows.append(run_one(bases["fcc_cu"], args.scratch_root, args.binary, "fcc_cu", 1, 2.0,
                        args.max_steps, args.output_dir, keep=args.keep_traces))

    # Reuse one completed fresh bcc run's generated Fe_out.nml to prove the
    # legacy restart route is explicit.  This is a provenance check, not a
    # second physics result.
    restart_base = args.scratch_root / "bcc_restart_source"
    shutil.copytree(bases["bcc_fe"], restart_base)
    prepare_case(bases["bcc_fe"], restart_base, "bcc_fe", 1, 0.0, 1, fresh_start=True)
    run_binary(args.binary, restart_base)
    patch_namelist(restart_base / "input.nml", {"self": {"nstep": 1}, "control": {"fresh_start": False}})
    run_binary(args.binary, restart_base)
    restart_result = summarize_run(restart_base, "bcc_fe", 1, 0.0, 0.0, 1, False)
    restart_result["case"] = "bcc_fe_restart_provenance"
    restart_result["status"] = "PROVENANCE_CHECK"
    rows.append(restart_result)

    fieldnames = [
        "case", "system", "xc", "functional", "backend", "seed_mub_nominal", "seed_bfsm_ry",
        "moment_at_seed_release_mub", "first_unconstrained_M_mub", "final_M_mub", "delta_Cd", "energy_ry",
        "band_energy_ry", "iterations", "status",
        "provenance", "fresh_start_requested", "diagnostic_rows", "residual_unweighted_first",
        "residual_integrated_first", "residual_unweighted_last", "residual_integrated_last",
    ]
    with (args.output_dir / "magnetic_scf_regression.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "campaign": "XCR-04 ordinary q=0 collinear magnetic SCF closure",
        "binary": str(args.binary.resolve()),
        "max_steps": args.max_steps,
        "fcc_fe": {"alat_angstrom": FCC_ALAT, "ws_radius_angstrom": FCC_WS_ANGSTROM,
                    "nominal_seed_mub": list(FCC_SEEDS), "seed_field_rule": "0.005 Ry per nominal mu_B label"},
        "functionals": FUNCTIONAL_NAMES,
        "temporary_seed": "B_fsm applied for first four outer iterations, then released; no constraint namelist is enabled",
        "rows": len(rows),
    }
    (args.output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"rows": len(rows), "output": str(args.output_dir / "magnetic_scf_regression.csv")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
