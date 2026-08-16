#!/usr/bin/env python3
"""VAL-08: damping and experimental magnetic-inertia convergence campaign.

The campaign consumes only production ``alldampings.out`` and
``allinertias.out`` records.  It does not reimplement either estimator.  The
small bcc-Fe fixture is metallic and has the supported on-site p/d SOC terms.
The existing recursion/Lehmann/Dyson damping triad remains the route oracle;
this campaign adds eta, k-mesh, energy-mesh, sign, and symmetry diagnostics.
Inertia is reported as the raw finite-difference torque kernel and is not
accepted against damping or exchange as an oracle.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FLOAT = r"[-+0-9.EeDd]+"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def patch_case(base: Path, workdir: Path, route: str, nk: int, eta: float,
               channels: int, damping: bool, inertia: bool, soc: bool) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)
    patch = {
        "calculation": {
            "gf_route": route,
            "do_damping": damping,
            "do_inertia": inertia,
        },
        "energy": {"channels_ldos": channels},
        "reciprocal": {
            "nk1": nk, "nk2": nk, "nk3": nk,
            "green_eta": eta, "green_backend": route,
            "use_symmetry_reduction": False,
        },
    }
    patched = workdir / "input.patched.nml"
    f90nml.patch(str(workdir / "input.nml"), patch, str(patched))
    text = patched.read_text()
    if "do_inertia" not in text.lower():
        text = re.sub(r"(?im)^(\s*do_damping\s*=.*)$", r"\1\n do_inertia = " + ("T" if inertia else "F"), text, count=1)
    patched.write_text(text)
    patched.replace(workdir / "input.nml")
    if not soc:
        fe = workdir / "Fe.nml"
        text = fe.read_text()
        text = re.sub(r"(?im)^\s*xi_p\s*=.*$", " xi_p = 0.0d0, 0.0d0", text, count=1)
        text = re.sub(r"(?im)^\s*xi_d\s*=.*$", " xi_d = 0.0d0, 0.0d0", text, count=1)
        fe.write_text(text)


def run_case(binary: Path, base: Path, workdir: Path, **kwargs: object) -> dict:
    patch_case(base, workdir, **kwargs)
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    runner = Path(__file__).resolve().parents[1] / "run_binary.sh"
    result = subprocess.run(
        ["/bin/bash", str(runner), str(binary)], cwd=workdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
    )
    log = (workdir / "testrun.log").read_text(errors="replace") if (workdir / "testrun.log").exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"VAL-08 run failed ({kwargs}):\n{log[-5000:]}")
    return {"parameters": kwargs, "damping": parse_damping(workdir) if kwargs["damping"] else None,
            "inertia": parse_inertia(workdir) if kwargs["inertia"] else None}


def parse_damping(workdir: Path) -> dict:
    rows = []
    for line in (workdir / "alldampings.out").read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 12 or not fields[0].lstrip("-").isdigit():
            continue
        values = [number(value) for value in fields[2:12]]
        rows.append(values)
    if not rows:
        raise RuntimeError(f"missing damping rows in {workdir}")
    matrix = [rows[0][0:3], rows[0][3:6], rows[0][6:9]]
    if not all(math.isfinite(value) for row in matrix for value in row):
        raise RuntimeError(f"non-finite damping tensor in {workdir}")
    symmetric = [[0.5 * (matrix[i][j] + matrix[j][i]) for j in range(3)] for i in range(3)]
    antisymmetric = max(abs(matrix[i][j] - matrix[j][i]) for i in range(3) for j in range(3))
    principal = [symmetric[0][0], symmetric[1][1], symmetric[2][2],
                 symmetric[0][0] * symmetric[1][1] - symmetric[0][1] ** 2,
                 symmetric[0][0] * symmetric[2][2] - symmetric[0][2] ** 2,
                 symmetric[1][1] * symmetric[2][2] - symmetric[1][2] ** 2,
                 symmetric[0][0] * (symmetric[1][1] * symmetric[2][2] - symmetric[1][2] ** 2)
                 - symmetric[0][1] * (symmetric[0][1] * symmetric[2][2] - symmetric[1][2] * symmetric[0][2])
                 + symmetric[0][2] * (symmetric[0][1] * symmetric[1][2] - symmetric[1][1] * symmetric[0][2])]
    return {"tensor": matrix, "alpha_xy": rows[0][9], "min_principal_symmetric": min(principal),
            "max_antisymmetric": antisymmetric}


def parse_inertia(workdir: Path) -> dict:
    rows = []
    for line in (workdir / "allinertias.out").read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) != 21 or not fields[0].lstrip("-").isdigit():
            continue
        values = [number(value) for value in fields[2:20]]
        rows.append(values)
    if not rows:
        raise RuntimeError(f"missing inertia rows in {workdir}")
    real = [rows[0][0:3], rows[0][3:6], rows[0][6:9]]
    imag = [rows[0][9:12], rows[0][12:15], rows[0][15:18]]
    if not all(math.isfinite(value) for row in real + imag for value in row):
        raise RuntimeError(f"non-finite inertia tensor in {workdir}")
    return {"real": real, "imag": imag,
            "max_abs": max(abs(value) for row in real + imag for value in row),
            "max_antisymmetric_real": max(abs(real[i][j] - real[j][i]) for i in range(3) for j in range(3))}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    base = root / "tests/regression/triad_bccFe_damping"
    scratch = args.scratch_root.resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    cases: dict[str, dict] = {}

    common = {"nk": 8, "eta": 0.02, "channels": 300, "soc": True}
    cases["damping_recursion"] = run_case(args.binary.resolve(), base, scratch / "damping_recursion",
                                           route="recursion", damping=True, inertia=False, **common)
    for eta in (0.01, 0.02, 0.04):
        cases[f"damping_lehmann_eta_{eta:g}"] = run_case(
            args.binary.resolve(), base, scratch / f"damping_lehmann_eta_{eta:g}",
            route="lehmann", damping=True, inertia=False, eta=eta, **{k: v for k, v in common.items() if k != "eta"})
    for nk in (4, 8, 12):
        cases[f"damping_lehmann_k_{nk}"] = run_case(
            args.binary.resolve(), base, scratch / f"damping_lehmann_k_{nk}",
            route="lehmann", damping=True, inertia=False, nk=nk, **{k: v for k, v in common.items() if k != "nk"})
    for channels in (120, 240, 360):
        cases[f"damping_lehmann_nE_{channels}"] = run_case(
            args.binary.resolve(), base, scratch / f"damping_lehmann_nE_{channels}",
            route="lehmann", damping=True, inertia=False, channels=channels,
            **{k: v for k, v in common.items() if k != "channels"})
    cases["damping_dyson_reference"] = run_case(
        args.binary.resolve(), base, scratch / "damping_dyson_reference",
        route="dyson", damping=True, inertia=False, **common)
    cases["damping_soc_zero"] = run_case(
        args.binary.resolve(), base, scratch / "damping_soc_zero",
        route="lehmann", damping=True, inertia=False, soc=False, **{k: v for k, v in common.items() if k != "soc"})

    inertia_cases = []
    for channels in (120, 240, 360):
        name = f"inertia_lehmann_nE_{channels}"
        cases[name] = run_case(args.binary.resolve(), base, scratch / name,
                               route="lehmann", damping=False, inertia=True,
                               nk=8, eta=0.02, channels=channels, soc=True)
        inertia_cases.append(cases[name]["inertia"])
    for eta in (0.01, 0.02, 0.04):
        name = f"inertia_lehmann_eta_{eta:g}"
        cases[name] = run_case(args.binary.resolve(), base, scratch / name,
                               route="lehmann", damping=False, inertia=True,
                               nk=8, eta=eta, channels=240, soc=True)
    cases["inertia_recursion_nE_240"] = run_case(
        args.binary.resolve(), base, scratch / "inertia_recursion_nE_240",
        route="recursion", damping=False, inertia=True, nk=8, eta=0.02, channels=240, soc=True)
    cases["inertia_soc_zero"] = run_case(
        args.binary.resolve(), base, scratch / "inertia_soc_zero",
        route="lehmann", damping=False, inertia=True, nk=8, eta=0.02, channels=240, soc=False)

    failures: list[str] = []
    for name, case in cases.items():
        if case["damping"] is not None:
            damping = case["damping"]
            if damping["min_principal_symmetric"] < -1.0e-8:
                failures.append(f"{name}: damping symmetric part is not positive semidefinite")
        if case["inertia"] is not None and case["inertia"]["max_abs"] > 1.0e300:
            failures.append(f"{name}: inertia overflow")
    if cases["damping_soc_zero"]["damping"]["alpha_xy"] > 1.0e-10 or max(
        abs(value) for row in cases["damping_soc_zero"]["damping"]["tensor"] for value in row) > 1.0e-10:
        failures.append("damping SOC-zero limit is nonzero")
    if cases["inertia_soc_zero"]["inertia"]["max_abs"] > 1.0e-10:
        failures.append("inertia SOC-zero limit is nonzero")

    lehmann = cases["damping_lehmann_eta_0.02"]["damping"]["alpha_xy"]
    dyson = cases["damping_dyson_reference"]["damping"]["alpha_xy"]
    if abs(lehmann - dyson) > 1.0e-5:
        failures.append(f"damping Lehmann/Dyson mismatch {abs(lehmann - dyson):.3e}")
    inertia_convergence = {
        "nE120_vs_nE240": abs(inertia_cases[0]["max_abs"] - inertia_cases[1]["max_abs"]),
        "nE240_vs_nE360": abs(inertia_cases[1]["max_abs"] - inertia_cases[2]["max_abs"]),
    }
    result = {
        "scope": "bcc-Fe metallic one-site pair; collinear nsp=2; supported on-site p/d SOC",
        "damping": {name: case["damping"] for name, case in cases.items() if case["damping"] is not None},
        "inertia": {name: case["inertia"] for name, case in cases.items() if case["inertia"] is not None},
        "inertia_convergence_absolute_max_abs": inertia_convergence,
        "interpretation": "damping has a route triad plus eta/k/energy/SOC diagnostics; inertia is raw finite-difference evidence only",
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if failures:
        print("VAL-08 FAIL")
        print("\n".join(failures))
        return 1
    print("VAL-08 PASS: damping triad diagnostics and finite inertia convergence records produced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
