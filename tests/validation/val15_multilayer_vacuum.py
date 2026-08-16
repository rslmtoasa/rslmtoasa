#!/usr/bin/env python3
"""VAL-15: validate the existing multilayer A|vacuum production route.

The fixture is deliberately assembled from the fcc-Cu(001) buildsurf example.
The short buildsurf probe is used to establish the symbolic layer/type order;
the converged parameter files already shipped with that same example are then
used for the interface run.  This keeps the validation from silently becoming
another hand-authored layer mapping.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path

FLOAT = r"[-+0-9.EeDd]+"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_input(source: Path, replacements: dict[str, str]) -> None:
    text = source.read_text()
    for old, new in replacements.items():
        require(old in text, f"input template no longer contains {old!r}")
        text = text.replace(old, new, 1)
    source.write_text(text)


def run_binary(binary: Path, root: Path, workdir: Path) -> str:
    result = subprocess.run(
        ["/bin/bash", str(root / "tests/run_binary.sh"), str(binary)],
        cwd=workdir, capture_output=True, text=True,
    )
    log_path = workdir / "testrun.log"
    log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout
    if result.returncode != 0:
        raise RuntimeError(f"VAL-15 run failed in {workdir}:\n{log[-8000:]}")
    return log


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("VAL-15 FAIL: " + message)


def buildsurf_mapping(binary: Path, root: Path, scratch: Path) -> dict:
    source = root / "example/surface/fccCu001"
    workdir = scratch / "buildsurf_nlay3"
    shutil.copytree(source, workdir)
    replace_input(workdir / "input.nml", {
        "nlay = 6": "nlay = 3",
        "ntype = 7": "ntype = 4",
        "ct(:) = 7*4.0d0": "ct(:) = 4*4.0d0",
        "label(7) = 'Cu-S-4'\n": "",
        "label(6) = 'Cu-S-3'\n": "",
        "label(5) = 'Cu-S-2'\n": "",
        "nstep = 500": "nstep = 1",
    })
    run_binary(binary, root, workdir)
    surfclu = (workdir / "surfclu.out").read_text(errors="replace")

    require("Maxtype:" in surfclu and re.search(r"Maxtype:\s+4", surfclu),
            "buildsurf did not construct the expected four symbolic types")
    require(re.search(r"Unique atoms type per layer:\s+1\s+1\s+1", surfclu),
            "buildsurf did not produce one active symbolic type per layer")
    rows = re.findall(
        r"^\s*-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+(\d+)\s+1\s+",
        surfclu, re.MULTILINE,
    )
    require(rows[:4] == ["1", "2", "3", "4"],
            f"unexpected buildsurf layer/type rows: {rows[:4]}")
    return {
        "symbolic_types": ["Cu", "ES", "Cu-S", "Cu-S-1"],
        "active_types": ["ES", "Cu-S", "Cu-S-1"],
        "source": str(source),
    }


def make_interface_case(root: Path, scratch: Path, buffer_width: int, nstep: int) -> Path:
    workdir = scratch / f"interface_buffer_{buffer_width}"
    workdir.mkdir(parents=True, exist_ok=True)
    source = root / "example/surface/fccCu001"
    vacuum = root / "example/interface/fccCu111_Avac/Vac.nml"
    shutil.copyfile(source / "Cu_out.nml", workdir / "A.nml")
    shutil.copyfile(source / "ES_out.nml", workdir / "Ac1.nml")
    shutil.copyfile(source / "Cu-S_out.nml", workdir / "Ac2.nml")
    shutil.copyfile(source / "Cu-S-1_out.nml", workdir / "Ac3.nml")
    shutil.copyfile(vacuum, workdir / "Vac.nml")
    shutil.copyfile(root / "example/interface/fccCu111_Avac/input.nml", workdir / "input.nml")
    replace_input(workdir / "input.nml", {
        "pre_processing = 'buildinterface'": "pre_processing = 'buildinterface'",
        "nlay = 1                  ! active layers": "nlay = 3                  ! active layers",
        "nlay_a = 1                ! frozen region A width": f"nlay_a = {buffer_width}                ! frozen region A width",
        "nlay_b = 1                ! vacuum width": f"nlay_b = {buffer_width}                ! vacuum width",
        "surftype = '1 1 1'": "surftype = '0 0 1'",
        "ntype = 3                 ! A + vacuum + one type per active layer": "ntype = 5                 ! A + vacuum + one type per active layer",
        "ct(:) = 3*4.0d0": "ct(:) = 5*4.0d0",
        "label(3) = 'Ac1'          ! active layer 1": "label(3) = 'Ac1'          ! buildsurf type 2: ES",
        "label(3) = 'Ac1'          ! buildsurf type 2: ES": "label(3) = 'Ac1'          ! buildsurf type 2: ES\nlabel(4) = 'Ac2'          ! buildsurf type 3: Cu-S\nlabel(5) = 'Ac3'          ! buildsurf type 4: Cu-S-1",
        "nstep = 30": f"nstep = {nstep}",
        "recur = 'chebyshev'": "recur = 'block'",
        "lld = 150": "lld = 21",
        "&charge\n": "&charge\nvmix = 0.2\n",
        "mixtype = 'broyden'": "mixtype = 'linear'",
    })
    return workdir


def parse_case(workdir: Path, log: str) -> dict:
    deck = (workdir / "surfclu.out").read_text(errors="replace")
    require(re.search(r"nlay_a, nlay, nlay_b:\s+\d+\s+3\s+\d+", deck),
            "interface deck did not retain three active layers")
    require(re.search(r"ntype, nbulk, nrec:\s+5\s+2\s+3", deck),
            "interface deck has the wrong frozen/active type counts")
    require(re.search(r"Sites per active layer:\s+1\s+1\s+1", deck),
            "interface deck has the wrong per-layer representative count")
    require(re.search(r"Active type -> charge-transfer reference type:\s+1\s+1\s+1", deck),
            "A|vacuum active types are not all referenced to region A")

    charge = re.findall(
        r"interfacepot: Q=\s*(%s) P=\s*(%s) residual=\s*(%s)" % (FLOAT, FLOAT, FLOAT), log
    )
    steps = re.findall(r"interfacepot: deep-probe potentials.*?step=\s*(%s)" % FLOAT, log)
    compensation = re.findall(
        r"compensation rows low/high=\s*(\d+)/(\d+).*?weights=\s*(%s)/\s*(%s)" % (FLOAT, FLOAT), log
    )
    vmad = re.findall(r"interfacepot: active type\s+(\d+) vmad=\s*(%s)" % FLOAT, log)
    require(charge and steps and compensation and vmad,
            "production diagnostics did not expose charge, compensation, step, and profile")

    q, p, residual = [number(value) for value in charge[-1]]
    step = number(steps[-1])
    low, high, wlow, whigh = compensation[-1]
    profile = [number(value) for _, value in vmad[-3:]]
    require(all(math.isfinite(value) for value in [q, p, residual, step] + profile),
            "non-finite electrostatic/profile quantity")
    require(abs(q) < 1.0e-8, f"compensated deviation charge is not neutral: Q={q}")
    require(number(whigh) < 1.0e-12 and abs(number(wlow) - 1.0) < 1.0e-6,
            f"vacuum received compensation: weights={wlow}/{whigh}")

    # The production Green/DOS path writes a local active-type DOS.  Vacuum is
    # frozen and is not a separately exposed local-GF observable on calctype=L.
    dos_path = workdir / "ES_dos.out"
    require(dos_path.exists(), "active local DOS was not emitted")
    dos = []
    for line in dos_path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            try:
                dos.append((float(fields[0]), float(fields[1])))
            except ValueError:
                pass
    require(dos and all(math.isfinite(x) and math.isfinite(y) for x, y in dos),
            "active local DOS contains non-finite values")
    require(min(y for _, y in dos) > -1.0e-5,
            "active local DOS violates the non-negative spectral-density check")
    return {
        "Q": q,
        "P": p,
        "residual": residual,
        "step_Ry": step,
        "active_vmad_Ry": profile,
        "compensation_rows": [int(low), int(high)],
        "compensation_weights": [number(wlow), number(whigh)],
        "active_dos_points": len(dos),
        "vacuum_local_gf": "not separately exposed; active local DOS finite/non-negative checked",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    scratch = args.scratch_root.resolve()
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)

    mapping = buildsurf_mapping(args.binary.resolve(), root, scratch)
    cases = {}
    for width in (1, 4):
        workdir = make_interface_case(root, scratch, width, nstep=6)
        log = run_binary(args.binary.resolve(), root, workdir)
        cases[str(width)] = parse_case(workdir, log)
        cases[str(width)]["buffer_width"] = width

    step_delta = abs(cases["1"]["step_Ry"] - cases["4"]["step_Ry"])
    require(step_delta < 0.25,
            f"buffer sensitivity is too large for this short damped campaign: {step_delta:.3e} Ry")
    print({"mapping": mapping, "cases": cases, "buffer_step_delta_Ry": step_delta})
    print("VAL-15 PASS: buildsurf-derived multilayer A|vacuum electrostatics checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
