#!/usr/bin/env python3
"""Validate internal/libXC XC equivalence through lean reciprocal SCF runs.

The production pair is legacy ``txc=5`` (internal PBE-LDA) and explicit
``txc=105`` (the existing ``LDA_X + LDA_C_PW`` libXC mapping).  The same
fixture, mesh, mixing, and SCF budget are used for both values.  This driver
is intentionally one compact workflow covering Si and collinear bcc Fe.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FLOAT = r"[-+0-9.EeDd]+"
REPORT_PATTERNS = {
    "total_energy": re.compile(rf"Total energy of system:\s*({FLOAT})"),
    "band_energy": re.compile(rf"Band energy of system:\s*({FLOAT})"),
    "fermi": re.compile(rf"Fermi energy:\s*({FLOAT})"),
    "occupation": re.compile(rf"Occupation at atom\s+1:\s*({FLOAT})"),
    "charge_transfer": re.compile(rf"Charge transfer at atom\s+1:\s*({FLOAT})"),
    "moment": re.compile(rf"Spin moment of atom\s+1:\s*({FLOAT})"),
    "rms": re.compile(rf"Total RMS Diff:\s*({FLOAT})"),
}
SCF_DIFF_RE = re.compile(rf"(?:Not converged! Diff=|Converged! Diff=)\s*({FLOAT})")


def as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def patch_case(case_root: Path, workdir: Path, txc: int, use_kspace: bool) -> None:
    shutil.copytree(case_root, workdir)
    for path in workdir.glob("*_out.nml"):
        path.unlink()
    patch = {
        "self": {"nstep": 24, "use_kspace": use_kspace, "conv_thr": 1.0e-8},
        "control": {"txc": txc},
        "lattice": {"strux_backend": "strux_lib"},
    }
    if use_kspace:
        patch["reciprocal"] = {
            "nk1": 4,
            "nk2": 4,
            "nk3": 4,
            "use_symmetry_reduction": False,
            "use_time_reversal": True,
            "use_shift": False,
            "dos_method": "gaussian",
            "gaussian_sigma": 0.01,
            "temperature": 300.0,
            "auto_find_fermi": True,
            "n_energy_points": 200,
        }
    input_path = workdir / "input.nml"
    patched_path = workdir / "input.patched.nml"
    f90nml.patch(str(input_path), patch, str(patched_path))
    patched_path.replace(input_path)


def first(pattern: re.Pattern[str], text: str, label: str) -> float:
    match = pattern.search(text)
    if match is None:
        raise RuntimeError(f"missing {label} in report.out")
    return as_float(match.group(1))


def parse_run(workdir: Path, system: str, txc: int, returncode: int) -> dict:
    if returncode != 0:
        log = (workdir / "testrun.log").read_text(errors="replace")
        raise RuntimeError(f"{system} txc={txc} failed:\n{log[-3000:]}")
    report = (workdir / "report.out").read_text(errors="replace")
    values = {key: first(pattern, report, key) for key, pattern in REPORT_PATTERNS.items()}
    diffs = [as_float(match.group(1)) for match in SCF_DIFF_RE.finditer(
        (workdir / "testrun.log").read_text(errors="replace")
    )]
    values["scf_diffs"] = diffs
    values["txc"] = txc
    values["system"] = system
    values["energy_components"] = {}
    for path in sorted(workdir.glob("*_out.nml")):
        parsed = f90nml.read(str(path))
        par = parsed.get("par", {})
        values["energy_components"][path.stem.removesuffix("_out")] = {
            key: float(par[key])
            for key in ("sumec", "sumev", "etot", "utot", "ekin", "rhoeps")
            if key in par
        }
    return values


def close(a: float, b: float, abs_tol: float, rel_tol: float = 0.0) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= abs_tol + rel_tol * max(abs(a), abs(b), 1.0)


def compare(internal: dict, libxc: dict) -> list[str]:
    failures: list[str] = []
    system = internal["system"]
    # Report text is six decimal places for site quantities; the tighter
    # energy/component checks use the full-precision *_out.nml records.
    report_tolerances = {
        "total_energy": 2.0e-5,
        "band_energy": 2.0e-5,
        "fermi": 2.0e-5,
        "occupation": 2.0e-5,
        "charge_transfer": 2.0e-5,
        "moment": 2.0e-5,
    }
    for key, tolerance in report_tolerances.items():
        if not close(internal[key], libxc[key], tolerance):
            failures.append(f"{system} {key}: internal={internal[key]:.12e} libXC={libxc[key]:.12e}")

    if len(internal["scf_diffs"]) != len(libxc["scf_diffs"]):
        failures.append(
            f"{system} SCF iteration count: internal={len(internal['scf_diffs'])} "
            f"libXC={len(libxc['scf_diffs'])}"
        )
    for index, (internal_diff, libxc_diff) in enumerate(zip(internal["scf_diffs"], libxc["scf_diffs"]), start=1):
        if not close(internal_diff, libxc_diff, 2.0e-4, 2.0e-2):
            failures.append(
                f"{system} SCF diff {index}: internal={internal_diff:.12e} libXC={libxc_diff:.12e}"
            )
    if not internal["scf_diffs"] or internal["scf_diffs"][-1] > 2.0e-5 or libxc["scf_diffs"][-1] > 2.0e-5:
        failures.append(f"{system} SCF did not reach the comparison residual")

    if set(internal["energy_components"]) != set(libxc["energy_components"]):
        failures.append(f"{system} output site set differs")
    for site in sorted(set(internal["energy_components"]) & set(libxc["energy_components"])):
        for key in ("sumec", "sumev", "etot", "utot", "ekin", "rhoeps"):
            a = internal["energy_components"][site].get(key)
            b = libxc["energy_components"][site].get(key)
            if a is None or b is None or not close(a, b, 2.0e-4, 2.0e-7):
                failures.append(f"{system} {site} {key}: internal={a} libXC={b}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument(
        "--kspace-only",
        action="store_true",
        help="Run only the four-run reciprocal CI subset; keep the RS check for validation campaigns.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    cases = [
        ("Si-kspace", root / "tests/scf/cases/bulk/diamondSi", True),
        ("bccFe-kspace", root / "tests/scf/cases/k_space_scf/bccFe", True),
        # One shared-XC-state check through the real-space recursion route.
        ("bccFe-rs", root / "tests/scf/cases/bulk/bccFe", False),
    ]
    if args.kspace_only:
        cases = [case for case in cases if case[2]]
    runner = root / "tests/run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[int, dict]] = {}
    failures: list[str] = []
    for system, case_root, use_kspace in cases:
        results[system] = {}
        for txc in (5, 105):
            workdir = args.scratch_root / f"{system}_txc{txc}"
            if workdir.exists():
                shutil.rmtree(workdir)
            patch_case(case_root, workdir, txc, use_kspace)
            completed = subprocess.run(
                ["bash", str(runner), str(args.binary.resolve())],
                cwd=workdir,
                text=True,
                capture_output=True,
                check=False,
            )
            results[system][txc] = parse_run(workdir, system, txc, completed.returncode)
            print(json.dumps(results[system][txc], sort_keys=True))
        failures.extend(compare(results[system][5], results[system][105]))

    if failures:
        print("VAL-03 FAIL")
        print("\n".join(failures))
        return 1
    print("VAL-03 PASS: txc=5 internal and txc=105 libXC agree for Si and bccFe")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
