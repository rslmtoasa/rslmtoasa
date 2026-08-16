#!/usr/bin/env python3
"""Run the compact VAL-02 reciprocal-SCF convergence campaign.

This is a campaign driver, not a CTest/quick test.  It runs the production
binary in isolated directories, parses production observables, and emits JSON
that can be summarized in the VAL-02 report.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

import f90nml


FLOAT = r"[-+0-9.EeDd]+"
CANONICAL_RE = re.compile(
    rf"Canonical k-space occupations: EF=\s*({FLOAT}).*?N=\s*({FLOAT}),"
    rf"\s*dN=\s*({FLOAT}),.*?weight_sum\(raw\)=\s*({FLOAT}),\s*EBAND=\s*({FLOAT})"
)
DOS_RE = re.compile(rf"DOS integrates to\s*({FLOAT})")
DOS_TETRA_RE = re.compile(rf"DOS integral\s*=\s*({FLOAT})")
DIFF_RE = re.compile(rf"(?:Not converged! Diff=|Converged! Diff=)\s*({FLOAT})")
TOTAL_RMS_RE = re.compile(rf"Total RMS Diff:\s*({FLOAT})")
ENERGY_RE = re.compile(rf"Total energy of system:\s*({FLOAT})")
BAND_RE = re.compile(rf"Band energy of system:\s*({FLOAT})")
MOM_RE = re.compile(rf"Spin moment of atom\s+1:\s*({FLOAT})")
MOM_PROJ_RE = re.compile(
    rf"Spin moment projections of atom\s+1:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})"
)
OCC_RE = re.compile(rf"Occupation at atom\s+1:\s*({FLOAT})")
CHARGE_RE = re.compile(rf"Charge transfer at atom\s+1:\s*({FLOAT})")
HERM_FATAL_RE = re.compile(r"non-Hermitian|non-Hermitian before eigensolution", re.I)


def as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def patch_case(case_root: Path, workdir: Path, patch: dict) -> None:
    shutil.copytree(case_root, workdir)
    for path in workdir.glob("*_out.nml"):
        path.unlink()
    input_path = workdir / "input.nml"
    patched_path = workdir / "input.patched.nml"
    input_patch = {key: value for key, value in patch.items() if not key.startswith("__")}
    f90nml.patch(str(input_path), input_patch, str(patched_path))
    patched_path.replace(input_path)
    if "__atom_mom" in patch:
        atom_path = workdir / "Fe.nml"
        atom_patched_path = workdir / "Fe.patched.nml"
        f90nml.patch(str(atom_path), {"par": {"mom": patch["__atom_mom"]}}, str(atom_patched_path))
        atom_patched_path.replace(atom_path)


def parse_first(pattern: re.Pattern[str], text: str) -> float | None:
    match = pattern.search(text)
    return as_float(match.group(1)) if match else None


def parse_vector(pattern: re.Pattern[str], text: str) -> list[float] | None:
    match = pattern.search(text)
    return [as_float(match.group(i)) for i in range(1, 4)] if match else None


def integrate_dos(path: Path) -> float | None:
    points: list[tuple[float, float]] = []
    if not path.exists():
        return None
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 2 or fields[0].startswith("#"):
            continue
        try:
            points.append((float(fields[0]), float(fields[1])))
        except ValueError:
            continue
    if len(points) < 2:
        return None
    return sum(
        0.5 * (right[1] + left[1]) * (right[0] - left[0])
        for left, right in zip(points, points[1:])
    )


def parse_run(workdir: Path, name: str, patch: dict, returncode: int) -> dict:
    log = (workdir / "testrun.log").read_text(errors="replace") if (workdir / "testrun.log").exists() else ""
    report = (workdir / "report.out").read_text(errors="replace") if (workdir / "report.out").exists() else ""
    canonical_matches = list(CANONICAL_RE.finditer(log))
    canonical = canonical_matches[-1] if canonical_matches else None
    method = patch["reciprocal"]["dos_method"]
    dos_pattern = DOS_RE if method == "gaussian" else DOS_TETRA_RE
    dos_matches = list(dos_pattern.finditer(log))
    dos_match = dos_matches[-1] if dos_matches else None
    diffs = [as_float(match.group(1)) for match in DIFF_RE.finditer(log)]
    result = {
        "name": name,
        "returncode": returncode,
        "mesh": patch["reciprocal"]["nk1"],
        "dos_method": method,
        "gaussian_sigma": patch["reciprocal"].get("gaussian_sigma"),
        "nstep": patch["self"]["nstep"],
        "nsp": patch["control"]["nsp"],
        "canonical": None,
        "dos_integral_log": as_float(dos_match.group(1)) if dos_match else None,
        "dos_integral_file": integrate_dos(workdir / "totaldos.out"),
        "scf_diffs": diffs,
        "total_rms_diff": parse_first(TOTAL_RMS_RE, report),
        "total_energy": parse_first(ENERGY_RE, report),
        "report_band_energy": parse_first(BAND_RE, report),
        "site_moment": parse_first(MOM_RE, report),
        "site_moment_vector": parse_vector(MOM_PROJ_RE, report),
        "site_occupation": parse_first(OCC_RE, report),
        "site_charge_transfer": parse_first(CHARGE_RE, report),
        "hermiticity_violation": bool(HERM_FATAL_RE.search(log)),
        "fatal_tail": log[-1200:] if returncode else None,
    }
    output_name = "Si1_out.nml" if patch["control"]["lmax"] == 1 else "Fe_out.nml"
    output_path = workdir / output_name
    if output_path.exists():
        output_nml = f90nml.read(str(output_path))
        par = output_nml.get("par", {})
        result["energy_components"] = {
            key: float(par[key])
            for key in ("sumec", "sumev", "etot", "utot", "ekin", "rhoeps")
            if key in par
        }
    if canonical:
        result["canonical"] = {
            "ef": as_float(canonical.group(1)),
            "electron_count": as_float(canonical.group(2)),
            "electron_residual": as_float(canonical.group(3)),
            "weight_sum": as_float(canonical.group(4)),
            "band_energy": as_float(canonical.group(5)),
        }
    return result


def base_patch(system: str, nsp: int, mesh: int, method: str, sigma: float, nstep: int) -> dict:
    is_si = system == "si"
    return {
        "lattice": {"strux_backend": "strux_lib"},
        "self": {"nstep": nstep, "use_kspace": True, "conv_thr": 1.0e-8},
        "energy": {
            "channels_ldos": 2001 if method == "tetrahedron" else 200,
            "energy_min": -1.5 if is_si else -2.0,
            "energy_max": 1.0 if is_si else 2.0,
        },
        "control": {
            "nsp": nsp,
            "lmax": 1 if is_si else 2,
            "recur": "block",
            "lld": 12,
        },
        "reciprocal": {
            "nk1": mesh,
            "nk2": mesh,
            "nk3": mesh,
            "use_symmetry_reduction": False,
            "use_time_reversal": True,
            "use_shift": False,
            "n_energy_points": 2001 if method == "tetrahedron" else 200,
            "dos_method": method,
            "gaussian_sigma": sigma,
            "temperature": 300.0,
            "auto_find_fermi": True,
        },
    }


def campaign_cases() -> list[tuple[str, str, dict]]:
    cases: list[tuple[str, str, dict]] = []
    for system in ("si", "fe"):
        nsp = 1 if system == "si" else 2
        for mesh in (2, 4, 6, 8, 10, 12, 16):
            patch = base_patch(system, nsp, mesh, "gaussian", 0.01, 12)
            cases.append((f"{system}_mesh{mesh}_gaussian_n12", system, patch))
        for method, sigma in (("gaussian", 0.01), ("gaussian", 0.02), ("tetrahedron", 0.01)):
            patch = base_patch(system, nsp, 4, method, sigma, 12)
            cases.append((f"{system}_dos_{method}_sigma{sigma:g}_mesh4_n12", system, patch))
        for nstep in (1, 3, 6, 12, 24):
            patch = base_patch(system, nsp, 4, "gaussian", 0.01, nstep)
            cases.append((f"{system}_scf_n{nstep}_mesh4_gaussian", system, patch))

    for nsp in (3, 4):
        patch = base_patch("fe", nsp, 2, "gaussian", 0.01, 12)
        cases.append((f"fe_nsp{nsp}_mesh2_gaussian_n12", "fe", patch))
        patch = base_patch("fe", nsp, 2, "gaussian", 0.01, 12)
        patch["hamiltonian"] = {"hoh": True}
        patch["__atom_mom"] = [1.0, 0.0, 0.0]
        patch["energy"] = {"channels_ldos": 200, "energy_min": -4.0, "energy_max": 4.0}
        cases.append((f"fe_nsp{nsp}_hoh_axisx_mesh2_gaussian_n12", "fe", patch))

    patch = base_patch("fe", 2, 8, "tetrahedron", 0.01, 12)
    cases.append(("fe_tetra_mesh8_n12", "fe", patch))
    return cases


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extended-only", action="store_true")
    parser.add_argument("--high-mesh-only", action="store_true")
    parser.add_argument("--dense-fe-only", action="store_true")
    parser.add_argument("--tetra-high-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    case_roots = {
        "si": root / "tests/scf/cases/bulk/diamondSi",
        "fe": root / "tests/scf/cases/k_space_scf/bccFe",
    }
    runner = root / "tests/run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    results = []
    cases = campaign_cases()
    if args.extended_only:
        cases = [case for case in cases if "hoh_axisx" in case[0]]
    if args.high_mesh_only:
        cases = [case for case in cases if "mesh8_" in case[0] or "mesh10_" in case[0]]
    if args.dense_fe_only:
        cases = [case for case in cases if case[0].startswith("fe_mesh12_") or case[0].startswith("fe_mesh16_")]
    if args.tetra_high_only:
        cases = [case for case in cases if case[0] == "fe_tetra_mesh8_n12"]
    for name, system, patch in cases:
        workdir = args.scratch_root / name
        if workdir.exists():
            shutil.rmtree(workdir)
        patch_case(case_roots[system], workdir, patch)
        completed = subprocess.run(
            ["bash", str(runner), str(args.binary.resolve())],
            cwd=workdir,
            text=True,
            capture_output=True,
            check=False,
        )
        results.append(parse_run(workdir, name, patch, completed.returncode))
        print(f"{name}: rc={completed.returncode}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    failures = [item for item in results if item["returncode"] != 0]
    print(f"completed={len(results)} failures={len(failures)} output={args.output}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
