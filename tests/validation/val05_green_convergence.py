#!/usr/bin/env python3
"""VAL-05: direct reciprocal Green-function convergence campaign.

The production ``kspace_green`` post-processing route is run on the fixed
bcc-Fe potential used by the existing route triads.  The Fortran driver writes
selected complex matrix elements of G_ii and G_ij, the full-block RS/Lehmann
differences, the trace DOS, and the Dyson(Sigma=0) invariant.  This script
only parses and compares those production records; it does not implement a
second Green-function evaluator.
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
SUMMARY_RE = {
    "dos_max": re.compile(rf"# C1  max\|dos_lehmann - dos_rs\|\s*=\s*({FLOAT})"),
    "dos_rms": re.compile(rf"#     rms dos_lehmann - dos_rs\s*=\s*({FLOAT})"),
    "block_max": re.compile(rf"# C2  max\|G_ii\^lehmann - G_ii\^rs\|\s*=\s*({FLOAT})"),
    "dyson_max": re.compile(rf"# B2\.4 max\|gij\^dyson - gij\^lehmann\|\s*=\s*({FLOAT})"),
    "weight_rs": re.compile(rf"# spectral weight\s+RS =\s*({FLOAT})"),
    "weight_lehmann": re.compile(rf"# spectral weight\s+RS =\s*{FLOAT}\s+Lehmann =\s*({FLOAT})"),
}
PAIR_BLOCK_RE = re.compile(rf"# direct G pair\s+(\d+)\s+\(\s*(\d+)\s+(\d+)\)\s+max\|G_RS-G_Lehmann\|\s*=\s*({FLOAT})")

PROBE_ENERGIES = (-0.8, -0.4, 0.0, 0.4, 0.8)
PAIR_KEYS = ("1_1", "1_335")


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def close(a: float, b: float, atol: float, rtol: float = 0.0) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= atol + rtol * max(abs(a), abs(b), 1.0)


def patch_case(base: Path, workdir: Path, nk: int, eta: float, energy_min: float = -1.2, energy_max: float = 1.2) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)
    patch = {
        "calculation": {"post_processing": "kspace_green", "gf_route": "recursion"},
        "self": {"nstep": 0},
        "energy": {
            "energy_min": energy_min,
            "energy_max": energy_max,
            "channels_ldos": 120,
        },
        "reciprocal": {
            "nk1": nk,
            "nk2": nk,
            "nk3": nk,
            "use_symmetry_reduction": False,
            "use_time_reversal": False,
            "green_backend": "lehmann",
            "green_eta": eta,
        },
        "lattice": {"njij": 2},
    }
    patched = workdir / "input.patched.nml"
    f90nml.patch(str(workdir / "input.nml"), patch, str(patched))
    text = patched.read_text()
    text = re.sub(r"(?im)^\s*ijpair\(1,\s*:\)\s*=.*$", "ijpair(1, :) = 1, 1", text)
    text = re.sub(r"(?im)^\s*ijpair\(2,\s*:\)\s*=.*$", "ijpair(2, :) = 1, 335", text)
    patched.write_text(text)
    patched.replace(workdir / "input.nml")


def run_case(binary: Path, base: Path, workdir: Path, nk: int, eta: float, energy_min: float = -1.2, energy_max: float = 1.2) -> dict:
    patch_case(base, workdir, nk, eta, energy_min, energy_max)
    runner = Path(__file__).resolve().parents[1] / "run_binary.sh"
    env = os.environ.copy()
    env["RSLMTO_OMP_THREADS_SERIAL"] = "1"
    result = subprocess.run(
        ["/bin/bash", str(runner), str(binary)], cwd=workdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
    )
    if result.returncode != 0:
        log = workdir / "testrun.log"
        detail = log.read_text(errors="replace")[-5000:] if log.exists() else result.stdout[-5000:]
        raise RuntimeError(f"VAL-05 run failed (nk={nk}, eta={eta}, window={energy_min}:{energy_max})\n{detail}")

    report = (workdir / "kspace_green_c1.dat").read_text(errors="replace")
    summary: dict[str, float] = {}
    for key, pattern in SUMMARY_RE.items():
        match = pattern.search(report)
        if match is None:
            raise RuntimeError(f"missing {key} in {workdir}/kspace_green_c1.dat")
        summary[key] = number(match.group(1))
    pair_block_max: dict[str, float] = {}
    for match in PAIR_BLOCK_RE.finditer(report):
        pair_block_max[f"{int(match.group(2))}_{int(match.group(3))}"] = number(match.group(4))
    for pair in PAIR_KEYS:
        if pair not in pair_block_max:
            raise RuntimeError(f"missing full-block direct-G metric for {pair}")

    rows: dict[str, list[dict[str, float]]] = {}
    for line in (workdir / "kspace_green_gf.dat").read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) != 13 or fields[0].startswith("#"):
            continue
        i, j = int(fields[0]), int(fields[1])
        values = [number(value) for value in fields[2:]]
        key = f"{i}_{j}"
        rows.setdefault(key, []).append({
            "energy": values[0],
            "rs": complex(values[1], values[2]),
            "lehmann": complex(values[3], values[4]),
            "dyson": complex(values[5], values[6]),
            "rs_down": complex(values[7], values[8]),
            "lehmann_down": complex(values[9], values[10]),
        })
    missing = [key for key in PAIR_KEYS if key not in rows]
    if missing:
        raise RuntimeError(f"missing selected direct-G pairs {missing} in {workdir}/kspace_green_gf.dat")
    for key in PAIR_KEYS:
        rows[key].sort(key=lambda row: row["energy"])
        if len(rows[key]) < 50:
            raise RuntimeError(f"too few direct-G rows for {key}: {len(rows[key])}")
    return {"nk": nk, "eta": eta, "energy_window": [energy_min, energy_max], "summary": summary,
            "pair_block_max": pair_block_max, "rows": rows}


def interpolate(rows: list[dict[str, float]], field: str, energy: float) -> complex:
    if energy < rows[0]["energy"] or energy > rows[-1]["energy"]:
        raise ValueError(f"energy {energy} outside direct-G grid")
    for left, right in zip(rows, rows[1:]):
        if left["energy"] <= energy <= right["energy"]:
            fraction = (energy - left["energy"]) / (right["energy"] - left["energy"])
            return left[field] + fraction * (right[field] - left[field])
    return rows[-1][field]


def direct_probe(case: dict, field: str, pair: str, energy: float) -> complex:
    return interpolate(case["rows"][pair], field, energy)


def max_probe_delta(left: dict, right: dict, field: str, pairs: tuple[str, ...] = PAIR_KEYS) -> float:
    return max(
        abs(direct_probe(left, field, pair, energy) - direct_probe(right, field, pair, energy))
        for pair in pairs for energy in PROBE_ENERGIES
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = Path(__file__).resolve().parents[2]
    base = root / "tests/regression/triad_bccFe_exchange"
    scratch = args.scratch_root.resolve()
    scratch.mkdir(parents=True, exist_ok=True)

    cases: dict[str, dict] = {}
    for nk in (4, 8, 12, 16):
        name = f"k{nk}_eta002"
        cases[name] = run_case(binary, base, scratch / name, nk, 0.02)
    for eta in (0.01, 0.04):
        name = f"k12_eta{int(eta * 1000):03d}"
        cases[name] = run_case(binary, base, scratch / name, 12, eta)
    # The recursion route uses this window for Chebyshev scaling.  A too-small
    # window (e.g. [-0.8,0.8]) is intentionally outside the supported regime
    # and can signal in the fast kernel; compare two safe windows instead.
    cases["k12_eta002_wide_window"] = run_case(binary, base, scratch / "k12_eta002_wide_window", 12, 0.02, -1.5, 1.5)

    failures: list[str] = []
    for name, case in cases.items():
        if case["summary"]["dyson_max"] > 1.0e-8:
            failures.append(f"{name}: Lehmann/Dyson Sigma=0 invariant {case['summary']['dyson_max']:.3e} > 1e-8")
        for pair in PAIR_KEYS:
            for row in case["rows"][pair]:
                for field in ("rs", "lehmann", "dyson", "rs_down", "lehmann_down"):
                    if not (math.isfinite(row[field].real) and math.isfinite(row[field].imag)):
                        failures.append(f"{name}/{pair}: non-finite {field}")

    k4, k8, k12, k16 = (cases["k4_eta002"], cases["k8_eta002"], cases["k12_eta002"], cases["k16_eta002"])
    k_convergence = {
        "k4_vs_k16_max_abs_G": max_probe_delta(k4, k16, "lehmann"),
        "k8_vs_k16_max_abs_G": max_probe_delta(k8, k16, "lehmann"),
        "k12_vs_k16_max_abs_G": max_probe_delta(k12, k16, "lehmann"),
        "k4_vs_k16_max_abs_G_down": max_probe_delta(k4, k16, "lehmann_down"),
        "k8_vs_k16_max_abs_G_down": max_probe_delta(k8, k16, "lehmann_down"),
        "k12_vs_k16_max_abs_G_down": max_probe_delta(k12, k16, "lehmann_down"),
    }

    eta_cases = [cases["k12_eta010"], cases["k12_eta002"], cases["k12_eta040"]]
    eta_convergence = {
        "eta001_vs_eta002_max_abs_G": max_probe_delta(eta_cases[0], eta_cases[1], "lehmann"),
        "eta002_vs_eta004_max_abs_G": max_probe_delta(eta_cases[1], eta_cases[2], "lehmann"),
        "eta001_vs_eta002_max_abs_G_down": max_probe_delta(eta_cases[0], eta_cases[1], "lehmann_down"),
        "eta002_vs_eta004_max_abs_G_down": max_probe_delta(eta_cases[1], eta_cases[2], "lehmann_down"),
    }
    window_delta = {
        "onsite_max_abs_G": max_probe_delta(cases["k12_eta002"], cases["k12_eta002_wide_window"], "lehmann", ("1_1",)),
        "intersite_max_abs_G": max_probe_delta(cases["k12_eta002"], cases["k12_eta002_wide_window"], "lehmann", ("1_335",)),
    }

    rs_comparison = {
        "eta002_k12_onsite_max_abs_G": max(cases["k12_eta002"]["summary"]["block_max"], 0.0),
        "eta002_k12_intersite_max_abs_G": cases["k12_eta002"]["pair_block_max"]["1_335"],
        "eta002_k12_dos_max_abs": cases["k12_eta002"]["summary"]["dos_max"],
        "eta002_k12_dos_rms": cases["k12_eta002"]["summary"]["dos_rms"],
        "eta002_k12_spectral_weight_rs": cases["k12_eta002"]["summary"]["weight_rs"],
        "eta002_k12_spectral_weight_lehmann": cases["k12_eta002"]["summary"]["weight_lehmann"],
    }
    result = {
        "scope": "bcc-Fe fixed potential; orthonormal full-BZ reciprocal GF; selected G11 and Gdd elements",
        "pairs": {"1_1": "onsite G_ii", "1_335": "first-neighbour intersite G_ij"},
        "probe_energies": PROBE_ENERGIES,
        "k_convergence": k_convergence,
        "eta_convergence": eta_convergence,
        "energy_window_sensitivity": window_delta,
        "rs_comparison_matched_eta002": rs_comparison,
        "runs": {name: {key: value for key, value in case.items() if key != "rows"} for name, case in cases.items()},
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if failures:
        print("VAL-05 FAIL")
        print("\n".join(failures))
        return 1
    print("VAL-05 PASS: direct G_ii/G_ij, k/eta/window convergence, and Sigma=0 D==E")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
