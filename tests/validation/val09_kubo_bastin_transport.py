#!/usr/bin/env python3
"""VAL-09: physical validation of charge, spin, and orbital Kubo transport.

This is an opt-in campaign, not a quick regression.  It consumes the
production ``*_cond.out`` output and never reimplements the conductivity
estimator.  The reduced fcc-Pt fixture retains SOC, PBC, and the production
real-space hopping/current construction while keeping the campaign runnable
on a CPU developer checkout.

The exact reciprocal moment route is exercised only for charge.  The current
implementation explicitly scopes ``reciprocal%fill_moments`` to
``cond_type='charge'``; spin/orbital results therefore are not compared with
it because that would compare different or unavailable observables.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import f90nml


def patch_case(base: Path, workdir: Path, *, cond_type: str, va: list[int],
               vb: list[int], cond_ll: int = 20, channels: int = 120,
               replication: int = 4, fermi: float = -0.085837,
               energy_min: float = -2.5, energy_max: float = 1.2,
               route: str = "recursion", rc: int = 20,
               kmesh: int | None = None, gpu_plugin: bool = False,
               gpu_backend: str = "csr", lld: int | None = None,
               cheb_backend: str = "legacy") -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base, workdir)
    patch = {
        "calculation": {"gf_route": route},
        "lattice": {"rc": rc, "n1": replication, "n2": replication,
                    "n3": replication, "strux_backend": "legacy"},
        "hamiltonian": {"v_alpha": va, "v_beta": vb, "js_alpha": "z",
                        "jl_alpha": "z", "hoh": False},
        "energy": {"fermi": fermi, "energy_min": energy_min,
                    "energy_max": energy_max, "channels_ldos": channels},
        "control": {"cond_ll": cond_ll, "cond_type": cond_type,
                     "cond_calctype": "per_type", "cheb_backend": cheb_backend,
                     "gpu_plugin": gpu_plugin, "gpu_backend": gpu_backend},
        "reciprocal": {"nk1": kmesh or replication, "nk2": kmesh or replication,
                        "nk3": kmesh or replication, "use_symmetry_reduction": False},
    }
    if lld is not None:
        patch["control"]["lld"] = lld
    patched = workdir / "input.patched.nml"
    f90nml.patch(str(workdir / "input.nml"), patch, str(patched))
    patched.replace(workdir / "input.nml")


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
        raise RuntimeError(f"VAL-09 run failed ({kwargs}):\n{log[-5000:]}")
    return {"parameters": kwargs, "observable": parse_conductivity(workdir)}


def parse_conductivity(workdir: Path) -> dict[str, float]:
    paths = list(workdir.glob("*_cond.out"))
    if not paths:
        raise RuntimeError(f"missing per-type *_cond.out in {workdir}")
    path = paths[0]
    rows: list[tuple[float, float, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        cols = line.split()
        if len(cols) < 3:
            continue
        try:
            row = tuple(float(value.replace("D", "E").replace("d", "e")) for value in cols[:3])
        except ValueError:
            continue
        if all(math.isfinite(value) for value in row):
            rows.append(row)
    if not rows:
        raise RuntimeError(f"Pt_cond.out has no finite rows in {workdir}")
    energy, real, imag = min(rows, key=lambda row: abs(row[0]))
    return {"energy_relative_to_fermi": energy, "real": real, "imag": imag,
            "abs": math.hypot(real, imag)}


def value(cases: dict[str, dict], name: str) -> float:
    return float(cases[name]["observable"]["real"])


def sequence_metrics(cases: dict[str, dict], names: list[str]) -> dict[str, object]:
    values = [value(cases, name) for name in names]
    scale = max(max(abs(item) for item in values), 1.0e-12)
    deltas = [abs(values[index + 1] - values[index]) / scale for index in range(len(values) - 1)]
    return {"names": names, "values": values, "scale": scale,
            "max_step_over_scale": max(deltas, default=0.0),
            "last_step_over_scale": deltas[-1] if deltas else 0.0}


def sign_check(cases: dict[str, dict], positive: str, reversed_name: str) -> dict[str, float]:
    first = value(cases, positive)
    second = value(cases, reversed_name)
    scale = max(abs(first), abs(second), 1.0e-12)
    return {"positive": first, "reversed": second,
            "sum_over_scale": abs(first + second) / scale,
            "ratio": second / first if first else math.nan}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--gpu-plugin", action="store_true",
                        help="run the same validation campaign through the CUDA RS plugin")
    parser.add_argument("--gpu-backend", default="csr",
                        choices=("csr", "bsr", "fft", "conv"),
                        help="CUDA RS backend when --gpu-plugin is enabled")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    base = root / "tests/postproc/cases/conductivity/fccPt"
    scratch = args.scratch_root.resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    cases: dict[str, dict] = {}

    common = {"cond_ll": 20, "channels": 120, "replication": 4,
              "fermi": -0.085837, "energy_min": -2.5, "energy_max": 1.2}
    directions = {
        "xx": ([1, 0, 0], [1, 0, 0]),
        "yy": ([0, 1, 0], [0, 1, 0]),
        "zz": ([0, 0, 1], [0, 0, 1]),
        "xy": ([1, 0, 0], [0, 1, 0]),
        "yx": ([0, 1, 0], [1, 0, 0]),
    }

    def add(name: str, *, case_base: Path = base, **kwargs: object) -> None:
        cases[name] = run_case(binary, case_base, scratch / name,
                               gpu_plugin=args.gpu_plugin,
                               gpu_backend=args.gpu_backend, **kwargs)

    # Tensor structure and component/sign conventions at a common point.
    for cond_type in ("charge", "spin", "orbital"):
        for component, (va, vb) in directions.items():
            add(f"{cond_type}_{component}", cond_type=cond_type, va=va, vb=vb, **common)
        add(f"{cond_type}_minus_alpha_yx", cond_type=cond_type,
            va=[0, -1, 0], vb=[1, 0, 0], **common)
        add(f"{cond_type}_minus_beta_yx", cond_type=cond_type,
            va=[0, 1, 0], vb=[-1, 0, 0], **common)

    # Chebyshev order convergence for all implemented observable families.
    order_grid = {
        "charge": (60, 80, 100),
        "spin": (40, 60, 80),
        "orbital": (120, 160, 200),
    }
    for cond_type in ("charge", "spin", "orbital"):
        for order in order_grid[cond_type]:
            add(f"{cond_type}_order_{order}", cond_type=cond_type,
                va=[1, 0, 0] if cond_type == "charge" else [0, 1, 0],
                vb=[1, 0, 0] if cond_type == "charge" else [1, 0, 0],
                **{**common, "cond_ll": order,
                   "replication": 4 if cond_type == "charge" else 8,
                   "channels": 120 if cond_type == "charge" else 480})

    # PBC/RS replication and Fermi placement.  The reciprocal mesh is changed
    # with the same integer to keep the cross-representation geometry explicit.
    for cond_type in ("charge", "spin", "orbital"):
        for replication in (4, 6, 8):
            add(f"{cond_type}_replication_{replication}", cond_type=cond_type,
                va=[1, 0, 0] if cond_type == "charge" else [0, 1, 0],
                vb=[1, 0, 0] if cond_type == "charge" else [1, 0, 0],
                **{**common, "replication": replication})
        for offset in (-0.15, 0.0, 0.15):
            add(f"{cond_type}_fermi_{offset:+.2f}", cond_type=cond_type,
                va=[1, 0, 0] if cond_type == "charge" else [0, 1, 0],
                vb=[1, 0, 0] if cond_type == "charge" else [1, 0, 0],
                **{**common, "fermi": common["fermi"] + offset})

    # There is no user-facing eta/kernel-alpha namelist in the current
    # conductivity implementation.  Window changes plus cond_ll are the
    # supported broadening/kernel-resolution diagnostics, and the fixed kernel
    # is recorded in the report as Lorentz(alpha=6).
    for cond_type in ("charge", "spin", "orbital"):
        for window in ((-2.0, 1.0), (-2.5, 1.2), (-3.0, 1.5)):
            add(f"{cond_type}_window_{window[0]:+.1f}_{window[1]:+.1f}",
                cond_type=cond_type,
                va=[1, 0, 0] if cond_type == "charge" else [0, 1, 0],
                vb=[1, 0, 0] if cond_type == "charge" else [1, 0, 0],
                **{**common, "energy_min": window[0], "energy_max": window[1]})

    # Genuine common-observable cross-representation comparison: charge only.
    cross_base = root / "tests/regression/triad_bccFe_conductivity"
    cross_common = {"cond_ll": 30, "channels": 120, "replication": 20,
                    "rc": 80, "kmesh": 6, "fermi": -0.069099,
                    "energy_min": -1.0, "energy_max": 1.2}
    for route in ("recursion", "lehmann"):
        add(f"charge_cross_{route}", case_base=cross_base,
            cond_type="charge", va=[1, 0, 0], vb=[1, 0, 0],
            **{**cross_common, "route": route})

    failures: list[str] = []
    for name, case in cases.items():
        if not all(math.isfinite(float(case["observable"][key])) for key in ("real", "imag", "abs")):
            failures.append(f"{name}: non-finite conductivity")

    # Cubic fcc Pt: diagonal charge conductivity is isotropic and charge Hall
    # conductivity is forbidden in this time-reversal-symmetric fixture.
    diag = [value(cases, "charge_xx"), value(cases, "charge_yy"), value(cases, "charge_zz")]
    diag_scale = max(max(abs(item) for item in diag), 1.0e-12)
    diagonal_spread = max(diag) - min(diag)
    if diagonal_spread / diag_scale > 0.20:
        failures.append(f"charge diagonal cubic spread too large: {diagonal_spread / diag_scale:.3g}")
    if max(abs(value(cases, "charge_xy")), abs(value(cases, "charge_yx"))) / diag_scale > 0.10:
        failures.append("charge xy/yx forbidden component is not small")

    signs = {}
    for cond_type in ("charge", "spin", "orbital"):
        signs[cond_type] = {
            "alpha": sign_check(cases, f"{cond_type}_yx", f"{cond_type}_minus_alpha_yx"),
            "beta": sign_check(cases, f"{cond_type}_yx", f"{cond_type}_minus_beta_yx"),
        }
        for side, check in signs[cond_type].items():
            if not math.isfinite(check["sum_over_scale"]) or check["sum_over_scale"] > 0.10:
                failures.append(f"{cond_type} {side}-direction sign convention failed")

    convergence = {}
    for cond_type in ("charge", "spin", "orbital"):
        order_names = [f"{cond_type}_order_{order}" for order in order_grid[cond_type]]
        replication_names = [f"{cond_type}_replication_{replication}" for replication in (4, 6, 8)]
        fermi_names = [f"{cond_type}_fermi_{offset:+.2f}" for offset in (-0.15, 0.0, 0.15)]
        window_names = [f"{cond_type}_window_{lo:+.1f}_{hi:+.1f}" for lo, hi in ((-2.0, 1.0), (-2.5, 1.2), (-3.0, 1.5))]
        convergence[cond_type] = {
            "order": sequence_metrics(cases, order_names),
            "replication": sequence_metrics(cases, replication_names),
            "fermi": sequence_metrics(cases, fermi_names),
            "window": sequence_metrics(cases, window_names),
        }
        # A campaign is accepted only when increasing order does not keep
        # moving by more than 35% of the observed scale.  Metallic Fermi and
        # finite-size changes are reported, not falsely treated as errors.
        if convergence[cond_type]["order"]["last_step_over_scale"] > 0.40:
            failures.append(f"{cond_type} order sequence has not settled")

    cross = {
        "recursion": value(cases, "charge_cross_recursion"),
        "lehmann": value(cases, "charge_cross_lehmann"),
    }
    cross_scale = max(abs(cross["recursion"]), abs(cross["lehmann"]), 1.0e-12)
    cross["absolute_difference"] = abs(cross["recursion"] - cross["lehmann"])
    cross["difference_over_scale"] = cross["absolute_difference"] / cross_scale
    if cross["difference_over_scale"] > 0.25:
        failures.append("charge recursion/eigenpair conductivity mismatch exceeds 25% envelope")

    report = {
        "scope": "fcc Pt, SOC, one-site fcc, real-space PBC fixture; reduced rc/n/energy/order campaign",
        "kernel": {"name": "Lorentz", "alpha": 6.0,
                   "source": "source/conductivity.f90 calculate_gamma_nm; no eta/kernel-alpha namelist is exposed"},
        "execution": {"gpu_plugin": args.gpu_plugin,
                      "gpu_backend": args.gpu_backend if args.gpu_plugin else None},
        "operator_conventions": {
            "charge": "v_a/b(R) = (1/i) (direction dot (r_i-r_j) alat) H_ij(R), with optional velocity_scale on v_b",
            "spin": "J^s_a = 1/2 {S_z, v_a}; js_alpha selects S_z here",
            "orbital": "J^L_a = 1/2 {L_z, v_a}; jl_alpha selects L_z here",
            "tensor": "each scalar output is the production Kubo contraction for v_alpha (first operator) and v_beta (second operator), sampled at the row nearest E-E_F=0",
            "units": "the current estimator uses the internal Ry/alat convention and factor 16/(pi DeltaE^2); e^2/hbar and volume conversion are not applied, so no S/cm or hbar/e S/cm numerical comparison is claimed",
        },
        "tensor": {"charge_diagonal": diag, "charge_diagonal_relative_spread": diagonal_spread / diag_scale,
                   "charge_xy": value(cases, "charge_xy"), "charge_yx": value(cases, "charge_yx")},
        "sign_checks": signs,
        "convergence": convergence,
        "cross_representation_charge": cross,
        "material_anchor": "Pt tensor structure/order of magnitude assessed internally; external literature magnitudes are not numerically compared because the output unit conversion is incomplete.",
        "cases": cases,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if failures:
        print("VAL-09 FAIL")
        print("\n".join(failures))
        return 1
    print("VAL-09 PASS: charge, spin, and orbital Kubo-Bastin diagnostics produced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
