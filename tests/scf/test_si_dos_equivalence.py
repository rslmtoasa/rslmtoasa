#!/usr/bin/env python3
"""Compare the production Si real-space and k-space DOS routes.

This is deliberately a comparison driver, not a DOS implementation. It runs
both production routes, reads their totaldos.out files, reconstructs absolute
energy from each route's reported Fermi level, and compares the spectra.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import time
from pathlib import Path

import f90nml


# Frozen after the one-time development study documented in the Si fixture
# README. These values must not be tuned during a test run.
CHEB_LLD = 200
KS_NK = 8
GAUSSIAN_SIGMA = 0.02
DOS_CHANNELS = 1000
ENERGY_MIN = -1.5
ENERGY_MAX = 1.0

EXPECTED_DOS_STATES = 16.0
EXPECTED_KS_OCCUPATIONS = 8.0
INTEGRAL_TOL = 0.02
WINDOW_COUNT_TOL = 0.10
DOS_RMS_REL_TOL = 0.20
DOS_MAX_ABS_TOL = 12.5
PEAK_POSITION_TOL = 0.08
WINDOW_MIN = -1.2
WINDOW_MAX = 0.8
PROBE_ENERGIES = (-0.8, -0.4, 0.0, 0.4, 0.6)


def patch_fixture(case_root: Path, workdir: Path, use_kspace: bool) -> None:
    shutil.copytree(case_root, workdir)
    patch = {
        "self": {"nstep": 1, "use_kspace": use_kspace},
        "energy": {
            "channels_ldos": DOS_CHANNELS,
            "energy_min": ENERGY_MIN,
            "energy_max": ENERGY_MAX,
        },
        "control": {
            "nsp": 1,
            "lmax": 1,
            # The k-space branch skips the recursion solver stage.  Keep the
            # real-space route on the stable functional backend; fast is
            # covered separately by the backend regression matrix.
            "lld": CHEB_LLD,
            "recur": "chebyshev",
            "cheb_backend": "legacy",
        },
    }
    if use_kspace:
        patch["reciprocal"] = {
            "nk1": KS_NK,
            "nk2": KS_NK,
            "nk3": KS_NK,
            "use_symmetry_reduction": False,
            "use_time_reversal": True,
            "use_shift": False,
            "n_energy_points": DOS_CHANNELS,
            "dos_method": "gaussian",
            "gaussian_sigma": GAUSSIAN_SIGMA,
        }
    patched = workdir / "input.patched.nml"
    f90nml.patch(str(workdir / "input.nml"), patch, out_path=str(patched))
    patched.replace(workdir / "input.nml")


def run_route(
    case_root: Path,
    workdir: Path,
    runner: Path,
    binary: Path,
    use_kspace: bool,
) -> float:
    patch_fixture(case_root, workdir, use_kspace)
    started = time.monotonic()
    # MSYS Python launches child processes through the native Windows API,
    # where a bare ``bash`` can resolve outside the configured MSYS2 runtime.
    # Match run_test.py and pass the MSYS-provided interpreter explicitly.
    bash = shutil.which("bash") if os.name == "nt" else "/bin/bash"
    if bash is None:
        raise RuntimeError("bash is required to run production routes on Windows/MSYS2")
    result = subprocess.run(
        [bash, str(runner), str(binary)],
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
    )
    elapsed = time.monotonic() - started
    if result.returncode != 0:
        log = workdir / "testrun.log"
        detail_parts = []
        if log.exists():
            detail_parts.append("--- testrun.log ---\n" + log.read_text(errors="replace")[-4000:])
        if result.stdout:
            detail_parts.append("--- runner stdout ---\n" + result.stdout[-4000:])
        if result.stderr:
            detail_parts.append("--- runner stderr ---\n" + result.stderr[-4000:])
        detail = "\n".join(detail_parts) or "(no output captured)"
        raise RuntimeError(f"production route failed in {workdir}\n{detail}")
    return elapsed


def parse_fermi(log_path: Path, use_kspace: bool) -> float:
    marker = "Canonical k-space occupations: EF=" if use_kspace else "Free Fermi energy:"
    for line in log_path.read_text(errors="replace").splitlines():
        if marker in line:
            value = line.split(marker, 1)[1].split()[0]
            return float(value.replace("D", "E").replace("d", "e"))
    raise RuntimeError(f"could not find {marker!r} in {log_path}")


def parse_ks_occupations(log_path: Path) -> tuple[float, float]:
    marker = "Canonical k-space occupations:"
    for line in log_path.read_text(errors="replace").splitlines():
        if marker in line:
            tail = line.split(marker, 1)[1]
            electron_count = float(tail.split("N=", 1)[1].split(",", 1)[0])
            band_energy = float(tail.split("EBAND=", 1)[1].split(",", 1)[0])
            return electron_count, band_energy
    raise RuntimeError(f"could not find {marker!r} in {log_path}")


def read_dos(path: Path, fermi: float) -> list[tuple[float, float]]:
    dos: list[tuple[float, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 2 or fields[0].startswith("#"):
            continue
        try:
            energy = float(fields[0]) + fermi
            value = float(fields[1])
        except ValueError:
            continue
        if not (math.isfinite(energy) and math.isfinite(value)):
            raise RuntimeError(f"non-finite DOS value in {path}: {line}")
        dos.append((energy, value))
    if len(dos) < 50:
        raise RuntimeError(f"too few DOS rows in {path}: {len(dos)}")
    if any(left[0] >= right[0] for left, right in zip(dos, dos[1:])):
        raise RuntimeError(f"DOS energy grid is not strictly increasing in {path}")
    return dos


def interpolate(dos: list[tuple[float, float]], energy: float) -> float:
    if energy < dos[0][0] or energy > dos[-1][0]:
        raise ValueError(f"energy {energy} is outside DOS grid")
    for (e0, d0), (e1, d1) in zip(dos, dos[1:]):
        if e0 <= energy <= e1:
            fraction = (energy - e0) / (e1 - e0)
            return d0 + fraction * (d1 - d0)
    return dos[-1][1]


def integrate(
    dos: list[tuple[float, float]], lower: float | None = None, upper: float | None = None
) -> float:
    if lower is None or upper is None:
        points = dos
    else:
        points = [(lower, interpolate(dos, lower))]
        points.extend((energy, value) for energy, value in dos if lower < energy < upper)
        points.append((upper, interpolate(dos, upper)))
    return sum(
        (e1 - e0) * (d0 + d1) / 2.0
        for (e0, d0), (e1, d1) in zip(points, points[1:])
    )


def central_comparison(
    rs_dos: list[tuple[float, float]], ks_dos: list[tuple[float, float]]
) -> tuple[float, float, float]:
    common = [
        (energy, interpolate(ks_dos, energy), value)
        for energy, value in rs_dos
        if WINDOW_MIN <= energy <= WINDOW_MAX
    ]
    if len(common) < 50:
        raise RuntimeError("too few points in the common central DOS window")
    differences = [ks_value - rs_value for _, ks_value, rs_value in common]
    rms = math.sqrt(sum(delta * delta for delta in differences) / len(differences))
    scale = math.sqrt(sum(rs_value * rs_value for _, _, rs_value in common) / len(common))
    maximum = max(abs(delta) for delta in differences)
    return rms / max(scale, 1.0e-12), maximum, scale


def peak_energy(dos: list[tuple[float, float]]) -> float:
    central = [point for point in dos if WINDOW_MIN <= point[0] <= WINDOW_MAX]
    return max(central, key=lambda point: point[1])[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.case_root = args.case_root.resolve()
    args.scratch_root = args.scratch_root.resolve()

    runner = Path(__file__).resolve().parents[1] / "run_binary.sh"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    rs_workdir = args.scratch_root / "real_space_chebyshev"
    ks_workdir = args.scratch_root / "k_space_gaussian"
    for workdir in (rs_workdir, ks_workdir):
        if workdir.exists():
            shutil.rmtree(workdir)

    total_started = time.monotonic()
    rs_seconds = run_route(args.case_root, rs_workdir, runner, args.binary, False)
    ks_seconds = run_route(args.case_root, ks_workdir, runner, args.binary, True)

    rs_fermi = parse_fermi(rs_workdir / "testrun.log", False)
    ks_fermi = parse_fermi(ks_workdir / "testrun.log", True)
    ks_electrons, ks_band_energy = parse_ks_occupations(ks_workdir / "testrun.log")
    if abs(ks_electrons - EXPECTED_KS_OCCUPATIONS) > 1.0e-8:
        raise RuntimeError(f"unexpected canonical k-space electron count: {ks_electrons}")
    if not math.isfinite(ks_band_energy):
        raise RuntimeError(f"non-finite canonical k-space band energy: {ks_band_energy}")

    rs_dos = read_dos(rs_workdir / "totaldos.out", rs_fermi)
    ks_dos = read_dos(ks_workdir / "totaldos.out", ks_fermi)
    rs_integral = integrate(rs_dos)
    ks_integral = integrate(ks_dos)
    if abs(rs_integral - EXPECTED_DOS_STATES) > INTEGRAL_TOL:
        raise RuntimeError(f"real-space DOS state count mismatch: {rs_integral}")
    if abs(ks_integral - EXPECTED_DOS_STATES) > INTEGRAL_TOL:
        raise RuntimeError(f"k-space DOS state count mismatch: {ks_integral}")
    if abs(rs_integral - ks_integral) > INTEGRAL_TOL:
        raise RuntimeError(f"DOS state-count mismatch between routes: {rs_integral} vs {ks_integral}")

    rs_window_count = integrate(rs_dos, WINDOW_MIN, WINDOW_MAX)
    ks_window_count = integrate(ks_dos, WINDOW_MIN, WINDOW_MAX)
    if abs(rs_window_count - ks_window_count) > WINDOW_COUNT_TOL:
        raise RuntimeError(f"central-window state-count mismatch: {rs_window_count} vs {ks_window_count}")

    rms_relative, maximum_difference, rms_scale = central_comparison(rs_dos, ks_dos)
    if rms_relative > DOS_RMS_REL_TOL:
        raise RuntimeError(f"central DOS RMS mismatch: {rms_relative} > {DOS_RMS_REL_TOL}")
    if maximum_difference > DOS_MAX_ABS_TOL:
        raise RuntimeError(f"central DOS maximum mismatch: {maximum_difference} > {DOS_MAX_ABS_TOL}")

    rs_peak = peak_energy(rs_dos)
    ks_peak = peak_energy(ks_dos)
    if abs(rs_peak - ks_peak) > PEAK_POSITION_TOL:
        raise RuntimeError(f"principal DOS peak mismatch: {rs_peak} vs {ks_peak}")

    print(
        f"fixed parameters: Chebyshev lld={CHEB_LLD}, KS mesh={KS_NK}x{KS_NK}x{KS_NK}, "
        f"Gaussian sigma={GAUSSIAN_SIGMA:.3f} Ry"
    )
    print(f"finite DOS rows: real-space={len(rs_dos)}, k-space={len(ks_dos)}")
    print(f"Fermi levels: real-space={rs_fermi:.8f}, k-space={ks_fermi:.8f} Ry")
    print(f"canonical k-space occupations: N={ks_electrons:.8f}, EBAND={ks_band_energy:.8f} Ry")
    print(f"integrated DOS states: real-space={rs_integral:.8f}, k-space={ks_integral:.8f}")
    print(f"central-window states: real-space={rs_window_count:.8f}, k-space={ks_window_count:.8f}")
    print(
        f"central DOS: relative RMS={rms_relative:.6f}, max abs={maximum_difference:.6f}, "
        f"scale={rms_scale:.6f}"
    )
    for energy in PROBE_ENERGIES:
        rs_value = interpolate(rs_dos, energy)
        ks_value = interpolate(ks_dos, energy)
        print(f"DOS probe {energy:+.2f} Ry: real-space={rs_value:.6f}, k-space={ks_value:.6f}")
    print(f"principal DOS peak: real-space={rs_peak:.6f}, k-space={ks_peak:.6f} Ry")
    print(
        f"route runtime: real-space={rs_seconds:.2f} s, k-space={ks_seconds:.2f} s, "
        f"total={time.monotonic() - total_started:.2f} s"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, OSError) as error:
        print(f"FAIL: {error}")
        raise SystemExit(1)
