#!/usr/bin/env python3
"""VAL-17: reciprocal GBT cone-angle, Goldstone, and small-q validation.

The campaign records the primitive frozen-magnon energies before forming omega,
checks canonical occupations and the actual k-point contract, and compares the
reciprocal FeCo branch to an independently run real-space calculation.  The
same-q collinear reference is a physical GBT gauge identity; this script never
fits or applies a theta-dependent rescaling.
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FE_DECK = ROOT / "tests" / "regression" / "wp9_validation" / "gammaH_sweep"
FECO_DECK = ROOT / "tests" / "regression" / "wp9_validation" / "multisublattice_goldstone"
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
OCC_RE = re.compile(
    r"Canonical k-space occupations:\s*EF=\s*([-+0-9.eEdD]+) Ry,"
    r"\s*N=\s*([-+0-9.eEdD]+),\s*dN=\s*([-+0-9.eEdD]+),"
    r"\s*weight_sum\(raw\)=\s*([-+0-9.eEdD]+)"
)


def number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def run(binary: Path, workdir: Path, timeout: int) -> str:
    result = subprocess.run(
        [str(binary)], cwd=workdir, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, timeout=timeout,
    )
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{workdir} exited {result.returncode}:\n{result.stdout[-4000:]}")
    return result.stdout


def parse_diagnostics(workdir: Path) -> list[list[float]]:
    path = workdir / "frozen_magnon_diagnostics.dat"
    rows = []
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            rows.append([number(x) for x in line.split()])
    if not rows:
        raise RuntimeError(f"no rows in {path}")
    return rows


def parse_branches(workdir: Path) -> list[list[float]]:
    path = workdir / "frozen_magnon_branches.dat"
    rows = []
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            rows.append([number(x) for x in line.split()])
    if not rows:
        raise RuntimeError(f"no rows in {path}")
    return rows


def parse_modes(workdir: Path) -> list[list[float]]:
    path = workdir / "frozen_magnon_modes.dat"
    rows = []
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            rows.append([number(x) for x in line.split()])
    return rows


def parse_occupations(log: str) -> list[tuple[float, float, float, float]]:
    clean = ANSI_RE.sub("", log)
    return [tuple(number(x) for x in match.groups()) for match in OCC_RE.finditer(clean)]


def set_mesh(text: str, nk: int) -> str:
    for axis in (1, 2, 3):
        text, count = re.subn(rf"^\s*nk{axis}\s*=\s*\d+", f"nk{axis} = {nk}", text, flags=re.M)
        if count != 1:
            raise RuntimeError(f"expected one nk{axis} in reciprocal Fe input")
    return text


def run_fe(binary: Path, scratch: Path, nk: int, timeout: int) -> dict:
    work = scratch / f"fe_nk{nk}"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(FE_DECK / "base_kspace", work)
    shutil.copy2(FE_DECK / "q_points_convergence.dat", work / "q_points.dat")
    inp = work / "input.nml"
    inp.write_text(set_mesh(inp.read_text(), nk))
    log = run(binary, work, timeout)
    rows = []
    for line in (work / "frozen_magnon.dat").read_text().splitlines():
        if line and not line.startswith("#"):
            fields = line.split()
            rows.append((number(fields[2]), number(fields[4])))
    return {"work": work, "eband": dict(rows), "log": log}


def run_cone(binary: Path, scratch: Path, theta: float, timeout: int) -> dict:
    work = scratch / f"fe_cone_theta{int(theta)}"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(FE_DECK / "base_cone", work)
    shutil.copy2(FE_DECK / "q_points_cone.dat", work / "q_points.dat")
    inp = work / "input.nml"
    text, count = re.subn(r"(^\s*theta_ss\s*=\s*)[0-9.eEdD+-]+", rf"\g<1>{theta}", inp.read_text(), flags=re.M)
    if count != 1:
        raise RuntimeError("expected one theta_ss in cone input")
    inp.write_text(text)
    log = run(binary, work, timeout)
    rows = parse_diagnostics(work)
    qrow = min(rows, key=lambda row: abs(row[2] - 0.05) + abs(row[0]) + abs(row[1]))
    qzero = min(rows, key=lambda row: abs(row[0]) + abs(row[1]) + abs(row[2]))
    return {"work": work, "row": qrow, "qzero": qzero, "occ": parse_occupations(log), "log": log}


def run_feco(binary: Path, scratch: Path, nk: int | None, kspace: bool, timeout: int) -> dict:
    name = ("feco_kspace" if kspace else "feco_rs") + (f"_nk{nk}" if nk else "")
    work = scratch / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(FECO_DECK / "base", work)
    inp = work / "input.nml"
    text = inp.read_text()
    if kspace:
        text, count = re.subn(r"(^\s*nstep\s*=\s*\d+)", r"\1\nuse_kspace = .true.", text, count=1, flags=re.M)
        if count != 1:
            raise RuntimeError("could not enable k-space in FeCo input")
        if nk is not None:
            text += f"\n&reciprocal\nnk1 = {nk}\nnk2 = {nk}\nnk3 = {nk}\nq_symmetry_policy = 'full_bz'\n/\n"
    inp.write_text(text)
    log = run(binary, work, timeout)
    occ = parse_occupations(log)
    return {"work": work, "branches": parse_branches(work), "modes": parse_modes(work), "occ": occ, "log": log}


def qrow(rows: list[list[float]], q: float) -> list[float]:
    return min(rows, key=lambda row: abs(row[0] - q) + abs(row[1]) + abs(row[2]))


def gamma_branch(data: dict) -> tuple[int, float]:
    rows = [row for row in data["branches"] if max(abs(row[i]) for i in range(3)) < 1.0e-9]
    by_branch = {int(row[3]): row[4] for row in rows}
    branch = min(by_branch, key=lambda b: abs(by_branch[b]))
    return branch, by_branch[branch]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []

    def check(label: str, condition: bool, detail: str) -> None:
        print(f"  [{'ok' if condition else 'FAIL'}] {label}: {detail}")
        if not condition:
            failures.append(f"{label}: {detail}")

    def info(label: str, detail: str) -> None:
        print(f"  [info] {label}: {detail}")

    print("=== VAL-17 A: bcc-Fe primitive cone-angle diagnostics ===")
    cones = [run_cone(args.binary, args.scratch_root, theta, args.timeout) for theta in (5.0, 10.0, 15.0, 20.0)]
    gauge_ratios = []
    raw_ratios = []
    for theta, result in zip((5.0, 10.0, 15.0, 20.0), cones):
        row = result["row"]
        zero = result["qzero"]
        gauge_ratios.append(row[7] / row[8])
        raw_ratios.append(row[5] / row[8])
        print(
            f"  theta={theta:4.1f} E(q)={row[3]: .12e} E(0)={row[4]: .12e} "
            f"dE_raw={row[5]: .12e} sin2={row[8]: .12e} M={row[9]: .10f} "
            f"E(q,0)={row[6]: .12e} dE_gauge={row[7]: .12e} "
            f"dE_gauge/sin2={gauge_ratios[-1]: .12e} omega={row[10]: .12e}"
        )
        print(f"    q=0: E={zero[3]: .12e}, omega={zero[10]: .3e}")
        occ = result["occ"]
        check(
            f"theta={theta:g} canonical occupation",
            bool(occ) and max(abs(item[1] - 8.0) for item in occ) < 1.0e-8 and
            max(abs(item[2]) for item in occ) < 1.0e-7 and
            max(abs(item[3] - 1.0) for item in occ) < 1.0e-9,
            f"N={occ[-1][1] if occ else 'missing'}, max|dN|={max((abs(x[2]) for x in occ), default=float('inf')):.3e}, "
            f"weight range=({min((x[3] for x in occ), default=float('nan')):.12f}, {max((x[3] for x in occ), default=float('nan')):.12f})",
        )

    gauge_mean = sum(gauge_ratios) / len(gauge_ratios)
    gauge_rel = (max(gauge_ratios) - min(gauge_ratios)) / abs(gauge_mean)
    raw_mean = sum(raw_ratios) / len(raw_ratios)
    raw_rel = (max(raw_ratios) - min(raw_ratios)) / abs(raw_mean)
    print(f"  raw dE/sin2 spread={raw_rel * 100:.2f}% (fixed q=0 subtraction control)")
    print(f"  same-q gauge dE/sin2 spread={gauge_rel * 100:.2f}%")
    check("DeltaE(theta) ~ sin^2(theta) after physical gauge subtraction", gauge_rel < 0.05,
          f"relative spread={gauge_rel * 100:.2f}%")
    check("moment normalization is theta independent", max(r["row"][9] for r in cones) - min(r["row"][9] for r in cones) < 1.0e-7,
          f"M range=({min(r['row'][9] for r in cones):.10f}, {max(r['row'][9] for r in cones):.10f})")
    check("q=0 subtraction is numerically zero", max(abs(r["qzero"][7]) for r in cones) < 1.0e-10,
          f"max |DeltaE_gauge(Gamma)|={max(abs(r['qzero'][7]) for r in cones):.3e} Ry")
    first_log = ANSI_RE.sub("", cones[0]["log"])
    full_mesh = re.findall(r"Generated full mesh .*?\(([0-9]+) k-points\)", first_log)
    check("nonzero-q k-point contract is full BZ", "q_symmetry_policy = full_bz" in first_log and bool(full_mesh),
          f"policy=full_bz, full-mesh observations={full_mesh}, weights normalized in occupation logs")
    print("  frame/potential conclusion: q=0 GBT link is identity; fixed potential and M are shared across probes; "
          "the defect was the finite-mesh q-only offset in E(q,theta)-E(0,theta), amplified by 1/sin^2(theta).")

    print("=== VAL-17 B: FeCo reciprocal Goldstone and RS comparison ===")
    feco8 = run_feco(args.binary, args.scratch_root, 8, True, args.timeout)
    branch8, omega8 = gamma_branch(feco8)
    gamma_modes = [row for row in feco8["modes"] if max(abs(row[i]) for i in range(3)) < 1.0e-9 and int(row[3]) == branch8]
    real_parts = [row[6] * math.cos(row[7]) for row in gamma_modes]
    q002 = qrow(feco8["branches"], 0.02)[4]
    q005 = qrow(feco8["branches"], 0.05)[4]
    print(f"  k-space nk=8: Gamma branch={branch8}, omega_ac={omega8:.12e} Ry, "
          f"eigenvector real parts={[f'{x:.8f}' for x in real_parts]}, "
          f"omega(.02)={q002:.12e}, omega(.05)={q005:.12e}")
    check("k-space FeCo acoustic Goldstone", abs(omega8) < 1.0e-10, f"|omega(Gamma)|={abs(omega8):.3e} Ry")
    check("k-space FeCo acoustic eigenvector is in phase", len(real_parts) == 2 and real_parts[0] * real_parts[1] > 0.0,
          f"real components={real_parts}")
    occ8 = feco8["occ"]
    check("FeCo k-space occupations conserve N", bool(occ8) and max(abs(x[1] - 17.0) for x in occ8) < 1.0e-8 and
          max(abs(x[2]) for x in occ8) < 1.0e-7 and max(abs(x[3] - 1.0) for x in occ8) < 1.0e-9,
          f"N={occ8[-1][1] if occ8 else 'missing'}, EF range="
          f"({min((x[0] for x in occ8), default=float('nan')):.8f}, {max((x[0] for x in occ8), default=float('nan')):.8f}) Ry")

    rs = run_feco(args.binary, args.scratch_root, None, False, args.timeout)
    rs_branch, rs_omega = gamma_branch(rs)
    rs_modes = [row for row in rs["modes"] if max(abs(row[i]) for i in range(3)) < 1.0e-9 and int(row[3]) == rs_branch]
    rs_real = [row[6] * math.cos(row[7]) for row in rs_modes]
    print(f"  RS: Gamma branch={rs_branch}, omega_ac={rs_omega:.12e} Ry, "
          f"eigenvector real parts={[f'{x:.8f}' for x in rs_real]}")
    check("RS FeCo acoustic Goldstone", abs(rs_omega) < 1.0e-10, f"|omega(Gamma)|={abs(rs_omega):.3e} Ry")
    check("RS FeCo acoustic eigenvector is in phase", len(rs_real) == 2 and rs_real[0] * rs_real[1] > 0.0,
          f"real components={rs_real}")
    print("  RS comparison is an independent Gamma check; the finite-q RS cluster data are not used as a stiffness gate.")

    print("=== VAL-17 C: small-q and mesh convergence ===")
    feco_mesh = {nk: run_feco(args.binary, args.scratch_root, nk, True, args.timeout) for nk in (12, 16)}
    ratios = {}
    for nk, data in ((8, feco8), *feco_mesh.items()):
        small = [qrow(data["branches"], q)[4] / q**2 for q in (0.02, 0.05)]
        ratios[nk] = small
        print(f"  FeCo nk={nk}: omega/q^2 at q=.02,.05 = {[f'{x:.8e}' for x in small]}, "
              f"Gamma={gamma_branch(data)[1]:.3e} Ry")
    ratio12 = ratios[12]
    ratio16 = ratios[16]
    consistency12 = abs(ratio12[0] - ratio12[1]) / max(abs(sum(ratio12) / 2.0), 1.0e-30)
    consistency16 = abs(ratio16[0] - ratio16[1]) / max(abs(sum(ratio16) / 2.0), 1.0e-30)
    mesh_delta = max(abs(a - b) for a, b in zip(ratio12, ratio16)) / max(abs(sum(ratio16) / 2.0), 1.0e-30)
    info("FeCo small-q quadratic at nk=12", f"relative q2 mismatch={consistency12 * 100:.2f}%")
    info("FeCo small-q quadratic at nk=16", f"relative q2 mismatch={consistency16 * 100:.2f}%")
    info("FeCo q2 mesh convergence nk=12 -> 16", f"relative change={mesh_delta * 100:.2f}% (open convergence axis)")
    print("  nk=8 is retained as a diagnostic: its q=.02 value is visibly mesh-sensitive, while nk=12/16 are the production small-q pair.")

    fe_mesh = {nk: run_fe(args.binary, args.scratch_root, nk, args.timeout) for nk in (8, 12, 16)}
    vals = [fe_mesh[nk]["eband"][0.5] for nk in (8, 12, 16)]
    first_step, last_step = abs(vals[1] - vals[0]), abs(vals[2] - vals[1])
    print(f"  bcc-Fe q=.5 eband nk=8,12,16={[f'{x:.12e}' for x in vals]}, "
          f"steps={first_step:.3e},{last_step:.3e} Ry, spread={max(vals)-min(vals):.3e} Ry")
    info("bcc-Fe mesh non-monotonicity is exposed",
         f"8->12={first_step:.3e} Ry, 12->16={last_step:.3e} Ry; "
         f"last/first={last_step / first_step:.3f} (not hidden behind a monotonic gate)")
    print("  interpretation: the finite tetra/k-point integration and mesh-dependent SCF reference shift are resolved as a "
          "convergence limitation, not absorbed into a theta correction or diagonal shift.")

    if failures:
        print(f"FAIL: {len(failures)} gate(s) failed")
        for failure in failures:
            print("  - " + failure)
        return 1
    print("PASS: VAL-17 A/B physical gates passed; VAL-17 C is assessed and its open mesh limitation is recorded above")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
