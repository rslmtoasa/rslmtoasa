#!/usr/bin/env python3
"""Read-only RA4 summary/comparison for completed RS-LMTO run directories.

This utility extracts registered scalar observables and d-projected LDOS.
It does not apply acceptance tolerances or infer scientific pass/fail.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def number(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def parse_report(path: Path, site: int) -> dict:
    txt = path.read_text(errors="replace")
    def scalar(label: str) -> float | None:
        m = re.search(rf"(?im)^\s*{re.escape(label)}\s*([-+\d.eEdD]+)", txt)
        return number(m.group(1)) if m else None

    def atom_scalar(label: str) -> float | None:
        m = re.search(rf"(?im)^\s*{re.escape(label)}\s*atom\s*{site}\s*:\s*({FLOAT})", txt)
        return number(m.group(1)) if m else None

    def orbital_occupation(spin: str) -> float | None:
        m = re.search(rf"(?im)^\s*{spin} orbital occupation at atom\s*{site}\s*:\s*(.*)$", txt)
        if not m:
            return None
        vals = [number(x) for x in re.findall(FLOAT, m.group(1))]
        return vals[2] if len(vals) >= 3 else None

    return {
        "site_index": site,
        "sphere_moment_muB": atom_scalar("Spin moment of"),
        "d_occupation_up_report": orbital_occupation("Up"),
        "d_occupation_down_report": orbital_occupation("Down"),
        "total_rms_diff_report_rounded": scalar("Total RMS Diff:"),
    }


def parse_potential(path: Path) -> dict:
    txt = path.read_text(errors="replace")

    def vector(pattern: str) -> list[float] | None:
        m = re.search(pattern + r"\s*=\s*([^\n]+)", txt, re.I)
        return [number(x) for x in re.findall(FLOAT, m.group(1))] if m else None

    # For l=2, flatten_ldm stores the 5x5 matrix row-major. Its diagonal
    # positions are 1, 7, 13, 19, 25 (Fortran 1-based indexing).
    traces = {}
    ql_d = {}
    for spin in (1, 2):
        # Full-precision ql(1,:,spin) is emitted in potential files; its last
        # component is the d occupation. Prefer it to report.out's 6 decimals.
        ql = vector(rf"ql\(\s*1\s*,\s*:\s*,\s*{spin}\s*\)")
        if ql:
            ql_d["up" if spin == 1 else "down"] = ql[-1]
        real = vector(rf"ldm_flatten\(\s*3\s*,\s*{spin}\s*,\s*:\s*\)")
        imag = vector(rf"ldm_flatten_imag\(\s*3\s*,\s*{spin}\s*,\s*:\s*\)")
        if real is not None and len(real) >= 25:
            diag = (0, 6, 12, 18, 24)
            traces["up" if spin == 1 else "down"] = {
                "real_trace": sum(real[i] for i in diag),
                "imag_trace": sum(imag[i] for i in diag) if imag and len(imag) >= 25 else None,
            }
    u, j = vector(r"hubbard_u"), vector(r"hubbard_j")
    active = bool((u and any(abs(x) > 1.0e-10 for x in u)) or (j and any(abs(x) > 1.0e-10 for x in j)))
    mode = re.search(r"(?im)^\s*hubbard_radial_moments\s*=\s*['\"]?([^,'\"\s]+)", txt)
    backend = re.search(r"(?im)^\s*hubbard_backend\s*=\s*['\"]?([^,'\"\s]+)", txt)
    return {
        "potential_file": str(path),
        "potential_sha256": sha256(path),
        "hubbard_u": u,
        "hubbard_j": j,
        "hubbard_active": active,
        "radial_mode": mode.group(1) if mode else None,
        "backend": backend.group(1) if backend else None,
        "d_occupation_up_ql": ql_d.get("up"),
        "d_occupation_down_ql": ql_d.get("down"),
        "ldm_d_trace_per_spin": traces,
        "ldm_trace_note": "inactive U/J; zero matrix is not an occupation oracle" if not active else "active U/J",
    }


def parse_ldos(path: Path) -> dict:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        vals = [number(x) for x in re.findall(FLOAT, line)]
        if len(vals) >= 19:
            # File columns: E, nine orbitals for spin up, nine for spin down.
            # The spd basis order has d orbitals at 5:9 in each spin block.
            rows.append((vals[0], sum(vals[5:10]), sum(vals[14:19])))
    if len(rows) < 2:
        raise ValueError(f"expected at least two 19-column LDOS rows in {path}")
    return {
        "path": str(path),
        "sha256": sha256(path),
        "energy_axis": "E - E_F (Ry), as emitted by bands.f90",
        "d_dos_units": "states/Ry/site per spin channel",
        "energy": [r[0] for r in rows],
        "d_up": [r[1] for r in rows],
        "d_down": [r[2] for r in rows],
    }


def status_from_log(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    clean = re.sub(r"\x1b\[[0-9;]*m", "", txt)
    records: list[dict] = []
    current = None
    for line in clean.splitlines():
        step = re.search(r"Perform recursion at step\s+(\d+)", line)
        if step:
            current = {"iteration": int(step.group(1)), "fixed_u_residual": None, "scf_marker": None, "scf_delta": None}
            records.append(current)
        fixed = re.search(
            rf"FIXED_U_RESIDUAL\s+rp_abs=\s*({FLOAT}).*?rv_abs_Ry=\s*({FLOAT})"
            rf".*?target_rp=\s*({FLOAT}).*?target_rv_Ry=\s*({FLOAT})"
            rf".*?physical=(\.true\.|\.false\.).*?streak=(\d+)", line, re.I
        )
        if fixed and current is not None:
            current["fixed_u_residual"] = {
                "rp_abs": number(fixed.group(1)), "rv_abs_Ry": number(fixed.group(2)),
                "target_rp": number(fixed.group(3)), "target_rv_Ry": number(fixed.group(4)),
                "physical": fixed.group(5).lower() == ".true.", "streak": int(fixed.group(6)),
            }
        marker = re.search(rf"(Not converged! Diff=|Converged!)\s*({FLOAT})?", line)
        if marker and current is not None:
            current["scf_marker"] = marker.group(1).strip()
            current["scf_delta"] = number(marker.group(2)) if marker.group(2) else None
    completed = [r for r in records if r["scf_marker"]]
    latest = completed[-1] if completed else None
    partial = records[-1] if records and (not latest or records[-1] is not latest) else None
    return {
        "run_log": str(path),
        "run_log_sha256": sha256(path),
        "scf_marker_last_completed_iteration": latest["scf_marker"] if latest else None,
        "scf_delta_last_completed_iteration": latest["scf_delta"] if latest else None,
        "scf_iterations_completed": latest["iteration"] if latest else None,
        "fixed_u_residual_same_iteration": latest["fixed_u_residual"] if latest else None,
        "trailing_iteration_in_progress": (
            {"iteration": partial["iteration"], "fixed_u_residual": partial["fixed_u_residual"]}
            if partial else None
        ),
    }


def choose_potential(run: Path, requested: str | None) -> Path:
    if requested:
        p = Path(requested)
        if not p.is_absolute():
            p = run / p
        return p
    candidates = sorted(run.glob("*_out.nml"))
    candidates = [p for p in candidates if not p.name.lower().endswith("-bulk_out.nml")]
    if len(candidates) != 1:
        raise ValueError(f"choose --potential; found {[p.name for p in candidates]} in {run}")
    return candidates[0]


def summarize(run: Path, site: int, potential: str | None, ldos: str | None) -> tuple[dict, dict]:
    run = run.resolve()
    report, log, inp = run / "report.out", run / "run.log", run / "input.nml"
    if not report.exists() or not log.exists():
        raise FileNotFoundError(f"completed run requires report.out and run.log: {run}")
    pot = choose_potential(run, potential)
    if ldos:
        lp = Path(ldos)
        if not lp.is_absolute():
            lp = run / lp
    else:
        choices = sorted(run.glob("*_orbital_dos.out"))
        if len(choices) != 1:
            raise ValueError(f"choose --ldos; found {[p.name for p in choices]} in {run}")
        lp = choices[0]
    files = {p.name: sha256(p) for p in (inp, report, log) if p.exists()}
    files.update({p.name: sha256(p) for p in sorted(run.glob("*.nml")) if p.name not in files})
    for name in ("clust", "map", "outnewclu", "str.out", "sbar", "view.sbar", "lattice.nml"):
        p = run / name
        if p.is_file():
            files[name] = sha256(p)
    obs = parse_report(report, site)
    potential_obs = parse_potential(pot)
    obs.update({k: potential_obs[k] for k in ("hubbard_u", "hubbard_j", "hubbard_active", "radial_mode", "backend", "ldm_d_trace_per_spin", "ldm_trace_note", "d_occupation_up_ql", "d_occupation_down_ql")})
    summary = {
        "run_dir": str(run),
        "files_sha256": files,
        "hash_note": "hashes describe current on-disk files; input launch hashes must be taken from immutable pre-run copies",
        "potential": {"path": str(pot), "sha256": sha256(pot)},
        "status": status_from_log(log),
        "observables": obs,
        "ldos": {k: v for k, v in parse_ldos(lp).items() if k not in ("energy", "d_up", "d_down")},
    }
    return summary, parse_ldos(lp)


def interp(xs: list[float], ys: list[float], x: float) -> float:
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    frac = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + frac * (ys[hi] - ys[lo])


def compare_ldos(a: dict, b: dict) -> dict:
    ea, eb = a["energy"], b["energy"]
    low, high = max(ea[0], eb[0]), min(ea[-1], eb[-1])
    if not low < high:
        raise ValueError("LDOS energy grids have no common interval")
    n = min(len(ea), len(eb))
    grid = [low + (high - low) * i / (n - 1) for i in range(n)]
    up = [abs(interp(ea, a["d_up"], x) - interp(eb, b["d_up"], x)) for x in grid]
    down = [abs(interp(ea, a["d_down"], x) - interp(eb, b["d_down"], x)) for x in grid]
    total = [
        abs(
            interp(ea, a["d_up"], x) + interp(ea, a["d_down"], x)
            - interp(eb, b["d_up"], x) - interp(eb, b["d_down"], x)
        )
        for x in grid
    ]
    return {
        "energy_axis": "E - E_F (Ry)",
        "d_dos_units": "states/Ry/site per spin channel",
        "energy_interval_common": [low, high],
        "comparison_points": n,
        "max_abs_d_ldos_up": max(up),
        "max_abs_d_ldos_down": max(down),
        "max_abs_d_ldos_spin_sum": max(total),
        "grid_definition": "uniform grid over interval intersection, linearly interpolated; no extrapolation",
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--run", action="append", required=True, help="run directory; repeat twice with --compare")
    p.add_argument("--site", action="append", type=int, help="1-based report site, default 1")
    p.add_argument("--potential", action="append", help="potential filename or path, aligned with --run")
    p.add_argument("--ldos", action="append", help="orbital DOS filename or path, aligned with --run")
    p.add_argument("--compare", action="store_true", help="compare two runs and add registered differences")
    args = p.parse_args()
    if args.compare and len(args.run) != 2:
        p.error("--compare requires exactly two --run arguments")
    n = len(args.run)
    sites = args.site or [1] * n
    pots = args.potential or [None] * n
    ldoses = args.ldos or [None] * n
    if not (len(sites) == len(pots) == len(ldoses) == n):
        p.error("if supplied, --site/--potential/--ldos must appear once per --run")
    summaries, dos = [], []
    for r, s, pot, ldos in zip(args.run, sites, pots, ldoses):
        sm, d = summarize(Path(r), s, pot, ldos)
        summaries.append(sm)
        dos.append(d)
    result = {"analyzer": "RA4 read-only observables v2", "tolerances_applied": False, "runs": summaries}
    if args.compare:
        a, b = [s["observables"] for s in summaries]
        result["comparison"] = {
            "scalar_differences_abs": {
                k: (abs(a[k] - b[k]) if a.get(k) is not None and b.get(k) is not None else None)
                for k in ("sphere_moment_muB", "d_occupation_up_ql", "d_occupation_down_ql")
            },
            "rounded_report_differences_abs": {
                k: (abs(a[k] - b[k]) if a.get(k) is not None and b.get(k) is not None else None)
                for k in ("d_occupation_up_report", "d_occupation_down_report")
            },
            "ldos": compare_ldos(dos[0], dos[1]),
        }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
