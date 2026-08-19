#!/usr/bin/env python3
"""B5.2 route-agnostic estimator triad runner.

Runs one bcc-Fe post-processing observable (J_ij or sigma) three ways --
recursion vs Lehmann vs Dyson(Sigma=0) -- on the SAME converged potential, and
checks the documented route-agreement envelopes:

  * lehmann == dyson  (tight): Sigma=0 makes backend D reproduce backend E, so the
    two k-space routes agree to solver tolerance. This is the strongest cross-route
    statement and the permanent invariant.
  * recursion vs lehmann (envelope): a DOCUMENTED band, not machine precision.
    For sigma it is tight (the exact k-space moments match the recursion moments up
    to k-mesh/Chebyshev truncation); for J_ij it is broad and shell-dependent (the
    two routes broaden G differently -- see docs/dev/reciprocal_green_convergence.md
    and docs/dev/route_agnostic_estimators.md).
  * golden reproducibility (informational): each route's value vs a stored golden,
    a loose guard against silent regressions (generous tol for cross-BLAS variance).

Mirrors run_matrix.py's mechanics (f90nml.patch, scratch copytree). Goldens live in
tests/regression/references_triad/<name>.json.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manifest_defaults import apply_manifest_defaults

ROUTE_PATCH = {
    "recursion": {"calculation": {"gf_route": "recursion"}},
    "lehmann": {"calculation": {"gf_route": "lehmann"}, "reciprocal": {"green_backend": "lehmann"}},
    "dyson": {"calculation": {"gf_route": "dyson"}, "reciprocal": {"green_backend": "dyson"}},
}

_JIJ_RE = re.compile(r"Jij between pair\s+(\d+)\s+and\s+(\d+)\s+is\s+([-\dEe.+]+)")


def load_case(cases_json: Path, case_name: str) -> dict:
    data = json.loads(cases_json.read_text())
    for case in data["cases"]:
        if case["name"] == case_name:
            return apply_manifest_defaults(data, case)
    raise KeyError(f"case {case_name!r} not found in {cases_json}")


def setup_workdir(base_dir: Path, workdir: Path) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base_dir, workdir)
    for pattern in ("*_out.nml", "cond_total.out", "run.log", "testrun.log"):
        for path in workdir.glob(pattern):
            path.unlink()


def run_route(binary: Path, base_dir: Path, workdir: Path, route: str, timeout: int,
               reciprocal_backend: str) -> str:
    setup_workdir(base_dir, workdir)
    tmp = workdir / "input.nml.tmp"
    patch = dict(ROUTE_PATCH[route])
    if route in ("lehmann", "dyson"):
        patch["reciprocal"] = {**patch.get("reciprocal", {}),
                                "reciprocal_backend": reciprocal_backend}
    f90nml.patch(str(workdir / "input.nml"), patch, str(tmp))
    tmp.replace(workdir / "input.nml")
    result = subprocess.run(
        [str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=timeout,
    )
    log = result.stdout
    (workdir / "run.log").write_text(log)
    if result.returncode != 0:
        print(log[-3000:])
        raise SystemExit(f"ERROR: route {route} exited {result.returncode}")
    return log


def extract_jij(workdir: Path, log: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for m in _JIJ_RE.finditer(log):
        out[f"{int(m.group(1))}_{int(m.group(2))}"] = float(m.group(3))
    if not out:
        raise SystemExit("ERROR: no 'Jij between pair' lines found in run log")
    return out


def extract_sigma(workdir: Path, fermi: float) -> dict[str, float]:
    path = workdir / "cond_total.out"
    if not path.exists():
        raise SystemExit("ERROR: cond_total.out not produced")
    best = None
    for line in path.read_text().splitlines():
        cols = line.split()
        if len(cols) < 2:
            continue
        try:
            e, s = float(cols[0]), float(cols[1])
        except ValueError:
            continue
        d = abs(e - fermi)
        if best is None or d < best[0]:
            best = (d, s)
    if best is None:
        raise SystemExit("ERROR: cond_total.out had no numeric rows")
    return {"sigma_xx": best[1]}


def extract_alpha(workdir: Path) -> dict[str, float]:
    # On-site Gilbert damping alpha = 0.5*(xx+yy), column 12 of alldampings.out
    # (cols: #i #j xx xy xz yx yy yz zx zy zz 0.5*(xx+yy) Dist rx ry rz).
    path = workdir / "alldampings.out"
    if not path.exists():
        raise SystemExit("ERROR: alldampings.out not produced (do_damping off, or nsp!=2?)")
    out: dict[str, float] = {}
    for line in path.read_text().splitlines():
        cols = line.split()
        if len(cols) < 12 or not cols[0].lstrip("-").isdigit():
            continue
        out[f"{int(cols[0])}_{int(cols[1])}"] = float(cols[11])
    if not out:
        raise SystemExit("ERROR: no numeric pair rows in alldampings.out")
    return out


def extract(case: dict, workdir: Path, log: str) -> dict[str, float]:
    if case["observable"] == "jij":
        return extract_jij(workdir, log)
    if case["observable"] == "sigma":
        return extract_sigma(workdir, float(case["fermi"]))
    if case["observable"] == "alpha":
        return extract_alpha(workdir)
    raise SystemExit(f"ERROR: unknown observable {case['observable']!r}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--cases-json", required=True)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--scratch-root", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--reciprocal-backend", choices=("lapack", "cuda"), default="lapack")
    parser.add_argument("--gen-ref", action="store_true")
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    cases_json = Path(args.cases_json).resolve()
    tests_dir = cases_json.parent
    case = load_case(cases_json, args.case_name)
    name = case["name"]
    timeout = int(case.get("timeout", 900))

    base_dir = tests_dir / case["base"]
    scratch = Path(args.scratch_root).resolve() / name

    values: dict[str, dict[str, float]] = {}
    for route in ("recursion", "lehmann", "dyson"):
        log = run_route(binary, base_dir, scratch / route, route, timeout, args.reciprocal_backend)
        values[route] = extract(case, scratch / route, log)

    ref_dir = Path(args.references).resolve()
    ref_path = ref_dir / f"{name}.json"
    if args.gen_ref:
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_path.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")
        print(f"WROTE [{name}]: {ref_path}")
        return 0

    failures: list[str] = []

    # (1) lehmann == dyson (Sigma=0), tight.
    ld_tol = float(case["lehmann_dyson_tol"])
    for key in values["lehmann"]:
        d = abs(values["lehmann"][key] - values["dyson"][key])
        if d > ld_tol:
            failures.append(f"lehmann!=dyson [{key}]: |{values['lehmann'][key]:.6g}-"
                            f"{values['dyson'][key]:.6g}|={d:.3g} > {ld_tol:g}")

    # (2) recursion vs lehmann envelope (documented band on the ratio; same sign).
    lo, hi = (float(x) for x in case["recursion_vs_lehmann_ratio"])
    for key in values["recursion"]:
        r, l = values["recursion"][key], values["lehmann"][key]
        if r == 0.0:
            continue
        ratio = l / r
        if (l < 0) != (r < 0) or not (lo <= ratio <= hi):
            failures.append(f"recursion/lehmann envelope [{key}]: lehmann/recursion="
                            f"{ratio:.4f} outside [{lo:g},{hi:g}] (rec={r:.6g} leh={l:.6g})")

    # (3) golden reproducibility (informational guard, generous tol).
    if ref_path.exists():
        golden = json.loads(ref_path.read_text())
        rtol = float(case.get("golden_rtol", 1.0e-2))
        atol = float(case.get("golden_atol", 1.0e-4))
        for route in values:
            for key, val in values[route].items():
                g = golden.get(route, {}).get(key)
                if g is None:
                    continue
                if abs(val - g) > atol + rtol*abs(g):
                    failures.append(f"golden drift [{route}/{key}]: {val:.6g} vs "
                                    f"{g:.6g} (rtol={rtol:g})")
    else:
        print(f"NOTE [{name}]: no golden {ref_path} (run with --gen-ref to create)")

    print(f"[{name}] values: " + json.dumps(values, sort_keys=True))
    if failures:
        for f in failures:
            print(f"  FAIL: {f}")
        print(f"FAIL [{name}]")
        return 1
    print(f"PASS [{name}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
