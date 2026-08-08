#!/usr/bin/env python3
"""WP9 Battery A -- generalized-Bloch-theorem commensurate-supercell known-answer
test: MFT (force theorem), full SCF, and constrained-SCF legs.

WHAT IS BEING TESTED
====================
For a *commensurate* spiral wavevector the generalized Bloch theorem is not an
approximation: the rotating-frame Hamiltonian K(dR) built in the small chemical
cell is the exact site-local gauge transform of the lab-frame Hamiltonian of an
explicit magnetic supercell,

    K = U^dagger H_lab U ,    U = diag_R D(q.R) ,   D(a) = diag(e^-ia/2, e^+ia/2)

(plan section 2.4, spec (*)).  U is unitary and site-diagonal with U = 1 at the
central site, so every *on-site* Green function -- hence per-site charge,
per-site moment and the band energy -- must come out identical on the two
sides.  This is an identity, not a convergence statement, and that is what
makes it a known-answer test.

This is a WP9 (physical validation) port of the B1r.5 known-answer test
originally built on commit `3fd21c0` (sibling branch `fable_v2_gbt`, not an
ancestor of the branch this file lives on). See the module docstring of the
staged reference, formerly at
`_b1r5_reference/run_gbt_supercell.py` (now removed -- ported here in full),
for the original derivation of supercell sizing and the antiperiodic-spinor
trap; that derivation is unchanged and summarized below. **This file adds two
legs the original did not have** (SCF and constrained-SCF) and re-derives all
tolerances fresh, because the underlying kernel changed materially between
`3fd21c0` and the branch this runs on (WP2-WP8 rewrote large parts of the
Hamiltonian assembly, k-space route, symmetry handling and density contract).

THREE LEGS
==========
1. **mft** -- the original comparison. Frozen potential, `nstep=1`, `beta=1.0`,
   `magbeta=0` (force theorem: one band evaluation, E_F re-found at fixed
   electron number, plan section 3.2). Decks: `<case>/gbt/`, `<case>/super/`.
2. **scf** -- new. Both sides converge to full self-consistency (`nstep=25`,
   `mixtype='broyden'`, `beta=0.15`, `magbeta` at its default of 1.0, i.e. the
   moment direction is allowed to relax). Decks: `<case>/gbt_scf/`,
   `<case>/super_scf/`.
3. **constrained** -- new. The `gbt_scf` deck re-run with
   `&control density_policy = 'constrained_spiral'` set explicitly (this is
   already the code default -- see `source/control.f90` restore_to_default --
   set explicitly here for documentation and to future-proof against the
   default changing), compared against the SAME `super_scf` result used by
   leg 2. The per-site transverse residual and constraint torque
   (`source/spin_density.f90`, `source/bands.f90:resolve_density_axes`) are
   extracted from the `DENSITY_POLICY constrained_spiral` log lines and
   reported as diagnostics (not asserted against a tolerance -- there is no
   known-answer target for them here). Decks: `<case>/gbt_constrained/`,
   reusing `<case>/super_scf/`.

Real-space recursion only (`recur='block'`) on all three legs, per the
project's own sign-off on B1r.5 ("Force theorem (band energy comparison) not
stable for k-space") and this task's explicit instruction to stay on that
route.

IMPORTANT MEASURED CAVEAT -- read before trusting a PASS
==========================================================
Porting these decks surfaced a real, large discrepancy that this runner
reports rather than hides: the MFT leg's frozen starting potential is carried
over unmodified from `3fd21c0`, which was self-consistent under a materially
different (older) GBT kernel. Evaluated once (`nstep=1`) under the current
kernel, the `gbt/` side reads a moment within noise of zero, versus the
`super/` side's real ~2.4 mu_B -- a large, reproducible gap, not a tolerance
miss. A full SCF re-relaxation from the same starting potential recovers most
but not all of the supercell's moment. See `tests/KNOWN_ISSUES.md`
("`magnetic_representation = 'gbt_single_q'`: a frozen potential carried over
from a different kernel evaluates to near-zero moment...") for the full
measured record, including what was verified directly versus what is
inferred. **This runner and its fixtures are the validation harness; the
numbers below are what the harness measured on the code as it stands, not a
claim that the identity holds.**

SUPERCELL SIZE, AND THE ANTIPERIODIC-SPINOR TRAP (plan section 3.3)
===================================================================
The spiral phase is alpha(R) = 2*pi*q_ss.R with R in alat units (locked by
tests/unit/test_qss_theta_conventions.f90).  bcc (001) layers are spaced
*half* an alat apart, so q_ss = (0,0,0.5) advances the moment by 90 deg per
layer, not 180 -- the repeat is FOUR layers, not two, and q_ss = (0,0,1/3)
repeats after SIX layers, not three.

The naive "L = 1/q_ss cells" sizing (2 and 3) is exactly the section-3.3 trap
in its real-space form: those cells have winding number N_w = q_ss.a_3 = 1/2,
i.e. the gauge D(q.R) returns to -1 rather than +1 after one period.  The
rotating-frame problem tolerates that (it is the antiperiodic spinor boundary
condition, handled by half-shifting the k-mesh); an *explicit lab-frame*
supercell does not -- a half-turn cell is simply not a periodic magnetic
structure and there is no set of hand-set moments that realises it.  The cells
used here therefore take the shortest z translation with an INTEGER winding
number:

    q_ss = (0,0,0.5)   -> a_3 = (0,0,2) alat, 4 atoms, alpha = 0, 90, 180, 270
    q_ss = (0,0,1/3)   -> a_3 = (0,0,3) alat, 6 atoms, alpha = 0, 60, ..., 300

Both have N_w = 1 -- the odd winding number section 3.3 warns about -- so the
odd-N_w case is the one actually exercised here, in the form where it is a
theorem rather than a nuisance.

WP9-SPECIFIC DECK ADAPTATIONS
==============================
The `3fd21c0` decks predate the WP3/WP4 representation-mode split and the
WP4 strux backend requirement. Every `gbt*/` deck here adds, relative to the
original:

    &lattice  strux_backend = 'strux_lib'   ! WP4: gbt_single_q is fatal on
                                             ! the legacy backend
    &hamiltonian  magnetic_representation = 'gbt_single_q'   ! WP3/WP5: plain
                                             ! bulk SCF (pre_processing=
                                             ! 'bravais') does not auto-force
                                             ! this the way
                                             ! post_processing_frozen_magnon
                                             ! does; it defaults to
                                             ! periodic_nc, under which q_ss
                                             ! is silently ignored.

The `super*/` decks reach their per-site moments through the ordinary
`periodic_nc` path (each atom is its own crystal type with its own fixed
`mom(:)`) and need neither key. Whether `super/` also needs
`strux_backend='strux_lib'` for a clean backend-matched comparison, versus
staying on `legacy`, was tested empirically -- see the WP9 battery report for
the measured backend sensitivity (or lack of it).

TOLERANCES
==========
Derived fresh from two independent noise measurements on THIS runner, THIS
binary, reported by every invocation:

(1) INTRA-RUN SYMMETRY SPREAD. Every site of the supercell is equivalent to
    every other by the spiral symmetry, and the GBT run has only one site.
    The spread of per-site charge/moment *within one run* is pure numerical
    noise of the recursion + DOS integration + Fermi search.
(2) PERTURBED STARTING GUESS. Each side is run twice, the second time from a
    perturbed starting guess (spin-up band centre/enu/c pushed by
    `PERTURB_DC` Ry, starting Fermi level shifted by `PERTURB_DEF` Ry,
    applied identically to both sides). The spread of the GBT-minus-supercell
    *difference* between the two starts is the second noise floor.

See the WP9 battery report for the actual measured numbers and the
tolerances derived from them (`charge_tol`/`moment_tol`/`eband_tol` in
`cases.json`).

USAGE
=====
    run_wp9_commensurate_supercell.py --binary <rslmto.x> --cases-json <json>
                         --case-name <name> --scratch-root <dir> [--report]

`--report` prints the measured numbers and skips the assertions; use it when
re-deriving tolerances. `cases.json` entries carry a `"leg"` key
(`"mft"|"scf"|"constrained"`) selecting which deck-pair convention applies.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

# --- perturbed starting guess (applied identically to both sides) -------------
PERTURB_DC = 0.01      # Ry added to the spin-up band centre / enu / c
PERTURB_DEF = -0.02    # Ry shift of the starting Fermi level

_BAND_RE = re.compile(r"^(\s*(?:center_band|enu|c)\(:, 1\)\s*=\s*)(.*?)(\s*)$", re.M)
_FERMI_RE = re.compile(r"^(\s*fermi\s*=\s*)([-\d.eEdD+]+)", re.M)
_EBAND_RE = re.compile(r"Band energy of system:\s*([-\d.]+)")
_PROJ_RE = re.compile(r"Spin moment projections of atom\s+(\d+):\s*(.*)")
_DENSITY_POLICY_RE = re.compile(
    r"DENSITY_POLICY constrained_spiral atom\s*(\d+)\s+m_long=\s*([-\d.]+)\s+"
    r"\|m_transverse\|=\s*([-\d.]+)\s+torque=\(\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\)"
)


def load_case(cases_json: Path, case_name: str) -> dict:
    data = json.loads(cases_json.read_text())
    for case in data["cases"]:
        if case["name"] == case_name:
            return case
    raise KeyError(f"case {case_name!r} not found in {cases_json}")


def setup_workdir(base_dir: Path, workdir: Path, perturb: bool) -> None:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base_dir, workdir)
    if not perturb:
        return
    for path in sorted(workdir.glob("Fe*.nml")):
        def bump(m: re.Match) -> str:
            vals = [float(v) + PERTURB_DC for v in m.group(2).split(",")]
            return m.group(1) + ", ".join("%.10f" % v for v in vals) + m.group(3)

        path.write_text(_BAND_RE.sub(bump, path.read_text()))
    inp = workdir / "input.nml"
    inp.write_text(_FERMI_RE.sub(
        lambda m: m.group(1) + "%.6f" % (float(m.group(2).replace("d", "e")) + PERTURB_DEF),
        inp.read_text()))


def run_side(binary: Path, base_dir: Path, workdir: Path, timeout: int, perturb: bool) -> dict:
    setup_workdir(base_dir, workdir, perturb)
    result = subprocess.run(
        [str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=timeout,
    )
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(f"ERROR: {base_dir.name} exited {result.returncode}")
    return extract(workdir)


def extract(workdir: Path) -> dict:
    """Per-site charge/moment from the output potentials, band energy from the
    report, and (if present) the constrained-spiral density-policy diagnostic.
    """
    report = (workdir / "report.out").read_text()
    charges, moments = [], []
    for path in sorted(workdir.glob("Fe*_out.nml"),
                       key=lambda p: int(re.search(r"Fe(\d*)_out", p.name).group(1) or 0)):
        rows = dict(re.findall(r"ql\(1, :, (\d)\) = (.*)", path.read_text()))
        up = sum(float(v) for v in rows["1"].split(","))
        down = sum(float(v) for v in rows["2"].split(","))
        charges.append(up + down)
        moments.append(up - down)
    directions = []
    for _, comps in sorted(_PROJ_RE.findall(report), key=lambda kv: int(kv[0])):
        directions.append([float(v) for v in comps.split()])
    natom = len(charges)
    log_path = workdir / "run.log"
    density_policy = []
    if log_path.exists():
        for m in _DENSITY_POLICY_RE.finditer(log_path.read_text()):
            density_policy.append({
                "atom": int(m.group(1)),
                "m_long": float(m.group(2)),
                "m_transverse": float(m.group(3)),
                "torque": [float(m.group(4)), float(m.group(5)), float(m.group(6))],
            })
    return {
        "natom": natom,
        "charge": charges,
        "moment": moments,
        "direction": directions,
        "eband_per_atom": float(_EBAND_RE.search(report).group(1)) / natom,
        "density_policy": density_policy,  # last-iteration entries are last in the list
    }


def spread(values: list[float]) -> float:
    return max(values) - min(values) if len(values) > 1 else 0.0


def imposed_angles(q_ss: float, natom: int) -> list[float]:
    """alpha_k = 2*pi*q_ss.R_k in degrees; R_k is the k-th (001) layer at z = k/2 alat."""
    return [math.degrees(2 * math.pi * q_ss * 0.5 * k) % 360.0 for k in range(natom)]


def angle_error(direction: list[float], expected_deg: float) -> float:
    got = math.degrees(math.atan2(direction[1], direction[0])) % 360.0
    return min(abs(got - expected_deg), 360.0 - abs(got - expected_deg))


def compare(gbt: dict, sup: dict) -> dict:
    """GBT has one site; every supercell site must match it."""
    return {
        "d_eband": abs(gbt["eband_per_atom"] - sup["eband_per_atom"]),
        "d_charge": max(abs(g - s) for g in gbt["charge"] for s in sup["charge"]),
        "d_moment": max(abs(g - s) for g in gbt["moment"] for s in sup["moment"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--cases-json", type=Path, required=True)
    ap.add_argument("--case-name", required=True)
    ap.add_argument("--scratch-root", type=Path, required=True)
    ap.add_argument("--report", action="store_true",
                    help="print measurements and skip the assertions")
    ap.add_argument("--starts", default="base,perturbed",
                    help="comma-separated starting guesses to run "
                         "(default both; pass 'base' alone to skip the "
                         "perturbed-restart noise-floor measurement)")
    args = ap.parse_args()

    case = load_case(args.cases_json, args.case_name)
    base = args.cases_json.parent / case["base"]
    scratch = args.scratch_root / case["name"]
    scratch.mkdir(parents=True, exist_ok=True)
    timeout = int(case.get("timeout", 900))
    gbt_dir = case.get("gbt_dir", "gbt")
    super_dir = case.get("super_dir", "super")
    starts = [s.strip() for s in args.starts.split(",") if s.strip()]

    runs: dict[tuple[str, str], dict] = {}
    for start in starts:
        perturb = start == "perturbed"
        for side, dirname in (("gbt", gbt_dir), ("super", super_dir)):
            runs[(side, start)] = run_side(args.binary, base / dirname,
                                           scratch / f"{side}_{start}", timeout, perturb)

    leg = case.get("leg", "mft")
    print(f"=== {case['name']} [{leg}]: q_ss = (0,0,{case['q_ss']}) 2*pi/alat, "
          f"{runs[('super', starts[0])]['natom']}-atom supercell, theta_ss = 90 deg")

    failures: list[str] = []

    # --- setup guard: the supercell really is the spiral we meant to build ----
    expected = imposed_angles(float(case["q_ss"]), runs[("super", starts[0])]["natom"])
    errs = [angle_error(d, a) for d, a in zip(runs[("super", starts[0])]["direction"], expected)]
    print("  supercell moment directions (deg): imposed "
          + " ".join("%.1f" % a for a in expected)
          + "  | max deviation of the computed moment %.2e deg" % max(errs))
    if max(errs) > case["angle_tol"]:
        failures.append("supercell moment directions deviate from the imposed spiral by "
                        f"{max(errs):.3e} deg > {case['angle_tol']} deg")

    # --- noise floor (1): intra-run symmetry spread ---------------------------
    sym_q = max(spread(runs[(s, t)]["charge"]) for s in ("gbt", "super") for t in starts)
    sym_m = max(spread(runs[(s, t)]["moment"]) for s in ("gbt", "super") for t in starts)
    print(f"  noise floor (1) intra-run symmetry spread: dq = {sym_q:.2e} e, dm = {sym_m:.2e} mu_B")

    # --- the comparison, at every starting guess run -----------------------------
    diffs = {}
    for start in starts:
        d = compare(runs[("gbt", start)], runs[("super", start)])
        diffs[start] = d
        g, s = runs[("gbt", start)], runs[("super", start)]
        print(f"  [{start:9s}] E_band/atom  gbt {g['eband_per_atom']:.12f}  "
              f"super {s['eband_per_atom']:.12f}   diff {d['d_eband']:.2e} Ry")
        print(f"  [{start:9s}] charge       gbt {g['charge'][0]:.10f}  "
              f"super {' '.join('%.10f' % c for c in s['charge'])}")
        print(f"  [{start:9s}] moment       gbt {g['moment'][0]:.10f}  "
              f"super {' '.join('%.10f' % m for m in s['moment'])}")
        print(f"  [{start:9s}] max |gbt - super|: dq = {d['d_charge']:.2e} e, "
              f"dm = {d['d_moment']:.2e} mu_B")
        if g["density_policy"]:
            last = g["density_policy"][-1]
            print(f"  [{start:9s}] gbt DENSITY_POLICY (last iteration): "
                  f"m_long={last['m_long']:.6f} |m_transverse|={last['m_transverse']:.6f} "
                  f"torque=({last['torque'][0]:.6f},{last['torque'][1]:.6f},{last['torque'][2]:.6f})")

    # --- noise floor (2): reproducibility of the difference itself ------------
    if len(starts) > 1:
        print("  noise floor (2) spread of the gbt-super difference over the starting guesses: "
              f"dq = {spread([d['d_charge'] for d in diffs.values()]):.2e} e, "
              f"dm = {spread([d['d_moment'] for d in diffs.values()]):.2e} mu_B, "
              f"dE = {spread([d['d_eband'] for d in diffs.values()]):.2e} Ry")

    for start, d in diffs.items():
        for key, tol_key, unit in (("d_charge", "charge_tol", "e"),
                                   ("d_moment", "moment_tol", "mu_B"),
                                   ("d_eband", "eband_tol", "Ry/atom")):
            if d[key] > case[tol_key]:
                failures.append(f"[{start}] {key} = {d[key]:.3e} {unit} "
                                f"exceeds {tol_key} = {case[tol_key]:.1e}")

    if args.report:
        print("REPORT mode: assertions skipped")
        return 0
    if failures:
        print("FAIL:")
        for f in failures:
            print("  - " + f)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
