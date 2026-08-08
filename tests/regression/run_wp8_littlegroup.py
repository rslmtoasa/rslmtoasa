#!/usr/bin/env python3
"""WP8 -- q little-group reduction and q-lifecycle, known-answer tests.

WHAT IS BEING TESTED
=====================
A finite-q GBT spin spiral lowers the symmetry of the electronic problem to
the little group of q_ss (docs/dev/plans/B1_gbt_frozen_magnons_v2.md section
3.1). The full chemical-cell BZ (&reciprocal q_symmetry_policy = 'full_bz',
the default) always remains a valid oracle. 'little_group' and
'little_group_common' are opt-in reductions
(source/reciprocal_bands.f90::generate_little_group_kpoint_mesh,
source/reciprocal_lifecycle.f90::force_full_bz_for_nonzero_q_gbt,
source/reciprocal.f90::ensure_kpoint_mesh's mesh cache key).

This script runs one bcc-Fe k-space frozen_magnon('mft') deck (wp8_littlegroup/base,
adapted from the B1r.6 fixture) under several q-lists/policies and checks:

  1. axial   -- q on the bcc Gamma-H line (little group C4v): 'little_group'
                agrees with the 'full_bz' oracle and the mesh is strictly
                reduced.
  2. generic -- q off every symmetry axis/plane (trivial little group):
                'little_group' agrees with the oracle and the mesh is (near)
                unreduced.
  3. shifted -- the same axial q on a shifted Monkhorst-Pack mesh: agreement
                and reduction both survive the shift.
  4. q_and_minus_q / common_subgroup -- a 4-row sweep (q=0 reference, +q_axial,
                -q_axial, q_generic) run under 'full_bz', 'little_group', and
                'little_group_common': every policy must agree row-for-row
                with the oracle, 'little_group' must rebuild the mesh once per
                distinct nonzero q, and 'little_group_common' must build
                exactly one mesh for the whole sweep (never row 1's mesh
                reused blindly for the others).
  5. reordered -- the same 4 q-points with rows 2-4 shuffled: per-q results
                (matched by q-vector, not row index) must be unchanged from
                case 4 under both 'little_group' and 'little_group_common'.

USAGE
=====
    run_wp8_littlegroup.py --binary <rslmto.x> --deck-root <dir>
                            --scratch-root <dir>
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

_POLICY_RE = re.compile(r"^(\s*q_symmetry_policy\s*=\s*)'[^']*'", re.M)
_SHIFT_RE = re.compile(r"^(\s*use_shift\s*=\s*)\.\w+\.", re.M)
_LITTLE_GROUP_RE = re.compile(
    r"reduced by the little group common to\s+(\d+)\s+q-point\(s\)\s+to\s+(\d+)\s+"
    r"irreducible k-points from\s+(\d+)\s+total points"
)
EBAND_TOL = 1.0e-6  # Ry; force-theorem probe on an identical Hamiltonian.


def run_variant(binary: Path, base_dir: Path, workdir: Path, q_points: Path, policy: str,
                 timeout: int, use_shift: bool = False) -> dict:
    if workdir.exists():
        shutil.rmtree(workdir)
    shutil.copytree(base_dir, workdir)
    shutil.copy(q_points, workdir / "q_points.dat")

    inp = workdir / "input.nml"
    text = inp.read_text()
    text, n = _POLICY_RE.subn(r"\g<1>'" + policy + "'", text)
    if n != 1:
        raise SystemExit(f"ERROR: expected exactly one q_symmetry_policy key in {inp}, found {n}")
    flag = ".true." if use_shift else ".false."
    text, n = _SHIFT_RE.subn(r"\g<1>" + flag, text)
    if n != 1:
        raise SystemExit(f"ERROR: expected exactly one use_shift key in {inp}, found {n}")
    inp.write_text(text)

    result = subprocess.run(
        [str(binary)], cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=timeout,
    )
    (workdir / "run.log").write_text(result.stdout)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        raise SystemExit(f"ERROR: {workdir.name} exited {result.returncode}")

    fm_dat = workdir / "frozen_magnon.dat"
    rows = [line.split() for line in fm_dat.read_text().splitlines() if line and not line.startswith("#")]
    if not rows:
        raise SystemExit(f"ERROR: no data rows in {fm_dat}")
    # q1 q2 q3 etot eband mtot_1..mtot_nrec omega -- eband is column index 4
    by_q = {(round(float(r[0]), 6), round(float(r[1]), 6), round(float(r[2]), 6)): float(r[4]) for r in rows}

    return {
        "eband_by_q": by_q,
        "eband_row1": float(rows[0][4]),
        "little_group_events": [(int(m[0]), int(m[1]), int(m[2])) for m in _LITTLE_GROUP_RE.findall(result.stdout)],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--deck-root", type=Path, required=True)
    ap.add_argument("--scratch-root", type=Path, required=True)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()

    base = args.deck_root / "base"
    scratch = args.scratch_root
    scratch.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []

    def check(label: str, cond: bool, detail: str) -> None:
        status = "ok" if cond else "FAIL"
        print(f"  [{status}] {label}: {detail}")
        if not cond:
            failures.append(f"{label}: {detail}")

    # --- 1/2/3: single-q axial / generic / shifted-mesh reduction vs oracle ---
    for case, qfile, expect_reduction, use_shift in [
        ("axial", "q_points_axial.dat", True, False),
        ("generic", "q_points_generic.dat", False, False),
        ("shifted", "q_points_axial.dat", True, True),
    ]:
        print(f"=== {case} ===")
        full = run_variant(args.binary, base, scratch / f"{case}_full", args.deck_root / qfile,
                           "full_bz", args.timeout, use_shift)
        lg = run_variant(args.binary, base, scratch / f"{case}_lg", args.deck_root / qfile,
                         "little_group", args.timeout, use_shift)
        d_eband = abs(full["eband_row1"] - lg["eband_row1"])
        check(f"{case}: eband agreement", d_eband < EBAND_TOL,
              f"full_bz {full['eband_row1']:.10f}, little_group {lg['eband_row1']:.10f}, diff {d_eband:.2e} Ry")
        check(f"{case}: little-group event logged", len(lg["little_group_events"]) >= 1,
              f"{len(lg['little_group_events'])} event(s)")
        if lg["little_group_events"]:
            nk_irred, nk_full = lg["little_group_events"][0][1], lg["little_group_events"][0][2]
            if expect_reduction:
                check(f"{case}: mesh strictly reduced", nk_irred < nk_full,
                      f"{nk_irred} irreducible from {nk_full} total")
            else:
                check(f"{case}: mesh (near-)unreduced (trivial little group)", nk_irred == nk_full,
                      f"{nk_irred} irreducible from {nk_full} total")

    # --- 4: multi-q sweep (q=0, +q_axial, -q_axial, q_generic) under all three policies ---
    print("=== multi-q common-subgroup / q-and-minus-q ===")
    order_a = {}
    for policy in ("full_bz", "little_group", "little_group_common"):
        order_a[policy] = run_variant(args.binary, base, scratch / f"orderA_{policy}",
                                      args.deck_root / "q_points_orderA.dat", policy, args.timeout)

    oracle_q = order_a["full_bz"]["eband_by_q"]
    for policy in ("little_group", "little_group_common"):
        for q, e_oracle in oracle_q.items():
            e_test = order_a[policy]["eband_by_q"].get(q)
            check(f"orderA[{policy}] q={q}", e_test is not None and abs(e_test - e_oracle) < EBAND_TOL,
                  f"oracle {e_oracle:.10f}, {policy} {e_test}, "
                  f"diff {abs((e_test or 0) - e_oracle):.2e} Ry" if e_test is not None else "row missing")

    n_lg_events = len(order_a["little_group"]["little_group_events"])
    n_common_events = len(order_a["little_group_common"]["little_group_events"])
    check("little_group rebuilds once per distinct nonzero q", n_lg_events == 3,
          f"{n_lg_events} event(s) (expected 3: +q_axial, -q_axial, q_generic)")
    check("little_group_common builds exactly one mesh for the whole sweep", n_common_events == 1,
          f"{n_common_events} event(s) (never row 1's mesh reused blindly, never a per-q rebuild either)")

    # --- 5: reordered q-list (rows 2-4 shuffled; row 1 = q=0 reference is unchanged) ---
    print("=== reordered q-list ===")
    order_b = {}
    for policy in ("little_group", "little_group_common"):
        order_b[policy] = run_variant(args.binary, base, scratch / f"orderB_{policy}",
                                      args.deck_root / "q_points_orderB.dat", policy, args.timeout)
        for q, e_a in order_a[policy]["eband_by_q"].items():
            e_b = order_b[policy]["eband_by_q"].get(q)
            check(f"reorder[{policy}] q={q} unchanged", e_b is not None and abs(e_b - e_a) < EBAND_TOL,
                  f"orderA {e_a:.10f}, orderB {e_b}, "
                  f"diff {abs((e_b or 0) - e_a):.2e} Ry" if e_b is not None else "row missing")

    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed:")
        for f in failures:
            print("  - " + f)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
