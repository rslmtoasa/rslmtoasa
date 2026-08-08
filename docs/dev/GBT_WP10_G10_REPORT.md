# WP10 / G10 report — legacy deletion and documentation

**Date:** 2026-08-08. **Branch:** `fable_v2_gbt_v2`.

## G9 status — explicitly not passing

WP10 is gated on G9 in `GBT_RS_LMTO_task_prompts.md`. G9 does **not** pass:
three items in `docs/dev/GBT_WP9_FAIL_WALKTHROUGH.md` remain open (a
commensurate-supercell band-energy residual 5-6x over tolerance, a
non-monotonic k-mesh refinement step at q3=0.5, and an unexplained ~776x
cone-angle scaling spread). Anders explicitly instructed proceeding with
WP10 anyway and skipping the G9 requirement. This report does not claim G9
passed, and none of WP10's work resolves or touches the WP9 physics
findings — they remain exactly as documented.

## Audit: most named deletion targets were already gone

Before touching any file, a repo-wide audit (grep-verified directly, not
taken on faith) checked WP10's six named deletion targets plus the two
"keep and document"/"add a repository check" requirements against the
actual state of `source/`, `tests/`, `docs/`, and `example/`. Result:

| WP10 checklist item | Finding |
| --- | --- |
| Absolute-position bulk GBT in `ham0m_nc` | **Already removed.** `ham0m_nc` fatal-errors if `magnetic_representation == gbt_single_q`; no `lattice%cr`/phase construction remains in it. Fixed in WP4 (moved to the ΔR-only recipe). |
| `fourier_transform_gbt*` | **Already removed** (WP5). Confirmed absent from `source/` by grep; only referenced by the test that checks it stays absent. |
| CCOR full-angle GBT rotations | **Already removed** (WP6b, commit `0188a6a`). Zero `cos(`/`sin(` calls remain in `hamiltonian_ccor.f90`; GBT CCOR routes through `build_ccor_pair_block_gbt`/`gbt_lift_orbital_block`, enforced by `tests/unit/test_gbt_wp6_ccor.py`'s token ban. |
| Mutually exclusive cone/sublattice routing | **Already fixed.** `gbt_endpoint_angles` selects the per-endpoint angle *source* (sublattice override vs. global), but `gbt_endpoint_link` is called unconditionally with `q_ss` regardless — the two mechanisms compose, they don't gate each other. |
| `gbt_kspace` deprecation/removal | **Still present** as a dead field (physics-inert since WP5). **Removed in this task** — see below. |
| Obsolete comments/known-wrong goldens | Source comments: clean. `tests/gbt_wp0/` (tracked, orphaned, one entry self-labeled `"diagnostic_known_wrong_gbt_not_golden"`, not wired into CMake, predates the representation-mode split) was the one real leftover — **deleted in this task**. |

So the actual remaining work was narrow: remove `gbt_kspace`, delete
`tests/gbt_wp0/`, add the q·bond repository check, fix two stale claims in
an untracked doc draft, and update three doc files. Everything else in the
table above is reported here as a verified "already clean" finding, not
silently assumed.

## Changes made

**Source (`gbt_kspace` removed outright, on Anders' explicit choice over
keeping it as a rejected-if-set field):**
- `source/hamiltonian.f90` — deleted the `logical :: gbt_kspace` type field.
- `source/include_codes/namelists/hamiltonian.f90` — deleted its namelist
  declaration and its entry in `namelist /hamiltonian/`.
- `source/hamiltonian_build.f90` — deleted the getter, setter, the
  deprecation-warning block that forced it to `.false.`, and the
  default-init site (four sites total).
- No example deck or test fixture set `gbt_kspace` (grep-confirmed before
  deleting), so this has no fixture fallout. An old deck that still sets it
  will now get a namelist parse error — a loud failure instead of a silent
  ignore, consistent with the project's existing fail-loud philosophy for
  unsupported GBT combinations.

**New repository check (WP10's "only `gbt_bond_phase` calculates q·bond"
requirement):**
- `tests/unit/test_gbt_wp10_source_contract.py` — scans every
  `source/*.f90` except `source/gbt_structure.f90` (which owns
  `gbt_bond_phase`/`gbt_endpoint_link`) for the component-indexed pattern
  `q_ss(1)`/`q_ss(2)`/`q_ss(3)`, a hand-rolled q·bond-style phase. The one
  legitimate exception — `explicit_texture`'s own absolute site-position
  phase in `prepare_explicit_texture_moments` (a different physical
  quantity: a site phase, not a bond gauge) — is named explicitly, and the
  test asserts the pattern is *still present* there (so the exemption
  itself can't silently go stale) before excluding just that routine's body
  from the scan. Registered in `CMakeLists.txt` as `UnitGbtWp10SourceContract`
  (labels `unit;gbt;wp10`).
- Verified the regex distinguishes real violations from clean code directly
  (`q_ss(1) * 2.0_rp` matches, `q_ss( 3 )` matches, unrelated text doesn't).

**Test/fixture cleanup:**
- `tests/unit/test_gbt_wp5_source_contract.py` — its old "`gbt_kspace` may
  only appear in these 3 files" check is now a straightforward "must not
  appear anywhere in `source/`" assertion; the now-dead `compatibility_block`
  slice check (which depended on code just deleted) was removed rather than
  left as vacuously-passing dead logic.
- `tests/gbt_wp6_matrix/feature_matrix.json` — removed the `gbt_kspace_flag`
  row; its cited guard message no longer exists in source. Confirmed safe by
  reading `matrix_contract()` in `test_gbt_wp6_matrix.py`: the anti-drift
  check only requires source-side GBT diagnostics be *covered* by some row,
  never the reverse, so deleting an inactive row cannot break it.
- `tests/unit/test_gbt_wp6_matrix.py` — updated one stale comment that
  referenced the now-deleted row as an example.
- `tests/gbt_wp0/` (`README.md`, `baselines.json`, `cases.json`,
  `cases/.gitkeep`) — deleted via `git rm` (tracked, fully recoverable from
  history). Confirmed not referenced by `CMakeLists.txt` or any other test
  before deleting.

**Docs:**
- `docs/source/keywords/hamiltonian_parameters.rst`,
  `docs/source/keywords/index.rst`, `docs/source/keywords/frozen_magnon.rst`
  — removed the `gbt_kspace` keyword entries; the first now notes it was
  removed in WP10 and points at `magnetic_representation='gbt_single_q'`.
- `docs/dev/GBT_OVERVIEW.md` (Anders' untracked draft) — two factual
  corrections only, not a rewrite: the K5 table row now separates the
  routing-exclusivity bug (fixed) from the unrelated, still-open
  auto-branch Goldstone-sum-rule gap; the K7 row (CCOR full-angle rotation)
  is now marked fixed (commit `0188a6a`) instead of "tracked, not yet
  revisited", with the same correction propagated to the §7 open-items list.
- `GBT_RS_LMTO_task_prompts.md` — WP10 checklist boxes checked off with an
  inline outcome note, matching how WP6c/WP6-integration/WP7 recorded their
  results.
- Examples: confirmed already clean — no deck references `gbt_kspace` or any
  other deleted mechanism.

**Not touched (explicitly out of scope):**
- `example/bulk/bccFe/q050/` — untracked, appears to be Anders' own
  in-progress run output; not safe to delete and not part of this task.
- Dated historical reports (`GBT_WP5_G5_REPORT.md`,
  `GBT_WP6_INTEGRATION_REPORT.md`, `PHASE3_STATUS.md`,
  `GBT_S_LEVEL_NEXT_SESSION_HANDOFF.md`, existing `tests/KNOWN_ISSUES.md`
  entries, `GBT_RS_LMTO_pre_B1r_migration_notes.md`) — point-in-time
  records; rewriting them would be revisionist.

## Test results

Rebuilt `build_13` (gfortran-13, Release) after the source/type changes —
compiled cleanly, no warnings related to the removed field.

- **`ctest -L unit`: 23/23 passed**, including the new
  `UnitGbtWp10SourceContract` and the updated `UnitGbtWp5SourceContract`/
  `UnitGbtWp6Matrix`.
- **`ctest -L regression`: 18/18 passed**, including `WP8LittleGroup` and
  `WP9MultisublatticeGoldstone` (the latter passes its loose real-space
  sanity bound — a measurement, not a strict gate, per the WP9 walkthrough).
- **`ctest -L example`: 37/45 passed.** The 8 failures are the exact
  pre-existing, already-documented set — verified against
  `tests/KNOWN_ISSUES.md`, not assumed:
  - `Example_bulk_bccFe_nsp4_block_spiral_qplus`,
    `Example_bulk_bccFe_nsp4_block_spiral_qminus`,
    `Example_frozen_magnon_bccFe`, `Example_frozen_magnon_bccFe_auto`,
    `Example_frozen_magnon_bccFe_auto_scf` — all fatal on
    `gbt_single_q requires strux_backend='strux_lib'; the legacy backend is
    unsupported.`, exactly the "Five `tests/scf/cases.json` GBT fixtures
    fatal on `strux_backend='legacy'`" issue in `tests/KNOWN_ISSUES.md`
    (reproduces identically on the pre-WP6c commit `0188a6a`, unrelated to
    this task).
  - `Example_bulk_bccFe_nsp2_block`, `Example_bulk_bccFe_nsp2_block_hoh` —
    a single `totaldos.out` row off by ~4-6e-5 (rel ~2e-6), matching the
    documented pre-existing gfortran-13 DOS-tolerance delta.
  - `Example_orbital_modern_bccFe` — times out at exactly 300.01s, matching
    the documented genuine timeout that reproduces identically on a
    pristine pre-WP7 tree at `b0e8e01`.

No new failure appeared anywhere in the three suites relative to the
documented pre-existing baseline. `grep -rn "gbt_kspace" source/` returns
nothing.

## G10: PASS

Every WP10 completion-checklist item is satisfied, with the specific
scoping noted above (five items were already done by earlier WPs; two were
done in this task). G9 remains FAIL and is untouched by this work — the
three WP9 physics-validation items are still open and are not claimed as
resolved anywhere in this report.
