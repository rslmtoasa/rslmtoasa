# GBT WP9 / Gate G9 report — physical validation

Date: 2026-08-07. Branch `fable_v2_gbt_v2`, on top of `936ea5a` (final G8:
PASS, `docs/dev/GBT_WP8_G8_REPORT.md`).

## Outcome

**Gate G9: FAIL**, with important qualifications. No production-code physics
was changed during this task — every finding below is a measurement, or (in
one case) a **test-fixture fix**, per the task's instruction ("do not change
physics during the first run... do not regenerate goldens or impose a
post-hoc Goldstone shift to pass").

- The algebraic/oracle layer (analytic Stoner, dense bond oracle),
  q-symmetry (q↔−q, smoothness through the zone boundary), and — the most
  encouraging result of this task — the multi-sublattice Goldstone check all
  measure clean.
- The commensurate-supercell test's headline failure (a near-zero reported
  moment) was traced to a **test-fixture bug**, not a kernel defect: the
  deck set the GBT-side `mom` input to the physical cone direction instead
  of the code's required collinear-rotating-frame convention. Fixed; the
  moment gap closed by 360–4100×. What remains is a smaller, real,
  still-open band-energy residual — the item still **fails** its tolerance,
  but for a much better-understood and smaller reason.
- Cone-angle (θ) scaling of ω(q) **fails** badly (~776× spread across
  5°–20°) and is not explained; time-reversal/BZ-reduction was checked and
  ruled out as the cause.
- The Γ–H sweep was reworked to **k-space only** per Anders' instruction
  (long-wavelength spirals can have real numerical problems on the
  real-space route from cluster-truncation/PBC artifacts) partway through
  this task. Under the new, k-space-only design, its mesh-refinement gate
  now fails mildly at one q-point (a real, small, uninvestigated
  non-monotonicity) — a different, better-targeted failure than the
  original mixed-route design's qualified pass.
- The LKAG stiffness comparison is **deferred** per Anders' explicit
  instruction, not attempted here.

Per the task's own ground rules, nothing here was tuned or hidden to pass.
Every finding is recorded in `tests/KNOWN_ISSUES.md` with what was verified
versus guessed, and registered as a `ctest` so it stays visible.

## 1. Execution model

WP9 was run as three independent, parallel validation batteries (matching
the task prompt's suggested delegation: "up to three ... runners; one ...
diagnosis owner"), followed by substantial integrator-led follow-up
diagnosis and fixture repair triggered by Anders' direct review of the
results:

- **Battery A** — commensurate-supercell known-answer test (MFT + SCF +
  constrained-SCF legs), `tests/regression/wp9_validation/commensurate_supercell/`.
- **Battery B** — bcc-Fe Γ–H frozen-magnon sweep, `tests/regression/wp9_validation/gammaH_sweep/`.
- **Battery C** — multi-sublattice Goldstone measurement (bcc FeCo),
  `tests/regression/wp9_validation/multisublattice_goldstone/`.

A mid-task editor/session restart killed Battery B's and Battery C's
in-flight numerical runs (setup work survived on disk; only the `rslmto.x`
runs were lost). Both were re-run directly by the integrator; Battery C's
re-run was seeded from an already-converged reference potential to keep the
restart cheap and avoid a mixing instability a cold start showed.

After the first full pass, Anders reviewed the results directly and asked
two pointed questions that materially changed this report's conclusions:

1. *Is the GBT formalism's `mom_i = (sinθcos(qr), sinθsin(qr), cosθ)`
   convention — constant projection along z — actually what the code
   implements, or could there be a local/global axis confusion?* This led
   to tracing the commensurate-supercell "near-zero moment" finding to a
   root cause (§3) and fixing the affected decks.
2. *The `gammaH_sweep` decks have `use_time_reversal = .true.` — GBT should
   require time-reversal off.* This led to explicitly checking (and, this
   time, confirming correct) the code's TR-handling for nonzero-q GBT (§3,
   §6).

Anders separately instructed that the Γ–H sweep be reworked to k-space only
(§6), which a delegated agent completed; the integrator reviewed and
independently reran that rework's `ctest`.

## 2. Item 1 — Analytic Stoner and dense bond oracle: **PASS**

No new work needed: `tests/unit/test_gbt_oracles.py` (`ctest`
`UnitGbtOracles`) already carries the one-orbital Stoner k∓q/2 identity and
the dense unequal-species pair oracle, established under WP1/G1 and
re-verified at every subsequent gate. Re-run on the current HEAD used for
this task: `ctest -L unit` **22/22**, including `UnitGbtOracles` and every
other GBT unit oracle. No regression.

## 3. Item 2 — Commensurate supercell, MFT then SCF: **FAIL** (root cause found, partially fixed)

`tests/regression/wp9_validation/commensurate_supercell/` — bcc Fe, θ=90°
flat spiral, q_ss=(0,0,0.5) → 4-atom supercell and q_ss=(0,0,1/3) → 6-atom
supercell (§3.3 sizing, integer winding number). Ported and adapted from a
fixture on the sibling branch `fable_v2_gbt` (commit `3fd21c0`, never merged
here — reused after review) built for an earlier, single-kernel GBT
architecture.

### 3.1 Root cause of the original "near-zero moment" failure — a test-fixture bug, now fixed

The MFT leg originally showed the GBT side's moment reading ≈0 against the
supercell's ≈2.4 μB — a six-orders-of-magnitude gap. Prompted by Anders'
question about local-vs-global axis handling, this was traced fully:

- In `gbt_single_q` (for these `hoh`-unset decks, where `build_obarm`/
  `build_enim` are never even called — confirmed empirically), the bond
  Hamiltonian and on-site correction are built **entirely** from
  `theta_ss`/`phi_ss_sublattice` plus the *sign* of `potential%mom(3)` (a
  binary up/down-sublattice selector, `source/hamiltonian_build.f90:1619`).
  `mom(1)`/`mom(2)` and the magnitude of `mom(3)` never enter the
  Hamiltonian construction.
- But `ql`'s up/down decomposition — what every moment readout (and SCF
  mixing) is built on — separately projects the density onto
  `potential%mom` as its quantization axis. For `gbt_single_q`, `mom` must
  therefore stay in the code's internal collinear-rotating-frame
  convention, nominally `(0,0,1)` **regardless of `theta_ss`** — exactly
  what the already-validated `tests/regression/wp8_littlegroup/base/Fe.nml`
  does even at `theta_ss=90`. The WP9 decks instead set `mom=(1,0,0)`,
  physically pre-rotated to the cone direction — the older architecture's
  convention, wrong for this one.
- **Direct proof:** with `mom=(1,0,0)`, `source/bands.f90`'s
  `DENSITY_POLICY` diagnostic (independent of `mom`) reads `m_long=0.000000
  |m_transverse|=2.395328` — the moment is present, magnitude ≈2.395 μB,
  matching the supercell almost exactly, just read out along the wrong
  axis. Editing only `mom` to `(0,0,1)` reproduces byte-identical
  `DENSITY_POLICY` physics (confirming the Hamiltonian doesn't depend on
  `mom`'s value at all here) but `ql` now correctly reads ≈2.395 μB.

**Fix applied** to all six `gbt_single_q` decks
(`gbt_supercell/{q050,q033}/gbt{,_scf,_constrained}/Fe.nml`):
`mom(:) = 0.0d0, 0.0d0, 1.0d0`. Effect on the registered ctests:

| case | moment gap before | moment gap after | eband gap (unaffected — frame-invariant) |
| --- | --- | --- | --- |
| q050 MFT | 2.396 μB | **5.78e-4 μB** (~4100× better) | 7.19e-4 Ry (unchanged) |
| q033 MFT | ~2.33 μB | **6.38e-3 μB** (~360× better) | 5.98e-3 Ry (unchanged) |
| q050 SCF | 1.29e-3 μB (already passing) | unchanged | 2.73e-4 Ry (unchanged, **now the sole failure**) |
| q050 constrained | (mirrors SCF) | unchanged | 2.73e-4 Ry (unchanged, **now the sole failure**) |

Time-reversal was checked at the same time (Anders' second question) and
ruled out as a contributor here: `force_full_bz_for_nonzero_q_gbt`
(`source/reciprocal_lifecycle.f90:526`) unconditionally forces both
`use_symmetry_reduction` and `use_time_reversal` to `.false.` for any
nonzero-q `gbt_single_q` k-space build regardless of the namelist —
confirmed by log inspection, not just code reading (see §6 for the same
check on the Γ–H sweep, where it's more directly relevant).

### 3.2 What remains — a real, smaller, still-open band-energy residual

All four `WP9CommensurateSupercell_*` cases **still fail** — now purely on
`ΔE`/atom, not moment direction:

- MFT leg: 7.19e-4 Ry (q050) / 5.98e-3 Ry (q033), against a 3.0e-6 Ry
  tolerance (that tolerance was inherited from the `3fd21c0` reference on a
  different kernel and was never re-derived for the current one — it may
  itself be too tight, but the gap is still ~6× the measured noise floor
  regardless).
- SCF leg: 2.73–2.81e-4 Ry, against 5.0e-5 Ry (~5-6× over, with the
  tolerance set at ~6× headroom over the measured 8.23e-6 Ry noise floor —
  this one is a real, reproducible gap, not a stale-tolerance artifact).

Leading hypothesis, not confirmed: the MFT leg's frozen starting potential
is carried over from the older `3fd21c0` kernel and may not be
self-consistent under the current one (the SCF leg, letting the same
starting potential relax, shows a smaller gap, consistent with this — but
an independent integrator q=0 control test using a *fresh*, non-carried-over
starting potential converged cleanly with no gap against `periodic_nc`,
which supports "stale starting potential" over "kernel defect" as the
explanation, without fully proving it). Not traced to a specific source
line. See `tests/KNOWN_ISSUES.md` for the full trace and a proposed
follow-up (regenerate the MFT leg's frozen potential under the current
kernel rather than reusing `3fd21c0`'s).

### 3.3 Constrained-SCF leg and backend note (unaffected by the fix, unchanged from the original pass)

Constrained-SCF leg: numerically identical to the SCF leg's base start
(`constrained_spiral` is already the code default). Transverse-residual
diagnostic reads `m_transverse=0.000000, torque≈0` — the single-sublattice
cell gives the constraint nothing to resist, matching WP7's own noted
limitation.

Backend note (checked per task instruction, not part of the gate):
`strux_backend='strux_lib'` badly breaks the explicit supercell's hand-built
`crystal_sym='file'` lattice — charge not conserved (7.972–8.347 e across
four translationally-equivalent sites, vs 8.0000000–8.0000010 e on
`legacy`) and a 15% site-to-site moment spread. `legacy` (used for the
`super/` decks) is correct. Recorded in `tests/KNOWN_ISSUES.md` as a
separate, `strux_lib`-specific finding, not investigated further.

**ctest:** `WP9CommensurateSupercell_*` (4 cases), labels `gbt;wp9`
(deliberately **not** `regression` — real, open, characterized gap).

## 4. Item 4 — bcc-Fe Γ–H and q-symmetry: **PASS**

`tests/regression/wp9_validation/gammaH_sweep/` (k-space route — see §6 for
why the real-space route was dropped). `magnetic_representation` forced to
`gbt_single_q` automatically by `post_processing_frozen_magnon_acoustic`.

**q↔−q:** exact `0.000e+00 Ry` over the full sweep. No SOC is present in
these runs, so this cleanly confirms no phase-convention asymmetry.

**Smoothness through H (q3=1.0):** clean, symmetric dispersion, no
discontinuity. Second differences (curvature) mirror exactly around
q3=1.0: `[(0.85,-1.18e-4), (0.90,-6.93e-4), (0.95,-1.15e-3), (1.00,-8.25e-4),
(1.05,-1.15e-3), (1.10,-6.93e-4), (1.15,-1.18e-4)]` — a textbook symmetric
zone-boundary crossing.

## 5. Item 3 — k-space EBAND convergence and RS/k-space agreement: **FAIL** (mild, reworked mid-task)

### 5.1 Why this changed from the first pass

The original design compared RS and k-space routes directly (a qualified
pass: it cleared its own RS-dominated, loose tolerance, but the raw
cross-route gap was ~7–12× WP7's non-GBT baseline). Anders instructed this
be reworked to k-space only: long-wavelength spirals can have real
numerical problems on the real-space recursion route, from finite cluster
truncation (`rc`) or periodic-boundary-condition artifacts, so validating
against a route that may itself be unreliable for this physics is the wrong
design. A delegated agent removed `base_realspace/` entirely and rewrote the
runner so the tolerance source is k-space mesh-refinement (`nk=8/12/16`)
self-consistency instead of a second route. The integrator reviewed the
diff and independently reran the resulting `ctest`.

### 5.2 Result under the new design

| q3 | k-space eband (nk=8,12,16) | spread | step(8→12) | step(12→16) | monotonic? |
| --- | --- | --- | --- | --- | --- |
| 0.5 | −1.98668388, −1.98714329, −1.98785909 | 1.175e-3 Ry | 4.594e-4 Ry | 7.158e-4 Ry | **no — fails** |
| 1.0 | −1.96445921, −1.96503119, −1.96555624 | 1.097e-3 Ry | 5.720e-4 Ry | 5.251e-4 Ry | yes |

The gate (refinement step must not grow between successive mesh
refinements) fails at q3=0.5: the step from nk=12→16 (7.16e-4 Ry) is larger
than nk=8→12 (4.59e-4 Ry). Both steps are the same order of magnitude and
well inside the overall 1.175e-3 Ry spread — not a divergence, plausibly
ordinary BZ-integration/tetrahedron discretization noise at that particular
q-point, but not diagnosed further. q3=1.0 is well-behaved.

**Time-reversal, checked directly at Anders' prompt:** confirmed by log
inspection (not just code reading) that for the nonzero-q probes in this
same deck, `force_full_bz_for_nonzero_q_gbt` correctly logs `"nonzero-q GBT
rebuilding the full chemical BZ mesh"` and builds the full, unreduced
1728-point (12³) mesh — no time-reversal, no spatial-symmetry reduction —
regardless of the deck's `use_symmetry_reduction=.true.`/
`use_time_reversal=.true.` namelist settings. Both settings are inert for
every nonzero-q build in this deck (and, separately, time-reversal is also
disabled generically for any noncollinear `nsp=3` build, GBT or not). Ruled
out as an explanation for either this finding or the cone-scaling finding
in §6. The decks were annotated with a comment recording this, left at
`.true.` to match the already-validated `wp8_littlegroup` convention.

**ctest:** `WP9GammaHSweep`, labels `gbt;wp9` (moved out of `regression` —
a real, open, mild, characterized gap).

## 6. Item 5 — 5–20° cone-angle scaling: **FAIL**

Same battery, fixed small q_ss=(0,0,0.05), k-space route, θ_ss swept
5°/10°/15°/20° — **unchanged by the k-space-only rework**, since this deck
(`base_cone/`) was already k-space-only from the first pass:

| θ_ss | ω(q=0.05) |
| --- | --- |
| 5° | −6.417e-3 Ry |
| 10° | −1.293e-3 Ry |
| 15° | −3.426e-4 Ry |
| 20° | −8.270e-6 Ry |

Not θ-independent — a ~776× spread across the window, decreasing
monotonically as θ grows, rather than the flat line the harmonic-regime
self-diagnostic (blueprint §2.10) predicts. A naive noise-amplification
argument (fixed absolute noise in ΔE, divided by `sin²θ`) predicts only a
~15× range across 5°→20°, an order of magnitude short of what's measured.
Time-reversal/BZ-reduction is ruled out as the cause (§5.2 — the same
per-build guard applies to every nonzero-q build in this run). Whether the
noise itself grows at small θ, or something else is at play, is unknown —
not diagnosed further, out of scope for this task. **"Small-cone scaling is
within the declared budget" is not met as measured.**

## 7. Item 6 — LKAG stiffness comparison: **DEFERRED**

Per explicit instruction (Anders: "we can skip the LKAG comparison for now.
I will do that later."), not attempted. Not scored as pass or fail.

## 8. Item 7 — Multi-sublattice Goldstone: measured, unshifted, **encouraging but incomplete**

`tests/regression/wp9_validation/multisublattice_goldstone/` — bcc FeCo
(`example/bulk/bccFeCo`), `nsp=3`, `branch_mode='auto'`, `mode='mft'`,
`theta_probe=20°`, **real-space route**. Reference SCF: two inequivalent
converged moments, 1.976 μB (Fe) / 1.657 μB (Co).

Measured `ω(Γ)`, acoustic branch: **`-2.15e-22 Ry`** (a second, independent
re-run through the registered ctest gave `-2.04e-15 Ry` — both zero to
numerical noise). This is **not** the ~0.28 Ry violation
`tests/KNOWN_ISSUES.md` records from the pre-WP1 architecture. Supporting
evidence this is genuine: the acoustic-branch eigenvector at Γ has both
sublattices in phase (a uniform rotation, the Goldstone-mode signature),
while the optical branch (finite gap, `ω(Γ)=7.04e-3 Ry`) has them out of
phase; small-q growth is smooth and roughly quadratic, not flat. No shift,
renormalization, or correction was applied to reach this number.

**Important caveat:** this measurement is on the **real-space** route. The
original ~0.28 Ry violation was measured on the **k-space** route — that
route was **not** re-tested here, so this is one new, positive,
real-space-only data point, not a closed finding.

**ctest:** `WP9MultisublatticeGoldstone`, labels `regression;gbt;wp9` — a
measurement with a loose sanity bound, not a Goldstone-satisfaction gate.

## 9. Regression evidence

| suite | before WP9 | after WP9 (incl. fixture fixes) |
| --- | --- | --- |
| `ctest -L unit` | 22/22 | 22/22 (unchanged) |
| `ctest -L "quick\|regression"` | 26/27 | 26/27 (identical single pre-existing failure, `Example_bulk_bccFe_nsp2_block_hoh`, the known gfortran-13 DOS-tolerance delta recorded in `docs/dev/GBT_WP0_G0_REPORT.md`) |

No production source files were changed in this task (verified via `git
status`/`git diff`): only `tests/regression/wp9_validation/` (new and
subsequently fixed), `tests/KNOWN_ISSUES.md`, and `CMakeLists.txt`'s
test-registration section changed. The full `regression`/`triad`/`backend`
multi-hour ctest suite was not re-run to completion (same rationale as WP8);
the fast `unit` and `quick`+`regression` (backend matrix) subsets were
re-run clean before and after, including after the mid-task fixture fixes.

## 10. Completion checklist

- [x] Analytic and dense-oracle tests pass.
- [ ] Commensurate MFT and SCF supercell tests pass — **FAIL**. The
      original near-zero-moment failure was a test-fixture bug (`mom`
      convention), traced and fixed (§3.1); a smaller, real band-energy
      residual remains and is still open (§3.2).
- [ ] RS/k-space observables and EBAND agree when converged — **FAIL**
      under the reworked k-space-only design (§5): mild mesh-refinement
      non-monotonicity at q3=0.5. (Under the original mixed-route design
      this was a qualified pass; the design itself was judged the wrong
      validation and replaced.)
- [x] bcc-Fe Gamma-H and q-symmetry tests pass — q↔−q and smoothness both
      clean (§4).
- [ ] Small-cone scaling is within the declared budget — **FAIL**, 776×
      spread across 5°–20°, time-reversal ruled out as the cause (§6).
- [ ] LKAG comparison includes shell convergence/error bars — **deferred**
      per Anders, not attempted (§7).
- [~] Multi-sublattice Goldstone test is unshifted and explained — measured
      gapless and unshifted on the real-space route; k-space route (where
      the original violation was measured) not re-checked (§8).
- [x] Inputs, raw data, convergence, uncertainty, and failures are
      archived — `tests/regression/wp9_validation/` (6 ctests),
      `tests/KNOWN_ISSUES.md` (5 entries), this report.
- [x] G9 PASS/FAIL is stated: **G9: FAIL** — one major finding was
      diagnosed to a test-fixture root cause and fixed, materially
      shrinking (but not closing) that item's gap; two items still fail on
      real, open, well-characterized findings; one is an incomplete
      positive result; LKAG is deferred. The algebraic/oracle layer and
      Γ–H/q-symmetry layer are solid. Nothing here was tuned, shifted, or
      hidden to pass; every finding is reproducible via its registered
      `ctest`.

## 11. Files changed

| file | change |
| --- | --- |
| `tests/regression/wp9_validation/commensurate_supercell/` | new, then fixed — decks (ported/adapted from `fable_v2_gbt`'s `3fd21c0`; `mom` convention corrected on all 6 `gbt_single_q` decks), `run_wp9_commensurate_supercell.py`, `cases.json` |
| `tests/regression/wp9_validation/gammaH_sweep/` | new, then reworked to k-space-only — `base_realspace/` deleted, `run_wp9_gammaH_sweep.py` rewritten, `use_time_reversal` annotated as inert |
| `tests/regression/wp9_validation/multisublattice_goldstone/` | new — deck (seeded from a converged potential), `run_wp9_goldstone.py` |
| `CMakeLists.txt` | 6 `ctest` registrations; `WP9CommensurateSupercell_*` and `WP9GammaHSweep` under `gbt;wp9` (not `regression` — real open gaps); `WP9MultisublatticeGoldstone` under `regression;gbt;wp9` |
| `tests/KNOWN_ISSUES.md` | 5 entries: commensurate-supercell `mom`-convention fixture bug (resolved) with the remaining band-energy residual as open follow-up; `strux_lib` supercell-backend finding; cone-angle scaling finding; k-space mesh non-monotonicity finding; multi-sublattice Goldstone entry updated with the fresh real-space measurement |

## 12. Unresolved risks / next task

- The commensurate-supercell band-energy residual (§3.2) is the main open
  physics question from this battery: 7e-4–6e-3 Ry/atom (MFT) and
  2.7-2.8e-4 Ry/atom (SCF), larger for the 6-atom q033 case. Leading
  hypothesis (stale starting potential) is plausible but unconfirmed.
  Suggested next step: regenerate the MFT leg's frozen starting potential
  under the *current* kernel instead of reusing `3fd21c0`'s, and re-measure.
- Cone-angle scaling's 776× spread (§6) needs either a dedicated
  DOS/mesh-convergence study at fixed small θ, or a documented reason it
  isn't expected to be constant in this regime, before the harmonic-regime
  self-diagnostic can be trusted for production frozen-magnon θ-sweeps.
  Time-reversal is ruled out; the mechanism is otherwise unknown.
- The k-space mesh-refinement non-monotonicity at q3=0.5 (§5.2) is mild and
  plausibly ordinary discretization noise, but unconfirmed.
- The multi-sublattice Goldstone result (§8) should be re-measured on the
  **k-space** route before treating `tests/KNOWN_ISSUES.md`'s original
  entry as resolved.
- `strux_backend='strux_lib'` breaking a hand-built `crystal_sym='file'`
  lattice (§3.3) is a finding independent of GBT — worth its own triage
  task.
- `WP9CommensurateSupercell_*` and `WP9GammaHSweep` are deliberately
  excluded from the `regression` label (labels `gbt;wp9` only) for their
  respective open gaps. Add `regression` back once each closes.

Next allowed task per the dependency graph: **WP10** (legacy deletion and
documentation) is gated on G9 passing — it does not, so WP10 should not
start until the remaining findings above are triaged (or Anders explicitly
overrides the gate). The LKAG comparison Anders deferred should also land
before any G9 re-evaluation.
