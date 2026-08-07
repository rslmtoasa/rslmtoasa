# GBT WP9 / Gate G9 report — physical validation

Date: 2026-08-07. Branch `fable_v2_gbt_v2`, on top of `936ea5a` (final G8:
PASS, `docs/dev/GBT_WP8_G8_REPORT.md`).

## Outcome

**Gate G9: FAIL**, with important qualifications. Physics was not changed
during this task — every finding below is a measurement, not a fix, per the
task's explicit instruction ("do not change physics during the first run...
do not regenerate goldens or impose a post-hoc Goldstone shift to pass").

- The algebraic/oracle layer (analytic Stoner, dense bond oracle),
  q-symmetry (q↔−q, smoothness through the zone boundary), and — the most
  encouraging result of this task — the multi-sublattice Goldstone check all
  measure clean.
- Two ordered validation items **fail**, with well-characterized, narrow
  causes rather than vague symptoms: the commensurate-supercell known-answer
  test, and cone-angle (θ) scaling of ω(q).
- The RS-vs-k-space EBAND agreement item is a **qualified pass**: it clears
  its own derived tolerance, but that tolerance is loose because the
  real-space route is not tightly converged in the range tested, and the
  production-resolution cross-route difference is measurably larger than
  the established non-GBT baseline (WP7).
- The LKAG stiffness comparison is **deferred** per Anders' explicit
  instruction, not attempted here.

Per the task's own ground rules, none of the failures below were fixed,
tuned around, or hidden. Each is recorded in `tests/KNOWN_ISSUES.md` with
what was verified versus guessed, and registered as a `ctest` so it stays
visible rather than silently regenerated away.

## 1. Execution model

WP9 was run as three independent, parallel validation batteries (matching
the task prompt's suggested delegation: "up to three ... runners; one ...
diagnosis owner") plus integrator verification:

- **Battery A** — commensurate-supercell known-answer test (MFT + SCF +
  constrained-SCF legs), `tests/regression/wp9_validation/commensurate_supercell/`.
- **Battery B** — bcc-Fe Γ–H frozen-magnon sweep (RS-vs-k-space EBAND,
  q↔−q, smoothness, cone scaling), `tests/regression/wp9_validation/gammaH_sweep/`.
- **Battery C** — multi-sublattice Goldstone measurement (bcc FeCo),
  `tests/regression/wp9_validation/multisublattice_goldstone/`.

A mid-task editor/session restart killed Battery B's and Battery C's
in-flight numerical runs (their setup work — decks, scripts — survived on
disk; only the actual `rslmto.x` runs were lost). Both were re-run directly
by the integrator rather than re-dispatching fresh subagents, since the
fixtures were already built and reviewed; Battery C's re-run was seeded
from an already-converged reference potential (see §4) to keep the restart
cheap and avoid the mixing instability a cold start showed. All numbers
reported below were either produced by the original battery runs (Battery
A, uninterrupted) or independently re-run and read directly by the
integrator (Battery B, Battery C) — nothing here is taken on a subagent's
word alone.

## 2. Item 1 — Analytic Stoner and dense bond oracle: **PASS**

No new work needed: `tests/unit/test_gbt_oracles.py` (`ctest`
`UnitGbtOracles`) already carries the one-orbital Stoner k∓q/2 identity and
the dense unequal-species pair oracle, established under WP1/G1 and
re-verified at every subsequent gate. Re-run on the current HEAD used for
this task: `ctest -L unit` **22/22**, including `UnitGbtOracles` and every
other GBT unit oracle (`UnitGbtWp5SourceContract`, `UnitGbtWp6Hoh`,
`UnitGbtWp6Ccor`, `UnitGbtWp6cVelocity`, `UnitGbtWp6Matrix`,
`UnitMagneticRepresentation`, `UnitGbtStructure`, `UnitSpinDensity`). No
regression.

## 3. Item 2 — Commensurate supercell, MFT then SCF: **FAIL**

`tests/regression/wp9_validation/commensurate_supercell/` — bcc Fe, θ=90°
flat spiral, q_ss=(0,0,0.5) → 4-atom supercell and q_ss=(0,0,1/3) → 6-atom
supercell (§3.3 sizing, integer winding number, not the naive 2-/3-atom
cells). Ported and adapted from a fixture on the sibling branch
`fable_v2_gbt` (commit `3fd21c0`, never merged here — reused after review,
per the delegation model's "reuse only after review" rule) built for an
earlier, single-kernel GBT architecture; adapting it to the current
representation-mode architecture (`magnetic_representation='gbt_single_q'`,
`strux_backend='strux_lib'`, both WP4 requirements not present in the
original) is itself part of what this battery measured.

**MFT leg (force theorem, `nstep=1`, frozen potential):** fails by a large,
reproducible margin — not noise:

| case | Δq | Δmoment | ΔE/atom |
| --- | --- | --- | --- |
| q050, base start | 9.8e-7 e | **2.396 μB** | **7.19e-4 Ry** |
| q050, perturbed start | 7.5e-7 e | **2.325 μB** | **1.43e-4 Ry** |

against derived noise floors of Δq≈1.8e-6 e, Δmoment≈1.1e-6 μB (intra-run
symmetry spread) — the moment gap is six orders of magnitude above noise.
**Root cause, narrowed but not fully diagnosed:** running the registered
ctest with `--output-on-failure` surfaces `source/bands.f90`'s
`DENSITY_POLICY` diagnostic for the `gbt_single_q` side:
`m_long=0.000000 |m_transverse|=2.395328
torque=(0.000000,-2.395328,0.000000)`. The moment is not actually missing —
its magnitude (2.395 μB) matches the supercell side almost exactly — it is
sitting entirely in the **transverse** channel instead of the longitudinal
(cone-axis) one the `ql`-based up/down extraction reads. This is a
frame/axis mismatch specific to the `nstep=1` force-theorem evaluation of a
potential carried over from a different kernel, not a magnitude-destroying
bug. A follow-up integrator control test (fresh SCF start, not the
carried-over potential, `q_ss=0, theta_ss=0`) showed `gbt_single_q` and
`periodic_nc` converge to matching moments (2.131 vs 2.128 μB, charge exact
to 7e-15 e) — so the core q=0 identity contract (§2.4 of the blueprint) is
**not** violated in general; the failure is scoped to the MFT leg's
frozen-potential/foreign-starting-point evaluation. Full detail, including
what was verified directly versus inferred, is in `tests/KNOWN_ISSUES.md`.

**SCF leg** (full self-consistency from the same starting potential,
`nstep=25`, broyden `beta=0.15`): much closer, but still a real,
above-noise gap:

| quantity | base start | perturbed start | derived tolerance |
| --- | --- | --- | --- |
| Δcharge | 2.04e-3 e | 5.69e-5 e | 1.0e-2 e (pass) |
| Δmoment | 1.29e-3 μB | 4.34e-4 μB | 5.0e-3 μB (pass) |
| ΔE/atom | **2.73e-4 Ry** | **2.81e-4 Ry** | 5.0e-5 Ry (**fail**, ~5-6×) |

Charge and moment pass; the band-energy leg fails a tolerance set with ~6×
headroom over its own measured noise floor (8.23e-6 Ry) — a real,
reproducible gap, not noise.

**Constrained-SCF leg:** numerically identical to the SCF leg's base start
(`constrained_spiral` is already the code default), confirming that. The
transverse-residual/torque diagnostic reads `m_transverse=0.000000,
torque≈0` at convergence — the single-sublattice cell gives the constraint
nothing to resist, matching WP7's own noted limitation that policy
distinctness needs a genuinely multi-sublattice constrained state.

**Backend note (checked per task instruction, not part of the gate):**
`strux_backend='strux_lib'` badly breaks the explicit supercell's
hand-built `crystal_sym='file'` lattice — charge not conserved
(7.972–8.347 e across four translationally-equivalent sites, vs
8.0000000–8.0000010 e on `legacy`) and a 15% site-to-site moment spread.
`legacy` (the default, used for the `super/` decks) is unaffected and
correct. This is recorded in `tests/KNOWN_ISSUES.md` as a separate,
`strux_lib`-specific finding, not investigated further (out of scope).

**ctest:** `WP9CommensurateSupercell_*` (4 cases), labels `gbt;wp9`
(deliberately **not** `regression` — see §7).

## 4. Item 3 — RS/k-space observables and EBAND agreement: qualified **PASS**

`tests/regression/wp9_validation/gammaH_sweep/` — bcc Fe,
`magnetic_representation` forced to `gbt_single_q` automatically by
`post_processing_frozen_magnon_acoustic` (confirmed by reading
`source/calculation.f90` directly), `strux_backend='strux_lib'`, `nsp=3`.
This is the first fixture in the repository to compare RS and k-space with
GBT actually active (nonzero `q_ss`) — WP7's own RS/k-space comparison used
a plain collinear ferromagnet.

Convergence spread (own-route, at q3=0.5 and q3=1.0/H):

| q3 | k-space (nk=8/12/16) spread | RS (lld=15/21/27) spread |
| --- | --- | --- |
| 0.5 | 1.175e-3 Ry | **1.149e-2 Ry** |
| 1.0 | 1.097e-3 Ry | **9.864e-3 Ry** |

The RS spread is dominated by `lld=15` being poorly converged (the
15→21 jump is ~0.011 Ry; 21→27 is ~0.0007 Ry) — `lld≥21` is needed for a
fair comparison here, consistent with but somewhat deeper than the WP7
precedent (`lld=15–27` was WP7's own range, without flagging `lld=15` as
insufficient — worth noting for future GBT convergence work).

Cross-route agreement at production resolution (k-space nk=12, RS lld=21),
checked against the (loose, RS-dominated) tolerance above:

| q3 | k-space eband | RS eband | diff | tol | result |
| --- | --- | --- | --- | --- | --- |
| 0.5 | −1.9871433 | −1.9861726 | 9.71e-4 Ry | 1.149e-2 Ry | pass |
| 1.0 | −1.9650312 | −1.9633936 | 1.64e-3 Ry | 9.864e-3 Ry | pass |

**Qualification:** both checks pass, but only against a tolerance inflated
by RS's own under-convergence at `lld=15`. The raw cross-route difference
itself (~1–1.6e-3 Ry) is roughly 7–12× larger than WP7's collinear (non-GBT)
RS/k-space etot agreement (1.33e-4 Ry, `docs/dev/GBT_WP7_G7_REPORT.md`).
This is not asserted as a failure — no non-GBT-comparable tolerance exists
to fail against — but it is a real, measured degradation worth flagging
for anyone converging a production GBT spiral: check `lld` convergence
explicitly rather than reusing WP7's non-GBT settings uncritically.

**ctest:** `WP9GammaHSweep`, labels `regression;gbt;wp9` (its own gated
checks pass cleanly).

## 5. Item 4 — bcc-Fe Γ–H and q-symmetry: **PASS**

Same deck/battery as §4.

**q↔−q:** k-space exact `0.000e+00 Ry` (over the full sweep); real-space
`5.89e-6 Ry` (at q3=±0.6) — both far inside the noise floor established in
§4. No SOC is present in these runs, so this cleanly confirms no
phase-convention asymmetry.

**Smoothness through H (q3=1.0):** both routes show a clean, symmetric
dispersion with no discontinuity. Second differences (curvature) mirror
exactly around q3=1.0 on both routes, e.g. real-space:
`[(0.85,-4.64e-4), (0.90,-6.68e-4), (0.95,-9.27e-4), (1.00,-1.05e-3),
(1.05,-9.27e-4), (1.10,-6.68e-4), (1.15,-4.64e-4)]` — a textbook symmetric
zone-boundary crossing.

## 6. Item 5 — 5–20° cone-angle scaling: **FAIL**

Same deck/battery, fixed small q_ss=(0,0,0.05), k-space route, θ_ss swept
5°/10°/15°/20°:

| θ_ss | ω(q=0.05) |
| --- | --- |
| 5° | −6.417e-3 Ry |
| 10° | −1.293e-3 Ry |
| 15° | −3.426e-4 Ry |
| 20° | −8.270e-6 Ry |

This is **not** θ-independent — a ~776× spread across the window, decreasing
monotonically as θ grows, rather than the flat line the harmonic-regime
self-diagnostic (blueprint §2.10, "ω(q) must come out independent of θ")
predicts. A naive noise-amplification argument (fixed absolute noise in
ΔE, divided by `sin²θ`) predicts only a ~15× range across 5°→20°, an order
of magnitude short of what's measured — so simple Fermi-search/DOS additive
noise does not, by itself, explain the full spread; whether the noise
itself grows at small θ (plausible: a smaller induced-moment signal is
proportionally noisier against the same absolute DOS/mesh discretization)
or something else is at play was not diagnosed further (out of this task's
scope — this is exactly the kind of finding to record and hand back, not
chase to a fix). The battery's runner script does not gate on this (it is
reported informationally per its own design, since no tolerance had been
derived for it going in); this report is where it is surfaced as a real,
open, checklist-relevant finding. **"Small-cone scaling is within the
declared budget" is not met as measured** — no budget was ever declared
that this satisfies, and 776× is not a defensible convergence spread by any
reasonable budget.

## 7. Item 6 — LKAG stiffness comparison: **DEFERRED**

Per explicit instruction (Anders: "we can skip the LKAG comparison for now.
I will do that later."), not attempted. Not scored as pass or fail.

## 8. Item 7 — Multi-sublattice Goldstone: measured, unshifted, **encouraging but incomplete**

`tests/regression/wp9_validation/multisublattice_goldstone/` — bcc FeCo
(`example/bulk/bccFeCo`), `nsp=3`, `branch_mode='auto'`, `mode='mft'`,
`theta_probe=20°`, **real-space route** (`recur='block'`, `lld=21`,
`strux_backend='strux_lib'`). Reference SCF: two inequivalent converged
moments, 1.976 μB (Fe) / 1.657 μB (Co) (seeded from an already-converged
starting potential after a cold start showed a real, if eventually
self-correcting, Broyden oscillation — see the script docstring).

Measured `ω(Γ)`, acoustic branch: **`-2.15e-22 Ry`** (a second, independent
re-run through the registered ctest gave `-2.04e-15 Ry` — both zero to
numerical noise by any reasonable standard, small run-to-run differences
reflecting the SCF reference not being driven to machine-precision
convergence). This is **not** the ~0.28 Ry violation
`tests/KNOWN_ISSUES.md` records from the pre-WP1 architecture. Supporting
evidence this is a genuine physical result, not a fluke: the acoustic-branch
eigenvector at Γ has both sublattices in phase (amplitudes 0.7375/0.6753,
both at phase −π — a uniform rotation, exactly the Goldstone-mode
signature), while the optical branch (finite gap, `ω(Γ)=7.04e-3 Ry`) has
them out of phase; small-q growth (`ω` at q3=0.02/0.05: 1.18e-5/7.87e-5 Ry)
is smooth and roughly quadratic, not flat (the old, broken measurement was
"nearly flat"). No shift, renormalization, or correction was applied to
reach this number — it is what the code reports.

**Important caveat, stated plainly per the task's honesty requirement:**
this measurement is on the **real-space** route. The original ~0.28 Ry
violation was measured on the **k-space** route
(`tests/KNOWN_ISSUES.md`'s original entry: "it appears on the k-space
path — the only one exercised so far for multi-sublattice"). This task did
**not** re-test the k-space route for the multi-sublattice case, so this
result does **not** by itself confirm the historical finding is resolved —
it is one new, positive, real-space-only data point on one system with one
θ_probe. `tests/KNOWN_ISSUES.md`'s entry is updated to record this
precisely (re-measured clean on RS, k-space not re-checked) rather than
marked resolved.

**ctest:** `WP9MultisublatticeGoldstone`, labels `regression;gbt;wp9` — a
measurement with a loose sanity bound (catches NaN/Inf or an unphysical
energy scale), not a Goldstone-satisfaction gate, per the task's explicit
"do not impose a post-hoc Goldstone shift to pass" instruction.

## 9. Regression evidence

| suite | before WP9 | after WP9 |
| --- | --- | --- |
| `ctest -L unit` | 22/22 | 22/22 (unchanged) |
| `ctest -L "quick\|regression"` | 26/27 | 26/27 (identical single pre-existing failure, `Example_bulk_bccFe_nsp2_block_hoh`, the known gfortran-13 DOS-tolerance delta recorded in `docs/dev/GBT_WP0_G0_REPORT.md`) |

No source files were changed in this task (verified: only
`tests/regression/wp9_validation/` (new), `tests/KNOWN_ISSUES.md`, and
`CMakeLists.txt`'s test-registration section changed). The full
`regression`/`triad`/`backend` multi-hour ctest suite was not re-run to
completion (same rationale as WP8: this task's diff cannot plausibly touch
real-space recursion/exchange machinery it doesn't call, and a full re-run
is outside a single task's "no marathons" budget); the fast `unit` and
`quick`+`regression` (backend matrix) subsets were re-run clean before and
after.

## 10. Completion checklist

- [x] Analytic and dense-oracle tests pass.
- [ ] Commensurate MFT and SCF supercell tests pass — **FAIL**, both legs,
      root cause narrowed to a frame/axis mismatch in the frozen-potential
      MFT evaluation (§3); recorded in `tests/KNOWN_ISSUES.md`.
- [~] RS/k-space observables and EBAND agree when converged — qualified
      pass; RS convergence at low `lld` is looser than assumed, and the
      cross-route gap is real though small (§4).
- [x] bcc-Fe Gamma-H and q-symmetry tests pass — q↔−q and smoothness both
      clean (§5).
- [ ] Small-cone scaling is within the declared budget — **FAIL**, 776×
      spread across 5°–20° (§6).
- [ ] LKAG comparison includes shell convergence/error bars — **deferred**
      per Anders, not attempted (§7).
- [~] Multi-sublattice Goldstone test is unshifted and explained — measured
      gapless and unshifted on the real-space route; k-space route (where
      the original violation was measured) not re-checked, so this is
      encouraging new evidence, not a closed finding (§8).
- [x] Inputs, raw data, convergence, uncertainty, and failures are
      archived — `tests/regression/wp9_validation/` (6 new ctests),
      `tests/KNOWN_ISSUES.md` (3 entries updated/added), this report.
- [x] G9 PASS/FAIL is stated: **G9: FAIL** — two ordered items fail with
      well-characterized causes, one is a qualified pass, one is
      encouraging-but-incomplete, and LKAG is deferred by instruction. The
      algebraic/oracle layer and q-symmetry layer are solid. Nothing here
      was tuned, shifted, or hidden to pass; every finding is reproducible
      via its registered `ctest`.

## 11. Files changed

| file | change |
| --- | --- |
| `tests/regression/wp9_validation/commensurate_supercell/` | new — decks (ported/adapted from `fable_v2_gbt`'s `3fd21c0`), `run_wp9_commensurate_supercell.py`, `cases.json` |
| `tests/regression/wp9_validation/gammaH_sweep/` | new — decks, `run_wp9_gammaH_sweep.py`, q-point files |
| `tests/regression/wp9_validation/multisublattice_goldstone/` | new — deck (seeded from a converged potential), `run_wp9_goldstone.py` |
| `CMakeLists.txt` | 6 new `ctest` registrations (`WP9CommensurateSupercell_*` ×4 under `gbt;wp9`; `WP9GammaHSweep`, `WP9MultisublatticeGoldstone` under `regression;gbt;wp9`) |
| `tests/KNOWN_ISSUES.md` | 2 new entries (frozen-potential/frame-mismatch finding; `strux_lib` supercell-backend finding) + 2 existing entries updated with fresh WP9 measurements (frozen-potential q=0 control test; multi-sublattice Goldstone real-space re-measurement) |

## 12. Unresolved risks / next task

- The commensurate-supercell MFT leg's frame/axis mismatch (§3) has a
  concrete next clue (`DENSITY_POLICY`'s `m_long`/`m_transverse` split) but
  is not diagnosed to a source line. Worth a dedicated follow-up before
  trusting any `gbt_single_q` force-theorem (`nstep=1`) evaluation seeded
  from a potential converged under a different code path.
- The SCF leg's ~2.7e-4 Ry/atom band-energy gap against the supercell
  (§3) — smaller than the MFT leg's, but still 5-6× over a
  noise-derived tolerance — is unexplained.
- Cone-angle scaling's 776× spread (§6) needs either a dedicated
  DOS/mesh-convergence study at fixed small θ, or a documented reason it
  isn't expected to be constant in this regime, before the harmonic-regime
  self-diagnostic can be trusted for production frozen-magnon θ-sweeps.
- The multi-sublattice Goldstone result (§8) should be re-measured on the
  **k-space** route before treating `tests/KNOWN_ISSUES.md`'s original
  entry as resolved — that is where the ~0.28 Ry violation was originally
  found, and it was not re-tested here.
- `strux_backend='strux_lib'` breaking a hand-built `crystal_sym='file'`
  lattice (§3, backend note) is a finding independent of GBT — it affects
  any `periodic_nc`/ordinary calculation using that combination, not just
  this battery's supercell decks. Worth its own triage task.
- The `WP9CommensurateSupercell_*` ctests are deliberately excluded from
  the `regression` label (labels `gbt;wp9` only) so the default
  `-L regression` run stays green for a known, documented, open gap
  instead of turning permanently red. Whoever picks up the frame-mismatch
  fix should add `regression` back to those labels once the gap closes.

Next allowed task per the dependency graph: **WP10** (legacy deletion and
documentation) is gated on G9 passing — it does not, so WP10 should not
start until the failures above are triaged (or Anders explicitly overrides
the gate). The LKAG comparison Anders deferred should also land before any
G9 re-evaluation.
