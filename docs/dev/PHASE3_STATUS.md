# Phase 3 status — feature sweep + documentation audit

**Scope.** Everything landed on `fable_v2` between commit `a8ec82e96`
("Fable v2 main merge", 2026-07-05) and `HEAD` (2026-07-28) — 59 commits.
Covers the blueprint items `B1`–`B7` from `BLUEPRINTS.md` /
`claude_files/plans/00_INDEX.md`, plus a doc-location audit and a proposed
homogenization of the three overlapping plan-file trees. `B8`–`B12` are not
started and out of scope here.

**How this was produced.** Full `git log a8ec82e96..HEAD` read commit by
commit; every named source file (`reciprocal_green.f90`, `reciprocal_bsf.f90`,
`electrostatics_multipole.f90`, `region_registry.f90`, `vacuum_lead.f90`,
`interfacepot`/`charge.f90`, etc.) read in its current state, not just diffed;
`tests/KNOWN_ISSUES.md`, `docs/dev/*.md`, and `claude_files/plans/*.md` cross
checked against the code; `docs/source/**/*.rst` checked against the actual
namelist keys and module list. B3's "is it really implemented" question was
spot-verified directly (not just via a subagent report) by reading
`calculation.f90:1488-1513` and `reciprocal_bsf.f90` and confirming it calls
the real `dyson_kspace_inverse` kernel, not a stub. B4 was spot-verified by
grepping the whole tree for `rsksp_gpu` (the planned plugin name) — zero
hits outside the unpromoted blueprint file.

---

## 1. Feature status by blueprint item

| # | Feature | Status | Test/example coverage |
|---|---|---|---|
| — | Pre-B1 groundwork (`3073adb`…`d86fe42`) | tetrahedron DOS integration fix, test-tolerance bump, first (broken) GBT/frozen_magnon attempt | superseded by B1; no separate coverage needed |
| **B1** | GBT spin-spiral fix + `frozen_magnon` post-processing | **done** | good |
| **B2** | `reciprocal_green` k-space GF (Lehmann + Dyson backends) | **done**, gates G-B2-1/G-B2-2 signed | orphaned (exists, not wired in) |
| **B3** | Bloch spectral functions (`post_processing='bsf'`) | **done** (verified directly, see above) | **none** |
| **B4** | GPU k-space (batched eigensolver plugin) | **not started** | n/a |
| **B5** | Route-agnostic Jij/damping/conductivity | **done** | orphaned (exists, not wired in) |
| **B6** | Surface/interface dipole electrostatics | **done**, narrower scope than blueprint | **none** |
| **B7** | Boundary self-energies (interface leads + vacuum) | **done**, bulk of recent work | best-covered |

### B1 — GBT spin-spiral fix + frozen magnons
Fixed `ham0m_nc`'s spiral branch: bond-vector phase (not absolute site
position) and half-angle SU(2) rotation, extended through HOH/CCOR
(`aabb842`, `130b20d`). Added `frozen_magnon` post-processing: E(q) sweep
and MFT-style adiabatic magnon dispersion, with the reference potential held
fixed across the q-sweep (`ce456e2`).

- **Coverage:** `tests/scf/cases.json` → `Example_bulk_bccFe_nsp4_block_spiral_q{plus,minus}`
  (E(q)=E(−q) regression) and `frozen_magnon/bccFe*` (E(q) sweep/dispersion).
- **Known-incomplete:** `branch_mode='auto'` multi-sublattice path is not
  gapless at Γ (~0.28 Ry Goldstone violation on two-sublattice bcc FeCo);
  root cause not found, deferred to B11. Documented in `tests/KNOWN_ISSUES.md`
  and `docs/DEVELOPER_MAP.md`.

### B2 — `reciprocal_green` (k-space Green's functions, flagship)
New submodule `source/reciprocal_green.f90` of `reciprocal_mod`, two
backends: **E** (strict Lehmann, Σ=0, one eigensolve per k) and **D** (direct
Dyson `[zS−H(k)−Σ]⁻¹`, batched `getrf/getrs`), dispatched via a new
`gf_route` key (`'recursion'|'lehmann'|'dyson'`). Fills the same
`green%gij/gji`, `gij_eta` ladder, and torque arrays the recursion route
fills — zero consumer changes. Supporting files: `sigma_provider.f90`
(abstract Σ interface), `lehmann_kernel.f90`, `dyson_kernel.f90`.

- **Gates:** G-B2-1 (contour convention) and G-B2-2 (`green_eta=0.02` Ry
  default, `nk=16³`) both signed. Convergence study in
  `docs/dev/reciprocal_green_convergence.md`.
- **Coverage:** validation exists only as a standalone, well-commented
  example — `example/exchange/bccFe_kspace_green/input.nml` — **not**
  registered in any `cases.json`, not documented in any rst page.

### B3 — Bloch spectral functions
`source/reciprocal_bsf.f90` + `bsf_kernel.f90`: `A(k,E) = −(1/π) Im Tr
G(k,E+iη)` along the canonical high-symmetry path, using the same
`dyson_kspace_inverse` primitive as B2 backend D. Wired into
`calculation.f90` (`post_processing_bsf`, dispatch case at line 276) —
confirmed real, not a stub.

- **Coverage:** **zero.** No `cases.json` entry, no `example/` directory, no
  rst page, no unit test beyond the low-level `tests/unit/test_bsf_sumrule.f90`
  (which pins the kernel's sum rule, not an end-to-end example).
- This is the easiest feature to mistake for "not started" — it isn't, it's
  just invisible to anyone not reading `calculation.f90` directly.

### B4 — GPU port of k-space functionality
**Not started.** No `rsksp_gpu` plugin, no batched eigensolver/`getrf`/GEMM
primitives anywhere in `source/`. `claude_files/plans/B4_gpu_kspace.md` is
the only artifact, and it is the original unpromoted blueprint text, not a
progress record.

### B5 — Route-agnostic G(E) post-processing
Three parts, all landed: (1) `reciprocal_moments.f90` — exact k-space
Chebyshev moment generator feeding the existing KPM conductivity pipeline
unchanged; (2) J_ij/σ triads pinned as regression comparing
recursion/lehmann/dyson on bcc Fe with documented tolerance envelopes
(`docs/dev/route_agnostic_estimators.md`); (3) Gilbert damping — audited
first (`docs/dev/B5.3_gilbert_damping_audit.md`: SOC-derivative torque
operator was already correct, no new physics needed), then wired via
`calculation%do_damping`. Along the way, fixed the real root cause of the
exchange-module NaN: a `simpson_f` one-element out-of-bounds read in
`math.f90` (`8b42928`) — unrelated to damping but blocking it.

- **Coverage:** `tests/regression/triad_bccFe_damping/` + `run_triad.py` is a
  solid, purpose-built harness with goldens — but it is a standalone script,
  **not** registered in `tests/regression/cases.json` or mentioned in
  `tests/README.md`.

### B6 — Surface/interface dipole electrostatics
Scope narrowed from the original blueprint's full 3D/2D Ewald l≤2
generalization to reusing the **existing** `charge%madl2d` dipole-monopole
(`dsz`) / dipole-dipole (`dzz`) Madelung matrices (already present, just
previously unfed) plus a new charge-side moment. New module
`source/electrostatics_multipole.f90` (`compute_dipole_moments`): l=1 (z)
dipole moment Q₁₀ per atom from on-site cross-orbital density-matrix
elements (s-pz, pz-dz2) and a new exported radial partial-wave amplitude
(`potential%q10`, `phi_amp` export hook). Gated by
`control%dipole_electrostatics`.

- **Validation performed:** bulk limit / feature-disabled reduces
  bit-identically to the monopole result (Q10→0).
- **Coverage:** **zero.** No test case, no example, no rst doc. Gate G-B6-1
  (which Skriver–Rosengaard reference tables define acceptance) is
  unsigned/not reached — no literature work-function comparison was done.
  The only trace of this work outside `source/` is the untracked
  `docs/dev/electrostatic_dipole_plan.md` (see §2 below).

### B7 — Boundary self-energies (interface leads + vacuum)
The bulk of the last three weeks' work, B7.0 through B7.6: Madelung
convention audit + real bug fix (`exh=-exh` stale value in `madl2d`,
`6e87860`); `source/region_registry.f90` (explicit per-site region
bookkeeping, replacing magic-offset arithmetic); `source/vacuum_lead.f90`
(frozen empty-sphere potential parameters, validated against an analytic
spherical-Bessel oracle); `interfacepot`/`align_regions` in `charge.f90`
(two-sided deviation-variable electrostatics reusing `surfmat`'s Madelung
matrices); an alignment solver (fixed-point V_r, gauge-anchored, gate
G-B7-2 signed); `calctype='L'` wired end-to-end (fixing two real bugs found
in the process — a missing `case('L')` in 8 select-cases that silently
skipped the Hamiltonian build, and a wrong representative-site pick); and
finally wiring `vacuum_lead` as a real consumer via `&lattice
region_b_kind='vacuum'` (fixing three index-space bugs that made interface
electrostatics identically zero).

- **Coverage:** best in the whole sweep. Cross-calctype fcc Cu oracle
  (`Example_bulk_fccCu_chebyshev` / `Example_impurity_fccCu_chebyshev` /
  `Example_interface_fccCu{001,111}_chebyshev`, same converged parameter set
  and pinned Fermi level so the four routes are directly comparable),
  documented at length in `tests/scf/README.md`. Dedicated rst page:
  `docs/source/user_guide/examples/interface_fcccu111.rst`, covering
  `region_b_kind='vacuum'`, the `nlay_a`/`nlay_b` dual-namelist trap, and the
  alignment solver with real numbers.
- **Known-incomplete, tracked in `tests/KNOWN_ISSUES.md`:**
  1. `(111)` interface DOS has an unexplained 2.05e-3 residual near the
     d-band peak (`(001)` is exact) — downstream of the Hamiltonian,
     unexplained, pinned rather than tolerated away.
  2. `A|vacuum` with >1 active layer diverges (step ≈ −334 Ry) —
     **untriaged**; a hand-built test deck is a live suspect, not confirmed
     as a code bug.
  3. Gate G-B7-3 (compensation-weight tolerance calibration) unsigned —
     needs buffer-thickness convergence data from a B7.7 that hasn't run.
  4. `A|vacuum-gap|B` needs no new code per maintainer decision (ordinary
     `A|active|B` with empty spheres in the active zone) but is untested.
  5. B7.7 (literature IEC validation, Fe/Cr and Co/Cu) is explicitly scoped
     as a separate, likely in-person/student task.

---

## 2. Plan/blueprint file location audit

Three separate trees currently hold overlapping copies of the same
per-feature (`B1`…`B12`) plan documents, at different stages of being kept
up to date. This is the actual mess, mapped out file by file.

### 2.1 `claude_files/` — entirely untracked, frozen at initial promotion

```
git ls-files claude_files/   → (empty — nothing here is in git at all)
```

This whole directory — `claude_files/00_INDEX.md`,
`claude_files/B1_gbt_frozen_magnons.md`, and everything under
`claude_files/plans/` (`00_INDEX.md`, `BLUEPRINTS.md`, `B1`…`B12`) — is a
**local working-tree artifact that was never committed**. It is the
original "promote all 12 blueprints into task plans in one pass" output.
Its checklists are still at their initial unchecked state
(`- [ ] B2.1 skeleton...`) even for items that have long since shipped —
compare `claude_files/plans/B2_reciprocal_green.md`'s checklist against
`docs/dev/plans/B2_reciprocal_green.md`'s (127-line diff, the latter has
per-task status prose and gate references the former doesn't).

Two files inside it are themselves stray duplicates of siblings one
directory up: `claude_files/00_INDEX.md` ≡ `claude_files/plans/00_INDEX.md`,
and `claude_files/B1_gbt_frozen_magnons.md` ≡
`claude_files/plans/B1_gbt_frozen_magnons.md` (both byte-identical).

### 2.2 `docs/dev/plans/` — the living fork, git-tracked, incomplete

```
tracked  docs/dev/plans/B2_reciprocal_green.md
tracked  docs/dev/plans/B3_bloch_spectral_functions.md
tracked  docs/dev/plans/B5_route_agnostic_postprocessing.md
tracked  docs/dev/plans/B7_interfaces_and_vacuum_leads.md
```

These four are real, actively-maintained implementation plans — updated
progress checklists, gate references, session handoff notes. B7's is even
under a **different filename** than its `claude_files/plans/` counterpart
(`B7_interfaces_and_vacuum_leads.md` vs `B7_boundary_self_energies.md`)
because the shipped design (region registry + frozen regions + vacuum lead)
diverged from the original López-Sancho-decimation blueprint — the old name
would actively mislead.

**Missing from this tree:** B1 and B6 have no `docs/dev/plans/` file at
all, despite both being done. `docs/DEVELOPER_MAP.md` line 50 references
`docs/dev/B1_GBT_SPIN_SPIRAL_PLAN.md` — **this file does not exist anywhere
in the repository.** It's a dead link (either never created, or created and
never committed, or renamed and the reference not updated).

### 2.3 `docs/dev/` (loose, non-`plans/`) — a different document genre, tracked

```
tracked    docs/dev/B2_GATE_G-B2-1.md
tracked    docs/dev/B2_GATE_G-B2-2.md
tracked    docs/dev/B2_RECIPROCAL_GREEN_HANDOVER.md
tracked    docs/dev/B5.3_gilbert_damping_audit.md
tracked    docs/dev/CONTRACT_FROZEN_REGION.md
tracked    docs/dev/CONVENTIONS_MADELUNG.md
tracked    docs/dev/MATH_AUDIT_T11.md
tracked    docs/dev/REFACTORING_PLAN.md
tracked    docs/dev/REVIEW_NOTES_T1-T6.md
tracked    docs/dev/reciprocal_green_convergence.md
tracked    docs/dev/route_agnostic_estimators.md
UNTRACKED  docs/dev/electrostatic_dipole_plan.md
```

These aren't plans — they're gates (sign-off records), conventions,
contracts, audits, and handover notes: a legitimate, distinct category from
the per-B task-plan documents in `plans/`. All correctly live directly
under `docs/dev/`. The one exception is `electrostatic_dipole_plan.md`
(the missing B6 plan, in substance) — it sits in the tracked directory but
was itself never `git add`ed, so it's invisible to anyone who clones the
repo fresh.

### 2.4 Root-level stray

```
UNTRACKED  BLUEPRINTS.md          (byte-identical to claude_files/plans/BLUEPRINTS.md)
```

`BLUEPRINTS.md` at repo root is untracked and is a plain copy of the one
inside `claude_files/plans/`. There is no tracked `BLUEPRINTS.md` anywhere
in the repository right now — the design document that `CLAUDE.md` and
`claude_files/plans/00_INDEX.md` both treat as load-bearing ("Phase-3
design document... promoted to a task plan when its turn comes") exists
only as an untracked working-tree file, in two copies, in the wrong
locations.

### 2.5 Summary table

| Current path | Tracked? | Kind | Verdict |
|---|---|---|---|
| `BLUEPRINTS.md` (root) | no | design doc | stray duplicate, should not exist at this path |
| `claude_files/00_INDEX.md` | no | index | stray duplicate of `claude_files/plans/00_INDEX.md` |
| `claude_files/B1_gbt_frozen_magnons.md` | no | plan | stray duplicate of `claude_files/plans/B1_gbt_frozen_magnons.md` |
| `claude_files/plans/00_INDEX.md` | no | index | authoritative *content*, wrong location (untracked) |
| `claude_files/plans/BLUEPRINTS.md` | no | design doc | authoritative *content*, wrong location (untracked) |
| `claude_files/plans/B1_gbt_frozen_magnons.md` | no | plan | only copy of the B1 plan; stale (no progress updates) |
| `claude_files/plans/B2/B3/B5/B7_*.md` | no | plan | **superseded** by the `docs/dev/plans/` versions — safe to drop |
| `claude_files/plans/B4/B6/B8/B9/B10/B11/B12_*.md` | no | plan | only copies of these plans; B6's is also superseded in substance by `docs/dev/electrostatic_dipole_plan.md` |
| `docs/dev/plans/B2/B3/B5/B7_*.md` | yes | plan | authoritative, current — keep as-is |
| `docs/dev/electrostatic_dipole_plan.md` | **no** | plan (B6) | authoritative content, needs `git add` + rename into `plans/` |
| `docs/dev/*.md` (gates/conventions/contracts/audits) | yes | reference | correct location already |

---

## 3. Homogenization — performed 2026-07-28

All 11 steps below have been executed (not just proposed). Nothing has been
`git add`ed yet — that's a deliberate pause so the diff can be reviewed
before staging; see the end of this section.

| # | Action | From | To | Done |
|---|---|---|---|---|
| 1 | Move | `claude_files/plans/BLUEPRINTS.md` (identical to root copy) | `docs/dev/BLUEPRINTS.md` | ✅ |
| 2 | Delete | root `BLUEPRINTS.md` | — | ✅ |
| 3 | Move + reconcile | `claude_files/plans/B1_gbt_frozen_magnons.md` | `docs/dev/plans/B1_gbt_frozen_magnons.md`, checklist annotated against verified reality (most sub-tasks' *named deliverables* — `test_gbt_block.f90`, `gbt_reference.py` — were not found in the repo, even though the kernel itself is shipped and in use; flagged honestly rather than guessed at) | ✅ |
| 4 | Move + rename + reconcile | `docs/dev/electrostatic_dipole_plan.md` (untracked) | `docs/dev/plans/B6_surface_electrostatics.md`, restructured as §1 "original blueprint spec" (from `claude_files/plans/B6_surface_electrostatics.md`) + §2 "as-shipped design" (the actual, narrower implementation) | ✅ |
| 5 | Move, unchanged | `claude_files/plans/{B4,B8,B9,B10,B11,B12}_*.md` | `docs/dev/plans/` verbatim (none started, nothing to reconcile) | ✅ |
| 6 | Drop (superseded) | `claude_files/plans/{B2,B3,B5,B7}_*.md` | — (never copied; the `docs/dev/plans/` versions were already the live ones) | ✅ (no-op by construction) |
| 7–9 | Delete stray dupes, preserve archive, remove directory | `claude_files/` (all of it, including `.DS_Store` and `rslmtoasa_phase3_plans.zip`) | zip copied to `/Users/andersb/Jobb/rslmto_fable_archive/rslmtoasa_phase3_plans.zip` (outside the repo, per explicit instruction to keep the original bundle rather than delete it); directory then removed | ✅ |
| 10 | Fix dead link | `docs/DEVELOPER_MAP.md:50` → `docs/dev/B1_GBT_SPIN_SPIRAL_PLAN.md` (never existed) | `docs/dev/plans/B1_gbt_frozen_magnons.md` §1.5 | ✅ |
| 11 | Reconcile | `docs/dev/BLUEPRINTS.md` §B7 | added a redirect note pointing at `docs/dev/plans/B7_interfaces_and_vacuum_leads.md`, kept the original design text below it as historical record | ✅ |

**Also updated in `docs/dev/plans/00_INDEX.md` while executing step 8**
(beyond a plain move): the gate registry now records actual sign-off status
for every gate (previously it only listed gate *existence* — this pass
checked each one against `docs/dev/*.md` gate/contract files and, for G-B1-*,
found no sign-off record at all, which the old table's blank "not yet"
framing didn't distinguish from "genuinely unreached"). G-B7-2 was added
(the table only had the obsolete G-B7-1). A new G-B7-3 row was added for the
compensation-weight-tolerance gate mentioned in `tests/KNOWN_ISSUES.md`
(previously untracked in any gate table). The file map gained a per-item
status column.

**Resulting layout:** `docs/dev/plans/` now holds one file per blueprint
item, B1–B12, all git-trackable (pending `git add`). `docs/dev/` (loose)
keeps its existing role — gates, conventions, contracts, audits, handovers —
unchanged. `docs/dev/BLUEPRINTS.md` is the one design-doc entry point.
`claude_files/` no longer exists. One item found during the sweep that was
**not** part of the 11-step plan and was left untouched: a root-level
`claude_files.zip` (an older, pre-`plans/`-subdirectory snapshot of the same
content, itself superseded by everything above) — flagged here rather than
deleted silently, since deleting untracked files beyond the agreed scope
wasn't authorized.

**Not yet done:** nothing in `docs/dev/plans/` or `docs/dev/BLUEPRINTS.md` is
`git add`ed yet. Recommend reviewing `git status`/`git diff --stat` for the
new tree before staging, since this touched ~20 files by move/copy/rewrite.

---

## 4. Documentation staleness (rst + top-level docs)

**Current and reliable** — verified accurate against the code:
- `docs/DEVELOPER_MAP.md` (aside from the one dead link in §3 item 10)
- `tests/KNOWN_ISSUES.md`
- `docs/source/user_guide/examples/interface_fcccu111.rst`
- `claude_files/plans/00_INDEX.md` is internally self-consistent (its gate
  table matches its own plan files) — it's just missing B7's later gate.

**Stale or silent on every B1–B7 addition:**
- `docs/source/theory/kspace_modes.rst` — no mention of `reciprocal_green`,
  `gf_route`, `green_backend`, `green_eta`, BSF, or `do_damping`.
- `docs/source/theory/green_functions.rst` — same; its `block_green_eta()`
  reference is an unrelated legacy real-space routine.
- `docs/source/theory/spin_dynamics.rst` — no mention of frozen-magnon
  post-processing or the GBT spin-spiral fix.
- `docs/source/keywords/control_parameters.rst:499-521` — documents
  `calctype` values `'bulk'/'surface'/'bands'/'exchange'/'conductivity'/'dos'`,
  which do not match the actual enum in `source/control.f90:635-638`
  (`'B'`, `'S'`, `'I'`, `'L'`) — predates this sweep, but `'L'` (new) is
  unsurprisingly absent too. No `gf_route`, `do_damping`, `green_backend`,
  `green_eta` documented anywhere in `keywords/`.
- `docs/source/keywords/hamiltonian_parameters.rst` — omits the new
  `gbt_kspace` namelist logical.
- No keyword page documents the new `&lattice` (`nlay_a`/`nlay_b`,
  `region_b_kind`), `&charge` (`fix_fermi_to_region`, `fermi_a`/`fermi_b`,
  `compensation_profile`, `bias`), `&control` (`dipole_electrostatics`,
  `dipole_mix`), or `&potential` (`q10`) additions. No `&frozen_magnon`
  keyword page exists at all.
- `docs/source/reference/module_overview.rst` — zero mentions of
  `reciprocal_green.f90`, `region_registry.f90`, `vacuum_lead.f90`,
  `sigma_provider.f90`, `electrostatics_multipole.f90`, `bsf_kernel.f90`,
  `dyson_kernel.f90`, `lehmann_kernel.f90`, `moment_kernel.f90`,
  `reciprocal_bsf.f90`, `reciprocal_moments.f90`, or `frozen_magnon.f90`.
- `docs/source/user_guide/examples.rst` — no page for frozen-magnon or BSF;
  the interface page exists but isn't cross-listed under a
  "post-processing" heading even though both are post-processing features
  with real passing regression cases.
- `tests/README.md` has zero references to any B1–B7 feature/test name —
  it's a generic CI doc, not a case index, so this is lower severity than
  the rst gaps but still worth a line each once the highlight examples in
  §5 exist.

---

## 5. Missing/orphaned highlight examples

Flagged for later action (not generated this session, per explicit
instruction to defer):

| Feature | Needs | Effort |
|---|---|---|
| **B3** (BSF) | new example built from scratch — pick a system, run, validate against band-structure eigenvalues at sharp δ | build + validate from zero |
| **B6** (dipole electrostatics) | new example built from scratch — a surface/interface system where the dipole shift is visible, ideally with a hand-checkable Q10 | build + validate from zero |
| **B2** (`reciprocal_green`) | promote `example/exchange/bccFe_kspace_green/` into `tests/postproc/cases.json` + write an rst page | promotion only, material already exists |
| **B5** (route-agnostic triads) | promote `tests/regression/triad_bccFe_damping/` + `run_triad.py` into `cases.json`/`tests/README.md` + rst page | promotion only, material already exists |

Note found while preparing to start this work (before it was deprioritized):
`Example_bulk_bccFe_nsp2_block_hoh` fails on a freshly rebuilt `build_13`
binary against current `HEAD`. This is a **pre-existing** failure per prior
session notes (`rslmto-test-baseline` memory), not introduced by anything in
this sweep — no source, test, or example files were touched this session.
