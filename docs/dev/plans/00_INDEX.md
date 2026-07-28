# RS-LMTO-ASA Phase 3 — Implementation Plans (B1–B12)

**What this is.** The promotion of `docs/dev/BLUEPRINTS.md` (Phase-3 design
document) into self-sufficient, delegable implementation plans. Each `B*.md`
file in this directory is written so that a coding agent given **only that
file, the repo-root `CLAUDE.md`, and the "session kit" file excerpts it
names** can execute its tasks without re-deriving physics or rediscovering
conventions.

**Repository state assumed:** branch `fable_v2` at or after commit `d86fe42`
("New effort on GBT. Still not working"). Phase 1 (structural refactor,
T1–T13) and Phase 2 (tests/CI/docs, P1–P10) complete; regression matrix and
example suites green.

**2026-07-28 homogenization note.** This directory previously existed in
three overlapping, partly-untracked copies: an untracked `claude_files/`
tree (the original bulk promotion, frozen at its initial state), this
tracked `docs/dev/plans/` tree (kept live for B2/B3/B5/B7 only), and a
stray root-level `BLUEPRINTS.md` duplicate. They have been reconciled into
this single tracked location — see `docs/dev/PHASE3_STATUS.md` for the full
audit and the reconciliation record. Two consequences worth knowing:

- `B7`'s plan file is named `B7_interfaces_and_vacuum_leads.md`, not
  `B7_boundary_self_energies.md` — the shipped design (region registry +
  frozen regions + vacuum lead) diverged from the original
  López-Sancho-decimation blueprint, and the old name would mislead.
- `B1` and `B6`'s plan files now carry "as-shipped vs originally-specified"
  reconciliation notes, since both features were implemented by a
  different, narrower route than their original text describes. Read their
  status headers before trusting the task lists below them.

---

## Dependency DAG and promotion order

```
B1 GBT fix ──────────────┐
                          ├─→ B11 LR-TDDFT (magnons)
B2 reciprocal_green ──┬──┤
   (two backends)     │   ├─→ B3 Bloch spectral functions
                      │   ├─→ B5 route-agnostic G(E) post-processing
B4 GPU k-space ───────┘   │      (Jij / damping / transport via Lehmann)
                          ├─→ B8 CPA + DLM (k-space)
                          └─→ B10 DMFT (Σ-provider API) ─→ B12 couplings
B6 surface electrostatics ─→ B7 boundary self-energies (leads + vacuum)
B7 ──→ B9 RS-CPA / RS-DLM (shares effective-medium machinery with B8)
```

**Campaign waves (from `docs/dev/BLUEPRINTS.md`, binding):**

| Wave | Items | Rationale |
|------|-------|-----------|
| 1 | B1, then B2 (+B4 in parallel) | keystone + conventions; B4 is independent mechanical work |
| 2 | B3, B5 | near-free payoffs of wave 1; establish the three-route regression triads |
| 3 | B6, then B7 | electrostatics before boundaries; B7 builds the per-energy RS driver |
| 4 | B8 (+ Σ-provider interface with B10 in mind) | CPA/DLM payoff; provider plumbing proven |
| 5 | B10, B11 | DMFT loop with external solvers; magnons |
| 6 | B9, B12, vertex-corrected CPA transport | research-adjacent tail |

Exception worth knowing: B12.3 (dJ_ij/du finite-displacement export for
UppASD) depends only on the existing exchange machinery and can promote
any time — see `B12_couplings.md` §4. B6 can run in parallel with B2–B5
(no shared files except `charge.f90`/SCF touchpoints).

**B1 first, and B1 is different.** B1 is a *bug fix with a completed fresh
audit* — the plan contains the exact diagnosis of why the current GBT
implementation (and the `d86fe42` repair attempt) is wrong, plus the derived
correct formulas in the code's own storage convention. Do not let any agent
"re-think" that derivation; it is maintainer-reviewed spec.

---

## Delegation policy (models and roles)

| Role | Model | Used for |
|---|---|---|
| **[OPUS]** | Claude Opus (or Fable) | Physics-critical kernels, convention-sensitive math, API design, anything a sign error would silently corrupt |
| **[SONNET]** | Claude Sonnet | Namelist plumbing, IO, CMake, test scaffolding, GPU boilerplate that follows an existing pattern, docstrings |
| **[GATE]** | Anders (maintainer) | Sign-off gates: convention pins, oracle numbers, physics decisions. No agent may proceed past an open gate. |

Every task in every plan carries one of these tags. Mixed tasks are split so
that the Opus part is a small, sharply-specified kernel and the Sonnet part
is the surrounding plumbing.

## Token-lean session protocol (applies to every plan)

1. **One task, one session, one commit** (Phase-1 ground rule 2 extended to
   sessions). A session receives: the plan file, repo `CLAUDE.md`, and the
   plan's *session kit* — an explicit list of files (with line ranges where
   given). Nothing else. Do not paste whole 1000+-line files when the kit
   names a range.
2. **The physics spec in the plan is authoritative.** If the agent believes
   the spec is wrong, it stops and writes a short escalation note for the
   maintainer instead of "fixing" the spec. This is the single most
   important rule: convention drift between sessions is how B1 failed
   the first time.
3. **Tests first.** Each task lists acceptance tests; the agent writes or
   extends the test *before* the implementation where feasible, and never
   deletes or loosens an existing tolerance to make a test pass.
4. **Run the regression matrix** (`tests/regression` + example suites)
   before starting and before committing. Bit-level behavior at q_ss = 0 /
   feature-off is the contract for every plan here.
5. **End-of-session handoff:** update the plan's progress checklist (each
   plan has one), note any surprises in ≤ 5 lines. No long narrative
   reports — the checklist plus the commit message is the record.
6. **Do not read the whole repo.** `docs/DEVELOPER_MAP.md` exists for
   orientation; consult it only when the session kit is insufficient, and
   prefer targeted `grep` over browsing.

## Standing conventions (pinned once, referenced everywhere)

- **Units:** energies in Ry internally (`rp` = real64); `q_ss` is given in
  the `&hamiltonian` namelist in **Cartesian units of 2π/alat** (q_ss = 0.5
  along a cubic axis is the zone boundary π/alat); bond phase
  α = 2π · q_ss·ΔR / alat with ΔR in absolute (alat-scaled) units.
  `theta_ss` stored internally in radians. Any plan that touches these
  re-asserts them in a doc-test/comment, never re-invents them.
- **Spin storage:** pair blocks are stored Pauli-decomposed in
  `hhmag(:,:,1:4)` = (x, y, z, charge) with assembly
  H = [[h4+h3, h1−i·h2],[h1+i·h2, h4−h3]] (see
  `hamiltonian_build.f90` lines ~1199–1302). All four components are complex
  orbital matrices; the basis {1, σ} with complex coefficients spans all
  2×2 spin structures — no storage change is ever needed to represent
  non-Hermitian bond blocks (Hermiticity holds pairwise, H_ji = H_ij†).
- **Architecture:** class-based Fortran (derived type + constructor +
  `restore_to_default` + type-bound procedures); new k-space code joins the
  `reciprocal_*` family; new plugins follow the C-API + `iso_c_binding` +
  CPU-fallback pattern established for CUDA/MKL kernels.
- **Fences:** `self.f90`, `symbolic_atom.f90` off-limits (edit around the
  edges via existing variables only — see B6 for the sanctioned pattern).
- **Known-answer discipline:** no estimator ships without a known-answer
  test; regression goldens are regenerated only by the sanctioned script
  path documented in `tests/README.md`.

## Gate registry (sign-off gates across all plans)

**Updated 2026-07-28** with actual sign-off status (previously this table
only listed gate existence, not whether each had been signed — see
`docs/dev/PHASE3_STATUS.md` for how each status below was verified).

| Gate | Plan | What Anders must sign off | Status |
|---|---|---|---|
| G-B1-1 | B1 | Numerical spot-check of one rotated bond block against the spec formulas | **unverified** — no sign-off record found; the named validation script (`tools/check_gbt_bond.py`) does not exist in the repo. See `B1_gbt_frozen_magnons.md` reconciliation note. |
| G-B1-2 | B1 | Frozen-magnon ω(q) normalization vs adiabatic-J(q) route on bcc Fe | **unverified** — workflow shipped and is in use, but no sign-off record found |
| G-B1-3 | B1 | SCF-mode policy for q_ss ≠ 0: constrained cone vs free-angle | **open** — no record of this decision being made |
| G-B2-1 | B2 | Energy-contour convention shared by both backends vs `green.f90` | **SIGNED** (Anders, 2026-07-13) — `docs/dev/B2_GATE_G-B2-1.md` |
| G-B2-2 | B2 | J_ij(R) k-mesh convergence acceptance thresholds | **SIGNED** (Anders, 2026-07-16) — `docs/dev/B2_GATE_G-B2-2.md`; default `green_eta=0.02` Ry |
| G-B4-1 | B4 | Precision policy for batched eigensolver (fp64-only vs mixed) | not reached — B4 not started |
| G-B6-1 | B6 | Which Skriver–Rosengaard reference tables define acceptance | **open** — B6 shipped at narrower scope without reaching this gate; no work-function validation performed |
| G-B7-1 | B7 | Principal-layer definition for the RS boundary recursion | **obsolete** — belonged to the original López-Sancho-decimation B7 design, which was superseded; see G-B7-2 |
| G-B7-2 | B7 | Frozen-region parameter contract (`vmad` persistence, required quantities, refuse-to-run policy) | **SIGNED** (Anders, 2026-07-26) — `docs/dev/CONTRACT_FROZEN_REGION.md`; note the refuse-to-run requirement was explicitly overridden in the same sign-off |
| G-B7-3 | B7 | Compensation-weight tolerance calibration | **open** — needs buffer-thickness convergence data from a B7.7 that has not run |
| G-B8-1 | B8 | CPA charge self-consistency scheme (SCF-CPA vs frozen-potential) | not reached — B8 not started |
| G-B10-1 | B10 | Σ(E) file/API format agreed with DMFT collaborators | not reached — B10 not started |

## File map

| File | Milestone | Effort | Lead model | Status (2026-07-28) |
|---|---|---|---|---|
| `B1_gbt_frozen_magnons.md` | GBT fix + frozen magnons | M | OPUS | done, gates unverified (see file) |
| `B2_reciprocal_green.md` | k-space GF, two backends | L | OPUS | done, both gates signed |
| `B3_bloch_spectral_functions.md` | BSF A(k,E) | S | SONNET | done, zero test/example/doc coverage |
| `B4_gpu_kspace.md` | Batched GPU eigensolver | M | SONNET | not started |
| `B5_route_agnostic_postprocessing.md` | Jij/damping/transport via Lehmann | M | OPUS | done |
| `B6_surface_electrostatics.md` | Monopole+dipole Madelung | M | OPUS | done at narrower scope; gate open; zero coverage |
| `B7_interfaces_and_vacuum_leads.md` | Leads + vacuum GF | L | OPUS | done, bulk of recent work; a few open items in `tests/KNOWN_ISSUES.md` |
| `B8_cpa_dlm_kspace.md` | k-space CPA + DLM | L | OPUS | not started |
| `B9_rs_cpa_dlm.md` | Real-space CPA/DLM | L | OPUS | not started |
| `B10_dmft_sigma_provider.md` | Σ-provider API (Hubbard-I first) | M | OPUS | not started |
| `B11_lr_tddft.md` | Transverse χ magnons | XL | OPUS | not started |
| `B12_couplings.md` | e-ph / e-magnon Σ roadmap | XL | OPUS | not started |
