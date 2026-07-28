# B1 — Generalized Bloch Theorem spin spirals: bug fix + frozen magnons

**Status:** fresh audit complete (July 2026, against `fable_v2` @ `d86fe42`).
**Effort:** M (1–2 weeks). **Lead:** OPUS for the kernel, SONNET for tests
and plumbing. **Blocks:** B2 (phase conventions), B11 (validation target).

---

## 1. Why the current implementation is wrong — the audit

This section is maintainer-reviewed spec. Agents implement it; they do not
re-derive or second-guess it. Escalate disagreements, do not "fix" them.

**⚠ Supersedes `BLUEPRINTS.md` §B1.** The blueprint's shorthand
"H_ij → U(q·ΔR/2) H_ij U†(q·ΔR/2)" is a *two-sided* transform — the same
structural error as the current code (see E1 below), differing only in
using the half angle. If you have both documents open, **this file wins.**
The correct transform is one-sided; see §1.3 eq. (★).

### 1.1 What the code does today

`source/hamiltonian_build.f90 :: ham0m_nc` (≈ lines 1365–1490):

1. When `q_ss ≠ 0` (and no per-sublattice angles), it sets **both** pair
   moments to the rotating-frame cone direction:
   `mom_ia = mom_ja = (sinθ, 0, cosθ)`.
2. It builds the Pauli-decomposed pair block `hhmag(:,:,1:4)` from these
   moments (so `dot = 1` for θ-independent parts, `cross = 0` always).
3. It then rotates the transverse components by the **full** bond angle
   α = 2π q·ΔR/alat:
   `h1' = h1 cosα − h2 sinα; h2' = h1 sinα + h2 cosα`,
   leaving `h3` (z) and `h4` (charge) untouched.

The `d86fe42` attempt correctly moved to the true bond vector `vet` (fixing
the earlier absolute-positions error) but kept steps 1 and 3 — and those are
the actual physics errors.

### 1.2 The two errors

**(E1) Wrong transformation type.** Rotating the (h1,h2) Pauli components by
the full angle α is the *adjoint* action, U(α) H U†(α) — the transformation
of an **on-site** object where both spin indices live at the same site. A
**bond** block transforms one-sidedly: its row spin index rotates with site
i, its column spin index with site j, and those angles differ by exactly
the spiral phase. The adjoint form corresponds to φ_i = φ_j, i.e. no spiral
at all in the spin-mixing channels; only a spurious transverse twiddle
survives, and the scalar/z channels (which carry most of E(q) for a small
cone) have **no q-dependence whatsoever**. This is why frozen-magnon E(q)
comes out flat/garbage.

**(E2) Wrong pair moments.** In the correct construction (below) the pair
block must be evaluated with the site-j moment at the *relative* spiral
azimuth: m_j = R_z(α) m₀, not m_j = m₀. With m_j = m₀ the bond-dependent
parts of the scalar exchange (`dot`) and the Dzyaloshinskii-like `cross`
term are killed identically.

**(E3, latent trap) `momc` aliasing for elemental systems.** Lines ~1416:

```fortran
momc(it, :) = cmplx(mom_ia, ...)
momc(jt, :) = cmplx(mom_ja, ...)
```

and the magnetic components then use `momc(it,m)` for the m_i term and
`momc(jt,m)` for the m_j term. **When it == jt** (every elemental crystal —
bcc Fe, the canonical test!), the second assignment overwrites the first
and both terms use m_j. Today this is masked because mom_ia == mom_ja; the
moment E2 is fixed, this aliasing silently re-breaks everything for
elemental systems while leaving multi-species tests green. This is very
plausibly a contributor to why previous "obviously correct" attempts still
failed. **Fix: pass the two pair moments as explicit local complex vectors
(`momc_i(3)`, `momc_j(3)`), never through the type-indexed array.**

### 1.3 The correct construction (spec)

Setup: flat or conical spiral with wavevector **q**, cone angle θ. Lab-frame
moments m_i = R_z(φ_i) m₀ with φ_i = q·R_i and m₀ = (sinθ, 0, cosθ).
Spinor z-rotation D(φ) = exp(−i σ_z φ/2) = diag(e^{−iφ/2}, e^{+iφ/2}),
which satisfies the covariance H[R_z(φ)m_i, R_z(φ)m_j] = D(φ) H[m_i,m_j] D†(φ).

The rotating-frame (GBT) bond block is K_ij = D†(φ_i) H_ij[lab] D(φ_j).
Using covariance with R_z(φ_i) and α ≡ φ_j − φ_i = 2π q·ΔR/alat (code
units; ΔR = `vet`):

```
K(ΔR) = H_pair[ m₀ , R_z(α) m₀ ] · D(α)              (★)
```

i.e. **two ingredients, both mandatory**:

1. Build the existing noncollinear pair block with
   `mom_ia = m₀ = (sinθ, 0, cosθ)` and
   `mom_ja = (sinθ·cosα, sinθ·sinα, cosθ)`.
   Then `dot = cos²θ + sin²θ·cosα` and
   `cross = (−sinθcosθ·sinα, sinθcosθ(cosα−1)·(−1)…)` — do NOT hand-code
   these; the existing `dot_product`/`cross_product` on the two moment
   vectors produces them.
2. Right-multiply the assembled 2×2 spin block by the **half-angle** spinor
   D(α). In the code's stored components (h4 = charge, h1,h2,h3 = x,y,z),
   with c = cos(α/2), s = sin(α/2), this is exactly:

```
h4' = c·h4 − i·s·h3
h3' = c·h3 − i·s·h4
h1' = c·h1 + s·h2
h2' = c·h2 − s·h1
```

Notes the implementer must respect:

- The charge and z components **mix with imaginary coefficients**. That is
  correct: K(ΔR) for ΔR ≠ 0 is not Hermitian by itself; Hermiticity holds
  pairwise, K(−ΔR) = K(ΔR)†. The complex Pauli-coefficient storage and the
  assembly at `hamiltonian_build.f90` lines ~1199–1302 already support this
  with zero structural change. Do not "symmetrize" or take real parts
  anywhere.
- The transverse components rotate by the **half** angle (and in the
  opposite sense to the current code). Do not blend spec (★) with the old
  full-angle rotation.
- On-site block: ΔR = 0 ⇒ α = 0 ⇒ D = 1, moments equal — reduces to the
  collinear-like on-site with the tilted cone moment m₀. The existing
  `vv ≤ 0.01` on-site branch (cx0/cex0, cx1/cex1 additions) is then correct
  as-is provided it uses m₀.
- q → 0 and θ arbitrary: K reduces to the current (correct) noncollinear
  FM block. q_ss = 0 paths must remain **bit-identical** (regression
  contract).
- Per-sublattice branch (`theta_ss_sublattice/phi_ss_sublattice`): the
  sublattice angles define per-type m₀(t) in the rotating frame; m_j must
  still be R_z(α) m₀(jt) and the D(α) factor still applies. Same formula
  (★), just per-type m₀.
- `hxc` array (exchange-field-only, lines ~1206–1209) drops h4 — under (★)
  the h4/h3 mixing makes that split convention-dependent for bonds. The
  frozen-magnon workflow uses **band energy**, not `hxc`; therefore: leave
  `hxc` computed as today, and add a fatal guard if any torque/`hxc`
  consumer runs with q_ss ≠ 0 until a dedicated decision is made
  (escalation, not silent output).

### 1.4 Why GBT pays off in a real-space code

The whole point of (★): K depends only on ΔR, so the rotating-frame
Hamiltonian is translationally invariant even for incommensurate q — one
recursion per inequivalent site describes the entire spiral. Frozen-magnon
E(q) on any q, including incommensurate ones, at the cost of one
collinear-like run. The same K(ΔR) Bloch-summed gives spiral bands in the
`reciprocal_*` family for free (task B1.6), and B2 inherits the phase
convention pinned here.

### 1.5 Frozen-magnon recipe (Essenberger-style MFT), already scaffolded

`source/frozen_magnon.f90` + `calculation.f90 ::
post_processing_frozen_magnon` sweep q over a list, mode `'mft'` (force
theorem: converge reference at q₁, single band-energy pass per q) or
`'scf'`. Magnon energy normalized ω(q) = ΔE(q)/Δm_z with
Δm_z ∝ M_tot sin²θ (small cone; see the code comment citing PRB 111). The
workflow is sound **given a correct K(ΔR)**; the plan keeps it and validates
its normalization at gate G-B1-2.

---

## 2. Tasks

### B1.1 [SONNET] — Pin conventions + failing-first tests
Write `tests/unit/test_gbt_block.f90` (or extend the existing unit-test
harness) that constructs `ham0m_nc` output for a single bond of a 2-atom
toy and asserts, against hand-evaluated spec-(★) references generated by the
Python helper in §3.1:
(a) q = 0 bit-identity with current output;
(b) the four component formulas of §1.3 for one nonzero (q, θ, ΔR);
(c) Hermiticity pairing K(−ΔR) = K(ΔR)† through the assembled 18×18 block.
These tests must FAIL against current code (except (a)). Commit them
red-guarded (skip-marked) with message `test(B1.1): GBT spec tests (red)`.
*Session kit:* this file §1; `hamiltonian_build.f90` lines 1100–1500;
`tests/README.md` conventions section.

### B1.2 [OPUS] — Rewrite the spiral branch of `ham0m_nc` to spec (★)
Implement §1.3 exactly: explicit `momc_i/momc_j` locals (kills E3), correct
pair moments, D(α) right-multiplication via the four component formulas.
Keep the collinear/q=0 code path byte-identical (guard the new math behind
the same `q_ss/θ` predicate that exists today). Un-skip B1.1 tests; all
green. Also produce `tools/check_gbt_bond.py` — a ~50-line standalone
script that reads one bond block dumped by a debug flag and verifies (★)
numerically. **Gate G-B1-1:** Anders runs the script on one bcc-Fe bond and
signs off before anything downstream proceeds.
*Session kit:* this file §1 (all of it); `hamiltonian_build.f90` lines
1100–1520; B1.1 test file. Nothing else.

### B1.3 [SONNET] — Convention doc-test for q_ss / theta_ss units
One place, asserted: q_ss in Cartesian 2π/alat; α = 2π q·ΔR/alat; θ in
radians; the namelist conversion path (`include_codes/namelists/
hamiltonian.f90`, cf. churn in commits `cd83d25`, `13f3d17`). Add a unit
test that feeds q_ss = 0.5 x̂ on a cubic cell and asserts the on-file α for
the nearest-neighbor bond equals π·(x-component of ΔR)/…, i.e. locks the
factor once and forever. Update the namelist doc comments.
*Session kit:* this file §1.3 conventions bullet; the two namelist files;
`hamiltonian_build.f90` lines 180–200.

### B1.4 [OPUS] — Supercell known-answer validation
The canonical GBT test: for **commensurate** q, the spiral described by
K(ΔR) in the chemical cell must reproduce an explicit noncollinear
supercell with the moments set to the lab-frame spiral by hand.
Implement as a regression-style test (not unit): bcc Fe, flat spiral
θ = 90°, q = (0,0,0.5)·2π/alat ⇒ 2-cell AFM-like spiral; run (i) GBT
1-atom cell and (ii) explicit 2-atom supercell with moments +x̂/−x̂ …
(the true spiral pattern), same recursion depth; assert per-site charge,
m, and band energy agree to the regression tolerances. Repeat for one
incommensurate-looking commensurate q (e.g. 1/3) with a 3-cell.
Acceptance: both pairs agree; document tolerances chosen and why.
*Session kit:* this file §1.3–1.4; an existing regression case dir as
template (e.g. the bcc Fe case under `tests/regression`); `tests/README.md`.

### B1.5 [SONNET] — Frozen-magnon workflow validation + G-B1-2
With B1.2–B1.4 green: run the `frozen_magnon` post-processing on bcc Fe
along Γ–H, mode `'mft'`, small cone (θ = 10°–30°, check θ-independence of
ω(q) in the small-cone regime as an internal consistency test). Compare
against the **adiabatic route already in the code**: J_ij from the existing
LKAG `exchange.f90` module on the FM state → ω_ad(q) = (4/M)[J(0) − J(q)]
(sum J_ij e^{iq·R}). The two routes are independent implementations of the
same adiabatic physics; agreement at long wavelength (say q < 0.3 of the
zone) within a few percent is the acceptance. Additional blueprint-mandated
checks in the same sweep: (i) E(q) = E(−q) to SCF tolerance (no SOC in
these runs — any asymmetry is a phase-convention bug); (ii) E(q) smooth
through the zone boundary (sample q slightly beyond H and assert
continuity); (iii) small-q: ΔE(q) quadratic, with the coefficient matching
the spin-wave stiffness D from the J_ij sum rule — the second independent
cross-check. **Gate G-B1-2:** Anders signs
the normalization (factor conventions between ΔE/Δm_z and 4[J(0)−J(q)]/M)
before the workflow is declared production-ready. Add the Γ–H sweep as a
`tests/postproc` example with a golden file.
*Session kit:* this file §1.5; `frozen_magnon.f90`; `calculation.f90`
lines 1180–1360; `exchange.f90` J(q) accessor points (grep `jij`).

### B1.6 [SONNET] — Spiral bands via the reciprocal family (cheap win)
Bloch-sum K(ΔR) in `reciprocal_fourier.f90`'s existing machinery and expose
spiral band structures/DOS. No new physics — K already carries everything;
the Bloch sum uses the ordinary e^{ik·ΔR} phase (the spiral phases live
inside K). Known-answer: q = 0 reproduces existing FM bands bit-level;
supercell band-folding check reusing B1.4's commensurate cases.
*Session kit:* this file §1.3–1.4; `reciprocal_fourier.f90`;
`reciprocal_bands.f90` lines 1–120.

### B1.7 [GATE G-B1-3, then OPUS if needed] — SCF policy for q_ss ≠ 0
Decision task. mode='scf' at fixed q needs a policy: (a) constrained cone —
fix θ, converge magnitudes only (cheapest, standard for E(q,θ) maps); or
(b) free-angle self-consistency (θ from the converged local moment
direction). Anders decides; if (a), implementation is a small projection in
the moment-update step; if (b), a dedicated task is scoped then. Until
decided, mode='scf' with q_ss ≠ 0 logs a prominent warning naming this
gate.

---

## 3. Shared assets

### 3.1 Reference generator (write once in B1.1, reuse everywhere)
`tools/gbt_reference.py`: NumPy, ~80 lines. Inputs: w0, w1 vectors (fake
orbital scalars are fine for unit tests), θ, q, ΔR. Builds the pair block
per the code's own formulas with the correct pair moments, applies D(α),
prints the four components and the assembled 2×2⊗orbital block. This is the
oracle for B1.1/B1.2 and Anders's G-B1-1 spot check.

### 3.2 Progress checklist

**Updated 2026-07-28, reconciled against the actual `fable_v2` state during
the docs/dev/plans homogenization pass.** The kernel rewrite and the
frozen-magnon workflow both shipped (commits `aabb842`, `130b20d`, `63a8856`,
`49fc075`, `cbbb320`, `ce456e2`, and the `Example_bulk_bccFe_nsp4_block_
spiral_q{plus,minus}` E(q)=E(-q) regression pair in `tests/scf/cases.json`,
which cites this file). However **the specific deliverables named in
B1.1/B1.2 below — `tests/unit/test_gbt_block.f90`,
`tools/gbt_reference.py`, `tools/check_gbt_bond.py` — do not exist in the
repository.** Either the actual implementation validated the spec by a
different route (the supercell/E(q) regression route, most likely) or that
specific validation step was skipped. This was **not independently
re-derived or re-verified** in this pass — flagging honestly rather than
guessing. Anyone resuming B1 work should treat B1.1/B1.2's checkboxes below
as genuinely open until that's resolved, even though the kernel itself is in
production use.

- [ ] B1.1 spec tests (red) — **no `test_gbt_block.f90` found in the repo.**
- [ ] B1.2 kernel rewrite + G-B1-1 signed — **kernel rewrite appears to have
      shipped** (spiral branch is live in `hamiltonian_build.f90` and feeds
      the regression pair above), **but no `tools/gbt_reference.py` /
      `check_gbt_bond.py` found, and no record of G-B1-1 being signed** was
      located in `docs/dev/`. Unverified either way — check with Anders
      before assuming this gate is open or closed.
- [ ] B1.3 unit-convention doc-test — q_ss/theta_ss unit fixes did land
      (commits `13f3d17` "restore pi/alat halving", `cd83d25` "Updated q_ss
      convention, again"), but no dedicated unit test for the convention was
      located.
- [ ] B1.4 supercell known-answer (×2 q-points) — not located; the landed
      regression instead pins E(q)=E(-q) on a single-cell flat spiral, which
      is a different (weaker) check than an explicit supercell comparison.
- [x] B1.5 frozen-magnon vs adiabatic-J(q) — **workflow shipped** and is
      documented in `docs/DEVELOPER_MAP.md` and exercised by
      `Example_frozen_magnon_bccFe*` in `tests/scf/cases.json`. The
      multi-sublattice `branch_mode='auto'` acoustic branch is **known to
      violate the Goldstone condition at Γ** (~0.28 Ry gap, unexplained,
      deferred to B11) — see `tests/KNOWN_ISSUES.md`. **Whether G-B1-2 was
      formally signed is not recorded here** — no sign-off note was found
      in `docs/dev/`.
- [x] B1.6 spiral bands (reciprocal route) — `source/reciprocal_fourier.f90`
      contains spiral-related code (`aabb842` "Add first-order k-space GBT
      transform", `130b20d` "Extend k-space GBT to HOH and CCOR"); not
      independently re-validated against B1.4's known-answer cases in this
      pass since those cases weren't located either.
- [ ] B1.7 SCF policy decided (G-B1-3) — no record found of this decision
      having been made.

### 3.3 Cost/efficiency notes
The fix adds four complex axpy-scale operations per bond block — noise
relative to `hmfind` and recursion. No new allocations in the hot path
(the current code already allocates two norb×norb temporaries per call for
the rotation; reuse them). Frozen-magnon MFT sweeps are embarrassingly
parallel over q — expose the existing MPI q-loop if not already (check
`post_processing_frozen_magnon_auto`); do NOT add threading inside
`ham0m_nc`.
