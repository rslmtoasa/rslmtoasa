# B8 — CPA + DLM in k-space

**Effort:** L. **Lead:** OPUS (single-site solver + conventions), SONNET
(plumbing, staging of observables). **Depends on:** B2 backend D (hard),
B3 (BSF payoff), B10 (Σ-provider interface co-designed — CPA is its first
real implementation and integration test).

## 1. Physics spec (Turek et al. TB-LMTO-CPA formulation)

Single-site CPA on each disordered sublattice: coherent potential function
𝒫_c(z) (equivalently coherent Σ_CPA(z) block) determined by the
self-consistency condition

```
Σ_c  x_c · t_c(z) = 0 ,   t_c = [ (P_c − 𝒫)⁻¹ − Ḡ_loc ]⁻¹
Ḡ_loc(z) = (1/N_k) Σ_k [ z − H_k(𝒫) ]⁻¹ |_site-block
```

i.e. the average single-site t-matrix of the real constituents embedded in
the coherent medium vanishes. Iteration: fixed-point with linear mixing
(mixing parameter namelist-controlled, default 0.2; Anderson acceleration
optional later — do not implement in v1). Ḡ_loc via **B2 backend D** with
`sigma_cpa` as the Σ provider — no new k-machinery.

**Convention pins (each gets a unit test):**
- CPA operates on the **potential-function / auxiliary-GF level** of the
  TB-LMTO representation (Turek's convention); the transform to physical G
  reuses the existing `auxiliary_gij`-family machinery. Do not mix
  representations mid-loop.
- Per-sublattice: only sites flagged disordered carry 𝒫; ordered sites
  keep their sharp P.
- **DLM = CPA with orientational disorder:** same atom, moment up/down at
  50/50 (collinear DLM v1; Lebedev-mesh/Weiss-weighted generalization is a
  named later task, not v1). Reuses the noncollinear machinery for the
  rotated single-site P.

**Charge self-consistency (Gate G-B8-1):** v1 ships frozen-potential CPA
(one-shot on converged constituent potentials, standard first step). The
SCF-CPA scheme (component-resolved charges via the CPA-averaged GF,
Poisson with net-charge handling per sublattice) is scoped as its own task
after Anders signs the scheme choice.

## 2. Observables staging — honesty about vertex corrections

Ship in this order; each stage has its own tests:
1. CPA-DOS, component-resolved DOS, total energy, moments.
2. Disorder-broadened **BSF** — B3 handles Σ automatically (zero new code;
   the payoff of the B2/B3 design).
3. Disorder-averaged J_ij (Lichtenstein-CPA, single-site vertex only) —
   documented as such.
4. **Transport with CPA is wrong without vertex corrections, by
   construction.** Until a Velický-ladder task exists, the code must
   refuse (fatal with message) on CPA + conductivity. This guard is part
   of stage 1, not an afterthought.

## 3. Validation (known-answer)
1. x → 0 and x → 1 limits reproduce the pure crystals (bit-level vs the
   clean backend-D run at Σ = `sigma_zero`).
2. **Dilute limit vs the code's own real-space single-impurity mode** — a
   self-contained cross-check most codes cannot do: 1% CPA vs one impurity
   embedded in the host, local DOS on the impurity site.
3. Published CuNi and FeV CPA-DOS benchmarks (semi-quantitative,
   figure-level agreement documented).
4. DLM bcc Fe: finite local moment, vanishing net magnetization; local
   moment magnitude vs literature (~1.9 μ_B ballpark — cite in the test
   doc).
5. Analyticity of 𝒫(z) (Im-part sign) asserted every iteration.

## 4. Tasks
- **B8.1 [OPUS]** `sigma_cpa` provider skeleton against B10's interface +
  single-site solver on a 1-band toy alloy (analytic 2-component
  Hubbard-binary check: CPA band edges/DOS vs published closed-form
  curves). Self-contained. *Kit:* this file; B10 interface file; nothing
  else.
- **B8.2 [OPUS]** Representation wiring: 𝒫 ↔ potential parameters,
  auxiliary↔physical transforms, per-sublattice bookkeeping; convention
  unit tests. *Kit:* B8.1; `green.f90` auxiliary transform routines;
  `potential.f90` P-function accessors (±40 lines each).
- **B8.3 [SONNET]** Loop driver + mixing + namelist + limits test (3.1) +
  the conductivity guard. *Kit:* B8.1/B8.2; B2.5 namelist pattern.
- **B8.4 [OPUS]** DLM mode + test 3.4. *Kit:* B8.3; the noncollinear
  moment-direction plumbing from `hamiltonian_build.f90` (moments region).
- **B8.5 [SONNET]** Dilute-limit impurity cross-check (3.2) + CuNi/FeV
  benchmark cases (3.3) + component-DOS/BSF examples. *Kit:* B8.3/B8.4;
  the impurity example (B2FeCo) as template; B3 interface.
- **B8.6 [OPUS]** Lichtenstein-CPA J_ij (stage-3 observable), documented
  single-site-vertex caveat. *Kit:* B8.3; `exchange.f90` contour loop;
  this file §2.3.

## 5. Checklist
- [ ] B8.1 solver + toy oracle
- [ ] B8.2 representation pins
- [ ] B8.3 loop + limits + guard (G-B8-1 decision recorded)
- [ ] B8.4 DLM
- [ ] B8.5 dilute cross-check + benchmarks
- [ ] B8.6 CPA-J_ij
