# B12 — Couplings roadmap: e–magnon, e–phonon, phonon–magnon

**Effort:** XL cumulative / research. **Lead:** OPUS. **Depends on:** B2,
B10, B11. **Status:** architecture guidance promoted piecewise — **never as
one task.** This plan defines the three sub-programs, their first
deliverables, and the format/interface decisions that must be made early so
B2/B10 choices don't foreclose them (they currently don't; keep it that
way).

## 1. Electron–magnon (first sub-program to promote)

**Design.** Σ_em(z) from the RPA/ALDA χ (B11) convolved with the electron
GF — one more `sigma_provider` implementation (B10 interface, unchanged):
```
Σ_em(k, iω_n) ∝ Σ_{q, iν_m}  g² · χ^{+−}(q, iν_m) · G(k−q, iω_n − iν_m)
```
with g the local exchange-splitting vertex in the ASA blocks (same B_xc/m
data as B11's kernel — one convention, one accessor). Matsubara-first
(backend D takes complex z natively); real-axis via analytic continuation
is out of scope v1.

**First target:** spin-polaron renormalization of the bcc Fe BSF — B3
displays it with zero new display code (the design payoff, again).

**Validation:** sum-rule checks on Σ_em (Im sign, high-frequency tail
∝ 1/z with the exact moment); weak-coupling limit vs second-order
perturbation theory on the chain fixture.

## 2. Electron–phonon (interface, not implementation)

**Design stance (binding, consistent with the parked-FCD decision):**
phonons and deformation potentials are **not** computed in-code — do not
rebuild what plane-wave codes do better. The deliverable is a Σ provider
that *accepts externally computed* coupling data:
- v1: isotropic Eliashberg — read α²F(ω) (or λ, ω_log), build
  Σ_ep(iω_n) on Matsubara, feed backend D. Standard, well-oracled
  (analytic weak-coupling mass enhancement 1+λ test).
- v2 (only on demonstrated need): full g(k,q) grids from an external code;
  format decision then, with the provider interface already proven.

## 3. Phonon–magnon (different character: boson–boson, adiabatic route)

The in-code contribution is **adiabatic ingredients**, not dynamics:
dJ_ij/du tensors from finite-displacement J_ij calculations — the exchange
machinery on displaced clusters, a genuine real-space specialty (k-space
codes need supercells; recursion just displaces atoms in the cluster).
Dynamical spin-lattice coupling lives downstream in UppASD-style
simulations.

**Format decision (early, cheap, important):** the dJ/du export format is
defined **with the UppASD/mope/spindle consumer conventions in mind** —
the Hellsvik et al. PRB 99 sign convention, the mRy unit constant, and the
normalization already audited across those packages. Reuse that
conventions contract verbatim; add a doc-test asserting the sign on a
2-site toy. (This is the cross-package factor-of-two/sign class of bug —
treat it with the same audit discipline.)

**First deliverable:** finite-displacement dJ_ij/du driver (displacement
stencil, central differences, symmetry reduction of the tensor) + UppASD
`llfile`-compatible export.

## 4. Promotion rules
- Each sub-program promotes as its own plan with its own gates when its
  dependencies are green (1 needs B11+B10; 2 needs B10; 3 needs only the
  existing exchange machinery — **3 can promote earliest** and has
  immediate synergy with the UppASD SLD pipeline).
- Any new Σ provider must pass B2's backend-D plumbing tests before any
  physics validation is attempted.
- No sub-program ships without its analytic-limit oracle (weak-coupling,
  1+λ, or 2-site sign test respectively).

## 5. Checklist (promotion tracker)
- [ ] 12.3 dJ/du promoted (earliest; only needs exchange machinery)
- [ ] 12.2 Eliashberg provider promoted (after B10)
- [ ] 12.1 e–magnon Σ promoted (after B11)
