# B11 — Linear-response TDDFT: transverse magnons

**Effort:** L–XL / medium-high risk (kernel + Goldstone hygiene).
**Lead:** OPUS throughout. **Depends on:** B1 (validation target +
conventions — hard), B2 backend E (χ₀ ingredients), B4 (χ₀ cost, GPU
batched GEMM). **Motivation:** true dynamical magnons ω(q) with Landau
damping, beyond the adiabatic picture; the flagship spectroscopy
deliverable.

## 1. Physics spec

**Bare transverse susceptibility** from eigenpair pairs across k and k+q
(B2 backend E products; the double-k structure is why B4's batched GEMM
matters):

```
χ₀^{+−}_{ij}(q, ω) = (1/N_k) Σ_k Σ_{n occ, m unocc}
   [ ψ†_{i,n↑}(k) ψ_{i,m↓}(k+q) ] [ ψ†_{j,m↓}(k+q) ψ_{j,n↑}(k) ]
   × ( f_{nk↑} − f_{mk+q↓} ) / ( ω − (ε_{mk+q↓} − ε_{nk↑}) + iη )
```

site/orbital-resolved in the ASA local blocks (i, j site indices; the
orbital structure kept at least site-diagonal for the kernel contraction —
v1: site-diagonal ALDA blocks, orbital-summed; orbital-resolved kernel is
a named extension).

**RPA/ALDA Dyson equation:**
```
χ(q,ω) = χ₀(q,ω) [ 1 − f_xc χ₀(q,ω) ]⁻¹
```
with the ASA-local transverse xc kernel f_xc = B_xc / m per site (the
adiabatic local approximation; extract B_xc and m from the converged
potential/moment data).

**Goldstone enforcement — document, don't hide.** The bare calculation
violates ω(q→0) = 0 by a finite amount (basis/kernel inconsistency);
enforce via the standard sum-rule kernel rescaling: scale f_xc per site by
λ chosen so χ(q→0, ω=0) has its zero mode exactly (equivalently the
static sum rule χ⁻¹(0,0)·m = 0). **Report the bare violation (meV) as a
quality metric in the output header** — it is a diagnostic of the
calculation, not an embarrassment to suppress.

**Observables:** Im χ^{+−}(q,ω) maps (magnon spectral function with Landau
damping visible as linewidth on entering the Stoner continuum); magnon
dispersion from peak positions; site-resolved composition for multi-
sublattice systems.

**Registered as** `post_processing_susceptibility` in `calculation.f90`
(insertion points recorded in `DEVELOPER_MAP.md`).

## 2. Validation — the three-way consistency (this is where B1 pays off)
1. **Goldstone:** ω(q→0) → 0 after enforcement; bare violation reported.
2. **Small-q slope/stiffness:** vs B1 frozen-magnon D and vs the J_ij
   sum-rule D — three independent routes to the same number; agreement
   within documented envelopes is the headline acceptance test.
3. bcc Fe and fcc Ni dispersions vs published LR-TDDFT and experiment
   (semi-quantitative; ASA-appropriate expectations documented).
4. Unit level: χ₀ on the 1-band chain fixture vs a brute-force
   supercell-response evaluation (small system, exact).

## 3. Efficiency spec
χ₀ cost: O(N_k · N_occ · N_unocc · N_site-blocks) per (q, ω) — batch over
ω (denominator-only change; numerators reused), GPU-GEMM the eigenvector
products (B4 primitive 3), MPI over q. Memory: never store the full
four-index χ₀; contract to site blocks on the fly. Frequency mesh from a
namelist with sane defaults (0 to ~1.5× the adiabatic bandwidth from B1's
ω(q)).

## 4. Tasks
- **B11.1 [OPUS]** χ₀ kernel on the chain fixture + brute-force oracle
  (test 2.4). Self-contained against B2 eigenpair access. *Kit:* this
  file; B2 eigenpair accessor interface; B4 GEMM primitive signature.
- **B11.2 [OPUS]** ASA f_xc extraction + Dyson inversion + Goldstone
  rescaling + quality-metric reporting. *Kit:* B11.1; `potential.f90`
  B_xc/moment accessors; this file §1.
- **B11.3 [SONNET]** `post_processing_susceptibility` wiring + namelist +
  ω/q mesh handling + Im χ map writer. *Kit:* B11.1/B11.2 interfaces;
  `calculation.f90` dispatch; B3's writer as format template.
- **B11.4 [OPUS]** bcc Fe campaign: tests 2.1–2.3, the three-way stiffness
  comparison, published-data comparison doc. *Kit:* B11.3; B1.5's
  frozen-magnon outputs; the J_ij stiffness from B5's triads.

## 5. Checklist
- [ ] B11.1 χ₀ + oracle
- [ ] B11.2 kernel + Goldstone
- [ ] B11.3 wiring
- [ ] B11.4 three-way validation + benchmarks
