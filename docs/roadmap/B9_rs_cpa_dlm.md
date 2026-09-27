# B9 — CPA / DLM for the real-space route

**Effort:** L–XL / high risk (research-adjacent). **Lead:** OPUS.
**Promote only after B7 + B8 exist** — B7 supplies the per-energy recursion
driver; B8 supplies the single-site solver, convention pins, and reference
results. **Motivation:** disordered *hosts* for the impurity machinery and
inhomogeneous disorder (layer-resolved concentration profiles at
interfaces) — exactly where the real-space identity shines and k-space CPA
cannot follow.

## 1. Design spec

**Effective-medium recursion.** Disordered sites carry Σ_CPA(E) (the
coherent 𝒫(E) from B8's single-site machinery), but the local GF entering
the single-site condition comes from **recursion on the effective medium**
instead of a k-sum:

```
iterate:  𝒫(E) → recursion on medium H[𝒫(E)] → Ḡ_loc(E) → single-site
          condition (B8's solver, unchanged) → 𝒫'(E) → mix
```

Energy-dependent on-site 𝒫(E) ⇒ the medium recursion runs **per energy**
— the exact cost pattern of B7's boundary recursion. **Reuse B7.3's
per-energy driver verbatim** (its callback mutates on-site blocks per
energy; that was designed for this). Do not write a second driver.

**Modes, staged:**
1. Homogeneous bulk RS-CPA (validation vehicle; k-space CPA exists, so
   this mode's only purpose is the two-route consistency proof).
2. **Impurity-in-alloy:** a real defect embedded in the converged CPA host
   via the existing Dyson impurity machinery — the scientific payoff.
3. **Inhomogeneous CPA:** layer-resolved concentrations at a surface/
   interface, 𝒫 per layer, composing with B7 boundary self-energies.

**Efficiency reality check (state in the namelist docs):** per-energy
recursion over the full contour times CPA iterations is expensive —
O(N_E × N_iter) recursions per medium site class. Mitigations, in scope:
coarse-to-fine energy sweeps (converge 𝒫 on a coarse mesh, refine),
warm-starting 𝒫 across adjacent energies, and the Chebyshev/GPU recursion
backends where applicable. Mitigations, out of scope v1: anything
smarter (rational interpolation of 𝒫(E), etc.) — note as future work.

## 2. Validation
1. **Two-route consistency:** homogeneous-bulk RS-CPA vs B8 k-space CPA,
   same alloy (CuNi), component DOS and total moments within documented
   recursion-depth/k-mesh envelopes.
2. Dilute limit vs impurity mode (mirror of B8's test 3.2, now from the
   RS side).
3. x → 0/1 limits ≡ clean recursion runs (bit-level).
4. Inhomogeneous smoke test: step concentration profile across an
   interface; 𝒫 far from the interface ≡ the two bulk-CPA values.

## 3. Tasks
- **B9.1 [OPUS]** Medium-recursion loop = B7.3 driver + B8.1 solver glue;
  homogeneous mode; tests 2.1/2.3 on a 1-band toy first, then CuNi. *Kit:*
  this file; B7.3 driver interface; B8.1 solver interface; `recursion.f90`
  entry points.
- **B9.2 [OPUS]** Impurity-in-alloy embedding via the existing Dyson
  impurity machinery; test 2.2. *Kit:* B9.1; the impurity-mode code path
  (locate via DEVELOPER_MAP, impurity workflow section).
- **B9.3 [OPUS]** Inhomogeneous/layered mode + B7 composition; test 2.4.
  *Kit:* B9.1/B9.2; B7 boundary interface.
- **B9.4 [SONNET]** Namelists, examples, cost-guidance docs (§1 reality
  check), regression cases.

## 4. Checklist
- [ ] B9.1 homogeneous + two-route consistency
- [ ] B9.2 impurity-in-alloy
- [ ] B9.3 inhomogeneous/layered
- [ ] B9.4 docs + regressions
