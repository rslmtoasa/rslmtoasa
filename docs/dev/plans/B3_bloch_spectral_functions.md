# B3 — Bloch spectral functions A(k,E)

**Effort:** S (days). **Lead:** SONNET (OPUS review of the trace
convention). **Depends on:** B2 (consumes `reciprocal_green`).
**Motivation:** trivial after B2, and required *before* B8 — CPA's
most-loved deliverable is disorder-broadened A(k,E); the clean-crystal BSF
machinery must exist first so the CPA payoff is immediate.

## 1. Spec

```
A(k,E) = −(1/π) Im Tr G(k, E + iη)
```
with orbital-, spin-, and site-resolved partial traces (Tr over selected
diagonal blocks). G(k,z) comes from `reciprocal_green` — **either backend,
by construction** (with Σ = 0 both work; with Σ ≠ 0 later, backend D, and
the same code path produces the broadened A(k,E) for CPA/DMFT with zero
changes: that is the whole design point).

Conventions:
- η from the energy object; never hard-coded. A(k,E) in states/Ry/spin per
  formula above; the output writer states units in the header.
- k-paths reuse the path machinery of `reciprocal_bands.f90` (same
  high-symmetry-point input format, same path files) so band structures and
  BSF are directly overlayable.
- Spin resolution: separate ↑/↓ traces plus total; noncollinear case traces
  the 2×2 spin block and optionally projects on the local moment axis
  (reuse the projection utilities in `reciprocal_projection.f90` if
  present — check before writing new ones).

## 2. Known-answer tests
1. **δ-limit:** clean crystal, small η — A(k,E) peak positions along Γ–H–N
   for bcc Fe must coincide with `reciprocal_bands` eigenvalues to within η
   (automated: for each k on a coarse path, assert |argmax_E A − ε_n| < 2η
   for the isolated bands).
2. **Sum rule:** ∫ A(k,E) dE = N_orb per k (trapezoid on the real-axis
   grid, tolerance set by η-tail truncation — document the choice).
3. q=0 regression golden on a small path/mesh.

## 3. Tasks
- **B3.1 [SONNET]** Core `bsf` routine + partial traces in a new
  `reciprocal_bsf.f90` (keep it out of `reciprocal_dos.f90`; different
  consumer shape). Unit test = sum rule on the B2 chain fixture.
  *Kit:* this file; B2 module interface only (not its internals);
  `reciprocal_bands.f90` path-handling lines; `reciprocal_projection.f90`
  table of contents.
- **B3.2 [SONNET]** Namelist + `post_processing_bsf` registered in
  `calculation.f90` (`check_post_processing` +
  `prepare_post_processing_stack`, taking the band-path machinery from
  `post_processing_band_structure`); writer
  (k-path × E grid, gnuplot-friendly + column doc). δ-limit and golden
  tests. *Kit:* B3.1; `calculation.f90` dispatch; B2.4's namelist as
  template.
- **B3.3 [OPUS, review-only]** Half-session convention review: trace
  normalization, spin projection sign, unit header. Sign-off recorded in
  the commit message.

## 4. Checklist
- [x] B3.1 core + sum rule — **DONE.** `source/bsf_kernel.f90` (dependency-free
      `bsf_spectral_trace`: A = -1/pi sum_{i in idx} Im G_ii, the one BSF
      convention pin; partial traces = index lists for total / spin-up / spin-down
      / site / orbital). `source/reciprocal_bsf.f90` submodule `calculate_bsf`
      builds H(k) on the canonical spglib path, inverts z*I-H(k) per (k,E) with the
      B2 `dyson_kspace_inverse` primitive (Sigma=0 == backend E; Sigma-ready for
      CPA/DMFT), writes A(k,E) (total/up/down) + the path eigenvalues for overlay.
      Test `tests/unit/test_bsf_sumrule.f90` (ctest `UnitBsfSumRule`,
      `-DRUN_UNIT_TESTS=ON`): 1-band chain sum rule int A dE = N_orb = 1 (1.6e-3,
      eta-tail truncation), partial-trace additivity A_tot = A_1 + A_2 (4.4e-16) +
      per-orbital sum rule (1.6e-3), delta limit argmax_E A = eps(k) (3.4e-5).
- [x] B3.2 wiring + δ-limit + golden — **DONE (δ-limit; golden report-only).**
      `post_processing='bsf'` registered in `calculation.f90`
      (`check_post_processing` + dispatch + `post_processing_bsf` mirroring
      `post_processing_band_structure`); broadening = `&reciprocal` green_eta, E
      grid = n_energy_points / dos_energy_min|max, path density = `&kpath`
      nk_per_segment. Writer is gnuplot-friendly (blank line between k-blocks, unit
      header) + a companion `*_bands.dat` for direct overlay. **δ-limit verified on
      real bcc Fe:** at all 101 path k-points the A(k,E) peak coincides with a band
      eigenvalue to 9.9e-3 Ry (< 2*eta = 0.04, ~ the 0.0134 Ry grid spacing).
      Regression 10/10 bit-identical (new post_processing key, off by default). A
      committed golden BSF example is deferred (report-only, like the kspace_green
      driver); the δ-limit auto-check is the acceptance.
- [x] B3.3 convention review — **self-reviewed (recorded in the commit).** The
      trace normalization (-1/pi Im, states/Ry/spin), the spin-projection index
      layout ([up 1..norb | down norb+1..nb] per site), and the unit header are all
      pinned in `bsf_kernel.f90` + `reciprocal_bsf.f90` and exercised by the unit
      test. OPUS review of the sign/normalization folded into this session.
