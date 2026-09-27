# Route-agnostic G(E) estimators — agreement envelopes (B5)

**Milestone:** B5 (route-agnostic post-processing). **Branch:** `fable_v2`.
**Status:** B5.1 (exact σ moments) and B5.2 (J_ij + σ triads) landed; the α
(Gilbert damping) triad and its η→0 protocol are reserved for B5.3.

This is the acceptance record for "the same G(E) observable produces the same
answer whether the Green's function / moments came from the real-space recursion
route or the k-space Lehmann/Dyson engine." It documents, per observable, the
**agreement envelope** within which the three routes are pinned as regression
cases (`tests/regression/run_triad.py`, `triad_cases.json`).

## The three routes

| Route | Producer | Σ | Notes |
|---|---|---|---|
| `recursion` | real-space Chebyshev/block recursion (legacy) | — | bit-identical default |
| `lehmann` | k-space backend E (strict Lehmann eigenpairs) | 0 | full unreduced BZ |
| `dyson` | k-space backend D (direct `[zS−H−Σ]⁻¹`) | 0 here | Σ≠0 ⇒ CPA/DMFT (B8/B10) |

Two qualitatively different agreements exist, and the triads test both:

1. **`lehmann` ≡ `dyson` (tight, permanent invariant).** With Σ = 0 backend D
   reproduces backend E to solver tolerance. This is the strongest cross-route
   statement and holds at machine/solver precision regardless of k-mesh or
   broadening. Pinned tightly (J_ij: `1e-5`; σ: `1e-6`).
2. **`recursion` vs `lehmann` (documented envelope).** These agree only within a
   **documented band**, because the two routes construct and broaden G(E)
   differently. The band is observable-specific (below). It is *documented*, not
   machine-precision — do not tighten it without a convergence argument.

## J_ij — exchange couplings

`post_processing='exchange'` + `gf_route` fills the SAME canonical `green` arrays
(`gij/gji`, the `gij_eta` ladder, torque families) all three routes fill, so
`calculate_exchange` runs unchanged. The full k-mesh / broadening convergence
study is `docs/dev/reciprocal_green_convergence.md` (gate G-B2-2); the triad here
pins route agreement + reproducibility on a CI-tractable mesh, not the converged
physical J.

**Envelope basis (from B2.6 / G-B2-2).** J_ij is a Fermi-surface-weighted
integral `∝ ∫^{E_F} Im Tr[dG_ij dG_ji] dE`. Unlike the DOS sum rule it is
**strongly broadening-dependent** for a metal: the recursion route broadens with
a near-sharp continued-fraction terminator (η ≈ 0), the Lehmann/Dyson routes with
the explicit retarded `green_eta`. So `recursion` and `lehmann` do **not** agree
tightly at fixed finite η — the mismatch is a *broadening/k-convergence artifact,
not a normalization error* (three independent lines of evidence in B2.6 §4; the
intersite kernel carries no normalization error, pinned `<1e-12` by
`UnitGammaSupercell`). The ratio is shell- and η-dependent (rules out any global
factor): `recursion/lehmann` ≈ 2.0 on the 1st-NN bond, ≈ 1.2 on the 2nd-NN bond.

**Triad case `triad_bccFe_jij`** (tracked `triad_bccFe_exchange` potential, 12³
mesh, η = 0.02 Ry, 1st- and 2nd-NN shells; ×10³/4π mRy):

| shell | pair | recursion (η≈0) | lehmann (η=0.02) | dyson | lehmann/recursion |
|---|---|---|---|---|---|
| 1st-NN | (1,335), a/2·(1,1,1) | 0.50788 | 0.25474 | 0.25474 | 0.50 |
| 2nd-NN | (1,336), a·(0,0,1)   | 0.38619 | 0.31323 | 0.31323 | 0.81 |

- `lehmann` ≡ `dyson` to ~3e-12 (tight pin).
- `recursion` vs `lehmann`: same sign, ratio in the documented band `[0.35, 1.15]`
  — the broadening envelope. To recover the near-sharp recursion J from the
  k-space routes, take η → 0 **jointly** with N_k → ∞ (B2.6 Tables A/B).

## σ — Kubo/KPM conductivity

`post_processing='conductivity'` + `gf_route`. The KPM pipeline consumes the
velocity-sandwiched double Chebyshev moment `μ_nm = <site|T_m(H̃) v_a T_n(H̃) v_b|
site>`, not `gij`. B5.1 added the **exact** k-space moment generator
(`reciprocal%fill_moments`, `moment_kernel.f90`): in the eigenbasis
`T_p(H̃(k)) = U diag(T_p(ε̃)) U†`, so the moments are exact from the eigenpairs and
fill the SAME `mu_nm_stochastic` array — `calculate_conductivity_tensor` is
untouched. The k-space velocity `v(k)=FT[v(R)]` is the recursion route's own
real-space velocity blocks Fourier-summed like `H(k)=FT[ee]`, so both routes
evaluate one operator in two representations.

**Envelope basis (the direct KPM error bound).** Because it is exact-vs-exact
(no broadening asymmetry like J_ij), `recursion` and `lehmann` σ agree **tightly**
— up to k-mesh sampling and the shared Chebyshev truncation. `UnitMomentKernel`
pins the eigenbasis moment identity to ~5e-16; the end-to-end σ difference is the
direct measurement of the KPM error on the same crystal.

**Scope pins.** Compare a **diagonal** component (σ_xx, `v_alpha == v_beta`): the
off-diagonal anomalous-Hall σ_xy is a symmetry-forced zero without SOC and is pure
coarse-mesh noise — not a route-agnostic observable. Exact generator scope:
`cond_type='charge'`, `cond_calctype='per_type'`, first-order H(k)/S=I; `dyson` ≡
`lehmann` (Σ=0 eigenpair path). Σ≠0 transport needs vertex corrections (B8) and is
out of scope.

**Triad case `triad_bccFe_sigma`** (tracked `triad_bccFe_conductivity` potential,
6³ mesh, cond_ll=30, σ_xx near E_F):

| route | σ_xx(E_F) | vs recursion |
|---|---|---|
| recursion | 3.23996 | — |
| lehmann   | 3.23505 | 0.15 % |
| dyson     | 3.23505 | 0.15 % |

- `lehmann` ≡ `dyson` bit-identical (same Σ=0 eigenpair moment path).
- `recursion` vs `lehmann`: ratio in `[0.97, 1.03]` (~0.15 % on this coarse mesh;
  tightens with N_k). This IS the KPM error bound.

## α — Gilbert damping (B5.3 / VAL-08)

The Kamberský torque-correlation damping consumes the canonical energy-grid
`green%gij/gji` arrays filled by all three route dispatches. The current
operator audit and convergence campaign are recorded in
[`VAL-08_DAMPING_INERTIA.md`](VAL-08_DAMPING_INERTIA.md).

- α is evaluated at the energy point nearest (E_F), with the damping tensor
  prefactor and SOC-derivative torque convention documented in VAL-08. It is
  broadening-defined, so recursion versus Lehmann is an eta/k-mesh envelope,
  not equality at fixed finite eta. The current on-site bcc-Fe measurements
  are (0.001341155) (recursion, (N_k=8^3)), (0.002527619) (Lehmann,
  (N_k=12^3)), and (0.002682953) (Lehmann/Dyson, (N_k=8^3,\eta=0.02)).
  These are convergence evidence, not an intrinsic single-eta value.

## Running the triads

```
python3 tests/regression/run_triad.py --binary build/bin/rslmto.x \
  --cases-json tests/regression/triad_cases.json --case-name triad_bccFe_jij \
  --scratch-root <scratch> --references tests/regression/references_triad
# --gen-ref (re)writes references_triad/<name>.json goldens.
```

The hard assertions are the cross-route ones (self-contained per run): `lehmann ≡
dyson` (tight) and `recursion` vs `lehmann` (documented band). The per-route
golden is a loose reproducibility guard (`golden_rtol=1e-2`) tolerant of
cross-BLAS eigensolver variance; regenerate it if the trusted numbers move for a
justified reason.
