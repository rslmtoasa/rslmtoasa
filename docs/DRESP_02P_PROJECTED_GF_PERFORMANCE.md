# DRESP-02P — Projected GF Performance Remediation

Status: `PASS — PERFORMANCE BLOCKER REMOVED`

DRESP-02P preserves the certified projected reciprocal real-axis GF/Kubo
response and makes the DRESP-02R bcc-Fe workload executable.  The material
GF/Lehmann closure itself remains unresolved; this document returns that gate
to DRESP-02 without changing its numerical controls or its physics status.

## 1. Locked DRESP-02R workload

The workload was taken from `docs/DRESP_02R_MATERIAL_GF_CLOSURE.md` and run
against the frozen accepted bcc-Fe state: full `4x4x4` mesh (64 k points),
`ham_only`, second-order Hamiltonian, orthogonal, collinear, no SOC, 300 K,
`EF=-0.0612124078445383 Ry`, and `response_lmax=4`.

| control | locked value |
| --- | ---: |
| response eta | `0.01 Ry` |
| base integration eta | `0.002 Ry` |
| fixed-window energy interval | `[-1.3638307295, 1.6141834982] Ry` |
| fixed-window mesh | `NE=1001, 2001, 4001` |
| controlled eta ladder | `0.008, 0.004, 0.002 Ry` |
| controlled `h/eta_int` | approximately `0.40` |
| controlled mesh | `NE=933, 1863, 3725` |
| window ladder | margins `0.30, 0.60, 1.00 Ry`, `NE=4001` |
| diagnostic point | `q=(0,0,0)`, `omega=0 Ry` |
| selectors | `d`, `spd` |
| quadrature | unchanged composite Simpson |

The unchanged projected GF implementation estimated the fine fixed-window
d sample at `103.303 s`.  The same sample on the optimized backend took
`0.803 s`; the complete driver rerun, including accepted-state setup and all
audit output, took `14.10 s` with one OpenMP thread.

## 2. P0 baseline profile

The baseline profile was measured from the unchanged implementation before
the optimized route was selected.  The material fine sample has
`NE*Nk = 4001*64 = 256064` k-energy points and one requested frequency.

| region | baseline evidence |
| --- | --- |
| GF/resolvent generation | `1,536,384` dense `build_weighted_resolvent` calls (`6*NE*Nk`); dominant repeated work |
| endpoint preparation | one site-component vertex tensor per request; the endpoint eigenbasis was not retained for the energy loop |
| DRESP-01 vertex construction | hoisted once before the energy loop, but consumed through repeated dense matrix products |
| transformation to site/operator representation | implicit in every dense bubble contraction; no explicit reusable site transition data |
| first Kubo term | `256064` bubble calls, each traversing all 16 affine component pairs |
| second Kubo term | `256064` bubble calls, with a second dense contraction for each component pair |
| energy integration | `256064` k-energy iterations, with fixed Simpson weights and the locked finite window |
| allocation/copy overhead | six `nbasis x nbasis x 3` GF work arrays plus per-component temporary matrices; the site result remained small |
| k loop | 64 k points per energy; no k parallelism |
| frequency loop | one physical frequency in the Fe workload; two additional dense resolvents per k-energy-frequency point |

The profile used wall-clock timing exposed by the response result.  `perf`
was also attempted on Ubuntu, but this environment has
`perf_event_paranoid=4`, so hardware counters are unavailable without an
external system permission change.

## 3. Optimization ladder

P1 was already present in the certified baseline: the DRESP-01 component
vertices were built outside the energy loop, and the two spectral resolvents
were reused across the frequency loop.  The implemented changes are P2 and
P3.

| stage | transformation | result |
| --- | --- | --- |
| P0 | unchanged dense-resolvent profile | reference retained as `dense-reference` |
| P1 | invariant vertex/spectral preparation audit | no physics change; existing hoists retained |
| P2 | transform each affine site vertex into the left/right endpoint eigenbases once per k; evaluate GF factors as scalar spectral/retarded/advanced denominators | zero dense resolvent-builder calls in production |
| P3 | combine the four affine components into one exact site transition matrix and form each Kubo term with vectorized site-space BLAS contractions | two site contractions per energy/k/frequency point |
| P4 | OpenMP/MPI/GPU expansion | not needed: the locked workload is already sub-second per fine projection; no MPI or GPU machinery was introduced |

The public `evaluate_projected_gf_chi0` seam now selects the eigenbasis
implementation.  The original implementation is retained as
`evaluate_projected_gf_chi0_reference` and is used by the unit certification.

## 4. Mathematical equivalence

For an already diagonalized endpoint Hamiltonian, the reference weighted
resolvent is

\[
 G^{(p)}(z)=U\,\mathrm{diag}\left(\frac{\epsilon_n^p}{z-\epsilon_n}\right)U^\dagger .
\]

For each site and endpoint k pair, the optimized backend first forms the exact
affine vertex matrix element

\[
 T_i^{nm}(k,q)=\sum_{p,q=0}^{1}
 \epsilon_{n,k}^{p}\epsilon_{m,k+q}^{q}
 \langle n,k|V_i^{(p,q)}|m,k+q\rangle .
\]

Linearity of the trace then gives each preserved Kubo term as a site-space
contraction over the same band-pair elements:

\[
 \chi^{(1)}_{ij}(E,\omega)=
 \sum_{nm} A_n(E)G_m^R(E+\omega)T_i^{nm}T_j^{nm*},
\]

\[
 \chi^{(2)}_{ij}(E,\omega)=
 \sum_{nm} A_m(E)G_n^A(E-\omega)T_i^{nm}T_j^{nm*}.
\]

The fixed Simpson energy integration, finite window, Fermi factor, physical
response eta, integration eta, both Kubo terms, q endpoint convention, and
the `d`/`spd` selector semantics are unchanged.  No integral was analytically
collapsed into a transition denominator and no Lehmann accumulator is called
by the GF backend.

## 5. Correctness evidence

The predetermined reference-comparison tolerance is `2e-10` absolute for
the complete site matrix and for each Kubo term.  The maximum differences
between the optimized and dense-reference projected GF results were:

| fixture / selector | full matrix | Kubo term 1 | Kubo term 2 |
| --- | ---: | ---: | ---: |
| one-site d, Gamma and finite q | `1.21e-38` | `3.54e-38` | `3.58e-38` |
| one-site spd, Gamma and finite q | `6.52e-18` | `5.20e-18` | `3.79e-18` |
| two-site d, Gamma and finite q | `2.51e-38` | `1.21e-38` | `1.47e-38` |
| two-site spd, Gamma and finite q | `5.75e-18` | `1.22e-17` | `1.74e-17` |

The complete two-site matrix, q=0, finite q, both selectors, positive
frequency sign, and q/−q covariance all pass.  The largest optimized GF
q/−q covariance residual in the focused unit fixture is `3.50e-19`, below
the existing `2e-8` covariance tolerance.  The independent product-GF oracle
also remains within the existing fixture limits; its largest reported
optimized residual is `4.39e-18`.

The focused regression passed all 8 tests:

```text
cmake --build build --target UnitLrProjectedReciprocalChi0 rslmto.x -j2
ctest --test-dir build -R 'UnitLrProjectedReciprocalChi0|UnitTddftProductionDriver|UnitLrProduct(KsSusceptibility|GfSusceptibility)$' --output-on-failure
```

## 6. Performance results

The focused unit profile before the change was `13.44 s` serially with
approximately `39.8 MB` RSS.  The optimized backend itself reduced the
projected fixture work to sub-second scale; the unit now also invokes the
dense reference for every correctness fixture, so its total test time is not
a direct optimized-only benchmark.

| workload | baseline wall | optimized wall | speedup | k-energy points/s, optimized |
| --- | ---: | ---: | ---: | ---: |
| Fe d, fixed eta, `NE=4001` | `103.303 s` | `0.803 s` | `128.7x` | `3.19e5` |
| Fe spd, fixed eta, `NE=4001` | not available in the partial baseline artifact | `0.803 s` | — | `3.19e5` |
| all 18 GF audit samples, d+spd | partial baseline only | `10.02 s` GF sample sum | — | — |

The optimized profile for the fine fixed-window samples records approximately
`0.075 s` scalar spectral/physical denominator work and `0.674 s` site
accumulation work, with zero dense resolver calls and `256064` site
accumulator calls.  Endpoint transforms and vertices are retained at site
resolution; no `232x232` susceptibility or point-response matrix is
allocated.  The production result remains `N_site x N_site`.

## 7. Exact Fe convergence rerun

The exact DRESP-02R campaigns were rerun in `/tmp/dresp02p_fe` with the
optimized backend.  The state fingerprints and all physical controls were
held fixed.

### Fixed integration eta

| selector | NE | h/eta_int | GF−Lehmann Frobenius difference | GF wall |
| --- | ---: | ---: | ---: | ---: |
| d | 1001 | 1.4890 | `2.99035e-1` | `0.205 s` |
| d | 2001 | 0.7445 | `3.77322e-1` | `0.405 s` |
| d | 4001 | 0.3723 | `3.18017e-1` | `0.803 s` |
| spd | 1001 | 1.4890 | `2.12807e0` | `0.207 s` |
| spd | 2001 | 0.7445 | `1.39168e0` | `0.404 s` |
| spd | 4001 | 0.3723 | `1.16484e0` | `0.803 s` |

### Controlled integration eta

| selector | eta_int | NE | h/eta_int | GF−Lehmann Frobenius difference |
| --- | ---: | ---: | ---: | ---: |
| d | `0.008` | 933 | 0.3994 | `1.35830e0` |
| d | `0.004` | 1863 | 0.3998 | `6.52224e-1` |
| d | `0.002` | 3725 | 0.3998 | `3.20268e-1` |
| spd | `0.008` | 933 | 0.3994 | `3.43760e0` |
| spd | `0.004` | 1863 | 0.3998 | `2.03628e0` |
| spd | `0.002` | 3725 | 0.3998 | `1.17844e0` |

### Energy window

With `eta_int=0.002 Ry` and `NE=4001` fixed:

| selector | margin `Ry` | GF−Lehmann Frobenius difference |
| --- | ---: | ---: |
| d | 0.30 | `3.20680e-1` |
| d | 0.60 | `3.18017e-1` |
| d | 1.00 | `3.13493e-1` |
| spd | 0.30 | `1.17120e0` |
| spd | 0.60 | `1.16484e0` |
| spd | 1.00 | `1.16375e0` |

The fine diagnostic sample reports the two Kubo terms closing internally at
`2.07e-12` Frobenius residual.  Its finite-window spectral residuals are
`1.47e-3` (zeroth), `2.36e-3` (first), and `7.90e-2` (Fermi weighted), for
both endpoint states.  These are the locked DRESP-02R observables and were
not hidden by the optimization.

## 8. Final gate status and return action

`DRESP-02P` is `PASS`: performance is no longer the reason the DRESP-02
material GF/Lehmann closure cannot be tested.  The exact controlled and
window ladders now execute, while their residuals remain macroscopic relative
to the already certified site/product oracle agreements.  Therefore this
remediation does not claim `PASS — numerical closure` and does not silently
reclassify the remaining material result as a new physics rung.

Return immediately to the DRESP-02 material GF↔Lehmann closure gate with the
optimized backend selected and the dense reference retained as its oracle.
The next DRESP-02 decision is whether the completed controlled observables
establish the remaining closure or identify a formulation blocker; no
integration control, tolerance, projection, or physical object is changed by
DRESP-02P.

## 9. Files changed

- `source/lr_projected_reciprocal_chi0.f90`: preserved dense reference and
  added the eigenbasis/vectorized projected GF backend plus timing counters;
- `source/tddft_production_driver.f90`: writes `gf_profile` timing/call-count
  rows beside each DRESP-02R audit sample;
- `tests/unit/test_lr_projected_reciprocal_chi0.f90`: compares full matrices
  and both Kubo terms against the dense reference for every d/spd Γ/finite-q
  fixture;
- `docs/DRESP_02P_PROJECTED_GF_PERFORMANCE.md`: this report.
