# DRESP-02 — Projected Reciprocal Bare Susceptibility

Status: implemented on `fable_v4` at the DRESP-01 baseline revision
`68cdc10d8f6f8fb8d561887fab401c07f4854bc8`.

The algebraic certification is PASS.  The material real-axis GF closure is
not certified, so the overall result is:

`ALGEBRAIC PASS — MATERIAL GF CLOSURE BLOCKED`

The remediation campaign and its material exit classification are recorded in
[`DRESP_02R_MATERIAL_GF_CLOSURE.md`](DRESP_02R_MATERIAL_GF_CLOSURE.md).

No interaction kernel, Goldstone correction, Dyson solve, loss matrix, or mode
fitting is part of this deliverable.

## Scope and representation

DRESP-02 evaluates the bare transverse circular response in the site basis

\[
  \chi^0_{ij}(\mathbf q,\omega),\qquad i,j=1,\ldots,N_{site},
\]

for the certified orthogonal, collinear, `ham_only`, second-order/HOH
reciprocal baseline.  The accepted endpoint is the exact folded `k+q`
state supplied by the reciprocal lifecycle.  DRESP-02 never adds a site
phase or rephases an endpoint eigenvector.

The supported selectors are the DRESP-01 selectors:

| selector | one-electron orbital content | response cutoff |
| --- | --- | --- |
| `d` | `l=2` only | complete `L=0..4` |
| `spd` | `l=0,1,2` | complete `L=0..4` |

For one site the complete product representation has the usual 232 compact
coordinates.  Production DRESP-02 retains only the `N_site x N_site`
susceptibility and the site-projected GF vertices; it does not allocate a
232-by-232 response matrix or a point-space response.

The implementation is `source/lr_projected_reciprocal_chi0.f90`, using the
unchanged DRESP-01 operator contract in
`source/lr_projected_site_spin.f90`.

## Lehmann backend

For each accepted pair of reciprocal endpoint eigenstates, DRESP-01 produces
the site transition vector

\[
  T_i^{nm}(\mathbf k,\mathbf q)=p_i^H z_{+}^{nm}(\mathbf k,\mathbf q).
\]

The production accumulator is directly

\[
 \chi^{0,+}_{ij}(\mathbf q,\omega)=
 \frac{2}{\sum_k w_k}\sum_{\mathbf k}w_k
 \sum_{nm}
 \frac{f_{n\mathbf k}-f_{m,\mathbf k+q}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+q}+i\eta_{response}}
 T_i^{nm}T_j^{nm*}.
\]

The minus channel uses the existing DRESP/product minus channel and is used
for the q/−q covariance audit.  Occupations are the immutable reciprocal
Fermi snapshot, not a separately solved occupation set.

The implementation path is therefore

`accepted eigenpairs -> DRESP-01 transition amplitudes -> site outer product`.

The compact product Lehmann service is not called by this production path.
It is used only in the unit test as an independent oracle after projection.

## Independent real-axis GF backend

The GF backend uses the existing weighted reciprocal resolvent builder and
the same DRESP-01 site functional.  The site vertex is formed from the
existing product component vertices as

\[
 V_{i,c}=\sum_a p^*_{ia}V_{a,c},
 \qquad c=1+p+2q,
\]

where `p` and `q` are the left and right affine endpoint powers.  No energy
denominator or response accumulation is embedded in the vertex.

With

\[
 A(E)=\frac{i}{2\pi}\left[G^R(E)-G^A(E)\right],
\]

the real-axis implementation evaluates the two existing Kubo terms with the
locked retarded convention:

\[
 \chi^{0,+}_{ij}(\omega)=\frac{2}{\sum_k w_k}\sum_k w_k
 \int dE\,f(E)\,
 \left\{
 \operatorname{Tr}[A_k(E)V_iG^R_{k+q}(E+\omega)V_j^\dagger]
 +\operatorname{Tr}[A_{k+q}(E)V_j^\dagger G^A_k(E-\omega)V_i]
 \right\}.
\]

The two backends have independent accumulators.  The GF routine never calls
the Lehmann routine and never consumes a compact susceptibility.

`eta_response` is the physical retarded pole width.  `integration_eta` is
the independent real-axis resolvent broadening and is required to satisfy

`0 <= integration_eta < eta_response`.

The GF energy interval is derived from both endpoint spectra and expanded by
the explicit `gf_energy_margin`.  Simpson quadrature requires an odd number
of points.

## Public API

The API is deliberately limited to a q-local bare site response:

```fortran
type(projected_chi0_request) :: request
type(projected_chi0_result) :: result

request%q = q
request%frequencies = omega
request%eta = eta_response
request%integration_points = 2001
request%integration_eta = integration_eta
request%energy_margin = energy_margin
request%channel = lr_channel_plus
request%contract => dresp_contract
request%product_basis => product_plus
request%electronic_state => accepted_state
request%q_endpoint_state => accepted_kq_state

call evaluate_projected_lehmann_chi0(request, result)
call evaluate_projected_gf_chi0(request, result)
```

The result records selector, channel, product dimension, q, frequency grid,
both eta values, energy interval, transition count, point-response flag,
vertex/susceptibility memory, timing, provenance, and the complex complete
`site x site x frequency` matrix.

## Unit certification

`tests/unit/test_lr_projected_reciprocal_chi0.f90` is registered as
`UnitLrProjectedReciprocalChi0`.  It uses numerically diagonalized complex
fixtures and checks:

- one-site and mandatory two-site cases;
- `d` and `spd` selectors;
- q=0 and finite q;
- Lehmann against an independently evaluated compact product Lehmann oracle;
- GF against an independently evaluated factorized compact product GF oracle;
- the complete two-site `chi11`, `chi12`, `chi21`, and `chi22` entries, all
  finite and nonzero;
- positive-frequency retarded pole sign;
- q/−q covariance for both selectors and both direct backends;
- a GF resolution/broadening ladder with decreasing GF-to-Lehmann residual;
- no production point-response allocation.

Observed maximum compact-oracle residuals were:

| fixture/selector | Lehmann oracle | GF product oracle |
| --- | ---: | ---: |
| one-site `d` | `2.35e-38` | `5.90e-38` |
| one-site `spd` | `3.18e-22` | `1.97e-18` |
| two-site `d` | `2.94e-38` | `1.03e-38` |
| two-site `spd` | `1.59e-22` | `9.91e-18` |

The GF-to-Lehmann Γ ladder decreased from `1.64e-23` to `9.46e-24`
for one-site `d`, from `1.26e-4` to `7.68e-5` for one-site `spd`,
from `1.51e-23` to `8.71e-24` for two-site `d`, and from `1.24e-4` to
`7.60e-5` for two-site `spd`.  The largest recorded q/−q Lehmann residual
was `7.18e-10` and the largest GF covariance residual was `1.39e-9`; both
are below the `2e-8` fixture tolerance.

The two-site Γ matrix is explicitly exercised at the second test frequency.
For example, the `spd` entries are approximately

```text
chi11 = -2.9205359425e-07 - 5.3123651813e-09 i
chi12 = -6.8837444758e-11 + 3.4614093855e-10 i
chi21 = -5.4628757547e-11 - 3.5744088222e-10 i
chi22 = -1.5733468156e-12 - 9.3554683085e-14 i
```

## Material gate

The driver accepts `backend = 'projected_chi0'`, reuses the accepted
reciprocal k-space SCF cache, and exits before any interaction or Dyson
object.  It writes a q-local site matrix artifact with the Lehmann/GF
entries and closure diagnostics.

The material execution used the accepted bcc-Fe `ham_only`, second-order
state on a 4x4x4 full k mesh, 64 accepted k points, 300 K, and
`eta_response = 0.01 Ry`.  The accepted state reported eight electrons and
an SCF residual of `1.8e-9`.  The production material gate reports the
accepted reciprocal band moments (the same m0 contract used by the SCF)
and checks the `spd` sum against the accepted state moment:

| selector | projected moment (`mu_B`) |
| --- | ---: |
| `d` | `1.970011` |
| `spd` | `1.950787` |

The matching accepted-state total was `1.950787628 mu_B`; the reported
`spd`/total residual was below `2e-5`.  The direct DRESP operator moment is
also written as a diagnostic in the artifact (`d=2.020544770`,
`spd=2.001437103`) and remains finite; it is not substituted for the
accepted SCF band-moment gate.

The same accepted state was also evaluated at Γ and the mesh-compatible
`q=(+/-0.25,0,0)` points, at `omega=0` and `0.015 Ry`.  The final scratch
artifact used Simpson quadrature with 101 energy points,
`integration_eta=0.002 Ry`, `gf_energy_margin=0.60 Ry`, and the derived GF
window `[-1.3638307295, 1.6141834982] Ry`.  The following rows show the
`chi_11` complex comparison; the `-0.25` rows agree with `+0.25` under the
q/−q check below.

| selector | q | omega (Ry) | Lehmann (Re, Im) | GF (Re, Im) | abs diff | rel diff |
| --- | ---: | ---: | --- | --- | ---: | ---: |
| `d` | `0` | `0` | `(-26.0971,-1.54443)` | `(-66.7336,-4.31107)` | `40.7305` | `1.59332` |
| `d` | `0` | `0.015` | `(-28.6335,-1.88439)` | `(-74.9168,-6.81628)` | `46.5454` | `1.59332` |
| `d` | `+0.25` | `0` | `(-24.7678,-0.930873)` | `(-46.1301,15.8644)` | `27.1741` | `1.74231` |
| `d` | `+0.25` | `0.015` | `(-26.4384,-1.31332)` | `(-79.6033,19.3510)` | `57.0395` | `1.74231` |
| `spd` | `0` | `0` | `(-30.5844,-3.50571)` | `(-127.626,-22.7752)` | `98.9367` | `3.59027` |
| `spd` | `0` | `0.015` | `(-35.4206,-8.77265)` | `(-161.622,-69.3080)` | `139.969` | `3.59027` |
| `spd` | `+0.25` | `0` | `(-26.4107,-0.967230)` | `(-51.6788,14.4489)` | `29.5996` | `1.66315` |
| `spd` | `+0.25` | `0.015` | `(-28.1447,-1.37691)` | `(-84.1858,9.18639)` | `57.0279` | `1.66315` |

At 1001, 2001, and 4001 GF points the q=0 site-matrix closure remained
unresolved/non-monotone at the material scale.  Representative Frobenius
differences were:

| GF points | `d` difference | `spd` difference |
| ---: | ---: | ---: |
| 1001 | `2.99e-1` | `2.13e0` |
| 2001 | `3.77e-1` | `1.39e0` |
| 4001 | `3.18e-1` | `1.16e0` |

The material state and both direct backends are available, but the GF ladder
does not establish closure at the material tolerance.  This is recorded as
the material GF closure blocker; it is not converted into an algebraic
pass by changing the tolerance or omitting the GF backend.

The pre-existing DRESP-01 accepted 8x8x8 Fe artifact remains a separate
ground-state moment reference.  A 101-point spot on that state was also
finite but under-resolved (`d` difference `2.51e1`, `spd` difference
`3.00e1`), so it is not used to claim material GF closure.

## Performance and storage

The projected production result retains only the site matrix and site-GF
vertices; it does not allocate the compact `232 x 232` susceptibility or a
point-response array.  In the final local unit run, the complete projected
fixture target took `14.25 s` with a maximum process resident set of about
`39,988 KB`.  The separate compact product-GF oracle target took `12.04 s`
and about `16,188 KB`; these are process-level unit observations, not a
cross-machine production benchmark.  These values are the pre-DRESP-02P
baseline.  The optimized eigenbasis GF backend and its measured Fe ladder are
reported in [`DRESP_02P_PROJECTED_GF_PERFORMANCE.md`](DRESP_02P_PROJECTED_GF_PERFORMANCE.md);
the result metadata also exposes the reference/optimized stage timings and
call counts for downstream measurements.

## Files changed for DRESP-02

- `source/lr_projected_reciprocal_chi0.f90`: direct site Lehmann and
  independent site GF services;
- `source/lr_projected_site_spin.f90`: DRESP-01 site-component vertex view;
- `source/lr_lmto_product_response_basis.f90`: optional selector-aware
  component/transition helpers, with the absent-selector route unchanged;
- `source/lr_ks_susceptibility.f90` and
  `source/lr_product_gf_susceptibility.f90`: selector-aware compact oracle
  inputs;
- `source/tddft_production_driver.f90`: material backend seam and output;
- `tests/unit/test_lr_projected_reciprocal_chi0.f90`: unit certification;
- `CMakeLists.txt` and `source/CMakeLists.txt`: target/source registration.

No existing DRESP-01 operator formula, selector meaning, core policy, or
accepted-state provenance was changed.

## Verification commands

```text
cmake --build build --target UnitLrProjectedReciprocalChi0 rslmto.x -j2
ctest --test-dir build -R '^UnitLrProjectedReciprocalChi0$' --output-on-failure
```

The focused unit test passes.  Existing product-Lehmann, product-GF, and
production-driver tests also passed in the focused regression run.

## Next API boundary

The next consumer may use `projected_chi0_result%susceptibility` as the
site-space bare input to a separately certified interaction/Dyson layer.
That layer is intentionally not started by DRESP-02.
