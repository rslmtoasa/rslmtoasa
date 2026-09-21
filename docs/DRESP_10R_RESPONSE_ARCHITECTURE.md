# DRESP-10R response architecture

DRESP-10R is the raw, static, full-spatial ALSDA Ward gate for the accepted
ham-only bcc-Fe state.  It stops before ALSDA interpretation when an upstream
representation gate is open.  It does not run Dyson dynamics, spectra, BES,
Halle, Goldstone repair, fitting, or kernel rescaling.

## Authoritative map

```text
P3 physical Pauli density
      |
      | Kxc = Bxc / P3
      v
SR physical field
      |
      | F_SR: six-branch 00/10/01/11/20/02 insertion
      v
coefficient-space delta H(k)
      |
      | exact static Frechet L_f(H)
      v
delta rho(k)
      |
      | Pauli six-branch observable
      v
fixed Pauli response
      |
      +---- augmentation/contact delta O
      v
total Pauli response
```

The compact coordinates are only an explicit 348-dimensional representation
of the physical radial source.  They are not an implicit conversion step.
Every production six-branch action goes through
`lmto_product_apply_branch_action`; its stored-dual mode applies the endpoint
orientation (H^q V H^p), with both powers looped independently.  The branch
metadata is `00=(0,0)`, `10=(1,0)`, `01=(0,1)`, `11=(1,1)`, `20=(2,0)`, and
`02=(0,2)`.

## Response objects

| object | physical meaning | representation | dimension | metric | units | producer | consumer | independent oracle |
|---|---|---|---:|---|---|---|---|---|
| P3 Pauli density | frozen Pauli spin polarization | physical radial profile | site × radial point | accepted radial volume metric | density | `pauli_density_from_moments` | `Kxc`, Pauli target | accepted Pauli projection / DRESP-09Y |
| SR physical field | full scalar-relativistic transverse source | site × (L,M) × radial point | 348 raw coordinates | response-space radial metric | field | `lr_full_spatial_source_components` | six-branch insertion | DRESP-09ZS physical vertices |
| compact density coordinates | Pauli/SR density response coefficients | live production product basis | 348 | production weighted modes | response coordinates | `lmto_product_response_basis` | `chi0`, static action | DRESP-10A-Q direct Lehmann/GF |
| compact field coordinates | covariant field coefficients | live production product basis | 348 | production weighted modes | field coordinates | product projection | source operator | raw/compact pairing and adjoint tests |
| coefficient-space (delta H) | perturbation of the accepted LMTO Hamiltonian | site-major spin/orbital matrix | 18 × 18 for one Fe site | Frobenius and band-action metrics | energy | six-branch field insertion | Frechet density | native (-i[G,H]) |
| delta rho | static density tangent | coefficient-space matrix | 18 × 18 | Frobenius / occupied action | density matrix | `lr_static_frechet_density` | endpoint observables | commutator (-i[G,rho]) |
| augmentation (delta O) | observable-side radial contact tangent | upper, lower-small, lower-angular(rank-0+rank-2), total | site × radial point | DRESP-09Y volume metric | observable density | `dresp09y_augmentation_density` | total Pauli/SR response | finite-angle DRESP-09Y tangent |
| static (chi_0) | exact finite-temperature compact response | dense product matrix | 348 × 348 | compact coordinate metric | response / energy | `lr_static_product_matrix` | static action gate | direct Frechet action |
| (K_{xc}) | declared ALSDA field kernel | pointwise radial ratio | site × radial point | positive-measure points only | field / density | direct `Bxc/P3` | ALSDA field insertion | pointwise (K_{xc}P3=Bxc) |
| mixed (K) | arbitrary-((L,M)) field-to-response operator | full source/contact composition | not assembled | fail-closed | operator | future `F_SR Kxc D_P` | static denominator | deterministic all-(L,M) action tests |
| contact (C) | augmentation/contact operator | full arbitrary-((L,M)) map | not assembled | fail-closed | response | future contact map | static denominator | L=0..4 contact-map actions |

The production Fe run currently reports `SR_SOURCE_SPACE_INCOMPLETE`: the
source span is tested using

\[
R = \frac{\lVert (I-P_{348})V_{\rm physical}\rVert}
              {\lVert V_{\rm physical}\rVert},
\]

for Pauli, SR upper, lower-small, lower rank-0, lower rank-2, SR total, and
augmentation delta-O at every (L=0\ldots4).  A nonzero generated operator is
not accepted as a completeness proof.

## Circular and Cartesian conventions

The production convention is

\[
V_+ = (V_x-iV_y)/2, \qquad V_-=(V_x+iV_y)/2,
\]

so (V_x=V_++V_-) and (V_y=i(V_+-V_-)).  This is tested against the native
SU(2) generator (G=\sigma_y/2), for which
\(-i[G,\sigma_z]=\sigma_x).  The six-branch source and observable paths use
the same circular convention, including the stored-dual branch swap.

## Gate order and current result

The gates are ordered as follows:

1. six-branch matrix action, including explicit 20/02 (H^2) tests;
2. exact static Fréchet matrix/direct action;
3. live 348-space source-span projection;
4. full SR field insertion versus the native tangent;
5. full DRESP-09Y contact and Pauli closure;
6. only then ALSDA field and raw Ward interpretation;
7. only after rigid closure, arbitrary-((L,M)) mixed/contact operators.

On the accepted Fe smoke state, the endpoint and static gates close at machine
precision.  The contact hierarchy is also reproduced: fixed observable versus
SR spin target is about `3.0095e-1`, while complete SR response is
`5.28e-15` and complete Pauli response versus P3 is `5.30e-15`.  The live
source-span gate remains open, with the largest physical SR residual about
`7.8010e-1` and the largest delta-O residual about `5.9902e-1`, both at
(L=4).  Therefore no ALSDA Ward defect is interpreted.

## Historical comparison paths

The following remain comparison oracles and are not promoted to production:

- DRESP-06A simple (U^H K U) Ward route;
- DRESP-09S scalar (L=0) source;
- DRESP-09Y rigid augmentation-frame response;
- the historical 232-mode four-branch basis.

The pre-consolidation DRESP-10 prototype residuals are preserved as historical
integration diagnostics, not ALSDA physics:

```text
field insertion = 7.870658e-1
static Frechet  = 8.282792e-2
native+contact = 9.587104e-1
ALSDA+contact  = 9.925727e-1
status         = SUPERSEDED INTEGRATION PROTOTYPE
```

`BES/Halle = OFF`, Goldstone correction is off, dynamic spectra are not run,
and the raw denominator remains unassembled while the source-span gate is
blocked.
