# DRESP-10A: six-branch LMTO product basis

## Scope

DRESP-10A promotes the live compact LMTO product basis from the historical
four affine endpoint branches to the certified second-order branch set

```text
00, 10, 01, 11, 20, 02
```

The scope is the LMTO radial/product representation, compact bare Lehmann
`chi0`, reciprocal-GF `chi0`, compact vertex contraction, legacy-subspace
embedding, and representation/memory metadata. Full ALSDA Ward closure, BES,
and Halle are deliberately outside this milestone and remain unclaimed.

The frozen files `source/exchange.f90`, `source/lr_lmto_turek_gf.f90`, and
`source/lr_lmto_turek_contour.f90` are not part of this migration.

## Authoritative contract

`source/lr_lmto_endpoint_branches.f90` is the single branch map and endpoint
power contract used by the compact product and reciprocal-GF paths:

| branch | left power `p` | right power `q` | radial coefficient |
|---|---:|---:|---|
| `00` | 0 | 0 | constant term plus `enu` shifts and second-order corrections |
| `10` | 1 | 0 | left first-order coefficient |
| `01` | 0 | 1 | right first-order coefficient |
| `11` | 1 | 1 | first-order cross coefficient |
| `20` | 2 | 0 | `0.5 * phiddot_left * phi_right` |
| `02` | 0 | 2 | `0.5 * phi_left * phiddot_right` |

For `dL = E_L - enu_L` and `dR = E_R - enu_R`, the positive-radius
polynomial is

```text
phiL*phiR + dL*dotL*phiR + phiL*dR*dotR + dL*dR*dotL*dotR
  + 0.5*dL**2*ddotL*phiR + 0.5*phiL*dR**2*ddotR.
```

The origin uses the established `s-s` `r**2` extrapolation and zero for
non-`s-s` channels. Endpoint powers are supported through two; the GF
resolvent/moment contract is supported through power four.

## Candidate inventory

The ordered candidate inventory is six branches for every allowed orbital pair:

| orbital selector | per-response-`L` candidate counts | total for `response_lmax=2` |
|---|---|---:|
| `s` | `6` | `6` |
| `sp` | `12, 12, 6` | `78` |
| `spd` | `18, 24, 24, 12, 6` | `348` |

The historical four-branch `spd` inventory is `232`. It is retained only as a
legacy subspace/shadow inventory. The live accepted `spd`, one-site,
`response_lmax=4` compact basis is expected to expose `unpruned_dimension=348`
and, for the accepted full-rank state, `product_dimension=348` in both circular
channels.

## Compact and GF paths

The compact product basis stores branch metadata (`second`, six branches, and
the ordered labels) and uses the same map in candidate coefficients,
component vertex tensors, and transition coordinates. The reciprocal-GF path
uses five moment slots (`0..4`) and supports all three contraction backends:

```text
scalar == optimized == factorized
```

The result metadata records radial order, branch count/labels, maximum GF
energy moment, and component/GF/susceptibility memory bytes. No point-space
transition allocation is required by the compact GF evaluator.

## Memory accounting

For the accepted one-site Fe `spd` coefficient space (`nbasis=18`, complex
double storage), the representation-only allocation comparison is:

| allocation | historical S4 | live S6 |
|---|---:|---:|
| susceptibility, per frequency | 861,184 B | 1,937,664 B |
| component vertices | 4,810,752 B | 10,824,192 B |
| six GF matrices, moments 0..2 vs 0..4 | 93,312 B | 155,520 B |

The susceptibility and vertex changes include both the dimension increase
(`232 -> 348`) and the component count (`4 -> 6`); the GF workspace grows from
three to five stored moments. Runtime result metadata reports the exact byte
counts for the active request.

## Legacy embedding and arbitrary-L audit

The compact-span audit keeps independent four-branch and six-branch weighted
SVDs. It reports principal overlaps, outside-span residuals, branch-20/02
projection residuals, candidate residuals, and direct arbitrary-`(L,M)` physical
source projections. The expected legacy relationship is that the old 232-mode
span is contained in the six-branch span, while the six-branch span can contain
additional modes. A production six-branch SVD is checked independently against
the direct six-branch radial polynomial and Gaunt assembly.

## Required gates

The focused gates are:

```text
UnitLrLmtoProductResponseBasis
UnitLrLmtoProductResponse
UnitLrProductKsSusceptibility
UnitLrProductGfSusceptibility
UnitLrProductGfQuadratureAudit
TddftDresp09ZSCompactSpan
Dresp09ZSFeArtifact
TddftDresp10AGfSpot
```

The unit fixtures set nonzero `phiddot` values so the `20` and `02` branches
are exercised. Expected independent checks include direct six-branch radial
polynomial values, compact-vs-point projection, scalar/optimized/factorized
GF equality, q/omega/eta/channel coverage, immutable-state checks, and the
legacy embedding diagnostics. A PASS-A classification is valid only when all
required gates and artifact checks pass; otherwise the campaign status is
`PARTIAL`, `BLOCKED`, or `REJECTED` with the failing metric named.

## Current evidence

The local focused unit evidence after the migration is:

| check | maximum relative residual |
|---|---:|
| compact GF component/R1 oracle | `2.43e-16` |
| GF transition factorization | `5.87e-20` |
| point-GF projection | `1.35e-15` |
| KS projection | `1.18e-15` |
| KS independent pair accumulation | `3.73e-16` |
| KS immutable-state mutation | `0` |

These numbers are unit-fixture evidence only. The accepted-Fe compact-span
classification and final campaign status must be taken from the generated
artifact and the final gate run, not inferred from these unit values.

The bounded accepted-Fe GF spot uses the same accepted 4x4x4 state and the
full `spd` product space: `Nprod=348`, 21 Simpson points, finite response,
Lehmann norm `3.8629`, GF norm `3.4957`, and relative Frobenius residual
`3.8038e-1`. This is a bounded execution/metadata check, not a quadrature
convergence claim. The long CMake-disabled quadrature campaign remains a
separate follow-up gate; its DRESP-10A-Q closure is recorded below. The fixed
6401-point audit was exercised manually and covered GF moments `p=0..4`.

## DRESP-10A-Q — six-branch reciprocal-GF quadrature gate

**Status: PASS for the independent mixed-eigenvector numerical oracle.**
The audit remains CMake-disabled and is run directly as a long-form numerical
campaign. No ALSDA/Ward, Dyson, spectra, BES, or Halle tests were run.

The mixed fixture is intentionally smaller than the accepted Fe production
space: it has `nbasis=8`, two k points, and `product_dimension=36`. The
accepted Fe representation remains `product_dimension=348`; the mixed fixture
is the independent dense-complex numerical oracle, not a replacement Fe
material fixture.

Before interpreting the quadrature, the rebuilt audit reported:

```text
product dimension       36 (mixed oracle); 348 (accepted Fe production)
branch count            6
branch labels           00,10,01,11,20,02
maximum GF moment       4
```

The transition-factorization oracle closed at:

| channel | maximum absolute residual | maximum relative residual |
|---|---:|---:|
| `chi_plus` | `6.56e-16` | `3.05e-16` |
| `chi_minus` | `1.52e-16` | `2.50e-16` |

The selected mixed eigenvector pairs have nonzero second-order transition
activity. The largest reported branch norms were `1.85e-2` (`20`) and
`8.48e-2` (`02`) in `chi_plus`, and `6.75e-3` (`20`) and `3.34e-3` (`02`)
in `chi_minus`. Thus the high-order branches are not present only as dead
metadata.

The independent weighted-resolvent oracle compared three complex energies for
every mixed-fixture k point and every power `p=0..4` against
`build_weighted_resolvent`; all five maximum relative residuals were zero to
the displayed precision. Independent unweighted moment checks against
Hamiltonian powers closed at `1.05e-15` or below. The fixture has nonzero
off-diagonal and imaginary contributions at `p=3` and `p=4`:

```text
p=3  max |M_ij|=1.10e-1  max |Im M_ij|=9.39e-3
p=4  max |M_ij|=9.55e-2  max |Im M_ij|=5.15e-3
```

The prescribed mixed eta ladder used `energy_margin=1.0 Ry` and remained
resolved at `h/eta=0.435981` for every point:

| `eta` (Ry) | `N_E` | `h` (Ry) | `h/eta` | `R_scalar` (`omega=0`) | `R_optimized` | `R_factorized` | `R_scalar` (`omega=0.17`) | `R_optimized` | `R_factorized` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0100 | 801 | 4.3598e-3 | 0.435981 | 6.226e-2 | 6.226e-2 | 6.226e-2 | 1.674e-1 | 1.674e-1 | 1.674e-1 |
| 0.0050 | 1601 | 2.1799e-3 | 0.435981 | 3.305e-2 | 3.305e-2 | 3.305e-2 | 9.293e-2 | 9.293e-2 | 9.293e-2 |
| 0.0025 | 3201 | 1.0900e-3 | 0.435981 | 1.703e-2 | 1.703e-2 | 1.703e-2 | 4.907e-2 | 4.907e-2 | 4.907e-2 |

The response residual decreases consistently with refinement for both
frequencies. At every ladder point the three GF backends agreed with each
other at approximately `1e-14` relative or better; a backend mismatch would
have been classified separately as `GF_CONTRACTION_BACKEND_MISMATCH`.
The per-element diagnostics were emitted for all three backends with a
nonzero reference-amplitude floor and worst-element indices.

For the required finite-q sample, `q=(0.23,0,0)`, `chi_minus`, `eta=0.0025`
Ry, and `N_E=3201` gave `h/eta=0.435981`,
`R_scalar=R_optimized=R_factorized=1.703e-2`, with pairwise backend
relative differences no larger than `1.12e-14`. The existing rotation oracle
also passed: the maximum eigenpair residual was `4.59e-16`, covariance
residuals were `1.38e-15`/`1.47e-15` for GF and `1.51e-15`/`1.43e-15` for
Lehmann, and the historical production-comparison thresholds were retained.

The direct commands were:

```text
build/bin/UnitLrProductGfQuadratureAudit mixed_one 1
build/bin/UnitLrProductGfQuadratureAudit mixed_one 2
build/bin/UnitLrProductGfQuadratureAudit mixed_one 3
build/bin/UnitLrProductGfQuadratureAudit mixed_q
build/bin/UnitLrProductGfQuadratureAudit mixed_rotation
```

The optional resolved accepted-Fe static Gamma sample was not run. The
independent mixed oracle closes the requested gate, while the existing Fe
21-point execution/metadata spot is retained unchanged and is not treated as
a convergence target.
