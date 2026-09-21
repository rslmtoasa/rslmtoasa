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
separate follow-up gate; the fixed 6401-point audit was exercised manually
and covered GF moments `p=0..4`.
