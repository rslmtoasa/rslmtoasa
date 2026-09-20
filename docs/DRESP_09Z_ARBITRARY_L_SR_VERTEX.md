# DRESP-09Z arbitrary-L scalar-relativistic response vertex

This document records the implementation boundary reached by DRESP-09Z.  It
does not enter ALSDA, Ward, Dyson, loss spectra, magnon dispersion, BES, or
Halle paths.

## Representation and convention

The live harmonic is the `hcpx` convention already owned by
`response_angular_basis.f90`:

```text
Y_live(l,m) = conjugate(Y_Condon--Shortley(l,m))
```

The production coefficient is therefore the existing covariant Gaunt

```text
G(lm, l'm'; LM) = integral Y_live(LM)* Y_live(lm)* Y_live(l'm') dOmega
```

`(L,M)` is always the physical source/density rank.  `(K,Q)` is reserved for
the orbital-product rank after lower-component recoupling.

The implementation is in:

```text
source/lr_sr_angular_vertex.f90
source/lr_sr_spatial_augmentation_tangent.f90
source/lr_dresp09z_bridge.f90
```

The angular production path uses analytic Gaunt recoupling.  The independent
oracle uses a 24-point Gauss--Legendre quadrature in `cos(theta)` and a
96-point periodic azimuthal quadrature; it does not call the analytic path.

## Angular algebra

For Cartesian component `i`, the lower operator is assembled as

```text
(sigma.n) sigma_i (sigma.n)
       = -sigma_i/3 + 2 sum(j,q) c(i,j,q) Y_live(2,q) sigma_j
```

The rank-2 part is recoupled in the explicit order

```text
(L,M) x (2,q) -> (K,Q)
```

and then contracted with the orbital Gaunt algebra.  No universal `-1/3`
factor is applied to the rank-2 contribution.  Terms with `K > 4` are passed
through the Gaunt selection rules and vanish after contraction with the spd
orbital product space; they are not hand-truncated in the analytic sum.

The lower rank map before the orbital cutoff is:

| source L | lower rank-2 product K retained for spd |
|---:|:---|
| 0 | 2 |
| 1 | 1, 3 |
| 2 | 0, 2, 4 |
| 3 | 1, 3 (`K=5` is forbidden by the spd product space) |
| 4 | 2, 4 (`K=6` is forbidden by the spd product space) |

The upper component remains the ordinary authoritative Gaunt coefficient with
`K=L`.

Cartesian matrices are assembled first.  Circular matrices use
`sigma_plus=(sigma_x+i sigma_y)/2` and
`sigma_minus=(sigma_x-i sigma_y)/2`.  With the live harmonic convention the
covariance tested by the unit oracle is

```text
O_minus(L,-M) = (-1)^M conjugate-transpose(O_plus(L,M))
```

## Six radial branches and augmentation tangent

The arbitrary spatial layer retains all six certified endpoint branches:

```text
00  10  01  11  20  02
```

For each branch it separately returns upper, lower-small, lower rank-0, lower
rank-2, and total matrices.  The augmentation tangent differentiates only the
spin-dependent radial matrices with

```text
delta A = -i [G,A]
```

The laboratory `Y_live(L,M)` coefficient is never rotated.  The finite-angle
oracle rotates the radial matrices and compares against the analytic tangent.

## Unit-oracle result

`UnitDresp09ZAngularVertex` exhaustively covers `l,l' <= 2`, all magnetic
indices, all `0 <= L <= 4`, all `M`, and all Cartesian components.

```text
upper analytic vs quadrature max error       9.8532e-15
lower KH analytic vs quadrature max error    7.3275e-15
forbidden upper-channel amplitude            0.0000e+00
circular covariance max error                1.7347e-18
L=4 rank-mixing fixture norm                 5.7797e-02
arbitrary-L augmentation FD max error        1.0910e-09
```

The branch finite-difference errors were:

```text
00  1.0910e-09
10  1.4540e-10
01  1.4540e-10
11  2.1136e-11
20  8.9156e-11
02  8.9156e-11
```

The nontrivial fixtures include `L=1,M=0`, `L=2,M=1`, `L=3,M=-2`, and
`L=4,M=4` in the augmentation finite-angle sweep.  The explicit L=4 fixture
shows nonzero lower rank mixing.

## DRESP-09ZR scalar L=0 projection

The direct arbitrary-orbital four-harmonic KH integral is independently
validated.  The historical comparison was not a scalar projection: it
compared a full orbital-resolved matrix element against the frozen scalar
observable pairwise.  Its measured residual is

```text
2.25675833e-01
```

For a physical `L=0,M=0` source, the accepted projector is applied shell by
shell:

```text
P0[O]_l = I_l/(2*l+1) * sum(m=-l:l) O_(lm,lm)
```

The full lower operator is retained as `K=0 + K=2`.  The independent
quadrature trace oracle confirms that the `K=2` diagonal shell trace vanishes
for `l=0,1,2`; it is zero identically for `l=0` and traceless, rather than
pairwise zero, for `p` and `d` shells.

The DRESP-09ZR unit oracle reports:

```text
scalar projector closure                     1.1102e-16
rank-2 trace closure l=0/1/2                5.5511e-16  5.6066e-15  3.3584e-15
projected upper/small/rank-0/rank-2/total   5.3291e-15  8.3267e-17  6.9389e-18  0  3.5527e-15
projected branches 00/10/01/11/20/02        5.3291e-15  4.4409e-16  4.4409e-16  5.5511e-17  5.5511e-17  5.5511e-17
```

The raw `2.25675833e-01` value is therefore
`EXPECTED_FULL_ORBITAL_vs_SCALAR_PROJECTED_DIFFERENCE`, not a correctness
residual.  The full rank-2 tensor remains present in the arbitrary-L path.

## Compact candidate audit boundary

The audit bridge reports the exact candidate inventories without modifying the
production compact basis:

```text
historical four-branch spd/response_lmax=4 candidates  232
shadow six-branch candidates                           348
new 20/02 candidate slots                              116
```

The old basis is still four-branch by construction.  The accepted-state
compact-span/SVD audit is the next audit boundary; no orthogonal complement is
discarded and no `chi0` or Dyson object is changed by DRESP-09ZR.

## DRESP-09ZR representation-gate verdict

```text
DRESP-09ZR scalar projection gate: CLOSED
Raw arbitrary-L SR vertex: CLOSED
Accepted-state compact-span audit: NOT RUN IN THIS COMMIT
```

The arbitrary-L analytic tensor, independent quadrature oracle, selection
rules, circular covariance, six-branch tangent layer, scalar projector, and
projected six-branch reduction all close.  This commit resolves the invalid
pairwise L=0 gate; it does not claim the final PASS-A/PASS-B compact-basis
classification because the accepted Fe six-branch SVD and field/density
duality audit have not been run here.

```text
ALSDA Ward: NOT RUN
BES/Halle: OFF
```

The next controlled milestone is to reconcile the direct nonspherical KH
orbital tensor with the certified scalar L=0 projection, then run the
accepted-state compact span/SVD audit.  Only after those two gates can a
DRESP-10 product-basis extension decision be made.
