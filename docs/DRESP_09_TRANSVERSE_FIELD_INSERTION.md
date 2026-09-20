# DRESP-09 transverse spatial-field insertion

## Starting point

The campaign was resumed at:

```text
branch: fable_v4
HEAD:   226e8add06b87261fec79774283a57638b3be901
commit: Validate native second-order LMTO field mapping
```

The worktree already contained unrelated generated smoke outputs. They were
left untouched.

## Representation result

The compact field algebra is implemented in
`source/lr_dresp09_field_insertion.f90`.

For a raw covariant field `B` and density `D`, the compact coordinates are

```text
b = U^H W^(1/2) B
d = U^H W^(1/2) D
```

and the hard pairing is

```text
b^H d = B^H W D.
```

The projection loops over the complete retained `(site,L,M,radial)` space;
there is no L=0 reduction and no second insertion of `W`. The origin remains
the exact zero-measure LR-04 null point.

The certified product transition tensor uses the live measurement-first
convention

```text
T(mu,n,m) = <n|O_mu|m>.
```

The compact metric-adjoint source is therefore

```text
V = sum_mu conjg(b_mu) O_mu^H,
```

while the physical retarded source form uses the same `O_mu^H` with the
supplied source coefficient un-conjugated. For linearized LMTO endpoint powers,
the adjoint reverses the powers exactly:

```text
(p,q) -> H^q O_mu^H H^p.
```

This is implemented independently of the raw-space projection helper.

## Validation

`UnitDresp09FieldInsertion` is a complex deterministic representation oracle.
It covers all 25 angular coordinates for `L=0..4`, raw/compact round trip,
compact-versus-raw energy pairing, the metric-adjoint trace identity, the
physical source trace identity, and endpoint-power reversal. The focused test
passed.

The DRESP-08 hard oracle remains the material rotation authority and was not
reopened or changed.

## Direct-route boundary

The requested independent direct radial/angular LMTO route cannot be certified
from the current accepted data contract. The production contract exposes:

- large and small scalar-relativistic radial numerators;
- same-spin `GFAC` normalization information in the radial solver;
- no certified mixed-up/down radial inner product for an off-diagonal spin
  operator;
- no certified lower-component spin-angular convention for nonspherical
  response harmonics and cross-`l` terms;
- no complete production snapshot of the small-component energy derivatives
  and the associated `GFAC`/potential provenance needed by the four
  `B**`, `B*dot`, `Bdot*`, and `Bdot*dot` branches.

The existing response mapping contract explicitly identifies this missing
lower-component spin-angular seam. `POTPAR` was not differentiated and no
average or geometric combination of spin-dependent `GFAC` values was inserted.

The direct API therefore fails closed with:

```text
BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED
```

`UnitDresp09DirectRouteBlocked` is the expected-failure regression for that
boundary. This is a capability block, not a numerical residual.

## ALSDA disposition

The raw ALSDA static Ward identity was not reevaluated as a certified DRESP-09
result. Doing so before the direct field vertex is closed would mix the
certified compact adjoint with an uncertified scalar-relativistic radial
metric. BES/Halle and kernel tuning remain off.

The next required input is an accepted augmentation contract that specifies the
mixed-spin large/small radial inner product, lower-component angular mapping,
and complete derivative provenance. Once supplied, Route A can be compared to
the compact adjoint before rerunning the raw ALSDA Ward observable.
