# DRESP-09W radial observable / ground-state provenance

DRESP-09W audits the radial seam after the DRESP-09V coefficient-space
`M0/M1/M2` and tangent closures. It is restricted to the accepted Fe
`ham_only`/HOH state, a zero-frequency global `L=0` rotation, and the accepted
radial mesh. It does not enter ALSDA/Ward, Dyson, spectra, BES, Halle, or the
DRESP-09U field mapping.

## Like-for-like objects

The production run keeps the Pauli and full scalar-relativistic observables
separate:

| object | valence/core | radial metric and source |
| --- | --- | --- |
| P1 Pauli target | valence | large component only; current `compute_accepted_pauli_magnetization` target |
| P2 Pauli direct | valence | independent occupied-state `phi + DeltaE*phidot`, large component only |
| P3 Pauli moments | valence | production `M0/M1/M2` polynomial, including `phi*phiddot` |
| SR1 accepted total | total | accepted `rho_weighted_up-rho_weighted_down` snapshot |
| SR2 SR moments | valence | live `lmto_density_from_moments` full-SR polynomial |
| SR3 SR direct | valence | independent occupied-state `lmto_state_density` oracle |
| SR4 SR+core total | total | SR2 plus the accepted frozen-core profile |

For every site, angular channel, and spin, the moment reconstruction uses the
accepted `enu_work` convention:

```text
Q0 = M0
Q1 = M1 - Enu_work*M0
Q2 = M2 - 2*Enu_work*M1 + Enu_work**2*M0
```

The Pauli production polynomial is

```text
Q0*phi**2 + 2*Q1*phi*phidot
  + Q2*(phidot**2 + phi*phiddot)
```

The full-SR polynomial delegates to the existing live implementation. Its
`A0/A1/A2` terms include `GFAC`, upper and lower components, and both first and
second energy derivatives. The independent SR occupied-state oracle uses the
same `lmto_state_density` formula state by state; it does not call the moment
reconstruction.

## Pauli target hierarchy and response

P1/P2 are the current accepted direct-state contract and agree in the Fe
artifact. P3 is intentionally not silently substituted for P1: the
production second-order polynomial contains `phi*phiddot`, while P1/P2 use
`|phi + DeltaE*phidot|**2`. DRESP-09W reports P1-P2, P1-P3, and P2-P3 both
totally and by `s/p/d` channel, so this target-definition seam remains
quantified.

The Pauli rigid response uses the exact DRESP-07/DRESP-09R coefficient-space
rotation and the large-component-only endpoint bilinear. It is compared
independently with P1, P2, and P3. A mismatch is an observable/target seam,
not permission to change the accepted Pauli target or the future `Kxc`.

## Full-SR hierarchy and core

SR2/SR3 are a hard direct-state versus moment-polynomial oracle. SR4 is
formed only after adding the separately captured frozen core. The total
accepted snapshot SR1 is never used as a valence target. The artifact also
reports the exact core decomposition identity and:

```text
core integrated spin moment
core weighted L2 norm
core/max-valence radial amplitude
core s-like integrated moment, L2 norm, and maximum amplitude
```

The core may have a tiny integrated spin moment but a large localized radial
amplitude, so it is not discarded based on its integral.

The full-SR rigid response remains the complete DRESP-09V endpoint tangent and
full-SR transverse bilinear. It is compared only with SR2. The historical
`0.305325` value is retained as:

```text
DRESP09V_0p305325 = CROSS-REPRESENTATION DIAGNOSTIC
```

It is the old full-SR-response versus Pauli-target comparison and is not a
radial correctness gate. The new `new_SR_to_SR_valence_residual` is the
representation-matched result.

## Origin and measure

Weighted radial observables are stored as

```text
weighted(r) = 4*pi*r**2*n(r)
```

The origin is handled by the existing radial-origin contract for the accepted
P1 target. For DRESP-09W residual norms and maxima, `r=0` is excluded because
it has zero volume measure; the physical pointwise conversion is performed
only for positive mesh points. Response-space `Y00` coefficients are
converted with:

```text
response coefficient -> sqrt(4*pi)*r**2*coefficient
```

All volume diagnostics use the accepted `r`, logarithmic mesh parameters `a`
and `b`, and the legacy Simpson Jacobian
`Simpson_weight*a*(r+b)`. No fitted normalization is applied.

## Classifications

The backend writes exactly one primary classification:

```text
RADIAL_OBSERVABLE_PROVENANCE_CLOSED
PAULI_TARGET_PROVENANCE_OPEN
PAULI_RESPONSE_OBSERVABLE_OPEN
SR_RADIAL_RECONSTRUCTION_OPEN
ACCEPTED_RADIAL_SNAPSHOT_PROVENANCE_OPEN
CORE_DECOMPOSITION_OPEN
MIXED
```

`PASS-A` requires all target hierarchies, the core split, and both
representation-specific response closures. `PASS-B` is reserved for a
bounded, explicitly localized target/projection seam after one representation
closes. A missing or unreproducible accepted radial snapshot is `BLOCKED`.

The backend writes `ALSDA_Ward = NOT RUN` and `BES_Halle = OFF` in every
artifact.
