# DRESP-09X scalar-relativistic physical-spin observable

DRESP-09X is the bounded `ham_only`, collinear, global `L=0` rigid-rotation
audit for the scalar-relativistic physical Pauli spin observable. It is a
radial/endpoint closure only. It does not enter ALSDA/Ward, Dyson, spectra,
BES, Halle, or arbitrary spatial-kernel generalization.

## Two spin quantities that must remain separate

The SCF channel polarization is the scalar-relativistic probability
difference. Its lower component has the ordinary positive probability weight:

```text
SCF channel polarization:
    scalar-relativistic probability density
    lower factor +1
    -> n_up - n_down  (SR2)
```

The transverse response is a physical Pauli-spin matrix element. After the
Kölling-Harmon angular average, the lower component carries `-1/3`:

```text
physical SR spin:
    Pauli spin operator <sigma>
    lower factor -1/3
    -> <sigma_z>       (SRspin2)
```

Therefore SR2 is not the physical-spin target of the transverse DRESP-09S
operator. The distinction is intentional: the SCF/XC radial field may use
`n_up-n_down`, while a transverse susceptibility measures a Pauli-spin
observable. DRESP-09X quantifies this seam and does not change the ALSDA
kernel or its denominator.

The Kölling-Harmon identity used by the oracle is

```text
(sigma . n) sigma_i (sigma . n) = 2 n_i (sigma . n) - sigma_i
average over n and i=x,y,z -> -sigma_i/3
```

The implementation tests all three Cartesian axes independently. The lower
factor is applied to both the small-component radial term and the
`l(l+1)/(TMC*r)^2` angular term. The zero-radius mesh point is a diagnostic
origin and is excluded from volume residuals; its angular term is set to zero.

## Complete second-order endpoint algebra

The accepted coefficient-space objects are

```text
M0 = rho
M1 = H rho
M2 = H^2 rho = rho H^2  (at the accepted equilibrium state)
```

The second-order endpoint response is explicitly ordered as:

```text
second-order endpoint response:
    00  10  01  11  20  02
    -> complete radial spin observable
```

The six branches are the product-rule derivatives of
`D_pq = H^p rho H^q`:

```text
d00 = d rho
d10 = dH rho + H d rho
d01 = d rho H + rho dH
d11 = dH rho H + H d rho H + H rho dH
d20 = dH H rho + H dH rho + H^2 d rho
d02 = d rho H^2 + rho dH H + rho H dH
```

The radial coefficients use the absolute-energy endpoint polynomial with
`phi`, `phidot`, and `phiddot`. The pure `20` and `02` branches contain the
explicit `1/2 phi*phiddot` terms. The equilibrium six-branch contractions
reproduce the production Pauli P3, SR2, and physical-spin SRspin2
polynomials. The historical four-branch routines remain unchanged and are
kept as regression/oracle paths.

## Response hierarchy and current classification

The Fe artifact reports:

```text
old four-branch result
six-branch fixed-observable result
physical-spin target SRspin2
```

It also reports the fixed-observable response by radial `s/p/d` channel and
the independent occupied-state SRspin3 oracle. A response mismatch against
SR2 is retained as the explicit `SCF_CHANNEL_POLARIZATION_vs_PHYSICAL_SPIN`
comparison; it is not silently used as evidence against SRspin2.

The fixed-observable coefficient term is only one term of a rigid rotation:

```text
delta m = Tr[delta rho O] + Tr[rho delta O]
```

The first term is implemented by the six-branch endpoint response. The second
term is the radial-observable/augmentation-frame tangent. DRESP-09X does not
define that missing term by subtracting the residual. The current Fe result is
therefore classified as:

```text
SECOND_ORDER_ENDPOINT_CLOSED_AUGMENTATION_ROTATION_REQUIRED
```

This is PASS-B: the six-branch algebra, target physics, endpoint tangents,
and independent physical-space rotation oracle are closed, while an analytic
augmentation-frame `delta O` remains the sole next milestone. PASS-A requires
adding that independently derived term and closing the total rigid response
against SRspin2.

## Scope and frozen paths

The accepted state is the converged Fe `ham_only` k-space SCF handoff. The
bridge requires the accepted reciprocal eigensystem, second-order Hamiltonian,
complete scalar-relativistic radial channels, and a global `L=0` rotation.

The following are deliberately not touched:

```text
source/exchange.f90
source/lr_lmto_turek_gf.f90
source/lr_lmto_turek_contour.f90
```

DRESP-09U representation-tangent behavior is frozen. No generated Fe output
files are part of the source change.

## Tests and artifact

`UnitDresp09xSrSpinObservable` covers:

1. the Pauli angular-average oracle for x, y, and z;
2. the `-1/3` physical-spin lower factor;
3. the second-order Taylor-product and explicit `20/02` fixtures;
4. Pauli, SR probability, and physical-spin six-branch equilibrium targets;
5. all six commutator tangents and a finite-angle endpoint oracle; and
6. an independent physical-space rigid-spin rotation.

`TddftDresp09XSrSpinObservable` runs the accepted Fe SCF handoff and writes
`/tmp/dresp09x_fe_4k.dat`. `Dresp09XFeArtifact` checks the KH axes, target
closures, six endpoint residuals, finite-angle ratios, PASS-B classification,
and the explicit `ALSDA_Ward = NOT RUN` / `BES_Halle = OFF` scope markers.
