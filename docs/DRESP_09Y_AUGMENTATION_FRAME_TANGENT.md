# DRESP-09Y scalar-relativistic augmentation-frame tangent

DRESP-09Y closes the missing observable-side term in the bounded
`ham_only`, collinear, global `L=0` rigid-rotation audit. It starts from the
accepted Fe reciprocal eigensystem and uses the existing DRESP-09X six-branch
endpoint algebra. It does not re-solve radial equations for rotated
directions and does not enter ALSDA/Ward, Dyson, spectra, BES, or Halle.

## Frame and tangent

For each local radial factor, the collinear diagonal matrix is promoted to a
global coefficient-frame matrix by

```text
A(theta) = U(theta) A_z U(theta)^dagger
U(theta) = exp(-i theta G),   G = sigma_y/2
delta A  = -i [G, A]
```

The observable is assembled before density contraction:

```text
O_sigma = A_L sigma A_R
delta O_sigma = delta A_L sigma A_R + A_L sigma delta A_R
```

The implementation applies this to all six absolute-energy branches
`00, 10, 01, 11, 20, 02`, and keeps the upper, lower-small, and
lower-angular terms separate. The lower angular factor is constructed from
`C_sigma = G_sigma/(TMC_sigma*r)` and retains the common outer `1/r^2` of the
scalar-relativistic radial observable. Thus the angular contribution carries
the same `1/r^4` pointwise normalization as the existing DRESP-09X
contract.

The coefficient-side contraction preserves the established DRESP-09X stored
dual ordering, `sum D(i,j)*O(i,j)`. The complete rigid response is therefore

```text
delta m = density-side six-branch term + augmentation-frame delta O term
        = sum delta D(i,j) O(i,j) + sum D(i,j) delta O(i,j)
```

No residual subtraction is used to construct the augmentation contribution.

## Production bridge and scope

`sr_augmentation_tangent` is a static Gamma-point backend requiring
`response_lmax=4`, the accepted k-space SCF handoff, and the second-order
Hamiltonian. The bridge reports:

- analytic-vs-finite-angle `delta O` errors for upper, lower-small,
  lower-angular, and total terms;
- all six branch errors and the circular Hermitian swap oracle;
- endpoint and `[H,rho]` commutator residuals;
- fixed, augmentation, and complete response norms and integrated values;
- `s/p/d` resolved contributions and the `2 Re` interference term; and
- Pauli hierarchy checks against P1 and P3.

The frozen source paths remain untouched:

```text
source/exchange.f90
source/lr_lmto_turek_gf.f90
source/lr_lmto_turek_contour.f90
```

## Fe result

Starting HEAD:

```text
9228e47f442a49888f362fadc904bc709a409292
```

The accepted Fe artifact is `/tmp/dresp09y_fe_4k.dat`. The production result
is:

```text
verdict = PASS-A
primary_classification = AUGMENTATION_FRAME_TANGENT_CLOSED
```

Representative closures are:

```text
max analytic delta-O branch error       2.6663390e-10
max circular Hermitian error            4.4408921e-16
complete response vs SRspin2            5.2796234e-15
integrated complete-target residual     3.3306691e-15
finite-angle theta^2 ratios             3.9999996, 3.9999988
Pauli complete vs P3                    5.3039529e-15
```

The fixed-observable residual remains reported as a reference quantity:

```text
fixed_observable_vs_SRspin2_relative = 3.0095449441e-01
```

It is closed only after adding the independently derived augmentation-frame
term. ALSDA/Ward is `NOT RUN`; BES/Halle is `OFF`.

## Checks

The focused unit test is:

```text
UnitDresp09YAugmentationTangent
```

It covers arbitrary-axis `delta A`, Cartesian and circular observables,
spin-independent controls, all six radial branches, upper/small/angular
components, circular Hermiticity, and an independent finite-dimensional
stored-dual trace oracle. The Fe integration and artifact validator are:

```text
TddftDresp09YAugmentationTangent
Dresp09YFeArtifact
```
