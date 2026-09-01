# TDDFT-06 K-space GF chi0 evidence

## Implemented contract

`source/tddft_chi0_green.f90` now supplies the validated K-space Green-function
`chi0` backend for the supported `ham_only`, collinear, SOC-free response
boundary.  The dynamic path evaluates

```text
chi_AB^R(q,w) = -1/(2*pi*i) integral dE f(E) Tr[
  A G_{k+q}^R(E+w) B (G_k^R(E)-G_k^A(E))
 +A (G_{k+q}^R(E)-G_{k+q}^A(E)) B G_k^A(E-w)
].
```

The energy integral is a direct real-axis trapezoid with explicit
`green_eta`.  The one-particle half-width defaults to `eta/2`, so the bubble
has the same effective denominator broadening as the explicit eigenpair
reference.  The provider uses `lehmann_kspace_resolvent`, which is the
validated orthogonal-basis `G(k,z)` primitive from TDDFT-05.

The left and right operators are constructed by the same
`site_projected_operator` used by the eigenpair response.  No XC/kernel
object is referenced by the GF bubble.  The source/measurement operators
therefore remain separate from the common Dyson kernel.

The K and K+q endpoint arrays are separate inputs.  The production driver
obtains K+q from `kpoint_workset%shifted(q)`, preserving the exact folded
endpoint while recording the requested direct-basis q.  The returned metadata
records the endpoint branches, q coordinate, energy window, energy-point
count, one-particle `green_eta`, effective response `eta`, Fermi level,
temperature, k mesh, and frequency grid.

The static API is intentionally separate from the dynamic near-real-axis
path.  The Lehmann source provides the exact zero-frequency divided
difference through the same endpoint spectral data.  A future Green-function
source that cannot provide an exact static limit fails explicitly; it is not
silently sampled at finite eta.

## Deterministic fixture evidence

The executable `UnitTddftGreenChiKS` contains:

- an analytic one-site, two-level spin-split fixture;
- a four-k-point periodic fixture at three exact k+q shifts;
- a two-site local-response fixture;
- an exact static `chi0(0,0)` comparison;
- a Ward residual comparison using the same circular response basis; and
- the retarded circular frequency/channel symmetry test.

The command is:

```text
cmake --build build --target UnitTddftGreenChiKS
ctest --test-dir build -R UnitTddftGreenChiKS --output-on-failure
```

The checked-in run passed.  The following pointwise rows are from the
periodic fixture; each row is the complex site `chi0(1,1)` element, not a
plotted spectrum.  The response broadening is `eta=0.001 Ry`, with
`green_eta=0.0005 Ry`, `T=2000 K`, and 16,001 real-axis energy nodes.

| q shift | omega (Ry) | Re GF | Im GF | Re eigenpair | Im eigenpair | absolute difference |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.025 | -41.248155623 | -0.080344180 | -41.400679120 | -0.552009055 | 0.495712791 |
| 0 | 0.110 | 307.339269900 | -32.098788411 | 307.485439900 | -30.748543994 | 1.358133155 |
| 1 | 0.075 | -0.306186751 | -1.201828162 | -0.638113965 | -0.779611058 | 0.537068858 |
| 1 | 0.110 | -23.900888966 | -0.807071140 | -24.040231620 | -0.809149626 | 0.139358156 |
| 2 | 0.110 | 265.128475710 | -27.540014914 | 265.324856510 | -27.251269103 | 0.349198457 |
| 2 | 0.170 | 26.610353589 | -0.802878082 | 26.712699620 | -0.783344318 | 0.104193464 |

The analytic two-level rows are:

| omega (Ry) | Re GF | Im GF | Re analytic | Im analytic | relative difference |
|---:|---:|---:|---:|---:|---:|
| 0.110 | -45.050811081 | -1.438617278 | -44.395116537 | -1.479837218 | 1.4790e-2 |
| 0.310 | 36.225616127 | -0.990395174 | 36.336609134 | -0.990998431 | 3.0535e-3 |
| 0.370 | 23.430754998 | -0.414694333 | 23.522086547 | -0.415095645 | 3.8822e-3 |

The finite residual in the analytic rows is the controlled finite-
`green_eta`/finite real-axis quadrature envelope; it is below the fixture
acceptance limit of 2%.  The static fixture uses the exact provider path and
therefore compares at machine precision; its expected circular static value
is `-20 Ry^-1`, and `chi0 B_xc = m` with `B_xc=-0.05 Ry` gives matching
eigenpair and K-GF Ward residuals.

## Scope boundary

Generalized-overlap/HOH, GBT, CCOR, Hubbard, SOC, and noncollinear response
remain rejected by the production TDDFT guards.  Native real-space
`G(R,z)->chi0(R,omega)` and mixed/complex-contour integration remain later
milestones.
