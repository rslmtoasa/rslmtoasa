# TDDFT Green-function energy integration

## Scope

TDDFT-08 adds a mixed Lounis/Mills contour path to the K-space Green-function
bubble while retaining the direct real-axis trapezoid as the numerical
reference.  The implementation is selected with:

```fortran
gf_integration = 'direct'          ! reference
gf_integration = 'mixed_contour'   ! RR/AA contour plus near-Fermi RA
```

The contour controls are `contour_points`,
`contour_subdivisions`, `near_fermi_points`, and the optional
`contour_height`.  Energies are in Rydberg.  The default height is selected
from the energy window, Green-function broadening, and temperature; an
explicit finite-temperature height is rejected when it approaches the first
Fermi-function pole.

## Decomposition and convention

The direct path evaluates the existing retarded bubble in the project
convention,

```text
chi0(w) = -1/(2 pi i) integral_W dE f(E) {
  Tr[A Gq^R(E+w) B (Gk^R(E)-Gk^A(E))]
  + Tr[A (Gq^R(E)-Gq^A(E)) B Gk^A(E-w)]
}.
```

Define

```text
RR(E) = Tr[A Gq^R(E+w) B Gk^R(E)]
RA(E) = Tr[A Gq^R(E+w) B Gk^A(E)]
AA(E) = Tr[A Gq^A(E) B Gk^A(E-w)].
```

Changing variables in the positive `Gq^R(E) Gk^A(E-w)` term gives

```text
chi0(w) = -1/(2 pi i) [ I_RR - I_RA - I_AA ],
```

where `RR` is integrated on an upper-half-plane contour, `AA` on the
corresponding lower contour, and the nonanalytic mixed term is

```text
I_RA = integral_W dE C_w(E) RA(E).
```

For an infinite window `C_w(E) = f(E)-f(E+w)`.  The implementation keeps the
exact finite-window correction.  For positive `w`, the coefficient is
`-f(E+w)` on `[Emin-w,Emin]`, `f(E)-f(E+w)` on the central interval, and
`f(E)` on `[Emax-w,Emax]`; the intervals are reversed with the corresponding
signs for negative `w`.  At finite temperature the central interval is
limited to 50 `kT` around the two Fermi edges, where the omitted occupation
tail is below the existing occupation cutoff.  The boundary intervals are
still retained, so the result does not depend on silently enlarging the
integration window.

The source contract is explicit: mixed integration calls
`green_function_provider%get_complex` at arbitrary complex energies.  The
eigenpair provider supplies this through the Lehmann resolvent.  A source
that only supplies sampled real-axis matrices receives an explicit failure;
real-axis interpolation is not treated as analytic continuation.

## Contour and quadrature

The RR and AA paths are three-segment polylines: vertical from the real-axis
endpoint, horizontal at `+H` or `-H`, and vertical back to the other endpoint.
Each segment uses Gauss-Legendre quadrature.  Horizontal segments are split
into `contour_subdivisions` pieces because a single global quadrature can
under-resolve narrow one-electron resolvent features even when the contour is
well away from the real axis.  The near-Fermi RA and finite-window boundary
pieces use `near_fermi_points` Gauss-Legendre nodes.

At zero temperature the same-contour interval ends at the Fermi level and
the occupation is exactly one there.  At finite temperature the full supplied
window is used, with the complex Fermi occupation evaluated on the contour.
The default finite-temperature height is below `pi*kT`, and the code rejects
contours that approach the first Fermi pole.  `eta` and `green_eta` remain
numerical broadenings; they are never reported or interpreted as a physical
magnon linewidth.

## Validation and measured performance

The CMake test `UnitTddftContour` compares both paths on a two-level circular
response fixture.  The direct reference uses 16,001 energy nodes at zero
temperature and 20,001 nodes at 2,000 K.  The following values are emitted by
the test on the development build; the error is the maximum absolute complex
difference divided by `max(1,max(abs(chi_direct)))`.

| case | contour nodes | horizontal subdivisions | near-Fermi nodes | normalized error | contour height | GF evaluations |
|---|---:|---:|---:|---:|---:|---:|
| zero T, fine | 64 | 8 | 256 | 1.67e-5 | 1.00 | 7,168 |
| zero T, coarse | 24 | 8 | 256 | 3.90e-5 | 1.00 | 3,968 |
| 2,000 K | 64 | 8 | 384 | 3.10e-6 | 9.95e-3 | 3,938 |

The GF count is obtained from result metadata; it is not a fixed algorithmic
constant because zero-weight quadrature intervals can be skipped.  The test also checks that the mixed
path's maximum complex energy is above `green_eta`, that the positive-frequency
circular response has the retarded sign, and that the mixed path uses fewer GF
evaluations than the direct reference.  On the recorded run the zero-
temperature direct/mixed CPU samples were approximately 5.50e-2 s and
1.17e-2 s, respectively, with 96,006 versus 7,168 GF evaluations.  These are
fixture measurements, not a production-material benchmark.

The direct reference remains available precisely because contour-node,
subdivision, height, and near-Fermi convergence must be checked for each
material and frequency window.  The backend reports all of these settings,
the maximum imaginary energy, evaluation count, integration window, and CPU
time in `tddft_chi0_metadata` and in text output headers.

## Native real-space status

TDDFT-R2-03 extends the native provider with a genuine arbitrary-complex-energy
source.  `tddft_native_green_source` is a thin adapter over the mature
`green%calculate_intersite_gf_complex` API; recursion, terminators, block
phase reconstruction, and Chebyshev reconstruction remain in `green_mod` and
its submodules.  The production chain is still

```text
G_ab(R,z), G_ba(-R,z) -> chi0_ab(R,omega) -> chi0_ab(q,omega).
```

With `gf_integration='mixed_contour'`, the native provider uses the same
RR/AA polyline plus finite-window RA decomposition as the K-space provider.
Each quadrature segment sends all required complex energies to the source in
one batch, traces every requested response component from that batch, and
Fourier-transforms the resulting `chi0(R,omega)` once for the complete q
batch.  Identical frequencies in a request reuse the already integrated
real-space response.  No `G(R,z)` to `G(k,z)` conversion or global analytic
continuation is performed.  The direct near-real-axis trapezoid remains
available as the reference route, and mixed integration without a complex
source fails explicitly.

The native mixed path reports the contour node counts, subdivisions,
near-Fermi mesh, contour height, maximum imaginary energy, source GF
evaluation count, and separate pair-integration/Fourier CPU samples in
`tddft_chi0_metadata`.  `UnitTddftRealspaceGF` checks the zero-frequency and
low-energy responses,
the same two-level poles inside the Stoner continuum, direct/contour
agreement, contour-node stability, complex-offset metadata, repeated-frequency
reuse, and reuse over a two-q batch.  Its development fixture produced a
direct/low/high contour pair-integration timing sample of approximately
`1.59e-3 / 1.41e-4 / 2.60e-3 s`; these are diagnostic fixture timings, not a
production benchmark.

## Analytic continuation and invariants

No Halle finite-complex-frequency or analytic-continuation mode is exposed by
this milestone.  It remains an explicitly deferred experimental follow-up;
therefore no continuation stability claim is made.  The mixed implementation
does not alter the response vertices, retarded denominator convention,
spectral sign, or normalization.  `UnitResponseConventions`, the existing
Green-function response regressions, and `UnitTddftContour` cover those
invariants at their respective abstraction levels.  Numerical `eta` remains
metadata-marked as nonphysical.
