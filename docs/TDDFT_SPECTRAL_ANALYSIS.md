# RS-LMTO TDDFT spectral analysis

This document defines the common post-`chi0` spectral products used by the
eigenpair, K-space Green-function, and native real-space Green-function
routes.  The implementation is in
[`source/tddft_dyson.f90`](../source/tddft_dyson.f90) and
[`source/tddft_modes.f90`](../source/tddft_modes.f90).  It applies to the
initial collinear, SOC-free, orthogonal (`ham_only`) response boundary.

## Loss matrix and Stoner reference

For a retarded susceptibility `chi(omega)`, the response-space loss matrix is

```text
L(omega) = -(chi(omega) - chi(omega)^H) / (2 i pi).
```

It is Hermitian by construction.  For a scalar diagonal response this reduces
to `L_ii = -Im chi_ii / pi`; with the project convention, positive-frequency
absorption therefore has positive loss.  Off-diagonal entries are not formed
with an elementwise `-Im`: they are the Hermitian anti-Hermitian-part product
and retain channel/site interference.

The Dyson result retains both response levels:

```text
chi_KS, L_KS       bare Kohn-Sham response and its Stoner spectrum
Xi = chi_KS K      enhancement product
chi, L             enhanced response and loss matrix
```

`chi_KS_site_loss` and `chi_KS_trace_loss` in a Dyson file are the Stoner
continuum reference.  `site_loss` and `trace_loss` are the corresponding
enhanced products.  The chi0 writer additionally emits `stoner_landau`
records: site-diagonal KS loss at a fixed q.  The per-q files make this a
q-resolved Landau-map-style diagnostic without requiring a separate runtime
representation.

Every Dyson frequency stores `max|L-L^H|` in
`loss_hermiticity_residual`.  When loss diagonalization is requested, the
Hermitian eigensystem is retained as `loss_eigenvalues` and
`loss_eigenvectors`; text output uses `loss_mode` and `loss_mode_vector`
records.  The eigenvalues are the spectral weight of orthogonal response
modes, not poles of `Xi`.

## Collective-mode branches

`Xi` is generally non-Hermitian.  At every `(q, omega)` the analyzer retains
right eigenvectors and left eigenvectors and normalizes safe pairs so
`l^H r = 1`.  A collective candidate is an interpolated sign-changing
crossing of

```text
Re(lambda_Xi(q, omega)) = 1,
```

with both bracketing eigensystems sufficiently conditioned.  A grid point that
is merely close to unity is retained as a diagnostic candidate but is never
reported as a collective crossing.

Branch labels are assigned in two stages.  Along frequency, and then along
the supplied q path, the analyzer solves a global maximum-weight assignment
for each adjacent slice.  The score combines:

- biorthogonal left/right eigenvector overlap (`0.55` weight);
- response-space character overlap from component intensities (`0.20`);
- eigenvalue continuity (`0.25`, with a bounded inverse-distance form).

The global assignment matters at crossings: a greedy first match can consume
the partner required by the next mode and silently exchange acoustic/optical
or sublattice branches.  The output retains the selected branch overlap,
character overlap, eigenvalue step, conditioning, and exceptional-point flag.

The same frequency/q continuation is applied to the Hermitian loss
eigenvectors when loss matrices are supplied.  Thus the loss-mode spectrum and
the non-Hermitian enhancement branches are both available for interpretation;
neither is reduced to a scalar trace.

## Quasiparticle-quality gate and linewidths

The analyzer only attempts a Lorentzian fit after these objective gates pass:

1. an interpolated `Re(lambda_Xi)=1` crossing exists;
2. the crossing is bracketed by well-conditioned eigensystems;
3. `abs(Im(lambda_Xi))` is below the configured damping threshold;
4. when a dressed loss matrix is supplied, the selected mode has positive
   projected enhanced loss;
5. the enhanced loss peak is isolated, half-height bracketed, and resolved by
   the frequency grid.

The fit checks a controlled relative residual against a constant local
background.  Its `center` is the sampled peak energy, and `fwhm`/`hwhm` are
reported explicitly with `fwhm = 2*hwhm`.  A single-run width is an observed
width: the numerical `eta` is not subtracted.  A separate multi-eta linear
extrapolation is available for controlled estimates of the zero-eta
intercept.

If any gate fails, no peak energy or linewidth is claimed.  The feature is
labelled `noncollective Stoner feature` when no unity crossing exists, or
`overdamped/continuum-like ...` when an enhanced crossing fails its quality
criteria.  The total KS and enhanced spectra remain in the output so a failed
fit is still useful evidence about the Stoner continuum and Landau damping.

## Backend independence and evidence

The analyzer consumes only `Xi`, the enhanced loss matrix, and the scalar
trace supplied by the common Dyson layer.  It has no backend-specific branches.
Consequently identical `chi0` inputs produce identical spectral products
regardless of whether those inputs came from explicit transitions,
`G(k,z)`, or native `G(R,z)`.

The focused evidence is `UnitTddftDysonModes`:

- direct loss sign, Hermiticity, and off-diagonal pairing;
- loss eigenvalue/eigenvector retention and orthonormality;
- synthetic two-mode avoided and q-order crossings;
- biorthogonal projection for a non-normal `Xi`;
- Stoner/no-crossing classification;
- pre-fit rejection of a strongly damped crossing;
- rejection of a non-Lorentzian structure;
- accepted broadened Lorentzian with explicit FWHM/HWHM;
- multi-eta zero-broadening extrapolation.

These tests establish the post-processing contract.  They do not claim that a
material calculation has passed Fe/Ni convergence, Goldstone, stiffness, or
backend-equivalence release gates.
