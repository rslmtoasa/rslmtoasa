# XCR-03 — radial GGA and libXC reconciliation

## Scope

This change repairs the spherical GGA derivative path for the legacy PBE
selector (`TXC=8`) and replaces the pointwise libXC GGA potential path used by
`TXC=108` with a complete radial functional derivative.  LDA selectors retain
their existing pointwise route.

## Mesh and derivative contract

The atomic radial mesh is

\[
  r_i=b\,[\exp(a(i-1))-1], \qquad x_i=a(i-1), \qquad
  \frac{dr}{dx}=r+b.
\]

`radgra` applies the existing finite-difference stencil in `x`, then divides
by `r+b`.  Its output is therefore the physical derivative `df/dr`, not a
logarithmic-mesh derivative.  A second call with the first derivative as its
input returns `d2f/dr2`.  The implementation is now in
`source/xc_radial.f90`, and `self.f90` uses this single implementation.

The radial density reconstructed by `VXC0SP` is

\[
  n_s(r)=\frac{R_s(r)}{4\pi r^2},
\]

with the origin value supplied by the existing regular extrapolation.  After
reconstruction, `RHOP` and `RHOPP` mean `dn_s/dr` and `d2n_s/dr2` explicitly.
No later `1/r` or `1/r2` transformation is applied.  This removes the former
double transformation in the GGA calls.

## Spin and libXC contracts

The historical `XCPOT_hybrid` interface has this ordering:

| quantity | slot 1 | slot 2 |
| --- | --- | --- |
| input `RHO1/RHO2` | down | up |
| output `V1/V2` | up | down |

The new mesh-level libXC helper uses explicit up/down arguments and performs
the conversion at the boundary.  libXC receives

```text
rho  = [rho_up, rho_down]
sigma = [grad_up**2, grad_up*grad_down, grad_down**2]
```

The `vrho` outputs are accumulated in up/down order.  For the GGA flux, with
`vsigma = [v_uu, v_ud, v_dd]`, the helper forms

```text
F_up   = 2*v_uu*grad_up + v_ud*grad_down
F_down = 2*v_dd*grad_down + v_ud*grad_up
```

and returns the Hartree-to-Ry-scaled potential.  The factor of two is applied
to the complete functional derivative, including the radial divergence, not
only to the pointwise libXC outputs.

For legacy PBE spin scaling, exchange evaluates the transformed density
`2*n_s`; the first and second derivatives are transformed too:
`2*n'_s` and `2*n''_s`.  Correlation uses the total-density spin variables as
before.

## Complete radial functional derivative

For a spherical spin channel, the multiplicative GGA potential is

\[
  v_s(r)=v_{\rho_s}(r)-\nabla\cdot F_s(r)\hat r,
\]

where the helper evaluates

\[
  \nabla\cdot[F(r)\hat r]=\frac{1}{r^2}\frac{d}{dr}[r^2F(r)].
\]

`radial_flux_divergence` differentiates `r2*F` over the complete mesh.  At the
origin it uses the regular limit `3 F'(0)` for `F(r)=F'(0)r+O(r3)`.  The
legacy PBE pointwise origin uses

\[
  \lim_{r\to0}\left(n''+\frac{2n'}r\right)=3n''(0),
\]

with `n'(0)=0`.  Neither path inserts an epsilon into `1/r`.

The old pointwise libXC GGA wrapper is no longer allowed to return an
incomplete potential: it reports a fatal diagnostic and is bypassed by
`VXC0SP`, which calls the mesh-level helper once for the full radial table.

## Validation

The following tests were added or run with GNU Fortran 13.3 and libXC 5.2.3:

| check | result |
| --- | --- |
| `UnitRadialGga` analytic `exp(-alpha*r)` and `r2*exp(-alpha*r)` derivatives | PASS; fine-mesh first-derivative error `7.02e-7`, second-derivative error `2.26e-5` |
| `UnitRadialGga` regular-origin derivative and flux-divergence limits | PASS |
| `UnitLegacyPbeGga` finite-difference functional derivative | PASS; final epsilon error `1.85e-8` |
| `UnitLibxcGgaRadial` sigma/channel mapping and finite-difference derivative | PASS; final epsilon error `8.93e-10` |
| `UnitPbeGgaComparison` on the same asymmetric smooth radial densities | PASS; energy-density difference below `9e-16`, maximum scalar/magnetic potential differences `5.06e-6`/`5.30e-6` |
| Existing LDA baseline, reconciliation, selector, and legacy-kernel tests | PASS |

The finite-difference checks use the radial volume integral
`4*pi*r2*dr*(rho_up+rho_down)*exc` and compare it with
`4*pi*r2*dr*(v_up*d rho_up + v_down*d rho_down)` over a decreasing epsilon
sweep.  The perturbations vanish at the mesh boundary and are regular at the
origin.

The production smoke cases below used four SCF steps, so their residuals are
not convergence claims.  They are crash/NaN and route-consistency checks.

| structure | selectors exercised | result |
| --- | --- | --- |
| bcc Fe | TXC=1, 8, 108 | all three exited successfully |
| fcc Fe | TXC=1, 8, 108 | all three exited successfully |

All six smoke runs exited successfully without `ERROR`, `FATAL`, or `NaN`; the
reported energies, moments, and residuals were finite.
