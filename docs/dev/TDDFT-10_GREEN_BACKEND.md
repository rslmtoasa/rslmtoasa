# TDDFT-10 — Green-function Kohn-Sham response backend

- [x] Implement `green_chi0_provider` matching the canonical KS-response provider interface.
- [x] Reuse the periodic one-electron resolvent through an eigenpair Green-function adapter.
- [x] Implement site projectors and charge/spin/circular vertices in full spinor form.
- [x] Support the same site response basis as the eigenpair backend.
- [x] Add real-axis energy-window, mesh, and broadening controls.
- [x] Keep complex-contour acceleration deferred until after real-axis validation.
- [x] Compare periodic finite-mesh Green and eigenpair responses over several q and omega values.
- [x] Demonstrate an inequivalent two-site open local-response geometry.

## Retarded bubble implemented

For a left response vertex `A` and right perturbing vertex `B`, the code
evaluates the equilibrium real-axis identity

```text
chi_AB^R(q,w) = -1/(2*pi*i) int dE f(E) Tr[
  A G_q^R(E+w) B (G_k^R(E)-G_k^A(E))
  + A (G_q^R(E)-G_q^A(E)) B G_k^A(E-w)
].
```

Here `G^A(E) = G^R(E)†`.  Inserting the spectral representation of the
one-electron Green function yields exactly the production Lehmann convention:

```text
(f_nk-f_m,k+q) <n,k|A|m,k+q><m,k+q|B|n,k>
------------------------------------------------ .
              w + e_nk - e_m,k+q + i*eta
```

The validated reference is a trapezoidal integral on a configurable real
energy mesh.  A one-particle Lorentzian half-width `green_eta` produces a
bubble linewidth of approximately `2*green_eta`; the default
`green_eta=0` therefore selects `eta/2` so that `eta` keeps its established
response-broadening meaning.  Finite mesh spacing, finite energy limits, and
finite one-particle broadening leave small differences from the analytic
eigenpair sum.  The unit comparison uses 16,001 points over `[-0.35,0.35] Ry`
and requires agreement within 2 percent for three commensurate q values and
four frequencies.  Production outputs record the actual energy window,
point count, one-particle width, and effective bubble width.

## Extension boundary and local geometry

`green_chi0_provider` consumes only the public `green_function_provider`
methods `get_retarded(branch,k,E,eta,G)` and `get_spectral_bounds`.  Its
periodic adapter builds resolvents from the same arbitrary-k eigenpairs used
by the reference backend.  A surface, impurity, recursion, or embedded Dyson
route can implement those two methods without exposing its private array
layout to TDDFT.  The two-site open test has inequivalent local levels and
checks its site-resolved response through this exact boundary.

For a frequency-dependent `Sigma(z)`, providing a dressed one-particle Green
function computes a dressed bubble only.  That is not automatically a
conserving many-body response: the vertex corrections required by the chosen
self-energy approximation are separate physics and are intentionally not
implied by this TDDFT-10 backend.

## Input

```fortran
&tddft
  chi0_backend = 'green'
  ! zero means eta/2; positive values set the one-particle half-width
  green_eta = 0.0
  ! Omit green_energy_min/max to infer a safe window from the source spectrum.
  green_energy_points = 2001
/
```
