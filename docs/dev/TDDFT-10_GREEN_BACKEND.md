# TDDFT-10 — eigenpair-resolvent reference backend

## Current production status

`chi0_backend = 'green'` is an eigenpair-resolvent reference implementation:
it constructs the one-electron resolvent from the same arbitrary-k eigenpairs
as the Lehmann backend. It is **not** connected to a native RS recursion,
surface, impurity, or reciprocal Green-function provider. Production
transverse runs require the real static Ward solver and therefore currently
require `chi0_backend = 'eigenpairs'`; the driver rejects `green` for that
route. Full response is separately capability-gated and is unavailable with
the shipped XC derivative providers.

The numerical reference and its tests remain useful, but this document is not
a supported production-input recipe until a native provider and a compatible
static-limit solver are validated.

- [x] Eigenpair-resolvent reference and real-axis integration tests.
- [x] Full-spinor site response vertices and Green/eigenpair comparison tests.
- [ ] Native RS Green-function provider integration.
- [ ] Static-limit and production-route validation for a native provider.
- [ ] Complex-contour acceleration (deferred).

## Reference retarded bubble

For a left response vertex `A` and right perturbing vertex `B`, the code
evaluates the equilibrium real-axis identity

```text
chi_AB^R(q,w) = -1/(2*pi*i) int dE f(E) Tr[
  A G_q^R(E+w) B (G_k^R(E)-G_k^A(E))
  + A (G_q^R(E)-G_q^A(E)) B G_k^A(E-w)
].
```

Here `G^A(E) = G^R(E)†`. Inserting the spectral representation of the
one-electron Green function yields the same convention as the production
Lehmann backend:

```text
(f_nk-f_m,k+q) <n,k|A|m,k+q><m,k+q|B|n,k>
------------------------------------------------ .
              w + e_nk - e_m,k+q + i*eta
```

The validated reference is a trapezoidal integral on a configurable real
energy mesh. A one-particle Lorentzian half-width `green_eta` produces a
bubble linewidth of approximately `2*green_eta`; the default
`green_eta=0` therefore selects `eta/2` so that `eta` keeps its established
response-broadening meaning. Finite mesh spacing, finite energy limits, and
finite one-particle broadening leave small differences from the analytic
eigenpair sum. The unit comparison uses 16,001 points over `[-0.35,0.35] Ry`
and requires agreement within 2 percent for three commensurate q values and
four frequencies.

## Native-provider extension boundary

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

## Reserved input

The controls remain parsed for the reference implementation, but do not use
this block for a production transverse calculation today. The production
metadata deliberately labels this choice `eigenpair-resolvent reference`.

```fortran
&tddft
  chi0_backend = 'green'
  ! zero means eta/2; positive values set the one-particle half-width
  green_eta = 0.0
  ! Omit green_energy_min/max to infer a safe window from the source spectrum.
  green_energy_points = 2001
/
```
