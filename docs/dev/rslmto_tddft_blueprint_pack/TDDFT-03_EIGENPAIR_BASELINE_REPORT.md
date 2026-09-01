# TDDFT-03 eigenpair baseline report

This report records the transparent TDDFT-03 reference path.  The numerical
values below are deterministic unit-fixture evidence, not fitted material
data.

## Implemented contract

- Dynamic `chi_KS` retains separate `(n,k)` and `(m,k+q)` eigenvalue and
  eigenvector endpoints through the tiled transition engine.
- `build_static_chi_ks_from_eigenpairs_at_q` evaluates the actual
  finite-temperature static divided difference for separate endpoints.  The
  existing `build_static_chi_ks_from_eigenpairs` API remains the explicit
  `q=0` wrapper.
- Equal and near-degenerate energies use the analytic finite-temperature Fermi
  derivative.  Static response does not consume dynamic `eta`; dynamic output
  labels `eta` as numerical broadening rather than a physical linewidth.
- Result and text-output metadata record the direct-basis q, endpoint
  provenance, k mesh, omega grid, Fermi level, temperature, projection, eta
  role, and static/dynamic status.
- The production entry point rejects non-`nsp=1`, nonzero-SOC, generalized
  overlap, HOH/second-order, GBT/explicit-texture, CCOR, orbital-polarization,
  Hubbard, and constrained-field branches until their response derivatives are
  derived.

The reciprocal q path remains in direct reciprocal coordinates.  The shifted
k-point workset folds Hamiltonian evaluation points into `[-0.5,0.5)` while
the requested, possibly unwrapped path coordinate is retained in provenance.
Gamma, boundary folding, negative coordinates, `q+G`, weight preservation, and
base-mesh immutability are covered by `UnitKpointWorkset`.

## Representative response values

`UnitTddftChiKS` reports the following fixed one-site, two-band fixture values:

| case | settings | representative result |
|---|---|---|
| dynamic `chi_KS(omega)` | `nk=4`, `T=100 K`, `eta=0.002 Ry`, `omega=0.100 Ry` | `Re chi=0`, `Im chi=-2.0000000000000000e3`, `-Im chi/pi=6.3661977236758139e2` |
| static `chi_KS(q,0)` | shifted endpoint, direct `q=(0,0,0.25)`, `nk=4`, `T=700 K`, `eta=0` | `Re chi=-1.6134844283791885e1`, `Im chi=-0`, static Stoner weight `-0` |

The dynamic result is independently checked against an explicit pair sum.  The
static shifted result is independently checked against the finite-temperature
divided-difference sum, including the no-eta path.  The same test also opens
the text output and verifies q, endpoint, eta-role, and omega-grid provenance.

## Verification and profiling

The TDDFT-focused CTest selection passed all 18 tests, including the response
conventions, q-workset/arbitrary-k mapping, dynamic and static chi0, transition
workspace, direct Xi, Green reference, Ward/Goldstone, Dyson/modes, config,
dispatch, and cross-milestone equivalence tests.

`UnitTddftCpuProfile` passed for both deterministic fixtures.  Its measured
TDDFT phase timings were:

| fixture | dimensions | dominant measured phases (seconds) |
|---|---|---|
| bccFe-labelled one-site fixture | `nk=16`, `nw=96` | GF energy integration `6.6044e-2`; response accumulation `1.8052e-2`; denominator generation `9.9850e-3` |
| fccNi-labelled two-site fixture | `nk=32`, `nw=192` | GF energy integration `4.9639e-1`; response accumulation `2.8603e-1`; denominator generation `1.5521e-1` |

These are profiling observations only; no material optimization or empirical
spectral shift was introduced.

## Remaining TDDFT-03 evidence

The repository contains a bcc-Fe production deck, but it currently requests
`nsp=2`; the new initial-boundary guard correctly rejects that relativistic
mode.  No converged real-material Fe/Ni spectra are therefore claimed by this
baseline report.  Small Fe/Ni production runs, k/R-cutoff convergence, and
material Ward/acoustic validation remain the explicit follow-up campaign once
SOC-free `nsp=1` inputs and the required restart/database state are supplied.
