# RS-LMTO libXC production contract

This document defines the production boundary between RS-LMTO-ASA and libXC.
It is an interface contract, not a replacement for the historical XC kernels.
The legacy `TXC=1–99` implementations remain available and are not silently
redirected to libXC.

## Supported family boundary

The production libXC interface supports:

| libXC family | Status | Required ASA data |
| --- | --- | --- |
| LDA | supported | spin density |
| GGA | supported | spin density and radial spin-density gradient |
| MGGA | not supported | kinetic-energy density and/or Laplacian data not supplied by ASA |
| HYB_LDA, HYB_GGA, HYB_MGGA | not supported | exact exchange and, for HYB_MGGA, additional orbital ingredients |
| LCA, OEP, orbital-dependent, kinetic-density-dependent families | not supported | orbital or nonlocal ingredients not supplied by ASA |

The family is queried from every active native libXC object at construction.
Only LDA and GGA objects may become active.  A combination with at least one
GGA component is routed as GGA; an all-LDA combination is routed as LDA.  The
family of the first component is never used as a proxy for the combination.

## Selector namespaces

`TXC` is an RS-LMTO selector, while `libxc_func_id(:)` contains native libXC
IDs only:

| Selector | Meaning |
| --- | --- |
| `1–99` | historical RS-LMTO implementation |
| `100–199` | predefined explicit libXC alias |
| `>=1000` | direct native libXC ID `TXC-1000` |

Examples:

```text
TXC=1001  -> native ID 1 only
TXC=1012  -> native ID 12 only
TXC=1101  -> native ID 101 only
```

The active native-ID list for a direct request always has exactly one entry.
No exchange or correlation partner is inferred from a functional name.  A
direct exchange-only or correlation-only request emits a warning, but is
evaluated exactly as requested.

Legacy reference mappings are documentation and comparison metadata only.
They are never inserted into `libxc_func_id(:)` for `TXC=1–99`.

## Predefined aliases

The aliases below are explicit native-ID combinations.  At runtime the
interface queries and reports each native libXC name, family, ID, and the
mapping quality retained in the selector metadata.

| TXC | Native IDs | Runtime functional names in libXC 5.2.3 | Mapping quality |
| ---: | --- | --- | --- |
| 101 | `[1,17]` | Slater exchange + von Barth & Hedin | `REFERENCE_EQUIVALENT` |
| 102 | `[1,24]` | Slater exchange + Gombas | `APPROXIMATE_ANALOGUE` |
| 103 | `[1]` | Slater exchange | `APPROXIMATE_ANALOGUE` |
| 104 | `[1,9]` | Slater exchange + Perdew & Zunger | `REFERENCE_EQUIVALENT` |
| 105 | `[1,12]` | Slater exchange + Perdew & Wang | `REFERENCE_EQUIVALENT` |
| 106 | `[1,7]` | Slater exchange + Vosko, Wilk & Nusair (VWN5) | `REFERENCE_EQUIVALENT` |
| 107 | `[1,5]` | Slater exchange + Gunnarson & Lundqvist | `APPROXIMATE_ANALOGUE` |
| 108 | `[101,130]` | Perdew, Burke & Ernzerhof + Perdew, Burke & Ernzerhof | `REFERENCE_EQUIVALENT` |
| 109 | `[117,130]` | Hammer, Hansen, and Norskov + Perdew, Burke & Ernzerhof | `APPROXIMATE_ANALOGUE` |

`REFERENCE_EQUIVALENT` means that the named parametrization is a controlled
comparison reference; it does not assert bitwise identity with a legacy
kernel. `APPROXIMATE_ANALOGUE` explicitly marks a related but non-identical
choice.  In particular, aliases do not imply that historical BH, BHJ, ASW,
or legacy GGA formulas have been replaced.

## Spin and density conventions

The physical density supplied by RS-LMTO is in electrons/bohr³ and the radial
coordinate is in bohr.  libXC receives the standard two-spin ordering:

```text
rho(:,1) = rho_up
rho(:,2) = rho_down
```

The pointwise historical wrapper has an explicit conversion boundary:

| Quantity | Historical wrapper slot |
| --- | --- |
| input `RHO1` | spin down |
| input `RHO2` | spin up |
| output `V1` | spin-down potential |
| output `V2` | spin-up potential |

Thus the wrapper passes `[RHO2,RHO1]` to libXC and maps returned `vrho(1)` to
`V2` and `vrho(2)` to `V1`.  `VXC0SP` performs the final swap required by its
potential-array convention: array channel 1 is up and channel 2 is down.
For the historical magnetic-field diagnostic used by this interface,

\[
B_{\rm xc}^{\rm historical} = (V_{\rm down}-V_{\rm up})/2.
\]

The response-provider field named `bxc_energy` is separately documented as
the coefficient of `sigma_z` in its `H=v\,I+B\,sigma_z` convention, namely
`(V_up-V_down)/2`.  These two named conventions must not be substituted for
one another.

## Energy and potential units

libXC returns the energy per particle, `exc`, and density derivatives in
Hartree-based atomic units.  RS-LMTO stores XC energies and multiplicative
potentials in Rydberg.  The sole conversion boundary is:

```text
1 Ha = 2 Ry
```

For LDA, `exc` and `vrho` are each multiplied by two after summing the active
native components.  For GGA, the radial functional derivative is assembled in
Hartree first and the complete result, including its divergence, is multiplied
by two once.  No later caller applies another XC unit conversion.

## Radial GGA formulation

The logarithmic ASA mesh is

\[
r_i=b[\exp(a(i-1))-1], \qquad dr/dx=r+b, \qquad x=a(i-1).
\]

The radial helper receives `rho_up`, `rho_down`, `d rho_up/dr`, and
`d rho_down/dr`.  For libXC it forms

```text
sigma = [(d rho_up/dr)^2,
         (d rho_up/dr)(d rho_down/dr),
         (d rho_down/dr)^2]
```

With `vsigma=[v_uu,v_ud,v_dd]`, it constructs the radial fluxes

```text
F_up   = 2*v_uu*d rho_up/dr + v_ud*d rho_down/dr
F_down = 2*v_dd*d rho_down/dr + v_ud*d rho_up/dr
```

and returns

\[
v_s(r)=v_{\rho_s}(r)-\frac{1}{r^2}\frac{d}{dr}[r^2F_s(r)].
\]

LDA contributions in a mixed LDA+GGA list are accumulated into `exc` and
`vrho` in the same radial call; GGA contributions additionally accumulate
`vsigma`.  This is why mixed combinations use the radial GGA route.

At the extrapolated origin, a regular radial flux uses
`div(F rhat)=3 F'(0)`.  The outer endpoint uses the five-point backward
stencil in logarithmic coordinate from `radgra`, then converts with `dr/dx`;
the endpoint flux is not silently forced to zero.  Interior and endpoint
derivative convergence are covered by `UnitRadialGga`.

## Density floor

`LIBXC_DENSITY_FLOOR = 1e-20 electrons/bohr³` is input protection for libXC
only.  It does not modify RS-LMTO density arrays, radial gradients, or
quadrature weights.  An exactly zero total density returns zero energy and
zero potential.  A positive total density with one zero spin channel is
evaluated as a finite polarized limit, with only the zero libXC input channel
regularized.

## Errors and fallback policy

* A libXC selector in a build without libXC fails during XC construction; it
  never falls back to a legacy selector.
* An undefined 100-series alias fails with a diagnostic directing the user to
  a direct native-ID selector.
* An invalid native ID, MGGA, hybrid, orbital-dependent, kinetic-density-
  dependent, or otherwise unsupported family fails before atomic SCF.
* The pointwise wrapper rejects GGA evaluation and directs callers to the
  full radial helper.
* Legacy selectors continue to use the historical implementation and do not
  require libXC.

Every temporary `xc_f03_func_t` used for metadata or pointwise evaluation is
ended before its metadata or object is reused.  The production `xc` object
owns only copied IDs and metadata; it does not retain libXC objects across
SCF iterations.

## Validation evidence

The compact regression suite is registered in CMake as
`UnitLibxcProductionContract` when libXC is enabled.  It independently
initializes native libXC objects and checks:

* every predefined alias, native ID, runtime name, component family, and
  aggregate route;
* direct IDs `1`, `12`, `101`, and `130` without implicit pairing;
* native LDA IDs `[1,5,7,9,12,17]` over unpolarized, polarized, near-full-
  polarization, and low-density points;
* the analytic spin-polarized LDA exchange oracle for `exc`, `v_up`, and
  `v_down`;
* standard libXC spin ordering, historical `V1/V2` conversion, the signed
  historical `B_xc`, and one Hartree-to-Rydberg conversion;
* PBE, RPBE, exchange-only GGA, correlation-only GGA, and an explicit mixed
  LDA+GGA evaluator list;
* zero density, one-spin-channel-zero, and density-floor behavior; and
  repeated construction and cleanup.

Additional registered tests provide the full radial evidence:

* `UnitLibxcGgaRadial` checks sigma ordering and the central finite-difference
  variational identity for smooth spherical PBE densities;
* `UnitPbeGgaComparison` compares the production PBE radial path against the
  legacy reference path on one mesh;
* `UnitRadialGga` checks analytic derivatives, the regular origin limit, and
  outer-boundary refinement/no-artifact behavior;
* `UnitXcLdaReconciliation` checks polarized LDA grids and the exchange oracle;
* `UnitXcSelectorSemantics` enforces the disjoint selector namespaces; and
  `UnitXcSelectorRequiresLibxc` is a negative test in a build without libXC.

Material calculations are smoke tests only.  No Fe moment, material-specific
threshold, or legacy XC formula is used as the primary certification oracle.
