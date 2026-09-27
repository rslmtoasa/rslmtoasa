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

For each component, initialization also retains its native ID, name, kind,
raw flags, 1D/2D/3D dimensionality, required ingredients, and EXC/VXC/FXC/KXC/LXC
availability. Metadata are copied before the temporary `xc_f03_func_t` is
ended; no native metadata pointer is used after that lifecycle boundary.

The compatibility policy is explicit. Kinetic functionals, non-3D
functionals, unsupported families, Laplacian/nonlocal requirements, invalid
IDs, and components without both EXC and VXC fail before SCF. Exchange-only,
correlation-only, multiple-X/C, combined-plus-additional, development-marked,
and unconventional but valid compositions warn and proceed without adding or
rewriting a component.

## Selector namespaces

`TXC` is an RS-LMTO selector, while `libxc_func_id(:)` contains native libXC
IDs only:

| Selector | Meaning |
| --- | --- |
| `1–99` | historical RS-LMTO implementation |
| `100–199` | predefined libXC XC bundle |
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

An arbitrary native component list is represented by `libxc_func_id(:)`; its
length and ordering are not semantically constrained. The supported internal
setter is `xc%set_libxc_component_ids([id1,id2,...])`, which revalidates and
rebuilds all metadata, marks the selection as explicit, and uses
`NO_EQUIVALENT` mapping quality. A user-facing namelist for this list remains
a future input-layer extension rather than another encoded TXC range.

## Predefined libXC XC bundles

The bundles below are explicit native-ID combinations.  At runtime the
interface queries and reports each native libXC name, family, ID, and the
mapping quality retained in the selector metadata.

| TXC | Native IDs | Runtime functional names in libXC 5.2.3 | Mapping quality |
| ---: | --- | --- | --- |
| 101 | `[1,17]` | Slater exchange + von Barth & Hedin | `REFERENCE_EQUIVALENT` |
| 102 | `[1,24]` | Slater exchange + Gombas | `APPROXIMATE_ANALOGUE` |
| 103 | `[1]` | Slater exchange | `APPROXIMATE_ANALOGUE` |
| 104 | `[1,9]` | Slater exchange + Perdew & Zunger | `REFERENCE_EQUIVALENT_UNPOLARIZED` |
| 105 | `[1,12]` | Slater exchange + Perdew & Wang | `REFERENCE_EQUIVALENT` |
| 106 | `[1,7]` | Slater exchange + Vosko, Wilk & Nusair (VWN5) | `REFERENCE_EQUIVALENT` |
| 107 | `[1,5]` | Slater exchange + Gunnarson & Lundqvist | `APPROXIMATE_ANALOGUE` |
| 108 | `[101,130]` | Perdew, Burke & Ernzerhof + Perdew, Burke & Ernzerhof | `REFERENCE_EQUIVALENT` |
| 109 | `[117,130]` | Hammer, Hansen, and Norskov + Perdew, Burke & Ernzerhof | `APPROXIMATE_ANALOGUE` |

`REFERENCE_EQUIVALENT` means that the named parametrization is a controlled
comparison reference; it does not assert bitwise identity with a legacy
kernel. `REFERENCE_EQUIVALENT_UNPOLARIZED` has the additional qualification
that the historical counterpart is unpolarized-only: TXC=104 is the libXC
reference counterpart to historical TXC=7 only in the unpolarized limit,
subject to the established normalization and unit conventions.
`APPROXIMATE_ANALOGUE` explicitly marks a related but non-identical choice.
In particular, bundles do not imply that historical BH, BHJ, ASW, or legacy
GGA formulas have been replaced.

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

`control%nsp` is not the number of XC channels. Its global modes are `1`
(collinear SR), `2` (collinear FR/SOC), `3` (noncollinear SR), and `4`
(noncollinear FR/SOC). The semantic queries in `control_mod` are
`is_collinear()`, `is_noncollinear()`, `has_soc()`,
`uses_spinor_representation()`, and `is_spin_polarized_mode()`.
The atomic radial/XC path uses two local spin-density channels even for global
`nsp=1`; noncollinear XC uses local eigenchannels
`n_±=(n±|m|)/2`. The current density, not the global mode integer, determines
whether a legacy unpolarized-only kernel is safe.

Legacy `TXC=6` (Wigner) and `TXC=7` (Perdew–Zunger) are explicitly
unpolarized-only. Equal local channels are accepted to the runtime tolerance;
an appreciable `rho_up-rho_down` is fatal with a spin-capable-functional
diagnostic. They are never silently evaluated with equal potentials for a
polarized density.

For arbitrary valid user-defined lists, initialization reports
`XC mode: custom libXC composition` and states that components are summed
exactly as requested. Unusual but compatible lists remain warnings, while
architecturally incompatible lists remain fatal.

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
* An undefined 100-series bundle fails with a diagnostic directing the user to
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

* every predefined bundle, native ID, runtime name, component family, and
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
  `UnitXcSelectorRequiresLibxc` is a negative test in a build without libXC;
* XCF-08 semantic tests cover all four `control%nsp` modes, the independent
  two-channel radial contract, explicit list composition, TXC=6/7 density
  guards, and kinetic/2D/MGGA/missing-VXC failures.

Material calculations are smoke tests only.  No Fe moment, material-specific
threshold, or legacy XC formula is used as the primary certification oracle.
