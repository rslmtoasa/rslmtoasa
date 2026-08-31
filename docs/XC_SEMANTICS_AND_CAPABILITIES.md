# XC semantics and libXC capabilities

This is the authoritative XCF-08 convention document for the exchange-
correlation path. It covers selection, composition, capability validation,
evaluation routing, and spin semantics. The historical LDA/GGA formulas are
not redefined here.

## 1. Three distinct XC concepts

The user-facing `TXC` value is an RS-LMTO selector. A native libXC component
is one libXC object, identified by a native integer ID such as `1`
(`XC_LDA_X`) or `101` (`XC_GGA_X_PBE`). An XC bundle, also called an XC
composition, is a list of native components whose energy and potential
contributions are summed.

These concepts use separate namespaces:

| `TXC` | Meaning |
| --- | --- |
| `1–99` | historical/internal RS-LMTO functional |
| `100–199` | predefined libXC XC bundle |
| `>=1000` | direct single native component, with `native_ID = TXC - 1000` |

For example, `TXC=101` is a predefined bundle `[1,17]`; it is not native
libXC ID 101. Conversely, `TXC=1101` selects exactly native ID 101. A direct
selection never adds an exchange or correlation partner:

```text
TXC=1001  -> [1]    XC_LDA_X only
TXC=1012  -> [12]   XC_LDA_C_PW only
TXC=1101  -> [101]  XC_GGA_X_PBE only
TXC=1130  -> [130]  XC_GGA_C_PBE only
```

The complete predefined table is maintained in
[`XC_SELECTOR_AND_LIBXC_MAPPING.md`](XC_SELECTOR_AND_LIBXC_MAPPING.md).

## 2. Predefined bundles and arbitrary compositions

The 100-series convenience interface is retained because the traditional
input model presents one XC choice while libXC commonly separates exchange
and correlation. The current predefined bundles are:

| `TXC` | Native IDs | Composition |
| ---: | --- | --- |
| 101 | `[1,17]` | LDA exchange + von Barth–Hedin correlation |
| 102 | `[1,24]` | LDA exchange + Gombas correlation |
| 103 | `[1]` | LDA exchange only |
| 104 | `[1,9]` | LDA exchange + Perdew–Zunger correlation |
| 105 | `[1,12]` | LDA exchange + Perdew–Wang correlation |
| 106 | `[1,7]` | LDA exchange + Vosko–Wilk–Nusair correlation |
| 107 | `[1,5]` | LDA exchange + Gunnarsson–Lundqvist correlation |
| 108 | `[101,130]` | PBE exchange + PBE correlation |
| 109 | `[117,130]` | RPBE exchange + PBE correlation |

Internally, `libxc_func_id(:)` stores an arbitrary-length list of native IDs.
There is no one-component/two-component assumption and no requirement that a
particular kind occupy a particular position. The setter
`xc%set_libxc_component_ids([id1,id2,...])` is the supported internal
extension point for a future explicit namelist such as
`libxc_components = 101,130`; no additional magic `TXC` ranges are used.
The setter validates and rebuilds all retained metadata, marks the selection
as an explicit composition, and does not claim a legacy mapping equivalence.
Components are never silently reordered, substituted, or completed.

Exchange-only, correlation-only, multiple-component, and unusual but valid
compositions are allowed when every component is individually compatible.

## 3. Initialization metadata and lifecycle

During XC initialization, each requested native component is temporarily
initialized through libXC. Before that object is ended, the code copies a
descriptor containing:

* native ID and human-readable name;
* family and kind (`exchange`, `correlation`, combined `exchange-correlation`,
  or `kinetic`);
* raw flags and derived dimensionality (`1D`, `2D`, `3D`);
* required ingredients;
* `EXC`, `VXC`, `FXC`, `KXC`, and `LXC` availability; and
* development, Laplacian, and nonlocal-ingredient markers.

The active `xc` object retains copied descriptors and parallel native-ID,
family, and kind arrays. It does not retain libXC objects or metadata pointers
after `xc_f03_func_end`. Evaluation creates short-lived native objects only
to evaluate the already validated IDs; it does not query native metadata for
routing or validation.

The initialization report prints the selector, backend, mode, every component
name/ID/kind/family/dimension, the complete classification, route, and
mapping quality once. It is provenance, not per-iteration output.

## 4. Compatibility policy

The policy has three levels:

* **FATAL:** reject before the SCF loop when the request is incompatible with
  the 3D ASA evaluator. This includes kinetic functionals, explicit 1D/2D or
  non-3D functionals, unsupported families such as MGGA/hybrid/LCA/OEP,
  Laplacian- or nonlocal-ingredient requirements, missing `EXC` or `VXC`,
  invalid native IDs, undefined 100-series selectors, and libXC selectors in
  a no-libXC build.
* **WARNING:** accept without changing the request for exchange-only or
  correlation-only compositions, more than one exchange/correlation
  component, a combined component plus additional components, unconventional
  valid X+C families, or libXC components marked as development.
* **OK:** proceed for compatible predefined bundles, ordinary compatible
  compositions, and intentionally selected single components.

The code never turns a direct component into a bundle and never silently adds
the missing partner. Warnings explicitly state that all requested components
will be summed exactly as requested.

The supported family boundary is:

| Family | Status | Required data |
| --- | --- | --- |
| LDA | supported | two local spin-density eigenchannels |
| GGA | supported | local spin densities and radial gradients |
| MGGA, hybrids, LCA, OEP, orbital/current/tau-dependent forms | rejected | ingredients unavailable in the present ASA path |

Kinetic components are rejected independently of family because a kinetic
functional must not contribute `delta T/delta n` to the KS XC potential.

## 5. Composition classification and evaluation route

After all components have been inspected, the code counts
`n_exchange`, `n_correlation`, `n_xc_combined`, and `n_kinetic`, and records
whether any component requires gradients. The route is determined by the
complete list:

```text
if any active component is GGA:
    radial GGA-capable evaluator
else:
    pointwise LDA evaluator
```

Thus LDA exchange plus LDA correlation is pointwise LDA, while GGA exchange
plus LDA correlation is radial GGA. In the radial evaluator, LDA components
contribute `exc` and `vrho`; GGA components contribute those terms plus
`vsigma`. The active list is summed exactly, without assumptions about its
length or ordering.

The libXC boundary uses standard channel order `[up, down]`. The historical
`XCPOT` boundary remains `RHO1/V1 = down` and `RHO2/V2 = up`; the wrapper maps
between these conventions. libXC Hartree outputs are converted to the
internal Rydberg convention once (`1 Ha = 2 Ry`).

## 6. Global spin mode versus local XC channels

The central distinction is:

```text
control%nsp != number of XC spin channels
```

`control%nsp` is the global electronic-structure mode:

| `control%nsp` | Collinear | Noncollinear | SOC |
| ---: | --- | --- | --- |
| 1 | yes | no | no |
| 2 | yes | no | yes |
| 3 | no | yes | no |
| 4 | no | yes | yes |

Use the semantic queries `is_collinear()`, `is_noncollinear()`, `has_soc()`,
`uses_spinor_representation()`, and `is_spin_polarized_mode()` for global
questions. The last query is a mode capability only; it does not claim that
the current density is magnetized.

The atomic radial/XC path uses two local density channels in its ordinary
spin-density representation, including the `control%nsp=1` scalar-relativistic
collinear mode. The explicit local count is named/documented as
`n_radial_spin_channels`; it is not inferred from the global mode integer.

For noncollinear calculations, the established production design uses the
local eigenchannels of the spin-density matrix,

\[
n_\pm(r)=\frac{1}{2}\left[n(r)\pm|\mathbf m(r)|\right],
\]

so a noncollinear global mode still supplies two local channels to LSDA/GGA.
Two local XC channels do not imply global collinearity, and `nsp > 2` never
means four XC channels. XCF-08 verifies that these remain two local XC
channels independent of the global `control%nsp` mode and removes the old
global/local `nsp` confusion. It is not a complete new end-to-end validation
of the reconstruction

\[
\rho_{\alpha\beta}\rightarrow n_\pm\rightarrow V^\pm_{\rm xc}
\rightarrow V^{\rm xc}_{\alpha\beta}
\]

for arbitrary noncollinear densities; that remains part of the broader
noncollinear magnetic validation suite.

## 7. Legacy TXC=6 and TXC=7

The historical Wigner (`TXC=6`) and Perdew–Zunger (`TXC=7`) implementations
are retained and explicitly marked `supports_spin_polarized_density = false`.
They provide equal spin potentials and must not be used to silently evaluate
a polarized density.

The runtime guard examines the actual density passed to `XCPOT`, not
`control%nsp`. For two local channels it checks

\[
\Delta n(r)=n_\uparrow(r)-n_\downarrow(r).
\]

An equal density is accepted to the documented tolerance, even in a mode that
permits magnetism. An appreciably unequal density terminates with a diagnostic
directing the user to a spin-capable functional or an unpolarized calculation.
This prevents a global mode value from being mistaken for the current
physical polarization and prevents silently using `V_up=V_down` for a
polarized system. The comparison is an absolute channel difference,
`abs(RHO1-RHO2)`, against `1e-10*max(abs(RHO),1e-20)`: the supplied `RHO` is
the total-density normalization and the second term only avoids a zero-scale
comparison.

## 8. Tests and evidence

The focused XCF-08 coverage is:

* `UnitControlSpinSemantics`: all four global-mode helper combinations and
  the independent two-channel radial contract;
* `UnitXcCompositionSemantics`: direct no-partner selection, retained
  arbitrary list order, complete-list GGA routing, warning classification, and
  exact requested-component sums;
* `UnitXcLegacyUnpolarizedDensity`: equal-density acceptance for TXC=6/7 with
  `control%nsp=1`;
* `UnitXcLegacyPolarizedRequiresFatal` and its TXC=7 invocation:
  unequal-density rejection for TXC=6/7 with `control%nsp=1`;
* `UnitXcLegacyPolarizationTolerance`: exact equality and `0.1*epsilon_spin`
  and `0.9*epsilon_spin` noise are accepted, while `1.1*epsilon_spin`,
  `10*epsilon_spin`, and a clearly polarized density are rejected for both
  TXC=6 and TXC=7;
* `UnitXcCapabilityRequiresFatal`, `UnitXcKineticRequiresFatal`,
  `UnitXcWrongDimensionalityRequiresFatal`, and
  `UnitXcMissingVxcRequiresFatal`: representative capability failures for
  native IDs 263, 50, 19, and 160; and
* `UnitLibxcProductionContract` and the established GGA/LDA tests: metadata
  lifecycle, exact component accumulation, units, radial flux, endpoint and
  origin behavior, and prior numerical formulas.

For arbitrary valid user-defined native lists, the implementation preserves
the supplied ordering and components, emits custom-composition provenance,
and only warns for scientifically unusual combinations. Architecturally
invalid components remain fatal.

The existing validated numerical routes remain unchanged for previously valid
legacy and libXC selections. Material-specific magnetic campaigns are outside
this semantic/capability scope.

Related authoritative documents:

* [`LIBXC_PRODUCTION_CONTRACT.md`](LIBXC_PRODUCTION_CONTRACT.md) — interface,
  units, radial formulation, and production evidence;
* [`XC_SELECTOR_AND_LIBXC_MAPPING.md`](XC_SELECTOR_AND_LIBXC_MAPPING.md) —
  selector table and reference-mapping labels; and
* [`XC_GENERAL_CLOSEOUT.md`](XC_GENERAL_CLOSEOUT.md) — prior numerical
  correctness closeout and regression evidence.
