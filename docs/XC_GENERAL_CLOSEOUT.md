# XCF-09 — General XC correctness closeout and semantic hardening

## Decision

The audited XC scope is accepted for production use:

* legacy RS-LMTO LDA and GGA routes remain the selected implementation for
  TXC=1–9, subject to the selector's declared spin restrictions;
* predefined libXC XC bundles TXC=101–109 and direct native selectors
  TXC=1000+ID select only validated LDA/GGA functionals;
* unsupported selectors, unsupported libXC families, and libXC requests in a
  no-libXC executable fail at construction/evaluation boundaries; and
* no libXC request falls back to legacy `XCPOT`.

No failed or unexplained result remains in the acceptance matrix. Converged
fcc-Fe physics and the GBT/Phase-II campaign are intentionally deferred; the
bounded magnetic-SCF findings remain documented in
[`XC_MAGNETIC_SCF_CLOSEOUT.md`](XC_MAGNETIC_SCF_CLOSEOUT.md).

The general RS-LMTO XC layer is closed for the supported scope: legacy
LDA/GGA paths at their documented validation level, and compatible 3D libXC
LDA/GGA components and bundles providing the required EXC/VXC capabilities.
This does not mean that every libXC functional is supported: the scope is
limited to compatible 3D LDA/GGA components whose required ingredients and
derivative capabilities are available to the ASA radial path.
Unsupported functional families or ingredients fail explicitly. Scientifically
unusual but mathematically valid libXC compositions are preserved exactly and
reported to the user. Material-specific magnetic behavior is not part of this
closeout.

## Interruption recovery

The interrupted worktree contained a large uncommitted replacement of
`source/xc.f90`. It was restored to `HEAD` as explicitly requested before the
closeout was resumed. The restored source then exposed one existing test-fixture
bug under the freshly configured Debug build: Fortran evaluates both operands
of `merge`, so the origin test indexed the absent second spin column for
`nsp=1`. The test now uses an explicitly guarded array assignment in
`test_vxc0sp_gga_origin.f90`.

The resumed source changes are narrow:

1. unknown legacy selectors now fail for both `nsp=1` and `nsp=2` instead of
   falling through to Barth–Hedin;
2. direct legacy `PBEGGA` and `LAGGGA` return exact zero for zero total
   density before any division by density; and
3. regressions cover the unknown-selector failure and zero-density behavior,
   plus the libXC GGA polarization/tail matrix.

## Final XC dispatch architecture

The production flow is:

```text
control%txc
    |
    v
xc constructor
    |
    +-- legacy TXC=1..99 ------------> XCPOT / PBEGGA / LAGGGA
    |
    +-- bundle TXC=101..199 ------------+
    |                                   v
    +-- direct TXC=1000+ID -------> init_libxc
                                        |
                                        +-- validate every native component
                                        +-- record family/kind/name metadata
                                        +-- LDA -> pointwise libXC wrapper
                                        +-- GGA -> radial libXC flux helper
                                        |
                                        v
                                   VXC0SP / projections
```

`init_libxc` is the single lifecycle boundary. It clears previous metadata,
builds the active native-ID list, validates every component, copies metadata
while each temporary libXC object is alive, and then ends the object. A mixed
LDA+GGA list is routed as aggregate GGA so its LDA pointwise terms and GGA
radial flux terms are both retained.

The ordinary production GGA path in `VXC0SP` reconstructs the physical density
from stored spherical radial density, obtains `dn/dr` and raw `d2n/dr2` from
the radial derivative routine, and passes raw derivatives to legacy `PBEGGA`.
At the regular origin, `PBEGGA` applies the spherical limit
`laplacian(n)=3*d2n/dr2` exactly once. The origin handoff and channel mapping
are covered by `UnitVxc0spGgaOrigin`.

## Selector and backend contract

| selector namespace | Meaning | Failure behavior |
| --- | --- | --- |
| `1–99` | historical RS-LMTO implementation | unknown values are fatal |
| `100–199` | predefined libXC XC bundle | undefined bundles are fatal |
| `>=1000` | exactly native libXC ID `TXC-1000` | incompatible/invalid family is fatal |
| `200–999` | no defined namespace | fatal |

The validated predefined mappings are:

| TXC | Native IDs | Route | Mapping status |
| ---: | --- | --- | --- |
| 101 | `[1,17]` | LDA X + VBH C | reference equivalent to TXC=1 |
| 102 | `[1,24]` | LDA X + Gombas C | approximate analogue |
| 103 | `[1]` | LDA X | approximate analogue |
| 104 | `[1,9]` | LDA X + PZ C | `REFERENCE_EQUIVALENT_UNPOLARIZED` to TXC=7 |
| 105 | `[1,12]` | LDA X + PW C | reference equivalent to TXC=5 |
| 106 | `[1,7]` | LDA X + VWN C | reference equivalent to TXC=4 |
| 107 | `[1,5]` | LDA X + GL C | approximate analogue |
| 108 | `[101,130]` | PBE X + PBE C | reference equivalent to TXC=8 |
| 109 | `[117,130]` | RPBE X + PBE C | approximate analogue |

The historical selectors do not all have the same evidence depth:

| Functional/path | Status |
| --- | --- |
| TXC=1 Barth–Hedin | VALIDATED / REPAIRED |
| TXC=2 Slater/X-alpha | RETAINED / analytic exchange-covered where applicable |
| TXC=3 BHJ | RETAINED / historical variant |
| TXC=4 VWN | CHARACTERIZED / VALIDATED TO HISTORICAL SCOPE |
| TXC=5 PW/PBE-LDA | CHARACTERIZED / VALIDATED TO HISTORICAL SCOPE |
| TXC=6 Wigner | LIMITED: UNPOLARIZED ONLY |
| TXC=7 PZ | LIMITED: UNPOLARIZED ONLY |
| TXC=8 PBE GGA | VALIDATED / REPAIRED |
| TXC=9 LAG | LEGACY RETAINED / NOT COMPREHENSIVELY CERTIFIED |

The rows above describe the historical RS-LMTO paths. The predefined
libXC bundles TXC=101–109 are separate explicit native-ID lists; arbitrary
compatible native libXC compositions are a third category and are preserved
and reported exactly as requested, without being presented as predefined
bundle validation.

Legacy-to-libXC labels describe controlled reference comparisons. They do not
claim bitwise identity. TXC=104 is the libXC reference counterpart to the
historical TXC=7 parametrization only in the unpolarized limit, subject to the
established normalization and unit conventions; legacy TXC=7 does not provide
a spin-polarized implementation. Direct native requests select exactly one
requested component; no exchange/correlation partner is inferred.

XCF-08 also makes the global/local spin distinction explicit. `control%nsp`
selects the global collinear/noncollinear and SR/FR(SOC) mode; it is not the
number of XC channels. The atomic radial/XC path uses two local channels even
when `control%nsp=1`. TXC=6/7 now inspect the actual density at `XCPOT`: equal
channels are accepted, while an appreciably polarized density is fatal.

The production noncollinear XC path is designed around the local eigenchannels
`n_±=(n±|m|)/2`. XCF-08 validates that these remain two local XC channels
independent of the global `control%nsp` mode and validates the global/local
semantic separation. It does not constitute a complete new end-to-end test of
`rho_(alpha beta) -> n_± -> V_xc^± -> V_xc_(alpha beta)` for arbitrary
noncollinear densities; that remains part of the broader noncollinear magnetic
validation suite.

## Density, channel, unit, and radial conventions

* At the historical `XCPOT` boundary, `RHO1/V1` are down-spin and `RHO2/V2`
  are up-spin. The libXC boundary receives standard `[up, down]` order and
  maps the result back to the historical boundary.
* libXC energies and derivatives are converted from Hartree to the internal
  Rydberg convention exactly once.
* GGA sigma is `[|grad up|², grad up·grad down, |grad down|²]`.
  The radial helper forms the corresponding spin fluxes and applies the
  spherical radial divergence; it does not treat a GGA as a pointwise LDA.
* `LIBXC_DENSITY_FLOOR=1e-20` is applied only to libXC input channels. It is
  never applied to RS-LMTO density arrays or quadrature weights.
* Exact zero total density returns exact zero. A positive total density with a
  zero channel is valid and is regularized only at the libXC boundary.

## Acceptance evidence

### Fixed-density LDA oracle

Command:

```bash
python3 tests/xc_reconciliation/run_xc_lda_reconciliation.py \
  --binary build/bin/UnitXcLdaReconciliation \
  --output-dir /tmp/xcf07_xc_lda_oracle --no-plots
```

The sweep contains 378 comparison rows: 351 regular rows and 27 exact fully
polarized endpoint rows. It includes `r_s=0.5,1,1.5,2,2.5,3,4,5,6`, the
polarizations `zeta=0, 1e-4, 1e-3, .01, .05, .1, .25, .5, .8, .9, .95,
.99, .999999`, and exact one-channel endpoints.

| comparison | max regular `|exc|` | max regular `|Vup|` | max regular `|Vdown|` | max regular `|Bxc|` | classification |
| --- | ---: | ---: | ---: | ---: | --- |
| TXC=1 / 101, BH/VBH | 8.88e-9 Ry | 1.18e-8 Ry | 9.16e-9 Ry | 6.46e-9 Ry | reconciled equivalent |
| TXC=5 / 105, PW/PW | 7.44e-7 Ry | 8.94e-7 Ry | 1.70e-5 Ry | 8.95e-6 Ry | expected variant difference |
| TXC=4 / 106, VWN/VWN | 5.08e-4 Ry | 8.36e-4 Ry | 1.17e-2 Ry | 5.83e-3 Ry | expected variant difference |

The analytic spin-polarized exchange oracle has maximum absolute error
`1.5543122344752192e-15` and maximum relative error
`1.6896093233077543e-14`. The finite fully polarized endpoint differences are
explained by the documented libXC input floor and remain finite; they are not
treated as unexplained failures. The largest regular fixed-density energy
gradient residuals are `8.45e-10` (BH), `9.48e-10` (PW), and `4.48e-9`
(VWN).

### GGA variational and cross-route checks

| check | result |
| --- | --- |
| `UnitLegacyPbeGga` fixed-density functional derivative | PASS; final residual `9.05e-10` |
| `UnitLibxcGgaRadial` radial functional derivative | PASS; final residual `9.05e-10` |
| `UnitPbeGgaComparison` TXC=8 vs TXC=108 | PASS; energy-density difference `8.88e-16`, scalar-potential difference `5.06e-6` Ry, magnetic-potential difference `5.30e-6` Ry |
| `UnitVxc0spGgaOrigin` | PASS for polarized and unpolarized production origin handoffs |
| `UnitLibxcProductionContract` | PASS; bundles, direct IDs, metadata, lifecycle, units, mixed components, radial flux, floor, ζ matrix, and zero tail |
| `UnitXcLegacyPolarizationTolerance` | PASS; exact equality, `0.1 epsilon_spin`, and `0.9 epsilon_spin` accepted; `1.1 epsilon_spin`, `10 epsilon_spin`, and clear polarization rejected for TXC=6 and TXC=7 |

The production-contract ζ matrix explicitly covers `0, .5, .9, .99, 1`, a
positive sub-floor total density, and exact zero density.

### Selector rejection and build variants

The final focused matrices are:

```bash
ctest --test-dir build --output-on-failure \
  -R '^(UnitXcSelectorSemantics|UnitXcUnknownSelectorRequiresFatal|UnitXcLdaLegacyKernel|UnitRadialGga|UnitLegacyPbeGga|UnitVxc0spGgaOrigin|UnitLibxcXcBaseline|UnitLibxcGgaRadial|UnitLibxcProductionContract|UnitXcUnsupportedFamilyRequiresFatal|UnitPbeGgaComparison|UnitXcLdaReconciliation)$'

ctest --test-dir build-xcf06-no-libxc --output-on-failure \
  -R '^(UnitXcSelectorSemantics|UnitXcUnknownSelectorRequiresFatal|UnitXcLdaLegacyKernel|UnitRadialGga|UnitLegacyPbeGga|UnitVxc0spGgaOrigin|UnitXcSelectorRequiresLibxc)$'
```

Results: libXC-enabled `12/12` passed; no-libXC `7/7` passed. The expected
failure tests prove that unknown legacy selectors, unsupported native MGGA
families, and libXC-only selectors in a no-libXC executable terminate rather
than silently selecting a legacy path.

Those earlier counts covered the pre-XCF-09 matrix. The authoritative XCF-09
execution evidence is below; expected-fatal tests are counted as passed by
CTest through their `WILL_FAIL` properties.

### XCF-08/XCF-09 execution evidence

Both builds used GNU Fortran 13.3.0, OpenMP enabled, MPI and CUDA disabled,
and standalone unit tests enabled:

```bash
cmake -S . -B build-xcf06-debug -DENABLE_LIBXC=ON -DRUN_UNIT_TESTS=ON -DCMAKE_BUILD_TYPE=Debug
cmake -S . -B build-xcf06-no-libxc -DENABLE_LIBXC=OFF -DRUN_UNIT_TESTS=ON -DCMAKE_BUILD_TYPE=Release
```

The enabled build found libXC/libXCF03 version 5.2.3 from `/usr/include` and
`/usr/lib/x86_64-linux-gnu`. The focused semantic/capability command was:

```bash
ctest --test-dir build-xcf06-debug --output-on-failure -j2 -R '^(UnitControlSpinSemantics|UnitXcSelectorSemantics|UnitXcCompositionSemantics|UnitXcLegacyUnpolarizedDensity|UnitXcLegacyPolarizedRequiresFatal|UnitXcLegacyPolarizedPzRequiresFatal|UnitLibxcProductionContract|UnitXcUnsupportedFamilyRequiresFatal|UnitXcCapabilityRequiresFatal|UnitXcKineticRequiresFatal|UnitXcWrongDimensionalityRequiresFatal|UnitXcMissingVxcRequiresFatal|UnitXcLegacyPolarizationTolerance(_txc[67]_case[0-5])?)$'
```

Result: `25` tests run, `25` passed, `0` failed. This includes the requested
XCF-08 semantic/capability tests, the separate TXC=7 fatal variant, and the
13 tolerance registrations (default accepted-case sweep plus six cases for
each of TXC=6 and TXC=7).

The corresponding no-libXC command was:

```bash
ctest --test-dir build-xcf06-no-libxc --output-on-failure -j2 -R '^(UnitControlSpinSemantics|UnitXcSelectorSemantics|UnitXcLegacyUnpolarizedDensity|UnitXcLegacyPolarizedRequiresFatal|UnitXcLegacyPolarizedPzRequiresFatal|UnitXcLegacyPolarizationTolerance(_txc[67]_case[0-5])?|UnitXcSelectorRequiresLibxc|UnitXcUnknownSelectorRequiresFatal|UnitXcLdaLegacyKernel|UnitRadialGga|UnitLegacyPbeGga|UnitVxc0spGgaOrigin)$'
```

Result: `24` tests run, `24` passed, `0` failed. The libXC-specific
composition, metadata, and native capability tests are intentionally absent
from this configuration; the no-libXC selector rejection and all applicable
legacy/tolerance sentinels are executed.

The broader enabled-build XC sentinel command was:

```bash
ctest --test-dir build-xcf06-debug --output-on-failure -j2 -R '^(UnitXcLdaLegacyKernel|UnitLibxcXcBaseline|UnitXcLdaReconciliation|UnitRadialGga|UnitLegacyPbeGga|UnitVxc0spGgaOrigin|UnitLibxcGgaRadial|UnitPbeGgaComparison|UnitLibxcProductionContract|UnitXcSelectorSemantics|UnitXcUnknownSelectorRequiresFatal|UnitXcUnsupportedFamilyRequiresFatal|UnitXcKineticRequiresFatal|UnitXcWrongDimensionalityRequiresFatal|UnitXcMissingVxcRequiresFatal|UnitXcCompositionSemantics|UnitControlSpinSemantics|UnitXcLegacyUnpolarizedDensity|UnitXcLegacyPolarizedRequiresFatal|UnitXcLegacyPolarizedPzRequiresFatal|UnitXcLegacyPolarizationTolerance(_txc[67]_case[0-5])?)$'
```

Result: `33` tests run, `33` passed, `0` failed. Previously valid XC
selections remained numerically unchanged in this sentinel run.

### Full production smoke

```bash
ctest --test-dir build --output-on-failure -j3 \
  -R '^(Example_bulk_fccCu_chebyshev|Example_bulk_bccFe_nsp2_block|Example_exchange_bccFe)$'
```

All three passed: nonmagnetic metallic fcc-Cu, robust magnetic bcc-Fe, and
the bcc-Fe exchange/post-processing path. A separate one-step scratch run
derived from `tests/scf/cases/bulk/bccFe`, with `nsp=2`, `txc=8`,
`recur='block'`, and `lld=20`, selected the legacy GGA route and completed
without error, fatal, NaN, or Infinity. Its expected one-step “Not
converged” diagnostic is a structural smoke result, not a convergence claim.

## Deferred scope

The following are deliberately not acceptance blockers for this closeout:

* converged high-spin fcc-Fe ordinary q=0 physics and any Phase-II material
  campaign;
* GBT validation and noncollinear extensions;
* unsupported libXC meta-GGA, hybrid, orbital-dependent, and
  long-range-corrected families; and
* unvalidated legacy-to-libXC mapping claims beyond the explicit bundles;
  compatible arbitrary native compositions remain supported subject to the
  capability contract.

No threshold retuning, material-specific XC tuning, implicit selector pairing,
or fallback behavior was introduced.

## Source and regression index

* [`source/xc.f90`](../source/xc.f90): selector construction, libXC lifecycle,
  LDA wrapper, radial GGA helper, and legacy zero-density guards.
* [`source/self.f90`](../source/self.f90): production `VXC0SP` density,
  derivative, origin, and channel handoff.
* [`source/xc_radial.f90`](../source/xc_radial.f90): physical radial
  derivatives and spherical flux divergence.
* [`tests/unit/test_libxc_production_contract.f90`](../tests/unit/test_libxc_production_contract.f90)
  and [`tests/unit/test_libxc_gga_radial.f90`](../tests/unit/test_libxc_gga_radial.f90):
  libXC contract and variational coverage.
* [`tests/unit/test_legacy_pbe_gga.f90`](../tests/unit/test_legacy_pbe_gga.f90),
  [`tests/unit/test_pbe_gga_comparison.f90`](../tests/unit/test_pbe_gga_comparison.f90),
  and [`tests/unit/test_vxc0sp_gga_origin.f90`](../tests/unit/test_vxc0sp_gga_origin.f90):
  legacy and production GGA coverage.
* [`tests/unit/test_xc_unknown_selector_requires_fatal.f90`](../tests/unit/test_xc_unknown_selector_requires_fatal.f90):
  selector rejection coverage.
* [`tests/unit/test_control_spin_semantics.f90`](../tests/unit/test_control_spin_semantics.f90),
  [`tests/unit/test_xc_composition_semantics.f90`](../tests/unit/test_xc_composition_semantics.f90),
  and the TXC=6/7/capability negative tests: XCF-08 semantic and
  pre-SCF compatibility coverage.
* [`docs/XC_LDA_RECONCILIATION.md`](XC_LDA_RECONCILIATION.md),
  [`docs/XC_GGA_RADIAL_RECONCILIATION.md`](XC_GGA_RADIAL_RECONCILIATION.md),
  and [`docs/LIBXC_PRODUCTION_CONTRACT.md`](LIBXC_PRODUCTION_CONTRACT.md):
  detailed prior evidence and contracts.
