# XC selector and libXC mapping contract

This document defines the exchange-correlation selector contract implemented
in `source/xc.f90`. It is about selection and provenance; it does not redefine
any legacy XC formula.

## Two independent namespaces

`TXC` is the RS-LMTO user-facing selector:

| TXC range | Meaning |
| --- | --- |
| `1–99` | Historical/internal RS-LMTO selectors |
| `100–199` | Predefined libXC XC bundles |
| `>=1000` | Direct native libXC request, `libXC_ID = TXC - 1000` |

`libxc_func_id(:)` contains native libXC identifiers only. Its value `17`, for
example, means native libXC `XC_LDA_C_VBH`; it never means legacy `TXC=17`.
Likewise, native libXC ID `101` is `XC_GGA_X_PBE`, while `TXC=101` is an
RS-LMTO bundle selector that expands to native IDs `[1,17]`.

## Backend dispatch

The dispatch invariant is:

* `TXC=1–99` uses the historical RS-LMTO implementation.
* `TXC=100–199` uses the libXC backend when the bundle is defined.
* `TXC>=1000` uses the libXC backend with exactly one native ID, `TXC-1000`.

An explicit libXC selector is never passed to legacy `XCPOT`. Unsupported
100-series bundles, invalid native IDs, and libXC selectors in a build without
libXC terminate with a descriptive error. No historical functional is silently
substituted.

## Historical RS-LMTO selectors

These are the production implementations in the legacy XC code. The native
IDs in the table are reference mappings only; they are not placed in the
active `libxc_func_id(:)` array for a legacy selector.

| TXC | Backend | RS-LMTO name | Reference native libXC IDs | Mapping quality | Notes |
| ---: | --- | --- | --- | --- | --- |
| 1 | legacy | Barth–Hedin | `[1,17]` | `REFERENCE_EQUIVALENT` | `XC_LDA_X + XC_LDA_C_VBH`; useful for direct comparison, but not an assertion of identical implementation details |
| 2 | legacy | Slater X-alpha | `[1]` | `APPROXIMATE_ANALOGUE` | Exchange-only reference; the historical X-alpha parameterization is not silently replaced |
| 3 | legacy | Barth–Hedin–Janak | `[1,17]` | `APPROXIMATE_ANALOGUE` | BHJ uses a different correlation parameter set from TXC=1 |
| 4 | legacy | Vosko–Wilk–Nusair | `[1,7]` | `REFERENCE_EQUIVALENT` | `XC_LDA_X + XC_LDA_C_VWN` reference |
| 5 | legacy | Perdew–Burke–Ernzerhof 96 LDA | `[1,12]` | `REFERENCE_EQUIVALENT` | `XC_LDA_X + XC_LDA_C_PW` reference |
| 6 | legacy | Wigner exchange | — | `NO_EQUIVALENT` | No validated ASA-equivalent libXC mapping is selected |
| 7 | legacy | Perdew–Zunger | `[1,9]` | `REFERENCE_EQUIVALENT_UNPOLARIZED` | `XC_LDA_X + XC_LDA_C_PZ` reference; historical TXC=7 is unpolarized-only |
| 8 | legacy | Perdew–Burke–Ernzerhof 96 GGA | `[101,130]` | `REFERENCE_EQUIVALENT` | `XC_GGA_X_PBE + XC_GGA_C_PBE` reference |
| 9 | legacy | Local Airy gas | — | `NO_EQUIVALENT` | No validated mapping is selected |
| 11 | legacy | Barth–Hedin ASW variant | `[1,17]` | `APPROXIMATE_ANALOGUE` | ASW uses a distinct parameter set; the VBH reference is not mathematical identity |

The three Barth–Hedin-family rows are intentionally separate. TXC=1, TXC=3,
and TXC=11 must not be collapsed into one claim merely because native
libXC `XC_LDA_C_VBH` is a useful comparison reference.

## Predefined libXC XC bundles

These bundles are the authoritative predefined entries in
`setup_libxc_functional_ids`. Every number shown is a native libXC ID.

| TXC | Backend | Functional names | Native libXC IDs | Mapping quality | Notes |
| ---: | --- | --- | --- | --- | --- |
| 101 | libXC | LDA exchange + von Barth–Hedin correlation | `[1,17]` | `REFERENCE_EQUIVALENT` | Explicit comparison route for legacy TXC=1 |
| 102 | libXC | LDA exchange + Gombas correlation | `[1,24]` | `APPROXIMATE_ANALOGUE` | Native ID 24 is `XC_LDA_C_GOMBAS` |
| 103 | libXC | LDA exchange | `[1]` | `APPROXIMATE_ANALOGUE` | Exchange-only route; it does not add correlation |
| 104 | libXC | LDA exchange + Perdew–Zunger correlation | `[1,9]` | `REFERENCE_EQUIVALENT_UNPOLARIZED` | Reference counterpart to TXC=7 in the unpolarized limit; libXC supports spin polarization |
| 105 | libXC | LDA exchange + Perdew–Wang correlation | `[1,12]` | `REFERENCE_EQUIVALENT` | Explicit reference route for legacy TXC=5 |
| 106 | libXC | LDA exchange + Vosko–Wilk–Nusair correlation | `[1,7]` | `REFERENCE_EQUIVALENT` | Explicit reference route for legacy TXC=4 |
| 107 | libXC | LDA exchange + Gunnarsson–Lundqvist correlation | `[1,5]` | `APPROXIMATE_ANALOGUE` | Native `XC_LDA_C_GL` reference |
| 108 | libXC | PBE exchange + PBE correlation | `[101,130]` | `REFERENCE_EQUIVALENT` | Explicit reference route for legacy TXC=8 |
| 109 | libXC | RPBE exchange + PBE correlation | `[117,130]` | `APPROXIMATE_ANALOGUE` | No legacy RS-LMTO RPBE implementation is claimed |

`TXC=100` and currently undefined values in the `100–199` range are still
libXC-style selectors, but they fail as unsupported bundles instead of
falling through to legacy code.

## Direct libXC selectors

`TXC=1000+ID` selects exactly native libXC functional ID `ID`; here `ID` is a
native libXC identifier, not another RS-LMTO `TXC` value. For example,
`TXC=1001` selects `[1]`, native `XC_LDA_X`, and remains exchange-only. The
selector layer does not infer or add a correlation/exchange partner from a
functional name. If a caller wants a pair, it must use a predefined bundle or
an explicitly defined pair in the interface.

Direct native requests do not claim a legacy-to-libXC equivalence; their
metadata quality is `NO_EQUIVALENT`. The selector, backend, native ID, and
native functional name provide the provenance for the direct request.

The selected native functional must be an ASA-compatible LDA or GGA. Invalid,
meta-GGA, hybrid, or otherwise unsupported native functionals fail during
initialization.

## Reference mappings and equivalence language

The quality labels mean:

* `EXACT_EQUIVALENT`: identity established for the stated implementation and
  conventions. No legacy-to-libXC row currently makes this strongest claim.
* `REFERENCE_EQUIVALENT`: same named/reference parametrization is useful for a
  controlled comparison; this does not promise bitwise or mathematical
  identity of legacy and libXC implementations.
* `REFERENCE_EQUIVALENT_UNPOLARIZED`: the same named/reference parametrization
  is a controlled counterpart only in the unpolarized limit. The historical
  implementation does not provide the spin-polarized extension, even when the
  libXC counterpart does.
* `APPROXIMATE_ANALOGUE`: a related or nearest practical reference, with a
  known parameterization or functional distinction.
* `NO_EQUIVALENT`: no validated mapping is available for this purpose.

A reference mapping is not itself a physics validation result.

In particular, TXC=104 is the libXC reference counterpart to the historical
TXC=7 Perdew–Zunger parametrization in the unpolarized limit, subject to the
established normalization and unit conventions. Legacy TXC=7 does not provide
a spin-polarized implementation.

## Runtime provenance and tests

XC initialization prints the selector, backend, native IDs when applicable,
native libXC names, and mapping quality. This makes it possible to identify
which backend produced a potential from the run log. An explicit native-ID
list is reported as `XC mode: custom libXC composition`, followed by an
informational provenance message stating that the requested components are
summed exactly as supplied. Valid but unusual lists retain their warning;
valid ordinary custom lists remain non-fatal and informative.

Focused selector coverage is provided by:

* `tests/unit/test_xc_selector_semantics.f90`, covering legacy selectors,
  bundles 101/105/108, and direct selector 1001;
* `tests/unit/test_xc_selector_requires_libxc.f90`, registered in a build
  without libXC and expecting a construction-time failure for a libXC-only
  selector; and
* `tests/unit/test_libxc_xc_baseline.f90`, which compares legacy TXC=5 with
  the explicitly selected libXC bundle TXC=105 without changing either kernel.
