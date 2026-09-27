# VAL-03 — libXC SCF equivalence

Date: 2026-08-16
Binary: `/tmp/rslmto-val03-libxc/bin/rslmto.x`
Configuration: GNU Debug, serial CPU, `strux_backend='strux_lib'`, libXC 5.2.3

## Scope and result

This validation establishes the narrow collinear PBE-LDA scope:

- `txc=5`: the existing internal PBE-LDA implementation;
- `txc=105`: the existing explicit libXC mapping `[XC_LDA_X, XC_LDA_C_PW]`.

No second XC implementation was added. `self%VXC0SP` now calls the existing
`XCPOT_hybrid` dispatcher. Legacy `txc=1-9` values remain on the internal
route; explicit 100-series and `1000+ID` values select the libXC wrapper.

The reciprocal workflow compares the same lean Si (`nsp=1`) and bcc Fe
(`nsp=2`) fixtures at a 4³ Gamma-centered mesh, Gaussian broadening 0.01 Ry,
300 K, and a common 24-step SCF budget. A separate bcc-Fe Block real-space
check is included in the validation driver. The CI registration runs only the
four reciprocal cases; the six-run driver is deliberately kept out of the
quick gate because the real-space pair takes about 3.5 minutes.

## Evidence

The fixed-density unit compares the internal and libXC energy density and both
spin-channel potentials at four spin densities. The largest reference-path
difference is about (8\times10^{-7}) Ry in the potential and (5.3\times10^{-7})
Ry in the energy density; all direct equivalence assertions pass.

The SCF comparisons use report values and the full-precision `*_out.nml`
component records:

| route / fixture | EF | site charge | moment | |Δ band energy| | |Δ total energy| | max component difference | SCF behavior |
|---|---:|---:|---:|---:|---:|---:|---|
| reciprocal / Si | 0 at printed precision | 0 at printed precision | 0 | 2.36×10⁻⁶ Ry | 1.72×10⁻⁵ Ry | 8.84×10⁻⁶ Ry | 24 iterations in both paths; final residuals ≤ 6×10⁻⁸ |
| reciprocal / bcc Fe | 0 at printed precision | 0 | 0 at printed precision | 2.51×10⁻⁶ Ry | 1.93×10⁻⁵ Ry | 2.00×10⁻⁵ Ry | 16 iterations in both paths; final residuals ≤ 2×10⁻⁸ |
| real space / bcc Fe | 1×10⁻⁶ Ry | 0 | 1×10⁻⁶ μB at printed precision | 5.20×10⁻⁶ Ry | 1.88×10⁻⁵ Ry | 2.16×10⁻⁵ Ry | 24-iteration budget; final residuals 1.0×10⁻⁶ and 4.9×10⁻⁷ |

The compared energy components are identically defined production `par`
fields: `sumec`, `sumev`, `etot`, `utot`, `ekin`, and `rhoeps`. Equality of
total energy was not assumed; it was checked only after those component
conventions matched. `rhoeps` is retained as the reported radial XC energy
term. The fixed-density potential assertions provide the pointwise XC
potential contract used by the SCF dispatcher.

## Systematic differences

The residual differences are systematic and expected: the internal legacy
PBE-LDA path and libXC 5.2.3's `LDA_X + LDA_C_PW` path use corresponding
mathematical functionals but not byte-identical implementations or constants.
The small XC energy/potential difference propagates through the nonlinear SCF
iteration into the reported band, total, and radial energy components. The
same iteration counts, charges, moments, EF values at report precision, and
decreasing residual sequences show no unexplained offset or route failure.
Tolerances were not loosened to hide a systematic mismatch.

## Checklist

- [x] Equivalent internal/libXC functional identified.
- [x] Nonmagnetic SCF compared (Si/sp).
- [x] Magnetic SCF compared (bcc Fe/spd).
- [x] XC potential compared in the fixed-density production-path unit gate.
- [x] Charge compared.
- [x] Moment compared.
- [x] Energy conventions verified before comparison.
- [x] SCF convergence compared.
- [x] No duplicate XC implementation created.
- [x] libXC maturity scope updated.

## Commands and files

```text
ctest --test-dir /tmp/rslmto-val03-libxc -R '^UnitLibxcXcBaseline$' --output-on-failure
ctest --test-dir /tmp/rslmto-val03-libxc -R '^Val03LibxcScfEquivalence$' --output-on-failure
python3 tests/validation/val03_libxc_scf.py ...       # reciprocal + RS campaign
```

Changed files are the XC dispatch, the existing conditional fixed-density
unit, the one conditional workflow registration and driver, this report, and
the Phase-II maturity ledger.
