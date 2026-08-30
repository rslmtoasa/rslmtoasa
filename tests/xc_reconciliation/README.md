# XCR-01 fixed-density XC reconciliation

This is a production-path diagnostic for the legacy `XCPOT` implementation
and the libXC wrapper. It uses identical spin densities generated from
`(r_s, zeta)` and does not enter the SCF loop.

Build with libXC and the standalone tests enabled, then run:

```text
python3 tests/xc_reconciliation/run_xc_lda_reconciliation.py \
  --binary build-xcr/bin/UnitXcLdaReconciliation \
  --output-dir tests/xc_reconciliation/results
```

The output directory contains the comparison CSV, Stoner-response CSV,
exchange-validation CSV, energy-gradient audit CSV, a JSON summary, a run log,
and two diagnostic PNGs.
The exact fully polarized rows are recorded separately with the legacy
zero-density guard status; they are not included in regular-range maxima.
