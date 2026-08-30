# XCR-04 magnetic SCF campaign

This directory contains the reproducible ordinary `nsp=2`, `q=0` collinear
campaign and its compact scalar results.  It uses fresh-start atomic data,
four outer iterations of a temporary `B_fsm` seed for nonzero labels, then
continues without the seed.  No constraint namelist is enabled.

Run it from the repository root with:

```bash
python3 tests/magnetic_scf/run_magnetic_scf_regression.py \
  --binary build-xcr/bin/rslmto.x \
  --output-dir tests/magnetic_scf/results \
  --max-steps 6 --keep-traces
```

The stored sweep is a bounded diagnostic campaign.  All rows retain their
convergence status; `BOUND_NOT_CONVERGED` is not a claim of a fixed point.
