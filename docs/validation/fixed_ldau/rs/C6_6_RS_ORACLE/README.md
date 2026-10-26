# C6.6 versioned RS no-U reference

This is a refreshed, current-runtime reference snapshot for the single-thread no-U RS identity
check. The original C5.2 reference remains unchanged and is retained as historical evidence.

The new reference was produced with the retained RD1 executable (SHA-256
`4bb76cc2bbd251d95cdc434eec0348b9082069af2deeb4f4859f6d8056779a4f`) in a fresh directory
containing only the archived `input.nml` and `Fe.nml`. Its startup record confirms the solver
loaded `./Fe.nml`; no `Fe_out.nml` existed before launch. The checked preflight passed, including
the host-E_F comparison at `-0.069587 Ry`. The current no-LibXC executable was run from the same
two initial files through the same wrapper. Its 27 output artifacts match this reference exactly
(normalizing only the runtime `cpu_s` field in `str.out`).

Against the historical C5.2 output, the current-runtime output changes 12 artifacts. The initial
deck and restart hashes still match the historical inputs. The solver rewrites `input.nml`'s
Fermi value from `-0.069587` to `-0.069588 Ry`; the logged free Fermi energy is also `-0.069588`
instead of the historical `-0.069587 Ry`. The d-occupation differences are below `1e-7`, the
reported sphere moment differs by `2.42e-5 μB`, and the orbital moment differs by `9.38e-6 μB`.
The DOS tables also show grid shifts up to `6.4e-4 Ry` and pointwise DOS changes up to about
`0.57`, so the full artifact drift is not adequately described as a small scalar delta alone.

The paired rerun supports refreshing the active reference while preserving the historical
snapshot. It does not establish that the Fermi rewrite alone caused every changed spectral
artifact. `ORACLE.json` records provenance and hashes; `SHA256SUMS` covers the 27 reference
artifacts. The legacy C5.7 identity record is left unchanged.

For a future checked identity run, pass this snapshot's `artifacts/` directory through
`tests/validation/fixed_ldau/c5_7_no_u_identity.py --rs-baseline` and write the result with
`--result` to a versioned C6 record. The default runner behavior continues to target the historical
C5.2 baseline.
