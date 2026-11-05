# RA1 result: impurity energy diagnostic identity

**Measured: pass.** The final B2FeCo `calctype='I'`, `nsp=1`, one-step run
reported impurity symbolic type 3 and its `E_int`, `E_DC`, and `E_Hub` agreed
with a separately invoked Fortran evaluator on the serialized d-occupation
matrix from `Fe_out.nml`. Maximum absolute difference was `8.891e-13 Ry`, below
the preregistered `1e-12 Ry` per-component threshold. The evaluator uses the
production `evaluate_hubbard_energy` functional; the matrix is independently
parsed from the emitted potential file, and the test explicitly distinguishes
symbolic type `nbulk+ia=3` from geometric representative `ib(1)=16`.

The exploratory mapping failure and unchanged tolerance are documented in
[`RA1_MAPPING_AMENDMENT.md`](RA1_MAPPING_AMENDMENT.md). Final run products are
under `/tmp/rslmto-ra1-verified/` and hashes are in
[`RA1_PROVENANCE.json`](RA1_PROVENANCE.json).
