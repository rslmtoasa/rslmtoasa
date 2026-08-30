# WP08 corrected constrained frozen-magnon MFT

This fixture documents the smallest runnable corrected-MFT input.  It uses the
ordinary Fe database entry from `example/bulk/bccFe`, the `strux_lib` GBT
backend, and `soc_scale = 0`/`fix_soc = .true.` because the current GBT audit
rejects active SOC.  The production input is
`example/bulk/bccFe/input_frozen_magnon_constrained.nml`; copy it to
`input.nml` beside `Fe.nml` before running.

The corrected path has three independent layers:

1. The ordinary reference potential is converged once with constraints disabled.
2. For every q, only the constraining field is iterated against the local GBT
   collinear axis.  The ordinary potential is snapshotted and restored around
   each Hamiltonian/DOS pass.
3. The converged field is frozen for the final one-shot occupied-eigenvalue sum.

`frozen_magnon.dat` retains the legacy row shape.  The corrected audit is in
`frozen_magnon_constrained_diagnostics.dat`; its rows include the physical
reference total energy, raw field-included band sum, controller residual,
controller penalty and field-coupling diagnostics, ordinary-potential checksum,
post-restore drift, per-sublattice field magnitudes, and the final normalized
Delta E.  The two controller energy values are explicitly diagnostics and are
not added to the physical DFT total energy.

Set `constraint_start_from_zero = .true.` for q-order-independent restarts.  A
false value deliberately enables field/controller continuation and is an
efficiency option rather than the reproducibility default.
