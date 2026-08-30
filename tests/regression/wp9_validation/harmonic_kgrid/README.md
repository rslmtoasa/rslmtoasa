# WP09 harmonic and k-grid fixture

This directory contains the small, repeatable q lists used by
`tests/validation/val19_gbt_harmonic_kgrid.py`:

- `q_points_cone.dat` is the fixed small-q cone-angle sweep (`Gamma` and
  `q=(0,0,0.05)`).
- `q_points_grid.dat` contains `Gamma`, the small-q probe, and `q=(0,0,0.5)`.
- `q_points_constrained_reference.dat` is the stable WP08 corrected-MFT
  reference (`Gamma` and `q=(0,0,+/-0.025)`).

The runner reuses the committed bcc-Fe frozen-potential deck from
`tests/regression/wp9_validation/gammaH_sweep/base_kspace` and patches only
the mode, cone angle, k mesh, and q-file. It emits one normalized
`wp09_sweep.csv` containing raw energies and per-q occupation/mesh metadata,
then invokes `tools/analyze_gbt_wp09.py` to produce `wp09_analysis.json`.

The q/2 mapping flag is exact for this cubic Cartesian fixture when
`q_i * nk_i / 2` is an integer. It is reported as a control, not used to
alter an energy or to declare convergence.

The corrected constrained-MFT route currently runs the original WP08
theta=90 degree reference convention. A small-angle corrected-MFT plateau is
not claimed unless the solver converges for that campaign; the runner keeps
the reference evidence separate from the bare-MFT harmonic sweep.
