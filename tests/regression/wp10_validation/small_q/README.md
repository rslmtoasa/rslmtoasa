# WP10 symmetric small-q fixture

`q_points_symmetric.dat` starts at Gamma and then contains matched positive
and negative q points along the cubic z direction.  The coordinate is the
Cartesian dimensionless value used by `hamiltonian%q_ss`, namely `2*pi/alat`.
The largest points are included to expose the `q^4` term; the smallest points
provide the long-wavelength window.

The fixture is consumed by `tests/validation/val20_gbt_small_q.py`, which
retains the raw WP09 rows in `wp10_sweep.csv` and writes the fit and
q-reversal analysis to `wp10_analysis.json`.
