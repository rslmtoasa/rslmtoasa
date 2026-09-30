# Validation: forced nonmagnetic atoms

Base: `rslmtoasa/rslmtoasa`, `main`, commit
`49c0031b6406a93cabc768a9d698796c1ed90fd9`.
Test branch: `feat/force-nonmagnetic`. No PR was opened.

## Build and regression results

- GNU Fortran 13.3.0, Release build; OpenBLAS 0.3.26; OpenMP enabled.
- Serial build: successful, including the projection-test executable.
- Open MPI 4.1.6 build (`ENABLE_MPI=ON`): successful.
- Final serial CTest suite: 3/3 passed (projection, nsp=1 SCF, nsp=4 SCF).
- Final two-rank MPI CTest: 1/1 passed.
- Both SCF tests run four iterations with Fe magnetic and Co constrained.
  They use linear mixing (beta=0.05), radius=12, recursion depth=10,
  800 energy channels, and `hoh=.true.`. The SOC input axis is (0.6,0,0.8).
  These are numerical regressions, not fully converged material predictions.
- Projection test checks every QL/PL spin sum and radial-density spin sum
  to absolute tolerance 1e-14, all four potential averages, `build_pot`
  differences, field removal, and the untouched magnetic atom.
- Omitted flags and explicit all-false flags produced exactly identical
  parsed potential files for both atoms, in serial and MPI runs.
- Separately compiled unmodified commit: default potential data agree within
  1e-7 absolute/relative tolerance. Largest absolute entry differences were
  4.05083611099144e-9 (collinear) and 1.1641532182693481e-9 (SOC).
- Zero-moment direction/angle outputs contain no NaNs or infinities.
- `git diff --check`: passed.

## Numerical results

Q values are valence occupations. Potential columns are maxima of absolute
values over all orbital components, evaluated by `build_pot` after the atomic
update. Moments are in Bohr magnetons. All listed zeros were exactly zero
in the runtime diagnostics, rather than values merely rounded to zero.

| Test | Co Q_up | Co Q_down | Q_up-Q_down | max abs cx1 | max abs wx1 | max abs cex1 | max abs obx1 | Co constrained moment | Fe moment |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Collinear, serial (nsp=1) | 4.590594089271 | 4.590594089271 | 0 | 0 | 0 | 0 | 0 | 0 | 2.839167 |
| Tilted SOC, serial (nsp=4) | 4.590252048665 | 4.590252048665 | 0 | 0 | 0 | 0 | 0 | 0 | 2.835225 |
| Tilted SOC, MPI (2 ranks) | 4.590252048543 | 4.590252048543 | 0 | 0 | 0 | 0 | 0 | 0 | 2.835225 |

The constrained moment is the projected SCF moment. Raw energy-resolved band
spin response is retained and can be finite through magnetic hybridization;
this implementation is not a zero-raw-moment Lagrange-multiplier solver.
Projection conserves charge at each application; charge remains free to
redistribute through the SCF cycle.

## Exact files changed

- `CMakeLists.txt`
- `source/include_codes/namelists/self.f90`
- `source/self.f90`
- `source/bands.f90`
- `source/mix.f90`
- `source/symbolic_atom.f90`
- `source/potential.f90`
- `tests/nonmagnetic/projection.f90` (new)
- `tests/nonmagnetic/run_scf.py` (new)
- `tests/nonmagnetic/README.md` (new)
- `tests/nonmagnetic/VALIDATION.md` (new; this report)

Hamiltonian routines were inspected. Their existing spin-sum/difference
construction already provides the required behavior once the four LMTO
parameter pairs are equal; no Hamiltonian source edit was required.

Portable build/test commands and the SCF application points are in README.md.
