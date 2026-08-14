# Diamond-Si sp Chebyshev fixture

This is a deterministic functional fixture, not a converged physical Si
benchmark. It is derived from `example/bulk/Si/SP/` and keeps its two-site
diamond primitive cell and production Si potential data.

The checked path is a periodic Bravais calculation with `n1=n2=n3=12`,
`lmax=1` (sp), `nsp=1`, `recur='chebyshev'`, `lld=200`,
`cheb_backend='legacy'`, and `strux_backend='strux_lib'`. It uses one minimal
SCF step and checks finite total/site DOS and energy-related output values.
The functional fixture deliberately uses the stable legacy Chebyshev route;
the FP32 `fast` route is covered separately by the backend regression matrix.

Measured provenance on 2026-08-13 using clean fixture directories: the
12x12x12, 200-moment fast case took 13.16 s through CTest in serial
`build-gc-debug`; its MKL-batch variant took 15.94 s in the
`ENABLE_MKL_KERNELS=ON` `build-gc-mkl-test`. The fixture is retained because it
gives the stable DOS shape needed by backend comparisons while remaining lean,
and is intentionally not marked `quick`.

The canonical functional reference is generated with the CI-equivalent runner.
The MKL-batch reference is generated and checked with an
`ENABLE_MKL_KERNELS=ON` build; both are clean production runs with the explicit
`strux_lib` path. The stored checks are deterministic
functional outputs (energy-related scalars, Wigner-Seitz radius, Madelung
potential, total DOS, Si-projected DOS, and Fermi level), not a converged Si
benchmark.

## Chebyshev/k-space DOS equivalence

`Example_si_chebyshev_kspace_dos_equivalence` runs this same fixture through
the real-space Chebyshev route and the production k-space eigensystem/DOS
route. Both routes use `strux_lib`, `lmax=1` (sp), `nsp=1`, the same
`-1.5 ... 1.0` Ry energy range, and 1000 DOS channels. The k-space route uses
an unreduced 8x8x8 mesh and Gaussian broadening with fixed sigma `0.02` Ry.
The real-space route uses the fixed Chebyshev order `lld=200`.
Here `nsp=1` retains the project's scalar-relativistic collinear spin
treatment; it is not a spin-disabled calculation.

These are deliberately lean functional settings, not converged material
physics. The 0.02 Ry broadening was frozen after a one-time comparison study:
with the rebuilt serial binary, the selected pair integrated to 16.000005 and
16.000002 states, had central-window relative DOS RMS 0.141 and maximum
absolute difference 11.60, and its principal DOS peaks differed by 0.005 Ry.
The CI driver uses fixed tolerances (relative RMS 0.20, maximum absolute
difference 12.5, and peak separation 0.08 Ry); it performs no runtime fitting
or DOS calculation. It only runs the two production routes, parses
`totaldos.out`, interpolates the shifted grids to a common absolute-energy
window `-1.2 ... 0.8` Ry, integrates, and reports differences. Route runtimes
are printed by the driver for each run; the final recorded serial run took
12.27 s for real-space Chebyshev, 2.15 s for k-space Gaussian, and 14.46 s
total (CTest wall time 14.52 s).
