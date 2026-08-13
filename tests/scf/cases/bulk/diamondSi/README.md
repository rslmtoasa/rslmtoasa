# Diamond-Si sp Chebyshev fixture

This is a deterministic functional fixture, not a converged physical Si
benchmark. It is derived from `example/bulk/Si/SP/` and keeps its two-site
diamond primitive cell and production Si potential data.

The checked path is a periodic Bravais calculation with `n1=n2=n3=12`,
`lmax=1` (sp), `nsp=1`, `recur='chebyshev'`, `lld=200`,
`cheb_backend='fast'`, and `strux_backend='strux_lib'`. It uses one minimal
SCF step and checks finite total/site DOS and energy-related output values.

Measured provenance on 2026-08-13 using clean fixture directories: the
12x12x12, 200-moment fast case took 13.16 s through CTest in serial
`build-gc-debug`; its MKL-batch variant took 15.94 s in the
`ENABLE_MKL_KERNELS=ON` `build-gc-mkl-test`. The fixture is retained because it
gives the stable DOS shape needed by backend comparisons while remaining lean,
and is intentionally not marked `quick`.

The fast reference was generated with `build-gc-debug`. The MKL-batch reference
was generated and checked with `build-gc-mkl-test`; both were clean production
runs with the explicit `strux_lib` path. The stored checks are deterministic
functional outputs (energy-related scalars, Wigner-Seitz radius, Madelung
potential, total DOS, Si-projected DOS, and Fermi level), not a converged Si
benchmark.
