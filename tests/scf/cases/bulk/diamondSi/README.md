# Diamond-Si sp Chebyshev fixture

This is a deterministic functional fixture, not a converged physical Si
benchmark. It is derived from `example/bulk/Si/SP/` and keeps its two-site
diamond primitive cell and production Si potential data.

The checked path is a periodic Bravais calculation with `n1=n2=n3=8`,
`lmax=1` (sp), `nsp=1`, `recur='chebyshev'`, `lld=200`,
`cheb_backend='fast'`, and `strux_backend='strux_lib'`. It uses one minimal
SCF step and checks finite total/site DOS and energy-related output values.

Measured provenance on 2026-08-13 using a clean fixture directory: serial
`build-gc-debug` run through `tests/run_binary.sh` completed in 5.25 s for the
selected 8x8x8, 200-moment case; the corresponding CTest run was 6.19 s
including harness setup. Clean 4x4x4 and 6x6x6 variants completed in 2.88 s
and 3.76 s, respectively, but 8x8x8 was retained for a more useful
deterministic DOS shape. The fixture is intentionally not marked `quick`.
