# Chebyshev backend matrix

The backend-only Chebyshev rows use the canonical nonmagnetic diamond-Si sp
fixture at `../scf/cases/bulk/diamondSi`: periodic 12x12x12 geometry,
`lmax=1`, `nsp=1`, `strux_lib`, and 200 moments. The regression runner accepts
the case's `output_file` so this shared fixture can be used without copying its
inputs into the regression tree.

The Si rows cover the production CPU dispatches `fast`, `batched`, `legacy`,
`mkl_batch`, and `mkl_sparse`. The two MKL rows remain gated by
`ENABLE_MKL_KERNELS`. References were generated on 2026-08-13 from clean
production runs: the first three with serial `build-gc-debug`, and the MKL rows
with `build-gc-mkl-test` configured with `ENABLE_MKL_KERNELS=ON` and Intel
oneMKL. The Si rows passed in 34.76 s in the debug build (MKL rows skipped) and
32.82 s in the MKL build with two CTest workers.

The magnetic/composition coverage remains on Fe deliberately:

- `bccFe_chebyshev_fast_hoh` and `bccFe_chebyshev_legacy_hoh` retain the HOH
  and legacy composition paths.
- `bccFe_chebyshev_fast_ccor_2c` retains derivative structure-constant/CCOR
  coverage.
- The SCF matrix retains the magnetic Fe Chebyshev nsp=2, nsp=3, and nsp=4
  cases, including their HOH variants. Fe remains the magnetic Chebyshev
  anchor; no Block or Lanczos case was moved to Si.

Before the change, the eight-row Fe regression Chebyshev matrix took 89.77 s
real time with two CTest workers in `build-gc-debug`; its plain backend rows
alone took 17.18 s (fast), 23.64 s (batched), and 22.05 s (legacy). Those rows
were replaced only after confirming that each Si manifest row selects the same
dispatch branch in `recursion_chebyshev.f90`.
