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

## Metallic Block and Lanczos coverage

Block and Lanczos remain metallic real-space recursion coverage. The canonical
Fe examples use 20 recursion levels, which is within the lean 15–20-level
target; their existing radius-cut cluster is retained because it is the
established faster production contract and does not map naturally to the Si
fixture's periodic replication controls. No Si Block or Lanczos case exists.
The DOS samples in these tests are coarse finite-depth recursion diagnostics,
not sharply resolved Kohn–Sham DOS benchmarks.

The Fe nsp=2/3/4 Block cases retain both HOH and non-HOH coverage. The direct
`Lanczos`, `Block`, and `Chebyshev` CTest baselines are intentionally retained:
they use separate Fe nsp=1 decks and therefore are not duplicates of the
manifest's nsp=2/3/4 cases. They check the legacy nsp=1 scalar contract
(`etot`, `ws_r`, `vmad`) and took 45.44 s real time together with two CTest
workers before TEST-08 (Block 7.68 s, Chebyshev 23.67 s, Lanczos 45.44 s).
The comparable 12-case manifest slice took 73.21 s before TEST-08; after
removing only the redundant Pt2MnGa Block+HOH point, the retained 11-case
slice took 54.98 s.

The nsp=2 Lanczos functional cases are labeled `known_issue` and are not
`quick`. They remain runnable scalar/moment reproducers only: the production
path currently emits NaN spectral DOS/`lmom`, so those values are deliberately
not referenced or described as validated. The root-cause repair is a separate
production task; see `tests/KNOWN_ISSUES.md`.

## Historical WP8/WP9 validation material

The WP8 little-group and WP9 validation decks/runners are retained here as
historical scientific evidence. They are intentionally outside the active
CTest harness: CMake does not register them, so they do not appear in
`ctest -N` or any ordinary label selector. Running one requires an explicit
validation decision and invoking its retained runner directly.

Pt2MnGa retains its multi-sublattice metallic Block path with serial/MPI
coverage and its Chebyshev paths. The redundant Pt2MnGa Block+HOH Cartesian
point was removed because Fe already covers Block+HOH and the retained Pt2
Block case covers the distinct multi-species path in both launch modes.
