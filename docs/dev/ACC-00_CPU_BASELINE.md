# ACC-00 CPU baseline report

This is the initial Phase III-A benchmark-harness evidence recorded on 2026-08-18.
The measurements are baseline evidence, not correctness references or timing
thresholds.

## Environment

| Field | Value |
|---|---|
| git commit | `b19ed5e9222e351b45309871204981867520317d` |
| compiler | GNU Fortran 13.3.0 (`/usr/bin/gfortran`) |
| build | Release; OpenMP enabled; MPI disabled |
| BLAS/LAPACK | Intel oneMKL 2026.1 |
| OpenMP threads | 2 |
| MPI ranks | 1 |
| CUDA/cuSOLVER | unavailable in this CPU campaign |
| GPU | unavailable |
| CPU model | reported by this host as `x86_64` |
| warm-ups / repetitions | 1 / 3 |

## End-to-end CPU baselines

Times are wall-clock seconds for the existing production runner, including
fixture setup and output checks. RS rows have no reciprocal k-point or tile;
those fields are recorded as not applicable by the harness.

| Benchmark | Fixture / route | Matrix dimension | k points | Median | Minimum | Spread |
|---|---|---:|---:|---:|---:|---:|
| reciprocal_si_sp | diamond Si / sp | 16 | 8 | 0.850508 | 0.850035 | 0.002333 |
| reciprocal_fe_spd | bcc Fe / spd magnetic | 18 | 8 | 2.606853 | 2.594572 | 0.021765 |
| rs_block_fe | bcc Fe / Block | 18 | n/a | 3.249646 | 3.232277 | 0.017640 |
| rs_lanczos_fe | bcc Fe / Lanczos | 18 | n/a | 15.141974 | 15.035657 | 0.144788 |
| rs_chebyshev_si | diamond Si / Chebyshev moments | 16 | n/a | 11.064130 | 10.998646 | 0.912704 |

## Reciprocal component profile

The existing `UnitTddftCpuProfile` supplies deterministic component records
for a one-site bcc-Fe model and a two-site fcc-Ni model. These are component
microbenchmarks; the production fixture rows above are the end-to-end physical
baseline.

| Profile | Fourier assembly median | Normal eigensolution median | Arbitrary-k eigensolution median |
|---|---:|---:|---:|
| bccFe_one_site | 0.002324 | 0.000009 | 0.000189 |
| fccNi_two_site | 0.004393 | 0.000005 | 0.003856 |

The raw machine-readable records were generated with
`tests/benchmarks/benchmark_harness.py`; they are intentionally kept outside
the source tree because machine-local timings are evidence, not correctness
fixtures. Re-run the campaign on the target machine before making crossover,
GPU, or batching decisions.

The production benchmark commands run smoke/output checks. Reference
comparison remains an explicit option because the current dirty validation
state does not match every committed reference; no reference was updated for
this task.
