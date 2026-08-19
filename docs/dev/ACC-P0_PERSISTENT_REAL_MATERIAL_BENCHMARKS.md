# ACC-P0 persistent real-material benchmark

This note records the methodology repair required by
`RS_LMTO_ASA_ACC_PERFORMANCE_RESCUE.md`. It is deliberately separate from the
legacy ACC-06 result: ACC-06 cold-process measurements remain useful as cold
wall evidence, but they cannot be interpreted as steady-state GPU throughput.

## Audit of ACC-06

The old path had all of the following properties:

- `tests/benchmarks/benchmark_harness.py` launched a fresh executable for every
  warm-up and repetition. Its wall timer surrounded that complete process.
- `tests/benchmarks/reciprocal_crossover.f90` constructed a hand-filled model
  in `setup_model`; the `Si_sp` and `bccFe_spd` labels were inferred from
  `sites/lmax`, not loaded from Si or Fe production inputs.
- The driver timed `solve_backend`, and `solve_backend` called
  `make_execution_backend` inside that interval. CUDA context, stream,
  cuSOLVER handle, and initial workspace setup therefore belonged to the
  measured cold solve.
- `source/cuda/reciprocal_cuda.cpp` kept cumulative H2D, solver, and D2H
  counters. The old Fortran adapter exposed snapshots but no interval reset or
  delta operation. ACC-06 consequently had no valid same-interval comparison
  between a one-call wall time and transfer/solver counters.
- The CUDA adapter still invokes the reference `cusolverDnZheevd` once per
  matrix. ACC-P0 does not change that algorithm; true batching is deferred to
  ACC-P1.

The old approximately 0.4-second floor is therefore a cold-process/context
measurement, not evidence about a persistent reciprocal solve.

## Implemented ACC-P0 path

`ReciprocalAccP0Benchmark` is an opt-in CMake target. It loads a production
input in the fixture directory, runs the normal lattice → structure constants
→ potential → H(R) → reciprocal Fourier path, creates one typed backend, calls
`prepare_operator` once, solves the first request, performs in-process warm-ups,
resets interval CUDA counters, and then repeats the same workload in one
process. Each measured repetition records H(k) assembly, backend solve wall
time, and total steady workload time. The outer Python harness records the
single process wall time as `cold_process_wall_s`.

The CUDA reset hook clears only H2D/solver/D2H interval accumulators and call
counts. It does not destroy the CUDA context, stream, cuSOLVER handle,
workspace, or lifetime request counters. Device free/total memory is queried
before the first solve; the driver rejects a CUDA case when its conservative
two-times dense tile estimate does not fit.

The machine-readable profile contains:

```text
fixture source workload L natom nmat nominal_mesh actual_unique_nk tile
vectors backend cold_process_wall_s cuda_context_backend_init_s first_solve_s
steady_solve_median_s metric_repetitions Hk_CPU_s H2D_s solver_s D2H_s total_steady_s
memory_estimate_mib memory_free_before_mib memory_total_mib
```

`tests/benchmarks/accp0_real_material.py` writes the complete CSV/JSON table,
keeps crossover and matched-density workload labels distinct, and uses
temporary campaign fixtures so production test directories are not modified.
The H2D/solver/D2H fields are the reset-to-end interval deltas normalized per
measured repetition; `metric_repetitions` records that normalization factor.
The driver intentionally solves the explicit full mesh rather than silently
using an irreducible integration workset; `actual_unique_nk` is therefore the
production reciprocal-coordinate mesh count after code-level reciprocal
folding/deduplication, with no unreported symmetry reduction.

## Corpus and physics oracle

The primitive fixtures are the canonical validated
`tests/scf/cases/bulk/diamondSi` and `tests/scf/cases/bulk/bccFe` inputs. The
Fe L=2,3,4,5 fixtures copy the converged primitive `Fe.nml` state to explicit
site labels and construct repeated bcc primitive vectors and coordinates in a
production `crystal_sym='file'` input. No large-supercell SCF is run.

Before any Fe supercell timing, the script solves primitive CPU LAPACK at the
reciprocal points compatible with the L×L×L translation and solves supercell
Gamma with the same production assembler. It compares state count, sorted
energy multiset, and degeneracy groups. A failed identity aborts the campaign;
GPU output is never used as the fixture oracle.

The campaign command is documented in `tests/benchmarks/README.md`. A quick
CPU-only smoke run is:

```bash
python3 tests/benchmarks/accp0_real_material.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda \
  --output-dir results/benchmarks/accp0-smoke \
  --quick --skip-cuda --warmups 1 --repetitions 2
```

The full CUDA campaign must be run on the intended host after preflight. Its
output is machine-local evidence, not a committed correctness reference. No
ACC-P0 conclusion should be drawn from nominal mesh sizes alone; use the
reported `actual_unique_nk` and the interval-consistent steady totals.

## CUDA-enabled validation evidence

On the available host (two NVIDIA RTX A4000 devices, 16,376 MiB each), the
three-repetition Gamma confirmation used two in-process warm-ups and one tile
per request. The bcc-Fe steady totals were:

| L | nmat | CPU steady | CUDA steady | CUDA/CPU |
|---:|---:|---:|---:|---:|
| 1 | 18 | 0.0356 ms | 0.6473 ms | 18.2x |
| 2 | 144 | 1.641 ms | 6.817 ms | 4.15x |
| 3 | 486 | 18.79 ms | 26.59 ms | 1.42x |
| 4 | 1152 | 158.2 ms | 101.8 ms | 0.64x |
| 5 | 2250 | 1.292 s | 482.6 ms | 0.37x |

Thus the current reference CUDA solver has a real-material crossover around
the L=4 regime on this host, despite losing for primitive through L=3. These
are performance observations for the unchanged one-matrix-at-a-time Zheevd
adapter, not a case for changing the default backend or claiming universal
GPU benefit. The vector-request smoke campaign also completed for Si and
bcc-Fe primitive/L=2 values-plus-eigenvectors workloads.

The canonical bcc-Fe input uses the legacy structure-constant backend. The
temporary supercell fixtures use that same backend so the CPU folding oracle
does not conflate structure-constant implementations with supercell physics.
