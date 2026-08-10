# TDDFT-11 CPU performance hardening

The measurements in this document are kernel-cost evidence only. They do not
lift the open WR-07/WR-08 mode-classification and material-physics gates in
`TDDFT_WARD_REPAIR_TRACKER.md`; GPU work remains deferred.

## Completed checklist

- [x] CPU timing is emitted with every `chi0`, Dyson, and mode output.
- [x] Profile fields separate arbitrary `k+q` eigensolutions, transition vertices, frequency denominators, response accumulation, Dyson solves/diagonalizations, mode analysis, and Green-function energy integration.
- [x] Folded duplicate `k+q` solutions remain reused by the arbitrary-k service.
- [x] χ_KS uses bounded transition tiles and BLAS `zgemm(N,T)` accumulation.
- [x] Frequencies are processed as a supplied batch; no response tensor over transition pairs is stored.
- [x] Exact default occupation handling is retained; occupation pruning remains opt-in and explicitly recorded.
- [x] MPI work partitioning remains over q points.
- [x] Scalar fixed-order and batched-GEMM χ_KS paths have a strict numerical-equivalence test.
- [x] TDDFT-03 through TDDFT-10 are covered by one CI equivalence gate.
- [x] GPU work is deferred: CUDA/HIP hardware is unavailable on this machine.

## CPU profile snapshot

`UnitTddftCpuProfile` was run in the local Debug build (GNU Fortran 16.1,
Apple Accelerate BLAS) on 2026-08-09.  It uses deterministic one-atom,
18-spinor `(s,p,d)` response fixtures matching the basis shape of the bcc Fe
and fcc Ni examples.  They are kernel-cost fixtures, not replacement material
validation calculations.

Times are process CPU seconds.  The columns are `vertices`, `denominators`,
`GEMM accumulation`, `k+q eigensolve`, `Dyson solve`, `Dyson diagonalization`,
`mode analysis`, and `GF integration`.

| Fixture | nk | nω | CPU profile (s) |
| --- | ---: | ---: | --- |
| bccFe | 16 | 96 | 0.005823, 0.011847, 0.011021, 0, 0.000175, 0, 0.000297, 0.081919 |
| fccNi | 32 | 192 | 0.011800, 0.047096, 0.042887, 0, 0.000233, 0, 0.000436, 0.081910 |

The fixture feeds precomputed eigenpairs, hence its `k+q eigensolve` column is
zero.  Production q runs populate it through
`profile_arbitrary_kq_eigensolve_cpu_s` in every χ₀ file.  The same output
headers expose all other fields, so an Fe/Ni production run needs no special
profiler build or changed physics input.

The result identifies Green-function integration, followed by denominator and
response accumulation, as the hotspots in these small-basis fixtures.  The
eigenpair route therefore uses GEMM tiles; no GPU path was added without
hardware evidence.

## Numerical and parallel policy

The optimized χ₀ result is compared with the retained scalar reduction path
on a complex, multi-channel fixture with a relative limit of
`256*epsilon(rp)`.  This permits only floating-point reduction-order effects;
the response formula, retarded denominator, and non-Hermitian right vertex are
unchanged.  `TddftCrossMilestoneEquivalence` runs the TDDFT-03…10 oracle suite
in CI.

The production scheduler uses MPI-over-q.  Each q rank calls BLAS for its
transition tiles; do not additionally OpenMP-parallelize that GEMM loop.
When the selected BLAS is threaded, set `OMP_NUM_THREADS=1` for the response
driver (or otherwise choose a single thread owner) to avoid nested-thread
oversubscription.  The scalar reference path keeps deterministic transition
ordering; batched BLAS reductions can differ only within the documented test
tolerance.

GPU acceleration, CPU/GPU equivalence testing, and GPU reduction-tolerance
reporting are intentionally deferred until CUDA or HIP hardware is available.
