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

`UnitTddftCpuProfile` uses deterministic completed reciprocal fixtures, not
an SCF material calculation: a one-site 18-spinor `(s,p,d)` model with
`(nk,nω)=(16,96)` and a two-site 36-spinor model with `(32,192)`.  They are
representative one- and multi-site workload shapes, not replacement material
validation calculations.

Times are informational process CPU seconds. `PROFILE_RECIPROCAL` labels
Fourier assembly, normal-k eigensolution, arbitrary-`k+q`
assembly/eigensolution, and finite-q LMTO pair-potential construction.
`PROFILE_TDDFT` labels vertex construction, denominator generation, response
accumulation, Green-function energy integration, Dyson, and mode analysis.
`PROFILE_MEMORY_MIB` reports analytical principal-array payloads: H(k), both
eigenpair sets, the pair-operator cache plus Q+ workspace, response tensors,
and their sum. It assumes FP64 reals (8 bytes) and complex FP64 values (16
bytes); it is deliberately not an OS-level RSS measurement.

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
