# TDDFT-14 performance and MPI report

This report records the measured optimization evidence for the validated
eigenpair, K-space Green-function, and native R-space Green-function paths.
The machine-readable campaign is [tddft14_cpu.json](../../results/benchmarks/tddft14_cpu.json).

## Campaign and backend coverage

`UnitTddftCpuProfile` exercises deterministic Fe/Ni-shaped fixtures with
`(nk,nω)=(16,96)` and `(32,192)`. It runs eigenpair transitions in batched and
scalar-reference modes, K-GF direct and mixed-contour integration, and native
R-GF direct integration with a three-q batch. The sampled real-axis native
source rejects mixed contour explicitly; this unsupported combination is
reported and guarded rather than approximated.

| fixture | eigenpair batched/scalar | K-GF direct/mixed | native R-GF | principal payload | native source + q batch |
|---|---:|---:|---:|---:|---:|
| bccFe one-site | 0.033961 / 0.520390 s | 0.028060 / 0.226710 s | 0.131230 + 0.000008 s | 1.6441 MiB | 1.3066 MiB |
| fccNi two-site | 0.496470 / 8.033100 s | 0.169200 / 1.473700 s | 0.135960 + 0.000009 s | 4.5315 MiB | 1.3066 MiB |

The native R-GF timing is split into energy/pair integration and the batched
R-to-q Fourier transform. K-GF's existing production timer covers its single
GF-contraction/energy-integration loop; it is reported as that combined stage.
Dyson and spectral timings are also emitted separately by the fixture: Fe
`0.000109 / 0.000424 s`, Ni `0.000126 / 0.000867 s`.

## Optimizations and scientific checks

The transition engine now samples `cpu_time` outside the `(n,m)` transition
loop. Native R-GF reuses interpolation, advanced-GF, and matrix-multiply
scratch across the entire pair/frequency build; the inner energy loop performs
no allocatable-matrix creation. The native path still performs one q-batched
Fourier transform after R-space response construction.

The scalar/batched eigenpair checksums are identical for Fe and differ by only
`2.75e-21` for Ni. The direct K-GF and mixed-contour checksums, plus the native
R-GF checksum, are retained in the JSON evidence. Existing backend-equivalence
and real-space Fourier tests remain the numerical gates.

Measured scalar-over-batched eigenpair speedups are `15.32x` (Fe) and `16.18x`
(Ni). The mixed-contour/direct K-GF ratios are `8.08x` and `8.71x`; these are
integration-cost observations, not claims that mixed contour is an
optimization.

## MPI decomposition

The planner selects q-outer ownership for reciprocal backends, with omega or k
fallbacks when q does not fill the worker set. For native R-GF it partitions
R blocks first, then energy nodes, reduces the partial `chi0(R,ω)` tensor, and
keeps the complete q batch on every rank for the final Fourier transform. Only
the normal q owner range writes q-labelled output. The four-rank MPI CTest
passed with 12 R blocks and 3 output q points per rank; the intentional q
duplication factor is 4 and is bounded to preserve R-space reuse.

The MPI unit contract also checks oversubscribed work ranges, actual-rank q
ownership, collective-reduction flags, and the shared-q invariant. Serial
backend equivalence, native R-GF, workspace, and performance tests all pass.

## GPU decision

GPU acceleration is explicitly deferred. No CUDA/HIP crossover measurement is
available on this host, and the measured fixture contains small matrix kernels
where transfer/setup overhead is material. The validated CPU BLAS path remains
the selected implementation until a sufficiently large, batched GPU campaign
can provide crossover and checksum evidence.

## Reproduction

```bash
python3 tests/benchmarks/tddft14_performance.py \
  --binary build/bin/UnitTddftCpuProfile \
  --output results/benchmarks/tddft14_cpu.json \
  --warmups 1 --repetitions 3
```

For the MPI contract, configure with `-DENABLE_MPI=ON` and run the registered
`UnitTddftPerformance_mpi` test with four ranks.
