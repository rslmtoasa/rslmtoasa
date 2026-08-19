# RS-LMTO-ASA Phase III-A Accelerator Blueprint and Luna Work Packages

**Target branch:** `fable_v3`  
**Phase:** III-A — Performance and accelerator establishment  
**Status entering this phase:** TEST, STAB, and VAL campaigns completed; GBT and TD-DFT remain separate Development tracks and are not on the accelerator critical path.

---

# 1. Purpose

Phase III-A introduces and establishes accelerator support for the parts of RS-LMTO-ASA that now have trustworthy CPU scientific reference paths.

The objective is **not** to "move RS-LMTO-ASA to the GPU".

The objective is to accelerate demonstrable bottlenecks while preserving:

- the real-space LMTO / Green-function identity of the code;
- the existing physical domain objects and public APIs;
- validated CPU production routes as correctness references;
- the established parent-module/submodule architecture;
- developer readability and maintainability;
- single-GPU usability as the primary design point;
- MPI rank-to-device extensibility without premature distributed-GPU complexity.

The first reciprocal accelerator target is deliberately narrow:

```text
short-range H(R)
      |
      | CPU Fourier assembly
      v
H(k) tile on host
      |
      | H2D
      v
GPU dense Hermitian eigensolver
      |
      | D2H
      v
eigenvalues/eigenvectors
      |
      v
existing validated CPU consumers
```

Only after that route is correct and profiled should the project decide whether device-resident eigenstates, GPU Lehmann contractions, GPU H(k) assembly, asynchronous overlap, or other accelerator extensions are justified.

---

# 2. Scientific and architectural constraints

## 2.1 Preserve the RS-LMTO-ASA center of gravity

RS-LMTO-ASA remains fundamentally a real-space LMTO / Green-function code.

The accelerator program must not reorganize the code around reciprocal space.

Existing real-space CUDA functionality remains a separate numerical implementation suited to the RS algorithms.

The reciprocal accelerator implementation should use the existing reciprocal execution-backend seam rather than creating a new global accelerator abstraction.

## 2.2 Protect the atomic LMTO core

The mature atomic/radial/potential routines in `self.f90` and related files are protected scientific infrastructure.

They should remain CPU-side and should not be modified for GPU convenience, memory-layout uniformity, stylistic modernization, kernel fusion, or device-residency schemes.

The accelerator boundary should start downstream of the established atomic LMTO physics.

## 2.3 Accelerators are execution technologies

The physical architecture remains:

```text
physical LMTO state
      |
      +---- real-space algorithms ---- RS CPU/CUDA execution
      |
      +---- reciprocal algorithms ---- CPU/LAPACK or CUDA execution
```

Do not introduce generic "GPU data providers", universal device arrays, or application-wide device ownership.

A new abstraction is acceptable only when at least two real production consumers require it and the performance evidence demonstrates its value.

---

# 3. Phase III-A non-goals

The following are explicitly outside the initial accelerator critical path:

- GBT-specific GPU kernels;
- TD-DFT-specific GPU kernels;
- GPU atomic/radial LMTO;
- GPU electrostatics;
- GPU surface/interface orchestration;
- distributed dense eigensolvers;
- NCCL;
- CUDA-aware MPI;
- mandatory multi-GPU execution;
- simultaneous CUDA and HIP implementations;
- speculative device-resident state frameworks;
- automatic GPU H(k) assembly before profiling;
- changing validated CPU numerical algorithms merely to resemble GPU code.

GBT and TD-DFT may later benefit automatically from shared reciprocal infrastructure after their scientific Development issues are resolved.

---

# 4. Correctness reference policy

The CPU implementation is the accelerator oracle only where Phase II established the CPU scope.

Accelerator validation must compare against the **same production physical route**, not a new Python physics implementation.

## 4.1 Reciprocal eigensolver contracts

Do not compare complex eigenvectors element-by-element.

For standard Hermitian problems test:

\[
|\epsilon_n^{GPU}-\epsilon_n^{CPU}|
\]

\[
r_n =
\frac{\|H v_n-\epsilon_n v_n\|}
     {\max(\|H\|,1)}
\]

\[
\|V^\dagger V-I\|
\]

For degenerate or nearly degenerate groups compare subspace projectors

\[
P = VV^\dagger
\]

rather than individual vector phases/orientations.

For generalized Hermitian problems, once supported:

\[
V^\dagger S V \simeq I
\]

and the generalized residual must be checked.

## 4.2 Downstream physical contracts

Use established functional fixtures to compare:

- electron count;
- Fermi energy;
- band energy;
- site charge;
- magnetic moment;
- selected eigenvalues;
- DOS/state count;
- selected Green-function entries;
- selected conductivity/exchange/damping quantities where the accelerated route feeds them.

## 4.3 Reference changes

Never update CPU scientific references merely because a GPU implementation differs.

First determine whether the discrepancy is due to legitimate floating-point ordering, degenerate-subspace freedom, tolerance choice, CPU bug, GPU bug, or unsupported algorithmic difference.

---

# 5. Benchmark-first policy

Performance decisions must be based on measured end-to-end cost.

A kernel is not considered successfully accelerated merely because the kernel itself is faster.

For reciprocal eigensolution measure:

\[
T_{\mathrm{total}} =
T_{H(k)} +
T_{\mathrm{H2D}} +
T_{\mathrm{eig}} +
T_{\mathrm{D2H}} +
T_{\mathrm{consumer}}
\]

as applicable.

For existing RS GPU paths measure comparable end-to-end production observables and timings.

Performance tests are not ordinary correctness gates. They should be reproducible, machine-described, and stored as benchmark evidence rather than flaky pass/fail timing assertions.

---

# 6. Benchmark fixture vocabulary

Use real LMTO Hamiltonians wherever possible.

| Fixture | Main performance role |
|---|---|
| diamond Si / sp | very small reciprocal matrix |
| bcc Fe / spd | small magnetic reciprocal matrix |
| B2 FeCo or another multisublattice metal | intermediate matrix / multi-site |
| larger existing multi-site primitive or controlled supercell | larger reciprocal matrix |
| existing metallic RS fixture | Block/Lanczos performance |
| existing Si/sp fixture | Chebyshev performance |
| existing Pt/SOC transport fixture | KPM/Kubo-Bastin performance |

Synthetic Hermitian matrices may be used for isolated eigensolver microbenchmarks, but every performance conclusion must also be checked on real assembled LMTO matrices.

---

# 7. Benchmark metadata

Every recorded benchmark should include enough environment metadata to interpret it later:

- git commit;
- compiler and version;
- build type;
- BLAS/LAPACK implementation;
- OpenMP thread count;
- MPI rank count;
- CUDA toolkit version;
- cuSOLVER version/API path if known;
- GPU model;
- GPU driver/runtime;
- CPU model;
- matrix dimension;
- number of k points;
- tile size;
- eigenvalues-only vs eigenvectors;
- standard vs generalized problem;
- transfer inclusion/exclusion;
- warm-up policy;
- number of repetitions.

Do not compare performance numbers from materially different environments as if they were controlled speedups.

---

# 8. Reciprocal backend design

The existing reciprocal execution backend is the correct seam.

The CUDA implementation should remain behind the same physical caller interface.

## 8.1 CUDA v1 supported scope

Initial CUDA reciprocal support should target:

- standard Hermitian problem;
- ordinary Hamiltonian-only representation;
- eigenvalues;
- eigenvectors;
- arbitrary-k API;
- ordinary reciprocal mesh;
- Si/sp nonmagnetic fixture;
- bcc-Fe/spd collinear magnetic fixture;
- single GPU;
- synchronous execution;
- CPU H(k) assembly.

## 8.2 Explicitly deferred from CUDA v1

Defer:

- generalized overlap;
- higher-order overlap paths not already scientifically established;
- noncollinear/SOC extensions until ordinary support is correct;
- GBT;
- TD-DFT;
- GPU H(k) assembly;
- device-resident spectral handoff;
- asynchronous multi-buffer execution;
- multi-GPU collectives.

---

# 9. Backend lifecycle and operator residency

The execution contract already contains operator generations.

A CUDA backend should use them to avoid repeated device uploads:

```text
if requested_operator_generation != prepared_generation:
    upload/rebuild persistent operator data
    prepared_generation = requested_operator_generation
else:
    reuse device state
```

Backend recreation or ordinary caller reuse must not invalidate a resident operator unless the physical operator actually changed.

Device state must have explicit ownership and release semantics. Avoid hidden global CUDA state.

---

# 10. H(R) -> H(k) policy

Initial policy:

> Keep reciprocal Fourier assembly on CPU.

After GPU eigensolution is working, profile H(k) assembly separately.

Only port it if it becomes a meaningful end-to-end bottleneck.

A Phase III-A work package is allowed to conclude:

> GPU H(k) assembly is not justified for the measured workload.

That is a successful performance result, not a failed task.

---

# 11. Eigenstate residency policy

CUDA v1 should return ordinary host eigenpairs to existing consumers.

Do not introduce a device-resident spectral handle before a real second GPU consumer exists.

If GPU Lehmann contractions later demonstrate that repeated D2H/H2D movement of eigenpairs is material, introduce the narrowest residency-aware mechanism required by the eigensolver + Lehmann pair.

Do not generalize it across the whole application.

---

# 12. MPI policy

Primary design points:

```text
single MPI rank + single GPU
```

and later:

```text
one MPI rank per GPU
```

Use the existing node-local rank/device-index machinery.

Do not add NCCL or CUDA-aware MPI during the first reciprocal implementation.

Multi-rank numerical equivalence should be established only after the one-rank/one-GPU route is stable.

---

# 13. CUDA/HIP policy

Implement CUDA first.

The Fortran reciprocal interface should remain vendor-neutral enough that a HIP backend can be introduced later if there is actual demand.

Do not introduce a portability framework merely to anticipate HIP.

The CUDA backend will require dense eigensolver support, so the build must add cuSOLVER cleanly and conditionally.

Do not impose CUDA/cuSOLVER dependencies on CPU-only builds.

---

# 14. Existing RS CUDA policy

The existing RS accelerator implementation is not to be rewritten.

Phase III-A should first determine exactly which existing RS GPU functionality is scientifically established.

Therefore:

1. inventory actual supported GPU kernels;
2. establish CPU<->GPU correctness for them;
3. profile them;
4. only then decide whether new RS kernels are needed.

This is especially important for KPM/Kubo-Bastin transport, where existing moment/velocity GPU machinery may already cover much of the expensive work.

---

# 15. CPU performance remains relevant

The LAPACK reciprocal backend is the permanent reference implementation.

Benchmark at least:

```text
A. serial k loop + threaded BLAS/LAPACK
B. parallel k loop + single-threaded BLAS/LAPACK
```

for representative small matrix sizes.

Do not add OpenMP parallelism until benchmarks show a clear benefit and nested-thread behavior is controlled.

GPU speedup should be reported against the best reasonable CPU reference, not an artificially slow configuration.

---

# 16. Tile-size policy

The backend already advertises preferred/maximum tile concepts.

Eventually the effective tile size should respect matrix dimension, eigenvectors vs eigenvalues only, backend preference, device memory, and solver API limitations.

Do not expose new user-facing tuning parameters until benchmarks demonstrate a need.

The first implementation may use a conservative fixed/internal tile policy.

---

# 17. Asynchrony policy

CUDA v1 should be synchronous and correct.

Do not begin with streams/double buffering.

Only add asynchronous overlap if profiling shows material benefit from overlapping CPU H(k) assembly, H2D, eigensolution, and D2H.

Correct synchronous code is the prerequisite.

---

# 18. Performance test tiers

Use three levels:

## Microbenchmark

Measures a narrow kernel or backend call.

Examples: batched eigensolver, H2D/D2H throughput, H(k) assembler, Lehmann contraction.

## Component benchmark

Measures a production component.

Examples: arbitrary-k eigenpair request, reciprocal SCF eigensystem phase, KPM moments, Lehmann GF for selected pairs/energies.

## End-to-end benchmark

Measures the real user workflow.

Examples: small/medium reciprocal SCF, representative transport calculation, representative Green-function post-processing.

Do not infer user speedup from microbenchmarks alone.

---

# 19. Definition of done for accelerator work

A work package is complete only when:

- [ ] CPU correctness baseline recorded;
- [ ] accelerator result satisfies the appropriate numerical contracts;
- [ ] invalid feature combinations were not hidden;
- [ ] no CPU scientific reference was silently changed;
- [ ] performance measured end-to-end where relevant;
- [ ] environment metadata recorded;
- [ ] memory/transfer overhead included where relevant;
- [ ] unsupported scope documented;
- [ ] no unrelated physics altered;
- [ ] focused tests pass;
- [ ] lean CPU test gate remains green;
- [ ] one focused commit created.

---

# 20. Work-package sequence

```text
ACC-00  benchmark harness
ACC-01  establish existing RS GPU coverage
ACC-02  reciprocal CUDA backend skeleton
ACC-03  isolated Hermitian GPU eigensolver
ACC-04  arbitrary-k integration
ACC-05  normal mesh / SCF integration
ACC-06  CPU-vs-GPU crossover and batching study
ACC-07  transfer/materialization reduction if justified
ACC-08  MPI rank-to-device establishment
ACC-09  validated reciprocal feature-scope extension
ACC-10  GPU Lehmann contractions
ACC-11  device-resident eigensystem handoff if justified
ACC-12  GPU H(k) assembly only if justified
ACC-13  RS KPM/transport GPU completeness and performance
ACC-14  accelerator support matrix and release gate
```

ACC-01 can run largely independently of ACC-02..12.

ACC-10 should wait until reciprocal eigensolution is established.

ACC-07, ACC-11, and ACC-12 are explicitly evidence-gated: they may finish with a documented decision not to add code.

---

# 21. Common Luna header

Prepend this section to every ACC prompt.

```text
You are working on the current HEAD of the fable_v3 branch of:

https://github.com/rslmtoasa/rslmtoasa

This task belongs to Phase III-A — Performance and Accelerator Establishment.

Before editing anything, inspect the CURRENT repository state, current CTest
labels, benchmark/test infrastructure, accelerator sources, and relevant
Phase-II validation evidence.

============================================================
GLOBAL ACCELERATOR RULES
============================================================

1. PRESERVE SCIENTIFIC CORRECTNESS

The validated CPU production route is the scientific reference for the
accelerated scope.

Do not change CPU reference values merely because GPU output differs.

2. ROOT CAUSES, NOT SYMPTOMS

Do not hide accelerator discrepancies by:
- loosening tolerances without analysis;
- silently falling back to CPU;
- disabling valid combinations;
- zeroing/clipping values;
- changing references;
- skipping failing data.

Identify the first violated numerical contract.

3. PROTECT LEGACY ATOMIC LMTO PHYSICS

Do not modify mature atomic/radial/potential routines in self.f90 and related
files for accelerator convenience.

4. USE EXISTING ARCHITECTURE

Prefer the existing reciprocal execution backend and existing RS CUDA plugin.

Do not create:
- a global GPU manager;
- generic device-data providers;
- new solver hierarchies;
- application-wide device arrays.

5. CPU-FIRST BASELINE

Before editing:
- run the smallest relevant correctness tests;
- record CPU observables;
- record performance baseline where relevant;
- record build/runtime environment.

6. PERFORMANCE CLAIMS

A kernel speedup is not an end-to-end speedup.

Include transfers, setup, synchronization, and materialization where relevant.

Compare against a reasonable CPU configuration.

7. TEST ORACLE POLICY

Do not implement a second physics algorithm in Python.

Python may:
- drive benchmarks;
- parse outputs;
- compute residual norms;
- compare timings;
- produce benchmark tables.

8. CUDA v1 SCOPE DISCIPLINE

Unless this task explicitly extends scope, do not add:
- GBT-specific support;
- TDDFT-specific support;
- GPU H(k) assembly;
- generalized eigensolvers;
- multi-GPU collectives;
- HIP;
- asynchronous stream pipelines.

9. GBT/TDDFT

GBT and TDDFT remain Development tracks.

Shared validated infrastructure may later serve them, but do not use them as
the accelerator correctness oracle in Phase III-A.

10. ONE FOCUSED COMMIT

When complete:
- tick every completed checkbox;
- leave genuinely incomplete items unticked;
- list files changed;
- list tests and benchmarks run;
- record correctness evidence;
- record performance evidence;
- state remaining unsupported scope;
- make one focused commit with the supplied one-line commit message.
```

---

# ACC-00 — Build the performance benchmark harness

## Objective

Create a stable benchmark harness that can measure CPU and later GPU execution without turning timing noise into a correctness gate.

This task should land **before new reciprocal CUDA kernels**.

## Scope

Inspect and extend the existing test/performance infrastructure rather than creating an unrelated benchmark framework.

The harness should support at minimum:

- metadata capture;
- repeated runs;
- warm-up runs;
- component timing records;
- machine-readable output;
- comparison of two benchmark runs;
- optional human-readable summary.

Use JSON or another simple repository-native format.

## Benchmark classes

Support labels/categories such as:

```text
performance
microbenchmark
component
end_to_end
rs
reciprocal
eigensolver
fourier
lehmann
transport
```

Do not make benchmark runtime thresholds normal CTest correctness assertions.

## Required initial CPU benchmarks

Add representative CPU baselines for reciprocal H(k) assembly where cleanly separable, eigensystem requests, arbitrary-k eigenpair requests, and the normal reciprocal eigensystem phase.

For RS add representative Block, Lanczos, Chebyshev, and KPM/moment timings where existing production entry points permit clean measurement.

Do not force all of these into one executable.

## Fixtures

At minimum cover:

- Si/sp;
- bcc Fe/spd;
- one multi-site metallic reciprocal fixture.

For RS use existing established metallic and Si fixtures as appropriate.

## Metrics

Record at minimum:

- wall time;
- repetition count;
- matrix dimension;
- number of k points;
- tile size;
- thread count;
- rank count.

Design the schema so later GPU runs can add H2D, D2H, eigensolver kernel time, device allocation/preparation time, and GPU model.

## Benchmark comparison

Provide a small script/tool that reports median, minimum, spread, speedup, and environment mismatch warnings.

Do not fail because performance differs by a fixed percentage.

## Checklist

- [x] Existing performance infrastructure audited
- [x] Benchmark schema defined
- [x] Environment metadata captured
- [x] Warm-up policy supported
- [x] Repetition/median supported
- [x] Machine-readable output implemented
- [x] Human-readable comparison implemented
- [x] Si reciprocal benchmark added
- [x] Fe reciprocal benchmark added
- [x] multi-site reciprocal benchmark added where practical
- [x] representative RS benchmark(s) added
- [x] no timing threshold added to ordinary correctness CI
- [x] benchmark usage documented
- [x] CPU baseline report recorded

## ACC-00 completion record

Implementation is in `tests/benchmarks/benchmark_harness.py`, with the initial
production inventory in `tests/benchmarks/manifest.json`, usage documentation
in `tests/benchmarks/README.md`, and the recorded CPU evidence in
`docs/dev/ACC-00_CPU_BASELINE.md`.

The focused tooling test is `BenchmarkHarnessSchema`; it validates profile
parsing, schema output, repetition capture, and comparison warnings without a
timing threshold. The completed CPU campaign covered the reciprocal profile,
Si/sp, bcc-Fe/spd, Block, Lanczos, and Chebyshev routes with one warm-up and
three repetitions. No GPU or cuSOLVER scope is claimed by ACC-00.

**Commit message:** `Add accelerator performance benchmark harness`

---

# ACC-01 — Establish existing RS CUDA correctness coverage

## Objective

Audit the existing RS CUDA plugin against the currently validated CPU RS functionality and close the gap between implemented GPU API and tested GPU API.

Do not rewrite the RS GPU architecture.

## Inventory

Map every public GPU route in the current RS accelerator layer, including where present:

- scalar Lanczos;
- Block Lanczos;
- Block DOS/GF;
- Chebyshev DOS/GF;
- Chebyshev moments;
- stochastic moments;
- orbital moments;
- velocity/current-related data;
- HOH/CCOR combinations actually supported.

For each record the CPU production counterpart, input limitations, current GPU tests, and missing correctness evidence.

## Validation

Use the strongest corresponding CPU contract.

Compare physical outputs, not only execution success: DOS, selected GF entries, recursion-derived moments, conductivity/KPM moments, and orbital moments where supported.

Use appropriate tolerances for GPU floating-point ordering. Do not demand bitwise equality.

## Unsupported combinations

If the plugin does not implement a CPU-supported combination, document it.

Do not add a CPU fallback and call the GPU feature supported.

## Performance

Run the ACC-00 harness on the established GPU routes and record which existing GPU kernels already deliver useful end-to-end benefit.

## Checklist

- [x] RS GPU API inventoried
- [x] CPU counterpart mapped for each route
- [x] existing GPU tests mapped
- [x] Chebyshev CPU/GPU contract established on GPU hardware
- [x] Block CPU/GPU contract established on GPU hardware
- [x] Lanczos CPU/GPU contract established on GPU hardware where supported
- [x] GF outputs compared where supported
- [x] KPM/moment outputs compared where supported
- [x] HOH/CCOR GPU scope documented
- [x] unsupported combinations not hidden by fallback
- [x] benchmark results recorded on GPU hardware
- [x] support matrix updated

## ACC-01 completion record

The inventory and evidence matrix is in
`docs/dev/ACC-01_RS_CUDA_COVERAGE.md`. The no-GPU source/API contract is
`tests/cuda/test_plugin_surface.py` and is registered as `CudaPluginSurface`
when standalone unit tests are enabled. The executable low-level contract is
`tests/cuda/rsrec_validate.py`, built and launched by
`tests/cuda/build_and_validate.sh`; it now binds every public C argument and
compares recursion, orbital, DOS, and eta-GF outputs with NumPy references.
ACC-01 also closes two dispatch issues found during the audit: structured
FFT/conv orbital moments are rejected before entering an implementation that
cannot support them, and Green reconstruction falls back when the executable
was built without the CUDA plugin.

ACC-01 is complete. The user-provided hardware
campaign passed all low-level CUDA routes (15 reported comparisons) at
approximately 1e-15 relative error, the CPU regression matrix at 10/10, and the production Chebyshev
CPU/GPU consistency matrix at 8/8. The ACC-00 harness then measured 2.872801
s for RS Block Fe and 1.922788 s for RS Chebyshev Si with `gpu_plugin=true`
on an NVIDIA RTX A4000. Compared with the same-build CPU medians, these are
1.14x and 5.73x speedups respectively. The build reported
`ENABLE_MKL_KERNELS=OFF`; optional MKL backend cases were skipped.

**Files changed:** `source/recursion_core.f90`, `source/green_chebyshev.f90`,
`source/green_block.f90`, `tests/cuda/test_plugin_surface.py`,
`tests/cuda/rsrec_validate.py`, `tests/cuda/build_and_validate.sh`,
`tests/run_test.py`, `tests/benchmarks/benchmark_harness.py`,
`tests/benchmarks/manifest.json`, `tests/benchmarks/README.md`,
`docs/dev/ACC-01_RS_CUDA_COVERAGE.md`, `docs/source/theory/gpu_acceleration.rst`,
and this blueprint.

**Checks run:** `python3 tests/cuda/test_plugin_surface.py` (pass), Python
syntax compilation (pass), CUDA compilation through the project CUDA build
(pass), standalone numerical validator (15 reported comparisons pass), CPU regression matrix
(10/10 pass), production CPU/GPU consistency matrix (8/8 pass), and the
GPU-specific ACC-00 Block/Chebyshev benchmark (3 repetitions each). The
optional MKL cases were skipped because `ENABLE_MKL_KERNELS=OFF`.

**Commit message:** `Establish RS CUDA correctness coverage`

---

# ACC-02 — Add reciprocal CUDA backend skeleton

## Objective

Introduce a CUDA reciprocal execution backend behind the existing reciprocal backend interface, without yet implementing the full eigensolver path.

This task establishes build, ownership, capability reporting, and lifecycle only.

## Build

Add CUDA/cuSOLVER dependencies conditionally.

CPU-only builds must remain unaffected.

Detect and document the minimum CUDA/cuSOLVER capabilities actually required.

ACC-02 adds the optional `ENABLE_CUDA_RECIPROCAL` build switch. Its CUDA
target requires the CUDA runtime and cuSOLVER libraries; no cuSOLVER solver
entry point is selected or called until ACC-03. The existing minimum toolkit
policy is unchanged.

Do not unnecessarily raise the toolkit minimum version before the eigensolver API decision in ACC-03 is complete.

## Backend

Create the smallest concrete backend implementation compatible with the existing reciprocal execution interface.

Implement initialize, capabilities, prepare_operator, synchronize, release, metrics, and a clear unsupported execution result if solve is not yet implemented.

Do not create a second public reciprocal API.

## Device ownership

Use explicit backend-owned state. No global CUDA singleton.

Use existing MPI/device-index information where available, but one rank/one GPU is sufficient in this task.

## Operator generation

Make persistent device preparation generation-aware.

Do not invalidate/re-upload an unchanged operator merely because a backend accessor is called again.

## Checklist

- [x] CPU-only build unchanged
- [x] CUDA build discovers cuSOLVER conditionally
- [x] reciprocal CUDA backend type created
- [x] existing reciprocal public API unchanged
- [x] capabilities implemented
- [x] initialize/release implemented
- [x] synchronize implemented
- [x] operator-generation reuse implemented
- [x] no hidden global CUDA state introduced
- [x] unsupported solve reports clearly
- [x] LAPACK backend unchanged
- [x] focused lifecycle tests pass

**Commit message:** `Add reciprocal CUDA backend skeleton`

## ACC-02 completion record

The optional `ENABLE_CUDA_RECIPROCAL` target links CUDA runtime and cuSOLVER
conditionally and leaves `ENABLE_CUDA_PLUGIN` independent. The concrete
`cuda_reciprocal_backend` uses the existing reciprocal execution interface;
its opaque backend-owned context owns the selected device and CUDA stream.
The node-local MPI rank is mapped through `g_parallel_context%device_index`.

Operator preparation stores the generation in the opaque context and returns
reuse for repeated preparation of the same generation. The execution method
clears the result and reports that solve is unsupported, with no CPU fallback.
Capabilities therefore advertise device residency but no numerical solve
feature until ACC-03.

**Files changed:** `CMakeLists.txt`, `source/CMakeLists.txt`,
`source/reciprocal.f90`, `source/reciprocal_backend.f90`,
`source/cuda/reciprocal_cuda.h`, `source/cuda/reciprocal_cuda.cpp`, and
`tests/unit/test_reciprocal_cuda_lifecycle.cpp`.

**Checks run:** CPU-only configure/build with `RUN_UNIT_TESTS=ON` (pass),
CUDA/cuSOLVER configure/build with `ENABLE_CUDA_RECIPROCAL=ON` (pass),
`UnitArbitraryKEigenpairs` (pass), `ReciprocalCudaLifecycle` (CTest pass;
hardware-aware skip without a CUDA driver), and `git diff --check` on the
ACC-02 files (pass).

No eigensolver correctness or performance result is claimed by ACC-02;
those belong to ACC-03 and later.

---

# ACC-03 — Implement standard Hermitian GPU eigensolution

## Objective

Implement the first actual reciprocal CUDA numerical kernel:

> many independent standard complex-Hermitian eigensystems.

Do not integrate it into every reciprocal workflow yet.

## Scope

Support:

- standard Hermitian `H v = e v`;
- eigenvalues;
- eigenvectors;
- host input Hamiltonian tile;
- host output eigenpairs;
- synchronous execution;
- single GPU.

Do not support generalized overlap in this task.

## Solver API selection

Use ACC-00 benchmark metadata and the installed CUDA toolkit to compare plausible cuSOLVER strategies for the small-matrix regime.

Evaluate the appropriate conventional dense Hermitian route and batched/uniform-batched route where available.

Do not choose solely from API novelty.

Use real LMTO matrices in addition to synthetic Hermitian microbenchmarks.

## Numerical contracts

For each matrix check eigenvalue agreement, residual norm, orthonormality, and projector comparison for degenerate groups.

Test dimensions representative of Si/sp, bcc Fe/spd, and an intermediate multi-site matrix.

## Workspace

Allocate reusable device workspace per backend/tile rather than malloc/free per k point.

Do not over-engineer a general memory pool.

## Checklist

- [x] standard Hermitian GPU solver implemented
- [x] eigenvalues supported
- [x] eigenvectors supported
- [x] reusable workspace implemented
- [x] no generalized solver added
- [x] Si-size matrix dimension covered by the focused synthetic test
- [x] Fe-size matrix dimension covered by the focused synthetic test
- [x] intermediate matrix dimension covered by the focused synthetic test
- [x] residual norms pass
- [x] orthogonality passes
- [x] degenerate-subspace comparison handled correctly
- [x] conventional-vs-batched API choice evaluated and documented
- [x] CPU reference unchanged
- [ ] real-LMTO GPU validation and performance results recorded

**Commit message:** `Implement reciprocal CUDA Hermitian eigensolver`

## ACC-03 implementation record

The CUDA backend now exposes a narrow standard-Hermitian tile solve through
`rslmto_reciprocal_cuda_solve_zheevd_batch`. It uses the installed CUDA 13.3
`cusolverDnZheevd` divide-and-conquer solver for each independent matrix in a
tile, while retaining one cuSOLVER handle, stream, device matrix/eigenvalue
buffers, solver workspace, and info buffer per backend context. H2D and D2H
transfers occur at the synchronous tile boundary. Eigenvalues-only requests
are supported; generalized overlap and H(k) assembly remain unsupported until
later work packages.

The uniform-batched Jacobi API (`cusolverDnZheevjBatched`) was reviewed but not
selected for the first implementation: `Zheevd` preserves the conventional
double-complex dense route and its workspace/lifecycle contract is clear for
the current tile seam. A GPU timing comparison remains deliberately pending
because this development host has CUDA 13.3 headers/libraries but no usable
NVIDIA driver. Consequently no real-LMTO GPU accuracy or speedup claim is made
by ACC-03; the hardware campaign should use the ACC-00 metadata and compare
this route with the batched API before adopting a batching policy.

The focused synthetic contract covers dimensions 8 (Si/sp-sized), 18
(bcc-Fe/spd-sized), and 36 (intermediate), plus a degenerate 4x4 case. It
checks known eigenvalues, residuals, orthonormality, eigenvalues-only output,
and the degenerate eigenspace projector. CPU-only and CUDA-enabled builds both
compile successfully; the CUDA CTest targets are hardware-aware and skip when
the driver is unavailable.

**Files changed:** `CMakeLists.txt`, `source/reciprocal.f90`,
`source/reciprocal_backend.f90`, `source/cuda/reciprocal_cuda.h`,
`source/cuda/reciprocal_cuda.cpp`, `tests/unit/test_reciprocal_cuda_eigensolver.cpp`,
and this blueprint.

**Checks run:** CUDA-enabled `rslmto` build (pass), CPU-only `rslmto` build
(pass), `ReciprocalCudaLifecycle` (pass/driver-aware skip),
`ReciprocalCudaEigensolver` (pass/driver-aware skip), and `git diff --check`
(pass). Real-GPU numerical and performance checks remain an explicit ACC-03
follow-up when CUDA hardware is available.

---

# ACC-04 — Integrate CUDA eigensolution into arbitrary-k requests

## Objective

Make the existing arbitrary-k reciprocal eigenpair API use the CUDA backend without changing its physics contract.

This is the preferred first end-to-end reciprocal integration because arbitrary-k requests already operate tile-wise and do not need the H(k) compatibility cache.

## Path

Preserve:

```text
k points
 -> existing CPU H(k) assembly
 -> host tile
 -> CUDA eigensolver
 -> host eigenpairs
 -> existing caller
```

Do not GPU-port Fourier assembly.

## Backend selection

Use the existing backend-selection mechanism.

Do not scatter CUDA conditionals through reciprocal physics routines.

If no CUDA backend is requested/available, LAPACK behavior must remain unchanged.

## Correctness fixtures

Use Si/sp arbitrary-k points, bcc Fe/spd arbitrary-k points, and duplicate/folded k-point behavior already tested by the CPU route.

Compare eigenvalues, residuals, and projectors/subspaces where required.

## Performance

Benchmark assembly, H2D, solve, D2H, and total arbitrary-k request. Record tile-size dependence.

## Checklist

- [x] arbitrary-k caller uses backend cleanly
- [x] CPU H(k) assembly retained
- [x] no H(k) compatibility return added unnecessarily
- [ ] Si arbitrary-k CPU/GPU matches
- [ ] Fe arbitrary-k CPU/GPU matches
- [x] duplicate/folded-k behavior preserved on the CPU reference route
- [ ] timings decomposed
- [x] LAPACK path unchanged
- [x] no CUDA conditionals leaked into physics code
- [ ] end-to-end benchmark recorded

**Commit message:** `Enable CUDA arbitrary-k eigensolution`

## ACC-04 implementation record

The arbitrary-k service now queries the typed execution-backend capabilities
before submitting each deduplicated tile. LAPACK retains its existing combined
CPU Fourier-assembly/eigensolution request. A CUDA v1 backend, which advertises
host-H(k)-input but no device Fourier assembly, receives a request-owned host
tile assembled by the existing `reciprocal_assembler`; no CUDA conditionals or
backend-name tests enter the reciprocal physics path. The service synchronizes
at the result boundary and rejects incomplete or unsupported feature results
instead of silently falling back to CPU.

The existing arbitrary-k CPU fixture continues to cover exact folded-point
duplicates, multi-tile scattering, and the Si-sized one-site route. The CPU
baseline profile recorded on this development host reports the arbitrary-k
assembly/eigensolution phase as 1.9700e-4 s for bcc Fe/spd-sized matrices
(18 x 18, 16 points) and 4.2840e-3 s for the two-site 36 x 36 profile (32
points). These are CPU reference measurements, not GPU speedup claims.

`ReciprocalCudaArbitraryK` now exercises the same arbitrary-k service with a
real CUDA context, two-point tiles, exact folded duplicates, CPU/GPU
eigenvalue comparison, projector comparison, residuals, and orthonormality.
On the external NVIDIA run (CUDA 13.3, driver 610.57.04, RTX A4000), the
worst eigenvalue, projector, and residual errors were respectively
5.551115e-17, 0.0, and 5.551115e-17. This is focused backend-integration
evidence; the real assembled Si/sp and bcc-Fe/spd production fixtures and the
H2D/solve/D2H timing decomposition remain pending.

**Files changed:** `CMakeLists.txt`, `source/reciprocal_backend.f90`,
`source/reciprocal_fourier.f90`, `tests/unit/test_acc04_arbitrary_k_source.py`,
`tests/unit/test_reciprocal_cuda_arbitrary_k.f90`, and this blueprint.

**Checks run:** CPU `UnitArbitraryKEigenpairs` and `Acc04ArbitraryKSource`
(pass), CUDA-enabled build with `ENABLE_CUDA_RECIPROCAL=ON` (pass), CUDA
`UnitArbitraryKEigenpairs`, `Acc04ArbitraryKSource`, and
`ReciprocalCudaArbitraryK` (pass on external hardware), and the
`ReciprocalCudaLifecycle`/`ReciprocalCudaEigensolver` checks (pass on external
hardware; driver-aware skip in the sandbox). `git diff --check` (pass).

---

# ACC-05 — Integrate CUDA eigensolution into the normal reciprocal mesh

## Objective

Use the CUDA backend for the ordinary reciprocal eigensystem path used by validated SCF/bands/DOS workflows while preserving current host semantics.

Do not optimize away H(k) materialization yet.

## Initial semantic rule

CUDA v1 must return the same host-side products currently expected by normal reciprocal workflows, including the compatibility H(k) cache where requested.

Correctness comes before transfer minimization.

## Validate

### Si/sp
Check electron count, EF, band energy, charge, band eigenvalues, and DOS/state count.

### bcc Fe/spd
Also check magnetic moment.

Run one controlled SCF/eigensystem path and the lean converged functional fixture where practical.

## Performance

Measure CPU H(k), H2D, GPU eigensolve, D2H eigenpairs, D2H H(k) compatibility transfer if applicable, total reciprocal solve phase, and end-to-end SCF workflow.

## Checklist

- [x] normal mesh CUDA execution wired
- [x] host semantics preserved
- [x] Si SCF CPU/GPU matches
- [x] Fe SCF CPU/GPU matches
- [x] electron count matches
- [x] EF matches
- [x] charge matches
- [x] Fe moment matches
- [x] band/DOS checks pass
- [x] transfer costs measured
- [x] LAPACK remains selectable
- [x] no reference values regenerated without cause

**Commit message:** `Enable CUDA reciprocal mesh eigensolution`

## ACC-05 implementation record

The normal-mesh execution seam now queries typed backend capabilities before
submitting each tile. LAPACK retains the existing combined host assembly and
eigensolution request. CUDA v1 instead receives a request-owned host H(k) tile,
solves the standard Hermitian systems on the device, and returns host
eigenpairs. The existing `hk_bulk`, `eigenvalues`, and `eigenvectors` caches are
populated exactly as before; no H(k) materialization was removed.

CUDA input-Hamiltonian solve requests are counted separately from combined
requests. A focused CUDA integration fixture now covers a tiled ordinary mesh,
CPU/CUDA cache and eigenpair comparisons, and the request-shape contract. The
fixture passes on the development GPU: NVIDIA RTX A4000, CUDA 13.3, driver
610.57.04, with worst eigenvalue, projector, and residual errors of
5.551115E-17, 0.0, and 5.551115E-17 respectively.

Production validation used the same spglib-enabled CPU baseline and CUDA
backend on 4x4x4 Gaussian-mesh Si/sp and bcc-Fe/spd cases, with both a
three-step controlled SCF and the nstep=12 converged fixture. The converged
Si results agree within 9.2E-7 Ry in EF, 1.6E-5 Ry in band energy, 1.3E-6
Ry in total energy, and 3E-6 in site charge; electron count is 8.0 and the
DOS state count is 16.0. bcc Fe agrees at reported precision for EF,
electron count, band energy, total energy, site charge, DOS state count, and
the 1.950775 Bohr-magneton site moment. The band-energy/eigensystem and DOS
checks therefore pass without changing scientific references.

CUDA phase timing was recorded on the RTX A4000 (CUDA 13.3, driver
610.57.04). At the final 48-tile point in the nstep=12 runs, the reciprocal
phase reported H2D 7.9E-4 s (Si) / 9.6E-4 s (Fe), GPU eigensolve 4.7E-1 s /
5.5E-1 s, D2H eigenpairs 1.1E-3 s / 1.4E-3 s, and total reciprocal phase
3.8E-2 s / 4.3E-2 s. The matched CPU baseline reported approximately
2.0E-3 s / 1.0E-3 s for host H(k) assembly. H(k) is assembled on the host and materialized directly
into the compatibility cache, so D2H H(k) is 0/not applicable. End-to-end
two-case campaign wall time was 3.71 s for CPU and 7.73 s for CUDA; this is
correctness evidence, not a claim of a speedup for these small matrices.

**Files changed:** `source/reciprocal.f90`,
`source/reciprocal_backend.f90`, `source/reciprocal_fourier.f90`,
`source/cuda/reciprocal_cuda.cpp`, `source/cuda/reciprocal_cuda.h`,
`source/include_codes/namelists/reciprocal.f90`,
`tests/unit/test_reciprocal_cuda_arbitrary_k.f90`,
`tests/validation/val02_reciprocal_scf.py`, `CMakeLists.txt`, and this
blueprint.

**Checks run:** CUDA-enabled build with `ENABLE_CUDA_RECIPROCAL=ON` (pass),
`UnitArbitraryKEigenpairs` and `Acc04ArbitraryKSource` in the CUDA-enabled
build (pass), `ReciprocalCudaLifecycle`, `ReciprocalCudaEigensolver`, and
`ReciprocalCudaArbitraryK` including the ACC-05 fixture (pass on the RTX
A4000), CPU-only build plus `UnitArbitraryKEigenpairs` and
`Acc04ArbitraryKSource` (pass), non-fused CPU compatibility build
`UnitArbitraryKEigenpairs` (pass), spglib-enabled CPU/CUDA production
controlled and nstep=12 Si/Fe campaigns (all four cases pass), CUDA phase
timing and end-to-end measurements (pass), and focused `git diff --check`
(pass).

---

# ACC-06 — Determine CPU/GPU crossover and tile policy

## Objective

Use the benchmark harness to establish when the reciprocal CUDA eigensolver is actually beneficial.

This task is primarily measurement and policy.

Do not add major new numerical functionality.

## Matrix-size study

Use real LMTO matrices spanning very small Si/sp, bcc Fe/spd, multi-site intermediate, and one larger existing cell.

For each sweep k count, tile size, eigenvalues-only vs eigenvectors where supported, CPU thread configuration, and GPU batching strategy.

## CPU strategies

Compare at least:

A. scalar k loop + threaded BLAS/LAPACK  
B. parallel k loop + controlled single-threaded BLAS/LAPACK, if a minimal benchmark implementation can be tested without destabilizing production code.

Do not merge OpenMP code into production solely for the benchmark unless evidence is already strong.

## Deliverable

Produce a crossover table such as:

```text
matrix size | Nk | best CPU | GPU | end-to-end speedup | recommended backend
```

Determine default tile-size policy, whether backend-reported preferred tile size should influence the caller, and whether automatic CPU/GPU dispatch is justified.

Do not add automatic dispatch merely because a crossover exists; explicit backend selection may remain preferable initially.

## Checklist

- [ ] representative real matrix sizes benchmarked
- [ ] Nk dependence measured
- [ ] tile-size dependence measured
- [ ] eigenvector cost separated
- [ ] best reasonable CPU baseline established
- [ ] GPU end-to-end cost recorded
- [ ] crossover documented
- [ ] preferred tile policy documented
- [ ] automatic-dispatch decision justified
- [ ] no unsupported performance claim made

**Commit message:** `Benchmark reciprocal CPU and CUDA crossover`

## ACC-06 implementation record

The ACC-06 driver is an opt-in benchmark target, not a CTest timing gate. It
uses the production reciprocal H(R)->H(k) assembler to create deterministic
Hermitian LMTO-shaped fixtures at the dimensions of the small Si/sp and bcc
Fe/spd cases, a two-site intermediate case, and a four-site larger case. The
fixtures exercise the validated dense standard-Hermitian request boundary;
they are performance fixtures, not regenerated physical references.

The generated full campaign report is committed in
[`docs/dev/ACC-06_CPU_GPU_CROSSOVER.md`](ACC-06_CPU_GPU_CROSSOVER.md).

`tests/benchmarks/acc06_crossover.py` drives the ACC-00 harness and records
one warm-up plus three samples for each combination of matrix size, Nk in
`{1, 8, 32}`, tile size in `{1, 8, 16}`, and eigenvalues-only/eigenvectors
mode. The full campaign contained 216 executable runs and produced 72
CPU-vs-GPU rows. GPU timing includes H2D, cuSOLVER, D2H, synchronization, and
the harness wall time includes executable/backend initialization. The driver
also reports assembly and backend solve phase times separately.

CPU strategy A is the production typed LAPACK backend with its serial k/tile
loop and threaded BLAS/LAPACK. Strategy B is a benchmark-only OpenMP loop
with one independent `zheev` workspace per k point. It is included as an
exploratory baseline only; no OpenMP loop was added to production reciprocal
code.

### Representative crossover table

The following rows are the full-campaign medians at Nk=8, tile=8, with
eigenvectors requested, on GNU Fortran 13.3.0, Intel oneMKL, two OpenMP
threads, one MPI rank, and two NVIDIA RTX A4000 devices visible. The CUDA
campaign used driver 610.57.04 and CUDA toolkit 13.3.

| matrix size | Nk | best CPU | GPU | end-to-end speedup | recommended backend |
|---:|---:|---:|---:|---:|---|
| 8 (Si/sp) | 8 | 0.009279 s | 0.399418 s | 0.02x | LAPACK |
| 18 (bcc Fe/spd) | 8 | 0.009149 s | 0.409600 s | 0.02x | LAPACK |
| 36 (two-site spd) | 8 | 0.013668 s | 0.424761 s | 0.03x | LAPACK |
| 72 (four-site spd) | 8 | 0.013675 s | 0.444491 s | 0.03x | LAPACK |

No CPU/GPU crossover was observed in the studied 8--72 matrix, Nk=1--32
range. The best GPU/CPU ratio in the complete table was 0.05x, so this
campaign provides no evidence for automatic CUDA dispatch or a production GPU
default. This is a bounded result for the measured fixtures and does not
claim that larger production matrices cannot cross over.

### Tile and eigenvector policy

Eigenvector cost was separated by running both request modes. At Nk=32 the
CUDA internal phase generally favored tile 16 for the tested sizes, while
Nk=1 and Nk=8 results moved between tile 1, 8, and 16 within the run-to-run
noise. The stable policy is therefore to retain the existing caller default
tile size of 16 and keep the CUDA capability's `preferred_tile_size=1`
neutral rather than advertise an unvalidated universal preference. No tile
autotuning or automatic CPU/GPU dispatch was added. Explicit backend
selection remains the policy until a larger, production-weighted workload
shows a material end-to-end crossover.

### Reproduction and checks

```bash
cmake -S . -B build-acc00
cmake --build build-acc00 --target ReciprocalCrossoverBenchmark --parallel
cmake -S . -B build-acc01-cuda
cmake --build build-acc01-cuda --target ReciprocalCrossoverBenchmark --parallel
python3 tests/benchmarks/test_benchmark_harness.py
python3 tests/benchmarks/acc06_crossover.py \
  --cpu-binary build-acc00/bin/ReciprocalCrossoverBenchmark \
  --cpu-build-dir build-acc00 \
  --cuda-binary build-acc01-cuda/bin/ReciprocalCrossoverBenchmark \
  --cuda-build-dir build-acc01-cuda \
  --output-dir /tmp/rslmto-acc06-full \
  --report /tmp/rslmto-acc06-full.md \
  --warmups 1 --repetitions 3
```

The CUDA executable and existing reciprocal correctness tests passed on the
RTX A4000 campaign host. The CPU benchmark target and CPU parser contract
also passed. Timing JSON and Markdown remain machine-local evidence under
`/tmp`; no performance result was promoted to a correctness reference.

- [x] representative production-assembler LMTO matrix dimensions benchmarked
- [x] Nk dependence measured
- [x] tile-size dependence measured
- [x] eigenvalue-only versus eigenvector cost separated
- [x] best reasonable CPU baseline established
- [x] GPU end-to-end cost recorded
- [x] crossover documented for the measured range
- [x] preferred tile policy documented
- [x] automatic-dispatch decision justified
- [x] no unsupported performance claim made

---

# ACC-07 — Reduce unnecessary host materialization if profiling justifies it

## Objective

Inspect the validated normal reciprocal workflows and determine whether returning/materializing host H(k) is a meaningful performance or memory cost.

This task is evidence-gated.

## First step

Profile current CUDA normal-mesh execution from ACC-05/06.

Determine the fraction of time spent transferring H(k) back, host memory consumed by the full H(k) cache, and actual downstream consumers requiring H(k).

Map every consumer before changing semantics.

## Decision A — not material

If H(k) transfer/materialization is not significant:

- do not change production code;
- document the result;
- close the task.

This is a successful outcome.

## Decision B — material

If significant, introduce the smallest explicit request distinction needed, for example:

```text
eigensystem-only
eigensystem + H(k) compatibility cache
```

or reuse an existing request flag if already present.

Do not make H(k) laziness implicit or surprising.

Preserve old behavior for callers that require it.

## Checklist

- [x] H(k) consumers mapped
- [x] H(k) transfer measured
- [x] H(k) host memory measured
- [x] optimization threshold justified
- [x] no code added because the measured benefit is insignificant
- [x] no request semantics changed because no optimization was justified
- [x] legacy callers preserved
- [x] SCF/bands/DOS unchanged
- [x] benchmark repeated after the Decision-A assessment (no code change)
- [x] net end-to-end benefit not applicable because no production code changed

**Commit message:** `Reduce reciprocal Hk materialization overhead`

## ACC-07 completion record

ACC-07 is closed with Decision A. The consumer map, transfer/memory
measurements, thresholds, fresh CPU/CUDA validation, and repeat benchmark are
recorded in [`docs/dev/ACC-07_HK_MATERIALIZATION.md`](ACC-07_HK_MATERIALIZATION.md).

The validated CUDA normal-mesh path reports zero H(k) device-to-host transfer:
host Fourier assembly produces the request tile and the compatibility cache is
populated from that same host data. H(k) remains required by the reciprocal
Dyson Green and BSF consumers, so removing it would change established caller
semantics. The measured 4x4x4 cache payloads are 0.0625 MiB for Si/sp and
0.3164 MiB for bcc-Fe/spd; the ACC-06 72x72, Nk=32 larger fixture is 2.5313
MiB. No production code or request semantics were changed.

Fresh current-build CUDA and LAPACK converged Si/Fe SCF/DOS runs matched at
reported precision, including EF, band energy, total energy, DOS integral, and
the Fe moment. The representative ACC-06 repeat remained CPU-favored across
all four fixtures, so no end-to-end materialization benefit was demonstrated
or claimed.

- [x] H(k) consumers mapped
- [x] H(k) transfer measured
- [x] H(k) host memory measured
- [x] optimization threshold justified
- [x] no code added if benefit is insignificant
- [x] if optimized, request semantics explicit (not applicable; Decision A)
- [x] legacy callers preserved
- [x] SCF/bands/DOS unchanged
- [x] benchmark repeated after decision
- [x] net end-to-end benefit demonstrated if code changed (not applicable; no code changed)

**Checks run:** fresh CUDA and current-build LAPACK `val02_reciprocal_scf.py
--converged-only` campaigns (2/2 cases pass each), representative ACC-06
repeat (24/24 benchmark runs pass), and the focused reciprocal CUDA/CPU tests
listed in the handoff. `git diff --check` passes.

---

# ACC-08 — Establish MPI rank-to-GPU execution

## Objective

Validate the existing MPI/device-index model for accelerator execution.

Do not implement distributed eigensolvers.

## Scope

Support and validate:

```text
one MPI rank -> one selected GPU
```

for nodes with one or more GPUs.

Use existing node-local rank/device-index logic.

Do not add NCCL or CUDA-aware MPI.

## Tests

At minimum:

- single rank / single GPU;
- two ranks mapped deterministically if two GPUs are available;
- explicit device override if the code supports one.

Validate device selection, no accidental same-device oversubscription unless explicitly configured, and numerical equality with one-rank CPU/GPU reference for decomposable workflows.

If hardware with multiple GPUs is unavailable, add unit-level mapping coverage and document the unverified hardware gate.

## Checklist

- [x] existing device mapping audited
- [x] reciprocal CUDA backend uses existing device index
- [x] one-rank/one-GPU validated
- [x] multi-rank mapping validated where hardware permits
- [x] explicit override validated where supported
- [x] no NCCL added
- [x] no CUDA-aware MPI dependency added
- [x] numerical equivalence checked
- [x] unsupported multi-GPU scope documented

## ACC-08 completion record

ACC-08 is closed with the existing node-local `MPI_Comm_split_type` context as
the single rank/device source for both CUDA paths. The RS CUDA recursion plugin
no longer passes device 0 unconditionally. It queries the CUDA-visible device
count, maps `local_rank -> device_index`, and logs
`CUDA_DEVICE_MAPPING world_rank=... local_rank=... local_size=...
visible_devices=... selected_device=...`. The reciprocal CUDA backend uses the
same context and mapping helper.

Production mapping rejects unconfigured local oversubscription: if more local
MPI ranks than visible GPUs are present, initialization fails instead of
silently wrapping ranks onto a shared device. An intentional override is
available through `RSLMTO_CUDA_DEVICE`, interpreted in the
`CUDA_VISIBLE_DEVICES` namespace; it is range-checked and explicitly permits
shared-device execution. The pure mapping helper and MPI unit tests also cover
invalid counts, four ranks/four devices, four ranks/two devices, and programmatic
override cases.

On the validation host, two NVIDIA RTX A4000 GPUs were visible with CUDA 13.3,
driver 610.57.04, and Open MPI 4.1.6. The clean MPI+CUDA build passed the
one-rank RS mapping test, the two-rank RS mapping test, the reciprocal
arbitrary-k/normal-mesh CPU-vs-CUDA test in one- and two-rank modes, and the
reciprocal CUDA lifecycle test. The two-rank logs selected device 0 for local
rank 0 and device 1 for local rank 1 in both RS and reciprocal backends; CPU/GPU
eigenvalue, residual, orthogonality, and projector comparisons passed at the
existing ACC-04 tolerances. `RSLMTO_CUDA_DEVICE=1` was validated in both
backends, including intentional two-rank sharing. A four-rank launch without
the override failed with the expected unconfigured-oversubscription diagnostic.

The supported multi-GPU scope remains one MPI rank per selected GPU on a shared
node. There is no distributed eigensolver, GPU-to-GPU collective, cross-rank
CUDA data exchange, NCCL, or CUDA-aware MPI dependency. Cross-node placement
and scheduler-specific GPU allocation remain outside this task; use
`CUDA_VISIBLE_DEVICES` to define each process's visible-device namespace.

- [x] mapping helper and local communicator audited
- [x] RS CUDA plugin uses the shared device index
- [x] reciprocal CUDA backend uses the shared device index
- [x] one-rank/one-GPU hardware test passed
- [x] two-rank/two-GPU deterministic mapping passed
- [x] explicit `RSLMTO_CUDA_DEVICE` override passed
- [x] unconfigured oversubscription rejected
- [x] CPU/GPU numerical equivalence passed
- [x] no NCCL or CUDA-aware MPI dependency added
- [x] unsupported distributed/multi-node scope documented

**Checks run:** CPU Open MPI `UnitParallelContext` at 1, 2, and 4 ranks plus
the override test; CUDA/Open MPI `CudaMpiDeviceMapping` at one and two ranks;
`ReciprocalCudaArbitraryK` in serial and two-rank modes;
`ReciprocalCudaLifecycle`; reciprocal two-rank mapping with and without the
explicit override; and the expected four-rank/two-GPU oversubscription failure.
The ACC-08 source diff passes `git diff --check`.

**Commit message:** `Validate MPI rank-to-GPU execution`

---

# ACC-09 — Extend reciprocal CUDA support to validated operator variants

## Objective

Extend the reciprocal CUDA eigensolver path beyond the initial ordinary Hamiltonian scope only for operator variants already scientifically established on CPU.

Do not broaden scope by assumption.

## Audit first

Inventory currently validated reciprocal combinations involving second-order/HOH, CCOR, and other ordinary standard-Hermitian variants.

Separate:

- variants that still reduce to the same standard Hermitian eigensystem;
- genuinely generalized-overlap problems;
- scientifically Development combinations.

## Implementation

For standard-Hermitian variants, reuse the existing CPU assembly and feed the same CUDA eigensolver.

Do not duplicate operator construction on GPU.

For generalized problems, defer unless Phase-II evidence clearly supports them and the current task is explicitly expanded.

## Validation

Use corresponding CPU functional fixtures and compare physical downstream outputs, not only eigenvalues.

## Checklist

- [ ] reciprocal operator variants inventoried
- [ ] CPU maturity checked for each
- [ ] standard-Hermitian validated variants enabled
- [ ] generalized variants not accidentally claimed
- [ ] HOH/CCOR outputs compared where supported
- [ ] no GPU operator duplication added
- [ ] unsupported combinations documented
- [ ] benchmark impact recorded

**Commit message:** `Extend CUDA reciprocal operator coverage`

---

# ACC-10 — Implement GPU Lehmann Green-function contractions

## Objective

Accelerate the validated Lehmann Green-function contraction while preserving canonical `green` outputs and existing downstream consumers.

Do not create a separate GPU Green-function object.

## Baseline

Use existing validated Lehmann/Dyson Sigma=0 equivalence, direct selected G(z) validation, and Jij/conductivity/damping route triads.

Record CPU performance for representative k counts, energy counts, and site/pair counts.

## Initial data path

The first GPU Lehmann implementation may consume host eigenpairs copied to device.

Do not require device-resident eigensolver handoff in this task.

Parallelize the expensive contraction dimensions as appropriate: k, band, energy, pair.

Preserve Green-function sign conventions, pair indexing, complex-energy conventions, and canonical output arrays.

## Validation

Compare GPU and CPU selected onsite G, selected intersite G, all requested pair outputs for a tiny fixture, and downstream Jij/damping/transport triads where relevant.

## Performance

Measure eigenpair H2D, contraction, output D2H, and total Lehmann request.

## Checklist

- [ ] CPU Lehmann baseline recorded
- [ ] GPU contraction implemented
- [ ] canonical green arrays preserved
- [ ] onsite G matches
- [ ] intersite G matches
- [ ] pair indexing preserved
- [ ] energy conventions preserved
- [ ] downstream triads remain valid
- [ ] transfer cost measured
- [ ] end-to-end speedup recorded
- [ ] no device-residency framework introduced yet

**Commit message:** `Accelerate Lehmann Green functions on CUDA`

---

# ACC-11 — Add device-resident eigensystem handoff only if justified

## Objective

Determine whether avoiding host round trips between the CUDA eigensolver and CUDA Lehmann contraction yields enough benefit to justify a narrow device-residency mechanism.

This task is explicitly conditional.

## Measure first

Using ACC-10, quantify eigensolver D2H eigenpair cost, Lehmann H2D eigenpair cost, frequency of eigensystem reuse, and memory footprint.

## Decision A — not significant

If transfer overhead is small relative to solve+contraction:

- make no new residency abstraction;
- document the result;
- close the task.

## Decision B — significant

Introduce the narrowest backend-owned handle/state required for:

```text
CUDA eigensolver -> CUDA Lehmann
```

Do not expose device pointers throughout physics modules.

Requirements:

- explicit validity generation;
- explicit release;
- host materialization remains available;
- CPU backend remains unaffected;
- no global device state.

## Checklist

- [ ] eigensystem transfer overhead measured
- [ ] reuse frequency measured
- [ ] memory cost measured
- [ ] no code added if benefit insignificant
- [ ] if added, residency scope limited to real consumers
- [ ] generation/lifecycle explicit
- [ ] host materialization retained
- [ ] CPU backend unchanged
- [ ] eigensolver and Lehmann correctness retained
- [ ] end-to-end benefit demonstrated if code changed

**Commit message:** `Add CUDA eigensystem residency when beneficial`

---

# ACC-12 — Port H(k) assembly only if profiling justifies it

## Objective

Determine whether CPU H(R)->H(k) assembly has become a material bottleneck after eigensolution acceleration.

Do not assume that it has.

## Measurement

Use ACC-06/10 representative workloads.

Measure assembly time, eigensolve time, transfer time, and total workflow time for small and larger reciprocal matrices.

## Decision A — assembly insignificant

If assembly remains a small fraction of total runtime:

- do not port it;
- document the result;
- close the task.

## Decision B — assembly material

Only then prototype GPU assembly.

Preserve the exact existing phase convention and short-range operator semantics.

Prefer a narrow backend/private kernel.

Do not replace the canonical CPU assembler.

Validate assembled H(k) directly using elementwise/norm comparison, Hermiticity, and downstream eigenvalues.

Do not use final SCF energy as the first oracle.

## Checklist

- [ ] H(k) assembly fraction measured
- [ ] representative workloads measured
- [ ] port decision justified quantitatively
- [ ] no GPU assembler added if unnecessary
- [ ] if added, phase convention preserved
- [ ] H(k) direct comparison passes
- [ ] Hermiticity passes
- [ ] eigenvalues unchanged
- [ ] CPU assembler retained
- [ ] net end-to-end benefit demonstrated if code changed

**Commit message:** `Accelerate reciprocal Hk assembly when justified`

---

# ACC-13 — Establish RS KPM/transport GPU completeness and performance

## Objective

Use the Phase-II validated charge/spin/orbital Kubo-Bastin contracts to establish how much of the transport workload is already accelerated by the existing RS CUDA plugin and fill only genuine high-value gaps.

Do not rewrite KPM transport.

## Audit

Map CPU KPM stages, current GPU moment kernels, velocity/current operator transfer, spin-current path, orbital-current path, and remaining CPU contractions/postprocessing.

Identify the true runtime hotspots using ACC-00.

## Correctness

For existing GPU-supported stages compare:

- KPM moments;
- charge conductivity;
- spin conductivity;
- orbital conductivity;
- symmetry-forbidden components.

Use the same conventions established in Phase II.

## Performance

Measure convergence-representative workloads over polynomial order, system size, and number of stochastic vectors if applicable.

Separate kernel time, data transfer, and CPU postprocessing.

## New kernels

Add a new GPU kernel only if profiling identifies a material CPU bottleneck, the scientific CPU path is validated, and the kernel has a clear narrow contract.

## Checklist

- [ ] transport CPU/GPU dataflow mapped
- [ ] existing GPU moment support validated
- [ ] charge conductivity CPU/GPU compared
- [ ] spin conductivity CPU/GPU compared
- [ ] orbital conductivity CPU/GPU compared
- [ ] symmetry constraints retained
- [ ] performance hotspots measured
- [ ] no redundant kernel added
- [ ] new kernels added only with benchmark evidence
- [ ] end-to-end transport speedup recorded
- [ ] support matrix updated

**Commit message:** `Establish CUDA KPM transport coverage`

---

# ACC-14 — Accelerator release gate and support matrix

## Objective

Close Phase III-A by documenting exactly what accelerator functionality is scientifically correct, performance-beneficial, and supported.

This is primarily a validation/documentation task.

Do not broaden GPU scope merely to make the support matrix look larger.

## Support matrix

Document at minimum:

### RS CUDA
For each relevant route:
- supported;
- validated;
- performance-tested;
- known limitations.

### Reciprocal CUDA
Record:
- standard Hermitian;
- eigenvalues/eigenvectors;
- arbitrary-k;
- normal mesh;
- Si/sp;
- Fe/spd;
- second-order/HOH/CCOR where established;
- generalized overlap status;
- MPI rank-device status;
- H(k) assembly location;
- eigenstate residency status;
- Lehmann status.

### Deferred
State clearly that GBT- and TD-DFT-specific GPU claims are not part of this release gate.

## Performance report

Summarize CPU/GPU crossover, recommended tile/batch policy, representative end-to-end speedups, cases where GPU is slower and CPU should be preferred, and transfer/memory limitations.

Do not report only best-case speedups.

## Correctness gate

Run lean CPU unit/quick gates, relevant RS GPU correctness matrix, reciprocal CPU/GPU fixtures, Lehmann CPU/GPU contracts if implemented, and transport CPU/GPU contracts if implemented.

## Build matrix

Verify:
- CPU-only build;
- CUDA build;
- MPI+CUDA build where available.

## Documentation

Update developer/accelerator documentation with build requirements, backend selection, supported feature table, benchmark methodology, and known limitations.

## Checklist

- [ ] RS CUDA support matrix completed
- [ ] reciprocal CUDA support matrix completed
- [ ] unsupported generalized scope documented
- [ ] H(k) CPU/GPU decision documented
- [ ] eigensystem residency decision documented
- [ ] Lehmann status documented
- [ ] transport status documented
- [ ] CPU/GPU crossover table published
- [ ] representative speedups published honestly
- [ ] slow-GPU regimes documented
- [ ] CPU-only build passes
- [ ] CUDA build passes
- [ ] MPI+CUDA checked where hardware permits
- [ ] lean CPU gate remains green
- [ ] accelerator correctness gate passes
- [ ] GBT/TDDFT accelerator scope explicitly deferred

**Commit message:** `Document validated accelerator support`
