# RS-LMTO-ASA `fable_v3` optimization and MPI audit

**Audit date:** 2026-08-11  
**Audited revision:** `c3856036dee703f2ca3ee5ce3a6dc05f71801e8f` (`fable_v3`)  
**Scope:** k-space SCF, TD-DFT, accelerator architecture, and MPI correctness/coverage. The legacy atomic LMTO/potential machinery is deliberately kept on the CPU.

## Executive assessment

The code already has a sensible CPU/GPU boundary: the atomic and real-space LMTO construction can remain on the CPU, while the reciprocal-space pipeline begins at Fourier assembly of \(H(\mathbf{k})\) and continues through diagonalization, occupations/projections, and TD-DFT response. That complete reciprocal pipeline—not an isolated eigensolver call—is the right GPU target.

The main findings are:

1. **GPU-porting k-space is feasible, but it must be batched and tiled.** Typical matrices are small enough that one GPU call per k point will be dominated by launch and transfer overhead. Fourier assembly, eigensolution, and the immediate TD-DFT contractions should share a device-resident tile.
2. **Start with `reciprocal_mode='ham_only'` in double precision.** This is already the supported mode for the pair-potential TD-DFT backend. Generalized overlap (`zhegv`) needs a separate design and should retain the CPU path initially.
3. **TD-DFT has both GPU and MPI opportunities.** It currently distributes only q points. Every rank redundantly solves the complete k mesh, and a q owner solves every k+q point serially. This scales poorly when `nq < nranks` and duplicates substantial memory.
4. **MPI cannot yet be called verified for the usual workflows.** The test matrix covers a few real-space SCF cases but no MPI k-space SCF, reciprocal post-processing, exchange, or end-to-end TD-DFT material case. Static inspection also found concrete initialization, build, output-file, and GPU-device-selection risks.
5. **Stabilization should precede the GPU port.** Fix the MPI correctness items, restore a green example-test workflow, and establish end-to-end timing fixtures before changing numerical backends.

## 1. Architecture and hot-path map

### 1.1 K-space SCF

The k-space SCF path is orchestrated by [`self.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/self.f90#L925-L1058):

1. generate/reuse the reciprocal mesh;
2. build \(H(\mathbf{k})\), and optionally \(S(\mathbf{k})\);
3. diagonalize each k point;
4. compute occupations, DOS, projected DOS, spin moments, and band energy;
5. map reciprocal observables back into the CPU-side self-consistency machinery.

This provides a clean offload seam after the CPU has produced real-space hopping/overlap blocks.

The principal compute regions are:

- [`reciprocal_fourier.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/reciprocal_fourier.f90): Fourier assembly over k points. The second-order path also forms \(H O H\) using `zgemm`, which maps naturally to strided-batched GEMM.
- [`reciprocal_bands.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/reciprocal_bands.f90#L11-L161): OpenMP loop over k points with LAPACK `zheev` or `zhegv`. Each thread owns matrix/workspace copies.
- `calculate_eigenpairs_at_kpoints` and `diagonalize_single_kpoint` in `reciprocal_fourier.f90`: arbitrary k/k+q service used by TD-DFT. It is serial across unique points, allocates work arrays per point, and performs a LAPACK workspace query per solve.
- [`reciprocal_dos.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/reciprocal_dos.f90): occupations, total/projected DOS, spin density, and tetrahedron/Gaussian accumulation.

Current MPI k distribution is conditional: the reciprocal Hamiltonian is distributed for the Gaussian DOS path, while tetrahedron/Blöchl and several Green-function/moment consumers require an undistributed mesh. Tetrahedron accumulation itself can be divided between ranks, but the preceding eigensystems remain replicated. Data ownership therefore needs to be made explicit before extending k-parallelism.

### 1.2 TD-DFT

The TD-DFT driver in [`calculation.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/calculation.f90#L1966-L2510) currently does the following:

- partitions q points between MPI ranks;
- computes the complete k eigensystem independently on every rank;
- for each owned q, computes the complete k+q eigensystem;
- builds transition vertices and accumulates \(\chi^0\) and, when selected, pair-potential \(\Xi\);
- performs Dyson enhancement and writes q-owned outputs;
- for mode analysis, allocates full q-frequency response tensors on every rank and all-reduces them, although only rank 0 consumes the result.

The eigenpair response kernels in [`tddft_chi0.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/tddft_chi0.f90) and [`tddft_xi.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/tddft_xi.f90) already batch transitions for matrix accumulation, but vertices and denominators are still generated in CPU loops. This is a good basis for GPU tiling.

The Green backend in [`tddft_chi0_green.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/tddft_chi0_green.f90) is explicitly an eigenpair-resolvent reference, not a native real-space Green-function provider. It contains deeper k-energy-frequency loops and per-call matrix allocations, but it should be optimized only after production profiling shows that this backend is strategically important.

## 2. Correctness and scalability findings

| Priority | Finding | Consequence | Recommended action |
|---|---|---|---|
| P0 | `rank`, `numprocs`, and `ierr` are not initialized in [`mpi.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/mpi.f90#L21-L28). Non-MPI code still reads these globals. | Non-MPI behavior relies on non-standard implementation initialization; `get_mpi_variables` can read undefined `numprocs`. | Initialize to `rank=0`, `numprocs=1`, `ierr=0`; add a serial unit test of range/mapping setup. |
| P0 | `report.out`, `minfo.out`, and `linfo.out` are opened with `status='replace'` by every rank before the root-only write section in [`self.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/self.f90#L1375-L1409). | Shared-file races and possible late truncation, especially on distributed filesystems. | Make rank 0 own open/write/close; keep collective calculations outside that guard. |
| P0 | `totaldos.out` and `magneticdos.out` are written by all ranks in the k-space SCF helper. | Distributed Gaussian runs race on shared output files. | Root-only output after reductions; test 1/2/4-rank file equivalence. |
| P0 | The CUDA recursion plugin always passes device 0 in [`rsrec_cuda_plugin.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/rsrec_cuda_plugin.f90#L217-L241). | Multiple local MPI ranks contend for GPU 0 unless the scheduler masks visibility per rank. | Add explicit local-rank-to-device mapping, device-count validation, user override, and mapping log. Prefer `MPI_Comm_split_type(...MPI_COMM_TYPE_SHARED...)` over launcher-specific environment variables. |
| P1 | `nomp` exists only when OpenMP is found, but the MPI banner uses it unconditionally in [`main.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/c3856036dee703f2ca3ee5ce3a6dc05f71801e8f/source/main.f90#L20-L52). | MPI-without-OpenMP configuration can fail to compile. | Declare/initialize `nomp=1` independently of OpenMP and overwrite it inside the OpenMP region. |
| P1 | CMake changes `CMAKE_Fortran_COMPILER` after the Fortran language is already enabled. | MPI configuration is brittle and toolchain dependent. | Link targets to `MPI::MPI_Fortran`; use target compile definitions. If a wrapper compiler is required, select it before `project()`. |
| P1 | The latest public downstream example-test run failed on Linux and macOS. | The current branch does not have a clean public validation signal. | Retrieve authenticated logs/artifacts, identify the failing cases, and restore the test gate before performance work. |
| P1 | TD-DFT distributes q only and replicates the full k solve, static work, and response memory on every rank. | Weak scaling for small q sets; idle ranks when `nq < nranks`; high memory and communication. | Introduce q groups with k workers; reduce/gather owned q slabs to root instead of all-reducing full tensors. |
| P1 | No MPI end-to-end coverage exists for k-space SCF or TD-DFT. | Collective/order, file ownership, and numerical equivalence regressions can pass CI unnoticed. | Add serial-vs-MPI paired tests and multi-rank output checks before refactoring ownership. |
| P2 | Arbitrary-k diagonalization allocates and queries LAPACK workspace for every point; duplicate folding is quadratic in the number of requested points. | Avoidable CPU overhead; awkward GPU batching boundary. | Introduce a reusable workspace/tile object and hash/index folded k points. |
| P2 | Mode output uses full-tensor `MPI_ALLREDUCE`, with full arrays allocated on all ranks. | Memory is \(O(n_r^2 n_\omega n_q)\) per rank and communication copies the entire result. | Use `MPI_REDUCE`/`Gatherv` to root or stream q slabs into the mode analyzer. |

No definite collective mismatch was found in the static scan, but that is not equivalent to runtime proof. The file races and uninitialized serial state are actionable independently of runtime reproduction.

## 3. GPU blueprint

### 3.1 Backend boundary

Retain these operations on the CPU:

- atomic potential and radial LMTO operations;
- construction of real-space Hamiltonian/overlap and legacy tangent data;
- SCF mixing and final text I/O;
- initially, generalized-overlap eigensolutions.

Move these operations behind a reciprocal accelerator API:

- upload immutable real-space block topology/values;
- generate Bloch phases and assemble tiled batches of \(H(\mathbf{k})\);
- second-order \(H O H\) with batched GEMM;
- batched standard Hermitian eigensolution;
- occupations, projections, and spin moments where this avoids a device-to-host round trip;
- TD-DFT transition vertices, denominators, and response accumulation.

The API should be backend-neutral even if CUDA is implemented first. Share device discovery, stream/error handling, and rank mapping with the existing recursion plugin, but use a separate reciprocal context: its lifetime and dimension fingerprint differ from the recursion context.

### 3.2 Numerical library choice

The existing B4 planning note should be corrected before implementation. Current cuSOLVER documentation provides single-matrix `cusolverDnZheevd` and batched Jacobi `cusolverDnZheevjBatched`; it does **not** provide `cusolverDnZheevdBatched`. The first prototype should benchmark `cusolverDnZheevjBatched` against a tile of independent `Zheevd` calls/streams for the actual 18/36/72-ish matrix sizes. See the official [cuSOLVER documentation](https://docs.nvidia.com/cuda/cusolver/index.html).

The current CUDA CMake target links cuBLAS but not cuSOLVER. Add `cusolver` discovery and `CUDA::cusolver` only when the reciprocal accelerator is enabled. Batched LU/solve routines in [cuBLAS](https://docs.nvidia.com/cuda/cublas/index.html) can later support Dyson/Green-function work.

All SCF and TD-DFT acceptance work should begin in complex FP64. Reduced precision is a separate, explicitly validated experiment; it should not inherit the recursion plugin's FP32 Chebyshev default.

### 3.3 Tiled execution model

A production tile should flow as:

1. CPU produces or updates the real-space LMTO blocks.
2. Upload the block data once per SCF iteration (or only when changed).
3. GPU assembles a tile of k Hamiltonians.
4. GPU diagonalizes the tile.
5. SCF path accumulates occupations/projections and returns only required observables; TD-DFT retains eigenvectors/eigenvalues for the response tile.
6. Double-buffer CPU-generated pair-potential operator tiles with GPU response contraction if legacy atomic evaluation remains on CPU.

Tile size should be selected from available device memory and benchmarked. Do not allocate complete host and device copies of every \(H(\mathbf{k})\) by default.

### 3.4 Generalized overlap

Phase 1 should support only `ham_only` on GPU and route `generalized_overlap_proxy` to the existing CPU implementation. A later path can transform the generalized Hermitian problem with batched Cholesky/triangular solves, but it needs independent conditioning, residual, and orthogonality validation. This restriction also aligns with the current pair-potential TD-DFT capability gate.

## 4. TD-DFT parallel blueprint

Use a two-level decomposition rather than q-only distribution:

- **Across q groups:** assign one or more q points to a rank group/GPU. The base k eigensystem is computed once per group and reused.
- **Within a q group:** distribute k tiles across workers for k+q eigenpairs and response accumulation, then reduce \(\chi^0\)/\(\Xi\) inside the group.
- **On a GPU:** batch k points, band transitions, and frequency tiles. Fuse occupation-difference/denominator weighting where practical and use GEMM for response accumulation.

This allows two useful operating modes:

- many q points: one q stream per GPU, maximizing k-eigensystem reuse;
- one/few q points: several ranks/GPUs cooperate over k, avoiding idle ranks.

Pair-potential operators can be large (`matrix × matrix × site × k`). Generate and consume them as tiles. The legacy atomic/tangent calculation remains CPU-side, while matrix contraction moves to the GPU. Pinned double buffers are appropriate only after profiling shows overlap is effective.

The current full-tensor mode-analysis arrays should be replaced by root-owned gathered q slabs or a streaming analysis interface. This is an MPI memory fix regardless of GPU work.

## 5. MPI validation matrix

Existing MPI coverage is narrow: six SCF cases carry `mpi_procs`, four of those explicitly register serial/MPI paired launches, and only one MPI k-space occupation unit test is present. The post-processing case matrix has no MPI process counts. There is no end-to-end TD-DFT material execution in the normal CTest matrix.

The minimum validation expansion should be:

| Workflow | Required configurations | Comparison |
|---|---|---|
| Real-space SCF | 1, 2, 4 ranks; OpenMP on/off | energy, charge, moments, recursion outputs, single-writer files |
| K-space SCF Gaussian | 1, 2, 4 ranks | Fermi level, band energy, DOS integrals/projected moments, output identity |
| K-space SCF tetra/Blöchl | 1, 2, 4 ranks | same, plus tetra ownership/reduction checks |
| Band path / reciprocal DOS | 1 and multi-rank | ordered k path, eigenvalues, root-owned files |
| Exchange/frozen magnon | 1 and multi-rank | exchange tensors/dispersion and deterministic output |
| Green/BSF/conductivity | supported rank counts | complex spectra and sum rules |
| TD-DFT eigenpair backend | 1, 2, 4 ranks; `nq<nranks` case | \(\chi^0\), \(\Xi\), loss, modes, Goldstone diagnostics, file manifest |
| CPU/GPU reciprocal | 1 rank and one rank/GPU | eigen residuals, orthogonality, SCF observables, TD-DFT spectra |
| Multi-GPU MPI | 1 rank/GPU and oversubscription rejection | device mapping, numerical equivalence, no shared-file races |

Use rank-count invariance tests, not merely “program exited successfully.” For output files, check both numerical content and ownership (one complete file, no per-rank truncation). Include sanitizer/debug builds where the compiler supports bounds checking, uninitialized-variable traps, and floating-point exception checks.

## 6. Recommended delivery order

### Gate 0 — correctness baseline

1. Restore green example-test CI and capture the exact failing cases.
2. Fix MPI global initialization, MPI-without-OpenMP compilation, and root-only file ownership.
3. Replace post-`project()` MPI compiler mutation with target-based MPI linkage.
4. Add k-space SCF and TD-DFT serial/MPI equivalence cases.
5. Add deterministic local-rank GPU device selection to the existing recursion plugin.

### Phase A — measurements and interfaces

1. Add timers for Fourier assembly, eigensolve, projection/DOS, k+q solve, pair-operator construction, vertex construction, denominator evaluation, accumulation, transfers, and MPI reductions.
2. Establish production-representative fixtures (bcc Fe/fcc Ni plus multi-atom cells), not only synthetic 16/32-k-point TD-DFT fixtures.
3. Record matrix size, k/q/frequency counts, memory high-water mark, and CPU thread/rank/GPU topology in performance metadata.
4. Introduce reusable k-point workspaces/tiles on CPU first; this reduces allocation overhead and defines the GPU API cleanly.

### Phase B — k-space GPU minimum viable path

1. CUDA reciprocal context with FP64 batched Fourier assembly and `ham_only` eigensolution.
2. CPU fallback and runtime capability query.
3. Gaussian k-space SCF integration, because its k ownership is already distributed.
4. Device-side projections/occupations where profiling shows transfers dominate.
5. Benchmark gate: require end-to-end benefit at production matrix/k sizes, including transfers, before making GPU the recommended backend.

### Phase C — TD-DFT accelerator and hybrid MPI

1. Reuse the base k eigensystem per q group.
2. GPU k+q assembly/eigensolution.
3. Batched vertex/denominator/\(\chi^0\) accumulation, followed by pair-potential \(\Xi\).
4. Two-level q-group/k-worker communicators and root-owned q-slab collection.
5. Validate spectra, sum rules, causality/spectral-weight diagnostics, and Goldstone behavior against FP64 CPU reference.

### Phase D — broaden reciprocal consumers

1. Make tetra/Blöchl data ownership compatible with distributed eigensystems.
2. Port or adapt band paths, projected moments, and reciprocal Green/BSF consumers.
3. Consider batched LU/solve for Dyson and Green paths only after profiling.
4. Evaluate a HIP backend once the API and CUDA path are stable, unless an immediate target system requires AMD support.

## 7. Verification performed for this audit

- Repository inspected at the pinned revision above; the worktree was not modified.
- Standalone Python TD-DFT dispatch test: **PASS**.
- TD-DFT validation campaign script: **PASS**.
- Native Fortran, MPI, and CUDA execution: **not run locally**, because this audit container has no Fortran compiler, MPI wrapper/runtime, or CMake executable.
- Public GitHub Actions evidence: the [branch build run](https://github.com/rslmtoasa/rslmtoasa/actions/runs/31479506589) completed successfully, while the latest downstream [example-test run](https://github.com/rslmtoasa/rslmtoasa/actions/runs/31482023846) failed on both Linux and macOS. Public annotations expose only an exit code, so authenticated logs/artifacts are needed to identify the failing case.

## Decision summary

The recommended first implementation target is a **CUDA, FP64, `ham_only`, tiled k-space pipeline** covering Fourier assembly + eigensolution, built behind a backend-neutral reciprocal API and validated first in Gaussian k-space SCF. In parallel, make MPI behavior trustworthy through the Gate 0 fixes and tests. Once the k-space tile remains device-resident, extend the same tile into TD-DFT k+q eigenpairs and response accumulation, with q groups plus k-level cooperation. This sequence keeps the legacy atomic LMTO code untouched while attacking the dominant reusable reciprocal-space costs.
