# RS-LMTO-ASA CPU-first reciprocal refactoring package

## Purpose

This package contains six sequential, copy-paste-ready implementation prompts for the `fable_v3` branch of RS-LMTO-ASA:

- **RF-01:** characterization, numerical oracles, and performance baseline;
- **RF-02:** MPI runtime state, build integration, rank-local device identity, and root-owned output;
- **RF-03:** explicit k-point worksets and ownership;
- **RF-04:** batched reciprocal assembly and reusable workspaces;
- **RF-05:** abstract reciprocal execution backend with a complete CPU/LAPACK implementation;
- **RF-06:** shared TD-DFT transition engine for \(\chi^0\) and pair-potential \(\Xi\).

All six tasks can be developed without a GPU. RF-02 and the distributed portions of RF-03 require an MPI compiler/runtime for full validation, but not accelerator hardware. CUDA/HIP implementation begins only after RF-05. TD-DFT GPU implementation begins only after RF-06.

The tasks are deliberately narrower than a general code modernization. They preserve the legacy atomic LMTO machinery and prepare one coherent reciprocal execution layer for later CPU, MPI, CUDA, and possibly HIP backends.

## Dependency order

```text
RF-01 → RF-02 → RF-03 → RF-04 → RF-05 → RF-06 → GPU implementation
```

Do not combine two RF tasks in one commit. Each stage must be green before starting the next.

---

# Rules shared by every RF task

Copy the following rules into the working context for every coding task.

## Repository discipline

1. Work from the latest `fable_v3` revision. Record `git rev-parse HEAD` before editing; do not reset to the revision used by the original audit.
2. Read, in order:
   - `CLAUDE.md`;
   - `docs/DEVELOPER_MAP.md`;
   - `docs/dev/REFACTORING_PLAN.md`;
   - `REFACTORING_PHASE2.md`;
   - `tests/README.md`;
   - `tests/KNOWN_ISSUES.md`;
   - the files named in the task-specific prompt.
3. Inspect `git status --short` before editing. Preserve all user changes and unrelated work.
4. One RF task equals one commit. Do not fix adjacent issues unless the prompt explicitly includes them.
5. If a new issue is discovered, record verified evidence in `tests/KNOWN_ISSUES.md`, mention it in the handoff, and continue only if it does not block the current task.
6. Never regenerate references simply because a test changed. Diagnose every difference first.
7. Keep `ENABLE_CUDA_PLUGIN=OFF` throughout RF-01–RF-06 validation. Do not add CUDA/HIP source, mock accelerator kernels, or hardware-dependent tests.

## Non-negotiable physics and scope constraints

1. Preserve the default FP64 CPU numerical path and all current physics conventions.
2. Preserve retarded TD-DFT denominators, occupation conventions, left/right circular-channel order, band-window semantics, and normalization by the explicit k-weight sum.
3. Preserve first-order, second-order/HOH, SOC, CCOR, GBT, and generalized-overlap behavior unless a task explicitly narrows its acceptance surface.
4. Do not refactor the atomic LMTO part.
5. `source/self.f90` and `source/symbolic_atom.f90` are not refactoring targets. At most, make a very small orchestration-level wiring edit when unavoidable; do not reorganize their bodies.
6. Do not change namelist defaults or introduce a new user-facing option unless explicitly requested.
7. Do not change precision, introduce approximate occupation pruning, or relax scientific tolerances for performance.
8. Retain a CPU numerical oracle whenever reduction ordering or batching changes.

## Required RS-LMTO-ASA implementation style

New code must look like it belongs beside `recursion.f90`, `hamiltonian.f90`, and `reciprocal.f90`:

1. Use a private-by-default module with explicitly public derived types and APIs.
2. Use concrete derived types for stateful concepts such as worksets and workspaces.
3. Use abstract types only at real substitution boundaries:
   - reciprocal execution backend in RF-05;
   - TD-DFT vertex provider in RF-06.
   Do not create an inheritance hierarchy for simple data containers.
4. Give stateful public types:
   - a generic constructor interface named after the type where construction needs arguments;
   - type-bound procedures;
   - `restore_to_default` when the object has meaningful reusable defaults;
   - a final destructor when the type owns allocatable or backend storage.
5. Constructors should associate pointers to longer-lived upstream objects declared with `target`, call `restore_to_default`, and then initialize owned state.
6. Owned arrays should be `allocatable`. Use pointers only for non-owning associations to objects whose lifetime is controlled elsewhere.
7. Put declarations and module-procedure interfaces in the parent module. Put substantial implementations in one or more submodules, following the `hamiltonian_mod`/`hamiltonian_build` and `reciprocal_mod`/`reciprocal_*` pattern.
8. Provide Doxygen-style `!>` documentation for every public type and procedure, including ownership, layouts, units, side effects, and distributed/local semantics.
9. Use `precision_mod, only: rp`; do not introduce literal-kind inconsistencies.
10. Use `g_logger` with `__FILE__` and `__LINE__` at application boundaries. Avoid long internal validation ladders; validate at constructors, public APIs, I/O boundaries, MPI/plugin boundaries, and LAPACK failures.
11. Use existing `g_timer` conventions where safe. Do not call a global timer concurrently from OpenMP workers unless its thread safety has been established. For parallel kernels, time outside the parallel region or reduce thread-local counters.
12. Follow current indentation, naming, continuation, preprocessor, and CMake source-list conventions.
13. Avoid allocation, deallocation, logging, and polymorphic dispatch inside the innermost k/band/transition loops.

## Live progress document

RF-01 should create `docs/dev/plans/RF_CPU_RECIPROCAL_REFACTOR.md` if it does not already exist. It should contain:

- the six-stage overview;
- the current baseline commit and toolchain;
- the checklist for every RF stage;
- validation commands and results;
- decisions that affect later GPU/MPI work;
- the commit SHA for each completed stage.

Every later RF task must update only its own section. Tick completed items from `[ ]` to `[x]`; leave blocked or unverified items unticked and explain them. Include the progress-document update in the same task commit.

## Minimum build matrices

Use clean, separate build directories. Adapt generator and compiler flags to the machine, but record the exact commands.

### Serial/OpenMP release build

```bash
cmake -S . -B build-rf-serial \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA_PLUGIN=OFF \
  -DENABLE_MPI=OFF \
  -DENABLE_OPENMP=ON \
  -DRUN_UNIT_TESTS=ON \
  -DRUN_REG_TESTS=ON \
  -DRUN_EXAMPLE_TESTS=ON
cmake --build build-rf-serial --parallel
ctest --test-dir build-rf-serial --output-on-failure
```

### MPI/OpenMP release build

Required from RF-02 onward when MPI is available:

```bash
cmake -S . -B build-rf-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA_PLUGIN=OFF \
  -DENABLE_MPI=ON \
  -DENABLE_OPENMP=ON \
  -DRUN_UNIT_TESTS=ON \
  -DRUN_REG_TESTS=ON \
  -DRUN_EXAMPLE_TESTS=ON
cmake --build build-rf-mpi --parallel
ctest --test-dir build-rf-mpi --output-on-failure
```

Also configure MPI with `ENABLE_OPENMP=OFF` in RF-02 to protect that supported build combination.

### Debug build

Use compiler-appropriate bounds, backtrace, floating-point exception, and uninitialized-value diagnostics. For GNU Fortran, a representative configuration is:

```bash
cmake -S . -B build-rf-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_CUDA_PLUGIN=OFF \
  -DENABLE_MPI=OFF \
  -DENABLE_OPENMP=OFF \
  -DRUN_UNIT_TESTS=ON \
  -DRUN_REG_TESTS=ON \
  -DRUN_EXAMPLE_TESTS=ON \
  -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan"
cmake --build build-rf-debug --parallel
ctest --test-dir build-rf-debug --output-on-failure
```

Do not paste GNU-only flags into non-GNU builds.

---

# RF-01 prompt — characterization and CPU baseline

## Copy-paste implementation prompt

You are implementing **RF-01: characterization, numerical oracles, and CPU performance baseline** in the RS-LMTO-ASA `fable_v3` branch.

This task must not introduce a new execution backend and must not reorganize production algorithms. Its purpose is to make the current reciprocal and TD-DFT behavior measurable before later ownership, batching, and abstract-interface changes.

Follow all shared rules from the RF package. In particular, read the repository guidance first, preserve bit-level/default-path behavior, keep CUDA disabled, and create or initialize `docs/dev/plans/RF_CPU_RECIPROCAL_REFACTOR.md` as the live checklist.

### Inspect before editing

Read at minimum:

- `source/reciprocal.f90`;
- `source/reciprocal_lifecycle.f90`;
- `source/reciprocal_fourier.f90`;
- `source/reciprocal_bands.f90`;
- `source/reciprocal_dos.f90`;
- `source/tddft_chi0.f90`;
- `source/tddft_xi.f90`;
- `source/calculation.f90`, limited to reciprocal SCF and TD-DFT orchestration;
- `tests/unit/test_arbitrary_k_eigenpairs.f90`;
- `tests/unit/test_kspace_occupations.f90`;
- `tests/unit/test_tddft_chi_ks.f90`;
- `tests/unit/test_tddft_direct_xi.f90`;
- `tests/unit/test_tddft_cpu_profile.f90`;
- relevant k-space SCF and post-processing cases and references.

Do not duplicate oracles that already exist. Extend them where coverage is missing.

### Required work

1. **Capture a trustworthy baseline.**
   - Record HEAD, compiler/version, BLAS/LAPACK, build type, OpenMP settings, CPU model, and test commands.
   - Run the existing serial unit, regression, SCF, and post-processing suites before editing.
   - If the starting revision is not green, identify the failing test names and stop numerical refactoring. Record the failure without changing references.

2. **Strengthen arbitrary-k and reciprocal oracles.**
   Extend the existing arbitrary-k test rather than replacing it. Cover, using small deterministic models:
   - `k`, folded `k`, and `k+G` equivalence;
   - repeated/permuted duplicate points and preservation of caller order;
   - an off-mesh point compared with direct Fourier assembly;
   - no mutation of the canonical reciprocal mesh by caller-owned arbitrary-k requests;
   - first-order and second-order/HOH assembly;
   - SOC where already supported by the fixture;
   - Hermiticity of assembled matrices;
   - eigenvalue agreement and gauge-invariant eigenspace comparison;
   - residual norm \(\|HU-U\Lambda\|\) and orthonormality \(\|U^\dagger U-I\|\);
   - generalized-overlap residual/metric orthogonality if a compact deterministic fixture can use the current `zhegv` route without expanding scope.

3. **Characterize SCF observables.**
   Ensure the k-space SCF bcc Fe case pins the observables that later GPU work must preserve:
   - Fermi level;
   - canonical electron count;
   - band energy;
   - charge and site spin moments;
   - total DOS normalization and selected DOS samples;
   - projected moment/DOS quantities already considered stable.
   Prefer extending comparison keys over creating a new physical input. Do not bless unexplained changes.

4. **Characterize TD-DFT equivalence.**
   Preserve the scalar path as numerical oracle and ensure tests compare scalar and batched paths for:
   - transverse circular channels with independently constructed left/right vertices;
   - finite temperature and zero-temperature limits;
   - exact pruning disabled;
   - a positive explicit pruning tolerance as an approximate opt-in;
   - restricted and full band windows;
   - static divided-difference behavior near degeneracy;
   - ordinary \(\chi^0\), static \(\chi^0\), constant pair-potential \(\Xi\), and k-dependent pair-potential \(\Xi\);
   - trace/site spectral products and metadata fields.

5. **Create a non-gating CPU profile.**
   Extend or replace `test_tddft_cpu_profile.f90` only as needed so that it reports separately:
   - reciprocal Fourier assembly;
   - k eigensolution;
   - arbitrary k+q assembly/eigensolution;
   - TD-DFT vertex construction;
   - denominator generation;
   - response accumulation;
   - pair-potential operator construction when available in a realistic driver;
   - peak host memory or analytically computed principal array sizes.
   Include at least representative matrix dimensions and k/frequency counts relevant to one-site and multi-site systems. Timings must be informational, not pass/fail CI thresholds.

6. **Document numerical contracts.**
   In the live progress document, state which comparisons are expected to be bitwise, which are tolerance-based because reduction order differs, and the actual tolerances with physical/numerical justification.

### Explicit non-goals

- No k-point workset type yet.
- No reciprocal backend interface yet.
- No MPI ownership redesign.
- No algorithm replacement.
- No GPU terminology in public input.
- No reference regeneration without a diagnosed physical reason.

### RF-01 checklist

- [ ] Repository guidance and relevant source/tests read.
- [ ] Clean starting status and baseline HEAD recorded.
- [ ] Serial release baseline run and results recorded.
- [ ] Debug baseline run and results recorded.
- [ ] Existing baseline failures, if any, isolated without reference churn.
- [ ] Arbitrary-k folding, ordering, and side-effect oracles complete.
- [ ] Hermiticity, residual, and orthogonality checks complete.
- [ ] K-space SCF observable contract documented and tested.
- [ ] TD-DFT scalar/batched equivalence surface complete.
- [ ] CPU profile reports the named phases without CI timing thresholds.
- [ ] Bitwise versus tolerance-based expectations documented.
- [ ] Full post-change serial release and debug suites pass.
- [ ] `git diff --check` passes and unrelated files are untouched.
- [ ] RF-01 section of the live progress document is fully updated.
- [ ] One RF-01 commit created.

### Required handoff

Report:

- files changed;
- baseline and final test commands/results;
- measured phase timings and array sizes, clearly labelled machine-specific;
- remaining characterization gaps;
- the completed checklist with `[x]`/`[ ]` state;
- the commit SHA.

### Exact one-line commit message

```text
test(RF-01): establish reciprocal CPU refactoring baselines
```

---

# RF-02 prompt — parallel runtime and output ownership

## Copy-paste implementation prompt

You are implementing **RF-02: MPI runtime state, build integration, local-rank identity, and root-owned output** in the RS-LMTO-ASA `fable_v3` branch.

RF-01 must already be complete and green. Follow all shared rules. This is a CPU/MPI task; keep CUDA disabled. The aim is to remove undefined serial state and shared-file races while introducing a small class-like parallel context that later reciprocal and GPU code can consume.

### Inspect before editing

Read at minimum:

- `source/mpi.f90`;
- `source/main.f90`;
- MPI branches in `source/self.f90`, without refactoring its body;
- `source/reciprocal_fourier.f90` and `source/reciprocal_dos.f90`;
- TD-DFT MPI sections in `source/calculation.f90`;
- root-only writers throughout `source/`;
- MPI configuration in `CMakeLists.txt` and `source/CMakeLists.txt`;
- the serial/MPI registration logic and `UnitKspaceOccupations_mpi` test.

### Required design

Introduce a concrete public type, preferably in the existing MPI module unless a small dedicated module is cleaner:

```fortran
type, public :: parallel_context
   integer :: rank
   integer :: size
   integer :: local_rank
   integer :: local_size
   logical :: mpi_enabled
contains
   procedure :: is_root
   procedure :: owns_work
   procedure :: split_range
   procedure :: barrier
   procedure :: restore_to_default
   final :: destructor
end type parallel_context
```

Use a generic constructor consistent with project style. Under MPI it should derive world rank/size and node-local rank/size. Prefer MPI-3 `MPI_Comm_split_type(..., MPI_COMM_TYPE_SHARED, ...)` for local identity. Under a non-MPI build it must deterministically produce rank 0, size 1, local rank 0, local size 1.

Do not put `MPI_Finalize` in a finalizer. MPI lifetime remains explicitly owned by the program. The context finalizer may free only communicators it created, and only while MPI is active.

For compatibility, existing module globals may temporarily remain, but initialize them at declaration (`ierr=0`, `rank=0`, `numprocs=1`) and synchronize them from the context in one place. Do not permit independent mutable sources of truth.

### Required work

1. **Fix deterministic serial state.**
   - Initialize MPI compatibility globals.
   - Test serial `get_mpi_range`/mapping, including zero items and more workers than items through an explicit-size test helper if needed.

2. **Fix MPI-without-OpenMP.**
   - Make `nomp=1` valid and declared regardless of OpenMP discovery.
   - Verify `ENABLE_MPI=ON, ENABLE_OPENMP=OFF` builds and starts correctly.

3. **Modernize CMake MPI linkage.**
   - Do not assign `CMAKE_Fortran_COMPILER` after `project()`/language enablement.
   - Link the main library/binary and unit tests to `MPI::MPI_Fortran`.
   - Use target compile definitions for `USE_MPI` where feasible instead of appending global flags.
   - Preserve supported compiler-wrapper workflows selected before configuration.

4. **Centralize local-rank identity.**
   - Store node-local rank/size even without a GPU.
   - Add a pure or side-effect-free helper that maps local rank onto an available-device count supplied as an argument. It must reject invalid counts and support an explicit user override supplied programmatically; do not add a namelist option yet.
   - Unit-test mappings such as 4 local ranks/4 devices, 4 ranks/2 devices, and an invalid override. Do not call CUDA.

5. **Fix output ownership included in the audit.**
   - Make `report.out`, `minfo.out`, and `linfo.out` root-owned across open/write/close.
   - Make k-space SCF `totaldos.out` and `magneticdos.out` root-owned after collective reductions.
   - Inspect nearby reciprocal/TD-DFT writers and add ownership tests where the task touches them, but do not widen into unrelated I/O cleanup.
   - Keep collective calculations outside root-only branches.

6. **Improve rank-count tests.**
   Add serial-vs-MPI paired coverage for the existing k-space SCF bcc Fe case, initially with Gaussian DOS because its k-space distribution is already implemented. Compare scientific values and assert that exactly one complete shared output exists. Exercise 1, 2, and—where CI resources permit—4 ranks.

### Explicit non-goals

- No q-group/k-worker communicator hierarchy yet.
- No replacement of every global `rank` use in the repository.
- No GPU runtime calls.
- No changes to physical algorithms or reduction formulas.
- No full TD-DFT MPI redesign.

### RF-02 checklist

- [ ] RF-01 commit present and baseline green.
- [ ] Parallel context follows constructor/type-bound/finalizer style.
- [ ] Serial defaults are deterministic.
- [ ] Local rank/size derived with MPI shared-memory communicator.
- [ ] Finalizer does not own `MPI_Finalize`.
- [ ] Compatibility globals have one synchronization point.
- [ ] MPI-without-OpenMP builds and reports one thread.
- [ ] CMake uses `MPI::MPI_Fortran` without late compiler mutation.
- [ ] Device-index mapping helper is hardware-independent and tested.
- [ ] Report and k-space DOS files are root-owned across open/write/close.
- [ ] Gaussian k-space SCF has serial/MPI paired coverage.
- [ ] 1/2/4-rank results agree at documented tolerances where available.
- [ ] Oversubscription/zero-work ranks complete collectives safely.
- [ ] Serial release, MPI release, MPI-no-OpenMP, and debug tests pass.
- [ ] `git diff --check` passes and no unrelated cleanup was included.
- [ ] RF-02 progress section updated and one commit created.

### Required handoff

Report the rank/device mapping policy, communicator ownership, files made root-owned, build matrices, serial/MPI numerical comparisons, unticked items, and commit SHA.

### Exact one-line commit message

```text
refactor(RF-02): centralize parallel state and output ownership
```

---

# RF-03 prompt — k-point worksets and ownership

## Copy-paste implementation prompt

You are implementing **RF-03: explicit k-point worksets and ownership** in the RS-LMTO-ASA `fable_v3` branch.

RF-01 and RF-02 must be complete and green. Follow all shared rules and keep CUDA disabled. The purpose is to replace implicit local/full-mesh assumptions with one concrete, testable object without yet changing reciprocal mathematics.

### Inspect before editing

Read all locations using:

- `k_points`, `k_weights`, `nk_total`, `nk_local`, `k_start`, `k_end`;
- `k_l2g_map`, `k_g2l_map`, and `k_mesh_distributed_active`;
- `get_mpi_range` for k or q partitions;
- full-mesh guards in reciprocal moments/Green/DOS code;
- arbitrary-k folding and duplicate handling;
- TD-DFT construction of `kq_points`.

### Required design

Introduce a concrete class-like type in a focused parent module and submodule, for example:

```fortran
type, public :: kpoint_workset
   integer :: nk_global
   integer :: nk_local
   integer :: global_start
   integer :: global_end
   logical :: distributed
   real(rp), allocatable :: points(:, :)
   real(rp), allocatable :: weights(:)
   integer, allocatable :: local_to_global(:)
   integer, allocatable :: global_to_local(:)
contains
   procedure :: restore_to_default
   procedure :: validate
   procedure :: local_index
   procedure :: global_index
   procedure :: weight_sum
   procedure :: fold
   procedure :: shifted
   procedure :: select_tile
   final :: destructor
end type kpoint_workset
```

Adapt exact names to the code, but document unambiguously whether `points` and `weights` are local or global. The preferred invariant is that a workset stores the points actually owned/processed by the current execution context; a replicated workset has `nk_local == nk_global`.

Do not add polymorphism here. This is a data/ownership type, not a backend.

### Required work

1. **Constructor and invariants.**
   - Provide construction from a complete point/weight list plus `parallel_context` and a `distributed` policy.
   - Provide construction of caller-owned arbitrary-point worksets that are replicated by default.
   - Support zero-local-work ranks with valid zero-length allocatables and `global_start=1`, `global_end=0` or another single documented empty-range invariant.
   - Validate shape and mapping consistency only at construction/public boundaries.

2. **Centralize folding and shifting.**
   - Move or delegate reciprocal folding through the workset API.
   - Implement q shifting without mutating the base workset.
   - Preserve coordinate convention exactly.
   - Preserve original caller order after folded duplicate reuse.

3. **Migrate reciprocal mesh ownership.**
   - Make `reciprocal` own or associate one authoritative workset.
   - Migrate Fourier assembly, diagonalization, occupations, and Gaussian DOS to consume it.
   - Do not maintain two independently mutable copies of points, weights, and mappings.
   - If transitional compatibility fields must remain for one task, populate them from the workset in one direction, document them as transitional/read-only, and add consistency assertions in debug/test code.

4. **Represent replicated requirements explicitly.**
   - Replace hidden assumptions such as “tetrahedron implies full mesh” with an explicit `distributed`/replicated requirement checked at the public consumer boundary.
   - Do not yet redesign tetrahedron, reciprocal Green, or moment algorithms for distributed input. They may request/require a replicated workset with a clear diagnostic.

5. **Use the workset for TD-DFT k and k+q preparation.**
   - Replace raw ad hoc q-shift array construction where possible.
   - Do not change q-level MPI distribution or TD-DFT response algorithms yet.

6. **Tests.**
   Add focused unit tests for:
   - replicated and distributed construction;
   - 1/2/4-rank mapping;
   - zero-work ranks;
   - local/global round trips;
   - weight preservation and local/global sums;
   - folded boundary points, negative coordinates, and `k+G`;
   - q shifts without base mutation;
   - duplicate/permuted points with restored output order;
   - explicit rejection of a distributed workset by a replicated-only consumer;
   - unchanged k-space SCF and arbitrary-k results.

### Efficiency requirement

Replace the current quadratic exact duplicate scan in arbitrary-k processing with a deterministic indexed strategy if this can be done without changing folding semantics. For example, construct a canonical integer key from the folded coordinate using a documented tolerance consistent with current reciprocal coordinates. Test boundary behavior carefully. If a robust key cannot be introduced without changing semantics, retain the current search in RF-03 and document it for a later isolated optimization; do not smuggle in a risky tolerance change.

### Explicit non-goals

- No batched Hamiltonian storage yet.
- No LAPACK/backend abstraction.
- No distributed tetrahedron rewrite.
- No q-group communicators.
- No GPU layout decisions beyond documenting contiguous tile-friendly arrays.

### RF-03 checklist

- [ ] RF-01/RF-02 commits present and green.
- [ ] All legacy k ownership fields/consumers inventoried.
- [ ] `kpoint_workset` is concrete, class-like, documented, and finalized.
- [ ] Local versus global storage semantics are unambiguous.
- [ ] Empty-rank invariant is defined and tested.
- [ ] Folding and q shifting are centralized without coordinate changes.
- [ ] Caller ordering survives duplicate folding/reuse.
- [ ] Reciprocal has one authoritative point/weight/mapping source.
- [ ] Any transitional fields are one-way and consistency-checked.
- [ ] Replicated-only consumers fail clearly at their boundary.
- [ ] Serial and MPI mapping tests pass at available rank counts.
- [ ] K-space SCF, arbitrary-k, GBT, DOS, and TD-DFT tests remain green.
- [ ] Performance/memory baseline re-run and differences explained.
- [ ] `git diff --check` passes.
- [ ] RF-03 progress section updated and one commit created.

### Required handoff

Include a compact ownership table for each major consumer, the empty-rank semantics, any transitional fields, test matrices, performance change, unticked items, and commit SHA.

### Exact one-line commit message

```text
refactor(RF-03): encapsulate reciprocal k-point ownership
```

---

# RF-04 prompt — batched assembly and reusable workspaces

## Copy-paste implementation prompt

You are implementing **RF-04: batched reciprocal assembly and reusable CPU workspaces** in the RS-LMTO-ASA `fable_v3` branch.

RF-01–RF-03 must be complete and green. Follow all shared rules and keep CUDA disabled. The task converts reciprocal execution to a batch-first CPU shape while retaining the existing public single-point behavior as wrappers and retaining current LAPACK algorithms.

### Inspect before editing

Read:

- all Fourier assembly routines and the second-order `HOH` path;
- overlap construction;
- `calculate_eigenpairs_at_kpoints` and its single-k helper;
- normal mesh construction and diagonalization;
- existing OpenMP regions and thread-private arrays;
- `operator_generation`/cache invalidation;
- existing safe-allocation and finalizer patterns.

### Required design

Introduce two concrete class-like types, with names adjusted to project conventions:

1. `reciprocal_workspace`
   - owns `H`, optional `S`, phase, second-order temporaries, eigenvalue/eigenvector, LAPACK work, real work, and `info` storage;
   - stores matrix dimension, tile capacity, active tile length, generalized/standard mode, and operator-generation fingerprint;
   - provides `ensure_capacity`, `clear`, `restore_to_default`, and a final destructor;
   - grows only when capacity is insufficient and reuses storage otherwise.

2. `reciprocal_assembler`
   - holds non-owning pointers to the Hamiltonian/lattice/control state needed for Fourier assembly;
   - is constructed like `recursion`/`hamiltonian` from longer-lived target objects;
   - provides `assemble_batch`, `assemble_overlap_batch`, and a batch-of-one wrapper;
   - owns only geometry/cache data that genuinely belongs to assembly.

Use the canonical layout `matrix(nmat,nmat,nbatch)` unless measurement demonstrates a compelling alternative. Document that the first two dimensions are Fortran column-major matrices and the third dimension is a contiguous strided batch.

### Required work

1. **Move to batch-first Fourier assembly.**
   - Implement first-order, second-order/HOH, SOC, overlap, and supported CCOR paths over a `kpoint_workset` or its tile.
   - Preserve the existing phase convention and accumulation order inside each matrix where feasible.
   - Keep current single-k procedures as thin batch-of-one wrappers for compatibility.

2. **Eliminate hot-loop allocation.**
   - No allocation/deallocation or LAPACK workspace query is permitted inside the k-point loop after workspace preparation.
   - Cache LAPACK workspace sizes by matrix dimension and solver mode.
   - Reuse Fourier and `HOH` temporaries across points/tiles.

3. **Convert arbitrary-k service to tiles.**
   - Fold/deduplicate once.
   - Assemble unique points in configurable internal CPU tiles.
   - Diagonalize using the existing LAPACK calls and preallocated workspace.
   - Scatter eigenpairs back into the caller's original order.
   - Do not expose a new user option for tile size yet; choose a conservative internal default and allow a test-only/programmatic override if needed.

4. **Keep OpenMP sound.**
   - Do not share mutable LAPACK work arrays between threads.
   - Choose either one workspace per OpenMP worker or tile-level parallelism with clearly exclusive slices.
   - Avoid nested oversubscription between OpenMP and threaded BLAS; document the expected setup but do not control external BLAS environment variables from the code.

5. **Preserve cache invalidation.**
   - Rebuild assembly-side cached data when `hamiltonian%operator_generation`, dimensions, order, SOC/overlap mode, or geometry changes.
   - Add direct tests showing unchanged state reuses capacity/cache and changed operator generation invalidates it.

6. **Tests and measurements.**
   Compare batch sizes 1, a partial tile, an exact full tile, and more than one tile for:
   - first-order Hamiltonians;
   - second-order/HOH Hamiltonians;
   - SOC;
   - generalized overlap;
   - duplicate/permuted arbitrary points;
   - single-site and multi-site dimensions;
   - OpenMP 1 thread and more than 1 thread.
   Check matrices before diagonalization as well as eigenpairs after it. Re-run RF-01 timings and report allocation/call-count changes separately from wall time.

### Explicit non-goals

- No abstract backend yet.
- No cuBLAS/cuSOLVER or C interoperability.
- No device/pinned memory.
- No new generalized-eigenproblem algorithm.
- No TD-DFT transition refactor.

### RF-04 checklist

- [ ] RF-01–RF-03 present and green.
- [ ] Workspace and assembler follow class-based project style.
- [ ] Batch memory layout and ownership documented.
- [ ] First-order, second-order/HOH, SOC, overlap, and supported CCOR assembly covered.
- [ ] Single-k APIs are thin batch-of-one wrappers.
- [ ] No allocation/deallocation in the prepared k loop.
- [ ] LAPACK workspace query is cached by relevant dimensions/mode.
- [ ] Arbitrary-k service tiles unique points and restores caller order.
- [ ] OpenMP workspaces/slices have exclusive ownership.
- [ ] Operator-generation invalidation is tested.
- [ ] Batch sizes 1/partial/full/multiple agree with RF-01 oracle.
- [ ] Serial, OpenMP, MPI, debug, SCF, DOS, GBT, and TD-DFT suites pass.
- [ ] CPU timings and allocation/call counts compared with RF-01.
- [ ] `git diff --check` passes.
- [ ] RF-04 progress section updated and one commit created.

### Required handoff

Report type/API layout, allocation sites removed, workspace-cache keys, OpenMP ownership model, batch-equivalence errors, performance changes, unticked items, and commit SHA.

### Exact one-line commit message

```text
refactor(RF-04): add batched reciprocal assembly workspaces
```

---

# RF-05 prompt — abstract reciprocal execution backend

## Copy-paste implementation prompt

You are implementing **RF-05: an abstract reciprocal execution backend and its complete CPU/LAPACK implementation** in the RS-LMTO-ASA `fable_v3` branch.

RF-01–RF-04 must be complete and green. Follow all shared rules and keep CUDA disabled. The abstract boundary must be useful for a later device-resident Fourier-plus-eigensolver backend; it must not be a cosmetic wrapper around one LAPACK call.

### Inspect before editing

Read:

- the RF-04 workset, assembler, and workspace types;
- `recursion.f90` backend dispatch and the existing CUDA plugin wrapper for style only;
- reciprocal SCF orchestration and arbitrary-k consumers;
- standard and generalized diagonalization paths;
- Hermiticity/overlap diagnostics and cache invalidation.

### Required design

Add a private-by-default module defining an abstract public backend type, for example:

```fortran
type, abstract, public :: reciprocal_execution_backend
contains
   procedure(backend_initialize_if), deferred :: initialize
   procedure(backend_capabilities_if), deferred :: capabilities
   procedure(backend_prepare_if), deferred :: prepare_operator
   procedure(backend_execute_if), deferred :: execute_batch
   procedure(backend_synchronize_if), deferred :: synchronize
   procedure(backend_release_if), deferred :: release
end type reciprocal_execution_backend
```

`execute_batch` should cover the meaningful residency boundary: assembly plus requested eigensolution, with results described by a request/result type. Do not force the caller to pull `H(k)` to host between assembly and solution. Allow tests/diagnostics to request assembled matrices explicitly.

Introduce small concrete request/capability types rather than strings scattered through code. Capabilities should include at least:

- standard Hermitian eigensystem;
- generalized Hermitian eigensystem;
- eigenvalues-only/eigenvectors;
- first/second-order assembly;
- overlap;
- maximum or preferred tile size;
- host/device residency indicator, initially host only.

Implement `lapack_reciprocal_backend` as a concrete extension. It should own/compose the RF-04 assembler and workspace, follow constructor/finalizer style, and reproduce all current CPU behavior.

### Required work

1. **Implement a real CPU backend.**
   - Standard path uses the existing `zheev` behavior.
   - Generalized path uses the existing `zhegv` behavior.
   - Preserve Hermiticity and overlap checks at the same logical boundary.
   - Preserve eigenvector layout, sorting, and failure diagnostics.

2. **Add a backend factory/constructor boundary.**
   - The current default selects only `lapack`/CPU.
   - Do not add a fake GPU backend.
   - A future CUDA factory addition must be possible without changing SCF or TD-DFT physics calls.

3. **Make reciprocal consumers use the abstraction.**
   - Normal k-space SCF construction/diagonalization and arbitrary k/k+q eigensystems should execute through the backend.
   - Keep compatibility type-bound procedures on `reciprocal`, but make them orchestration/adapters rather than duplicated numerical implementations.
   - Avoid modifications to `self.f90`; if a call-site change is unavoidable, keep it minimal.

4. **Keep polymorphism out of hot loops.**
   - Perform one deferred call per tile or large operation, never per matrix element, band pair, or frequency.
   - The concrete CPU backend handles its tile internally.

5. **Capability enforcement.**
   - Validate unsupported requests at the public backend boundary with a concise diagnostic.
   - The CPU backend must report generalized overlap support.
   - Request/result types should make ownership and validity explicit: assembled matrix requested or not, eigenvectors requested or not, local point count, and operator generation.

6. **Tests.**
   Add backend contract tests for:
   - construction/destruction/reconstruction;
   - capability reporting;
   - standard eigenvalues/eigenvectors;
   - generalized eigenvalues and metric orthogonality;
   - eigenvalues-only and eigenvector requests;
   - batch sizes across RF-04 cases;
   - operator-generation refresh;
   - residual, Hermiticity, and orthogonality checks;
   - compatibility APIs returning RF-01 results;
   - deterministic cleanup under debug allocation checks.

### Architecture acceptance criterion

At completion, adding a CUDA backend should require:

- one new concrete backend implementation;
- factory/selection wiring;
- build-system additions;

It should not require rewriting reciprocal SCF, arbitrary-k TD-DFT orchestration, k-point ownership, or physics formulas. Demonstrate this in the progress document with a short projected call flow, but do not implement CUDA.

### Explicit non-goals

- No CUDA/HIP code or C API.
- No automatic backend selection from hardware.
- No precision changes.
- No generalized-problem reformulation.
- No device-side DOS/projections yet.
- No TD-DFT transition backend in this task.

### RF-05 checklist

- [ ] RF-01–RF-04 present and green.
- [ ] Abstract type exists only at the true execution substitution boundary.
- [ ] Request, result, and capability types document ownership/layout.
- [ ] Concrete LAPACK backend follows constructor/finalizer/submodule style.
- [ ] CPU backend supports current standard and generalized modes.
- [ ] SCF and arbitrary-k services dispatch once per tile, not inside hot loops.
- [ ] Existing reciprocal public procedures are adapters without duplicate algorithms.
- [ ] Capability failures occur at the public boundary.
- [ ] Backend construction/destruction and cache refresh tests pass.
- [ ] Matrix/eigenpair/residual/orthogonality contracts match RF-01.
- [ ] Serial/OpenMP/MPI/debug and all affected scientific suites pass.
- [ ] CPU performance is compared with RF-04 and regressions explained.
- [ ] Future CUDA call flow documented without implementing it.
- [ ] `git diff --check` passes.
- [ ] RF-05 progress section updated and one commit created.

### Required handoff

Report the public backend contract, concrete CPU ownership, dispatch granularity, compatibility adapters, capability table, numerical/performance comparisons, unticked items, and commit SHA.

### Exact one-line commit message

```text
refactor(RF-05): introduce reciprocal execution backends
```

---

# RF-06 prompt — shared TD-DFT transition engine

## Copy-paste implementation prompt

You are implementing **RF-06: a shared class-like TD-DFT transition engine for \(\chi^0\) and pair-potential \(\Xi\)** in the RS-LMTO-ASA `fable_v3` branch.

RF-01–RF-05 must be complete and green. Follow all shared rules and keep CUDA disabled. The objective is to factor duplicated transition batching, denominator generation, and matrix accumulation without weakening the separate physics definitions of ordinary response and pair-potential response.

### Inspect before editing

Read:

- `source/tddft_chi0.f90`;
- `source/tddft_xi.f90`;
- `source/response_vertices.f90` or the current vertex module;
- `source/tddft_four_component.f90`;
- `source/tddft_chi0_green.f90` only to preserve its separate reference status;
- TD-DFT orchestration in `source/calculation.f90`;
- all TD-DFT unit tests and validation campaign fixtures;
- RF-01 TD-DFT numerical contract.

### Required design

Introduce:

1. A concrete `tddft_transition_engine` type that owns a reusable transition workspace and implements static/dynamic response accumulation.
2. A concrete `tddft_transition_workspace` containing bounded tiles of:
   - left and right vertices;
   - occupation factors;
   - transition energies;
   - complex denominators;
   - weighted left vertices;
   - any GEMM scratch/result tile.
3. An abstract `tddft_vertex_provider` with deferred **tile-level** operations, not per-transition scalar dispatch.
4. Concrete providers for:
   - site/channel response vertices used by \(\chi^0\);
   - pair-potential weighted right vertices used by \(\Xi\), including k-dependent operators.

The provider should fill a whole transition tile or a whole k/band tile per deferred call. The engine should not make a polymorphic call for every `(k,n,m,frequency)` element.

Follow project style: constructor, type-bound procedures, `restore_to_default`, final destructor, parent-module declarations, submodule implementations, Doxygen ownership/layout documentation.

### Required work

1. **Extract shared transition selection.**
   One implementation should own:
   - band-window normalization;
   - k-weight normalization;
   - Fermi occupations;
   - optional occupation pruning;
   - transition-energy generation;
   - deterministic transition ordering;
   - transition-tile flushing.

2. **Extract shared denominator and accumulation.**
   One implementation should generate

   \[
   \omega + \epsilon_{n\mathbf{k}}-\epsilon_{m,\mathbf{k+q}}+i\eta
   \]

   and perform the tiled matrix accumulation for both \(\chi^0\) and \(\Xi\). Preserve the existing no-conjugation/transposition convention: the right circular-channel factor is computed independently as \(\langle m|B|n\rangle\), not assumed to be the conjugate of the left factor.

3. **Preserve static response.**
   Route static \(q=0,\omega=0\) response through a shared static mode that uses the existing divided-difference limit. Do not emulate it with finite eta or a small dynamic frequency.

4. **Convert current public routines into adapters.**
   Keep public entry points and result types used by callers/tests where possible:
   - `build_chi_ks_from_eigenpairs`;
   - `build_static_chi_ks_from_eigenpairs`;
   - `build_direct_xi_from_eigenpairs`;
   - `build_direct_xi_from_k_dependent_eigenpairs`;
   - static direct Xi variants.
   These should configure a provider and engine rather than contain duplicated band/denominator/GEMM loops.

5. **Prepare pair-operator streaming without rewriting atomic LMTO.**
   - The pair provider must accept an operator tile or k-slice and consume it without requiring a second complete copy.
   - Retain adapters for existing full-array callers.
   - Do not refactor `build_lmto_pair_potential_at_kpoint` or legacy tangent physics.
   - Document the later producer/consumer path in which CPU LMTO produces a tile and a GPU response backend consumes it.

6. **Preserve scalar oracle.**
   Keep a selectable/internal scalar reference implementation for tests. It may live in the engine, but it must remain independent enough to catch batching mistakes.

7. **Metadata and timing.**
   Preserve all current metadata fields and meanings. Time at least vertex fill, occupation/transition preparation, denominator generation, and accumulation. Ensure timing does not introduce calls inside the innermost loops when profiling is disabled/default.

8. **Tests.**
   Compare old/scalar and shared-engine paths for:
   - ordinary transverse \(\chi^{+-}\);
   - multi-channel/four-component response where it uses the same eigenpair machinery;
   - static divided-difference response;
   - constant and k-dependent pair-potential \(\Xi\);
   - batch sizes 1, partial, exact, and multiple tiles;
   - finite and zero temperature;
   - eta ladder and near-degenerate transitions;
   - complete and restricted band windows;
   - pruning disabled and explicit approximate pruning;
   - spectral products, Ward/Goldstone diagnostics, Dyson/mode downstream tests;
   - existing validation campaign outputs.

9. **Do not absorb the Green reference backend.**
   `tddft_chi0_green.f90` remains a separate eigenpair-resolvent reference. Only share small convention helpers if doing so is clearly safe and within scope.

### Architecture acceptance criterion

At completion, a future GPU transition engine should be able to consume eigenpair/operator tiles from the RF-05 backend and produce partial response matrices without changing vertex physics, public TD-DFT result semantics, or `calculation.f90`'s scientific decisions. Document the projected interface; do not implement GPU code.

### Explicit non-goals

- No GPU implementation.
- No q-group/k-worker MPI redesign.
- No changes to Goldstone correction physics.
- No longitudinal capability expansion.
- No replacement of the Green reference backend.
- No changes to eta, frequency, occupation, or spectral sign conventions.

### RF-06 checklist

- [ ] RF-01–RF-05 present and green.
- [ ] Transition engine/workspace follow class-based project style.
- [ ] Vertex abstraction dispatches per tile, not per scalar transition.
- [ ] Site/channel and pair-potential providers are implemented.
- [ ] Band selection, occupations, pruning, and transition order are shared.
- [ ] Dynamic denominator and GEMM accumulation are shared.
- [ ] Static divided-difference response remains mathematically distinct.
- [ ] Existing public χ⁰/Xi APIs are thin adapters.
- [ ] Pair provider accepts k/operator tiles without duplicating the full tensor.
- [ ] Legacy LMTO/tangent construction is untouched.
- [ ] Scalar reference path remains available to tests.
- [ ] Metadata and timing meanings are preserved.
- [ ] χ⁰, Xi, four-component, Ward, Goldstone, Dyson, modes, and validation tests pass.
- [ ] Batch/tile variants agree at RF-01 tolerances.
- [ ] Serial/OpenMP/MPI/debug builds pass.
- [ ] CPU performance and peak-array-size changes are reported.
- [ ] Future device-resident TD-DFT call flow documented.
- [ ] `git diff --check` passes.
- [ ] RF-06 progress section updated and one commit created.

### Required handoff

Report the engine/provider/workspace API, dispatch granularity, formulas/conventions explicitly preserved, public adapters, pair-operator tile semantics, numerical/performance comparisons, unticked items, and commit SHA.

### Exact one-line commit message

```text
refactor(RF-06): unify TDDFT transition response kernels
```

---

# Final gate before GPU implementation

The reciprocal GPU backend may begin only when all of the following are true:

- [x] RF-01 through RF-06 each exist as a separate reviewed commit.
- [x] Serial release, OpenMP, debug, MPI+OpenMP, and MPI-only builds are green.
- [x] K-space SCF has serial/MPI equivalence coverage.
- [x] Reciprocal matrices and eigenpairs meet residual/orthogonality contracts.
- [x] TD-DFT scalar and tiled engines agree across the full response test surface.
- [x] No hot-loop allocation remains in arbitrary-k eigensolution or transition accumulation.
- [x] K-point ownership and empty-rank semantics are explicit.
- [x] The CPU backend is selected through one backend factory boundary.
- [x] A future backend can keep Fourier assembly and eigensolution in one resident tile.
- [x] Pair-potential operators can be consumed as tiles.
- [x] Timings identify the measured CPU bottlenecks on production-relevant fixtures.
- [x] Known CI failures have been resolved or explicitly demonstrated unrelated.

## Final gate review — 2026-08-12

The initial adversarial review correctly blocked the gate on strict-Debug
defects, per-k/transition-tile allocations, absent resident normal-mesh
execution, and a non-streamed pair-operator contract. GC-01 through GC-05
closed those findings in focused reviewed commits. GC-06 qualified the result
on source commit `464cbe7f7193ea433f154542c80fb0c2e7f86bc7`.

The complete CPU-only qualification matrix passed with CUDA disabled:

- Release serial/OpenMP: 105/105 enabled tests passed.
- Strict Debug with GNU runtime checks and invalid/zero/overflow traps:
  105/105 enabled tests passed.
- Release MPI/OpenMP: 116/116 enabled tests passed.
- Release MPI-only: 116/116 enabled tests passed.

The MPI matrices include bcc-Fe k-space SCF at serial and MPI ranks 1, 2, and
4. `UnitTddftCpuProfile` records the one-site bcc-Fe and two-site fcc-Ni
payloads and phases. Seven repository-disabled WP8/WP9 tests are explicitly
excluded from all enabled counts; they are not GPU-reciprocal gate failures.

The full evidence and independent-review record are in
[`GPU_RECIPROCAL_GATE_CLOSURE_PLAN.md`](plans/GPU_RECIPROCAL_GATE_CLOSURE_PLAN.md)
at commits `44e958b` and `79ce8a1`.

**Decision: GPU cleared.** The next task is one concrete FP64, `ham_only`
reciprocal backend; it must preserve the RF numerical contracts above rather
than modify the established scientific orchestration.
