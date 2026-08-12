# Reciprocal GPU gate-closure plan

## Objective

Close the blockers identified by the 2026-08-12 adversarial review before any
CUDA source, plugin, or GPU-specific build option is introduced.  This plan
deliberately follows RF-01 through RF-06; it does not amend or combine those
six commits.

| Stage | Scope | Depends on | Commit |
| --- | --- | --- | --- |
| GC-01 | Debug-safe TD-DFT workspace lifecycle | — | `fix(GC-01): make TDDFT workspace capacity guard debug-safe` |
| GC-02 | Reusable overlap-validation scratch | GC-01 | `refactor(GC-02): reuse reciprocal overlap validation scratch` |
| GC-03 | Allocation-free TD-DFT vertex tiles | GC-01 | `refactor(GC-03): reuse TDDFT vertex tile scratch` |
| GC-04 | Resident normal-mesh reciprocal tile path | GC-02 | `refactor(GC-04): fuse normal-mesh reciprocal tile execution` |
| GC-05 | Streamed pair-potential operator tiles | GC-03, GC-04 | `refactor(GC-05): stream pair-potential operator tiles` |
| GC-06 | Independent qualification and final review record | GC-01–GC-05 | `test(GC-06): qualify reciprocal GPU gate closure` |

GC-02 and GC-03 can be developed independently after GC-01, but they must be
rebased and reviewed as separate commits.  GC-04 and GC-05 must preserve all
RF numerical contracts.  Only GC-06 may change the final-gate checklist.

## Rules for every stage

- Start from the current gate-closure HEAD; record `git rev-parse HEAD` and
  `git status --short` before editing.  Preserve user-owned files and all
  unrelated changes.
- Keep `ENABLE_CUDA_PLUGIN=OFF`.  Do not add CUDA/HIP files, GPU mocks,
  user-facing backend options, changed physical defaults, or reference churn.
- Preserve FP64, retarded TD-DFT denominators, occupation/window semantics,
  explicit k-weight normalization, left/right vertex order, and the scalar
  TD-DFT oracle.
- One GC stage is one focused commit.  Do not fold validation evidence for a
  later stage into an earlier commit.
- Update this plan's matching stage section in the same commit.  Record exact
  commands, compiler/runtime, pass counts, and any skipped test with its
  reason.
- Treat a Debug failure as a correctness defect unless a minimal standalone
  reproducer proves it is outside this codebase.  Do not call an unrecorded
  failure "environmental" merely because Release passes.

## Baseline facts to retain

At `5f1ff71`, Release/OpenMP units pass 42/42.  MPI+OpenMP and MPI-only units
pass 46/46 when OpenMPI may create local PMIx sockets.  The bcc-Fe k-space SCF
oracle passes serial and MPI ranks 1, 2, and 4.  The current Debug unit suite
fails seven TD-DFT tests at `tddft_transition_engine.f90:329` because an
unallocated allocatable component is passed to `size` in a compound logical
expression.  This is the first fix, not a test waiver.

---

## GC-01 — Debug-safe TD-DFT workspace lifecycle

### Outcome

`tddft_transition_workspace%ensure_capacity` must work for an uninitialized
workspace under GNU runtime checking and remain safe for reuse and resizing.
The fix must not change the response formula or hide the test with disabled
checks.

### Copy-paste implementation prompt

> You are implementing **GC-01: Debug-safe TD-DFT workspace lifecycle**.
> Work only on the unsafe capacity/reuse guard and its direct tests.  Do not
> change TD-DFT physics, tiling strategy, allocation ownership, public input,
> or backend architecture.
>
> First reproduce the failure with the configured Debug build and
> `ctest --test-dir build-rf-debug -R 'Unit(TddftChiKS|TddftDirectXi)' --output-on-failure`.
> Inspect `source/tddft_transition_engine.f90`, `tests/unit/test_tddft_chi_ks.f90`,
> and `tests/unit/test_tddft_direct_xi.f90` before editing.
>
> Replace the compound `allocated(...) .and. size(...)` reuse condition with
> control flow that never evaluates `size`, `shape`, or an array section until
> the corresponding allocatable component is known to be allocated.  Remember
> that Fortran does not guarantee short-circuit evaluation of `.and.`.  The
> valid reuse condition must require: sufficient capacity; all required
> workspace components allocated; and matching left/right extents.  A partial
> workspace must be cleared and rebuilt rather than reused.
>
> Add a focused unit test that exercises all of: first use of a default
> workspace, same-shape reuse, changed capacity, changed response dimensions,
> `clear`, and reuse after `clear`.  Test observable allocation/shape state;
> do not make a source-text test the only proof.  Keep normal numerical tests
> unchanged except where they can directly use the new lifecycle test.
>
> Required validation:
>
> ```bash
> cmake --build build-rf-debug --parallel
> ctest --test-dir build-rf-debug -L unit --output-on-failure -j 1
> cmake --build build-rf-serial --parallel
> ctest --test-dir build-rf-serial -R 'Unit(TddftChiKS|TddftDirectXi|TddftCpuProfile|TddftGreenChiKS|TddftFourComponent)|TddftCrossMilestoneEquivalence' --output-on-failure
> git diff --check
> ```
>
> Record the original runtime-check failure, the safe guard invariant, and
> exact test results in the GC-01 section of this document.  Commit exactly:
>
> ```text
> fix(GC-01): make TDDFT workspace capacity guard debug-safe
> ```

### Acceptance checklist

#### Evidence — 2026-08-12, `5f1ff71e57de997f32c6813c6c79305c52fbd0b2` + GC-01 worktree

The original reproducer failed `UnitTddftChiKS` and `UnitTddftDirectXi` in
`tddft_transition_engine.f90:329` under GNU Fortran 13.3.0 Debug runtime
checking: the compound allocation/`size` guard evaluated `size` on an
unallocated component.  The replacement checks allocation of every component
before entering a separate shape/capacity check; a partial workspace falls
through to `clear` and complete reallocation.

`UnitTddftTransitionWorkspace` covers default first use, same-shape reuse,
capacity and response-dimension changes, an intentionally partial workspace,
`clear`, and reuse after `clear` through observable allocation and shape
state.  `cmake --build build-rf-debug --parallel` followed by
`ctest --test-dir build-rf-debug -L unit --output-on-failure -j 1` passed
43/43 unit tests.  `cmake --build build-rf-serial --parallel` followed by
`ctest --test-dir build-rf-serial -R 'Unit(TddftChiKS|TddftDirectXi|TddftCpuProfile|TddftGreenChiKS|TddftFourComponent)|TddftCrossMilestoneEquivalence' --output-on-failure`
passed 6/6 Release TD-DFT tests.  CUDA remained disabled; no references or
runtime checks changed.

- [x] Default, cleared, partially initialized, reused, and resized workspaces are safe.
- [x] All seven previously failing Debug TD-DFT tests pass.
- [x] Release TD-DFT numerical results remain unchanged within RF-01 tolerances.
- [x] No check is weakened and no reference is regenerated.
- [ ] GC-01 evidence is recorded and the commit is independently reviewed.

---

## GC-02 — Reusable overlap-validation scratch

### Outcome

The generalized arbitrary-k eigensolution must not allocate a Cholesky copy
once per k point.  It must retain validation: `S(k)` remains Hermitian and
positive definite before `zhegv`.

### Copy-paste implementation prompt

> You are implementing **GC-02: reusable reciprocal overlap-validation
> scratch**.  Do not change generalized-overlap mathematics, the LAPACK
> routine, tolerance, arbitrary-k folding/deduplication, or TD-DFT code.
>
> Inspect `source/reciprocal_backend.f90`, `source/reciprocal.f90`,
> `source/reciprocal_fourier.f90`, and
> `tests/unit/test_arbitrary_k_eigenpairs.f90`.  Confirm that
> `lapack_backend_execute_batch` calls `assert_overlap` in its per-k solve
> loop and that `assert_overlap` currently allocates `chol`.
>
> Make the Cholesky copy owned by `reciprocal_workspace` (or an equivalently
> long-lived backend workspace), sized in `ensure_capacity`, cleared in its
> lifecycle method, and reused for every overlap validation in a prepared
> tile.  Pass that storage explicitly to the overlap validator; never mutate
> the source overlap matrix merely to validate it.  The per-k path may copy
> values into preallocated scratch but may not allocate/deallocate or query
> LAPACK workspace.
>
> Extend the generalized two-site arbitrary-k fixture to execute multiple
> prepared tiles and assert stable workspace allocation/query counters before
> and after the repeated solve.  Add a narrow source-contract test only if the
> existing runtime counters cannot observe the Cholesky scratch; it must target
> the validator body rather than ban legitimate allocation elsewhere.  Keep
> the existing residual and metric-orthogonality checks.
>
> Required validation:
>
> ```bash
> cmake --build build-rf-serial --parallel
> ctest --test-dir build-rf-serial -R 'UnitArbitraryKEigenpairs' --output-on-failure
> cmake --build build-rf-debug --parallel
> ctest --test-dir build-rf-debug -R 'UnitArbitraryKEigenpairs' --output-on-failure
> bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -R "UnitArbitraryKEigenpairs|UnitKspaceOccupations_mpi" --output-on-failure'
> git diff --check
> ```
>
> Record storage ownership, the validation invariant, and counter evidence in
> GC-02.  Commit exactly:
>
> ```text
> refactor(GC-02): reuse reciprocal overlap validation scratch
> ```

### Acceptance checklist

#### Evidence — 2026-08-12, `8b3f8082232a751c69851e9b0f6d68ada76e7a88` + GC-02 worktree

`reciprocal_workspace` now owns `overlap_cholesky(nmat,nmat)`, allocated only
while preparing a generalized tile and released by `clear`.  The backend passes
that scratch explicitly to `assert_overlap`; it first copies `S(k)` into the
scratch, verifies Hermiticity, and applies `ZPOTRF` there.  Thus the source
overlap matrix remains intact for the following `ZHEGV`, and the per-k solve
loop has neither allocation/deallocation nor a workspace query.

The generalized two-site arbitrary-k fixture now checks the prepared scratch
shape, generalized residual and metric orthogonality, and stable
`storage_allocations`/`lapack_workspace_queries` across a repeated prepared
tile sequence.  GNU Fortran 13.3.0 validation passed:

```text
cmake --build build-rf-serial --parallel
ctest --test-dir build-rf-serial -R 'UnitArbitraryKEigenpairs' --output-on-failure
# 1/1 passed

cmake --build build-rf-debug --parallel
ctest --test-dir build-rf-debug -R 'UnitArbitraryKEigenpairs' --output-on-failure
# 1/1 passed

bash -lc 'source env/openmpi.sh; cmake --build build-rf-mpi --target UnitArbitraryKEigenpairs --parallel'
bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -R "UnitArbitraryKEigenpairs|UnitKspaceOccupations_mpi" --output-on-failure'
# 2/2 passed outside the sandbox
```

The same MPI ctest inside the sandbox failed only when OpenMPI could not create
its PMIx listener (`opal_ifinit: socket() failed with errno=1`); after the
approved unsandboxed rerun, both rebuilt tests passed.  CUDA remained disabled;
no overlap tolerance, LAPACK routine, reference, or TD-DFT code changed.

- [x] Generalized per-k solve path has no allocate/deallocate or workspace query.
- [x] `S(k)` validation remains a copied Cholesky test and preserves `S(k)`.
- [x] Generalized residual and metric-orthogonality tests pass.
- [x] Repeated prepared arbitrary-k tiles show no new workspace allocation.
- [ ] GC-02 evidence is recorded and the commit is independently reviewed.

---

## GC-03 — Allocation-free TD-DFT vertex tiles

### Outcome

Dynamic and static transition accumulation allocate all scratch once per
engine/workspace preparation, not once per transition tile.  Vertex providers
continue to dispatch per completed tile, never per scalar transition.

### Copy-paste implementation prompt

> You are implementing **GC-03: allocation-free TD-DFT vertex tiles**.  Do
> not change TD-DFT physics, scalar-oracle order, batching semantics, BLAS
> formulas, public chi_KS/Xi APIs, or pair-operator representation yet.
>
> Inspect `source/tddft_transition_engine.f90`, `source/tddft_chi0.f90`,
> `source/tddft_xi.f90`, and the chi_KS/direct-Xi/four-component unit tests.
> The site-channel and pair-operator `fill_vertex_tile` implementations
> currently allocate `bra` and `ket` on every flushed transition tile.
>
> Extend `tddft_transition_workspace` with reusable coefficient-space
> `bra`/`ket` tile storage.  Give the transition engine a safe way to learn
> the provider coefficient dimension once per accumulation call (for example,
> a tile-level provider capability/query), then size all response- and
> coefficient-space storage in `ensure_capacity`.  Change the tile interface
> to receive caller-owned scratch slices; providers must fill those slices and
> must not allocate/deallocate.  Any added dispatch may occur once per engine
> call or completed transition tile, never per `(k,n,m,omega)` element.
>
> Keep the scalar reduction as an independent numerical oracle.  Add explicit
> workspace allocation/reuse counters if necessary and test dynamic and
> static chi_KS plus constant and k-dependent Xi with more than one k point,
> a partial final tile, and repeated identical calls.  The test must prove no
> new workspace allocation after preparation and preserve scalar/GEMM equality.
>
> Required validation:
>
> ```bash
> cmake --build build-rf-serial --parallel
> ctest --test-dir build-rf-serial -R 'Unit(TddftChiKS|TddftDirectXi|TddftFourComponent|TddftCpuProfile)|TddftCrossMilestoneEquivalence' --output-on-failure
> cmake --build build-rf-debug --parallel
> ctest --test-dir build-rf-debug -L unit --output-on-failure -j 1
> bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -L unit --output-on-failure -j 1'
> git diff --check
> ```
>
> Record workspace layout, dispatch granularity, allocation evidence, and
> numerical tolerances in GC-03.  Commit exactly:
>
> ```text
> refactor(GC-03): reuse TDDFT vertex tile scratch
> ```

### Acceptance checklist

#### Evidence — 2026-08-12, `128e66a75efb1555767410a858de50dc52f621f7` + GC-03 worktree

`tddft_transition_workspace` now owns coefficient-space `bra` and `ket`
tiles with the response tile arrays.  Each provider exposes its coefficient
dimension once per `accumulate_dynamic`/`accumulate_static` call; the engine
then prepares all storage together.  The tile interface receives caller-owned
`bra`/`ket` buffers and both site-channel and pair-operator providers only fill
their active slices.  They perform no allocation/deallocation in a tile flush;
provider dispatch remains once per completed transition tile.

`UnitTddftTransitionWorkspace` now uses real site-channel chi_KS and constant
pair-operator Xi providers over two k points with batch size 3 (including the
partial final tile).  It verifies prepared `(ncoefficient,capacity)` scratch,
unchanged allocation counters and stable numerical results for repeated dynamic
and static calls.  Existing scalar-vs-GEMM coverage remains in the chi_KS and
direct-Xi tests.  GNU Fortran 13.3.0 validation passed:

```text
cmake --build build-rf-serial --parallel
ctest --test-dir build-rf-serial -R 'Unit(TddftChiKS|TddftDirectXi|TddftFourComponent|TddftCpuProfile)|TddftCrossMilestoneEquivalence' --output-on-failure
# 5/5 passed

cmake --build build-rf-debug --parallel
ctest --test-dir build-rf-debug -L unit --output-on-failure -j 1
# 43/43 passed

bash -lc 'source env/openmpi.sh; cmake --build build-rf-mpi --parallel'
bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -L unit --output-on-failure -j 1'
# 47/47 passed outside the sandbox
```

The MPI run used the approved unsandboxed OpenMPI runtime because the sandbox
blocks PMIx listener socket creation (recorded in GC-02).  CUDA remained
disabled; no TD-DFT response formulas, scalar reduction, batching semantics,
or reference results changed.

- [x] Site and pair providers allocate no `bra`/`ket` scratch in a tile flush.
- [x] Repeated dynamic/static chi_KS and Xi calls reuse prepared storage.
- [x] Scalar and GEMM/tiled paths agree over the RF-01 response surface.
- [x] Debug TD-DFT suite passes with checks enabled.
- [ ] GC-03 evidence is recorded and the commit is independently reviewed.

---

## GC-04 — Resident normal-mesh reciprocal tile path

### Outcome

The ordinary k-space SCF mesh crosses the backend boundary once per tile with
both Fourier assembly and eigensolution requested.  A future device backend
can retain H/S and eigensolver data on one resident tile; host compatibility
caches are materialized only after that operation, never between assembly and
solution.

### Copy-paste implementation prompt

> You are implementing **GC-04: resident normal-mesh reciprocal tile
> execution**.  Do not add a GPU backend or a new user-facing option.  Keep
> `ham_only` FP64 as the default and preserve all legacy public reciprocal
> procedures as compatibility adapters.
>
> Inspect `source/reciprocal.f90`, `source/reciprocal_backend.f90`,
> `source/reciprocal_fourier.f90`, `source/reciprocal_bands.f90`,
> `source/reciprocal_lifecycle.f90`, and the bcc-Fe k-space SCF/arbitrary-k
> tests.  Today normal mesh construction requests an assembled host H tile and
> a later call allocates/copies it back as `input_hamiltonian` for eigensolve.
>
> Introduce one internal normal-mesh service that, for each owned workset tile,
> makes a single `execute_batch` request with assembly and eigensolution both
> enabled.  It must populate the existing `hk_bulk`, eigenvalue, and
> eigenvector compatibility caches with the same local/global ownership and
> operator-generation rules as before.  Make `build_kspace_hamiltonian` and
> `diagonalize_hamiltonian` adapters to that service, or use an equivalently
> clear cache-validity state, so unchanged scientific callers cannot trigger a
> second solve.  A stale operator generation, k-path, generalized mode, empty
> local workset, or request for a diagnostic matrix must have explicit,
> tested behavior.
>
> Preserve the abstract request/result boundary.  The CPU backend may copy
> final compatibility output to host; the contract must not require a host H
> copy between Fourier assembly and eigensolution.  Add a backend/test
> observable proving one combined request per normal tile and no
> assemble-only-then-input-H sequence on the normal SCF route.  Compare fused
> and pre-GC-04 oracle H/eigenpairs, residuals, orthogonality, occupations, and
> bcc-Fe SCF observables.
>
> Required validation:
>
> ```bash
> cmake --build build-rf-serial --parallel
> ctest --test-dir build-rf-serial -R 'Unit(ArbitraryKEigenpairs|KspaceOccupations)|Example_k_space_scf_bccFe' --output-on-failure
> cmake --build build-rf-debug --parallel
> ctest --test-dir build-rf-debug -R 'Unit(ArbitraryKEigenpairs|KspaceOccupations)' --output-on-failure
> bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -R "Unit(KspaceOccupations_mpi|ParallelContext_mpi_[124])|Example_k_space_scf_bccFe(_mpi(_[14])?)?" --output-on-failure -j 1'
> git diff --check
> ```
>
> Record the combined request semantics, cache-generation invariant, host-copy
> boundary, and numerical comparison in GC-04.  Commit exactly:
>
> ```text
> refactor(GC-04): fuse normal-mesh reciprocal tile execution
> ```

### Acceptance checklist

- [ ] Normal mesh requests assembly and eigensolution together once per tile.
- [ ] No host H tile is required between those two backend operations.
- [ ] Compatibility H/eigenpair caches and invalidation behavior are unchanged.
- [ ] Serial/MPI bcc-Fe SCF and reciprocal contracts match the RF-01 oracle.
- [ ] GC-04 evidence is recorded and the commit is independently reviewed.

---

## GC-05 — Streamed pair-potential operator tiles

### Outcome

The TD-DFT pair provider accepts a one-k operator tile from a producer without
requiring a complete `(nmat,nmat,nright,nk)` tensor.  Existing constant and
full-tensor APIs remain thin compatibility adapters.

### Copy-paste implementation prompt

> You are implementing **GC-05: streamed pair-potential operator tiles**.
> Do not change pair-potential physics, endpoint phases, q conventions,
> left/right vertex order, public user input, or the LMTO construction itself.
>
> Inspect `source/tddft_transition_engine.f90`, `source/tddft_xi.f90`,
> `source/reciprocal_fourier.f90`, `source/lmto_pair_potential.f90`, and
> `tests/unit/test_tddft_direct_xi.f90`.  The current provider retains a
> pointer to the complete rank-four operator tensor and selects one k slice;
> that is not a streamed tile contract.
>
> Define a narrow, documented pair-operator tile source boundary.  It must
> fill a caller-owned `(nmat,nmat,nright)` operator tile for one local k index
> and may report its coefficient/channel dimensions.  The pair provider must
> prepare/fetch that tile once per k point, retain it for every transition
> batch at that k, and consume it without copying a full k-resolved tensor.
> Keep source/provider dispatch outside the scalar transition and frequency
> loops.  A constant source and a legacy rank-four cached source must adapt the
> existing direct-Xi public APIs so their results and metadata do not change.
>
> Provide a deterministic mock streaming source in unit tests.  It must prove
> that: (1) its output is fetched exactly once per k point, not once per
> transition batch; (2) partial final transition tiles see the same operator;
> (3) the streamed result equals the current k-resolved Xi oracle; and (4) a
> constant source still equals the constant-operator oracle.  Test finite q,
> complex operators, static Xi, dynamic Xi, and scalar/GEMM equivalence.
>
> Required validation:
>
> ```bash
> cmake --build build-rf-serial --parallel
> ctest --test-dir build-rf-serial -R 'Unit(TddftDirectXi|TddftChiKS|LmtoPairPotential|TddftCpuProfile)|TddftCrossMilestoneEquivalence' --output-on-failure
> cmake --build build-rf-debug --parallel
> ctest --test-dir build-rf-debug -L unit --output-on-failure -j 1
> bash -lc 'source env/openmpi.sh; ctest --test-dir build-rf-mpi -L unit --output-on-failure -j 1'
> git diff --check
> ```
>
> Record source ownership/lifetime, fetch granularity, legacy-adapter
> semantics, memory bound, and numerical comparisons in GC-05.  Commit
> exactly:
>
> ```text
> refactor(GC-05): stream pair-potential operator tiles
> ```

### Acceptance checklist

- [ ] Pair provider consumes one caller-owned operator tile per k point.
- [ ] No complete k-resolved operator tensor is required by the new path.
- [ ] Constant and legacy tensor inputs remain compatible adapters.
- [ ] Dynamic/static Xi and scalar/tiled tests cover streamed tiles.
- [ ] GC-05 evidence is recorded and the commit is independently reviewed.

---

## GC-06 — Final qualification and review record

### Outcome

Produce auditable, current evidence for every final-gate box.  This stage
contains no scientific implementation changes.  If any prerequisite fails,
leave the relevant box unchecked and fix it in a new focused task.

### Copy-paste implementation prompt

> You are implementing **GC-06: reciprocal GPU gate qualification**.  Do not
> refactor production code, regenerate references, or add GPU code.  Audit
> GC-01 through GC-05 and update only their evidence plus the final gate.
>
> Confirm each GC commit is separate, has an independent review record, and
> has no mixed unrelated change.  Review `tests/KNOWN_ISSUES.md` and current
> CI configuration.  For every remaining known failure, record its exact
> reproducer, present status, and evidence that it is unrelated to the
> reciprocal/TD-DFT gate; otherwise leave the final CI box unchecked.
>
> Configure clean, separate directories and run the complete matrix with CUDA
> disabled:
>
> ```bash
> cmake -S . -B build-gc-serial -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA_PLUGIN=OFF -DENABLE_MPI=OFF -DENABLE_OPENMP=ON -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON -DRUN_EXAMPLE_TESTS=ON
> cmake --build build-gc-serial --parallel
> OMP_NUM_THREADS=4 ctest --test-dir build-gc-serial --output-on-failure
>
> cmake -S . -B build-gc-debug -DCMAKE_BUILD_TYPE=Debug -DENABLE_CUDA_PLUGIN=OFF -DENABLE_MPI=OFF -DENABLE_OPENMP=OFF -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON -DRUN_EXAMPLE_TESTS=ON -DCMAKE_Fortran_FLAGS_DEBUG='-O0 -g -fcheck=all,no-recursion -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan'
> cmake --build build-gc-debug --parallel
> ctest --test-dir build-gc-debug --output-on-failure -j 1
>
> cmake -S . -B build-gc-mpi -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA_PLUGIN=OFF -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON -DRUN_EXAMPLE_TESTS=ON
> cmake --build build-gc-mpi --parallel
> OMP_NUM_THREADS=4 ctest --test-dir build-gc-mpi --output-on-failure -j 1
>
> cmake -S . -B build-gc-mpi-noomp -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA_PLUGIN=OFF -DENABLE_MPI=ON -DENABLE_OPENMP=OFF -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON -DRUN_EXAMPLE_TESTS=ON
> cmake --build build-gc-mpi-noomp --parallel
> ctest --test-dir build-gc-mpi-noomp --output-on-failure -j 1
> ```
>
> If sandboxing prevents OpenMPI from creating PMIx sockets, repeat only the
> MPI commands in an approved unsandboxed runtime and record both the failure
> signature and successful external result.  Explicitly rerun bcc-Fe k-space
> SCF serial plus MPI ranks 1, 2, and 4; arbitrary-k standard/generalized
> residual and orthogonality tests; all TD-DFT scalar/tile tests; and the CPU
> profile.  Record the dominant phases, matrix/k/frequency dimensions, and
> principal array payload.
>
> Tick a final-gate box only with a command, current commit SHA, and reviewer
> evidence.  The final decision must be either "GPU cleared: implement one
> FP64 ham_only reciprocal backend" or "GPU blocked", followed by the exact
> unchecked item.  Commit exactly:
>
> ```text
> test(GC-06): qualify reciprocal GPU gate closure
> ```

### Acceptance checklist

- [ ] Four complete build/test matrices pass on the same qualified revision.
- [ ] Every no-allocation and tile-residency claim has direct test evidence.
- [ ] K-space SCF serial/MPI, reciprocal, TD-DFT, and profile evidence is current.
- [ ] CI failures are closed or rigorously classified unrelated.
- [ ] Final-gate boxes and reviewer record are complete and truthful.

## Review record template

Use this template for each GC commit in the handoff and in the commit review:

| Field | Record |
| --- | --- |
| Commit SHA | |
| Reviewer and date | |
| Scope checked | |
| Numerical contracts rechecked | |
| Debug result | |
| MPI result | |
| Allocation/residency evidence | |
| Open findings / follow-up issue | |
| Approved for next GC stage | Yes / No |
