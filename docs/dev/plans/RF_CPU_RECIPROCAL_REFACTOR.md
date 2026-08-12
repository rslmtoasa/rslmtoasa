# RF CPU reciprocal refactor progress

## Six-stage overview

| Stage | Scope | Status |
| --- | --- | --- |
| RF-01 | Characterization, numerical oracles, and CPU baseline | Complete |
| RF-02 | MPI runtime state, build integration, and output ownership | Complete |
| RF-03 | Explicit k-point worksets and ownership | Complete (MPI runtime validation pending environment repair) |
| RF-04 | Batched reciprocal assembly and reusable workspaces | Pending RF-03 |
| RF-05 | Reciprocal execution backend with CPU/LAPACK implementation | Pending RF-04 |
| RF-06 | Shared TD-DFT transition engine for chi0 and Xi | Pending RF-05 |

## Baseline record

- **Baseline commit:** `c3856036dee703f2ca3ee5ce3a6dc05f71801e8f` (`fable_v3`), recorded before edits on 2026-08-11.
- **Starting worktree:** two user-owned untracked planning documents, preserved unchanged.
- **Host:** Intel(R) Core(TM) i9-11900 @ 2.50 GHz; 16 logical CPUs.
- **Toolchain:** GNU Fortran 13.3.0; CMake 3.28.3; GNU OpenMP 4.5; Intel oneMKL 2026.1 LP64 BLAS/LAPACK; Python 3.12.3 and pytest 9.1.1. CUDA plugin disabled.

### Validation commands and results

Serial Release configuration:

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

The Release build passed. The first CTest attempt stopped at `Lanczos` because `pytest` was not on `PATH`; pytest 9.1.1 was then installed in the active Python environment. The focused RF-01 command then passed 5/5 tests:

```bash
ctest --test-dir build-rf-serial \
  -R 'Unit(ArbitraryKEigenpairs|KspaceOccupations|TddftChiKS|TddftDirectXi|TddftCpuProfile)' \
  --output-on-failure
```

GNU 14.2.0 Release/OpenMP validation was completed in `build-rf-release14`.
The unit suite passed 40/40; the regression suite passed 17/17 after excluding
the user-waived WP8/WP9 tests; and the post-processing suite passed 12/12.
The SCF suite initially exposed five invalid GBT fixtures. They were repaired
on 2026-08-11 (`strux_lib`; `nsp=3`; direct spirals explicitly select
`gbt_single_q`) and their reviewed goldens were regenerated. A focused CTest
rerun passed all five. The previously reported fcc-Cu "failures" were not
CTest failures: the `abs_tol`/`rel_tol` fields stored in `cases.json` are not
consumed by the CMake registration, which correctly supplies the configured
suite tolerances (`1e-4` absolute, `1e-5` relative). The observed printed-DOS
deltas are within those configured tolerances, so no fcc-Cu golden or
tolerance change was made.

Debug configuration requested by RF-01:

```bash
cmake -S . -B build-rf-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_CUDA_PLUGIN=OFF \
  -DENABLE_MPI=OFF \
  -DENABLE_OPENMP=OFF \
  -DRUN_UNIT_TESTS=ON \
  -DRUN_REG_TESTS=ON \
  -DRUN_EXAMPLE_TESTS=ON \
  -DCMAKE_Fortran_FLAGS_DEBUG='-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan'
cmake --build build-rf-debug --parallel
```

GNU 13.3.0 and 14.2.0 both expose a GNU recursion-check linker defect when Debug uses `-fcheck=all`: unresolved compiler-generated `is_recursive.*` symbols. The GNU Debug configuration now appends `-fcheck=all,no-recursion`, retaining the other requested checks. A clean GNU 14.2.0 Debug build links all production and unit targets and passes `ctest -L unit` 40/40. The Dyson test temporarily masks only the divide-by-zero emitted internally by oneMKL `zheev`; the missing WP6 Python test files were restored as independent oracles. ifx 2025.3.3 Debug compilation also completes after redundant PUBLIC declarations were removed in the pair-potential and magnetic-tangent modules and the screening-alpha constant was made parent-module scoped. The evidence is recorded in `tests/KNOWN_ISSUES.md`.

## RF-01 numerical contracts and measurements

The existing arbitrary-k test uses exact equality for canonical mesh side-effect checks and folded duplicate coordinates, and `1.0e-12` for eigenvalue and gauge-overlap comparisons. These are FP64, single-process operations with unchanged LAPACK/reduction order. Future batched or parallel comparisons must be tolerance-based and document their reduction order; no bitwise claim is made for them.

The arbitrary-k oracle also checks first- and second-order/SOC assembled-matrix
Hermiticity, eigenpair residual `max|HU-UΛ|`, and metric-free orthonormality
`max|U†U-I|` at `1.0e-12`. The repaired k-space bcc-Fe SCF case explicitly
sets `self.use_kspace=true` and pins 19 values: canonical EF, electron count,
and EBAND; site valence/charge/spin; cumulative DOS state count and selected
DOS samples; and stable namelist moments. Its `[-2,2] Ry` cumulative DOS state
count is 18 as expected. The raw sampled-DOS integral (13.507092) is also
pinned but differs from 18; this pre-existing diagnostic discrepancy is
recorded in `tests/KNOWN_ISSUES.md` rather than treated as a numerical oracle
for the canonical occupations.

The non-gating `UnitTddftCpuProfile` now uses completed, deterministic reciprocal fixtures rather than synthetic eigenpairs: `bccFe_one_site` has `nmat=18`, `(nk,nomega)=(16,96)`, and `fccNi_two_site` has `nmat=36`, `(nk,nomega)=(32,192)`. It separately emits reciprocal Fourier assembly, normal-k eigensolution, arbitrary-`k+q` assembly/eigensolution, and finite-q LMTO pair-potential construction. It also emits TD-DFT vertex construction, denominator generation, response accumulation, Green-function integration, Dyson, and mode-analysis timings.

`PROFILE_MEMORY_MIB` gives analytical resident principal-array payloads (complex FP64 matrices counted as 16 bytes; real FP64 arrays as 8 bytes): H(k), normal and arbitrary-k+q eigenpairs, pair-operator cache plus one Q+ workspace, response tensors, and their sum. Timings are informational process CPU seconds only; CI asserts successful execution, never a performance threshold.

Current GNU 13.3 / oneMKL Release measurement on the recorded i9-11900 host
(2026-08-12):

| Fixture | Reciprocal phases (s) | TD-DFT phases (s) | Principal payload (MiB) |
| --- | --- | --- | ---: |
| bccFe_one_site | Fourier 2.70e-4; k solve 1.2271e-2; k+q 2.36e-4; pair operator 1.3435e-2 | vertices 6.315e-3; denominators 1.2986e-2; accumulation 3.6587e-2; Green integration 6.6798e-2 | 0.33746 |
| fccNi_two_site | Fourier 4.66e-4; k solve 7.482e-3; k+q 1.1739e-2; pair operator 4.207e-2 | vertices 4.038e-2; denominators 1.3727e-1; accumulation 4.2724e-1; Green integration 4.0008e-1 | 3.2249 |

## Stage checklists

### RF-01

- [x] Repository guidance and relevant source/tests read.
- [x] Baseline HEAD and toolchain recorded; user-owned changes preserved.
- [x] Serial Release build and focused reciprocal/TD-DFT baseline passed.
- [x] Full serial Release unit, regression, SCF, and post-processing suites passed (WP8/WP9 disabled by user direction).
- [x] Debug baseline failure isolated and repaired without reference churn.
- [x] Arbitrary-k folding, ordering, side-effect, direct-Fourier, HOH/SOC spectrum, and gauge oracles complete.
- [x] Hermiticity, residual, and orthogonality checks complete.
- [x] K-space SCF observable contract documented and tested.
- [x] TD-DFT scalar/batched equivalence surface complete for chi_KS and direct Xi fixtures.
- [x] CPU profile reports every named phase and principal memory sizes.
- [x] Full post-change serial Release and debug suites pass (Debug units 40/40; Release unit/regression/SCF/post-processing green with WP8/WP9 disabled).
- [x] Verified baseline issue recorded in `tests/KNOWN_ISSUES.md`.
- [x] `git diff --check` passes and RF-01 commit is created.

WP8/WP9 are explicitly disabled at CTest registration by user direction for
RF-01 validation; they remain independently documented executable tests.

### RF-02

- [x] RF-01 commit `28960f1` present and baseline green.
- [x] `parallel_context` has deterministic serial defaults, generic construction,
  type-bound range/root/barrier helpers, local rank/size, and a finalizer that
  frees only its duplicated shared-memory communicator; it never finalizes MPI.
- [x] Legacy MPI globals are initialized at declaration and synchronized only by
  context initialization.
- [x] `ENABLE_MPI=ON, ENABLE_OPENMP=OFF` builds, passes all 45 unit tests, and
  reports one OpenMP thread per MPI process.
- [x] MPI target definitions/linkage use `MPI::MPI_Fortran`; CMake no longer
  changes the Fortran compiler after enabling the language.
- [x] Device assignment is hardware-independent: default
  `mod(local_rank, device_count)`, with a validated programmatic override.
- [x] `report.out`, `minfo.out`, `linfo.out`, `totaldos.out`, and
  `magneticdos.out` are root-owned through open/write/close after collectives.
- [x] Gaussian bcc-Fe k-space SCF is paired against one serial reference at
  serial and MPI ranks 1, 2, and 4; each leaves exactly one complete shared
  output file set.
- [x] Zero-work and oversubscribed ranges are deterministic and collective-safe.
- [x] Serial Release and Debug units each pass 41/41; MPI/OpenMP and
  MPI/no-OpenMP units each pass 45/45.
- [x] `git diff --check` passed; RF-02 progress is updated in this commit.

#### RF-02 validation record

Configured with GNU Fortran 13.3.0, OpenMPI 4, and CUDA disabled. The MPI
wrapper was selected before configuration; the host has an Intel oneAPI MPI
runtime on `LD_LIBRARY_PATH`, so the OpenMPI library directory was put first in
`LIBRARY_PATH` while linking validation binaries to avoid mixing MPI ABIs.

```bash
# Serial/OpenMP Release: 41/41 unit tests
cmake --build build-rf-serial --parallel 4
ctest --test-dir build-rf-serial -L unit -j 4 --output-on-failure

# MPI/OpenMP Release: 45/45 unit tests
cmake -S . -B build-rf-mpi -DCMAKE_Fortran_COMPILER=/bin/mpifort \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA_PLUGIN=OFF \
  -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DMPIEXEC_EXECUTABLE=/usr/bin/mpirun
LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib cmake --build build-rf-mpi --parallel 8
ctest --test-dir build-rf-mpi -L unit -j 4 --output-on-failure

# MPI without OpenMP: 45/45 unit tests
ctest --test-dir build-rf-mpi-noomp -L unit -j 4 --output-on-failure

# Debug build: 41/41 unit tests
cmake --build build-rf-debug --parallel 1
ctest --test-dir build-rf-debug -L unit -j 1 --output-on-failure
```

The bcc-Fe Gaussian SCF registration runs serial plus MPI ranks 1, 2, and 4
against `tests/scf/references/Example_k_space_scf_bccFe`. All four passed with
the CTest tolerances (`abs=1e-4`, `rel=1e-5`) over 18 checked values. Its
canonical values are unchanged by rank count: `EF=-0.0463582057 Ry`,
`N=8`, and `EBAND=-1.8920685 Ry`; Gaussian DOS integrates to 18. The MPI run
also verifies exactly one each of the root-owned report/minfo/linfo and DOS
output files.

### RF-03

- [x] RF-02 is present and green (`41fc840` plus documentation correction `23b853c`).
- [x] Added concrete `kpoint_workset`: local point/weight storage, contiguous ownership maps, finalizer, validation, folding, and non-mutating q shift.
- [x] Defined the empty-rank invariant as allocated zero-length arrays with global range `[1,0]`.
- [x] Reciprocal construction creates a replicated workset; Gaussian k assembly selects distributed local ownership. Fourier assembly and canonical occupations consume local workset arrays.
- [x] TD-DFT obtains k+q from `k_workset%shifted`, retaining the base mesh unchanged.
- [x] Tetrahedron DOS, reciprocal Green, and reciprocal moments reject distributed worksets at their public boundaries.
- [x] Added `UnitKpointWorkset`, covering simulated 1/2/4 rank ownership, empty ranks, map round-trips, weight preservation, folding boundaries/negative/k+G, duplicate order, and q shifts.
- [x] Serial Release build and the complete unit-labelled CTest suite pass (42/42), including `UnitKpointWorkset`, occupations, arbitrary-k, and TD-DFT coverage.
- [x] `git diff --check` passes.
- [x] Open MPI runtime validation passes: the executable resolves only `libmpi.so.40`/`libmpi_mpifh.so.40`, and `Unit(KpointWorkset|KspaceOccupations|ParallelContext)` passes 7/7, including 1/2/4-rank parallel-context and MPI-occupation tests. Validation used `env/openmpi.sh` with `RSLMTO_FC=/usr/bin/mpifort` on the host shell.
- [ ] Full SCF/GBT/DOS/TD-DFT regression matrix and a new performance baseline remain to be rerun; no reciprocal mathematical or storage-layout batching was introduced, so no material performance change is expected.

**Ownership table.** `kpoint_workset` owns the authoritative local points,
weights, and mappings. Reciprocal mesh construction creates it; Fourier and
occupation consumers read it; TD-DFT creates a caller-local shifted workset.
Tetrahedron/Green/moment consumers explicitly require a replicated workset.
Legacy `reciprocal%k_points`, `k_weights`, and mapping fields are transitional
read-only full-list compatibility views during the remaining migration; their
map metadata is asserted against the workset at setup. The arbitrary-k exact
duplicate scan deliberately remains quadratic: its bit-exact comparison has no
safe tolerance-key replacement under the existing coordinate convention.

### RF-04

- [ ] Start only after RF-03 is green and committed.

### RF-05

- [ ] Start only after RF-04 is green and committed.

### RF-06

- [ ] Start only after RF-05 is green and committed.

## Decisions affecting later MPI/GPU work

- The default numerical path remains FP64 CPU with Intel oneMKL LAPACK.
- `ENABLE_CUDA_PLUGIN=OFF` remains mandatory through RF-06 validation.
- No workset, backend, MPI ownership, or production-algorithm change has begun while RF-01 baseline status is blocked.

## Completed-stage commits

| Stage | Commit |
| --- | --- |
| RF-01 | Created (exact SHA reported in the handoff) |
| RF-02 | Created in this commit (exact SHA reported in the handoff) |
| RF-03 | Created in this commit (exact SHA reported in the handoff) |
| RF-04 | Not started |
| RF-05 | Not started |
| RF-06 | Not started |
