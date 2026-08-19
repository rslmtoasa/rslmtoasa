# RS-LMTO-ASA Phase III-A Performance Rescue
## ACC-P0 ... ACC-P4 — Real-material, persistent, batched GPU performance campaign

**Target branch:** `fable_v3`
**Purpose:** determine honestly whether reciprocal CUDA eigensolution can become performance-beneficial after the ACC-00...ACC-09 correctness/infrastructure campaign.

---

# Why this rescue campaign exists

The first ACC campaign established a correct reciprocal CUDA backend but did not establish a performant algorithm.

The current measured GPU path is dominated by two effects that must be separated before making a scientific/engineering decision:

1. small Hermitian matrices are currently passed to cuSOLVER one at a time rather than exposed as a true same-size batch;
2. the ACC-06 crossover campaign includes fresh-process CUDA/runtime/backend startup in every executable invocation.

Therefore the current ~0.4 s GPU floor and 20–100x slowdown for matrices of dimension 8–72 do **not yet answer** the intended question:

> Can persistent GPU execution accelerate either:
>
> - many small real LMTO H(k) matrices; or
> - fewer much larger real LMTO H(k) matrices?

The rescue campaign answers that question using real production Si and Fe Hamiltonians, including periodic bcc-Fe supercells.

---

# Mandatory benchmark plane: matrix dimension x number of k points

Do not study only the primitive-cell "many k points" limit.

For collinear bcc Fe with an spd basis, the ordinary spin Hamiltonian dimension is expected to scale as:

```text
1x1x1 : 1 atom    -> ~18
2x2x2 : 8 atoms   -> ~144
3x3x3 : 27 atoms  -> ~486
4x4x4 : 64 atoms  -> ~1152
5x5x5 : 125 atoms -> ~2250
```

The exact dimension must be read from the production `basis`/Hamiltonian state and asserted rather than hard-coded blindly.

The benchmark campaign should cover both axes:

```text
                    Nk
                    ^
                    |
 many tiny matrices|  Si primitive / Fe primitive
                    |
                    |
                    |
                    +--------------------------> matrix dimension
                                            Fe 2^3 ... Fe 5^3
                                            few larger matrices
```

A GPU may lose badly in the upper-left corner and still win in the lower-right corner.

---

# Physics construction of the large Fe benchmark

The large benchmark must not be a random Hermitian matrix.

Use the current validated production LMTO path.

Recommended construction:

1. start from the current canonical converged collinear bcc-Fe/spd reference state;
2. use the same physical lattice and atom type;
3. build LxLxL periodic supercells for L=2,3,4,5;
4. replicate the same converged primitive-cell Fe potential / potential parameters / collinear moment state over translationally equivalent sites;
5. construct structure constants and H(R) through the normal production code;
6. assemble H(k) through the normal reciprocal Fourier assembler;
7. use `strux_lib` unless the current support matrix explicitly requires another backend for the benchmark;
8. do not run a separate unconverged 125-atom SCF just to create a large eigensystem.

Use an existing frozen/restart/reference-potential path if one exists.

Do not add a new production "benchmark-only Hamiltonian" path.

If no clean frozen-potential route exists, use the smallest existing production route that loads the validated primitive potential and performs only the minimum required setup. Document it.

---

# Supercell band-folding physics oracle

For a periodic LxLxL supercell built by exact repetition of the primitive bcc-Fe state, the supercell Gamma spectrum must reproduce the union of primitive-cell spectra at the primitive reciprocal vectors compatible with the supercell periodicity.

Do this using the actual lattice transformation and the code's reciprocal-coordinate convention.

Do not hard-code Cartesian cubic k points if the production lattice basis differs.

For a diagonal LxLxL repetition, the compatible primitive reciprocal coordinates are conceptually:

```text
(m1/L, m2/L, m3/L)
```

for the appropriate integer ranges, modulo the primitive reciprocal lattice.

Use the actual supercell transformation matrix to generate them robustly.

Compare:

- total state count;
- sorted eigenvalue multiset;
- degeneracy groups;
- spin/channel structure where applicable.

This is a physics/symmetry oracle for the construction, independent of GPU timing.

A large supercell is not accepted as a benchmark fixture until this folding identity is satisfied on CPU.

---

# Two benchmark modes

## Mode A — hardware crossover map

Purpose: isolate dependence on matrix size and number of matrices.

Use real production H(k) matrices but deliberately sweep:

- matrix dimension / supercell size;
- actual unique Nk;
- tile/batch size;
- eigenvalues only;
- eigenvalues + eigenvectors.

The k points must be real production reciprocal points, but the workload does not have to correspond to identical physical k-density across every supercell.

## Mode B — physically matched workload

Purpose: estimate realistic user performance.

Choose a primitive k density, then reduce the supercell k mesh consistently with supercell enlargement so comparable Brillouin-zone sampling density is maintained.

Record the actual unique k-point count after folding/symmetry reduction.

Do not compare nominal mesh dimensions if the production worksets contain different numbers of unique points.

---

# Timing policy

Every benchmark must separate:

```text
T_context_init
T_backend_init
T_operator_prepare
T_Hk_CPU
T_H2D
T_solver
T_D2H
T_post
T_total_steady
```

where applicable.

Cold-start and steady-state timings are both useful but must never be mixed.

Steady-state timing must:

1. create the process once;
2. create the CUDA context/backend once;
3. prepare persistent operator state once;
4. perform warm-up solves;
5. reset interval metrics;
6. perform multiple measured repetitions;
7. report median and spread.

Do not launch a fresh executable for every measured repetition.

Do not compare cumulative CUDA counters with per-call wall times.

---

# Solver strategy

The current one-matrix-at-a-time Zheevd path must remain as a correctness/reference GPU solver.

Add a real same-size batched strategy for the small/mid-size regime.

Candidate:

```text
cusolverDnZheevjBatched
```

The batch strategy must be benchmarked rather than assumed to be faster.

For large matrices such as 1152 and 2250, a single/few dense Zheevd solves may be more appropriate than Jacobi batching.

Therefore represent solver strategy explicitly:

```text
zheevd_serial
zheevj_batched
```

and delay the final dispatch threshold until the measured campaign.

Do not hard-code "all CUDA uses batched Jacobi".

---

# Decision gate

After ACC-P0...P3, classify reciprocal GPU eigensolution by measured workload region:

```text
CPU-preferred
GPU-preferred
inconclusive / unsupported
```

A reasonable performance promotion criterion is:

- >=1.5x steady-state eigensolver speedup on at least one realistic real-material production regime;
- and >=1.2x reciprocal-phase/end-to-end benefit in that regime;
- with full numerical correctness retained.

These are engineering decision gates, not scientific tolerances.

If no realistic regime crosses over after true batching and persistent execution, stop investing in primitive/small-supercell GPU eigensolution and pivot accelerator effort toward higher-arithmetic-intensity workloads such as Lehmann contractions and KPM/Kubo-Bastin transport.

---

# COMMON LUNA HEADER

Prepend this to every ACC-P task:

```text
You are working on the CURRENT HEAD of the fable_v3 branch:

https://github.com/rslmtoasa/rslmtoasa

This task belongs to the Phase III-A reciprocal GPU PERFORMANCE RESCUE campaign.

The ACC-00...ACC-09 campaign established substantial correctness and backend
infrastructure, but it DID NOT establish a performant reciprocal GPU algorithm.

Your job in ACC-P is explicitly performance-first, while preserving all
validated CPU/GPU numerical contracts.

============================================================
NON-NEGOTIABLE RULES
============================================================

1. INSPECT CURRENT CODE BEFORE EDITING

Read the current:
- reciprocal backend Fortran code;
- CUDA reciprocal implementation;
- ACC benchmark harness;
- ACC reports/docs;
- Phase-II/ACC support documentation;
- current CPU/GPU tests.

Do not rely on old task wording if current HEAD differs.

2. REAL LMTO WORKLOADS ARE MANDATORY

Synthetic/deterministic matrices may be retained for microbenchmarks, but
performance conclusions MUST include real production Si and bcc-Fe Hamiltonians.

For large-matrix tests use real periodic Fe supercells constructed through the
normal LMTO -> H(R) -> H(k) production path.

Do not label a hand-filled benchmark matrix "Si" or "Fe".

3. PRESERVE PHYSICS

Do not change:
- LMTO Hamiltonian conventions;
- structure constants;
- spin ordering;
- H(R)->H(k) phase convention;
- occupations;
- SCF physics;
- reference energies

to improve GPU timing.

4. CPU IS THE SCIENTIFIC ORACLE

Before optimizing, record CPU eigenvalues/physical outputs.

GPU correctness requires:
- eigenvalue agreement;
- eigenpair residuals;
- orthogonality;
- projector/subspace comparison for degeneracies;
- downstream physical agreement where applicable.

Do not compare raw eigenvector phases.

5. NO SILENT CPU FALLBACK

If a CUDA solver strategy is unsupported or fails, report it explicitly.

Do not execute CPU LAPACK behind a "CUDA" benchmark result.

6. COLD AND STEADY STATE MUST BE SEPARATE

Do not benchmark GPU performance by spawning a new process for every measured
solve and then call that steady-state performance.

Measure context/backend initialization separately from repeated solves within
one persistent process.

7. METRICS MUST SHARE THE SAME INTERVAL

Never compare cumulative H2D/solver/D2H counters to a one-call wall time.

Reset/snapshot/difference metrics around the same measured region.

8. TRUE BATCHING MEANS ONE BATCH SOLVER CALL

A host loop that invokes one cuSOLVER dense eigensolver per matrix is NOT a
batched GPU algorithm.

When this task requests batching, expose the tile as a true same-size batch to
a batched solver API.

9. LARGE MATRICES MATTER

The campaign must include the real matrix-dimension axis.

For bcc Fe/spd/nsp=2 inspect/verify approximately:

    primitive : 18
    2^3       : 144
    3^3       : 486
    4^3       : 1152
    5^3       : 2250

Do not hard-code these without checking the actual production dimension.

10. MEMORY SAFETY

Before large GPU solves:
- estimate required matrix/eigenvector/workspace storage;
- query available GPU memory where possible;
- avoid OOM;
- skip a case explicitly with recorded reason if it cannot fit.

Never reduce matrix dimension silently.

11. PERFORMANCE CLAIMS ARE END-TO-END

Record kernel/solver speedup AND reciprocal-phase/end-to-end speedup.

A faster cuSOLVER call does not establish application acceleration.

12. KEEP THE LAPACK AND EXISTING ZHEEVD PATHS

The CPU LAPACK backend remains the trusted reference.

The current CUDA Zheevd path remains a GPU reference/fallback while the batched
strategy is evaluated.

Do not delete them merely because another strategy is faster.

13. DO NOT START ACC-10/11/12 OPTIMIZATIONS INSIDE ACC-P

No GPU Lehmann implementation.
No device-resident spectral framework.
No GPU H(k) assembler.

ACC-P determines whether/how the eigensolver itself should continue.

14. ONE FOCUSED COMMIT PER TASK

When done:
- tick every completed checkbox;
- leave genuinely incomplete boxes unticked;
- state exactly what was measured on real hardware;
- list numerical evidence;
- list performance evidence;
- list unsupported/skipped cases and why;
- make one focused commit with the supplied single-line message.
```

---

# ACC-P0 — Repair the benchmark methodology and build the real-material corpus

## Objective

Replace the misleading cold-process/synthetic-only interpretation of ACC-06
with a benchmark harness that measures:

1. cold CUDA startup;
2. persistent steady-state backend performance;
3. real Si/Fe production matrices;
4. the large-matrix Fe-supercell axis.

Do NOT optimize the CUDA solver in this task.

This task establishes trustworthy evidence before algorithm changes.

## A. Audit current benchmark implementation

Inspect the current ACC-06 benchmark driver and document exactly:

- where the process starts/stops;
- where `make_execution_backend` is called;
- whether CUDA context creation is inside the timed region;
- whether cuSOLVER handle/stream creation is inside the timed region;
- whether each repetition launches a fresh executable;
- whether timing counters are cumulative or interval-local;
- whether matrices are synthetic or come from production H(R).

Add these findings to the benchmark documentation.

Do not proceed until the old ~0.4 s floor is decomposed.

## B. Add persistent-process benchmark mode

Add a mode that performs, within ONE executable invocation:

```text
initialize physical system
build reciprocal object
create selected backend once
prepare operator once

warmup solve 1
warmup solve 2

reset/snapshot metrics

repeat measured solve N times:
    assemble/request same class of workload
    solve
    synchronize at required API boundary only

report median/min/spread
report metric deltas
```

Recommended measured repetitions:
- at least 5 for fast small cases;
- at least 3 for expensive 4^3/5^3 cases.

Do not force five repetitions of a very expensive 2250x2250 solve if runtime is
unreasonable; document the chosen count.

Report separately:

```text
cold_process_wall
cuda_context_backend_init
first_solve
steady_solve_median
```

## C. Make timing counters interval-consistent

Inspect CUDA backend metrics.

If current H2D/solver/D2H counters are cumulative, add one of:

- explicit reset before benchmark region; or
- snapshot-before / snapshot-after and difference.

Do not destroy useful lifetime counters if production diagnostics rely on them.

The benchmark report must never combine lifetime GPU counters with one-call wall time.

## D. Real diamond-Si primitive fixture

Locate the CURRENT canonical validated Si/sp reciprocal fixture.

Use it directly or construct a benchmark wrapper around the same production
input.

Requirements:
- real Si lattice;
- real sp basis;
- real production structure constants;
- real H(R);
- real H(k) Fourier assembly;
- normal reciprocal backend.

Record:
- atom count;
- basis dimension;
- actual Hamiltonian dimension;
- actual unique Nk.

Do not fill a deterministic matrix by hand and call it Si.

## E. Real bcc-Fe primitive fixture

Likewise use the current validated collinear bcc-Fe/spd reciprocal reference.

Prefer a converged/frozen reference state so benchmark repetitions do not
measure unrelated SCF drift.

Record the actual Hamiltonian dimension.

## F. Build real bcc-Fe L^3 supercell benchmark fixtures

Construct:
- 2x2x2;
- 3x3x3;
- 4x4x4;
- 5x5x5

periodic bcc-Fe supercells from the canonical primitive reference.

Use:
- identical Fe type;
- identical primitive converged potential/potential parameters;
- identical collinear spin direction;
- same Hamiltonian options;
- production structure constants;
- production H(R);
- production reciprocal H(k).

Do NOT run separate unconverged large-supercell physics merely to obtain a
matrix.

If an existing restart/frozen-potential route exists, use it.

If the repository has a validated supercell generator/helper, reuse it.

Do not create a benchmark-only alternative Hamiltonian builder.

## G. Validate every supercell by band folding BEFORE timing

For each L:

1. generate the primitive reciprocal k set compatible with the supercell
   translation matrix;
2. solve those primitive H(k) matrices using CPU LAPACK;
3. solve supercell Gamma using CPU LAPACK;
4. compare the sorted eigenvalue multiset.

Check:
- number of primitive folded states equals supercell matrix dimension;
- eigenvalues agree within a tight production numerical tolerance;
- degeneracies are handled as eigenvalue multisets/subspaces;
- no site-dependent moment/potential drift entered the construction.

Do not use GPU results to validate the fixture.

If folding fails:
STOP.
The large benchmark fixture is not valid.
Fix/understand the construction before performance work.

## H. Benchmark modes

Implement two selectable modes:

### crossover mode
Sweep actual production k points over useful Nk values independent of strict
equal physical k density.

### matched-density mode
Choose comparable BZ sampling density between primitive and supercell cases.

Record both nominal mesh and actual unique Nk after production workset creation.

## I. Suggested initial campaign grid

Primitive Si/Fe:
- actual Nk near 1;
- ~32/64;
- ~128/216;
- ~512;
- ~1000 if affordable.

Use real mesh/workset generation rather than fabricating duplicate k points.

Fe 2^3:
- Nk = 1, ~8, ~32, and one matched-density case.

Fe 3^3:
- Nk = 1, ~8, and matched-density.

Fe 4^3:
- Nk = 1 and one small/matched-density case.

Fe 5^3:
- Gamma first;
- add more k points only if memory/runtime permit.

For every case test eigenvalues-only and eigenvectors where production supports both.

Do not require all combinations if memory/runtime is prohibitive; record skips.

## Deliverables

Produce a machine-readable table with:

```text
fixture
source physical fixture
L
Natom
nmat
nominal k mesh
actual unique Nk
tile
vectors yes/no
backend
cold init
first solve
steady median
Hk CPU
H2D
solver
D2H
total
memory estimate
```

## Checklist

- [x] Old ACC-06 timing interval audited
- [x] CUDA startup isolated
- [x] persistent-process benchmark implemented
- [x] metric reset/delta implemented
- [x] real Si production fixture benchmarked
- [x] real primitive Fe production fixture benchmarked
- [x] Fe 2^3 production supercell constructed
- [x] Fe 3^3 production supercell constructed
- [x] Fe 4^3 production supercell constructed
- [x] Fe 5^3 production supercell constructed
- [x] actual matrix dimensions recorded
- [x] every supercell CPU band-folding identity validated
- [x] crossover mode implemented
- [x] matched-density mode implemented
- [x] large-case memory preflight implemented
- [x] cold and steady timings reported separately
- [x] no CUDA solver algorithm changed
- [x] benchmark methodology documented honestly

ACC-P0 evidence from the CUDA-enabled host: the three-repetition Gamma
confirmation measured bcc-Fe steady GPU/CPU speedups of 0.05x (L=1, nmat=18),
0.24x (L=2, nmat=144), 0.71x (L=3, nmat=486), 1.55x (L=4, nmat=1152), and
2.68x (L=5, nmat=2250). The L=2,3,4,5 CPU folding checks passed with maximum
sorted-eigenvalue errors of 1.55e-13, 5.92e-13, 5.16e-13, and 4.28e-13,
respectively, with matching degeneracy groups. The quick Si/Fe campaign also
covered values-only and eigenvector requests, crossover and matched-density
labels, tiles 1/8, and nominal 1x1x1 and 2x2x2 meshes. GPU preflight reported
approximately 15.7 GiB free on the selected RTX A4000; no case was skipped for
memory. Raw JSON/CSV evidence was kept machine-local under `/tmp`.

The full suggested Nk grid was intentionally not forced for L=4/5: one dense
2250x2250 solve is already materially more expensive than the small cases.
The campaign driver supports the broader grid and records the chosen
repetition policy, but the committed evidence uses the safe Gamma-heavy
campaign described above.

**Commit message:** `Build persistent real-material GPU benchmarks`

---

# ACC-P1 — Implement and validate a true batched small/mid-size CUDA eigensolver

## Objective

Replace the current host loop over individual Zheevd calls as the only
small-matrix CUDA strategy with a TRUE same-size batched eigensolver strategy.

Keep the current Zheevd implementation as a reference/fallback.

The goal is to expose independent H(k) matrices concurrently to the GPU.

## A. Confirm current behavior

Before editing, inspect the current C++/CUDA solver and verify whether it still
contains conceptually:

```cpp
for (ibatch = 0; ibatch < batch_size; ++ibatch) {
    cusolverDnZheevd(... matrix ibatch ...);
}
```

Record this as the pre-change GPU strategy.

If current HEAD already contains real batching, do not duplicate it; audit what
it actually does and adapt the task.

## B. Add explicit internal solver strategies

Introduce a narrow internal strategy selection such as:

```text
zheevd_serial
zheevj_batched
```

Do not expose a broad new public solver framework.

For benchmark/debug selection, a small backend option/environment/test control
is acceptable if it does not pollute normal scientific input.

Do not set final automatic thresholds in this task.

## C. Implement `cusolverDnZheevjBatched` correctly

For a same-size tile:

- matrices must be contiguous;
- each matrix is column-major;
- eigenvalues are contiguous per matrix;
- one `info` value must exist per matrix;
- eigenvectors overwrite the device matrix buffer when requested.

Verify that the existing Fortran/C interoperability layout maps correctly.

Do not transpose/copy merely by assumption.

Add a layout test with known matrices.

Create and retain the Jacobi parameter object for backend/context lifetime
rather than every solve.

Query required workspace using the proper batched buffer-size call.

Reuse/grow workspace across solves.

Allocate `device_info` for `batch_size`, not one integer.

Copy/check every `info[i]`.

If any matrix does not converge:
- fail the CUDA request explicitly;
- report which batch matrix failed;
- do not silently rerun it on CPU.

## D. Accuracy settings

Start from cuSOLVER's documented/default Jacobi parameters.

Then measure production residuals.

Only change tolerance/max-sweeps if the default fails the same numerical
contract used by CPU/GPU validation.

Document any non-default settings.

Do not loosen the RS-LMTO validation tolerance to accommodate a weak GPU solve.

## E. Numerical validation

For real Si and primitive Fe tiles, and at least one Fe 2^3 / 3^3 case if the
batched solver can handle it reasonably:

Compare against CPU LAPACK and existing CUDA Zheevd:

1. eigenvalues;
2. residual:
       ||H v - e v|| / max(||H||,1)
3. orthogonality:
       ||V^H V - I||
4. degenerate-subspace projector comparison;
5. per-matrix `info`.

Because `ZheevjBatched` does not provide a usable batched residual report,
compute the RS-LMTO residual explicitly from the original H and returned
eigenvectors.

Keep an untouched host copy of H for validation.

Do not judge correctness from solver `info=0` alone.

## F. Performance comparison

Using the persistent ACC-P0 harness compare:

```text
CPU LAPACK
CUDA zheevd_serial
CUDA zheevj_batched
```

on REAL matrices.

Mandatory regions:

### many small matrices
Si primitive and Fe primitive:
- Nk increasing to several hundred or more.

### medium matrices
Fe 2^3 (~144)
Fe 3^3 (~486)

Try meaningful batch sizes constrained by memory.

### large matrices
Fe 4^3 / 5^3:
do NOT force Jacobi batching if memory/API behavior is unsuitable.

The existing Zheevd path is the expected candidate there.

Record where `zheevj_batched` ceases to be beneficial or feasible.

## G. No premature dispatch

Do not add:

```text
if n < X use Jacobi else use Zheevd
```

to production yet.

ACC-P3 will establish the threshold from the complete campaign.

Keep strategy explicitly selectable for measurement.

## Checklist

- [ ] Existing serial GPU solver behavior confirmed
- [ ] narrow solver-strategy selection added
- [ ] true `ZheevjBatched` path implemented
- [ ] matrices passed as one contiguous batch
- [ ] batched workspace queried/reused
- [ ] Jacobi parameter object reused
- [ ] per-matrix device info implemented
- [ ] no silent CPU fallback exists
- [ ] Si real-matrix residuals pass
- [ ] primitive Fe residuals pass
- [ ] Fe 2^3 tested
- [ ] Fe 3^3 tested where feasible
- [ ] eigenvalue agreement passes
- [ ] orthogonality passes
- [ ] degenerate subspaces compared correctly
- [ ] CPU/Zheevd/Zheevj performance compared
- [ ] no final automatic threshold hard-coded
- [ ] benchmark report records non-convergence/unsupported cases honestly

**Commit message:** `Add true batched CUDA Hermitian eigensolver`

---

# ACC-P2 — Remove hot-path overhead and make GPU timing trustworthy

## Objective

After true batching exists, remove avoidable overhead that distorts the
small-matrix path and improve staging for larger matrices ONLY where measured.

This is a profiling-driven optimization task.

Do not alter the LMTO/Fourier physics.

## A. Reuse timing events

Inspect whether CUDA events are currently created/destroyed for every tile.

If so:

- create reusable timing events in backend/context initialization; or
- disable fine-grained event timing outside benchmark/profile mode.

Do not pay repeated event-construction overhead in the production hot path.

Preserve enough instrumentation for ACC-P3.

## B. Remove redundant synchronization

Map the synchronization boundary exactly.

If the C++ solver already performs `cudaStreamSynchronize()` before returning
and the Fortran caller immediately synchronizes the same backend again:

- choose ONE required synchronous API boundary;
- remove the redundant synchronization.

Do not remove synchronization needed before host output is consumed.

Validate correctness under compute-sanitizer/debug tooling if available.

## C. Reuse staging buffers

Profile current host allocation/copy behavior for:

- request-owned H(k) tile;
- H2D staging;
- eigenvalue/eigenvector D2H staging.

If repeated Fortran allocation/copy is material:
- reuse an existing reciprocal workspace or backend staging allocation;
- preserve clear ownership;
- avoid introducing a global memory pool.

Do not optimize a copy that is below timing noise.

## D. Evaluate pinned host memory for LARGE transfers

This is specifically relevant to Fe 4^3/5^3, where a dense complex matrix is
large.

Benchmark pageable versus pinned host staging for:
- H2D matrix;
- D2H eigenvectors.

Only retain pinned allocations if they yield a measurable end-to-end benefit.

Pinned memory must:
- have explicit ownership/lifecycle;
- not leak;
- not be allocated per matrix.

Do not pin large application arrays globally.

## E. Persistent workspace policy

For each solver strategy:

- query required workspace for `(n,batch,vectors)` state;
- reuse if capacity is sufficient;
- grow only when needed;
- do not cudaMalloc/cudaFree inside the per-matrix loop;
- do not shrink repeatedly.

Record allocation count in benchmark diagnostics if practical.

## F. Metric semantics

Define metrics clearly:

```text
lifetime_*
last_request_*
```

or provide snapshot/delta functionality.

Do not leave ambiguous counters.

The ACC-P3 report must be able to state a mathematically consistent relation
between component and total time.

## G. Rebenchmark after each optimization

For each accepted change, show before/after on at least:

- primitive Fe many-k;
- Fe 3^3;
- Fe 5^3 Gamma if feasible.

Do not batch many unrelated micro-optimizations into one unexplained speedup.

## Checklist

- [ ] event creation overhead measured
- [ ] timing events reused or profiling-gated
- [ ] synchronization path mapped
- [ ] redundant synchronization removed if confirmed
- [ ] host request/staging copies profiled
- [ ] staging reused only where material
- [ ] pinned memory benchmarked for large cases
- [ ] pinned memory retained only if beneficial
- [ ] solver workspace persistent
- [ ] hot-path cudaMalloc/cudaFree removed
- [ ] metric lifetime/interval semantics explicit
- [ ] primitive Fe benchmark improved/reported
- [ ] Fe 3^3 benchmark improved/reported
- [ ] Fe 5^3 benchmark improved/reported where feasible
- [ ] no physics or CPU references changed

**Commit message:** `Reduce reciprocal CUDA hot-path overhead`

---

# ACC-P3 — Run the decisive real-material 2D crossover campaign

## Objective

Run the performance campaign that ACC-06 was intended to provide:

> determine GPU crossover over BOTH number of k points and actual LMTO matrix
> dimension, using persistent execution and real Si/Fe Hamiltonians.

This task should make minimal production-code changes.

Its main deliverable is evidence.

## A. Freeze code before campaign

Before timing:

- run lean CPU correctness gate;
- run reciprocal CPU/GPU correctness fixtures;
- run band-folding checks for every Fe supercell;
- record git commit;
- record compiler, MKL, CUDA, driver, GPU, CPU, thread count.

Do not optimize while collecting the final table.

## B. CPU reference configurations

Benchmark a reasonable CPU reference.

At minimum:

1. production LAPACK backend with the normal threading configuration;
2. best established benchmark-only k-parallel CPU strategy if ACC-06 showed it
   can beat production for small matrices.

Control MKL/OpenMP thread counts explicitly.

Avoid oversubscription.

Report which CPU strategy wins each row.

GPU speedup must be computed against the best reasonable CPU row.

## C. Mandatory real fixtures

### Si primitive/sp
Use several real k meshes up to the practical high-hundreds/~1000 unique-k
regime.

### bcc Fe primitive/spd
Same.

### bcc Fe 2^3
nmat expected near 144.

### bcc Fe 3^3
nmat expected near 486.

### bcc Fe 4^3
nmat expected near 1152.

### bcc Fe 5^3
nmat expected near 2250.

Record actual dimensions from runtime state.

Do not substitute synthetic matrices for any mandatory row.

## D. Mandatory solver strategies

For rows where supported, compare:

- CPU best;
- CUDA Zheevd serial;
- CUDA ZheevjBatched.

Do not omit the old Zheevd GPU path merely because Jacobi is faster.

For large matrices, record Jacobi as unsupported/not competitive if that is the
measured result.

## E. Two campaign slices

### Slice 1 — crossover plane

Choose sufficient `(nmat, Nk)` points to expose the boundary.

Suggested shape, adjusted for actual memory/runtime:

```text
primitive: Nk ~ 1, 32/64, 128/216, 512, ~1000
2^3:       Nk ~ 1, 8, 32, 128 if feasible
3^3:       Nk ~ 1, 8, 32 if feasible
4^3:       Nk ~ 1, 4/8 if feasible
5^3:       Nk = 1 mandatory; >1 only if feasible
```

### Slice 2 — matched physical sampling

Pick at least one physically comparable k-density sequence across primitive and
supercells.

Record nominal and actual unique k counts.

## F. Eigenvalues-only versus eigenvectors

For every feasible row compare both modes.

This matters because:
- SCF normally needs eigenvectors;
- bands/spectral reconnaissance may sometimes need only eigenvalues.

Do not claim a speedup for SCF based only on eigenvalues-only results.

## G. Cold versus steady state

For every fixture class report:

- backend/context init;
- first solve;
- steady-state median.

For application recommendations use steady-state for repeated SCF-style work,
but retain cold cost for short one-shot workflows.

## H. End-to-end production workflows

In addition to isolated backend timing run actual production calculations:

Mandatory:
- real Si primitive reciprocal workflow;
- real Fe primitive reciprocal workflow.

At least one larger end-to-end case:
- Fe 2^3 or Fe 3^3, chosen by practical runtime.

For Fe 4^3/5^3, a validated frozen-potential eigensystem benchmark is sufficient
if full SCF is unreasonably expensive.

Report:

```text
solver speedup
reciprocal-phase speedup
whole executable speedup
```

Do not infer whole-code speedup from backend-only timing.

## I. Memory report

For large cases record:
- GPU free memory before solve;
- allocated matrix bytes;
- eigenvector bytes;
- workspace bytes;
- peak approximate device usage.

If Fe 5^3 cannot fit:
- record exact reason;
- do not silently lower L;
- the skip is valid evidence.

## J. Required output

Produce a 2D crossover/support map.

At minimum:

```text
fixture | L | nmat | Nk | vectors | CPU strategy | CPU steady |
GPU strategy | GPU steady | solver speedup | reciprocal speedup |
whole-run speedup | memory | recommendation
```

Use categories:

```text
CPU-preferred
GPU-preferred
no clear benefit
unsupported
```

## Checklist

- [ ] correctness gate passed before campaign
- [ ] environment frozen/recorded
- [ ] best CPU baseline used
- [ ] real Si campaign complete
- [ ] real Fe primitive campaign complete
- [ ] Fe 2^3 campaign complete
- [ ] Fe 3^3 campaign complete
- [ ] Fe 4^3 campaign complete
- [ ] Fe 5^3 Gamma measured or explicitly memory/runtime-skipped
- [ ] actual nmat recorded for every row
- [ ] actual unique Nk recorded
- [ ] Zheevd GPU retained in comparison
- [ ] ZheevjBatched compared where supported
- [ ] eigenvalues-only measured
- [ ] eigenvectors measured
- [ ] cold startup separated
- [ ] steady-state separated
- [ ] matched-density slice completed
- [ ] actual Si production workflow timed
- [ ] actual Fe production workflow timed
- [ ] at least one larger production/frozen-potential workflow timed
- [ ] GPU memory reported for large cases
- [ ] 2D crossover map produced
- [ ] no result omitted because it was unfavorable to GPU

**Commit message:** `Measure real-material reciprocal GPU crossover`

---

# ACC-P4 — Make the hard reciprocal-GPU go/no-go decision and revise ACC-10...14

## Objective

Use ACC-P3 evidence to decide where reciprocal GPU eigensolution should be a
real supported performance path.

This task must be willing to conclude that some or all matrix-size regimes
should remain CPU-only.

Do not optimize further merely to avoid an unfavorable conclusion.

## A. Classify workload regions

Build a dispatch/support table in terms of:

- nmat;
- Nk/tile;
- eigenvalues-only vs eigenvectors;
- solver strategy;
- cold vs repeated workload.

Example conceptual output:

```text
nmat <= 72, Nk small:
    CPU

nmat <= 72, Nk >= ...:
    batched GPU if measured crossover exists

nmat 144-486:
    measured policy

nmat >= 1152:
    Zheevd GPU / CPU according to measured result
```

Do not use these example boundaries without ACC-P3 evidence.

## B. Performance promotion gate

A workload region may be called `GPU-preferred` only if:

1. numerical correctness is fully retained;
2. steady-state eigensolver speedup is at least about 1.5x;
3. reciprocal-phase or meaningful end-to-end speedup is at least about 1.2x;
4. the result is from real LMTO matrices;
5. it is reproducible across repeated runs on the campaign hardware.

If a different threshold is scientifically/engineering-wise justified, explain
it explicitly rather than quietly changing the gate.

Do not promote a path because a microkernel alone is faster.

## C. Decide backend-selection policy

Choose one:

### Policy 1 — explicit only
Keep LAPACK default; user explicitly selects CUDA.

### Policy 2 — optional auto
Provide a conservative `auto` mode based on measured matrix/batch thresholds,
while retaining explicit LAPACK/CUDA choices.

Do NOT make CUDA the global default based on one RTX A4000 campaign.

If `auto` is implemented:
- unknown hardware must default conservatively;
- threshold logic must be documented;
- no hidden CPU fallback inside explicit CUDA mode.

## D. Decide the fate of small-matrix reciprocal CUDA

If true batching still loses for real primitive Si/Fe:

state explicitly:

```text
primitive-cell small-matrix reciprocal eigensolution remains CPU-preferred
```

Do not keep optimizing it indefinitely.

The CUDA backend may still be retained for:
- correctness/reference;
- larger matrices;
- future consumers;
- other hardware.

## E. Decide the fate of large-matrix reciprocal CUDA

The Fe 3^3/4^3/5^3 results are critical.

If GPU crosses over there:
- define that as the primary reciprocal CUDA production regime;
- document matrix/memory limits;
- keep large-matrix solver strategy explicit.

If it does not cross over:
- state that dense reciprocal eigensolution is not presently a useful RTX A4000
  acceleration target.

## F. Revise ACC-10...ACC-14 plan

### If eigensolver GPU has a real winning regime

Proceed to ACC-10 GPU Lehmann.

ACC-11 device-resident handoff becomes meaningful only if:
- eigensolver and Lehmann both execute on GPU;
- transfer cost is material.

ACC-12 H(k) assembly remains profile-gated.

### If eigensolver has no useful winning regime

Do NOT abandon GPU acceleration generally.

Revise the plan:

- ACC-10 should accelerate Lehmann contractions starting from CPU eigenpairs;
- ACC-11 device-resident eigensolver handoff should be skipped/deferred;
- ACC-12 H(k) GPU assembly should normally be skipped;
- ACC-13 RS KPM/Kubo-Bastin becomes the highest-priority accelerator path.

This is a valid successful outcome.

## G. Document Luna task status honestly

Update the ACC report to distinguish:

```text
correctness fulfilled
infrastructure fulfilled
performance fulfilled
performance incomplete
```

for ACC-00...ACC-09.

In particular do not retroactively describe ACC-03/ACC-06 as performance-complete
unless ACC-P evidence justifies it.

## H. No more wishful checkboxes

A checkbox requiring demonstrated speedup is ticked only if the measured
real-material gate passed.

If it did not, leave it unticked and record the negative result.

## Checklist

- [ ] workload regions classified
- [ ] real-material speedup gate applied
- [ ] small-matrix policy decided
- [ ] medium-matrix policy decided
- [ ] large-matrix policy decided
- [ ] explicit/auto backend policy decided
- [ ] CUDA not made global default from one machine
- [ ] no hidden explicit-CUDA CPU fallback introduced
- [ ] ACC-00...09 completion status revised honestly
- [ ] ACC-10 plan revised from evidence
- [ ] ACC-11 go/defer decision made
- [ ] ACC-12 go/defer decision made
- [ ] ACC-13 priority revised if necessary
- [ ] negative performance results documented, not hidden
- [ ] final crossover/support map committed

**Commit message:** `Define reciprocal GPU performance policy`
