# SCF-B0C — Consolidate and complete the reciprocal/SCF CPU/GPU benchmark harness

You are working on the **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task follows the completed reciprocal accelerator work from:

```text
ACC-00 ... ACC-09
ACC-P0 ... ACC-P4
ACC-10 ... ACC-14 where applicable
```

and should also reuse lessons from the now-completed KPM benchmark consolidation:

```text
KPM-B0C
KPM-B1
```

Do **not** restart the ACC optimization campaign.

Do **not** reimplement the reciprocal backend.

Do **not** optimize kernels unless this task uncovers a benchmark-integrity defect
that prevents fair measurement.

The purpose of SCF-B0C is to turn the existing reciprocal/eigensolver benchmark
infrastructure into one coherent, strict, reusable benchmark harness for a final
SCF-B1 campaign.

---

# 0. Current architectural context

RS-LMTO-ASA is primarily a real-space LMTO-ASA code, but reciprocal-space
workflows are used for k-space Hamiltonians, diagonalization, SCF, spectral
methods, and TD-DFT-related functionality.

The existing reciprocal accelerator campaign already established much of the
hard infrastructure:

```text
real Si and bcc-Fe production fixtures
real H(R) -> H(k) construction
primitive and Fe supercell scaling
persistent CUDA execution
cold-start versus steady-state timing
CPU LAPACK reference
CUDA Zheevd
CUDA batched Jacobi where supported
eigenvalues-only versus eigenvectors
best-CPU comparison
real matrix-size x Nk crossover studies
GPU memory accounting
eigenvalue/eigenvector correctness checks
band-folding checks for Fe supercells
large-matrix CUDA crossover evidence
```

The benchmark plane already covered conceptually:

```text
                  number of k points
                         ^
                         |
many small matrices     |  Si primitive / Fe primitive
                         |
                         |
                         +----------------------------> matrix dimension
                                             Fe 2^3 ... Fe 5^3
```

with approximate Fe/spd/nsp=2 dimensions:

```text
1x1x1 : ~18
2x2x2 : ~144
3x3x3 : ~486
4x4x4 : ~1152
5x5x5 : ~2250
```

Actual dimensions must always be read from runtime state.

The missing piece is not a new eigensolver campaign.

The missing piece is a **single end-to-end SCF benchmark contract** that connects:

```text
eigensolver
    ->
reciprocal phase
    ->
SCF iteration
    ->
SCF convergence
```

and produces fair, reproducible CPU/GPU evidence.

---

# 1. Objective

At completion, there must be one canonical SCF benchmark workflow that can be
used directly by a later SCF-B1 frozen campaign.

It must answer, for a given real material/workload:

1. how long H(R)->H(k) assembly takes;
2. how long diagonalization takes;
3. how much time is spent in H2D/D2H transfers;
4. how much time is spent in occupations / Fermi-level handling;
5. how much time is spent reconstructing charge/spin density;
6. how much time is spent updating/mixing the potential;
7. how long one SCF iteration takes;
8. how many iterations are required to converge;
9. how long the complete SCF run takes;
10. whether CPU and GPU converge to the same physical result.

The benchmark must preserve the already-established reciprocal microbenchmark
capability, but add the missing SCF-level accounting.

---

# 2. Non-negotiable rules

## 2.1 Inspect CURRENT HEAD first

Read at minimum:

```text
current reciprocal execution backend
current CUDA reciprocal backend
current H(R)->H(k) assembly path
current SCF driver / SCF loop
current diagonalization interface
current density / occupation / Fermi-level code
current mixing / potential update path
existing ACC benchmark harness
existing ACC-P reports
existing reciprocal correctness tests
existing KPM B0C/B1 harness code and reports
```

Do not assume historical task text still matches CURRENT HEAD.

## 2.2 Preserve physics

Do not change:

```text
LMTO Hamiltonian conventions
structure constants
spin ordering
H(R)->H(k) phase convention
occupation scheme
Fermi-level algorithm
density construction
mixing algorithm
SCF convergence criterion
reference energies
potential update equations
SOC treatment
```

to improve benchmark numbers.

## 2.3 CPU remains the scientific oracle

A GPU row is not a valid performance result unless the corresponding physical
outputs agree with the CPU reference within the established scientific
tolerances.

## 2.4 No silent fallback

Explicit CUDA mode must never execute CPU LAPACK behind the scenes and still be
reported as a GPU row.

Unsupported routes must be marked:

```text
UNSUPPORTED
```

or:

```text
SKIPPED
```

with reason.

## 2.5 No benchmark-only physics path

Use real production Hamiltonians and SCF states.

Do not introduce synthetic SCF states merely to make timing easier.

## 2.6 B0C is infrastructure, not optimization

If a stage is slow, measure and report it.

Do not start a new optimization WP inside SCF-B0C unless the harness itself is
incorrect or unfair.

---

# 3. SCF-B0C-A — Audit and consolidate the existing ACC harness

Before editing, produce a short gap table:

```text
requirement | existing | partial | missing | action
```

Cover at least:

```text
real Si fixture
real bcc-Fe fixture
Fe supercell fixtures
persistent process
cold timing
steady timing
best CPU baseline
OMP/BLAS control
eigenvalues-only
eigenvectors
Zheevd
batched Jacobi
matrix dimension
actual Nk
GPU memory
eigensolver correctness
band-folding correctness
H(k) timing
reciprocal-phase timing
whole-executable timing
SCF-iteration timing
SCF-to-convergence timing
physical convergence comparison
JSON
CSV
Markdown
raw logs
strict pairing
headline eligibility
```

Do not duplicate mechanisms that already exist.

Prefer extending the strongest existing benchmark driver.

---

# 4. SCF-B0C-B — Define benchmark levels explicitly

Every benchmark row must have:

```text
benchmark_level
```

with one of:

```text
eigensolver
reciprocal_phase
scf_iteration
scf_convergence
```

## eigensolver

Measures the diagonalization backend in isolation using real H(k).

## reciprocal_phase

Measures:

```text
H(k) assembly
+
diagonalization
+
required reciprocal postprocessing
```

for one production reciprocal request.

## scf_iteration

Measures one complete SCF iteration from the established iteration boundary to
the next.

## scf_convergence

Measures a complete SCF run from a defined initial/restart state to the normal
production convergence criterion.

Do not mix these levels in one speedup number.

---

# 5. SCF-B0C-C — Instrument the complete SCF iteration

Add exclusive top-level profiling around the SCF loop.

Use names consistent with current code style, but provide the equivalent of:

```text
P_scf_iteration_total

P_hamiltonian_prepare
P_hk_assembly
P_eigensolver
P_eigenpair_transfer
P_occupations_fermi
P_density_build
P_charge_spin_accumulate
P_potential_update
P_mixing
P_scf_io
P_scf_misc
```

If some stages are already profiled elsewhere, reuse them.

Do not double-count nested timers.

Require:

```text
sum(exclusive P_*) ~= P_scf_iteration_total
```

with:

```text
profile_closure_error <= 0.03
```

for rows used in performance conclusions.

Any:

```text
P_scf_misc > 5%
```

must be decomposed or explained.

---

# 6. SCF-B0C-D — Preserve reciprocal detail timers

Keep or expose the existing reciprocal timing decomposition:

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

These are nested/detail timers relative to SCF iteration timing.

Do not sum inclusive child timers into the SCF parent total.

The final harness must support both:

```text
SCF-level view
```

and:

```text
reciprocal-backend view
```

without ambiguity.

---

# 7. SCF-B0C-E — Cold start versus repeated SCF execution

Keep the existing ACC-P distinction:

```text
cold_process_wall
cuda_context_backend_init
first_reciprocal_solve
steady_reciprocal_median
```

and add:

```text
first_scf_iteration
steady_scf_iteration_median
```

for repeated iterations in a real SCF run.

Do not count CUDA context creation as a repeated SCF iteration cost.

For complete SCF convergence also report:

```text
full_scf_wall
```

from process start to convergence.

Therefore the harness should allow a reader to distinguish:

```text
one-shot overhead
steady iteration cost
total convergence cost
```

---

# 8. SCF-B0C-F — Define strict physical pairing

Create one canonical SCF comparison key.

At minimum include:

```text
material / fixture
supercell / replication
Natom
Hamiltonian dimension
nsp
SOC state
basis / lmax information
structure-constant backend
k-mesh definition
actual unique Nk
smearing / occupation parameters
electron count
Fermi-level policy
mixing method
mixing parameters
SCF convergence thresholds
starting-state/restart identity
potential / reference-state identity
eigenvectors required
numeric mode
```

Where appropriate also include:

```text
GBT
CCOR
Hubbard
overlap/generalized eigenproblem
other feature flags
```

Do not compare two SCF rows that differ in physics and call the result a backend
speedup.

Generate a stable row/pairing fingerprint.

---

# 9. SCF-B0C-G — Define numerical modes

Every row must record independently:

```text
Hamiltonian precision
eigensolver precision
eigenvector precision
density-accumulation precision
SCF canonical precision
```

Then classify:

```text
numeric_mode = fp64 / fp32 / mixed
```

according to the actual route.

Do not infer "FP32" merely because the CUDA eigensolver is FP32 if the rest of
the SCF iteration remains FP64.

If a genuine end-to-end FP32 SCF path does not exist, state that clearly.

It is acceptable for SCF-B1 to compare:

```text
FP64 CPU vs FP64 GPU
```

and separately:

```text
mixed GPU production route
```

if that is the real implementation.

Do not manufacture a full-FP32 SCF path solely for benchmark symmetry.

---

# 10. SCF-B0C-H — CPU fairness contract

For every headline SCF workload run CPU with:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

unless the current host makes one configuration invalid.

Control:

```text
MKL_NUM_THREADS
OPENBLAS_NUM_THREADS
OMP_PROC_BIND
OMP_PLACES
```

and avoid oversubscription.

For each physics/numeric-mode group define:

```text
best_cpu_same_physics
=
minimum median valid SCF iteration time
```

and separately:

```text
best_cpu_full_scf_wall
```

for full convergence campaigns.

Retain all thread rows.

---

# 11. SCF-B0C-I — GPU strategy metadata

Record explicitly:

```text
gpu_solver_strategy
```

for example:

```text
zheevd
zheevj_batched
other current supported route
```

Also record:

```text
batch/tile size
eigenvalues_only or eigenvectors
device index
workspace bytes
matrix bytes
eigenvector bytes
peak approximate GPU bytes
```

For SCF headline rows:

```text
eigenvectors = required
```

unless CURRENT HEAD has a scientifically valid SCF route not requiring them.

Do not claim SCF acceleration from eigenvalues-only measurements.

---

# 12. SCF-B0C-J — Correctness contract

Every headline SCF pair must carry structured correctness evidence.

At minimum compare:

## Eigensystem

```text
eigenvalues
residual norms
orthogonality
subspace/projector comparison for degeneracies
```

Do not compare raw eigenvector phases.

## Per-iteration physical outputs

Where practical:

```text
Fermi energy
total band energy / relevant energy terms
total charge
site charges
spin moments
density norm
```

## Converged SCF outputs

At minimum:

```text
final total energy
Fermi energy
total charge
site-resolved charge
site-resolved magnetic moments where applicable
number of SCF iterations
final convergence residual
```

If GPU and CPU converge in different iteration counts but to the same physical
solution, record that explicitly.

Do not automatically declare FAIL solely because iteration counts differ.

A full-SCF pair is headline-eligible only if the final physical state passes the
scientific tolerance contract.

---

# 13. SCF-B0C-K — Starting-state contract for convergence benchmarks

SCF-to-convergence comparisons are only fair if CPU/GPU begin from the same
well-defined state.

Support at least two modes:

## Restart / near-converged mode

Start from the same saved/reference potential.

Purpose:

```text
steady production iteration performance
```

## Cold / standard initial mode

Start from the same normal production initial state.

Purpose:

```text
total time to convergence
```

Do not compare CPU cold-start convergence against GPU restart convergence.

Record:

```text
starting_state_id
```

or equivalent provenance.

---

# 14. SCF-B0C-L — Material fixtures

Consolidate at least:

## diamond Si

Use the existing validated real Si reciprocal/SCF fixture.

This is useful for:

```text
small matrix
nonmagnetic
many-k regime
```

## bcc Fe

Use the existing validated collinear Fe fixture.

This is useful for:

```text
magnetic
spd basis
spin-polarized SCF
```

Support:

```text
primitive
2x2x2
3x3x3
```

for SCF-level work where practical.

Keep:

```text
4x4x4
5x5x5
```

available for eigensolver/frozen-potential benchmark levels even if full SCF is
too expensive.

Do not require a 5x5x5 full SCF run merely to satisfy a checklist.

---

# 15. SCF-B0C-M — Physically matched k-point sampling

Retain the ACC-P distinction between:

## crossover-plane benchmarks

Deliberately vary:

```text
matrix dimension
Nk
solver strategy
```

to map hardware behavior.

## physically matched SCF benchmarks

When comparing supercells, use physically consistent k-point density as far as
practical.

Record:

```text
nominal k mesh
actual unique Nk
supercell transformation
```

Do not compare nominal mesh dimensions if the actual unique workloads differ.

---

# 16. SCF-B0C-N — Full iteration and full convergence speedups

Define separately:

```text
S_solver =
    best CPU eigensolver time
    /
    GPU eigensolver time
```

```text
S_reciprocal =
    best CPU reciprocal-phase time
    /
    GPU reciprocal-phase time
```

```text
S_iteration =
    best CPU steady SCF iteration median
    /
    GPU steady SCF iteration median
```

```text
S_convergence =
    best CPU full SCF wall
    /
    GPU full SCF wall
```

Do not call `S_solver` an SCF speedup.

Do not call `S_iteration` a convergence speedup.

All four should be available where the benchmark level supports them.

---

# 17. SCF-B0C-O — Iteration-history output

For full SCF runs, record iteration-by-iteration:

```text
iteration number
iteration wall
H(k) time
solver time
density time
mixing time
SCF residual
total energy
Fermi energy
magnetic moment summary
```

This enables later plots such as:

```text
iteration time vs iteration
SCF residual vs iteration
energy convergence vs iteration
```

Do not store huge arrays per iteration unless necessary.

Prefer compact scalar diagnostics.

---

# 18. SCF-B0C-P — Strict headline eligibility

A CPU/GPU pair may become a headline result only if:

```text
comparison key matches
numeric mode is valid for the comparison
profile closure PASS
eigensystem correctness PASS
physical output correctness PASS
normal production path used
no silent fallback
```

For `scf_convergence` additionally require:

```text
same starting state
same convergence criterion
final physical-state correctness PASS
```

Store:

```text
headline_speedup_eligible = true/false
```

plus a reason.

Example rejection reasons:

```text
physics_mismatch
numeric_mode_mismatch
starting_state_mismatch
correctness_failed
profile_failed
unsupported_gpu_route
silent_fallback_detected
```

---

# 19. SCF-B0C-Q — Canonical output package

As with KPM-B0C, produce one canonical dataset and derive:

```text
campaign.json
campaign.csv
campaign.md
```

plus:

```text
raw/
correctness/
iteration_history/
```

JSON is the full-fidelity canonical format.

CSV should contain one summarized benchmark configuration per row.

Markdown should present a concise human-readable closure report.

Do not calculate speedups independently in several output writers.

---

# 20. SCF-B0C-R — Required CSV fields

At minimum:

```text
row_id
benchmark_level
material
supercell
Natom
nmat
nnz if meaningful
nsp
SOC
basis
Nk_nominal
Nk_unique
starting_state
numeric_mode
solver_strategy
vectors
OMP_threads
BLAS_threads
GPU_device
profile_status
correctness_status
P_hk_assembly
P_eigensolver
P_occupations_fermi
P_density_build
P_potential_update
P_mixing
P_scf_io
P_scf_misc
P_scf_iteration_total
steady_iteration_median
n_scf_iterations
full_scf_wall
S_solver
S_reciprocal
S_iteration
S_convergence
headline_speedup_eligible
best_cpu_row_id
```

Use null/not-applicable rather than inventing values for benchmark levels where a
field does not apply.

---

# 21. SCF-B0C-S — Environment and provenance

Capture:

```text
git commit
git dirty state
compiler
compiler version
build type
optimization flags where available
BLAS/LAPACK vendor
MPI build/runtime status
MPI ranks
OMP settings
CPU model
physical/logical cores
RAM
OS/kernel
CUDA toolkit
CUDA driver
GPU model
selected GPU
GPU VRAM
compute capability
CUDA_VISIBLE_DEVICES
```

Optional where available:

```text
NUMA topology
CPU governor
GPU persistence mode
```

Do not make optional metadata fatal.

---

# 22. SCF-B0C-T — Required harness contract tests

Add focused tests that do not require long real SCF runs.

Test at least:

## pairing

```text
same physics -> eligible
different material -> reject
different Nk -> reject
different SOC -> reject
different starting state -> reject for convergence comparison
different numeric mode -> reject equal-precision headline
```

## correctness

```text
missing correctness -> reject
failed eigensystem correctness -> reject
failed final SCF correctness -> reject
```

## profile

```text
closure failure -> reject
```

## CPU selection

```text
best OMP row selected only among valid same-physics rows
```

## fallback

```text
explicit CUDA row with CPU fallback -> reject
```

## outputs

```text
JSON schema
CSV required fields
Markdown required sections
```

## material staging

```text
Si fixture selection
Fe fixture selection
Fe supercell selection
```

---

# 23. SCF-B0C-U — Closure campaign

Run a manageable real-material campaign to prove the harness.

This is **not** SCF-B1 yet.

Minimum recommended closure set:

## Si primitive

Run one normal reciprocal SCF case:

```text
CPU OMP 1/2/4/8
GPU FP64 or current canonical GPU precision
```

with:

```text
full eigenvectors
real k mesh
SCF iteration timing
full convergence timing
```

## Fe primitive

Same:

```text
CPU OMP 1/2/4/8
GPU
```

with magnetic correctness.

## Fe 2x2x2

At minimum:

```text
steady SCF iteration
```

and full convergence if practical.

## Fe 4x4x4 or 5x5x5

Do not require full SCF.

Run one frozen-potential eigensolver/reciprocal-phase smoke row to verify the
existing large-matrix benchmark remains integrated with the consolidated
harness.

---

# 24. SCF-B0C-V — Do not overinterpret closure timings

SCF-B0C performance rows are harness validation only.

Do not write final application recommendations from B0C.

The later SCF-B1 campaign will freeze the code and run the full benchmark matrix.

B0C should conclude only:

```text
harness ready
```

or:

```text
harness not ready because ...
```

---

# 25. SCF-B0C-W — Reconcile accelerator documentation

Update the appropriate master/current steering document so it states clearly:

```text
ACC reciprocal infrastructure/correctness campaign complete
ACC-P reciprocal performance rescue complete
KPM B0C/B1 complete
SCF-B0C complete after this task
SCF-B1 next
```

Do not rewrite historical reports.

Add a concise current-status block.

---

# 26. Acceptance criteria

SCF-B0C is complete only when:

```text
1. one canonical reciprocal/SCF benchmark workflow exists;

2. existing ACC reciprocal microbenchmark capability is preserved;

3. complete SCF iteration timing is instrumented;

4. full SCF convergence timing is supported;

5. Si and Fe production fixtures are supported;

6. Fe supercell scaling remains supported;

7. cold versus steady timing is explicit;

8. solver / reciprocal / iteration / convergence speedups are distinct;

9. CPU OMP 1/2/4/8 fairness is enforced;

10. strict physical pairing is implemented;

11. starting-state matching is enforced for convergence comparisons;

12. correctness is attached to performance rows;

13. eigenvector-phase-safe correctness checks are retained;

14. profile closure is enforced;

15. no silent CUDA fallback is possible in valid benchmark rows;

16. JSON is generated;

17. CSV is generated;

18. Markdown is generated;

19. raw logs and iteration histories are retained;

20. environment/build provenance is captured;

21. focused harness tests pass;

22. Si real closure run passes;

23. Fe real closure run passes;

24. at least one larger Fe reciprocal benchmark passes;

25. SCF-B1 can be run without another harness redesign.
```

---

# 27. Checklist

## Audit

- [ ] CURRENT HEAD inspected
- [ ] ACC benchmark drivers inspected
- [ ] ACC-P reports inspected
- [ ] SCF loop inspected
- [ ] KPM-B0C lessons inspected
- [ ] gap table written
- [ ] no duplicate harness introduced

## Benchmark levels

- [ ] eigensolver level
- [ ] reciprocal-phase level
- [ ] SCF-iteration level
- [ ] SCF-convergence level

## SCF profiling

- [ ] SCF iteration parent timer
- [ ] Hamiltonian preparation timer
- [ ] H(k) assembly timer
- [ ] eigensolver timer
- [ ] eigenpair transfer timer
- [ ] occupations/Fermi timer
- [ ] density-build timer
- [ ] charge/spin accumulation timer
- [ ] potential-update timer
- [ ] mixing timer
- [ ] SCF I/O timer
- [ ] SCF misc timer
- [ ] profile closure <=3%
- [ ] misc <5% or explained

## Reciprocal profiling

- [ ] cold process separated
- [ ] backend/context init separated
- [ ] first solve separated
- [ ] steady solve median retained
- [ ] H2D retained
- [ ] solver retained
- [ ] D2H retained
- [ ] GPU workspace/memory retained

## Physics pairing

- [ ] material included
- [ ] supercell included
- [ ] basis included
- [ ] nsp included
- [ ] SOC included
- [ ] k mesh included
- [ ] actual unique Nk included
- [ ] occupations included
- [ ] mixing settings included
- [ ] convergence thresholds included
- [ ] starting state included
- [ ] feature flags included where relevant
- [ ] stable comparison key generated

## Precision

- [ ] Hamiltonian precision recorded
- [ ] solver precision recorded
- [ ] density precision recorded
- [ ] canonical SCF precision recorded
- [ ] numeric mode classified honestly
- [ ] unsupported full-FP32 SCF not fabricated

## CPU fairness

- [ ] OMP=1
- [ ] OMP=2
- [ ] OMP=4
- [ ] OMP=8
- [ ] BLAS threads controlled
- [ ] oversubscription prevented
- [ ] best CPU iteration row selected fairly
- [ ] best CPU convergence row selected fairly

## Correctness

- [ ] eigenvalues checked
- [ ] residuals checked
- [ ] orthogonality checked
- [ ] degenerate-subspace comparison retained
- [ ] Fermi energy checked
- [ ] total charge checked
- [ ] site charge checked
- [ ] magnetic moments checked where applicable
- [ ] total energy checked
- [ ] final convergence residual checked
- [ ] iteration-count differences reported rather than blindly failed
- [ ] final physical-state correctness gate implemented

## Starting states

- [ ] restart/near-converged mode supported
- [ ] normal/cold initial mode supported
- [ ] starting-state ID recorded
- [ ] mismatched starting states rejected for convergence speedup

## Materials

- [ ] Si primitive staged
- [ ] Fe primitive staged
- [ ] Fe 2x2x2 staged
- [ ] Fe 3x3x3 staged where practical
- [ ] Fe 4x4x4/5x5x5 reciprocal benchmark retained
- [ ] real production H(R)->H(k) used
- [ ] no synthetic performance fixture substituted

## Outputs

- [ ] canonical JSON
- [ ] canonical CSV
- [ ] canonical Markdown
- [ ] raw logs
- [ ] correctness evidence
- [ ] iteration history
- [ ] row IDs
- [ ] sample IDs
- [ ] speedup eligibility reasons retained

## Environment

- [ ] git commit
- [ ] dirty state
- [ ] compiler/version
- [ ] build flags
- [ ] BLAS/LAPACK
- [ ] MPI status
- [ ] OMP settings
- [ ] CPU model
- [ ] RAM
- [ ] OS/kernel
- [ ] CUDA toolkit
- [ ] driver
- [ ] GPU model
- [ ] selected GPU
- [ ] VRAM
- [ ] compute capability

## Harness tests

- [ ] physics mismatch rejection
- [ ] Nk mismatch rejection
- [ ] starting-state mismatch rejection
- [ ] precision mismatch rejection
- [ ] correctness missing/fail rejection
- [ ] profile-fail rejection
- [ ] silent-fallback rejection
- [ ] best-CPU selection test
- [ ] JSON schema test
- [ ] CSV schema test
- [ ] Markdown test
- [ ] Si fixture test
- [ ] Fe fixture test

## Closure campaign

- [ ] Si CPU OMP sweep
- [ ] Si GPU run
- [ ] Si full convergence comparison
- [ ] Fe primitive CPU OMP sweep
- [ ] Fe primitive GPU run
- [ ] Fe primitive full convergence comparison
- [ ] Fe 2x2x2 steady iteration
- [ ] larger Fe reciprocal smoke row
- [ ] all headline closure rows correctness PASS
- [ ] SCF-B1 declared next

---

# 28. Required closeout report

Write:

```text
docs/dev/RS_LMTO_ASA_SCF_B0C_REPORT.md
```

Recommended structure:

## 1. Existing infrastructure

Summarize what ACC/ACC-P already provided.

## 2. Gap analysis

State what was missing before SCF-B0C.

## 3. SCF profiling model

Show the exclusive iteration decomposition.

## 4. Pairing and correctness contract

Explain the comparison key and final-state validation.

## 5. Material fixtures

Document Si and Fe.

## 6. Output/provenance structure

Document JSON/CSV/Markdown/raw/iteration-history output.

## 7. Closure evidence

Show the small Si/Fe validation campaign.

## 8. Remaining limitations

For example:

```text
full SCF not practical for largest Fe supercells
no multi-GPU
no MPI scaling claim
full FP32 SCF unsupported if applicable
```

## 9. Next task

State explicitly:

```text
SCF-B1 — frozen final CPU/GPU SCF benchmark campaign
```

---

# 29. Commit message

Use:

```text
Consolidate SCF CPU GPU benchmark harness
```

---

# 30. Final mindset

This task is not:

```text
make the GPU faster
```

It is:

```text
make the SCF performance evidence complete, fair, and reproducible
```

The existing ACC work already tells us a great deal about reciprocal
eigensolution.

SCF-B0C must connect that evidence to the scientific workflow users actually
care about:

```text
How much faster is an SCF iteration?
How much faster is convergence?
Does it converge to the same physical solution?
```

At completion, SCF-B1 should be mostly:

```text
freeze
run
collect
plot
interpret
```

not another benchmark-harness development cycle.
