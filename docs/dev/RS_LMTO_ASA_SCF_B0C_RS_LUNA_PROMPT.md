# SCF-B0C-RS — Add the real-space SCF lane to the consolidated CPU/GPU benchmark harness

You are working on the **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task is a **targeted addendum** to the completed:

```text
SCF-B0C
```

and must reuse the infrastructure already established there.

The existing SCF-B0C benchmark entry point is:

```text
tests/benchmarks/scf_b0c.py
```

or whatever equivalent/common code CURRENT HEAD now uses after SCF-B0C.

The reciprocal/k-space lane already measures:

```text
eigensolver
    ->
reciprocal phase
    ->
SCF iteration
    ->
SCF convergence
```

with:

```text
strict physical pairing
CPU OMP 1/2/4/8 fairness
precision metadata
correctness gates
profile closure
cold/steady timing
full convergence timing
JSON / CSV / Markdown
raw logs
iteration histories
environment provenance
```

The purpose of SCF-B0C-RS is to add the **real-space SCF route** to that same
framework so that the later SCF-B1 campaign can benchmark both:

```text
real-space SCF
and
reciprocal/k-space SCF
```

without another harness redesign.

---

# 1. Core objective

At completion, the common SCF benchmark harness must support:

```text
scf_route = reciprocal
scf_route = real_space
```

with the same high-level contract.

For the real-space route, the benchmark hierarchy is:

```text
RS solver kernel
    ->
RS electronic-structure phase
    ->
SCF iteration
    ->
SCF convergence
```

The final harness must answer:

1. How long does the production RS solver kernel take?
2. How much GPU acceleration is obtained for that kernel?
3. How much of that acceleration survives the complete RS electronic-structure phase?
4. How much survives one full SCF iteration?
5. How much survives to total SCF convergence?
6. Does CPU and GPU RS SCF converge to the same physical state?
7. Which RS stage becomes dominant after kernel acceleration?
8. Is the result ready to be compared later with the reciprocal/k-space SCF route?

SCF-B0C-RS is **not** the final performance campaign.

It is harness completion and closure evidence only.

---

# 2. Non-negotiable rules

## 2.1 Inspect CURRENT HEAD before editing

Read at minimum:

```text
SCF-B0C benchmark driver and common utilities
SCF-B0C closeout report
real-space SCF driver
real-space Hamiltonian construction
production block-recursion implementation
production GPU block-recursion implementation
production Chebyshev electronic-structure implementation
production GPU Chebyshev implementation if present
Lanczos/scalar recursion path if still production-relevant
Green-function / DOS reconstruction code
Fermi-level / energy-integration path
charge/spin-density construction
potential update and mixing
existing RS correctness/regression tests
existing RS CPU/GPU benchmark drivers
ACC/Phase-III accelerator documentation
```

Do not rely on historical task descriptions if CURRENT HEAD differs.

## 2.2 Do not create a second independent benchmark framework

Extend:

```text
tests/benchmarks/scf_b0c.py
```

or factor shared code cleanly underneath it.

Do not create:

```text
an unrelated RS benchmark stack
with different row schemas
different correctness rules
different provenance
or different speedup definitions
```

The RS and reciprocal lanes must share the same canonical campaign model.

## 2.3 Preserve physics

Do not change:

```text
LMTO Hamiltonian conventions
structure constants
basis conventions
spin ordering
recursion equations
terminators
Chebyshev scaling
kernel polynomial definitions
Fermi-level calculation
density construction
SCF mixing
SCF convergence criteria
potential equations
reference energies
SOC treatment
```

to improve timing.

## 2.4 Do not modify protected legacy atomic LMTO routines unless strictly necessary

The older atomic/self-consistency LMTO routines are scientifically sensitive.

Prefer instrumentation and benchmark integration around existing interfaces.

If a protected legacy routine must be touched:

- document why;
- make the minimum possible change;
- do not refactor it opportunistically.

## 2.5 CPU remains the scientific oracle

GPU performance evidence is valid only if the corresponding CPU/GPU physics
comparison passes.

## 2.6 No silent fallback

A row labelled:

```text
backend = cuda
```

must not execute the CPU solver kernel silently.

Fallback/unsupported conditions must become:

```text
UNSUPPORTED
SKIPPED
FAIL
```

with explicit reason.

## 2.7 This task is not a new optimization campaign

Do not rewrite recursion kernels merely because profiling shows them to be slow.

Do not introduce a new memory layout, algorithm, terminator, or GPU kernel in
SCF-B0C-RS unless a benchmark-integrity defect makes measurement impossible.

Measure first.

Optimization decisions belong after SCF-B1.

---

# 3. Reuse the common SCF-B0C row contract

Every row must continue to contain the common physical metadata already
established by SCF-B0C, including where applicable:

```text
material
fixture revision
supercell
Natom
Hamiltonian dimension
nsp
SOC
basis / lmax
structure-constant backend
smearing / electronic temperature
electron count
Fermi-level policy
mixing method
mixing parameters
SCF convergence thresholds
starting-state identity
potential / reference-state identity
numeric mode
feature flags
git/build/environment provenance
```

Add:

```text
scf_route
```

with:

```text
reciprocal
real_space
```

For RS rows add route-specific metadata such as:

```text
rs_solver
rs_backend
recursion_depth / nrec
block_size where applicable
terminator
Chebyshev_order where applicable
Chebyshev_scaling / spectral bounds where applicable
RS boundary/PBC mode
stochastic settings if genuinely part of SCF
```

Do not fill irrelevant fields with invented values.

Use null / not-applicable.

---

# 4. Real-space solver taxonomy

Inspect CURRENT HEAD and classify only genuine production routes.

Expected categories may include:

```text
block_recursion
lanczos / scalar_recursion
chebyshev
```

but do **not** force all three into the campaign if CURRENT HEAD does not use them
as production SCF solvers.

Mandatory policy:

## Block recursion

Treat as a primary RS SCF route if it is the normal metallic production solver.

It is especially important for:

```text
bcc Fe
other metallic magnetic systems
```

## Chebyshev electronic-structure SCF

Treat as a separate primary route if it is a genuine production electronic-
structure solver.

Keep its benchmark identity clearly separate from:

```text
KPM / Kubo-Bastin transport
```

The latter was already benchmarked in the KPM campaign.

## Scalar/Lanczos recursion

Include only if it remains a supported/meaningful production or reference route.

Do not multiply the B0C-RS closure matrix merely because an old code path exists.

---

# 5. Benchmark levels for the RS route

Extend the common benchmark-level taxonomy.

For `scf_route = real_space`, support:

```text
rs_kernel
rs_electronic_structure
scf_iteration
scf_convergence
```

## rs_kernel

Measures the repeated numerical solver kernel itself.

Examples:

```text
block-recursion recurrence
Lanczos recurrence
Chebyshev recurrence
```

Use real production Hamiltonians.

## rs_electronic_structure

Measures the complete RS electronic-structure phase needed by one SCF step,
including as appropriate:

```text
solver recurrence
Green-function / spectral reconstruction
DOS / LDOS reconstruction
energy integration
Fermi-level determination
charge/spin-density generation
```

Do not include potential update/mixing unless the production architecture makes
that boundary impossible to separate.

## scf_iteration

Measures one complete real-space SCF iteration.

## scf_convergence

Measures the complete real-space SCF calculation to the normal production
convergence criterion.

Do not call an `rs_kernel` speedup an SCF speedup.

---

# 6. RS-specific timing decomposition

Add an exclusive real-space profiling decomposition consistent with CURRENT HEAD.

Use names appropriate to the code, but provide the equivalent of:

```text
P_rs_hamiltonian_prepare
P_rs_solver_kernel
P_rs_green_function
P_rs_spectral_reconstruct
P_rs_energy_integration
P_rs_fermi
P_rs_density_build
P_rs_charge_spin_accumulate
```

plus the common SCF stages:

```text
P_potential_update
P_mixing
P_scf_io
P_scf_misc
P_scf_iteration_total
```

Not every RS solver will use every stage.

For example:

- Chebyshev may not have a Green-function stage;
- block recursion may combine reconstruction and integration;
- Fermi-level handling may be embedded in an existing production routine.

Do not invent fake stage boundaries.

Instead:

1. inspect the real production call graph;
2. define the narrowest meaningful exclusive stages;
3. document the mapping in the closeout report.

Require:

```text
profile_closure_error <= 0.03
```

for headline-eligible closure rows.

Require:

```text
P_scf_misc <= 5%
```

or explain the residual.

---

# 7. Route-neutral SCF timing must remain identical

Do not create separate definitions of:

```text
SCF iteration
full SCF wall
starting state
convergence residual
```

for RS and reciprocal routes.

The common meanings must remain:

```text
P_scf_iteration_total
first_scf_iteration
steady_scf_iteration_median
full_scf_wall
n_scf_iterations
final_convergence_residual
```

This is necessary for later RS-versus-k-space comparison.

---

# 8. GPU detail timing for RS kernels

Where CUDA RS kernels already expose detail metrics, retain or add benchmark-safe
snapshot/delta reporting for:

```text
T_rs_H2D
T_rs_kernel
T_rs_D2H
T_rs_sync
T_rs_setup
```

and:

```text
H2D bytes
D2H bytes
workspace bytes
allocation count
workspace reuse count
```

where practical.

Do not force every CUDA diagnostic into the SCF parent timer.

GPU detail timers are nested diagnostics.

The exclusive SCF parent sum must not double-count them.

---

# 9. CPU/GPU kernel timer-equivalence gate

Before running the closure campaign, verify for each mandatory RS solver that:

```text
CPU kernel timer
```

and:

```text
GPU kernel timer
```

cover semantically equivalent numerical work.

For block recursion, explicitly check whether either timing includes:

```text
coefficient initialization
Hamiltonian packing
host/device staging
recurrence
orthogonalization
terminator/reconstruction
postprocessing
```

For Chebyshev, explicitly check whether either timing includes:

```text
spectral scaling
moment initialization
recurrence
moment accumulation
reconstruction
```

If CPU/GPU boundaries differ:

- fix reporting/timer boundaries;
- do not alter scientific work;
- document the final contract.

Do not compare pure GPU kernel execution against a larger CPU phase and label it
kernel speedup.

---

# 10. Precision taxonomy

Do not infer precision from the name of a GPU routine.

Record independently:

```text
Hamiltonian precision
RS kernel precision
recursion coefficient / moment precision
spectral reconstruction precision
density precision
canonical SCF precision
```

Classify:

```text
fp64
fp32
mixed
```

from the actual data path.

If the GPU RS solver uses FP32 internally but returns into an FP64 SCF path:

```text
numeric_mode = mixed
```

unless the full route is genuinely FP32.

Do not fabricate a full-FP32 RS SCF path for symmetry with KPM transport.

---

# 11. CPU fairness

For every closure workload that will later support a CPU/GPU conclusion, run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

where meaningful.

Control:

```text
MKL_NUM_THREADS
OPENBLAS_NUM_THREADS
OMP_PROC_BIND
OMP_PLACES
```

Avoid oversubscription.

For sparse/recursion CPU kernels, record whether:

```text
OpenMP
BLAS threading
both
or neither
```

actually influence the stage.

Do not assume that more threads are faster.

Retain all rows.

---

# 12. RS correctness contract

Every CPU/GPU RS pair must carry structured correctness evidence.

## 12.1 Kernel-level correctness

Use solver-appropriate invariants.

### Block recursion / Lanczos

Where practical compare:

```text
selected recursion coefficients
selected block norms
selected Green-function values
selected DOS/LDOS values
integrated spectral weight
```

Do not require huge full-array dumps.

### Chebyshev electronic structure

Where practical compare:

```text
selected moments
DOS/LDOS
integrated spectral weight
electron count
```

The comparison should detect a broken recurrence without making the benchmark
output enormous.

## 12.2 SCF-level correctness

For every real-space SCF closure row compare at minimum:

```text
Fermi energy
total charge
site charges
magnetic moments where applicable
total energy / stable production energy quantities
final SCF residual
```

If stable and already available, include:

```text
density norm
band energy
DOS-integrated electron count
```

## 12.3 Full convergence

CPU/GPU may require different iteration counts.

That is not automatically a failure.

Headline eligibility requires:

```text
both satisfy the normal convergence criterion
and
final physical observables agree within the scientific tolerance contract
```

Retain the iteration-count difference.

---

# 13. Starting-state contract

Reuse SCF-B0C exactly.

For full-convergence comparisons:

```text
CPU and GPU must start from the same state
```

Support:

```text
standard initial state
same restart / near-converged state
```

Record:

```text
starting_state_id
```

Do not compare CPU cold-start with GPU restart.

---

# 14. Material / solver closure fixtures

SCF-B0C-RS must prove the harness on physically sensible RS workloads.

Do not use synthetic matrices as the main closure evidence.

## 14.1 Mandatory metallic fixture — bcc Fe

Use a validated periodic bcc-Fe RS fixture.

Primary purpose:

```text
block recursion / metallic RS SCF
magnetic correctness
CPU/GPU recurrence timing
```

Choose a system large enough to exercise the production RS kernels meaningfully
but small enough for closure work.

Prefer an existing canonical fixture.

If no canonical closure size exists, a reasonable starting point is a periodic:

```text
6x6x6
```

or:

```text
8x8x8
```

Fe supercell.

Do not hard-code that choice if CURRENT HEAD already has a better validated
fixture.

Run:

```text
CPU OMP 1/2/4/8
GPU current production route
```

with:

```text
steady SCF iteration
```

and preferably:

```text
full SCF convergence
```

if practical.

Also retain one larger kernel/steady-iteration smoke workload if an existing
10x10x10-class Fe benchmark is already available.

Do not require full convergence of the largest smoke case.

## 14.2 Mandatory insulating/semiconducting fixture — Si

Use the validated periodic Si RS fixture for the Chebyshev electronic-structure
route if CURRENT HEAD supports it as production SCF.

A typical useful Chebyshev test may use approximately:

```text
M ~ 200
```

and a periodic supercell large enough to exercise the sparse kernel, e.g.
an existing 8x8x8–12x12x12-class fixture.

Use CURRENT validated fixture/settings rather than blindly adopting these
numbers.

Run:

```text
CPU
GPU if a production GPU Chebyshev SCF path exists
```

with:

```text
steady iteration
```

and full convergence if practical.

If there is no production GPU Chebyshev SCF route:

```text
record UNSUPPORTED
```

and do not create one in B0C-RS.

## 14.3 Optional scalar/Lanczos fixture

Include only if CURRENT HEAD still treats it as a supported production/reference
RS SCF path.

Do not make it mandatory merely for completeness.

---

# 15. Recursion-depth / polynomial-order metadata

For block/Lanczos recursion record:

```text
nrec
terminator
block dimension
```

where applicable.

For Chebyshev record:

```text
M
kernel/damping scheme if applicable
spectral bounds / scaling policy
```

These are part of the physical/numerical comparison key.

CPU and GPU rows with different recursion depth or Chebyshev order are not
valid backend pairs.

---

# 16. Distinguish RS Chebyshev SCF from KPM transport

This must be explicit in code, metadata, and documentation.

Use distinct labels such as:

```text
solver_family = chebyshev_electronic_structure
```

versus the already completed:

```text
KPM_Kubo_Bastin_transport
```

Do not reuse KPM transport timing fields or call a transport benchmark an SCF
benchmark.

It is fine to reuse generic utility code if appropriate.

---

# 17. Prepare for later RS-versus-k-space comparison

SCF-B0C-RS should make future cross-route comparison possible without claiming
it yet.

Add a route-neutral field such as:

```text
physical_case_id
```

or:

```text
cross_route_case_id
```

that can identify physically comparable RS and reciprocal calculations.

A future RS-versus-k-space pair should require agreement in:

```text
material
lattice / structure
basis
spin/SOC state
electronic temperature / smearing
electron count
SCF functional/potential model
convergence thresholds
starting-state class
```

while allowing route-specific numerical controls to differ:

```text
RS supercell size
recursion depth / M
k-space k mesh
```

because those are different convergence parameters.

Do **not** define cross-route speedup merely from equal nominal input sizes.

Later SCF-B1 must compare routes only after physical accuracy is matched.

---

# 18. Route-specific convergence metadata

To support that later comparison, retain convergence-control metadata.

For RS:

```text
supercell size
nrec / M
terminator
spectral broadening / kernel
```

For reciprocal:

```text
k mesh
Nk_unique
DOS/integration method
```

These must remain visible in the canonical dataset.

---

# 19. Speedup definitions for the RS lane

Define:

```text
S_rs_kernel =
    best valid CPU RS kernel time
    /
    valid GPU RS kernel time
```

```text
S_rs_phase =
    best valid CPU RS electronic-structure phase time
    /
    valid GPU RS electronic-structure phase time
```

```text
S_iteration =
    best valid CPU steady SCF iteration median
    /
    valid GPU steady SCF iteration median
```

```text
S_convergence =
    best valid CPU full SCF wall
    /
    valid GPU full SCF wall
```

Reuse the existing common:

```text
S_iteration
S_convergence
```

fields where possible.

Do not create duplicate meanings.

For reciprocal rows:

```text
S_solver
S_reciprocal
```

remain as already defined.

For route-neutral reporting, optionally add:

```text
S_route_kernel
S_route_phase
```

with clear route-specific interpretation.

---

# 20. Headline eligibility

Reuse the SCF-B0C eligibility mechanism.

An RS CPU/GPU pair may be performance-eligible only if:

```text
scf_route matches
physical comparison key matches
solver family matches
recursion depth / M matches
numeric mode is a valid comparison
profile closure PASS
kernel correctness PASS
SCF physical correctness PASS
normal production path used
no silent fallback
```

For full convergence additionally require:

```text
same starting state
same convergence criterion
final physical state PASS
```

Store explicit rejection reasons.

Suggested additions:

```text
rs_solver_mismatch
recursion_depth_mismatch
chebyshev_order_mismatch
terminator_mismatch
rs_kernel_correctness_failed
```

---

# 21. Canonical output package

Do not create a separate RS output format.

The same campaign package must support both routes:

```text
campaign.json
campaign.csv
campaign.md
raw/
correctness/
iteration_history/
```

If useful, also include:

```text
route_detail/
```

for compact RS kernel diagnostics.

JSON remains the full-fidelity source of truth.

CSV remains one summarized configuration per row.

Markdown is derived from JSON.

---

# 22. Required schema additions

Add at minimum:

```text
scf_route
rs_solver
rs_backend
recursion_depth
block_size
terminator
chebyshev_order
chebyshev_kernel
spectral_bounds_policy
P_rs_hamiltonian_prepare
P_rs_solver_kernel
P_rs_green_function
P_rs_spectral_reconstruct
P_rs_energy_integration
P_rs_fermi
P_rs_density_build
T_rs_H2D
T_rs_kernel
T_rs_D2H
S_rs_kernel
S_rs_phase
cross_route_case_id
```

where applicable.

Do not require every field to be populated for every solver.

For reciprocal rows these should normally be null.

For RS rows, reciprocal-only fields should normally be null.

---

# 23. Harness contract tests

Add focused tests for the RS lane.

At minimum:

## Routing

```text
scf_route=real_space selects RS benchmark path
scf_route=reciprocal preserves existing behavior
```

## Pairing

Reject:

```text
block_recursion vs chebyshev
different nrec
different M
different terminator
different numeric mode for equal-precision headline
different physical fixture
different starting state for convergence
```

## Timing

```text
RS profile closure failure rejects headline
nested CUDA timers are not double-counted
```

## Correctness

```text
missing RS kernel correctness rejects row
failed final SCF correctness rejects row
```

## Fallback

```text
CUDA-labelled RS row executing CPU kernel is rejected
```

## Outputs

Ensure:

```text
JSON round trip
CSV required RS fields
Markdown route label
```

## Regression

Existing reciprocal SCF-B0C tests must continue to pass.

---

# 24. Closure campaign

Run a bounded real-material closure campaign.

This is not the final B1 performance matrix.

## Closure A — Fe block recursion

Mandatory if production block recursion exists.

At minimum:

```text
one validated Fe RS fixture
CPU OMP 1/2/4/8
GPU production route
same nrec
same terminator
same starting state
```

Require:

```text
kernel correctness PASS
profile closure PASS
SCF physical correctness PASS
```

Obtain:

```text
steady SCF iteration
```

and full convergence if practical.

## Closure B — Si Chebyshev electronic-structure SCF

Mandatory if both CPU and GPU production routes exist.

At minimum:

```text
one validated Si RS fixture
same Chebyshev M
same scaling
same kernel
```

Require correctness/profile closure.

If GPU Chebyshev SCF is unsupported:

```text
record this explicitly
```

and keep the CPU path staged for later RS-vs-k-space scientific comparison.

## Closure C — larger RS kernel smoke

Run one existing production-sized RS workload to prove:

```text
large sparse state
GPU memory/accounting
kernel timing
no hidden allocation/fallback
```

No full SCF convergence is required.

---

# 25. B0C-RS interpretation discipline

The closeout report may state:

```text
RS harness ready
```

or:

```text
RS harness not ready because ...
```

It must **not** make the final production recommendation:

```text
RS GPU is faster than CPU
RS is faster than k-space
```

from closure rows alone.

Those are SCF-B1 conclusions.

Preserve negative closure measurements but label them correctly as harness
validation evidence.

---

# 26. Required closeout report

Write:

```text
docs/dev/RS_LMTO_ASA_SCF_B0C_RS_REPORT.md
```

Recommended structure:

## 1. Purpose

Explain that this is an addendum to SCF-B0C.

## 2. Existing RS infrastructure

Document:

```text
production block recursion
GPU recursion
Chebyshev electronic structure
other supported RS solvers
existing RS tests
```

as actually found in CURRENT HEAD.

## 3. Harness integration

Explain:

```text
scf_route = reciprocal / real_space
```

and shared contracts.

## 4. RS profiling model

Show actual production call boundaries and the exclusive timing decomposition.

## 5. CPU/GPU kernel timer contract

Document like-for-like timing boundaries.

## 6. Precision model

State actual FP64/FP32/mixed routes.

## 7. Correctness model

Document solver-level and SCF-level validation.

## 8. Fixtures

Document Fe and Si closure fixtures.

## 9. Closure evidence

Show:

```text
CPU/GPU correctness
profile closure
steady iteration timing
full convergence where available
```

without turning closure numbers into final performance policy.

## 10. Remaining limitations

Examples:

```text
GPU Chebyshev SCF unsupported
full convergence not run on largest RS fixture
single GPU only
single MPI rank
some legacy solver retained only as reference
```

## 11. SCF-B1 readiness

State explicitly whether the unified final campaign can proceed.

---

# 27. Update the steering/status documentation

Update the current accelerator/SCF steering document to show:

```text
SCF-B0C       DONE
SCF-B0C-RS    DONE
SCF-B1        NEXT
```

and describe SCF-B1 as a unified:

```text
RS CPU vs RS GPU
+
k-space CPU vs k-space GPU
+
selected physically matched RS vs k-space comparisons
```

campaign.

Do not rewrite historical ACC reports.

---

# 28. Acceptance criteria

SCF-B0C-RS is complete only when:

1. `scf_route=real_space` is integrated into the existing SCF benchmark framework.
2. The reciprocal B0C route still works unchanged.
3. Production RS solver families are identified from CURRENT HEAD.
4. Block recursion is supported if it is a production solver.
5. Chebyshev electronic-structure SCF is supported if it is a production solver.
6. KPM transport remains clearly separate.
7. RS kernel timing is instrumented.
8. RS electronic-structure phase timing is instrumented.
9. Complete SCF iteration timing remains route-neutral.
10. Full convergence timing remains route-neutral.
11. CPU/GPU RS kernel timer boundaries are verified equivalent.
12. Precision metadata are explicit.
13. No fake full-FP32 SCF route is introduced.
14. RS solver-level correctness is attached to rows.
15. Final SCF physical correctness is attached to rows.
16. CPU OMP fairness is preserved.
17. Profile closure <=3% is enforced.
18. Silent fallback is rejected.
19. JSON/CSV/Markdown support RS rows.
20. Raw logs and iteration histories are retained.
21. Fe RS closure passes.
22. Si Chebyshev closure passes or unsupported GPU status is documented.
23. One larger RS smoke workload passes.
24. Cross-route metadata are sufficient for later accuracy-matched comparison.
25. SCF-B1 can start without another harness redesign.

---

# 29. Checklist

## Audit

- [ ] CURRENT HEAD inspected
- [ ] SCF-B0C code inspected
- [ ] SCF-B0C report inspected
- [ ] RS SCF driver inspected
- [ ] block recursion inspected
- [ ] GPU block recursion inspected
- [ ] Chebyshev SCF inspected
- [ ] GPU Chebyshev support status determined
- [ ] scalar/Lanczos status determined
- [ ] RS correctness tests inspected
- [ ] no duplicate harness created

## Shared architecture

- [ ] `scf_route` added
- [ ] reciprocal route preserved
- [ ] real-space route integrated
- [ ] shared comparison-key machinery reused
- [ ] shared starting-state machinery reused
- [ ] shared correctness machinery reused
- [ ] shared provenance reused
- [ ] shared JSON/CSV/Markdown reused

## RS solver metadata

- [ ] solver family recorded
- [ ] backend recorded
- [ ] nrec recorded where applicable
- [ ] block size recorded where applicable
- [ ] terminator recorded where applicable
- [ ] Chebyshev M recorded where applicable
- [ ] Chebyshev kernel recorded where applicable
- [ ] spectral scaling/bounds recorded where applicable

## Profiling

- [ ] RS Hamiltonian/setup timing
- [ ] RS solver-kernel timing
- [ ] GF timing where applicable
- [ ] spectral reconstruction timing where applicable
- [ ] energy-integration timing
- [ ] Fermi timing
- [ ] density timing
- [ ] common potential-update timing
- [ ] common mixing timing
- [ ] common iteration total
- [ ] profile closure <=3%
- [ ] misc <=5% or explained

## GPU detail

- [ ] H2D timing
- [ ] kernel timing
- [ ] D2H timing
- [ ] synchronization timing where meaningful
- [ ] bytes recorded
- [ ] workspace recorded
- [ ] allocation/reuse diagnostics retained where practical
- [ ] nested GPU timers not double-counted

## Timer equivalence

- [ ] block recursion CPU/GPU boundaries checked
- [ ] Chebyshev CPU/GPU boundaries checked if supported
- [ ] setup not hidden asymmetrically
- [ ] transfer not hidden asymmetrically
- [ ] reconstruction not hidden asymmetrically
- [ ] final timer contract documented

## Precision

- [ ] Hamiltonian precision
- [ ] RS kernel precision
- [ ] coefficient/moment precision
- [ ] reconstruction precision
- [ ] density precision
- [ ] canonical SCF precision
- [ ] numeric mode classified honestly
- [ ] no fabricated full-FP32 SCF

## Correctness

- [ ] solver-specific kernel correctness
- [ ] integrated spectral weight
- [ ] Fermi energy
- [ ] total charge
- [ ] site charges
- [ ] magnetic moments where applicable
- [ ] total energy / stable energy quantities
- [ ] final residual
- [ ] iteration-count differences retained
- [ ] final physical-state gate

## CPU fairness

- [ ] OMP 1
- [ ] OMP 2
- [ ] OMP 4
- [ ] OMP 8
- [ ] BLAS threads controlled
- [ ] oversubscription prevented
- [ ] all rows retained

## Fe closure

- [ ] validated production fixture selected
- [ ] block recursion route selected
- [ ] CPU sweep run
- [ ] GPU run
- [ ] kernel correctness PASS
- [ ] profile closure PASS
- [ ] steady iteration measured
- [ ] full convergence obtained if practical
- [ ] magnetic final-state correctness PASS where converged

## Si closure

- [ ] validated production fixture selected
- [ ] Chebyshev electronic-structure route selected
- [ ] CPU row run
- [ ] GPU row run if production GPU route exists
- [ ] unsupported GPU status recorded otherwise
- [ ] M/scaling matched
- [ ] correctness PASS
- [ ] profile closure PASS
- [ ] steady iteration measured
- [ ] full convergence obtained if practical

## Larger RS smoke

- [ ] real production-sized fixture
- [ ] kernel executes
- [ ] memory diagnostics
- [ ] no silent fallback
- [ ] timing retained
- [ ] no full-convergence requirement fabricated

## Cross-route preparation

- [ ] cross-route case ID supported
- [ ] shared physical metadata sufficient
- [ ] RS convergence controls retained
- [ ] reciprocal convergence controls retained
- [ ] no invalid equal-size RS/k-space comparison introduced

## Tests

- [ ] route-selection test
- [ ] solver-family mismatch rejection
- [ ] nrec mismatch rejection
- [ ] M mismatch rejection
- [ ] terminator mismatch rejection
- [ ] precision mismatch rejection
- [ ] correctness failure rejection
- [ ] profile failure rejection
- [ ] fallback rejection
- [ ] JSON test
- [ ] CSV test
- [ ] Markdown test
- [ ] reciprocal regression tests still pass

## Closeout

- [ ] SCF-B0C-RS report written
- [ ] steering document updated
- [ ] limitations stated
- [ ] B1 readiness stated
- [ ] no optimization campaign opened inside B0C-RS
- [ ] checkboxes ticked honestly

---

# 30. Commit message

Use:

```text
Add real-space SCF benchmark lane
```

---

# 31. Final mindset

SCF-B0C-RS is not:

```text
benchmark every recursion algorithm
```

and it is not:

```text
make RS CUDA faster
```

It is:

```text
connect the real-space production SCF route to the same trustworthy benchmark
contract already built for reciprocal SCF
```

At completion, the later SCF-B1 campaign must be able to run:

```text
A. RS CPU vs RS GPU
   RS kernel
      ->
   RS electronic-structure phase
      ->
   SCF iteration
      ->
   SCF convergence

B. k-space CPU vs k-space GPU
   eigensolver
      ->
   reciprocal phase
      ->
   SCF iteration
      ->
   SCF convergence

C. selected RS vs k-space comparisons
   only after physical accuracy is matched
```

without redesigning the harness again.
