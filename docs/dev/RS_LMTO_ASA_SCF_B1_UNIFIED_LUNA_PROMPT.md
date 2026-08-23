# SCF-B1 — Frozen unified real-space + reciprocal CPU/GPU SCF benchmark campaign

You are working on the **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task follows the completed benchmark-infrastructure work:

```text
ACC-00 ... ACC-09
ACC-P0 ... ACC-P4
SCF-B0C
SCF-B0C-RS
```

SCF-B0C and SCF-B0C-RS have now established one shared benchmark framework for:

```text
scf_route = reciprocal
scf_route = real_space
```

The common harness provides:

```text
real production fixtures
strict physics/numerical pairing
CPU OMP 1/2/4/8 fairness
SCF iteration profiling
full convergence timing
route-specific kernel/phase timing
profile closure checks
final-state correctness
fallback detection
JSON / CSV / Markdown outputs
raw logs
iteration histories
hardware/build provenance
```

The real-space route supports the production solver identities found by SCF-B0C-RS:

```text
block      -> block Haydock recursion + block Green reconstruction
chebyshev  -> production Chebyshev recursion + spectral reconstruction
lanczos    -> scalar/Lanczos recursion + scalar Green path
```

The reciprocal route retains the production:

```text
H(R) -> H(k) -> dense eigensystem -> occupations/density -> SCF
```

path.

SCF-B1 is **not another harness-development task**.

Its purpose is to:

```text
freeze the implementation
run the final real-material CPU/GPU campaign
measure the RS and reciprocal routes honestly
perform a small physically matched RS-vs-k-space comparison
write the final workload/performance recommendation
```

---

# 1. The three B1 lanes

SCF-B1 contains three distinct evidence lanes.

## Lane A — Real-space CPU vs GPU

Measure:

```text
RS solver kernel
    ->
RS electronic-structure phase
    ->
complete SCF iteration
    ->
SCF convergence
```

for production RS workloads.

## Lane B — Reciprocal/k-space CPU vs GPU

Measure:

```text
dense eigensolver
    ->
reciprocal electronic-structure phase
    ->
complete SCF iteration
    ->
SCF convergence
```

for production reciprocal workloads.

## Lane C — Physically matched RS vs reciprocal

For a small number of selected materials, compare the two formulations only
after demonstrating comparable physical accuracy.

Do not mix the three questions.

A kernel speedup is not an SCF speedup.

A CPU/GPU route comparison is not automatically an RS-vs-k-space comparison.

An RS-vs-k-space timing comparison is invalid unless the physical accuracy
contract is satisfied.

---

# 2. Freeze policy

Before any headline timing:

```bash
git status
git rev-parse HEAD
```

Require a clean tracked implementation state.

Record at minimum:

```text
frozen_git_commit
git_dirty
compiler
compiler version
build type
optimization flags
BLAS/LAPACK vendor
MPI build/runtime status
MPI ranks
OMP settings
CPU model
RAM
OS/kernel
CUDA toolkit
CUDA driver
GPU model
selected GPU
GPU VRAM
compute capability
```

Do not modify solver algorithms once the headline campaign starts.

If B1 uncovers a genuine correctness or benchmark-contract defect:

1. stop affected rows;
2. document the defect;
3. make only the minimum required repair;
4. freeze a new commit;
5. rerun all affected paired rows.

Do not mix headline rows from different frozen commits.

---

# 3. Critical pre-campaign gate A — SCF convergence semantics

SCF-B0C-RS observed a 40-iteration Fe block-recursion campaign with final
reported residuals approximately:

```text
CPU  ~8.64e-9
CUDA ~1.67e-8
```

while the executable convergence flag remained false.

Before any full-convergence B1 claim, determine exactly why.

Inspect the actual production convergence criterion and record:

```text
which residual/quantity controls convergence
which threshold is used
whether more than one criterion must pass
whether a minimum/maximum iteration rule applies
whether the harness field being called "residual" is the actual production criterion
```

Do not change the convergence criterion.

Do not reinterpret a small diagnostic residual as convergence if the production
SCF says otherwise.

If the harness reports the wrong convergence quantity, repair only the
instrumentation and rerun affected closure checks before B1.

B1 full-convergence rows are valid only when the normal production convergence
flag/criterion is satisfied.

---

# 4. Critical pre-campaign gate B — precision policy

This campaign must preserve a strict distinction between arithmetic precision
and physical correctness.

## 4.1 Reciprocal route

The primary reciprocal headline remains:

```text
CPU FP64
vs
GPU FP64
```

for equal-precision claims.

Existing FP32 eigensolver routes embedded in FP64 SCF state remain:

```text
numeric_mode = mixed
```

unless CURRENT HEAD genuinely provides a complete FP32 SCF route.

## 4.2 Real-space route

SCF-B0C-RS established that the current CUDA RS plugin uses:

```text
FP32 working paths
+
FP64 Hamiltonian/density/canonical SCF state
```

therefore the current production RS CUDA route is:

```text
numeric_mode = mixed
```

Do **not** label:

```text
CPU FP64 vs RS GPU mixed
```

as an equal-precision GPU speedup.

Before timing, inspect whether CURRENT HEAD already provides a genuinely
precision-matched CPU route using the same working precision for the RS kernel.

### If a matching CPU mixed/FP32-working route already exists

Benchmark it.

Then:

```text
CPU matched mixed
vs
GPU matched mixed
```

may produce an equal-numerical-mode RS accelerator speedup.

### If no matching CPU route exists

Do not implement a new solver merely for benchmark symmetry.

Report:

```text
CPU FP64 scientific reference
GPU mixed production route
```

with a clearly labelled production-route timing ratio, **not** an equal-precision
speedup.

The ratio may still be useful if physical correctness passes, but the report
must state explicitly that arithmetic precision differs.

## 4.3 No false FP32 story

Do not reuse the KPM result to imply that FP32 is automatically safe for SCF.

The KPM transport campaign is separate.

---

# 5. Speedup/ratio terminology

Use unambiguous terminology.

## Equal-precision / equal-mode CPU-GPU speedups

Only when the CPU and GPU numeric-mode contract matches:

```text
S_rs_kernel
S_rs_phase
S_solver
S_reciprocal
S_iteration
S_convergence
```

## Non-equal-precision production comparison

If RS CPU is FP64 and RS GPU is mixed, use clearly labelled fields such as:

```text
R_rs_kernel_production
R_rs_phase_production
R_iteration_production
R_convergence_production
```

Never put these ratios in the equal-precision headline table.

Every table must identify:

```text
CPU numeric mode
GPU numeric mode
equal_precision_eligible = true/false
```

---

# 6. Lane A — Real-space campaign

The RS campaign should focus on the production solver families that matter,
not mechanically benchmark every historical path.

## 6.1 Block recursion — mandatory primary RS GPU case

Use bcc Fe.

This is the central RS accelerator workload because SCF-B0C-RS verified the
production CUDA block path and no silent fallback.

Measure:

```text
P_rs_solver_kernel
RS electronic-structure phase
P_scf_iteration_total
full_scf_wall where converged
```

with the existing block Green reconstruction.

Record:

```text
rs_backend
nrec
block_size
terminator
system size
Natom
Hamiltonian dimension / sparse size metadata
nnz where meaningful
PBC
numeric mode
```

## 6.2 Chebyshev electronic-structure SCF

Use the validated Si fixture.

Benchmark CPU production Chebyshev SCF.

Benchmark GPU only if CURRENT HEAD has a production GPU Chebyshev SCF route.

If it does not:

```text
GPU = UNSUPPORTED
```

Do not implement one inside B1.

Record:

```text
Chebyshev order M
kernel/damping scheme
spectral bounds/scaling
system size
```

Keep this distinct from KPM/Kubo-Bastin transport.

## 6.3 Scalar/Lanczos recursion

Retain as a secondary CPU/reference route if it remains production-relevant.

Do not make it a major GPU axis unless a real production GPU route already
exists.

---

# 7. RS system-size ladder

Use existing validated production fixtures where possible.

For block-recursion Fe, benchmark a size ladder that exercises the sparse
kernel meaningfully.

Suggested shape, adjusted to CURRENT validated fixtures:

```text
small closure-sized Fe RS cell
4x4x4 or similar
6x6x6
8x8x8
10x10x10 if already practical/validated
```

Do not force exactly these sizes if the existing RS corpus uses a better ladder.

For each size record actual:

```text
Natom
Hamiltonian/sparse dimension
nnz
nrec
block size
```

Full SCF convergence is not mandatory at every large size.

Require:

```text
full convergence on at least one representative Fe RS workload
steady SCF iteration on several sizes
kernel/phase scaling on the largest practical workload
```

---

# 8. RS CPU fairness

For each important RS workload run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

where the production kernel supports OpenMP meaningfully.

Control BLAS/OpenMP settings and avoid oversubscription.

Retain all rows.

Select the best valid CPU row separately for:

```text
RS kernel
RS phase
SCF iteration
full convergence
```

because the best thread count may differ by level.

---

# 9. RS correctness requirements

For block recursion / Lanczos, retain structured kernel invariants and add
scientifically meaningful comparisons where already available:

```text
finite recursion coefficients
selected recursion coefficients
selected block/scalar Green values
DOS/LDOS samples
integrated spectral weight
electron count
```

For Chebyshev:

```text
selected moments
DOS/LDOS
integrated spectral weight
electron count
```

For converged RS SCF compare:

```text
total energy
Fermi energy
total charge
site charges
site magnetic moments where applicable
final production convergence quantity
```

A CUDA row is invalid if:

```text
rs_gpu_used != true
fallback_detected != false
kernel correctness != PASS
```

---

# 10. RS timing interpretation

SCF-B0C-RS established that currently:

```text
T_rs_kernel = P_rs_solver_kernel
```

and separate:

```text
T_rs_H2D
T_rs_D2H
T_rs_sync
T_rs_setup
```

are not exposed by the backend and remain zero with:

```text
rs_detail_timers_status = not_exposed_by_backend
```

Preserve this honest state.

Do not infer transfer time from a wall-time remainder.

Do not open a B1 optimization task merely to expose lower-level CUDA counters.

The outer production RS solver boundary is sufficient for application-level
performance conclusions.

---

# 11. Lane B — Reciprocal/k-space campaign

Retain the already-designed reciprocal B1 workload ladder.

## 11.1 Si primitive

Purpose:

```text
small matrix
nonmagnetic
many-k
GPU-overhead baseline
```

Run full SCF convergence.

## 11.2 Fe primitive

Purpose:

```text
small magnetic matrix
many-k
GPU-overhead baseline
```

Run full SCF convergence.

## 11.3 Fe 2x2x2

Approximate matrix dimension:

```text
~144
```

This is a key application crossover workload.

Attempt full CPU/GPU FP64 SCF convergence.

## 11.4 Fe 3x3x3

Approximate matrix dimension:

```text
~486
```

Attempt full convergence.

If impractical:

```text
same-starting-state steady SCF iteration campaign
```

with formal justification.

## 11.5 Fe 4x4x4

Approximate matrix dimension:

```text
~1152
```

Frozen-potential reciprocal/eigensolver benchmark is sufficient.

Optional SCF iterations only if practical.

## 11.6 Fe 5x5x5

Approximate matrix dimension:

```text
~2250
```

Use the validated frozen-potential reciprocal/eigensystem route.

Do not require full SCF.

Record memory limitations honestly.

Actual runtime dimensions are authoritative.

---

# 12. Reciprocal CPU fairness

For all headline reciprocal SCF workloads run:

```text
OMP 1
OMP 2
OMP 4
OMP 8
```

with controlled BLAS threading.

Compare GPU to the best valid CPU configuration.

Do not use OMP=1 as the sole CPU reference.

---

# 13. Reciprocal GPU strategies

Compare supported FP64 strategies where meaningful:

```text
Zheevd
ZheevjBatched
```

or their CURRENT HEAD names.

Do not force one solver strategy across all matrix sizes.

Use the best valid measured GPU strategy in the final workload map.

Retain unfavorable strategies in the evidence package.

---

# 14. Reciprocal timer-equivalence sanity check

Before the final timing matrix, verify on:

```text
Si primitive
Fe 2x2x2
```

that:

```text
CPU P_eigensolver
```

and:

```text
GPU T_solver
```

cover semantically equivalent numerical work.

Check:

```text
packing
workspace setup
H2D
D2H
postprocessing
```

If boundaries differ, repair reporting only.

Document the final timer contract.

---

# 15. Shared SCF timing

For both routes preserve the same route-neutral quantities:

```text
first_scf_iteration
steady_scf_iteration_median
P_scf_iteration_total
n_scf_iterations
full_scf_wall
final_convergence_residual/criterion
```

Require:

```text
profile closure <= 3%
```

for performance-eligible rows.

If:

```text
P_scf_misc > 5%
```

explain it.

---

# 16. Cold vs steady state

For both routes distinguish:

```text
cold process wall
backend/plugin initialization
first solver request
steady solver/phase
first SCF iteration
steady SCF iteration
full SCF convergence wall
```

Do not mix fresh-process CUDA startup with repeated SCF iteration cost.

Do not hide cold-start cost for one-shot workflows.

---

# 17. Full convergence protocol

For all full-convergence CPU/GPU comparisons use identical:

```text
material/structure
basis
spin/SOC state
starting state
mixing method
mixing parameters
smearing/electronic temperature
SCF convergence criterion
maximum-step policy
solver numerical controls
```

Do not loosen GPU convergence tolerances.

Different iteration counts are allowed if both converge to the same physical
solution.

---

# 18. Lane C — physically matched RS vs reciprocal comparison

This lane must be deliberately small.

Its purpose is not to prove one formulation universally superior.

It should answer:

> For selected periodic materials, at comparable physical accuracy, how do the
> real-space and reciprocal formulations compare in computational cost?

Use the `cross_route_case_id` / equivalent metadata introduced by SCF-B0C-RS.

## 18.1 Mandatory matched material — bcc Fe

Use one converged RS Fe block-recursion case and one converged reciprocal Fe
case.

Because the numerical convergence controls differ, do not require:

```text
same supercell size
same nominal matrix dimension
```

Instead establish an accuracy contract.

At minimum compare:

```text
total energy per atom
Fermi energy
magnetic moment per atom
total/site charge
SCF final residual
```

Optionally compare:

```text
DOS near EF
```

if already available robustly.

## 18.2 Optional matched material — Si

If the production RS Chebyshev SCF route is validated sufficiently, compare
with reciprocal Si.

This may remain CPU-vs-CPU if GPU Chebyshev SCF is unsupported.

---

# 19. Cross-route accuracy matching

Do not compare arbitrary RS and reciprocal settings.

For each cross-route case, construct a small numerical convergence study.

## RS controls may include

```text
supercell size
nrec
terminator
Chebyshev M
spectral broadening/kernel
```

## Reciprocal controls may include

```text
k mesh
Nk_unique
DOS/integration method
smearing/broadening
```

Choose production configurations such that both routes meet a common physical
accuracy target relative to the best practical reference.

Use existing scientific tolerances where already established.

If no defensible accuracy-matched pair can be constructed:

```text
cross-route performance comparison = INCONCLUSIVE
```

Do not force a winner.

---

# 20. Two tiers of cross-route performance comparison

Because RS GPU may be mixed precision, present cross-route comparisons in two
tiers.

## Tier 1 — arithmetic-clean scientific comparison

Prefer:

```text
RS CPU FP64
vs
reciprocal CPU FP64
```

and where a genuine FP64 GPU route exists:

```text
FP64 vs FP64
```

This establishes formulation-level cost without arithmetic-mode ambiguity.

## Tier 2 — best validated production configuration

Optionally compare:

```text
best validated RS production route
vs
best validated reciprocal production route
```

even if one route is mixed precision.

This must be clearly labelled:

```text
best_validated_production_configuration
```

and must include numerical-mode labels and physical correctness evidence.

Do not call this an equal-precision comparison.

---

# 21. Mandatory performance quantities

## RS lane

Where equal-mode pairing exists:

```text
S_rs_kernel
S_rs_phase
S_iteration
S_convergence
```

Where only FP64 CPU vs mixed GPU exists:

```text
R_rs_kernel_production
R_rs_phase_production
R_iteration_production
R_convergence_production
```

## Reciprocal lane

```text
S_solver
S_reciprocal
S_iteration
S_convergence
```

## Cross-route lane

Use direct wall-time comparisons only after accuracy matching.

Examples:

```text
T_RS_CPU / T_K_CPU
T_RS_best_validated / T_K_best_validated
```

Do not call these GPU speedups.

---

# 22. Timing protocol

For steady-state kernel/phase/iteration measurements:

```text
2 warmups
5 measured repetitions
```

where practical.

For expensive cases:

```text
minimum 3 measured repetitions
```

is acceptable with explicit annotation.

Use:

```text
median
min
max
MAD or IQR
```

Do not headline the fastest single sample.

For full SCF convergence:

- repeat key small/moderate cases at least 3 times if affordable;
- otherwise retain one validated convergence wall plus repeated steady-iteration timing.

---

# 23. Required RS tables

## Table RS-1 — block-recursion production performance

```text
material
system size
Natom
nnz/dimension
nrec
CPU mode
GPU mode
best CPU OMP
CPU kernel
GPU kernel
equal-mode S_rs_kernel OR production ratio
CPU phase
GPU phase
equal-mode S_rs_phase OR production ratio
correctness
```

## Table RS-2 — RS SCF iteration/convergence

```text
material
system size
solver
numeric modes
CPU iteration
GPU iteration
S_iteration or production ratio
CPU iterations-to-convergence
GPU iterations-to-convergence
CPU full wall
GPU full wall
S_convergence or production ratio
final-state correctness
```

## Table RS-3 — Chebyshev

```text
material
system size
M
CPU timing
GPU timing/status
correctness
notes
```

## Table RS-4 — scalar/Lanczos secondary evidence

Only if meaningful.

---

# 24. Required reciprocal tables

## Table K-1 — FP64 full SCF convergence

```text
material
supercell
Natom
nmat
Nk_unique
best CPU OMP
GPU strategy
CPU iterations
GPU iterations
CPU full wall
GPU full wall
S_convergence
correctness
```

## Table K-2 — steady SCF iteration

```text
material
supercell
nmat
Nk_unique
CPU iteration
GPU iteration
S_iteration
```

## Table K-3 — reciprocal phase

```text
material
nmat
Nk
CPU reciprocal
GPU reciprocal
S_reciprocal
GPU strategy
```

## Table K-4 — eigensolver

```text
material
nmat
Nk
CPU solver
GPU solver
S_solver
GPU strategy
memory
```

---

# 25. Required cross-route table

For each valid `cross_route_case_id`:

```text
material
RS solver
RS convergence controls
reciprocal convergence controls
RS numeric mode
reciprocal numeric mode
energy difference
Fermi difference
moment difference
charge difference
accuracy status
RS SCF iteration wall
reciprocal SCF iteration wall
RS full convergence wall
reciprocal full convergence wall
comparison tier
interpretation
```

No cross-route performance number is valid if:

```text
accuracy status != PASS
```

---

# 26. Required plots

Generate at least:

## RS plots

1. RS kernel time / production ratio vs system size.
2. RS SCF iteration time vs system size.
3. RS SCF stage fractions for representative Fe workloads.
4. CPU OMP scaling for the block-recursion route.

## Reciprocal plots

5. `S_solver` vs matrix dimension.
6. `S_reciprocal` vs matrix dimension.
7. `S_iteration` vs matrix dimension.
8. `S_convergence` vs matrix dimension where available.
9. reciprocal SCF stage fractions.

## Shared / comparison plots

10. SCF residual/criterion vs iteration for representative CPU/GPU RS and reciprocal cases.
11. Optional `nmat x Nk` reciprocal crossover map.
12. Accuracy-vs-wall-time plot for the selected RS-vs-k-space comparison if enough convergence-study points exist.

Do not force a plot with too few meaningful points.

---

# 27. Required interpretation questions — RS

Answer explicitly:

1. Does the production CUDA block-recursion route reduce kernel wall time?
2. Is that an equal-precision result or a mixed-vs-FP64 production ratio?
3. How much of the kernel reduction survives the RS electronic-structure phase?
4. How much survives the complete SCF iteration?
5. Does the GPU RS route reach the same converged physical state?
6. Does the GPU alter the iteration count?
7. Which RS phase dominates after acceleration?
8. How does CPU OpenMP scale for block recursion?
9. Is Chebyshev GPU SCF supported?
10. Is scalar/Lanczos worth retaining as a performance route or mainly reference functionality?
11. What RS workload sizes are clearly CPU-preferred, GPU-beneficial, or unresolved?
12. Is further RS kernel optimization justified?

---

# 28. Required interpretation questions — reciprocal

Answer explicitly:

1. Does GPU FP64 help primitive Si?
2. Does GPU FP64 help primitive Fe?
3. Where does eigensolver crossover occur?
4. Where does reciprocal-phase crossover occur?
5. Where does complete SCF-iteration crossover occur?
6. Does Fe 2x2x2 achieve full-SCF crossover?
7. Does Fe 3x3x3 strengthen the GPU case?
8. How much eigensolver acceleration survives at SCF level?
9. Which non-eigensolver stage becomes dominant?
10. Which GPU solver strategy is preferred in each regime?
11. How important is cold CUDA startup?
12. Is further reciprocal SCF optimization justified?

---

# 29. Required interpretation questions — RS vs reciprocal

Answer explicitly:

1. Can a defensible physically matched Fe comparison be constructed?
2. At matched accuracy, which CPU formulation is faster?
3. At matched accuracy, which validated production configuration is faster?
4. Is the answer sensitive to the chosen accuracy target?
5. Does RS gain relative advantage as the physical system becomes larger or less translationally simple?
6. Does reciprocal space remain superior for small perfect primitive crystals?
7. Which conclusions are measured, and which remain architectural expectations rather than benchmark evidence?

Do not generalize beyond measured fixtures.

---

# 30. No cherry-picking

Preserve:

```text
GPU losses
non-convergence
mixed-precision limitations
unsupported GPU routes
poor OpenMP scaling
large-memory skips
batched-solver losses
cross-route inconclusive cases
```

A negative result is valid evidence.

Do not remove a row because it weakens the desired narrative.

---

# 31. Recommended-backend policy

Do not change production defaults in B1.

The report should end with an evidence-based support table, for example:

```text
workload regime
RS CPU
RS GPU
k-space CPU
k-space GPU
best measured route
confidence/caveat
```

Possible regime descriptors may include:

```text
small perfect primitive cell
many-k small matrix
moderate periodic supercell
large sparse real-space system
large dense reciprocal eigensystem
broken translational symmetry / disorder
```

Only mark:

```text
best measured route
```

where direct evidence supports it.

For architectural cases not directly benchmarked, label:

```text
expected / not measured
```

rather than presenting them as performance results.

---

# 32. Output package

Write under:

```text
results/benchmarks/scf_b1/
```

using the shared canonical pipeline:

```text
campaign.json
campaign.csv
campaign.md
raw/
correctness/
iteration_history/
plots/
```

All summary tables and plots must derive from canonical data.

Do not manually transcribe timing results into Markdown.

---

# 33. Final report

Write:

```text
docs/dev/RS_LMTO_ASA_SCF_B1_REPORT.md
```

Recommended structure:

## 1. Executive summary

Include:

```text
frozen commit
hardware
main RS conclusion
main reciprocal conclusion
cross-route conclusion
precision caveat
```

## 2. Methodology

Document:

```text
three-lane structure
CPU fairness
precision policy
correctness
starting states
timing protocol
convergence semantics
```

## 3. Real-space workloads and solver taxonomy

## 4. Real-space kernel results

## 5. Real-space SCF iteration/convergence results

## 6. Reciprocal workloads

## 7. Reciprocal eigensolver results

## 8. Reciprocal SCF iteration/convergence results

## 9. CPU OpenMP scaling

## 10. Precision / mixed-route evidence

Explicitly distinguish:

```text
equal-precision speedups
production mixed-vs-FP64 ratios
```

## 11. Physically matched RS vs reciprocal comparison

## 12. Stage/Amdahl analysis

## 13. GPU memory / startup limitations

## 14. Recommended workload map

## 15. Remaining optimization opportunities

## 16. Final performance conclusion

State whether the SCF accelerator campaign can close.

---

# 34. Stop condition

SCF-B1 is complete when:

```text
implementation is frozen;

convergence semantics are verified;

RS block-recursion CPU/GPU production evidence exists;

RS precision status is handled honestly;

at least one RS workload reaches full valid convergence on CPU and the relevant
GPU route if possible;

Si Chebyshev CPU evidence exists and GPU support status is explicit;

reciprocal Si full convergence is benchmarked;

reciprocal Fe primitive full convergence is benchmarked;

reciprocal Fe 2x2x2 full convergence is benchmarked or explicitly impractical;

Fe 3x3x3 has full convergence or justified steady-iteration evidence;

Fe 4x4x4 large reciprocal evidence exists;

Fe 5x5x5 evidence exists or has an explicit skip;

CPU OMP fairness is complete;

all equal-precision headline rows pass correctness;

mixed/non-equal-precision comparisons are separated;

a physically matched Fe RS-vs-k-space comparison is attempted;

cross-route claims are made only when accuracy PASS;

plots and tables are generated;

backend/workload recommendations are written;

negative results are retained;

remaining limitations are explicit.
```

Do not invent SCF-B2 automatically.

A follow-up should be proposed only if B1 identifies one concrete, high-value
correctness or performance problem.

---

# 35. Checklist

## Freeze and preflight

- [ ] clean tracked implementation
- [ ] frozen commit recorded
- [ ] environment recorded
- [ ] production convergence criterion identified
- [ ] false/non-obvious convergence-flag issue explained
- [ ] harness convergence field verified
- [ ] RS precision path audited
- [ ] reciprocal precision path audited
- [ ] equal-precision eligibility rules verified

## RS route

- [ ] block recursion production fixture selected
- [ ] Fe size ladder selected
- [ ] CPU OMP 1/2/4/8
- [ ] CUDA block route
- [ ] no fallback
- [ ] kernel correctness PASS
- [ ] RS phase correctness PASS
- [ ] steady iteration measured
- [ ] full convergence attempted
- [ ] final physical state checked
- [ ] equal-mode speedup OR production ratio labelled correctly
- [ ] stage fractions generated

## Chebyshev RS

- [ ] Si production fixture
- [ ] M recorded
- [ ] scaling/bounds recorded
- [ ] CPU route measured
- [ ] GPU route measured if supported
- [ ] UNSUPPORTED retained otherwise
- [ ] correctness PASS
- [ ] distinct from KPM transport

## Lanczos/scalar

- [ ] production/reference status stated
- [ ] secondary evidence retained if useful
- [ ] no artificial GPU task added

## Reciprocal Si

- [ ] CPU OMP sweep
- [ ] GPU FP64
- [ ] full convergence
- [ ] solver/phase/iteration/convergence metrics
- [ ] correctness PASS

## Reciprocal Fe primitive

- [ ] CPU OMP sweep
- [ ] GPU FP64
- [ ] full convergence
- [ ] magnetic correctness
- [ ] all four metrics

## Reciprocal Fe2

- [ ] actual nmat/Nk
- [ ] CPU OMP sweep
- [ ] GPU FP64
- [ ] full convergence attempted
- [ ] correctness
- [ ] all available metrics

## Reciprocal Fe3

- [ ] actual nmat/Nk
- [ ] CPU OMP sweep
- [ ] GPU FP64
- [ ] full convergence attempted
- [ ] convergence OR formal steady substitute
- [ ] correctness

## Reciprocal Fe4/Fe5

- [ ] real production frozen-potential rows
- [ ] actual dimensions
- [ ] CPU solver
- [ ] GPU solver
- [ ] reciprocal phase
- [ ] memory
- [ ] correctness
- [ ] explicit skips if necessary

## Cross-route Fe

- [ ] shared physical case ID
- [ ] RS convergence study
- [ ] reciprocal convergence study
- [ ] common accuracy target
- [ ] total energy/atom comparison
- [ ] Fermi comparison
- [ ] moment comparison
- [ ] charge comparison
- [ ] accuracy PASS/FAIL
- [ ] CPU RS vs CPU reciprocal timing if PASS
- [ ] best validated production comparison if PASS
- [ ] arithmetic modes shown
- [ ] no forced conclusion if INCONCLUSIVE

## Timing/profiling

- [ ] RS kernel boundary fair
- [ ] reciprocal solver boundary fair
- [ ] cold vs steady separated
- [ ] profile closure <=3%
- [ ] misc <=5% or explained
- [ ] no inferred unexposed CUDA RS transfer timings

## Correctness

- [ ] solver/kernel correctness
- [ ] Fermi energy
- [ ] total charge
- [ ] site charges
- [ ] magnetic moments
- [ ] total energy
- [ ] production convergence quantity
- [ ] iteration-count differences retained
- [ ] no silent fallback

## Outputs

- [ ] canonical JSON
- [ ] CSV
- [ ] Markdown
- [ ] raw logs
- [ ] correctness evidence
- [ ] iteration histories
- [ ] plots
- [ ] final unified B1 report
- [ ] steering/master status updated

## Final interpretation

- [ ] RS recommendation
- [ ] reciprocal recommendation
- [ ] precision caveat
- [ ] cross-route conclusion
- [ ] workload map
- [ ] remaining bottleneck
- [ ] further-work decision
- [ ] SCF campaign closed if justified

---

# 36. Commit message

Use:

```text
Publish unified RS and reciprocal SCF benchmark campaign
```

---

# 37. Final mindset

SCF-B1 is not:

```text
make CUDA win
```

and it is not:

```text
prove RS is better than k-space
```

It is:

```text
measure the production formulations honestly
```

The final report must distinguish four separate questions:

```text
1. How fast is the RS numerical kernel?
2. How fast is the reciprocal eigensolver?
3. How much of either acceleration reaches the complete SCF calculation?
4. At comparable physical accuracy, which formulation is preferable for the
   measured workload?
```

The strongest result is not necessarily the largest speedup.

The strongest result is a defensible workload map showing where each method
belongs, with numerical correctness and precision limitations made explicit.
