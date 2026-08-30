# SCF-B1R2 — Chebyshev GPU scaling and final SCF performance closure

You are working on the **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This is the final targeted performance-closeout task following:

```text
SCF-B0C
SCF-B0C-RS
SCF-B1
SCF-B1R
```

It is **not** a new optimization campaign.

It is also **not** the RS-vs-k-space scientific-equivalence task. The failed
fixed-potential RS-vs-k-space accuracy comparison will be handled separately
later as a validation task.

The purpose of SCF-B1R2 is to close the remaining performance gaps:

```text
A. properly benchmark RS Chebyshev CPU vs GPU, including size scaling;
B. complete common-state correctness for reciprocal Fe2/Fe3 iteration timing;
C. suppress or minimally repair the invalid K1 CPU eigensolver component timer;
D. freeze clean source provenance for the final campaign;
E. reconcile the final SCF performance report.
```

The central performance metric remains:

```text
steady production SCF iteration wall time
```

with route-specific kernel/phase decomposition.

SCF iteration count is a mixing/minimization/convergence diagnostic, **not** an
accelerator performance metric.

---

# 1. Current evidence and why this task exists

The reconciled B1R campaign established:

## Real-space block recursion

For bcc-Fe 3x3x3 (`nmat=486`, `nrec=20`), the fair CPU OMP sweep gave
approximately:

```text
CPU OMP1: kernel 198.1 s, iteration 222.2 s
CPU OMP2: kernel 121.0 s, iteration 140.4 s
CPU OMP4: kernel  89.1 s, iteration 100.1 s
CPU OMP8: kernel  69.7 s, iteration  78.3 s
CUDA mixed: kernel 73.4 s, iteration 78.0 s
```

Thus the current block-recursion CUDA route is approximately at parity with the
best measured threaded CPU at this size.

Because the CUDA RS route uses FP32 working paths with FP64 canonical SCF state,
it remains:

```text
numeric_mode = mixed
```

and must not be presented as an equal-precision FP64 speedup.

## Large reciprocal FP64 eigensystems

Current-build frozen-potential/eigenvector-enabled rows established approximately:

```text
nmat=486  -> 1.49x CPU/GPU
nmat=1152 -> 3.24x
nmat=2250 -> 4.28x
```

for FP64 dense eigensolution.

These are valid current-build solver results.

## Reciprocal full iteration

The reconciled Fe2/Fe3 rows showed approximately:

```text
Fe2 nmat=144, Nk=512:
CPU iteration ~3.43 s
GPU iteration ~5.63 s

Fe3 nmat=486, Nk=512:
CPU iteration ~45.74 s
GPU iteration ~27.12 s
```

but the GPU rows do not yet have the common-state correctness pairing needed to
promote these application-level iteration ratios.

## Chebyshev

Earlier B1 evidence contained a very promising RS Chebyshev result for diamond
Si:

```text
CPU OMP1 iteration ~6.56 s
CPU OMP2 iteration ~3.60 s
CPU OMP4 iteration ~2.27 s
CPU OMP8 iteration ~1.70 s
CUDA mixed iteration ~0.527 s
```

but the CUDA row failed the correctness/convergence eligibility gate.

Therefore the raw approximately 3.2x ratio versus CPU OMP8 is **not yet a valid
performance claim**.

Chebyshev is expected architecturally to be more GPU-friendly than block
recursion because its dominant recurrence can have a more regular sparse
matrix-vector / polynomial structure with less complex block-recursion
bookkeeping. Treat this only as a **hypothesis to test**, not as a conclusion.

SCF-B1R2 must establish whether that expected scaling is actually present in the
production implementation.

---

# 2. User-run execution workflow — mandatory

The user does not want Luna spending tokens monitoring long benchmark jobs.

Use the same two-phase pattern as B1R.

## PHASE A — prepare

Luna should:

1. inspect the current B1R canonical evidence and harness;
2. implement only the small benchmark/report changes required below;
3. create one resume-safe bench-all command;
4. run only short unit/parser/smoke tests;
5. print the exact command for the user;
6. STOP.

Do not launch the long benchmark matrix.

Do not poll or monitor it.

Preferred command:

```bash
bash tests/benchmarks/run_scf_b1r2_all.sh
```

## PHASE B — after the user reports completion

Luna should:

1. ingest the completed result directory;
2. verify all expected cases;
3. rerun nothing expensive unless evidence is missing/corrupt;
4. compute canonical ratios;
5. regenerate affected plots/tables;
6. reconcile the final SCF report;
7. write a short B1R2 closeout;
8. state whether the SCF **performance** campaign can close.

---

# 3. Clean provenance before benchmarking

B1R verified that the actual Release compile commands use effective `-O3`, but
the tested source state was still reported as dirty.

Before the B1R2 long campaign:

1. inspect tracked changes;
2. commit the source/build-system fixes and benchmark-harness changes required
   for B1R2;
3. require:

```bash
git diff --quiet
git diff --cached --quiet
```

for tracked source before building the final benchmark executable;

4. record the exact new commit;
5. build the benchmark executable from that clean tracked commit.

Untracked machine-local build/results directories do not need to be committed,
but distinguish:

```text
git_dirty_tracked
git_untracked_present
```

in provenance.

Do not benchmark a numerically different dirty working tree and later pretend a
commit is equivalent.

Use the clean B1R2 commit as the final performance provenance.

---

# 4. Lane A — RS Chebyshev CPU vs GPU is mandatory

This is the main new performance task.

Use the production:

```text
real_space --rs-solvers chebyshev
```

route and the validated diamond-Si fixture family.

Do not confuse this with KPM/Kubo-Bastin transport.

The benchmark target is:

```text
Chebyshev recurrence/kernel
    ->
Chebyshev electronic-structure phase
    ->
steady SCF iteration
```

The full number of SCF cycles to convergence is not the performance metric.

---

# 5. Chebyshev numerical configuration

Inspect CURRENT HEAD and use the validated production settings.

The earlier validated campaign used approximately:

```text
M = 200
```

for Si.

Use the same production:

```text
Chebyshev order M
kernel/damping scheme
spectral bounds/scaling policy
energy grid/integration policy
smearing/broadening
```

for paired CPU/GPU rows.

These controls are part of the comparison key.

Do not alter M or spectral resolution differently on CPU and GPU to improve
timing.

If B1R2 includes an M-scaling diagnostic, it must be a **separate secondary
study**, not mixed into the backend speedup comparison.

---

# 6. Chebyshev size-scaling ladder

The earlier result is insufficient because it is essentially one fixture.

Benchmark a real periodic Si size ladder large enough to reveal CPU/GPU scaling.

Prefer existing validated fixtures.

If no canonical ladder already exists, stage a modest production ladder such
as:

```text
Si 4x4x4
Si 6x6x6
Si 8x8x8
Si 10x10x10 or 12x12x12 if practical
```

The exact sizes should be chosen from CURRENT production capability and runtime
cost.

Do not mechanically force the largest size if memory/runtime is unreasonable.

For every row record actual:

```text
replication
Natom
Hamiltonian/sparse dimension
nnz
M
spectral bounds
backend
numeric mode
```

The objective is to establish:

```text
CPU-preferred regime
crossover/parity regime
GPU-beneficial regime
```

if those regimes are observed.

Do not assume they must exist.

---

# 7. Chebyshev CPU fairness

For every important size run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

where production OpenMP applies.

Control BLAS threading and avoid oversubscription.

Retain every row.

Select independently:

```text
best CPU Chebyshev kernel
best CPU Chebyshev phase
best CPU SCF iteration
```

because the optimal OMP count can differ.

Do not compare GPU only against CPU OMP1.

---

# 8. Chebyshev GPU precision policy

Audit the actual CUDA Chebyshev numerical path before interpreting results.

Record independently:

```text
Hamiltonian precision
Chebyshev recurrence precision
moment/coefficient precision
spectral reconstruction precision
density precision
canonical SCF precision
```

If, as expected, the GPU path is FP32-working + FP64 canonical state:

```text
numeric_mode = mixed
```

and CPU FP64 vs GPU mixed is a:

```text
production-route ratio
```

not an equal-precision speedup.

If a matching CPU mixed/FP32-working Chebyshev route already exists, benchmark
it as an additional same-mode comparison.

Do not implement a new precision path merely for benchmark symmetry.

---

# 9. Chebyshev correctness gate — use a common state, not SCF cycle count

The earlier CUDA Chebyshev row failed eligibility. Resolve whether this was:

```text
a numerical electronic-structure mismatch
or
merely failure to satisfy the full SCF convergence budget
```

without turning B1R2 into a mixing study.

## Preferred correctness protocol

Use the same validated Si:

```text
pot.par
restart
or fixed potential/state
```

for CPU and GPU.

From that exact common state, run either:

```text
one fixed-potential electronic-structure evaluation
```

or, if the production architecture requires it:

```text
one / a small fixed number of identical SCF iterations
```

before mixing history can diverge substantially.

Compare:

```text
Fermi energy
total electron count / charge
site charge
stable electronic/band energy quantity
DOS / LDOS where available
selected Chebyshev moments where available
density norm where available
```

For nonmagnetic Si, magnetic moment is not a primary discriminator.

## DOS comparison

Strongly prefer a DOS correctness check using the same:

```text
energy grid
energy window
broadening/kernel
normalization
```

Report:

```text
relative L2 DOS difference
max absolute DOS difference
integrated DOS/electron-count difference
```

Do not compare curves with different broadening.

## Eligibility

If CPU and GPU common-state electronic observables pass the established
scientific tolerance contract, the steady-iteration production ratio may be
reported even if a long SCF run would require different cycle counts.

If the electronic observables fail:

```text
Chebyshev GPU performance = measured but scientifically ineligible
```

and the discrepancy becomes a separate correctness/precision follow-up.

Do not fix the algorithm inside B1R2.

---

# 10. Chebyshev timing quantities

Measure separately:

```text
P_rs_solver_kernel
P_rs_spectral_reconstruct
P_rs_fermi
P_rs_density_build
RS electronic-structure phase
P_scf_iteration_total
```

using existing production boundaries.

If lower-level CUDA transfer/setup timers remain unexposed:

```text
rs_detail_timers_status = not_exposed_by_backend
```

and do not infer them.

For every valid paired size report either:

```text
S_cheb_kernel
S_cheb_phase
S_iteration
```

for equal-mode comparisons,

or:

```text
R_cheb_kernel_production
R_cheb_phase_production
R_iteration_production
```

for FP64 CPU vs mixed GPU.

---

# 11. Chebyshev timing protocol

For small/moderate rows:

```text
2 warmups
5 measured repetitions
```

For the largest expensive rows:

```text
1-2 warmups
3 measured repetitions
```

is acceptable.

Use:

```text
median
min
max
MAD or IQR
```

Do not use the fastest sample as the headline.

---

# 12. Required Chebyshev tables

## Table C1 — CPU scaling

```text
Si size
Natom
nnz/dimension
M
OMP
kernel s
phase s
iteration s
status
```

## Table C2 — CPU/GPU production comparison

```text
Si size
Natom
M
CPU mode
GPU mode
best CPU OMP
CPU kernel
GPU kernel
kernel ratio
CPU phase
GPU phase
phase ratio
CPU iteration
GPU iteration
iteration ratio
correctness status
equal_precision_eligible
reason
```

## Table C3 — correctness

```text
Si size
common state ID
Fermi difference
charge difference
energy difference
DOS relative L2
DOS max abs
integrated DOS difference
status
reason
```

---

# 13. Required Chebyshev plots

Generate:

## Plot C1

```text
Chebyshev kernel wall time vs system size
```

CPU best-threaded and GPU.

## Plot C2

```text
Chebyshev SCF iteration wall time vs system size
```

CPU best-threaded and GPU.

## Plot C3

```text
Chebyshev production ratio vs system size
```

with a horizontal parity reference at 1.

Label clearly that it is a mixed-vs-FP64 production ratio if that is the actual
precision configuration.

## Plot C4

For one representative correctness case:

```text
CPU vs GPU DOS
```

only if the grids/broadening match.

---

# 14. Lane B — reciprocal Fe2/Fe3 common-state correctness

Do not rerun a broad reciprocal campaign.

Use the existing Fe2 and Fe3 production configurations:

```text
Fe2: nmat=144, Nk_unique=512
Fe3: nmat=486, Nk_unique=512
```

with:

```text
CPU FP64
GPU FP64
```

from the exact same potential/restart state.

Run a small fixed number of iterations sufficient to validate per-iteration
electronic correctness.

Compare at minimum:

```text
Fermi energy
total charge
site charge
magnetic moment
stable total/band-energy quantity
density/residual quantity
```

Do not require the CPU and GPU to take the same number of cycles to eventual
SCF convergence.

If correctness passes, promote the existing/new steady iteration comparison to
an eligible application-level result.

The expected current timing values are approximately:

```text
Fe2 CPU ~3.43 s, GPU ~5.63 s
Fe3 CPU ~45.74 s, GPU ~27.12 s
```

but B1R2 must use the clean frozen B1R2 build and canonical results rather than
copying these historical values.

---

# 15. K1 eigensolver component timer — suppress or minimally repair

The B1R K1 CPU `P_eigensolver` values are not physically credible for:

```text
nmat=18
nmat=144
nmat=486
```

because they remain at the ~10^-5 s scale even though the reciprocal basis is
genuinely enlarged.

Do not spend a large development effort on this.

## Preferred action

Inspect the accumulator/boundary once.

If the bug is trivial, e.g.:

```text
timer overwritten instead of accumulated
timer scoped inside only the final k point
wrong accumulator field
```

repair it and add a focused regression test.

## Otherwise

Suppress:

```text
K1 S_solver
```

and mark:

```text
solver_component_timer_status = INVALID_FOR_K1
```

with reason.

Use the already validated K2 frozen-potential campaign as the authoritative
solver-level evidence.

K1 remains authoritative for:

```text
complete reciprocal electronic-structure / SCF iteration timing
```

once common-state correctness passes.

Do not block B1R2 closure on repairing a redundant component timer.

---

# 16. Keep K2 large reciprocal evidence as authoritative solver benchmark

Retain the current-build FP64/eigenvector K2 rows:

```text
nmat=486
nmat=1152
nmat=2250
```

as the authoritative dense-eigensolver result.

Do not rerun them unless:

```text
the clean B1R2 source commit changes numerical/solver code
or
the build differs materially from the validated B1R build
```

If only harness/report code changes, preserve the current measured data with
explicit source/build provenance as allowed by the repository's evidence
policy.

If source/build changes require rerun, include K2 in the bench-all script.

---

# 17. RS block-recursion result is already closed

Do not rerun the broad block-recursion campaign unless the clean-source rebuild
invalidates the old B1R timing evidence.

The current scientific performance conclusion should remain unless new evidence
changes it:

```text
At Fe3, mixed CUDA block recursion is approximately at parity with the best
measured 8-thread FP64 CPU route.
```

This is not an equal-precision speedup.

The purpose of B1R2 is to determine whether **Chebyshev**, as expected, shows a
more favorable GPU scaling trend.

---

# 18. RS-vs-k-space formulation equivalence is explicitly deferred

Do not attempt to repair the failed Fe RS-vs-k-space fixed-potential accuracy
contract in this task.

Record the existing conclusion:

```text
RS-vs-k-space formulation comparison = INCONCLUSIVE
```

The differences in energy/Fermi/moment require a dedicated validation study.

That later task should investigate:

```text
RS cluster size
recursion depth
terminator
Chebyshev order where relevant
k mesh
DOS convergence
common pot.par
```

but that is outside SCF-B1R2.

---

# 19. Bench-all script

Create:

```text
tests/benchmarks/run_scf_b1r2_all.sh
```

or an equivalent single-command Python driver.

It must be:

```text
resume-safe
self-logging
case-manifest driven
non-interactive
```

and should print progress like:

```text
[B1R2 01/N] Chebyshev Si4 CPU OMP1
[B1R2 02/N] Chebyshev Si4 CPU OMP2
...
[B1R2 XX/N] Reciprocal Fe3 GPU common-state
```

Recommended result directory:

```text
results/benchmarks/scf_b1r2/
```

Support:

```text
--resume
--force
--dry-run
```

if easy.

At completion print:

```text
B1R2_RUN_COMPLETE
```

and summarize completed/failed/skipped cases.

Luna must give the user exactly one command and then stop.

---

# 20. Canonical outputs

Write:

```text
results/benchmarks/scf_b1r2/
```

with:

```text
manifest.json
campaign.json
campaign.csv
campaign.md
raw/
correctness/
iteration_history/
chebyshev/
reciprocal_common_state/
plots/
```

Reuse the B1/B1R schema rather than creating an incompatible format.

---

# 21. Final report reconciliation

Update:

```text
docs/dev/RS_LMTO_ASA_SCF_B1_REPORT.md
```

and write a concise addendum:

```text
docs/dev/RS_LMTO_ASA_SCF_B1R2_REPORT.md
```

The final SCF performance report should contain distinct sections:

```text
1. build/provenance
2. RS block recursion
3. RS Chebyshev
4. reciprocal SCF iteration
5. large reciprocal FP64 eigensolver
6. CPU OpenMP scaling
7. precision caveats
8. final workload map
9. deferred RS-vs-k-space validation issue
```

---

# 22. Final workload map

The final report should be able to state, from measured evidence:

```text
RS block recursion:
    CPU vs mixed GPU scaling

RS Chebyshev:
    CPU vs mixed/equal-mode GPU scaling

reciprocal small matrices:
    CPU/GPU iteration result

reciprocal moderate matrix (Fe3):
    CPU/GPU iteration result if common-state correctness passes

large reciprocal dense eigensystem:
    FP64 GPU solver scaling
```

Do not extrapolate beyond measured sizes.

For each route identify:

```text
CPU-preferred
parity
GPU-beneficial
inconclusive
unsupported
```

with precision caveats.

---

# 23. Required final interpretation questions

## Chebyshev

1. Does Chebyshev GPU outperform the best-threaded CPU at the smallest tested
   size?
2. Where is the crossover?
3. How does the production ratio evolve with system size?
4. Does Chebyshev scale more favorably on GPU than block recursion in the
   measured range?
5. Is the gain mainly in the recurrence kernel or does it survive spectral
   reconstruction and the full SCF iteration?
6. Does the common-state CPU/GPU electronic structure pass correctness?
7. Is the current CUDA route mixed precision?
8. Is Chebyshev GPU ready for production SCF iteration work?

## Reciprocal

9. Does Fe2 common-state correctness pass?
10. Does Fe3 common-state correctness pass?
11. If yes, what are the valid steady-iteration CPU/GPU ratios?
12. Does Fe3 show an application-level GPU crossover consistent with the K2
    solver crossover?
13. Is the K1 solver component timer fixed or explicitly suppressed?

## Overall

14. Which current SCF route has the strongest measured GPU scaling:
    block recursion, Chebyshev, or large reciprocal diagonalization?
15. Which conclusions are equal-precision and which are mixed-vs-FP64
    production comparisons?
16. Can the scoped SCF performance campaign now close?

---

# 24. Acceptance criteria

SCF-B1R2 is complete only when:

- [ ] final benchmark source is tied to a clean tracked commit;
- [ ] one-command bench-all driver exists;
- [ ] Luna does not monitor long jobs;
- [ ] Chebyshev Si size ladder is measured;
- [ ] CPU OMP1/2/4/8 is measured for important Chebyshev sizes;
- [ ] CUDA Chebyshev route is measured where supported;
- [ ] Chebyshev precision is classified honestly;
- [ ] CPU/GPU common-state Chebyshev correctness is evaluated;
- [ ] DOS correctness is evaluated on at least one representative paired case
      if production output supports it;
- [ ] Chebyshev kernel, phase, and iteration ratios are reported;
- [ ] Chebyshev size scaling is plotted;
- [ ] block-recursion result is retained as comparison context;
- [ ] reciprocal Fe2 common-state correctness is evaluated;
- [ ] reciprocal Fe3 common-state correctness is evaluated;
- [ ] eligible reciprocal iteration ratios are promoted only after correctness
      PASS;
- [ ] invalid K1 solver timer is fixed trivially or explicitly suppressed;
- [ ] K2 remains the authoritative dense-solver benchmark;
- [ ] RS-vs-k-space formulation equivalence remains explicitly deferred;
- [ ] final workload map is updated;
- [ ] final report distinguishes equal-precision from mixed production ratios;
- [ ] no new optimization campaign is opened.

---

# 25. PHASE A checklist

Before giving the command to the user:

- [ ] inspect B1R canonical results
- [ ] inspect Chebyshev production CPU/GPU code paths
- [ ] determine actual GPU precision
- [ ] select validated Si size ladder
- [ ] prepare common-state correctness fixture
- [ ] inspect K1 eigensolver timer once
- [ ] fix or suppress K1 component timer
- [ ] prepare reciprocal Fe2/Fe3 common-state cases
- [ ] commit benchmark/source changes
- [ ] verify tracked source clean
- [ ] build clean Release executable
- [ ] verify effective optimization
- [ ] create bench-all script
- [ ] dry-run script
- [ ] run only short tests
- [ ] print exact one-line user command
- [ ] STOP

---

# 26. PHASE B checklist

After the user says the script completed:

- [ ] load manifest
- [ ] verify all expected cases
- [ ] select best CPU OMP rows
- [ ] compute Chebyshev production ratios
- [ ] evaluate Chebyshev correctness
- [ ] evaluate DOS agreement
- [ ] evaluate Fe2/Fe3 reciprocal correctness
- [ ] promote only eligible iteration ratios
- [ ] confirm K1 solver timer status
- [ ] regenerate affected plots
- [ ] update final B1 report
- [ ] write B1R2 closeout
- [ ] state final scoped performance conclusions
- [ ] state whether SCF performance campaign closes

---

# 27. Commit messages

Phase A:

```text
Prepare final SCF performance closure
```

Phase B:

```text
Close final SCF performance campaign
```

If one combined final commit is preferred:

```text
Finalize SCF CPU GPU benchmark campaign
```

---

# 28. Final mindset

The missing benchmark is not more block recursion.

The missing benchmark is:

```text
RS Chebyshev CPU vs GPU
```

with:

```text
fair CPU threading
real size scaling
common-state correctness
honest precision labels
kernel -> phase -> iteration accounting
```

The working hypothesis is that Chebyshev should expose more favorable GPU
scaling than block recursion because its dominant recurrence is more regular.

B1R2 must test that hypothesis.

Do not encode the expected answer into the report.

The final performance picture should be able to distinguish:

```text
block recursion:
    approximately CPU/GPU parity at the largest measured Fe3 workload

Chebyshev:
    measured size-dependent CPU/GPU behavior from B1R2

reciprocal dense eigensolution:
    clear FP64 GPU advantage for sufficiently large matrices

complete reciprocal SCF iteration:
    promoted only after common-state correctness passes
```

After that, the SCF **performance** campaign should be closed.

The separate fixed-potential RS-vs-k-space electronic-structure mismatch belongs
to a later scientific validation task, not another performance benchmark.
