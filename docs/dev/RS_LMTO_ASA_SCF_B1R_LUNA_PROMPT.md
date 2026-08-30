# SCF-B1R — Targeted reconciliation of the unified SCF benchmark campaign

Repository / branch:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This is a **small benchmark-reconciliation task**, not a new optimization work package.
It follows:

```text
SCF-B0C
SCF-B0C-RS
SCF-B1
```

The current final-campaign report is:

```text
docs/dev/RS_LMTO_ASA_SCF_B1_REPORT.md
```

and the canonical B1 evidence is under:

```text
results/benchmarks/scf_b1/
```

The purpose of B1R is to resolve seven specific issues before closing the SCF performance campaign:

1. verify whether the suspicious trailing `-O0` actually contaminated the benchmark build;
2. resolve the Fe2/Fe3 reciprocal `nmat`, `Nk`, and eigensolver-timer semantics;
3. complete the missing Fe3 real-space CPU OpenMP sweep;
4. establish one trustworthy RS Fe correctness/reference case, while treating **SCF iteration time rather than number of SCF cycles** as the performance quantity;
5. rerun the Fe 3^3 / 4^3 / 5^3 reciprocal frozen-potential eigensystem rows on the current validated build;
6. attempt one careful FP64 CPU RS-vs-k-space Fe comparison using a common potential / `pot.par` and, preferably, a common DOS comparison;
7. repair the final report tables so solver strategy and failure/ineligibility reasons are explicit.

Do not start a new solver/kernel optimization campaign inside B1R.

---

# 1. Performance philosophy

For this project the primary SCF accelerator metric is:

```text
steady production SCF iteration wall time
```

and its route-specific decomposition.

The number of SCF iterations required to reach convergence is primarily governed by:

```text
mixing
minimization
starting state
physical convergence behavior
```

and is **not** to be treated as a GPU/CPU optimization metric.

Therefore distinguish clearly:

```text
PERFORMANCE
  kernel / solver time
  electronic-structure phase time
  steady SCF iteration time

CORRECTNESS / STABILITY
  whether the route can converge
  whether CPU/GPU give the same observables
  whether a common restart/potential remains physically consistent
```

Do not tune mixing merely to make one backend converge in fewer iterations.

Full convergence may be used as a correctness anchor where practical, but it is not the main performance headline.

---

# 2. Mandatory two-phase workflow — Luna must not monitor long runs

The user does not want Luna to spend tokens waiting for long benchmark jobs.

## PHASE A — prepare only

In the first session:

1. inspect CURRENT HEAD, the B1 harness, B1 report, and canonical results;
2. make the minimal benchmark-integrity/reporting changes required below;
3. create **one resume-safe bench-all driver** for the complete B1R matrix;
4. run only short unit/parser tests and tiny smoke checks;
5. print the exact **single command** the user should run;
6. STOP.

Do **not** launch the long benchmark matrix yourself.
Do **not** poll it.
Do **not** wait for it.

Preferred user command:

```bash
bash tests/benchmarks/run_scf_b1r_all.sh
```

or, if a Python driver is clearly cleaner:

```bash
python3 tests/benchmarks/scf_b1r.py --all
```

Prefer a wrapper so the user needs only one command.

## PHASE B — after the user returns and says the command is done

Then:

1. read the generated B1R evidence;
2. verify completeness;
3. do not rerun expensive rows unless something is genuinely missing/corrupt;
4. rebuild canonical JSON/CSV/Markdown;
5. regenerate affected plots/tables;
6. reconcile `RS_LMTO_ASA_SCF_B1_REPORT.md`;
7. write a short B1R closeout report;
8. state whether SCF-B1 can now close.

The bench-all script must archive enough evidence that the user does not need to manually transcribe timings.

---

# 3. Bench-all driver requirements

Create preferably:

```text
tests/benchmarks/run_scf_b1r_all.sh
```

Use robust shell behavior such as:

```bash
set -euo pipefail
```

or an equivalent Python implementation.

The driver must:

- locate the repo root automatically;
- use a dedicated result directory:
  `results/benchmarks/scf_b1r/`;
- print/log every command before execution;
- archive stdout/stderr per case;
- record timestamps and exit status;
- be resume-safe;
- skip already-completed valid cases unless `--force` is requested;
- abort before performance timing if build preflight fails;
- retain ordinary benchmark FAIL / UNSUPPORTED / SKIPPED cases without destroying the rest of the matrix;
- write a machine-readable manifest;
- print a compact final summary;
- finish with a clear sentinel such as:

```text
B1R_RUN_COMPLETE
```

Useful options if simple to implement:

```text
--resume
--force
--dry-run
```

Do not create an elaborate new framework merely for these options.

---

# 4. Step 1 — verify whether `-O0` contaminated the Release build

The current B1 provenance records a suspicious Release flag string equivalent to:

```text
-O3 ... -O3 ... -O0
```

Because compiler flag ordering matters, this could make the supposedly Release benchmark effectively unoptimized.

## Required check

Do not rely on CMake cache metadata alone.

Inspect actual compile commands using whatever is available:

```text
compile_commands.json
CMakeFiles/*/flags.make
ninja -t commands
make VERBOSE=1
saved build logs
```

Inspect representative performance-critical Fortran objects from at least:

```text
SCF driver
real-space recursion
Hamiltonian construction
reciprocal/k-space path
```

Determine the **effective optimization flags actually passed to gfortran**.

Write:

```text
results/benchmarks/scf_b1r/build_preflight.json
```

with at least:

```text
build_dir
compiler
build_type
representative compile commands
effective optimization conclusion
status = PASS/FAIL
reason
timing_reuse_from_B1 = allowed/forbidden
```

### If `-O0` is not actually effective

Document why the B1 provenance string was misleading and fix provenance reporting if needed.

### If `-O0` is actually effective

This invalidates previous performance timings.

Do not benchmark with that build.

Create/configure a clean Release build, for example:

```text
build-b1r-release-cuda/
```

Verify the actual compile commands again before running benchmarks.

If the source CMake configuration itself incorrectly injects `-O0` into Release mode, make the smallest possible build-system fix.

Do not change numerical algorithms.

If the build was contaminated, the one-command B1R runner must rerun every B1R performance row with the corrected build.

---

# 5. Step 2 — resolve Fe2/Fe3 reciprocal `nmat`, `Nk`, and eigensolver timing

The current B1 report has reciprocal rows labelled `Fe 2x2x2` and `Fe 3x3x3` with missing `nmat`/`Nk` and nearly primitive-like CPU eigensolver timing.

Do **not** assume this is automatically wrong.

There are two possibilities:

## Case A — true reciprocal supercell

The reciprocal Hamiltonian really has a larger basis, e.g. approximately:

```text
2x2x2 -> nmat ~144
3x3x3 -> nmat ~486
```

Then a primitive-like dense eigensolver time would indeed be suspicious.

## Case B — the label denotes an RS/physical cluster replication while the reciprocal basis remains primitive

Then the dense k-space problem may legitimately remain around:

```text
nmat ~18
```

and nearly constant CPU eigensolver time is reasonable.

In that case the problem is metadata/report semantics, not performance.

## Required action

Trace the actual production path and determine which case applies.

For each relevant reciprocal row record runtime-authoritative:

```text
physical fixture / replication label
reciprocal basis replication
nmat
nominal k mesh
Nk_unique
solver strategy
vectors yes/no
```

Do not infer `nmat` from a fixture name.
Do not force 144 or 486.

If useful, distinguish explicitly:

```text
physical_replication
reciprocal_basis_replication
```

### Timer validation

For at least one primitive case and every genuinely enlarged reciprocal-basis case, verify that:

```text
CPU P_eigensolver
GPU T_solver
```

cover comparable solver work.

If the reciprocal basis remains primitive, retain the constant solver time and explain it.

If the basis is genuinely enlarged and the CPU timing boundary is wrong, repair only the instrumentation/reporting and rerun affected rows.

**Do not treat constant-time single-site real-space recursion as a k-space timing bug.**

---

# 6. Step 3 — complete the Fe3 RS CPU OMP sweep

The B1 block-recursion data contain Fe3 CPU OMP=1 and CUDA mixed, but not the fair CPU OMP sweep.

Run the Fe3 RS block workload with:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

with controlled BLAS/OpenMP settings and no oversubscription.

For this expensive case it is acceptable to use:

```text
1-2 warmups
3 measured repetitions
```

provided the reduced protocol is recorded.

Report separately the best CPU configuration for:

```text
RS solver kernel
RS electronic-structure phase
steady SCF iteration
```

The best OMP count may differ by stage.

The current RS CUDA route remains:

```text
mixed precision
```

so CPU FP64 vs CUDA mixed must remain a **production-route ratio**, not an equal-precision speedup.

The purpose of this step is to determine whether the apparent size trend:

```text
small: CPU preferred
Fe2: near parity
Fe3: GPU beneficial
```

survives comparison against the best-threaded CPU.

---

# 7. Step 4 — establish one RS Fe correctness/reference anchor

Do not turn this into a mixing-performance campaign.

The aim is simply to obtain a trustworthy common state for RS CPU/GPU correctness and for the RS-vs-k-space comparison.

## Preferred option A — existing validated `pot.par` / restart

Search the existing validated Fe examples/tests/results for a scientifically trustworthy converged or near-converged:

```text
pot.par
restart state
or equivalent potential state
```

Use it if available.

## Option B — generate one CPU FP64 reference state

If no suitable state exists, run one CPU FP64 calculation long enough to create a reference potential.

This is a correctness/reference preparation run.

Do not headline:

```text
number of cycles
full wall to convergence
```

as accelerator metrics.

Do not tune mixing to favor a backend.

## Common-state RS CPU/GPU check

Starting from exactly the same reference state, run a small fixed number of production iterations with:

```text
RS CPU FP64
RS CUDA mixed
```

Compare per-iteration/final observables such as:

```text
Fermi level
total charge
site charge
magnetic moment
stable energy quantities
production density/mixing residual
```

If CUDA also converges normally, retain that as useful correctness evidence.
If not, report it honestly.

The performance quantity remains:

```text
steady iteration time
```

not iteration count.

---

# 8. Step 5 — rerun current-build Fe 3^3 / 4^3 / 5^3 reciprocal eigensystem rows

The B1 report imported useful ACC-P4 FP64 evidence from an older commit.

Reproduce the essential large reciprocal rows on the current validated B1R build.

Target the real frozen-potential reciprocal workloads that actually produce approximately:

```text
nmat ~486
nmat ~1152
nmat ~2250
```

**only if those are the true runtime reciprocal matrix sizes.**

Runtime values are authoritative.

For each supported row run:

```text
CPU FP64
GPU FP64
eigenvectors = yes
```

Use the production solver, normally `fp64_zheevd`, unless CURRENT HEAD says otherwise.

Timing protocol:

```text
2 warmups + 3 measured repetitions
```

is sufficient for these expensive rows.

Record:

```text
frozen/current commit
nmat
Nk
vectors
CPU solver time
GPU solver time
CPU reciprocal/frozen-potential phase
GPU reciprocal/frozen-potential phase
GPU memory
correctness
```

If Fe5 is unsafe/unavailable, emit:

```text
SKIPPED
```

with exact reason.

After this step, old-commit ACC-P4 rows may remain as historical context but must not be the current B1 headline evidence.

---

# 9. Step 6 — one FP64 CPU RS-vs-k-space Fe comparison using a common potential

The primary formulation comparison should avoid GPU precision complications.

Compare:

```text
RS CPU FP64
vs
reciprocal CPU FP64
```

for bcc Fe.

## 9.1 Use the same potential state

Prefer the exact same validated:

```text
pot.par
```

or equivalent restart/potential state as input to both routes.

This removes mixing-history differences from the electronic-structure comparison.

## 9.2 Prefer fixed-potential evaluation

If the production code supports a true frozen/fixed-potential electronic-structure evaluation, use it.

Otherwise use the cleanest production equivalent, such as one iteration from the exact same restart and compare the observables at a common pre-/post-update boundary.

Do not invent a non-production code path.

## 9.3 Compare physical observables

At minimum compare:

```text
Fermi energy
integrated electron count / total charge
magnetic moment
stable directly comparable band/electronic energy quantity
```

Strongly prefer a DOS comparison as the main spectral check.

### DOS comparison requirements

Before comparing, put both DOS results on the same:

```text
energy window
energy grid
broadening/resolution convention
spin normalization
per-atom/per-cell normalization
```

Then report, where meaningful:

```text
relative L2 DOS difference
maximum absolute DOS difference
integrated DOS/electron-count difference
```

Do not directly compare differently broadened DOS curves.

## 9.4 Optional `pot.par` output comparison

If both routes perform one identical update from the common input potential, compare the resulting `pot.par`/potential-parameter fields as a **secondary diagnostic**.

Do not require byte identity.

Compare scientifically meaningful fields/norms.

## 9.5 Accuracy decision

Reuse existing validation tolerances where available.

Do not invent overly tight tolerances simply to force PASS/FAIL.

If no established tolerance exists for an observable, report the measured difference explicitly.

The cross-route result may be:

```text
PASS
FAIL
INCONCLUSIVE
```

Only `PASS` permits a formulation timing comparison.

## 9.6 Timing comparison

If accuracy passes, compare:

```text
RS CPU FP64 electronic-structure phase
reciprocal CPU FP64 electronic-structure phase
```

and, if the iteration boundary is scientifically equivalent:

```text
RS CPU FP64 steady SCF iteration
reciprocal CPU FP64 steady SCF iteration
```

Do not compare iteration counts.
Do not generalize one Fe result into a universal RS-vs-k-space claim.

---

# 10. Step 7 — repair report tables and status semantics

Update the report generator / B1 report so reciprocal SCF iteration rows visibly include:

```text
solver strategy
numeric mode
nmat
Nk_unique
status
failure/ineligibility reason
```

This is mandatory because the current human-readable table contains multiple CUDA rows with different solver strategies but does not identify the strategy.

For RS rows retain/show:

```text
solver
RS backend
nrec or M where applicable
numeric mode
OMP
status
failure/ineligibility reason
equal_precision_eligible
```

For all performance tables:

- do not show a speedup for an ineligible pair;
- show null/`-` plus reason;
- distinguish explicitly:
  `PASS`, `NOT_CONVERGED`, `NOT_APPLICABLE`, `UNSUPPORTED`, `FAIL`, `SKIPPED`;
- preserve the reason in canonical JSON/CSV.

Add a compact build-integrity summary near the methodology section:

```text
effective optimization
frozen/current commit
tracked dirty state
timing reuse from old B1 allowed?
```

---

# 11. Exact B1R bench-all matrix

The single command should run only what is needed for reconciliation.

## PRE — build preflight

Verify actual optimization. Abort performance timing if invalid.

## RS1 — Fe3 block recursion

```text
CPU OMP 1/2/4/8
CUDA mixed production route
```

Reuse a prior CUDA row only if same commit/build/physics and build preflight proves timing reuse valid; otherwise rerun.

## RS2 — Fe common-reference correctness

```text
common pot.par/restart
CPU FP64
CUDA mixed
small fixed iteration count
```

Optional full convergence only as correctness evidence.

## K1 — reciprocal Fe2/Fe3 metadata/timer validation

Run the minimum set required to determine whether these are:

```text
true enlarged reciprocal matrices
or
primitive reciprocal matrices attached to larger physical/RS fixture labels
```

and validate `P_eigensolver` / `T_solver` semantics.

## K2 — current large reciprocal rows

```text
Fe 3^3
Fe 4^3
Fe 5^3
CPU FP64
GPU FP64
vectors=yes
```

as supported.

## X1 — RS vs k-space Fe fixed-potential comparison

```text
same pot.par/restart
RS CPU FP64
reciprocal CPU FP64
common observables
DOS comparison where available
```

No broad cross-route sweep is required.

---

# 12. Runtime-conscious behavior

The whole point of the bench-all driver is that the user can launch it directly and return after it finishes.

The script should print progress such as:

```text
[B1R 1/12] build preflight
[B1R 2/12] Fe3 RS CPU OMP=1
[B1R 3/12] Fe3 RS CPU OMP=2
...
```

The script itself handles sequencing, logging, resume, and failures.

Luna must not monitor the long run.

---

# 13. Canonical outputs

Write under:

```text
results/benchmarks/scf_b1r/
```

with at least:

```text
build_preflight.json
manifest.json
campaign.json
campaign.csv
campaign.md
raw/
correctness/
iteration_history/
cross_route/
plots/
```

Reuse the B1 schema wherever possible.
Do not create an incompatible benchmark data model.

---

# 14. Plots to regenerate

Only regenerate plots affected by B1R.

Mandatory:

```text
Fe block CPU OMP scaling including Fe3
RS kernel/phase/iteration production ratio vs size
current-build large reciprocal FP64 solver speedup vs nmat
```

If the Fe cross-route accuracy test passes:

```text
RS vs k-space DOS overlay
```

and optionally a DOS-difference plot.

Do not plot a cross-route speed comparison if accuracy is `INCONCLUSIVE`.

---

# 15. Questions the final B1R report must answer

## Build

1. Was the B1 executable actually optimized?
2. If not, which measurements were rerun?

## Reciprocal metadata/timers

3. Do Fe2/Fe3 labels represent true larger reciprocal matrices or larger physical/RS fixtures with a primitive reciprocal basis?
4. What are actual runtime `nmat` and `Nk_unique`?
5. Are CPU/GPU eigensolver timer boundaries valid?

## RS

6. What is the best-threaded Fe3 CPU block-recursion time?
7. What is the corrected Fe3 mixed-GPU production ratio?
8. Does the size trend CPU-preferred -> parity -> GPU-beneficial survive the fair CPU sweep?
9. From a common potential/restart, do RS CPU and CUDA preserve the important physical observables per iteration?

## Large reciprocal

10. Are the historical large-matrix FP64 GPU gains reproduced on the current build?
11. Where is the current measured crossover?

## RS vs k-space

12. Using the same Fe potential and matched spectral resolution, do RS and reciprocal electronic structures agree?
13. What does the DOS comparison show?
14. If accuracy passes, which FP64 CPU formulation is faster for this measured case?
15. If accuracy does not pass, state `INCONCLUSIVE`.

## Reporting

16. Can a reader identify solver strategy and failure/ineligibility reason for every important row?

---

# 16. Acceptance checklist

## Build integrity

- [ ] actual compile commands inspected
- [ ] effective optimization established
- [ ] `-O0` contamination ruled out or corrected
- [ ] timing reuse policy recorded

## Reciprocal semantics

- [ ] Fe2 actual reciprocal basis size known
- [ ] Fe3 actual reciprocal basis size known
- [ ] `nmat` populated
- [ ] `Nk_unique` populated where meaningful
- [ ] fixture replication distinguished from reciprocal-basis replication if needed
- [ ] solver timers validated
- [ ] no RS constant-time behavior incorrectly called a k-space bug

## RS Fe3 fairness

- [ ] OMP1
- [ ] OMP2
- [ ] OMP4
- [ ] OMP8
- [ ] best CPU kernel selected
- [ ] best CPU phase selected
- [ ] best CPU iteration selected
- [ ] CUDA route remains labelled mixed
- [ ] production ratio recomputed

## RS correctness anchor

- [ ] validated common Fe potential/restart found or generated
- [ ] CPU FP64 run from common state
- [ ] CUDA mixed run from common state
- [ ] charge checked
- [ ] moment checked
- [ ] Fermi checked
- [ ] stable energy quantity checked
- [ ] residual/density change checked
- [ ] iteration count not used as accelerator metric

## Large reciprocal rerun

- [ ] Fe3 true large reciprocal case
- [ ] Fe4 true large reciprocal case
- [ ] Fe5 true large reciprocal case or explicit skip
- [ ] CPU FP64
- [ ] GPU FP64
- [ ] eigenvectors=yes
- [ ] correctness
- [ ] current commit provenance

## Cross-route Fe

- [ ] same input `pot.par` / potential state
- [ ] RS CPU FP64
- [ ] reciprocal CPU FP64
- [ ] fixed-potential/common-boundary comparison
- [ ] common DOS grid/broadening
- [ ] Fermi difference
- [ ] charge difference
- [ ] moment difference
- [ ] energy difference where comparable
- [ ] DOS difference metrics
- [ ] PASS/FAIL/INCONCLUSIVE
- [ ] timing comparison only if PASS

## Report

- [ ] reciprocal strategy visible
- [ ] numeric mode visible
- [ ] nmat/Nk visible
- [ ] failure reason visible
- [ ] ineligible speedups suppressed
- [ ] build-integrity summary visible
- [ ] affected plots regenerated
- [ ] B1 report reconciled

---

# 17. PHASE A stop checklist

Before handing control to the user:

- [ ] B1 harness/report inspected
- [ ] build preflight implemented
- [ ] reciprocal metadata/timer validation implemented
- [ ] report schema/table fix implemented
- [ ] cross-route Fe comparison driver implemented
- [ ] bench-all script created
- [ ] resume behavior checked
- [ ] dry-run checked
- [ ] only short tests/smokes run by Luna
- [ ] long campaign NOT launched by Luna
- [ ] exact single command printed for user

Then **STOP** and wait for the user to return after the command finishes.

---

# 18. PHASE B checklist

After the user says the bench-all run is complete:

- [ ] load manifest
- [ ] require build preflight PASS
- [ ] account for every expected case
- [ ] avoid unnecessary expensive reruns
- [ ] rebuild canonical campaign outputs
- [ ] select fair Fe3 CPU baselines
- [ ] resolve reciprocal basis semantics in report
- [ ] summarize current large-matrix FP64 evidence
- [ ] evaluate RS common-state correctness
- [ ] evaluate Fe RS-vs-k-space accuracy
- [ ] evaluate DOS comparison
- [ ] regenerate affected plots
- [ ] repair final tables
- [ ] update final B1 report
- [ ] write B1R closeout
- [ ] state close / not-close decision

---

# 19. B1R closeout report

Write:

```text
docs/dev/RS_LMTO_ASA_SCF_B1R_REPORT.md
```

Recommended sections:

1. Executive result
2. Build-integrity check
3. Reciprocal Fe2/Fe3 basis/timer resolution
4. Fe3 RS CPU OMP reconciliation
5. RS common-potential correctness
6. Current-commit large reciprocal FP64 rerun
7. Fixed-potential RS vs k-space Fe comparison
8. Report/schema corrections
9. Final workload conclusions
10. Closeout decision

---

# 20. Commit messages

Phase A harness/preparation commit:

```text
Prepare SCF B1 benchmark reconciliation
```

Phase B final reconciliation commit:

```text
Reconcile final SCF benchmark campaign
```

If only one final commit is desired, use the latter.

Do not commit machine-local build directories or transient raw artifacts unless existing repository policy explicitly requires them.

---

# 21. Final mindset

B1R is not:

```text
make SCF converge in fewer cycles
make CUDA win
optimize the kernels again
```

It is:

```text
make the final SCF performance evidence trustworthy
```

The main performance quantity is:

```text
time per production SCF iteration
```

The cross-route scientific comparison is anchored by:

```text
same Fe potential / pot.par
+ matched numerical resolution
+ matching DOS / charge / moment / Fermi observables
```

before comparing RS and reciprocal timing.

Operationally:

```text
Luna prepares one bench-all command
user runs it without Luna monitoring
user returns when complete
Luna gathers the stored evidence and reports
```
