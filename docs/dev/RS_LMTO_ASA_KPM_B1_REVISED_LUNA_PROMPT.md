# KPM-B1 — Frozen systematic real-material CPU/GPU benchmark campaign

You are working on **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task follows the completed KPM optimization and benchmark-harness campaign:

```text
G0      DONE
G1      DONE
G1.1    DONE
G1.2    DONE
G1.3    DONE
G1.4    DONE
G2      DONE
G3      SUPERSEDED
B0C     DONE
B1      THIS TASK
```

B0C is assumed to have produced the canonical benchmark workflow with:

```text
strict physics pairing
strict fp64/fp32/mixed numeric-mode classification
CPU OMP 1/2/4/8 fairness
deterministic random_vec pairing
correctness attached to rows
profile closure enforcement
best-same-precision CPU selection
JSON / CSV / Markdown output
raw-log provenance
Pt and magnetic bcc-Fe material support
headline-speedup eligibility enforcement
```

Do **not** redesign that harness in B1.

---

# 1. Objective

The purpose of B1 is to generate the final systematic performance evidence for
the KPM/Kubo-Bastin CPU/GPU implementation.

This is a **measurement and reporting campaign**, not an optimization task.

The final result must answer:

1. How much faster is GPU KPM transport than the best fair CPU implementation?
2. How does that speedup depend on:
   - real-space system size;
   - Chebyshev order M=`cond_ll`;
   - numerical precision;
   - conductivity type;
   - stochastic estimator / trace count?
3. Where is the CPU/GPU crossover?
4. Does the GPU advantage grow for the large systems where acceleration matters most?
5. What is the performance difference between:
   - FP64;
   - FP32;
   - mixed practical routes?
6. Does the optimized GPU moment kernel remain the dominant scalable advantage?
7. Does the G2 block-stochastic path help in practice?
8. What backend/precision should a production user select?

The final campaign must be credible to a skeptical technical audience.

Do not optimize for a favorable headline.

---

# 2. Freeze the implementation first

Before running the campaign:

```text
git status
git rev-parse HEAD
```

Require:

```text
working tree clean
```

Record the exact frozen commit.

Do not modify transport algorithms after timing begins.

If B1 reveals a genuine correctness or benchmark-contract defect:

1. stop the affected campaign;
2. document the defect;
3. make the smallest necessary fix;
4. create a new frozen commit;
5. rerun all affected paired rows.

Do not mix rows from different implementation commits in one headline table.

---

# 3. Primary hardware policy

Use one primary GPU for all paired headline measurements.

For the current campaign this should be one RTX A4000 unless the current test
host has changed.

Record:

```text
GPU model
device index
VRAM
driver
CUDA toolkit
compute capability
CUDA_VISIBLE_DEVICES
```

Do not mix A4000 device 0 and device 1 in the same paired performance series
unless explicitly running a secondary device-variance check.

Use:

```text
1 MPI rank
```

for the primary CPU/GPU campaign unless a separate MPI scaling study is
explicitly added later.

Multi-GPU and multi-node scaling are out of scope.

---

# 4. CPU fairness policy

For every CPU algorithm / numerical mode used in a headline comparison, run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

Use the B0C controlled threading policy.

Record:

```text
OMP_NUM_THREADS
BLAS_NUM_THREADS
OMP_PROC_BIND
OMP_PLACES
```

Prevent nested oversubscription.

For each valid comparison group define:

```text
best_cpu_same_precision
=
minimum median T_transport_total
```

over valid OMP=1/2/4/8 CPU rows.

Retain all CPU thread rows.

Do not compare GPU only against OMP=1.

---

# 5. Numeric modes

The final report must separate exactly:

## A. FP64 headline

```text
CPU FP64 moments + CPU FP64 reconstruction
vs
GPU FP64 moments + GPU FP64 reconstruction
```

Both rows must satisfy B0C:

```text
numeric_mode = fp64
```

## B. FP32 headline

```text
CPU FP32 moments + CPU FP32 reconstruction
vs
GPU FP32 moments + GPU FP32 reconstruction
```

Both rows must satisfy:

```text
numeric_mode = fp32
```

Canonical final host/output arrays may be widened to FP64 after the reduced
FP32 reconstruction result.

That does not change the numerical-mode classification.

## C. Mixed practical route

Keep separate:

```text
numeric_mode = mixed
```

For example:

```text
FP32 moments + FP64 reconstruction
```

may be practically useful.

Do not include MIXED rows in equal-FP32 or equal-FP64 headline tables.

---

# 6. Mandatory real-material systems

Synthetic matrices may appear only in existing correctness/unit tests.

All performance conclusions require real RS-LMTO Hamiltonians.

## 6.1 SOC Pt

This is the primary scaling material.

Use the validated real Pt SOC transport fixture.

Mandatory real-space sizes:

```text
r4
r6
r8
```

Report actual:

```text
N
Natom if available separately
Hamiltonian dimension
nnz
```

Do not report only replication labels.

Add one larger Pt system only if:
- memory permits safely;
- runtime remains practical;
- it provides meaningful extension of the crossover trend.

Do not make a larger system mandatory if it risks destabilizing the final campaign.

## 6.2 Magnetic bcc Fe

Run at least one real magnetic metallic bcc-Fe transport workload using the
validated material staging consolidated in B0C.

Prefer:

```text
spin transport
realistic M
same benchmark schema
```

Use at least one medium/large real-space size that is meaningful for transport.

If multiple Fe sizes are already clean and affordable, run a short size series.

Do not invent a toy magnetic system.

The purpose is to show that the performance conclusions are not unique to Pt.

---

# 7. Mandatory production anchor

The core production anchor remains:

```fortran
&control
calctype = 'B'
nsp = 2
cond_ll = 500
lld = 150
recur = 'chebyshev'
cond_type = 'spin'
cond_calctype = 'per_type'
/
```

with:

```text
M  = 500
NE = 2510
```

or the exact current fixture values if CURRENT HEAD records them slightly
differently.

Run this anchor for:

```text
Pt r4
Pt r6
Pt r8
```

in:

```text
CPU FP64 OMP 1/2/4/8
GPU FP64

CPU FP32 OMP 1/2/4/8
GPU FP32

mixed route(s) separately
```

This is the central final comparison.

---

# 8. Moment-order scaling

For at least:

```text
one medium Pt size
one large Pt size
```

run:

```text
M = 250
M = 500
M = 750
```

at fixed:

```text
material
real-space size
lld
energy range
estimator
cond_type
precision
```

Add:

```text
M = 1000
```

if memory and runtime are practical.

If M=1000 is skipped, record:

```text
SKIPPED
reason
```

Do not silently omit it.

This series must isolate the scaling with Chebyshev order.

Report:

```text
T_moments
T_gamma
T_reconstruction
T_transport
whole_wall
```

and speedup versus M.

---

# 9. Real-space size scaling

At fixed:

```text
M = 500
lld = 150
spin
per_type
```

run:

```text
Pt r4
Pt r6
Pt r8
```

for all headline numerical modes.

This is the primary plot/table for convincing collaborators.

Produce:

```text
system size / N
vs
CPU transport
GPU transport
S_moments
S_transport
S_whole
```

The final interpretation must state whether GPU advantage grows with system size.

Do not fit a crossover model unless the data support it.

If interpolation is used, label it as approximate.

---

# 10. Conductivity-type coverage

Mandatory:

```text
spin
```

At Pt r4 and preferably one larger Pt size, also run:

```text
charge
orbital
```

using the same physical state where supported.

The purpose is to verify that the performance behavior is generic across the
shared numerical transport machinery.

Do not create different physical fixtures merely to fill these rows.

Report any mode-specific difference honestly.

---

# 11. Estimator coverage

## 11.1 per_type

This is the primary production workflow.

Run the full mandatory Pt campaign with:

```text
cond_calctype = per_type
```

## 11.2 random_vec

Use the existing G2/B0C random-vector route.

Run a systematic but bounded stochastic series.

Recommended:

```text
Ntrace = 1, 2, 4, 8, 16
```

where scientifically meaningful and affordable.

Use identical:

```text
seed
random-vector distribution
normalization
Ntrace
physics
numeric_mode
```

for CPU/GPU pairs.

Report:

```text
total moment time
time per trace
traces/s
transport total
```

## 11.3 G2 block width

Do not assume block stochastic processing is beneficial.

Carry the G2 evidence honestly.

For representative GPU random_vec rows evaluate the already-supported block
widths, for example:

```text
B = 1, 2, 4, 8
```

and B=16 only if memory preflight allows.

Do not rerun an exhaustive block-width matrix if G2 already established the
trend convincingly.

The final report must include negative or neutral block-width scaling if that
is what the evidence shows.

---

# 12. Correctness is mandatory for headline rows

Every headline CPU/GPU pair must satisfy B0C:

```text
profile_status = PASS
correctness_status = PASS
headline_speedup_eligible = true
```

The correctness comparison must use the same physical workload.

For random_vec:

```text
same stochastic realization
```

is mandatory.

At minimum record:

```text
conductivity spectrum max abs
conductivity spectrum relative L2
integrated/tensor difference
moment comparison/evidence where available
tolerance set
```

Do not contaminate optimized GPU timing with diagnostic full-moment downloads.

Use separate correctness companion runs where needed.

Do not publish a speedup for a numerically failed row.

---

# 13. Production-output policy

Headline timing must use:

```text
normal production output
```

not benchmark-no-write mode.

The no-write diagnostic may be included in a secondary attribution section
only.

Do not use no-write timings as application performance claims.

Record output I/O time separately.

---

# 14. Timing protocol

For each benchmark configuration:

## Warmup

Use:

```text
2 warmups
```

where practical.

## Measured repetitions

Default:

```text
5 measured repetitions
```

For genuinely expensive large CPU rows:

```text
minimum 3 measured repetitions
```

is acceptable.

Document every reduced-repetition exception.

Report:

```text
median
min
max
MAD or IQR
```

Use medians for headline speedups.

Do not select the fastest individual sample.

---

# 15. Required timing quantities

For every summarized row retain:

```text
P_operator_setup
P_trace_setup
P_moments_total
P_gamma_basis_setup
P_gamma_generation
P_reconstruction_total
P_result_unpack
P_energy_integration
P_tensor_postprocess
P_output_prepare
P_output_io
P_misc
T_transport_total
whole_wall
```

and available nested GPU details such as:

```text
moment H2D
moment GPU kernel
moment D2H
GPU mu pack
Gamma generation
cuBLAS reconstruction
reduced-result D2H
```

Profiler closure must PASS.

Do not mix inclusive child timers into top-level totals.

---

# 16. Speedup definitions

Use exactly:

```text
S_moments =
    best_same_precision_CPU(P_moments_total)
    /
    GPU(P_moments_total)

S_transport =
    best_same_precision_CPU(T_transport_total)
    /
    GPU(T_transport_total)

S_whole =
    best_same_precision_CPU(whole_wall)
    /
    GPU(whole_wall)
```

For MIXED practical rows use a separate label:

```text
S_practical_mixed
```

Do not call `S_moments` the application speedup.

Do not call `S_whole` the kernel speedup.

---

# 17. Required campaign matrix

The final campaign should contain at least the following.

## Matrix A — Pt production size scaling

```text
material = Pt SOC
M = 500
lld = 150
cond_type = spin
cond_calctype = per_type

sizes:
    r4
    r6
    r8

routes:
    CPU FP64 OMP=1,2,4,8
    GPU FP64
    CPU FP32 OMP=1,2,4,8
    GPU FP32
    selected mixed route(s)
```

## Matrix B — M scaling

At one medium and one large Pt size:

```text
M = 250, 500, 750
M = 1000 if practical
```

Use FP64 and FP32 headline routes.

## Matrix C — conductivity modes

At Pt r4:

```text
charge
spin
orbital
```

At least FP64 and FP32 GPU plus the valid CPU comparators.

## Matrix D — stochastic estimator

At one representative Pt size:

```text
random_vec
Ntrace = 1,2,4,8,16
```

with valid CPU/GPU pairing.

Use selected block widths based on G2 evidence.

## Matrix E — magnetic material

At least one bcc-Fe workload:

```text
real magnetic Hamiltonian
realistic M
spin transport
per_type
```

with:

```text
CPU FP64 OMP=1,2,4,8
GPU FP64
CPU FP32 OMP=1,2,4,8
GPU FP32
```

where scientifically supported.

If a specific precision route is unsupported, report:

```text
UNSUPPORTED
```

rather than silently substituting another mode.

---

# 18. Required headline tables

## Table 1 — equal-precision FP64

Columns:

```text
material
size
N
nnz
M
estimator
Ntrace
best CPU OMP
CPU moment
GPU moment
S_moments
CPU transport
GPU transport
S_transport
CPU wall
GPU wall
S_whole
correctness
```

Only:

```text
numeric_mode = fp64
headline_speedup_eligible = true
```

rows.

## Table 2 — equal-precision FP32

Same schema.

Only:

```text
numeric_mode = fp32
headline_speedup_eligible = true
```

rows.

## Table 3 — mixed practical routes

Separate table.

Do not merge with equal-precision tables.

## Table 4 — CPU thread scaling

For representative r4/r6/r8 rows show:

```text
OMP=1
OMP=2
OMP=4
OMP=8
```

for FP64 and FP32.

## Table 5 — stochastic scaling

Show:

```text
Ntrace
block width
time per trace
traces/s
transport time
```

for CPU/GPU.

---

# 19. Required plots

Generate at least:

## Plot 1

```text
S_transport vs real-space system size
```

Separate:

```text
FP64
FP32
```

## Plot 2

```text
S_moments vs system size
```

## Plot 3

```text
S_transport vs M
```

for medium and large Pt systems.

## Plot 4

```text
CPU transport vs OMP threads
```

for representative systems.

## Plot 5

```text
GPU/CPU stochastic traces per second vs Ntrace
```

if random_vec coverage is sufficient.

## Plot 6

Stage-fraction plot or equivalent data presentation:

```text
moments
Gamma
reconstruction
postprocess
I/O
other
```

for representative CPU FP64 / GPU FP64 / CPU FP32 / GPU FP32 rows.

The purpose is to show why total speedup differs from moment-kernel speedup.

---

# 20. Required interpretation questions

The final B1 report must explicitly answer:

1. What is the FP64 CPU/GPU crossover?
2. What is the FP32 CPU/GPU crossover?
3. How does GPU advantage change from r4 to r8?
4. How does speedup change with M?
5. How much faster is the GPU moment stage than the full transport workflow?
6. What stages limit FP64 GPU speedup?
7. What stages limit FP32 GPU speedup?
8. Does CPU OpenMP scaling materially shift the crossover?
9. Is FP32 scientifically acceptable for the validated transport workflows?
10. Is FP64 GPU worthwhile for large systems?
11. Are charge/spin/orbital performance trends similar?
12. Does random_vec scale differently from per_type?
13. Does block stochastic processing help?
14. Are Pt and magnetic Fe performance conclusions qualitatively consistent?
15. What backend/precision should users select for small/medium/large systems?
16. Is any further KPM GPU optimization justified after this campaign?

Do not force a universal GPU recommendation.

---

# 21. No cherry-picking

Publish unfavorable rows.

If CPU wins for small systems, OMP=8 beats GPU in a regime, FP64 GPU crossover
is weak, FP32 fails a correctness tolerance, M=1000 cannot fit, block
stochastic B=8 is slower than B=1, or Fe behaves differently from Pt, state
that explicitly.

The purpose is a trustworthy performance map.

---

# 22. Output package

Use B0C's canonical output pipeline.

Produce:

```text
campaign.json
campaign.csv
campaign.md
raw logs
correctness evidence
plots
```

Use one stable results directory, preferably:

```text
results/benchmarks/kpm_b1/
```

or the repository's established equivalent.

Every final table/plot must be reproducible from the canonical dataset.

Do not manually edit numerical values into Markdown after generation.

---

# 23. Final report

Write:

```text
docs/dev/RS_LMTO_ASA_KPM_B1_REPORT.md
```

Recommended structure:

1. Executive summary
2. Benchmark methodology
3. Workloads
4. FP64 results
5. FP32 results
6. Mixed practical route
7. CPU thread scaling
8. Stochastic estimator / G2 results
9. Stage-level analysis
10. Magnetic Fe cross-check
11. Production recommendation
12. Remaining limitations
13. Final campaign conclusion

The production recommendation should be an evidence-based table such as:

```text
workload regime | preferred backend | preferred precision | caveat
```

Do not enable automatic backend switching in B1.

---

# 24. Production recommendation constraints

Do not change the application's backend defaults automatically in B1.

Hardware-specific thresholds from one RTX A4000 are not portable universal
constants.

Recommendations are evidence, not hard-coded policy.

---

# 25. Stop condition

B1 is complete when all mandatory campaign matrices are run or explicitly
skipped with reason; all headline pairs pass profile and correctness contracts;
FP64, FP32, and mixed tables exist; CPU OMP scaling is documented; Pt size
scaling, M scaling, random_vec/G2 evidence, and magnetic Fe evidence are
documented; plots are generated; backend recommendations are written; and
remaining limitations are explicit.

At that point:

```text
close the KPM GPU performance campaign.
```

Do not invent KPM-B2 merely because another small optimization might exist.

---

# 26. Checklist

## Final execution record

The corrected campaign was frozen at commit `e23fab86f10dda9ddc1d5dd9452f9586c2fc428a` after the initial run exposed and stopped on an OMP-wrapper propagation defect. The smallest harness fix was committed before timing, and all final rows below were rerun from that corrected freeze. The canonical package is `results/benchmarks/kpm_b1_final/`; the generated report is `docs/dev/RS_LMTO_ASA_KPM_B1_REPORT.md`.

The timed build is `build-b1-frozen-cuda/bin/rslmto.x`, Release, gfortran 13.3, oneMKL, serial/one-rank, CUDA plugin enabled. The primary device is CUDA device 0, NVIDIA RTX A4000, 16376 MiB, compute capability 8.6, CUDA 13.3. GPU timing used `OMP_NUM_THREADS=1`; CPU anchor sweeps used OMP 1/2/4/8 with BLAS threads 1 and `OMP_PROC_BIND=close`, `OMP_PLACES=cores`.

The final dataset retains failed or non-headline rows. In particular, strict FP64 correctness excluded Pt charge, Pt orbital, Fe spin, and Pt r6 M=750 GPU rows from headline speedups; the generated report records these failures. Pt r6 M=1000, Pt r8 M=750/M=1000, and random_vec Ntrace=16 are explicit bounded-campaign skips with reasons. The initial pre-fix artifacts remain outside the canonical final dataset.

## Freeze
- [x] freeze state recorded: tracked tree was clean at corrected freeze; pre-existing untracked artifacts preserved and disclosed
- [x] frozen commit recorded
- [x] build recorded
- [x] primary GPU fixed
- [x] CPU/hardware metadata recorded

## Pt production anchor
- [x] r4 FP64 CPU OMP 1/2/4/8
- [x] r4 GPU FP64
- [x] r4 FP32 CPU OMP 1/2/4/8
- [x] r4 GPU FP32
- [x] r6 same FP64 matrix complete
- [x] r6 same FP32 matrix complete
- [x] r8 same FP64 matrix complete
- [x] r8 same FP32 matrix complete
- [x] mixed route separately recorded (CPU route; mixed GPU route unavailable under current CUDA precision coupling)

## M scaling
- [x] medium Pt M=250/500/750
- [x] medium Pt M=1000 or explicit skip
- [x] large Pt M=250/500 measured; M=750 explicitly skipped with reason
- [x] large Pt M=1000 or explicit skip
- [x] FP64 M-scaling complete for affordable points; failed M=750 and skips retained
- [x] FP32 M-scaling complete for affordable points; large M=750/1000 bounded explicitly

## Conductivity modes
- [x] Pt charge
- [x] Pt spin
- [x] Pt orbital
- [x] equal-precision comparisons valid where correctness PASS; failed rows are ineligible and retained
- [x] mode-specific correctness assessed (FP32 charge/orbital PASS; strict FP64 charge/orbital failures disclosed)

## random_vec / G2
- [x] Ntrace=1/2/4/8
- [x] Ntrace=16 or explicit skip
- [x] identical CPU/GPU seeds
- [x] block width evidence included
- [x] traces/s reported
- [x] time/trace reported
- [x] neutral/negative G2 results retained

## Magnetic Fe
- [x] real bcc-Fe fixture used
- [x] realistic M used
- [x] FP64 CPU OMP sweep
- [x] GPU FP64
- [x] FP32 CPU OMP sweep
- [x] GPU FP32
- [x] correctness assessed; FP32 headline PASS and strict FP64 failure retained as ineligible
- [x] Fe/Pt comparison interpreted

## Correctness
- [x] every headline FP64 pair PASS
- [x] every headline FP32 pair PASS
- [x] random_vec same-realization checks PASS
- [x] failed/ineligible rows retained
- [x] no diagnostic full-moment D2H contaminates timing

## Timing
- [x] 2 warmups used where practical
- [x] 5 measurements used by default
- [x] expensive-row exceptions documented (r6/r8 and secondary M/G2 rows use 3 measurements)
- [x] medians used
- [x] spread reported
- [x] profile closure PASS
- [x] transport and whole-wall kept separate

## Tables
- [x] equal-FP64 table
- [x] equal-FP32 table
- [x] mixed table
- [x] CPU OMP table
- [x] stochastic table

## Plots
- [x] transport speedup vs size
- [x] moment speedup vs size
- [x] speedup vs M
- [x] CPU thread scaling
- [x] stochastic throughput
- [x] stage-fraction plot/data

## Interpretation
- [x] FP64 crossover stated
- [x] FP32 crossover stated
- [x] size scaling interpreted
- [x] M scaling interpreted
- [x] OpenMP effect interpreted
- [x] FP32 numerical acceptability stated
- [x] FP64 usefulness stated
- [x] random_vec behavior stated
- [x] block stochastic result stated
- [x] Fe generalization stated
- [x] production recommendation written
- [x] remaining limitations written
- [x] further KPM GPU work recommendation stated

## Outputs
- [x] campaign JSON
- [x] campaign CSV
- [x] campaign Markdown
- [x] raw logs archived
- [x] correctness evidence archived
- [x] plots generated
- [x] final B1 report written
- [x] master tracker updated
- [x] KPM GPU performance campaign closed if criteria met

---

# 27. Commit message

Use:

```text
Publish final KPM CPU GPU benchmark campaign
```

---

# 28. Final mindset

B1 is not:

```text
find another optimization
```

It is:

```text
freeze the code and produce the final fair evidence
```

The value of B1 is credibility:

```text
real systems
realistic M
best CPU OMP comparison
same precision
same stochastic realization
correctness attached
negative results retained
stage-level explanation
one frozen commit
```

At completion, collaborators should be able to inspect the dataset and
understand exactly when KPM GPU acceleration is useful, by how much, and under
what numerical assumptions.
