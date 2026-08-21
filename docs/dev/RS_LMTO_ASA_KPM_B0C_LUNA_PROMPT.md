# KPM-B0C — Close and consolidate the fair KPM CPU/GPU benchmark harness

You are working on **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task is the benchmark-harness closure step after the completed KPM GPU
optimization campaign.

The implementation sequence has effectively become:

```text
KPM-G0      baseline / ownership profiling                 DONE
KPM-G1      precision routes / persistent CUDA workspaces  DONE
KPM-G1.1    CPU BLAS Gamma*mu reconstruction               DONE
KPM-G1.2    exclusive profiler + CPU fairness campaign     DONE
KPM-G1.3    resident GPU moments + tiled Gamma + cuBLAS    DONE
KPM-G1.4    fixed-overhead cleanup                         DONE
KPM-G2      stochastic block-vector study                  DONE
KPM-G3      old planned reconstruction task                SUPERSEDED
                                                               |
                                                               v
KPM-B0C     THIS TASK
                                                               |
                                                               v
KPM-B1      frozen final real-material campaign
```

Do **not** reopen G0-G2.

Do **not** optimize the transport kernels in B0C.

The purpose of B0C is to turn the benchmark infrastructure already built
during G1.2/G2 into a strict, publication-grade harness that makes unfair
CPU/GPU comparisons difficult or impossible.

---

# 0. Current state that must be preserved

Inspect CURRENT HEAD before editing.

At minimum read:

```text
tests/benchmarks/benchmark_harness.py
tests/benchmarks/kpm_g12_transport.py
tests/benchmarks/kpm_g0_transport.py
tests/validation/val09_kubo_bastin_transport.py
source/kpm_profile.f90
source/conductivity.f90
source/recursion_transport.f90
docs/dev/RS_LMTO_ASA_KPM_TRANSPORT_GPU_FOLLOWUP.md
all G1.2/G1.3/G1.4/G2 reports currently present in the tree
```

The present infrastructure already provides substantial functionality.
Preserve and consolidate it rather than rewriting it.

In particular, CURRENT HEAD already has or has recently established:

```text
real SOC Pt r4/r6/r8 workload staging
M=500 / lld=150 / NE=2510 production-scale defaults
charge/spin/orbital selection
per_type and random_vec estimators
fixed stochastic seed for CPU/GPU random_vec pairing
GPU stochastic block-width controls
CPU cheb backends legacy / fast / fast_dp
GPU fp32 / fp64 controls
OMP thread sweep
BLAS thread controls
OMP_PROC_BIND / OMP_PLACES controls
warmups and repetitions
median / min / max / MAD / IQR
exclusive KPM phase parsing
profile-closure validation
nested child<=parent validation
moment time per trace
traces/second
best-CPU speedup machinery
JSON campaign output
memory-limit skip records for stochastic GPU blocks
benchmark-only no-write mode
```

The older G0 driver additionally knows how to stage:

```text
real SOC fcc Pt
magnetic bcc Fe
```

Do not lose that material coverage when consolidating the final harness.

---

# 1. B0C objective

At completion there must be one canonical KPM transport benchmark workflow
that can generate the B1 dataset with:

1. strict same-physics pairing;
2. strict numerical-precision labeling;
3. fair CPU thread selection;
4. deterministic stochastic pairing;
5. correctness evidence attached to performance rows;
6. full environment/build provenance;
7. JSON + CSV + Markdown outputs;
8. archived raw logs;
9. Pt and magnetic-Fe real-material support;
10. speedup claims emitted only when the pairing/correctness contract passes.

The final B0C harness must make it difficult to accidentally publish:

```text
CPU FP64 vs GPU FP32
different stochastic seeds
different M/lld/energy grids
different current operators
failed correctness rows
```

as an "equal-precision GPU speedup".

---

# 2. Non-negotiable rules

## 2.1 Do not redesign the physics

Do not change:

```text
Hamiltonian construction
current operators
SOC treatment
charge/spin/orbital conventions
Chebyshev scaling
kernel
Gamma definition
Kubo-Bastin formula
energy integration
per_type estimator
random_vec estimator
normalization
scientific output meaning
```

## 2.2 Do not optimize kernels

B0C is benchmark infrastructure.

The one permitted numerical addition is a **minimal CPU FP32 reconstruction
route** if CURRENT HEAD still lacks one and it is required for a truly
FP32-vs-FP32 benchmark.

Do not otherwise change performance kernels.

## 2.3 Preserve the optimized GPU production path

Do not regress:

```text
resident GPU moments
zero full-moment D2H in optimized execution
device-tiled Gamma
cuBLAS FP64/FP32 reconstruction
reduced-result D2H
G1.4 fixed-overhead cleanup
```

## 2.4 Do not weaken correctness tolerances

Reuse established validation tolerances from current G1/G1.3/VAL-09 evidence.

If no suitable tolerance exists for a new paired comparison, derive and
document one from the established numerical precision/error contract before
using it.

Do not choose a tolerance merely because a row would otherwise fail.

## 2.5 Unfavorable results are first-class evidence

Do not hide:

```text
CPU-faster rows
FP32 failures
memory-limited block widths
poor OpenMP scaling
negative G2 block-vector results
```

---

# 3. B0C-A — Audit the current harness and produce a gap table

Before implementation, produce a short table:

```text
B0 requirement | already present | partial | missing | action
```

Cover all original B0 requirements:

```text
benchmark dimensions
OMP 1/2/4/8
BLAS control
GPU metadata
warmups/repetitions
same-precision pairing
same stochastic states
stage timing
correctness per row
CPU-best reference
CSV
JSON
Markdown
raw logs
git/build/hardware metadata
```

Do not implement a second mechanism where the current one is already adequate.

The expected major gaps are likely to include:

```text
strict FP32 CPU reconstruction parity
correctness_status still being a placeholder
pairing key too permissive
CSV output
Markdown summary
raw-log archival
dirty-tree / fuller hardware metadata
material generalization of kpm_g12_transport.py
```

Confirm against CURRENT HEAD rather than assuming this list is exact.

---

# 4. B0C-B — Establish an explicit numerical-mode taxonomy

Every benchmark row must carry independent metadata for:

```text
moment_backend
moment_precision
reconstruction_backend
reconstruction_precision
canonical_output_precision
```

From these derive one machine-readable field:

```text
numeric_mode
```

with exactly:

```text
fp64
fp32
mixed
```

Suggested rule:

```text
fp64:
    moment_precision == fp64
    AND reconstruction_precision == fp64

fp32:
    moment_precision == fp32
    AND reconstruction_precision == fp32

mixed:
    all other supported combinations
```

Canonical host/output arrays may remain FP64 after widening.

Do not call:

```text
FP32 moments + FP64 reconstruction
```

an FP32 end-to-end numerical route.

It is:

```text
numeric_mode = mixed
```

This classification must appear in:

```text
JSON
CSV
Markdown
pairing logic
headline-table logic
```

---

# 5. B0C-C — Add the minimal CPU FP32 reconstruction route if still missing

## Why

The current optimized CPU `fast` route has historically meant approximately:

```text
FP32 Chebyshev/stochastic moments
+
FP64 host reconstruction
```

while the optimized GPU FP32 path after G1.3 is:

```text
FP32 moments
+
FP32 cuBLAS reconstruction
+
widen reduced C result to canonical FP64 host output
```

Those routes are useful practical comparisons, but they are not strict
FP32-vs-FP32 end-to-end numerical comparisons.

B1 requires a defensible:

```text
CPU FP32 vs GPU FP32
```

headline.

## Required implementation

First inspect CURRENT HEAD.

If a CPU FP32 reconstruction already exists and is correctly exposed, use it.

Otherwise implement the smallest possible CPU analogue of the validated G1.1
flattened reconstruction:

```text
G(NE,M*M) * U(M*M,nb) -> C(NE,nb)
```

using the project's CPU BLAS:

```text
CGEMM
```

for the FP32 route.

Requirements:

```text
FP32 Gamma/reconstruction arithmetic
FP32 packed moments
FP32 C result
widen only reduced C to canonical FP64 host integration/output
same q = n + (m-1)*M convention as G1.1/G1.3
same factor
same per_type/random_vec semantics
no conjugation change
no scalar fallback hidden under "FP32"
```

Keep the existing FP64 ZGEMM path unchanged.

Do not introduce a global precision framework.

Use a narrow transport/reconstruction strategy consistent with current style.

## Validation

Compare:

```text
CPU FP64 reference
CPU FP32
GPU FP32
GPU FP64
```

on focused deterministic contraction tests and real Pt outputs.

The CPU FP32 path must not become the default production route merely because
it exists for benchmarking.

---

# 6. B0C-D — Define the canonical comparison key

The current speedup grouping is not strict enough if it groups only by size,
conductivity type and numerical precision.

Create one canonical `comparison_key` / pairing fingerprint.

It must include all physics/numerical fields that must be identical for a
fair CPU/GPU comparison.

At minimum include:

```text
material / fixture identity
fixture revision or stable fixture identifier
replication / supercell
N
nnz
nsp
SOC state
cond_type
current-operator specification / va / vb / spin-orbital selector as relevant
cond_calctype
Ntrace / projector count
random seed or deterministic projector identifier
M = cond_ll
lld
NE / channels
rc
energy_min
energy_max
fermi
kernel and kernel parameters
Chebyshev scaling contract
numeric_mode
```

Do not include performance-strategy choices that are legitimately allowed to
differ, such as:

```text
OMP thread count
GPU stochastic block width
CPU backend implementation name
GPU energy tile width
```

Those belong in row metadata, not the physics pairing key.

For `random_vec`, pairing must require:

```text
same seed
same Ntrace
same vector distribution/normalization contract
```

For `per_type`, pairing must require the same deterministic type/projector
semantics.

Generate a stable hash/string representation of the comparison key in each row.

---

# 7. B0C-E — Enforce pairing, do not merely annotate it

Create a strict pairing validator.

A CPU/GPU pair is eligible for an equal-precision speedup only if:

```text
comparison_key identical
numeric_mode identical and in {fp32, fp64}
correctness_status == PASS
profile_status == PASS
normal production output mode used
```

For random-vector performance:

```text
same stochastic seed/realization required
```

If any condition fails:

```text
headline_speedup_eligible = false
```

and store a clear reason, for example:

```text
precision_mismatch
physics_mismatch
seed_mismatch
correctness_missing
correctness_failed
profile_failed
benchmark_no_write
```

Do not silently fall back to a nearby CPU row.

Do not compute an "equal precision" speedup for MIXED rows.

Mixed rows may have their own explicitly labeled practical comparison table.

---

# 8. B0C-F — Best CPU selection must happen only inside a valid pair group

For each valid physics + precision comparison group:

run/retain CPU:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

with controlled BLAS threading.

Then select:

```text
best_cpu_same_precision
=
minimum median T_transport_total
```

from CPU rows that:

```text
share the comparison key
share the numeric mode
PASS correctness
PASS profile validation
```

Retain all CPU thread rows in the output.

The selected row must record:

```text
moment backend
reconstruction backend
OMP threads
BLAS threads
transport median
wall median
```

Headline GPU speedup must use this selected CPU reference.

Do not compare GPU against OMP=1 merely because it is convenient.

---

# 9. B0C-G — Attach correctness evidence to performance rows

The current placeholder:

```text
correctness_status =
"attach validated moment/conductivity comparison before speedup claim"
```

must disappear from final B0C rows.

## Required correctness object

Every performance row must contain or reference a structured object such as:

```text
correctness:
    status: PASS / FAIL / NOT_APPLICABLE
    reference_row_id:
    tolerance_set:
    moment:
        available:
        max_abs:
        rel_l2:
    conductivity_spectrum:
        max_abs:
        rel_l2:
    integrated_or_tensor:
        max_abs:
        rel:
    validation_evidence_id:
```

Do not require full-moment download during every timed optimized GPU run.

That would invalidate the production timing path.

Instead:

### Performance run

Use the normal optimized path:

```text
full moment D2H = 0
```

where applicable.

### Correctness companion

For each unique backend/precision/workload pair, run or reuse a deterministic
validation companion that is outside the timed performance samples.

It may enable diagnostic full-moment download where necessary.

At minimum compare the actual production conductivity output generated from
the same physical workload.

For moment-level correctness, attach:

```text
direct moment comparison where available
OR
validated low-level/diagnostic evidence ID
```

Do not fabricate a moment error when the optimized route intentionally does not
download the full tensor.

## Same stochastic realization

For `random_vec`, correctness companion CPU/GPU runs must consume exactly the
same random vectors/seed.

The harness must fail the pairing if this cannot be established.

---

# 10. B0C-H — Use production output for numerical comparison

Reuse the existing VAL-09 philosophy:

```text
compare the production *_cond.out outputs;
do not reimplement the Kubo-Bastin estimator in Python.
```

The benchmark correctness layer should read and compare outputs, not compute a
parallel "reference conductivity" formula.

For each paired run compare whichever production files are relevant:

```text
cond_total.out
orbital-resolved outputs
per-type outputs
other established charge/spin/orbital files
```

Preserve current output semantics.

For charge/spin/orbital use the actual files produced by that route.

Record:

```text
max absolute difference
relative L2 difference
selected integrated/final-value difference
```

using established tolerances.

---

# 11. B0C-I — Raw-log and output provenance

Every measured sample must have a stable row/sample ID.

Archive raw execution evidence under the results tree, for example:

```text
results/benchmarks/kpm_b0c/
    campaign.json
    campaign.csv
    campaign.md
    raw/
        <row_id>/
            warmup_01.log
            warmup_02.log
            sample_01.log
            sample_02.log
            ...
    correctness/
        <pair_id>/
            cpu/
            gpu/
            comparison.json
```

Exact directory names may follow current repository conventions.

Do not rely on a single `testrun.log` that gets overwritten between samples.

Raw logs must be referenced from the JSON row/sample metadata.

For correctness companion runs, retain the compared output files or stable
checksums plus sufficient provenance to reproduce them.

Do not unnecessarily archive giant diagnostic moment tensors if a small
comparison summary and reproducible test is sufficient.

---

# 12. B0C-J — Produce JSON, CSV, and Markdown from one canonical dataset

JSON remains the canonical full-fidelity format.

Generate:

```text
campaign.json
campaign.csv
campaign.md
```

from the same in-memory row data.

Do not independently calculate results in three separate code paths.

## JSON

Contains:

```text
full environment
full physics metadata
all samples
statistics
correctness objects
pairing keys
speedup eligibility
speedups
skipped rows
raw-log references
```

## CSV

One row per summarized benchmark configuration.

Include at minimum:

```text
row_id
material
size
N
nnz
M
NE
lld
cond_type
cond_calctype
Ntrace
seed
moment_backend
moment_precision
reconstruction_backend
reconstruction_precision
numeric_mode
OMP threads
BLAS threads
GPU strategy
block width
profile status
correctness status
transport median
whole-wall median
moment median
Gamma median
reconstruction median
output median
S_moments
S_transport
S_whole
headline_speedup_eligible
best_cpu_row_id
```

Keep individual repetition samples in JSON rather than exploding the main CSV,
unless a separate sample CSV is useful.

## Markdown

Generate a human-readable summary that includes:

```text
environment
physics/workload
CPU thread winner by group
equal-FP64 pairs
equal-FP32 pairs
mixed practical rows
failed/ineligible rows
memory-limited rows
correctness summary
definitions of S_moments / S_transport / S_whole
```

Do not cherry-pick favorable rows.

---

# 13. B0C-K — Expand environment and build provenance

The existing `capture_environment()` already records useful core metadata.

Extend it conservatively.

Required where available:

```text
git commit
git dirty state
compiler path/name/version
CMake build type
relevant optimization flags
BLAS/LAPACK vendor
OMP threads
BLAS threads
OMP_PROC_BIND
OMP_PLACES
MPI ranks
CPU model
physical/logical CPU count
RAM
OS/kernel
CUDA toolkit
CUDA driver
GPU model
selected GPU/device index
GPU VRAM
GPU compute capability
CUDA_VISIBLE_DEVICES
```

Useful optional fields if available without fragile dependencies:

```text
CPU governor
NUMA topology/binding
GPU persistence mode
```

If a field is unavailable:

```text
null
```

plus optionally a short reason.

Do not make the benchmark fail because a laptop or CI host cannot provide CPU
governor/NUMA/GPU metadata.

## Selected GPU

Do not record only a comma-separated list of all visible GPUs.

Record which device the benchmark actually selected.

For the final B1 campaign, one GPU must be the primary device for paired rows.

---

# 14. B0C-L — Record the exact physical fixture contract

The current G1.2 driver is Pt-centric.

Consolidate material staging so the final harness can use at least:

```text
fccPt_SOC
bccFe_magnetic
```

without maintaining separate incompatible benchmark schemas.

Reuse the proven G0 material staging rather than inventing a new Fe fixture.

A material specification should provide:

```text
material name
base fixture path
SOC state
nsp
default cond_type
va
vb
minimum/current rc
replication semantics
any material-specific required patch
```

Keep this simple.

Do not create a general materials database/framework.

The final B1 harness must be able to run:

```text
--material pt
--material fe
```

or equivalent explicit material selections.

B0C needs only smoke/contract validation for Fe.

The full magnetic performance campaign belongs to B1.

---

# 15. B0C-M — Preserve estimator semantics and G2 evidence

For:

```text
cond_calctype = per_type
```

retain the optimized production route.

For:

```text
cond_calctype = random_vec
```

retain:

```text
fixed seed pairing
Ntrace metadata
block width metadata
memory-limit skip behavior
moment time per trace
traces/second
```

Do not change the G2 algorithm.

Do not conceal the negative/neutral block-width scaling result.

For headline CPU/GPU stochastic comparisons, CPU and GPU must share:

```text
same seed
same Ntrace
same physics
same precision mode
```

GPU block width may differ as an implementation strategy.

---

# 16. B0C-N — Clarify stage timing versus whole-process timing

The current benchmark driver launches the application as a process for each
measured sample.

That is acceptable for whole-application evidence.

Do not rewrite RS-LMTO-ASA into an in-process persistent benchmark solely for
B0C.

Instead clearly distinguish:

```text
T_transport_total
    internal exclusive transport phase

whole_wall
    complete process invocation measured by the driver
```

Warmup runs should remain explicit.

GPU context/startup cost that lies outside `T_transport_total` is therefore
visible in `whole_wall`.

If CURRENT HEAD already has a clean persistent benchmark entry point, it may
be retained as an additional metric, but do not create a large new application
API in this task.

B1 must report both transport and whole-wall speedups with their definitions.

---

# 17. B0C-O — Benchmark-no-write rows are diagnostic only

The G1.4 no-write mode is useful for attributing I/O.

Keep it.

However:

```text
output_mode = benchmark_no_write
```

must automatically make:

```text
headline_speedup_eligible = false
```

for final production-performance tables.

No-write rows may appear in diagnostic sections only.

Normal B1 headline rows must use production output behavior.

---

# 18. B0C-P — Add hard contract/unit tests for the harness

Do not require real multi-minute benchmark execution in normal unit tests.

Add focused tests for:

## Pairing

```text
same physics + same precision + PASS -> eligible

fp64 vs fp32 -> rejected

mixed vs fp32 -> rejected

different M -> rejected

different lld -> rejected

different NE -> rejected

different material -> rejected

different seed for random_vec -> rejected

correctness missing -> rejected

correctness FAIL -> rejected

profile FAIL -> rejected

benchmark_no_write -> rejected for headline
```

## Best CPU

Verify selection among:

```text
OMP=1,2,4,8
```

uses the minimum median from valid same-precision rows only.

## Output

Verify:

```text
JSON schema
CSV required columns
Markdown required sections
```

## Environment

Verify missing optional metadata is represented safely rather than crashing.

## Material staging

Smoke-test Pt and Fe fixture selection without requiring a full B1 campaign.

---

# 19. B0C-Q — Produce a dry-run B0C evidence package

After implementation, run a manageable real-material closure campaign.

It does NOT need to be the final B1 matrix.

Minimum:

```text
Pt r4
M=500
NE=2510
lld=150
spin
per_type
```

Run:

```text
CPU FP64 OMP 1/2/4/8
GPU FP64

CPU FP32 OMP 1/2/4/8
GPU FP32

current mixed CPU route
```

Use:

```text
2 warmups
>=3 measured repetitions
```

for this closure evidence if five runs are unnecessarily expensive.

Also run one small/controlled:

```text
random_vec
```

pair to prove seed matching and pairing logic.

Run one bcc-Fe smoke/short evidence row to prove the material-generalized
harness.

This is harness validation, not the final published campaign.

---

# 20. B0C-R — Speedup definitions

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

Do not mix these quantities.

Do not call `S_moments` the application speedup.

Do not call `S_whole` the KPM-kernel speedup.

For MIXED rows, use an explicitly different label such as:

```text
practical_mixed_speedup
```

rather than an equal-precision headline.

---

# 21. B0C-S — Reconcile the master tracker

Update:

```text
docs/dev/RS_LMTO_ASA_KPM_TRANSPORT_GPU_FOLLOWUP.md
```

so the task state matches the implementation.

Do not rewrite the historical scientific discussion.

Add/repair a compact status section stating:

```text
G0    DONE
G1    DONE
G1.1  DONE
G1.2  DONE
G1.3  DONE
G1.4  DONE
G2    DONE
G3    SUPERSEDED by G1.3/G1.4
B0C   DONE after this task
B1    NEXT
```

Tick stale B0-equivalent checkboxes only when B0C actually satisfies them.

Mark old G3 as superseded rather than pretending its original checklist was
performed verbatim.

State that B1 is the next/final substantive campaign.

---

# 22. B0C acceptance criteria

B0C is complete only when all of the following are true:

```text
1. one canonical KPM benchmark driver/workflow exists for B1;

2. Pt and bcc-Fe material selection are supported;

3. CPU FP64 / GPU FP64 exact numerical-mode pairing is available;

4. CPU FP32 / GPU FP32 exact numerical-mode pairing is available,
   OR a technically justified "unsupported" status is recorded;
   do not fake FP32 parity with a mixed CPU row;

5. MIXED routes are explicitly classified;

6. pairing uses the full physics/numerical comparison key;

7. random_vec pairing enforces identical seed/Ntrace;

8. correctness is structured and attached;
   no placeholder correctness string remains;

9. speedup eligibility is enforced programmatically;

10. best CPU selection uses valid OMP 1/2/4/8 rows at the same precision;

11. JSON is generated;

12. CSV is generated;

13. Markdown is generated;

14. raw logs are archived and referenced;

15. git dirty state and fuller hardware/build metadata are captured;

16. benchmark-no-write rows cannot become headline results;

17. profile closure remains mandatory;

18. focused harness contract tests pass;

19. Pt r4 real closure campaign succeeds;

20. Fe harness smoke test succeeds;

21. master follow-up document is reconciled;

22. B1 can be run without further benchmark-harness redesign.
```

---

# 23. Checklist

## Current-harness audit

- [x] CURRENT HEAD inspected
- [x] existing G1.2/G2 features preserved
- [x] B0 gap table written
- [x] no duplicate benchmark framework introduced

## Precision

- [x] moment precision recorded
- [x] reconstruction precision recorded
- [x] canonical output precision recorded
- [x] numeric_mode={fp64,fp32,mixed} implemented
- [x] current mixed CPU fast route labelled MIXED
- [x] CPU FP32 reconstruction route inspected
- [x] minimal CPU CGEMM route added if required
- [x] CPU FP32 correctness validated
- [x] GPU FP32 correctness validated
- [x] FP64 routes unchanged

## Pairing

- [x] canonical comparison key implemented
- [x] material included in key
- [x] fixture/physics identity included
- [x] nsp/SOC included
- [x] current-operator contract included
- [x] M/lld/NE included
- [x] energy window/Fermi included
- [x] kernel included
- [x] estimator included
- [x] Ntrace included
- [x] random seed included for random_vec
- [x] numeric mode included
- [x] strict eligibility validator implemented
- [x] mismatched precision rejected
- [x] mismatched seed rejected
- [x] mismatched physics rejected
- [x] no-write headline rejected
- [x] failed correctness headline rejected
- [x] failed profile headline rejected

## Correctness

- [x] placeholder correctness string removed
- [x] structured correctness object implemented
- [x] production conductivity outputs compared
- [x] max-abs metric recorded
- [x] relative-L2 metric recorded
- [x] integrated/tensor metric recorded
- [x] established tolerance set recorded
- [x] diagnostic moment evidence attached where available
- [x] optimized GPU performance path keeps full-moment D2H disabled
- [x] random_vec correctness uses identical realization

## CPU fairness

- [x] OMP=1 automated
- [x] OMP=2 automated
- [x] OMP=4 automated
- [x] OMP=8 automated
- [x] BLAS thread policy controlled
- [x] oversubscription prevented
- [x] best CPU selected only from valid same-precision rows
- [x] all CPU thread rows retained

## Materials

- [x] Pt staging retained
- [x] bcc-Fe staging consolidated from G0
- [x] material metadata recorded
- [x] Pt smoke/closure run passes
- [x] Fe smoke run passes
- [x] no toy/synthetic material substituted for B1 support

## Estimator

- [x] per_type retained
- [x] random_vec retained
- [x] fixed seed recorded
- [x] Ntrace recorded
- [x] GPU block width recorded
- [x] memory-limited skips retained
- [x] G2 negative results not hidden

## Environment

- [x] git commit recorded
- [x] dirty state recorded
- [x] compiler/version recorded
- [x] build type recorded
- [x] optimization flags recorded where available
- [x] BLAS/LAPACK recorded
- [x] OMP/BLAS settings recorded
- [x] CPU model recorded
- [x] CPU count recorded
- [x] RAM recorded
- [x] OS/kernel recorded
- [x] CUDA toolkit recorded
- [x] CUDA driver recorded
- [x] selected GPU recorded
- [x] GPU model recorded
- [x] GPU VRAM recorded
- [x] compute capability recorded
- [x] CUDA_VISIBLE_DEVICES recorded
- [x] optional unavailable fields handled safely

## Output/provenance

- [x] stable row IDs implemented
- [x] stable sample IDs implemented
- [x] raw warmup logs retained
- [x] raw measurement logs retained
- [x] correctness outputs/provenance retained
- [x] JSON generated
- [x] CSV generated
- [x] Markdown generated
- [x] all three derived from same canonical row data
- [x] speedup definitions included
- [x] ineligible rows preserved with reasons

## Tests and closure

- [x] precision mismatch contract test
- [x] seed mismatch contract test
- [x] physics mismatch contract test
- [x] correctness missing/fail contract tests
- [x] profile fail contract test
- [x] no-write ineligibility test
- [x] best-CPU-selection test
- [x] JSON schema test
- [x] CSV schema test
- [x] Markdown-section test
- [x] missing-metadata robustness test
- [x] Pt material-selection test
- [x] Fe material-selection test
- [x] real Pt B0C closure campaign run
- [x] random_vec pairing closure run
- [x] bcc-Fe smoke run
- [x] master tracker reconciled
- [x] B1 declared next

---

# 24. Required B0C closeout report

Write a concise development report, preferably:

```text
docs/dev/RS_LMTO_ASA_KPM_B0C_REPORT.md
```

It must state:

## What already existed

Describe what G1.2/G2 had already built.

## What B0C added

List only new closure/consolidation work.

## Precision contract

State exactly what:

```text
FP64
FP32
MIXED
```

mean.

## Pairing contract

List the comparison-key fields and rejection rules.

## Correctness contract

State how timed rows receive correctness evidence without contaminating
optimized timing with diagnostic full-moment transfers.

## Material support

State Pt and Fe fixtures.

## Environment/provenance

Show example captured metadata.

## Output package

Show JSON/CSV/Markdown/raw-log layout.

## Closure evidence

Give the small B0C real Pt campaign results only as harness evidence.

Do not turn B0C into the final performance publication.

## Remaining work

State explicitly:

```text
KPM-B1 is next.
No further benchmark-harness redesign is expected before B1 unless B0C finds
a correctness or pairing defect.
```

---

# 25. One-line commit message

Use:

```text
Close fair KPM CPU GPU benchmark harness
```

---

# 26. Final mindset for Luna

This task is not:

```text
"build another benchmark script"
```

It is:

```text
"turn the benchmark machinery we already have into one strict evidence
pipeline that cannot accidentally make an unfair speedup claim."
```

At completion, B1 should be mostly:

```text
freeze commit
choose campaign matrix
run
collect
plot
interpret
publish
```

not another round of harness engineering.
