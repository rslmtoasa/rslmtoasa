# RS-LMTO-ASA Phase III-A — Current GPU Performance Steering
## Post ACC-P0 / ACC-P1 / ACC-P2 / ACC-P3 update — active entry point: ACC-P4

**Target branch:** `fable_v3`
**Purpose:** keep the accelerator campaign focused on measured GPU performance after the original ACC-00...ACC-09 work and the subsequent ACC-P rescue campaign.

---

# 1. Executive steering decision

The accelerator effort is no longer in the state assumed by the original rescue document.

The current evidence has established three important facts:

1. **The original ACC-06 result was misleading as a global performance conclusion.**
   Fresh-process CUDA startup and tiny matrices made the GPU appear roughly 20--100x slower.

2. **True batching helps the very-small-matrix CUDA route, but does not make primitive Si/Fe competitive with optimized CPU LAPACK on the tested RTX A4000.**

3. **Large real LMTO Hamiltonians already have a genuine FP64 CUDA winning regime.**
   The real bcc-Fe supercell campaign shows the current dense `Zheevd` CUDA path crossing the CPU route between approximately `nmat=486` and `nmat=1152`, with a strong win by `nmat=2250`.

Therefore the current Phase III-A goal is no longer:

> "Can reciprocal CUDA work at all?"

It is now:

> **Determine the best precision/solver/transfer strategy as a function of real LMTO matrix size and workload, optimize the already-useful medium/large-matrix CUDA path, and then decide which downstream GPU work packages are justified.**

## Current KPM interlude

The established production Pt anchor has exposed a separate CPU bottleneck:
the host `Gamma_nm * mu_nm` reconstruction, not the CUDA stochastic-moment
kernel. KPM-G1.1 is therefore an active scoped interlude under ACC-13:
replace the scalar `O(NE*M^2*nb*Ntrace)` accumulation with a CPU FP64
`ZGEMM` baseline before any stochastic-vector batching or GPU reconstruction.

This interlude preserves the Kubo-Bastin algebra, current conventions,
per-type/random-vector semantics, normalization, and canonical outputs. It is
not a reason to reopen the reciprocal P4 evidence or to add a CUDA Gamma*mu
kernel. Production-scale timing and real-Pt observable validation are now
recorded in the KPM follow-up document; the CPU reconstruction is validated
for the production `M=500/lld=150` workload.

---

# 2. Current work-package status

Use this table as the primary steering reference.

| Work package | Status | Current interpretation | Action |
|---|---|---|---|
| ACC-00 | Completed | Benchmark infrastructure established | Archive / reuse |
| ACC-01 | Completed | Existing RS CUDA correctness scope established | Reuse in ACC-13 |
| ACC-02 | Completed | Reciprocal CUDA backend/lifecycle established | Keep |
| ACC-03 | Completed for correctness/infrastructure | FP64 reciprocal CUDA eigensolver works; original performance expectation was incomplete | Superseded by P0/P1/P3 evidence |
| ACC-04 | Completed | Arbitrary-k CUDA integration established | Keep |
| ACC-05 | Completed | Normal reciprocal mesh/SCF CUDA correctness established | Keep |
| ACC-06 | Completed as historical evidence, **not** final performance policy | Tiny/cold-start campaign showed no crossover | Superseded by P0/P1/P3 |
| ACC-07 | Completed | H(k) host-materialization optimization was not justified then | Revisit only if P2/P4 profiling says so |
| ACC-08 | Completed | MPI rank-to-device support established | Keep |
| ACC-09 | Completed | Validated standard-Hermitian operator variants use CUDA path | Keep |
| **ACC-P0** | **Completed** | Persistent real Si/Fe + Fe-supercell benchmark corpus established | Permanent benchmark foundation |
| **ACC-P1** | **Completed** | True batched Jacobi implemented for supported small matrices; large matrices remain dense Zheevd | Permanent solver evidence |
| **ACC-P1b** | **Evidence input** | FP32/mixed-precision study; consume its current result if present on HEAD | Do not restart; include result in P4, or record FP32 as not yet established |
| **ACC-P2** | **Completed** | Hot-path optimization has been performed | Use its post-optimization timings as P4 evidence |
| **ACC-P3** | **Completed** | Decisive real-material 2D crossover campaign is done | Reuse its corpus/results; only rerun targeted rows if P2 changed their timings materially |
| **ACC-P4** | **NEXT ACTIVE WP** | Convert current P1/P1b/P2/P3 evidence into production CPU/GPU/precision policy | Start here |
| ACC-10 | Pending, likely relevant | GPU Lehmann contractions | Proceed after P4, with route chosen from P4 evidence |
| ACC-11 | Conditional | Device-resident eigensystem handoff | Only if GPU eigensolver + GPU Lehmann and transfer cost justify it |
| ACC-12 | Conditional / low priority | GPU H(R)->H(k) assembly | Only if P2/P4 show H(k) assembly has become material |
| **ACC-13 / KPM-G1.1** | **CPU VALIDATION COMPLETE; TUNING REMAINS** | Establish the production CPU `Gamma*mu` reconstruction baseline with zero-copy Gamma storage and reusable diagonal-moment packing | Preserve the validated CPU path; defer KPM-G2/GPU reconstruction until separately justified |
| ACC-14 | Final gate | Accelerator support/performance matrix | Run last |

---

# 3. Current validated benchmark corpus

The performance campaign must continue using **real production LMTO workloads**, not hand-filled matrices labelled as materials.

The established real-material corpus includes:

- diamond Si primitive / sp;
- bcc Fe primitive / spd;
- explicit bcc-Fe `2x2x2`;
- explicit bcc-Fe `3x3x3`;
- explicit bcc-Fe `4x4x4`;
- explicit bcc-Fe `5x5x5`.

The Fe supercells use the production:

```text
explicit lattice
    -> structure constants
    -> H(R)
    -> reciprocal H(k)
    -> eigensolver
```

path.

The repeated-supercell Hamiltonians are physically checked by band folding before being trusted as benchmark inputs.

For the current collinear Fe/spd state the runtime dimensions are:

```text
bcc Fe 1^3 : 18
bcc Fe 2^3 : 144
bcc Fe 3^3 : 486
bcc Fe 4^3 : 1152
bcc Fe 5^3 : 2250
```

Si primitive is currently reported as `nmat=16` in the real-material ACC-P campaign.

These actual runtime dimensions, not hand-coded assumptions, are the benchmark dimensions of record.

---

# 4. Established FP64 performance picture

The ACC-P0/P1 evidence already demonstrates a real matrix-size crossover.

Representative matched-density **values-only** rows from ACC-P1 were:

| Fixture | nmat | CPU LAPACK | CUDA `zheevd_serial` | CUDA `zheevj_batched` |
|---|---:|---:|---:|---:|
| Si primitive | 16 | 2.44e-5 s | 5.88e-4 s | 3.70e-4 s |
| bcc Fe 1^3 | 18 | 3.49e-5 s | 6.44e-4 s | 3.21e-4 s |
| bcc Fe 2^3 | 144 | 1.83e-3 s | 6.82e-3 s | unsupported for current batched API |
| bcc Fe 3^3 | 486 | 1.93e-2 s | 2.60e-2 s | unsupported for current batched API |
| bcc Fe 4^3 | 1152 | 1.62e-1 s | 1.02e-1 s | unsupported for current batched API |
| bcc Fe 5^3 | 2250 | 1.30 s | 4.82e-1 s | unsupported for current batched API |

The current interpretation is:

```text
n ~ 16--18:
    CPU strongly preferred

n ~ 144:
    CPU preferred

n ~ 486:
    crossover / near-parity region

n ~ 1152:
    GPU preferred

n ~ 2250:
    GPU strongly preferred
```

This is already sufficient to reject two incorrect global conclusions:

```text
"CUDA reciprocal eigensolution is always slower"
```

and

```text
"CUDA should be the default for all reciprocal calculations"
```

Neither is supported.

---

# 5. What ACC-P1 established

ACC-P1 is considered **successfully completed**.

The current small-matrix CUDA path now contains a true same-size batched strategy rather than a host loop pretending to be a batch.

The resulting evidence is scientifically useful even though it is negative for primitive cells:

- true batching improves the primitive CUDA path materially;
- optimized CPU LAPACK remains much faster for primitive Si/Fe;
- the current documented cuSOLVER batched route is limited to the small-matrix regime in this implementation/API path;
- the large-matrix GPU wins come from the retained dense `Zheevd` route;
- unsupported larger batched cases are reported explicitly rather than silently falling back.

Therefore:

> **Do not continue spending major engineering effort trying to make `n=16--18` primitive eigensolution beat MKL on the current hardware.**

Primitive cells remain:
- correctness fixtures;
- overhead probes;
- possible future high-`Nk` batching tests;
- not the main GPU performance target.

---

# 6. Numerical correctness status entering P1b/P2

The FP64 CPU/GPU numerical agreement is excellent.

Representative maximum errors reported by ACC-P1 are:

```text
max |delta eigenvalue|        ~ 4.2e-15
max residual                  ~ 4.9e-14
max orthogonality error       ~ 4.5e-15
max degenerate-projector err  ~ 1.2e-13
```

CPU supercell band folding passed for all Fe supercells, with reported sorted-multiset errors of order `1e-13` in the focused P1 report.

Therefore the FP64 reciprocal CUDA route is not merely "fast enough to benchmark"; it is a strong scientific reference against which FP32 can now be evaluated.

---

# 7. ACC-P1b — precision evidence

## Status

ACC-P1b is no longer an entry-point task in this steering document.

The assistant starting from current HEAD should **read and use whatever ACC-P1b result is present**, but should not restart the task.

If the P1b evidence is present, ACC-P4 must incorporate its scoped conclusion, e.g.:

```text
FP32 acceptable for ordinary reciprocal SCF
FP32 acceptable only for selected workloads
FP32 scientifically inadequate
```

If P1b evidence is absent or incomplete, do **not** create a new gating loop before P4. ACC-P4 should simply record:

```text
FP32 policy: not yet established
```

and define the production policy from the established FP64 evidence.

The application-level scientific arrays remain FP64 unless a scoped FP32 mode has actually passed its own correctness evidence.

---

# 8. Optional P1c — only if already justified by existing P1b evidence

P1c is **not part of the default active sequence**.

Do not stop before P4 to invent a P1c.

Only formulate a narrow P1c later if the existing P1b report already demonstrates all of the following:

- a large FP32 performance gain;
- otherwise-good eigenvectors/subspaces and SCF behavior;
- one narrow, understood precision deficit;
- a cheap correction such as FP64 Rayleigh-quotient reevaluation or FP64 cleanup iterations is strongly motivated.

Otherwise leave P1c unscheduled.

---

# 9. ACC-P2 — completed

ACC-P2 has already been performed.

Do **not** ask the next assistant to repeat its optimization work.

The current assistant should read the P2 report and extract, at minimum:

- the post-P2 `vectors=yes` timings for real Fe `3^3`, `4^3`, and `5^3`;
- the final breakdown of `H(k)`, staging, H2D, solver, D2H eigenvectors, synchronization, and total time;
- which hot-path changes were retained;
- which attempted optimizations were rejected as non-beneficial;
- whether eigenvector transfer is now material enough to motivate ACC-11;
- whether CPU `H(k)` assembly is now material enough to motivate ACC-12;
- whether the `n~486` crossover moved;
- whether the `n~1152` and `n~2250` GPU advantages improved.

These are inputs to ACC-P4, not reasons to reopen ACC-P2.

# 10. ACC-P3 — completed evidence corpus

ACC-P3 is complete.

Its real-material crossover corpus is an **evidence source**, not a pending work package.

After ACC-P2, only rerun the specific established P3 rows needed to quantify the effect of P2 on current HEAD. Do not recreate the full P3 development task.

Use the same real fixtures and benchmark conventions so pre-P2 and post-P2 results remain comparable.

The current assistant should treat the latest valid P3/P2 combination as the performance dataset feeding ACC-P4.

---

# 11. ACC-P4 — NEXT ACTIVE WORK PACKAGE

ACC-P4 is now the correct entry point for a new assistant.

There is no additional "gate" before it.

## What the assistant should do first

Before changing production code, perform a short read-only evidence consolidation on current HEAD:

1. read the ACC-P1 solver/correctness report;
2. read the current ACC-P1b result if present;
3. read the completed ACC-P2 optimization report;
4. read the completed ACC-P3 crossover report/corpus;
5. identify the latest post-P2 benchmark rows for real Si and bcc-Fe `1^3...5^3`;
6. confirm which rows are `vectors=yes` and therefore relevant to SCF;
7. record the current CPU, FP64 CUDA, and—if established—FP32 CUDA winners.

This is preparation **inside ACC-P4**, not a separate work package.

## ACC-P4 objective

Freeze an honest production policy for reciprocal eigensolution based on the evidence already generated.

Classify by the dimensions that actually matter:

```text
nmat
actual Nk / tile
vectors yes/no
precision
solver strategy
cold vs persistent workload
```

Allowed outcomes include:

```text
CPU-preferred
FP64 CUDA-preferred
FP32 CUDA-preferred, scoped
no clear benefit
unsupported
```

Do not force one backend to win everywhere.

## Required policy decisions

ACC-P4 must decide:

- small/primitive matrix policy;
- medium-matrix policy around the measured crossover;
- large-matrix policy;
- FP32 status, if P1b evidence exists;
- whether LAPACK remains the global safe default;
- whether an optional conservative `auto` selection is justified;
- whether ACC-11 is justified by measured eigenpair transfer;
- whether ACC-12 is justified by measured H(k)-assembly cost;
- how ACC-10 should receive eigenpairs;
- how ACC-13 should be prioritized.

Do not silently execute CPU inside an explicitly selected CUDA mode.

Do not make CUDA globally default based on one RTX A4000 campaign.

## Performance wording

Use real-material, persistent, post-P2 results.

For SCF policy, use `vectors=yes` evidence.

Values-only results may define a separate spectral/bands policy but may not be presented as SCF speedup.

# 12. Revised relevance of original ACC-10...ACC-14

This is the main steering section for what happens after P4.

---

## ACC-10 — GPU Lehmann Green-function contractions

### Status
**Still relevant and likely high value.**

### Why

Lehmann contractions have much higher arithmetic intensity than primitive-cell eigensolution.

The current large-matrix eigensolver evidence makes the attractive path:

```text
CPU H(k)
    -> GPU eigensolve for large n
    -> GPU Lehmann contraction
    -> small canonical G output
```

If P4 concludes small/medium eigensystems remain CPU-preferred, ACC-10 is still valid:

```text
CPU eigensolver
    -> copy eigenpairs once
    -> GPU Lehmann contraction
```

Therefore ACC-10 should **not depend on GPU eigensolution winning everywhere**.

### Action
Proceed after P4.

Use P4 to select whether ACC-10 starts from:
- CPU eigenpairs;
- GPU eigenpairs;
- both, depending on regime.

---

## ACC-11 — device-resident eigensystem handoff

### Status
**Conditional.**

### Proceed only if

All of these are true:

1. ACC-10 GPU Lehmann is worthwhile;
2. a workload region uses GPU eigensolution;
3. P2 shows D2H eigenvectors + ACC-10 H2D eigenvectors are material;
4. reuse frequency/memory justify lifecycle complexity.

### Skip/defer if

- GPU eigensolution is CPU-preferred in the targeted Lehmann regime; or
- eigenpair transfer is a small fraction of solve+contraction.

Do not add residency merely because GPU eigensolution exists.

---

## ACC-12 — GPU H(R)->H(k) assembly

### Status
**Conditional and currently lower priority.**

P2 must report the new H(k)-assembly fraction after solver optimization.

Proceed only if CPU H(k) assembly has become a material bottleneck.

If it remains small:
- close ACC-12 with "not justified";
- keep the canonical CPU Fourier assembler.

Do not port Fourier assembly preemptively.

---

## ACC-13 — RS KPM/Kubo-Bastin GPU completeness and performance

### Status
**CPU validation complete; retain the measured baseline for future KPM scope decisions.**

### KPM-G1.1 — CPU `Gamma*mu` reconstruction baseline

The current scalar contraction is mathematically:

```text
I_l^t(E_i) = factor * sum_(n,m) Gamma_nm(E_i) * mu_nm_stochastic(l,l,n,m,t)
```

For `M=cond_ll`, `NE=channels_ldos+10`, and `K=M*M`, the implementation now
forms `C_t = factor * G(NE,K) * U_t(K,nb)` with:

```text
q        = n + (m-1)*M
G(i,q)   = gamma_nm(i,n,m)
U_t(q,l) = mu_nm_stochastic(l,l,n,m,t)
```

The existing Fortran storage order makes `G` a zero-copy rank-2 view of the
allocated Gamma tensor. Only the diagonal `mu` blocks are packed into a
reusable workspace; for `M=500`, `nb=18`, and complex FP64 this is about 72 MB.
`ZGEMM('N','N',...)` is intentional: the existing scalar expression has no
complex conjugation. The scalar `gamma_mu_reference` helper remains available
for focused validation, and no user-facing physics option was added.

The profiler now distinguishes `T_mu_pack` from `T_gamma_mu` and emits exact
`bytes_gamma`/`bytes_mu_pack` reconstruction workspace sizes. The focused
`UnitKpmTransport` test covers the encoded storage layout and scalar/BLAS
agreement for both `per_type` and `random_vec` accumulation semantics.

Production closure is recorded in
`docs/dev/RS_LMTO_ASA_KPM_TRANSPORT_GPU_FOLLOWUP.md`: real SOC-Pt `M=500/lld=150`
spin, charge, and orbital runs completed; valid larger Pt replications
`N=3888` and `N=9216` completed; the optimized r4 spin output agrees with the
pre-interlude scalar output to `1.12e-13` relative L2 or better; and the
direct `Gamma*mu` stage fell from `490.636 s` to `1.874 s`. The optimized
production path is therefore ready as the CPU reconstruction baseline.

The remaining work is tuning and scope selection, not correctness closure:
the measured host thread sweep is complete for `OPENBLAS_NUM_THREADS=1`, the
full `random_vec` production campaign is not yet required, and no GPU
`Gamma*mu` kernel should be added without a new measured justification.

This task is largely orthogonal to reciprocal eigensolver policy.

It should remain a major Phase III-A target because:
- RS large-system algorithms are natural GPU workloads;
- existing CUDA machinery already covers substantial recursion/moment functionality;
- Phase II established charge/spin/orbital Kubo-Bastin CPU scientific contracts.

If P4 concludes reciprocal GPU has only a narrow useful regime, ACC-13 becomes even more important.

---

## ACC-14 — accelerator release gate

### Status
**Still required as the final task.**

It must publish an honest support/performance map, not merely a feature list.

The final documentation should distinguish:

```text
supported
validated
performance-beneficial
performance-tested but CPU-preferred
experimental
unsupported
```

for both RS and reciprocal GPU paths.

---

# 13. Recommended sequence from the current point

The live sequence is now simple:

```text
CURRENT HEAD
    |
    v
ACC-P4  <-- START HERE
    |
    +-----------------------------+
    |                             |
    v                             v
ACC-10                         ACC-13
GPU Lehmann                    RS KPM/transport
    |
    +--> ACC-11 only if P2/P4 transfer evidence justifies residency

ACC-12 only if P2/P4 show H(k) assembly is material

    |
    v
ACC-14 final accelerator support/performance gate
```

ACC-P1b, ACC-P2, and ACC-P3 are **inputs to P4**, not preliminary gates that the next assistant should execute again.

If the current P1b evidence is incomplete, ACC-P4 records FP32 as unestablished and proceeds with the FP64 policy. It does not block the rest of Phase III-A.

# 14. Work packages that should NOT be restarted

Do not restart these merely because older steering documents list them earlier in the chain:

```text
ACC-00 ... ACC-09
ACC-P0
ACC-P1
ACC-P3
```

Their code/evidence should be reused.

Only reopen one if:
- P2 exposes a concrete correctness defect in it;
- a benchmark contract is shown to be invalid;
- P1b/P4 proves a specific implementation assumption wrong.

Do not rerun whole development campaigns for documentation symmetry.

---

# 15. Performance-first interpretation of current results

The current GPU story is already useful:

```text
primitive/small reciprocal matrices:
    CPU/MKL is the likely production choice

medium reciprocal matrices:
    performance boundary region

large reciprocal matrices:
    FP64 CUDA dense eigensolution already wins

FP32:
    potentially shifts the crossover lower and/or increases large-n speedup,
    but only if the scientific precision gate passes

RS recursion/KPM:
    separate GPU performance domain

Lehmann:
    promising higher-arithmetic-intensity reciprocal GPU consumer
```

The objective is not to force one backend to win everywhere.

The objective is to produce a **measured heterogeneous execution policy**.

---

# 16. Rules for all remaining tasks

These rules remain mandatory.

## 16.1 Real workloads

Performance claims require real LMTO Si/Fe/production Hamiltonians.

Synthetic matrices are supplementary microbenchmarks only.

## 16.2 No silent fallback

Explicit CUDA means CUDA.

Unsupported CUDA must report unsupported/failure.

Do not make a CPU result look like a CUDA result.

## 16.3 CPU scientific oracle

CPU FP64 remains the reference unless a scoped mixed-precision mode is explicitly validated.

## 16.4 Eigenvectors matter

For SCF performance, `vectors=yes` is the primary workload.

Do not advertise values-only speedups as SCF speedups.

## 16.5 Cold vs persistent

Keep startup and steady-state timings separate.

## 16.6 Same timing interval

Never compare lifetime component counters to one-call wall times.

## 16.7 Physics stays unchanged

Do not alter:
- structure constants;
- LMTO Hamiltonian conventions;
- reciprocal phase conventions;
- SCF tolerances;
- occupations;
- reference physics

to create a GPU speedup.

## 16.8 Negative results are valid

A task may correctly conclude:
- CPU remains faster;
- FP32 is not accurate enough;
- pinned memory is not worth it;
- H(k) GPU assembly is unnecessary;
- device residency is not worth the complexity.

Do not keep adding complexity to avoid a negative result.

## 16.9 Checkbox discipline

Tick a performance checkbox only when the measured requirement actually passed.

Do not tick:
- "optimized" because code changed;
- "GPU preferred" because CUDA ran;
- "validated" because no exception occurred.

## 16.10 One focused commit

Every remaining implementation work package should end with:
- completed checklist;
- exact files changed;
- tests run;
- real-hardware benchmark evidence;
- unsupported/negative results;
- one focused single-line commit.

---

# 17. START HERE — instruction for the next assistant

**Do not start at the top of this document.**

The completed history is:

```text
ACC-00...09  completed
ACC-P0       completed
ACC-P1       completed
ACC-P2       completed
ACC-P3       completed
```

ACC-P1b is an evidence input if its result is present on current HEAD.

## Your first active work package is ACC-P4.

Begin ACC-P4 with a short read-only consolidation of the existing P1/P1b/P2/P3 evidence, then make the reciprocal CPU/GPU/precision policy decisions described in Section 11.

Do not repeat:
- benchmark-harness construction;
- supercell construction;
- true-batching implementation;
- P2 hot-path optimization;
- P3 crossover campaign.

Only rerun established benchmark rows when needed to obtain current post-P2 numbers on the same corpus.

After ACC-P4:

```text
ACC-10 and ACC-13 are the main implementation tracks.
ACC-11 is conditional on transfer evidence.
ACC-12 is conditional on H(k)-assembly evidence.
ACC-14 closes the accelerator campaign.
```

The current steering principle is:

> **The benchmarking/solver-rescue phase is complete. Freeze the measured execution policy in ACC-P4, then spend engineering effort on the downstream workloads that can deliver additional real application speedup.**
