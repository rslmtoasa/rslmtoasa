# RS-LMTO-ASA KPM GPU Follow-up
## KPM-G1.2 — trustworthy benchmark/profiling semantics
## KPM-G1.3 — GPU reconstruction on the optimized BLAS track

**Target branch:** `fable_v3`

**Starting point:** current HEAD after KPM-G1 and KPM-G1.1.

**Important established evidence:**

A realistic production-like Pt spin transport workload has:

```text
N = 1152
M = cond_ll = 500
NE = 2510
lld = 150
cond_type = spin
cond_calctype = per_type
```

KPM-G1.1 replaced the old scalar Gamma*mu accumulation by a BLAS formulation.

Representative measured results supplied after G1.1:

```text
before scalar:
    T_cheb_moments ~522.057 s
    T_gamma         ~8.024 s
    T_gamma_mu    ~490.636 s
    transport     ~522.526 s
    wall          ~522.816 s

after BLAS, spin r4:
    T_cheb_moments  ~33.359 s
    T_gamma          ~8.298 s
    T_mu_pack        ~0.095 s
    T_gamma_mu       ~1.874 s
    transport       ~33.822 s
    wall            ~34.153 s
```

The BLAS reconstruction itself is a clear success:

```text
Gamma*mu:
    ~490.6 s -> ~1.87 s
```

However the reported timing categories overlap:

```text
T_cheb_moments + T_gamma + T_gamma_mu > T_transport_total
```

so the current timing labels cannot yet be interpreted as mutually exclusive
top-level stages.

KPM-G1.2 fixes this measurement problem and establishes fair CPU/GPU baselines.

KPM-G1.3 then makes the GPU reconstruction follow the **same optimized algebra**
as the G1.1 CPU track. It must not accelerate the obsolete scalar loops.

---

# COMMON RULES FOR G1.2 AND G1.3

```text
You are working on CURRENT HEAD of:

https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3

This is a performance/scientific-correctness task.

Do not optimize for a favorable headline.
Do not change the physical Kubo-Bastin formula.
Do not change the stochastic estimator.
Do not change the Chebyshev normalization/scaling.
Do not loosen numerical tolerances merely to accept a new backend.

============================================================
MANDATORY ENGINEERING RULES
============================================================

1. INSPECT CURRENT HEAD BEFORE EDITING

Read at least:
- source/kpm_profile.f90
- source/conductivity.f90
- source/recursion_transport.f90
- current CUDA stochastic-moment implementation
- KPM-G1 report/evidence
- KPM-G1.1 report/evidence
- current transport benchmark scripts

Do not assume the source still matches an older prompt.

2. PRESERVE PHYSICS

Do not change:
- H;
- current operators;
- charge/spin/orbital-current conventions;
- Lorentz/Jackson kernel choice;
- energy scaling;
- Gamma_nm definition;
- Kubo-Bastin prefactors;
- energy grid;
- Simpson/integration semantics;
- per_type/random_vec normalization.

3. SAME-PRECISION COMPARISONS

Keep separate:

    CPU FP64 vs GPU FP64
    CPU FP32 vs GPU FP32

Mixed precision is labelled MIXED and never presented as an equal-precision
speedup.

4. CPU THREAD FAIRNESS

For final CPU baselines run:

    OMP_NUM_THREADS = 1, 2, 4, 8

and control BLAS threading explicitly.

Do not accidentally run:
    OpenMP 8 threads * MKL 8 threads

unless nested execution is deliberately being studied.

5. REAL WORKLOADS

The mandatory anchor is real Pt:

    M=500
    lld=150
    spin
    per_type

Use r4, r6, r8 production real-space sizes where established.

Do not replace it with M=20/40.

6. NO SILENT FALLBACK

An explicit CUDA route must not silently execute the CPU reconstruction.

7. NEGATIVE RESULTS ARE VALID

If the optimized CPU reconstruction is already faster than a GPU
reconstruction in a regime, record that and retain CPU.

8. ONE FOCUSED COMMIT PER WP

Tick only completed checkboxes.
Record exact tests and hardware.
Record unfavorable results.
Use the supplied one-line commit message.
```

---

# KPM-G1.2 — Repair KPM profiling and establish fair post-G1.1 baselines

## Objective

Make the KPM performance data mathematically interpretable.

The current profiler contains useful stage timers, but the reported values mix
inclusive/nested scopes and exclusive CUDA component times. As a result,
numbers such as:

```text
T_cheb_moments
T_gamma
T_gamma_mu
T_transport_total
```

cannot currently be summed or used directly to infer the dominant exclusive
stage.

G1.2 must establish a timing model where:

1. top-level phase timers are disjoint;
2. nested detail timers are clearly labelled as nested;
3. backend and precision are explicit;
4. CPU/GPU comparisons use the same physical workload;
5. CPU threading is fairly swept;
6. the final benchmark report identifies the actual new bottleneck after G1.1.

Do NOT perform another major kernel optimization in this WP.

---

## G1.2-A — Audit the timing scopes

Before changing code, produce a short table for every existing profiler key.

For each key record:

```text
timer name
start location
stop location
CPU path / CUDA path / both
inclusive or exclusive
may overlap with which other timer
wall-clock source or CUDA-event source
accumulated over how many calls
```

At minimum inspect:

```text
T_operator
T_trace_setup
T_cheb_moments
T_H2D
T_D2H
T_gamma
T_mu_pack
T_gamma_mu
T_energy_integral
T_transport_total
```

If G1.1 introduced additional timer names, include them.

Explicitly explain why the old G1.1 table could show:

```text
T_cheb_moments + T_gamma + T_gamma_mu > T_transport_total
```

Do not merely rename fields to hide the overlap.

---

## G1.2-B — Define two timing namespaces

Use a simple distinction:

### 1. Exclusive top-level phases

These phases must not overlap and should approximately partition the profiled
transport workflow:

```text
P_operator
P_trace_setup
P_moments_total
P_gamma
P_reconstruction_total
P_energy_integration
P_output_io
P_other
```

The exact names may follow existing style, but the semantics must be explicit.

For one complete measured transport calculation require:

```text
P_operator
+ P_trace_setup
+ P_moments_total
+ P_gamma
+ P_reconstruction_total
+ P_energy_integration
+ P_output_io
+ P_other
~= T_transport_total
```

within timer resolution and explicitly documented setup outside the transport
scope.

A small tolerance such as 2-5% is acceptable for reporting/clock overhead,
but the harness must FLAG larger discrepancies.

Do not fake `P_other=0`.
If something is uninstrumented, measure the residual honestly.

### 2. Nested/detail timers

Retain useful implementation detail:

```text
D_moment_H2D
D_moment_GPU_kernel
D_moment_D2H

D_mu_pack
D_reconstruction_BLAS

D_gamma_basis
D_gamma_fill
```

where applicable.

These detail timers are children of a top-level phase.

They are NOT added again when calculating the total.

For example:

```text
P_moments_total

contains:
    D_moment_H2D
    D_moment_GPU_kernel
    D_moment_D2H
    host conversion / setup / other
```

and:

```text
P_reconstruction_total

contains:
    D_mu_pack
    D_reconstruction_BLAS
```

This avoids the present parent/child double counting.

Do not build a generic hierarchical profiler framework.
A small fixed transport-specific representation is sufficient.

---

## G1.2-C — Give CPU and CUDA moment timers identical semantic boundaries

The CPU and GPU `P_moments_total` measurements must cover the same conceptual
work:

```text
input projector/reference state ready
        ->
mu_nm_stochastic available in the representation required by the next
reconstruction stage
```

For the current CPU route that means the complete CPU moment calculation.

For the current GPU route decide explicitly whether the boundary ends with:

```text
A. mu moments downloaded to host
```

or:

```text
B. mu moments valid/resident on GPU
```

For G1.2, if the production CUDA path still downloads moments, use A.

G1.3 will introduce a resident reconstruction route.

Within the GPU parent timer report separately:

```text
H2D/setup
CUDA moment kernel
D2H
FP32->FP64 widening if applicable
```

Do not use the CUDA-event kernel timer as the CPU-vs-GPU "moment total" if the
CPU timer includes extra work.

---

## G1.2-D — Backend and precision metadata are mandatory

Every emitted profile row must contain machine-readable fields for:

```text
moment_backend:
    cpu_legacy
    cpu_fast
    cuda

moment_precision:
    fp64
    fp32
    mixed

reconstruction_backend:
    scalar_reference
    cpu_blas
    cuda_blas   [future G1.3]

reconstruction_precision:
    fp64
    fp32
    mixed

OMP_NUM_THREADS
BLAS_NUM_THREADS
```

Use actual detected/selected values where feasible.

At minimum the benchmark harness must record the environment values used.

The old ambiguous statement:

```text
T_cheb_moments = 33 s
```

is not sufficient.

The new report must say, for example:

```text
moment_backend=cpu_fast
moment_precision=fp32
P_moments_total=...
```

or:

```text
moment_backend=cuda
moment_precision=fp32
D_moment_GPU_kernel=...
P_moments_total=...
```

---

## G1.2-E — Separate output/file I/O from numerical integration

The conductivity path performs energy integration and writes several output
files.

Do not hide file I/O inside a numerical timer.

Measure separately:

```text
P_energy_integration
P_output_io
```

If integration and writing are currently interleaved and cannot be separated
without invasive refactoring, perform the smallest safe split.

Do not change output content/order/format.

This distinction matters because the final GPU route may make numerical stages
so fast that file output becomes visible.

---

## G1.2-F — Correctly profile Gamma generation and G1.1 reconstruction

For G1.1 record:

```text
P_gamma
P_reconstruction_total

D_mu_pack
D_reconstruction_BLAS
```

Require:

```text
D_mu_pack + D_reconstruction_BLAS <= P_reconstruction_total
```

up to measurement noise.

Record:

```text
Gamma bytes
mu packed bytes
integrand/result bytes
```

For the mandatory M=500 / NE=2510 case, report the actual allocated Gamma
memory rather than a hand-estimated value only.

Do not redesign Gamma storage in G1.2.

---

## G1.2-G — Add profiler consistency checks

The benchmark parser should calculate:

```text
exclusive_sum =
    P_operator
  + P_trace_setup
  + P_moments_total
  + P_gamma
  + P_reconstruction_total
  + P_energy_integration
  + P_output_io
  + P_other
```

and:

```text
profile_closure_error =
    abs(exclusive_sum - T_transport_total) / T_transport_total
```

Report it for every row.

If closure error exceeds the chosen tolerance:

```text
PROFILE_STATUS=FAIL
```

and do not use that row for bottleneck conclusions.

Also sanity-check children:

```text
sum(moment children) <= P_moments_total + tolerance
sum(reconstruction children) <= P_reconstruction_total + tolerance
```

A row with impossible stage accounting must not silently enter the final CSV.

---

## G1.2-H — Fair CPU thread sweep

For each realistic CPU path run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

For the BLAS reconstruction, control the actual BLAS library threads.

Recommended baseline policy:

### CPU legacy / OpenMP-heavy moment path

```text
OMP_NUM_THREADS = requested
BLAS threads = 1
```

unless current implementation requires otherwise.

### CPU fast path

Use the threading model actually intended by that routine and prevent nested
oversubscription.

### CPU BLAS reconstruction

Run a consistent BLAS-thread sweep.

If reconstruction is called outside an OpenMP region, it may use the requested
1/2/4/8 BLAS threads.

Document:

```text
OMP_NUM_THREADS
MKL_NUM_THREADS or OPENBLAS_NUM_THREADS
OMP_PROC_BIND
OMP_PLACES
```

Do not quote GPU speedup only versus OMP=1.

For each precision/workload define:

```text
best_cpu = minimum median transport time over 1/2/4/8
```

---

## G1.2-I — Precision-matched benchmark matrix

Run at least these routes where current code supports them:

```text
CPU FP64 moments + CPU FP64 BLAS reconstruction

CPU FP32 moments + CPU FP32 BLAS reconstruction

GPU FP64 moments + CPU FP64 BLAS reconstruction

GPU FP32 moments + CPU FP32 BLAS reconstruction
```

If the application canonical output is widened to FP64 after an FP32
reconstruction, label the numerical reconstruction itself FP32.

If CPU FP32 BLAS reconstruction does not yet exist, implement only the minimal
`CGEMM` analogue necessary for a fair benchmark. Do not redesign transport.

Mixed routes such as:

```text
GPU FP32 moments -> widen -> CPU ZGEMM FP64
```

may also be measured, but label them:

```text
MIXED
```

and keep them out of equal-precision headline speedups.

---

## G1.2-J — Mandatory workload corpus

Use the same real Pt input/physics for all paired rows.

Mandatory:

```text
r4:
    N ~1152
    M=500
    NE=2510
    lld=150

r6:
    N ~3888
    M=500
    NE=2510
    lld=150

r8:
    N ~9216
    M=500
    NE=2510
    lld=150
```

Mandatory physics:

```text
spin
per_type
```

At r4 also run:

```text
charge
orbital
```

where already validated.

Do not change the physical parameters between CPU/GPU paired rows.

---

## G1.2-K — Repetitions and statistics

For moderate rows use:

```text
2 warmups
5 measured repetitions
```

in a persistent process where practical.

For very expensive CPU FP64 r8 runs, three measured repetitions are acceptable
if documented.

Report:

```text
median
min
max
MAD or IQR
```

Do not publish a single lucky run.

GPU context initialization must be reported separately from persistent
transport timing.

---

## G1.2-L — Correctness attached to performance

For each precision-matched pair record:

```text
moment tensor:
    max abs
    relative L2

conductivity spectrum:
    max abs
    relative L2

integrated/final conductivity:
    absolute/relative differences
```

Use identical `per_type` projectors.

For later `random_vec`, use identical seeds/vectors.

No speedup claim for a failed numerical row.

---

## G1.2-M — Required report

Produce a table like:

```text
material
size
N
nnz
M
NE
lld
cond_type
estimator
moment_backend
moment_precision
reconstruction_backend
reconstruction_precision
OMP threads
BLAS threads

P_operator
P_trace_setup
P_moments_total
  D_H2D
  D_GPU_kernel
  D_D2H
  D_conversion
P_gamma
P_reconstruction_total
  D_mu_pack
  D_BLAS
P_energy_integration
P_output_io
P_other
T_transport_total
whole_wall

profile_closure_error
correctness_status
```

Then publish three speedups:

```text
S_moments
S_transport
S_whole
```

against the best same-precision CPU row.

---

## G1.2-N — Decision output for G1.3

The final G1.2 report must answer:

1. What is the true CPU FP64 moment time?
2. What is the true CPU FP32 moment time?
3. What is the true GPU FP64 moment total and kernel-only time?
4. What is the true GPU FP32 moment total and kernel-only time?
5. What is the best CPU thread count for r4/r6/r8?
6. What fraction is Gamma generation?
7. What fraction is BLAS reconstruction?
8. What fraction is integration/output?
9. What is the equal-precision CPU/GPU end-to-end speedup?
10. What data movement would be removed by a GPU-resident reconstruction?

Do not proceed to a GPU reconstruction based on impossible/overlapping timing
rows.

---

## G1.2 checklist

- [ ] current timer scopes audited
- [ ] reason for old overlapping totals documented
- [ ] exclusive top-level phase timers implemented
- [ ] nested/detail timers clearly identified
- [ ] CPU/GPU moment parent boundaries made equivalent
- [ ] backend metadata emitted
- [ ] moment precision emitted
- [ ] reconstruction backend emitted
- [ ] reconstruction precision emitted
- [ ] OMP thread count recorded
- [ ] BLAS thread count recorded
- [ ] integration separated from output I/O
- [ ] Gamma bytes recorded
- [ ] packed-mu bytes recorded
- [ ] profiler closure check implemented
- [ ] child<=parent consistency checks implemented
- [ ] CPU OMP 1/2/4/8 sweep run
- [ ] FP64 equal-precision rows run
- [ ] FP32 equal-precision rows run where supported
- [ ] mixed routes separately labelled
- [ ] r4 M=500 benchmark complete
- [ ] r6 M=500 benchmark complete
- [ ] r8 M=500 benchmark complete
- [ ] charge/spin/orbital r4 coverage complete
- [ ] repeated-run statistics recorded
- [ ] correctness attached to each performance pair
- [ ] new bottleneck identified from valid exclusive timing

**Commit message:**

`Make KPM transport profiling precision-fair`

---

# KPM-G1.3 — GPU reconstruction using the optimized Gamma*mu GEMM track

## Objective

Make the GPU transport path follow the **same optimized mathematical
reconstruction introduced by G1.1**, rather than the obsolete scalar
Gamma*mu loops.

The CPU reference now evaluates the contraction as a flattened matrix product:

```text
G(NE, M*M) * U(M*M, nb) -> C(NE, nb)
```

G1.3 must implement the GPU analogue using cuBLAS and, where practical, keep
the stochastic moments resident so that the large `mu_nm` object is not
downloaded merely to be multiplied on the CPU.

This is NOT a task to invent a new conductivity formula.

This is NOT a task to port the old scalar Fortran loops to CUDA.

This is NOT G2 stochastic-vector batching.

---

## G1.3-A — Preconditions

Before editing:

1. read the G1.2 final report;
2. verify profiler closure PASS for the benchmark rows used to motivate G1.3;
3. record the best CPU FP64 and FP32 BLAS reconstruction baselines;
4. record current GPU FP64/FP32 moment performance;
5. record D2H bytes/time for `mu_nm`.

If G1.2 shows that GPU reconstruction cannot materially improve the target
workload, document that before implementing.

The user nevertheless wants the GPU path to follow the optimized CPU track, so
a scoped implementation is still allowed for architectural consistency, but
its performance status must remain evidence-based.

---

## G1.3-B — Preserve the G1.1 algebra exactly

Use the same flattening convention established and tested in G1.1.

Let:

```text
K = M*M
q = flatten(n,m) using the exact validated Fortran storage convention

G(i,q)   = Gamma_nm(i,n,m)
U(q,l)   = mu_nm_stochastic(l,l,n,m,t)

C(i,l)   = factor * sum_q G(i,q) * U(q,l)
```

The GPU implementation must produce this same C.

No unintended conjugation.

Do not use a Hermitian/symmetric BLAS routine unless the mathematical identity
is independently proven and tested.

The primary implementation should use:

```text
cuBLAS ZGEMM-equivalent for FP64
cuBLAS CGEMM-equivalent for FP32
```

through the existing project CUDA style.

Do not expose raw CUDA pointers to Fortran physics code.

---

## G1.3-C — Do not upload/materialize the full 10-GB Gamma tensor on the GPU

For:

```text
NE = 2510
M  = 500
```

the full complex-FP64 Gamma tensor is approximately 10 GB.

Do NOT implement:

```text
CPU constructs full Gamma
-> upload all 10 GB
-> one GPU GEMM
```

as the production solution.

That would consume excessive device memory and transfer bandwidth and would
not scale well.

Use **energy tiling**.

Choose an energy block size:

```text
BE = number of energies in one Gamma block
```

and process:

```text
for each energy block:
    construct Gamma_block(BE,K)
    GEMM Gamma_block * U
    produce C_block(BE,nb)
```

Benchmark several BE values such as:

```text
16, 32, 64, 128
```

subject to memory.

Do not hard-code the first tested value as universal.

---

## G1.3-D — Generate Gamma blocks from the exact production formula

Use the existing production Gamma definition exactly.

The current algebra is built from:

```text
wscale(E)
sqrt(1-wscale^2)
acos(wscale)
Chebyshev T_n(wscale)
cn(E,n)
cm(E,m)
kernel(n)
weights(n)
```

and the established expression:

```text
Gamma(E,n,m) =
    [ cn(E,n) * T_m(E)
    + cm(E,m) * T_n(E) ]
    / (1-wscale(E)^2)^2
    * kernel(n) * kernel(m)
    * weights(n) * weights(m)
```

Verify this against CURRENT HEAD before coding.

For G1.3 implement one of these, in order of preference:

### Preferred: GPU block generation

Precompute/upload only the small 1-D/2-D basis data needed for the current
energy block and generate `Gamma_block` directly on device.

This avoids:
- the full host Gamma allocation for the GPU route;
- ~10 GB H2D transfer;
- a second full Gamma copy.

### Acceptable first diagnostic route

CPU-generate a contiguous Gamma block and upload that block.

This may be used to validate the cuBLAS reconstruction independently.

It is NOT the final preferred production route if repeated host block
generation/H2D dominates.

Do not change the Gamma mathematics for speed.

Do not apply the analytic separability/factorization beyond block generation
in this task. A more radical factorized-Gamma algorithm is a separate future
WP if profiling justifies it.

---

## G1.3-E — Keep GPU stochastic moments resident

The current GPU stochastic-moment engine already computes `mu_nm` on device.

Avoid:

```text
GPU computes mu
-> D2H full mu
-> host packs U
-> H2D U
-> GPU GEMM
```

That would be architecturally pointless.

Add a narrow backend-owned residency contract for the moment result.

The CUDA backend/context may retain:

```text
device_mu
precision
M
nb
trace/type index
system generation/token
valid flag
```

Fortran must not receive a raw device pointer.

Provide backend methods conceptually like:

```text
stochastic_moments(..., retain_device=.true.)

reconstruct_conductivity_from_resident_moments(...)
```

or an equivalent design consistent with current backend style.

A new Hamiltonian/operator/precision/shape must invalidate stale residency.

Do not create a general-purpose GPU object registry.

---

## G1.3-F — Pack only the needed diagonal moments on GPU

The G1.1 BLAS contraction uses:

```text
mu(l,l,n,m,t)
```

for the reconstruction.

Do not transfer or unnecessarily duplicate the complete orbital matrix tensor
if only the diagonal blocks are required by the current conductivity formula.

On device, build/reuse:

```text
U(K,nb)
```

with the same flattening as CPU G1.1.

For:

```text
M=500
nb=18
```

this is modest relative to the full Gamma tensor.

Pack on GPU from the resident moment representation.

Profile:

```text
D_gpu_mu_pack
```

separately.

If the CUDA moment engine can naturally emit the required packed diagonal form
without changing canonical host semantics, that may be used, but do not
contort G2/vector batching into this task.

---

## G1.3-G — Precision routes

Implement and validate two explicit reconstruction routes.

### GPU FP64

```text
FP64 moments
FP64 Gamma block
FP64 packed U
cuBLAS complex-double GEMM
FP64 C_block
```

Return canonical FP64 host result.

### GPU FP32

```text
FP32 moments
FP32 Gamma block
FP32 packed U
cuBLAS complex-single GEMM
FP32 C_block
```

Widen the small final `C(NE,nb)` result to FP64 on host for the existing
integration/output layer.

Label this:

```text
FP32 reconstruction with FP64 canonical output
```

not "FP64 reconstruction".

Optional mixed route:

```text
FP32 moments -> FP64 reconstruction
```

may be retained for scientific testing only if useful, but must be labelled
MIXED and not replace the two like-for-like routes.

---

## G1.3-H — Return only the reduced reconstruction result by default

For `per_type`, after the GPU GEMM the host needs only the reduced result
needed to populate the existing integrand semantics:

```text
C(NE,nb)
```

or type-resolved equivalents.

For FP64 r4 this is orders of magnitude smaller than full `mu_nm`.

The normal optimized GPU route should therefore transfer:

```text
C_block / C
```

rather than full moments.

However:

- preserve an explicit diagnostic/download path for full `mu_nm`;
- current validation tests must still be able to compare moments;
- do not delete canonical host-moment support.

Record D2H bytes before and after.

---

## G1.3-I — Preserve per_type semantics exactly

For each type `t`, the GPU result must reproduce the G1.1 CPU semantics:

```text
integrand_at(:, :, :, t)
```

and total:

```text
integrand = sum_t integrand_at(...,t)
```

with the same factor and later normalization.

Do not change the order/meaning of type-resolved output files.

For `random_vec`, if the existing code simply accumulates trace contributions,
support the existing scalar trace semantics.

Do not implement G2 block stochastic vectors here.

---

## G1.3-J — Keep energy integration/output on CPU for this WP

Once `C(NE,nb)` is downloaded, use the existing canonical CPU:

```text
energy integration
orbital-resolved integration
file output
```

unchanged.

Do not GPU-port Simpson integration or file/output work in G1.3.

The purpose is to isolate the benefit of:
- resident moments;
- tiled Gamma generation;
- cuBLAS reconstruction.

---

## G1.3-K — Correctness ladder

### 1. Tiny deterministic reconstruction test

Use small:

```text
NE
M
nb
```

with arbitrary complex basis data/moments.

Compare:

```text
scalar CPU reference
CPU BLAS G1.1
GPU FP64 tiled reconstruction
GPU FP32 tiled reconstruction
```

Verify:
- flattening;
- no conjugation error;
- block-boundary correctness.

Choose NE not divisible by BE to test the final partial block.

### 2. Gamma-block correctness

For small/moderate M compare GPU-generated Gamma blocks with the production CPU
`calculate_gamma_nm` output.

Report:
- max abs error;
- relative L2.

FP64 should agree at FP64 numerical accuracy.

FP32 should be judged against an FP32-rounded/reference contract and also
reported against FP64 for scientific precision.

### 3. Pt M=500

Mandatory:

```text
r4
M=500
NE=2510
lld=150
spin
per_type
```

Compare full conductivity spectrum against:
- CPU FP64 BLAS reference;
- CPU FP32 BLAS reference for FP32 route.

### 4. Larger real systems

Run r6 and r8.

The reconstruction itself should be largely N-independent; this verifies that
the integrated GPU route does not regress with system size/memory pressure.

### 5. Charge/orbital

At r4 validate charge and orbital conductivity where established.

---

## G1.3-L — Benchmark energy block size

For the mandatory r4 case benchmark:

```text
BE = 16
BE = 32
BE = 64
BE = 128
```

and a larger value if clearly memory-safe/useful.

Report:

```text
Gamma block generation
Gamma block memory
cuBLAS GEMM
result D2H
total reconstruction
```

Choose a conservative default from evidence.

Do not benchmark only one tile and call it optimal.

---

## G1.3-M — Fair CPU/GPU performance comparison

Use G1.2's best same-precision CPU results.

Headline comparisons:

```text
best CPU FP64 BLAS transport
vs
GPU FP64 moments + GPU FP64 reconstruction

best CPU FP32 BLAS transport
vs
GPU FP32 moments + GPU FP32 reconstruction
```

Also report the intermediate route:

```text
GPU moments + CPU BLAS reconstruction
```

This isolates the incremental benefit of G1.3.

For every row report:

```text
P_moments_total
P_gamma / Gamma basis preparation
P_reconstruction_total

GPU detail:
    moment H2D
    moment kernel
    moment D2H avoided
    GPU mu pack
    Gamma block generation
    Gamma block H2D if any
    cuBLAS GEMM
    reduced result D2H

P_energy_integration
P_output_io
T_transport_total
whole wall
```

Profiler closure must PASS.

---

## G1.3-N — Performance acceptance

Do not require a predetermined speedup.

Retain the GPU reconstruction as production-beneficial only if it improves
realistic end-to-end transport or provides a clear memory-transfer advantage
without material regression.

The task should explicitly calculate:

```text
S_reconstruction =
    CPU_BLAS_reconstruction / GPU_reconstruction

S_transport =
    best_same_precision_CPU_transport / GPU_transport

S_whole =
    best_same_precision_CPU_wall / GPU_wall
```

Also record:

```text
full-mu D2H bytes before
reduced-result D2H bytes after
```

If GPU reconstruction only saves milliseconds because CPU BLAS is already
~1.9 s and Gamma generation/integration dominate, state that honestly.

Architectural parity alone is not enough to call it a major speedup.

---

## G1.3-O — Memory accounting

For M=500 report peak device allocations for:

```text
Hamiltonian/operators
resident moments
packed U
Gamma basis arrays
Gamma_block
C_block
cuBLAS workspace if applicable
other stochastic workspace
```

Do not approach device OOM without preflight.

Do not materialize the full FP64 Gamma tensor on a 16-GB A4000.

---

## G1.3-P — Do not implement factorized Gamma algebra yet

Although the production Gamma formula may admit a lower-rank/separable
reformulation, G1.3 must stay on the already validated G1.1 algebraic track:

```text
Gamma block
*
packed moment matrix
```

This keeps CPU and GPU mathematically aligned.

If profiling after G1.3 shows:

```text
Gamma block generation dominates
```

then propose a later focused task to derive and validate a factorized
Gamma contraction.

Do not mix that derivation into this WP.

---

## G1.3 checklist

- [ ] G1.2 profiler closure passes before optimization
- [ ] CPU G1.1 flattening contract reused exactly
- [ ] no scalar GPU Gamma*mu kernel introduced
- [ ] cuBLAS FP64 GEMM path implemented
- [ ] cuBLAS FP32 GEMM path implemented
- [ ] full 10-GB Gamma upload avoided
- [ ] energy tiling implemented
- [ ] BE=16 benchmarked
- [ ] BE=32 benchmarked
- [ ] BE=64 benchmarked
- [ ] BE=128 benchmarked
- [ ] final partial energy block tested
- [ ] Gamma block formula matches production CPU
- [ ] GPU moment residency implemented narrowly
- [ ] residency invalidation is explicit
- [ ] no raw device pointer leaks into Fortran physics
- [ ] packed diagonal U generated on GPU
- [ ] full mu D2H avoided in optimized route
- [ ] diagnostic full-moment download retained
- [ ] reduced C result returned to canonical CPU integration
- [ ] per_type total semantics preserved
- [ ] per_type resolved semantics preserved
- [ ] random_vec scalar semantics preserved
- [ ] FP64 tiny deterministic test passes
- [ ] FP32 tiny deterministic test passes
- [ ] FP64 Gamma-block test passes
- [ ] FP32 Gamma-block error quantified
- [ ] Pt r4 M=500 spin FP64 passes
- [ ] Pt r4 M=500 spin FP32 passes
- [ ] Pt r6 tested
- [ ] Pt r8 tested
- [ ] charge/orbital r4 validated where supported
- [ ] same-precision CPU/GPU performance table produced
- [ ] GPU-moments+CPU-BLAS intermediate row retained
- [ ] D2H byte reduction reported
- [ ] peak device memory reported
- [ ] profiler closure passes after G1.3
- [ ] no energy integration GPU port introduced
- [ ] no factorized-Gamma redesign introduced
- [ ] unfavorable results retained

**Commit message:**

`Add GPU BLAS Kubo-Bastin reconstruction`

---

# Recommended sequence after G1.3

After G1.2 and G1.3:

```text
G1.2
    trustworthy timing + fair CPU/GPU baselines
        |
        v
G1.3
    resident GPU moments
    + tiled Gamma
    + cuBLAS reconstruction
        |
        v
reprofile M=500 r4/r6/r8
        |
        +-- Gamma generation dominates
        |      -> design a focused factorized/tiled-Gamma optimization
        |
        +-- stochastic trace throughput is next concern
        |      -> KPM-G2 block stochastic traces
        |
        +-- integration/output dominates
               -> optimize only if realistic wall fraction justifies it
```

Do not decide the next bottleneck from the pre-G1.2 overlapping timers.
