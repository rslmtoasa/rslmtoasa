# KPM-G1.4 — Verify and minimize the fixed host-side floor of the optimized GPU transport path

You are working on the **CURRENT HEAD** of:

```text
https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3
```

This task follows KPM-G1, G1.1, G1.2, and G1.3.

Do **not** restart those work packages.

---

# 0. Established evidence entering G1.4

The current realistic production-like benchmark is real SOC Pt transport with:

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

and real-space system sizes:

```text
r4: N = 1152
r6: N = 3888
r8: N = 9216
```

## CPU post-G1.1 evidence

The optimized CPU FP64 route has:

| Size | P_moments_total | P_gamma | P_reconstruction_total | P_energy_integration | P_output_io | P_other |
|---|---:|---:|---:|---:|---:|---:|
| r4 | 10.386 s | 8.198 s | 1.939 s | 0.001 s | 2.870 s | 0.584 s |
| r6 | 34.359 s | 8.202 s | 1.944 s | 0.001 s | 2.852 s | 0.716 s |
| r8 | 82.080 s | 8.211 s | 1.943 s | 0.001 s | 2.858 s | 1.110 s |

The important conclusion is already established:

```text
CPU moment generation is the scaling bottleneck.

Gamma generation, BLAS reconstruction, and output I/O are approximately
independent of real-space system size N at fixed M and NE.
```

## GPU post-G1.3 evidence

The optimized GPU path now uses:

```text
GPU stochastic Chebyshev moments
        ->
resident diagonal-packed U(M*M,nb)
        ->
GPU tiled Gamma generation
        ->
cuBLAS reconstruction
        ->
reduced C(NE,nb) D2H
        ->
existing CPU integration/output
```

Representative FP32 results:

| Size | P_moments_total | P_reconstruction_total | T_transport_total |
|---|---:|---:|---:|
| r4 | 0.439 s | 0.021 s | 5.225 s |
| r6 | 1.022 s | 0.020 s | 6.190 s |
| r8 | 2.248 s | 0.021 s | 7.751 s |

GPU tiled Gamma generation is approximately:

```text
~0.7 s
```

for the production-order workload.

Therefore the remaining GPU transport floor is no longer an unknown
large numerical kernel.

A plausible current decomposition is approximately:

```text
GPU Gamma                  ~0.7 s
output I/O                 ~2.8-2.9 s
setup/postprocess/misc     ~0.5-1.0 s
```

which already explains most or all of the residual 4-5 s.

This must now be **verified**, not rediscovered from scratch.

---

# 1. Objective

The purpose of G1.4 is:

> **Verify and minimize the fixed host-side floor of the optimized GPU KPM transport route.**

Specifically:

1. prove where the remaining ~4-5 s comes from;
2. verify that no obsolete CPU work survives on the CUDA route;
3. remove only demonstrably redundant or materially optimizable work;
4. quantify the unavoidable output/I/O floor;
5. stop when the remaining numerical stages are too small to justify further pre-G2 work.

This is **not** another broad GPU-porting campaign.

This is **not** KPM-G2.

---

# 2. Non-negotiable rules

## 2.1 Inspect current HEAD first

Read at minimum:

```text
source/kpm_profile.f90
source/conductivity.f90
source/recursion_transport.f90
current CUDA KPM/transport implementation
G1.2 profiling report/evidence
G1.3 report/evidence
current benchmark harness for r4/r6/r8
```

Do not assume the source still matches old task text.

## 2.2 Preserve physics

Do not change:

```text
Hamiltonian
current operators
charge/spin/orbital conventions
Chebyshev scaling
Gamma definition
Kubo-Bastin prefactors
energy grid
integration rule
normalization
per_type semantics
random_vec semantics
output meaning
```

## 2.3 G1.2 timing semantics remain mandatory

Use only exclusive top-level phase timings for bottleneck decisions.

Nested detail timers must remain clearly marked as children.

Require profiler closure before interpreting a row.

## 2.4 Precision must be explicit

Every benchmark row must identify:

```text
moment backend
moment precision
reconstruction backend
reconstruction precision
```

Do not attribute a speedup to GPU when it arises from changing precision.

## 2.5 CPU threading must remain fair

For any material host-side numerical stage benchmark:

```text
OMP_NUM_THREADS = 1, 2, 4, 8
```

Control BLAS threading and avoid oversubscription.

## 2.6 Do not start G2

Do not implement block stochastic/random-vector processing here.

## 2.7 Negative results are valid

If the residual is mostly required I/O plus small host setup, conclude that
G1.4 is complete and proceed to G2.

Do not invent more work merely to reduce the wall time further.

---

# 3. G1.4-A — Confirm the current exclusive GPU time budget

Run the current optimized GPU route on:

```text
Pt r4
Pt r6
Pt r8

M=500
NE=2510
lld=150
spin
per_type
```

Primary route:

```text
GPU FP32 moments
GPU FP32 reconstruction
```

Also run:

```text
GPU FP64 r4
```

as a reference.

For r4 also run charge and orbital if already validated and available without
changing the fixture.

Use the G1.2 exclusive timing model.

At minimum report:

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

Require:

```text
profile_closure_error <= 0.03
```

for every row used in conclusions.

Any `P_misc > 5%` of `T_transport_total` must be decomposed further.

Do not accept a generic 4-second `P_misc`.

---

# 4. G1.4-B — Verify the expected fixed-floor model

For r4/r6/r8, explicitly test whether:

```text
fixed_floor ~=
    P_gamma_basis_setup
  + P_gamma_generation
  + P_reconstruction_total
  + P_result_unpack
  + P_energy_integration
  + P_tensor_postprocess
  + P_output_prepare
  + P_output_io
  + P_misc
```

is approximately independent of real-space system size N.

Produce:

| Size | GPU moments | GPU Gamma | GPU reconstruction | host postprocess | output I/O | other | total |
|---|---:|---:|---:|---:|---:|---:|---:|

Also produce:

| Size | CPU fast FP64 moments | CPU post-moment fixed work | GPU FP32 moments | GPU post-moment fixed work |
|---|---:|---:|---:|---:|

The purpose is to establish clearly:

```text
CPU runtime increasingly follows real-space moment generation.

GPU runtime increasingly approaches a nearly N-independent host-side floor.
```

If the data do not support this, explain why.

---

# 5. G1.4-C — Mandatory audit for obsolete legacy CPU work

This is the highest-priority correctness/performance audit.

For the optimized CUDA reconstruction path, verify that none of the following
are executed unnecessarily:

```text
full host gamma_nm(NE,M,M) allocation
full host gamma_nm fill
CPU calculate_gamma_nm for an unused tensor
CPU mu packing
CPU ZGEMM reconstruction
full moment download
legacy scalar Gamma*mu loop
legacy integrand reconstruction pass
unused large temporary zero/fill
```

For each item state explicitly:

```text
executed? yes/no
required? yes/no
time
memory
action
```

## Especially important: full host Gamma

At:

```text
NE = 2510
M  = 500
```

a full complex-FP64 Gamma tensor is about 10 GB.

The G1.3 CUDA route generates Gamma in tiles on device.

Therefore the optimized GPU route must not allocate/fill a full legacy host
Gamma tensor unless there is a real remaining consumer.

If such work still exists:

1. prove it is redundant;
2. bypass it only on the CUDA tiled-reconstruction route;
3. preserve it for CPU/reference routes where needed;
4. remeasure immediately.

Do not remove reference functionality globally.

---

# 6. G1.4-D — Reconfirm resident-moment behavior

For the optimized CUDA route verify:

```text
full moment D2H bytes = 0
```

unless explicit diagnostic mode is requested.

Confirm that:

```text
gpu_moment_download = false
```

or its current equivalent does not cause an accidental host copy.

Also verify that no host moment buffer is allocated/zeroed at full production
size when the optimized resident route does not need it.

Report:

```text
resident device moment bytes
host full-moment bytes allocated
full-moment D2H bytes
reduced-result D2H bytes
```

The diagnostic full-moment download path must remain available for validation.

---

# 7. G1.4-E — Quantify output I/O as a real performance floor

The optimized CPU profile already shows:

```text
P_output_io ~2.85 s
```

almost independently of r4/r6/r8.

This is potentially the largest single fixed stage of the FP32 GPU route.

Separate:

```text
P_output_prepare
P_output_io
```

and, if needed for diagnosis, further split:

```text
formatting
buffer preparation
open/close
write
flush
```

## Diagnostic no-write timing

Add a **benchmark-only diagnostic mode** that executes all numerical work and
output preparation but suppresses the actual filesystem write.

This must NOT become a silent production shortcut.

Its purpose is only to measure:

```text
T_numerical_transport
```

versus:

```text
T_numerical_plus_IO
```

The normal production path must continue to write identical output.

Report:

```text
wall with normal output
wall with benchmark no-write
I/O-attributable difference
```

for r4/r6/r8 FP32.

This is important for interpreting GPU speedup.

Do not quote no-write wall time as normal application runtime.

---

# 8. G1.4-F — Optimize I/O only if the implementation is inefficient

If I/O is ~2.8-2.9 s because of genuinely required data volume, record it as an
I/O floor.

If instead the time comes from inefficient mechanics such as:

```text
many tiny formatted writes
repeated open/close
repeated format conversion
unnecessary flushes
temporary string construction
```

then optimize the implementation without changing:

```text
file names
file contents
numerical precision
ordering
format
```

Preferred allowed changes:

```text
buffer larger output blocks
reuse open units when semantics permit
reduce repeated formatting work
avoid unnecessary flush/open/close cycles
```

Benchmark before/after.

Do not reduce the amount of scientific output to claim a speedup.

---

# 9. G1.4-G — Inspect the remaining ~0.5-1 s host overhead

After Gamma, reconstruction, and I/O are accounted for, inspect the remaining
host-side fixed work.

Likely categories:

```text
operator/setup
trace/projector setup
result unpack
tensor assembly
normalization
allocation/deallocation
large zero/fill
precision conversion
CUDA synchronization
energy/Chebyshev basis preparation
```

For each stage above 5% of total:

1. measure it explicitly;
2. classify whether it depends on N, M, NE, conductivity type, or only the call;
3. determine whether it is repeated unnecessarily.

Do not spend engineering effort on sub-5% stages unless the fix is trivial.

---

# 10. G1.4-H — Audit invariant table reuse

Check whether these or equivalent quantities are rebuilt every transport call:

```text
energy-grid transformation
Chebyshev T_n(E) tables
kernel weights
quadrature/Simpson weights
Gamma basis arrays
normalization factors
type/projector maps
```

Classify their validity dependencies:

```text
M
NE
energy range
kernel choice
weights
Hamiltonian scaling bounds
conductivity type
```

If a table is safely reusable across repeated transport calls with identical
parameters, implement a **narrow persistent cache/reuse mechanism** with
explicit invalidation.

Do not add a generic caching framework.

Do not retain large data that no longer provides a measurable benefit.

---

# 11. G1.4-I — CPU optimization policy for remaining host numerical stages

If a remaining host numerical stage is material, optimize CPU first unless
the GPU handoff is clearly superior.

For each candidate run:

```text
OMP_NUM_THREADS = 1
OMP_NUM_THREADS = 2
OMP_NUM_THREADS = 4
OMP_NUM_THREADS = 8
```

Use controlled:

```text
OMP_PROC_BIND
OMP_PLACES
BLAS_NUM_THREADS
```

Prefer:

```text
better memory traversal
SIMD/vectorization
reuse of invariant values
OpenMP over sufficiently large independent loops
```

before adding another GPU kernel.

Do not parallelize tiny loops.

---

# 12. G1.4-J — Optimization thresholds

Use these decision thresholds.

## High priority

A stage is a strong G1.4 optimization target if:

```text
>20% of T_transport_total
```

for at least one realistic FP32 r4/r6/r8 workload.

## Investigate if cheap

```text
10-20%
```

especially if the stage grows badly with M or NE.

## Normally stop

```text
<10%
```

unless:
- the fix is trivial;
- work is demonstrably redundant;
- it blocks G2/B0.

After each retained optimization:

1. rerun r4;
2. validate numerically;
3. reprofile;
4. reconsider the ranking.

Do not implement a predetermined list of optimizations.

---

# 13. G1.4-K — Explicit stop condition

G1.4 should terminate once all of the following are true:

```text
1. no obsolete CPU Gamma/reconstruction work remains on the CUDA route;

2. full moment D2H remains zero in normal optimized execution;

3. profiler closure passes;

4. P_misc <5% or is fully explained;

5. output I/O is separately quantified;

6. no remaining individual numerical host stage exceeds ~10-15% of realistic
   transport total, OR such a stage has been shown not worth optimizing;

7. the residual fixed floor is explainable.
```

In particular:

> If output I/O accounts for roughly half or more of the remaining fixed floor
> and all remaining numerical stages are individually small, STOP G1.4.

Do not GPU-port more postprocessing merely because a few seconds remain.

Proceed to G2.

---

# 14. G1.4-L — Correctness requirements

Every retained code change must preserve:

```text
charge conductivity
spin conductivity
orbital conductivity
type-resolved outputs
energy-resolved spectra
integrated outputs
```

where applicable.

Mandatory checks:

## Unit/surface tests

All existing KPM transport and CUDA surface tests remain green.

## Pt r4 M=500

Compare before/after G1.4:

```text
max absolute spectrum difference
relative L2 spectrum difference
integrated conductivity difference
```

For pure allocation/I/O/reuse optimizations, expect no scientific change beyond
roundoff.

## r6/r8

Verify selected integrated values/checksums and no drift.

Do not loosen tolerances.

---

# 15. G1.4-M — Benchmark protocol

Primary route:

```text
GPU FP32
Pt
spin
per_type
M=500
NE=2510
lld=150
```

Sizes:

```text
r4
r6
r8
```

Also run at r4 where supported:

```text
charge
orbital
```

Reference:

```text
GPU FP64 r4
```

and one larger FP64 row if practical.

For r4:

```text
2 warmups
5 measured repetitions
```

For r6/r8:

```text
2 warmups
>=3 measured repetitions
```

Use the same A4000 device for paired before/after comparisons.

Do not mix device 0 and device 1 within a microbenchmark comparison unless
explicitly studying device-to-device variance.

Report:

```text
median
min
max
MAD or IQR
```

---

# 16. G1.4-N — Required final tables

## Table A — exclusive GPU phase profile

```text
Size
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
profile_closure_error
```

for r4/r6/r8 FP32.

## Table B — fixed-floor interpretation

```text
Size
GPU moments
GPU fixed post-moment work
fraction moments
fraction fixed floor
```

## Table C — CPU versus GPU scaling context

```text
Size
best CPU fast FP64 moments
best CPU fast FP64 total
GPU FP32 moments
GPU FP32 total
```

This table is contextual only unless precision-matched validation permits a
headline speedup.

## Table D — I/O attribution

```text
Size
normal production wall
benchmark no-write wall
output attributable time
output fraction
```

Do not use no-write runtime as a production headline.

---

# 17. G1.4-O — Required interpretation

The final report must answer explicitly:

1. Was the remaining ~4-5 s GPU floor fully explained?
2. Was any obsolete CPU Gamma/reconstruction work still executing?
3. Was any full host Gamma tensor still allocated/fill on the CUDA route?
4. Was any full moment D2H still occurring?
5. How much time is genuine output I/O?
6. How much time is output preparation/formatting?
7. What is the largest remaining numerical host-side stage?
8. Does any remaining host numerical stage justify another optimization WP?
9. Is a factorized-Gamma task still justified?
10. Should the project now proceed to G2?

The expected likely outcome is:

```text
GPU moment stage scales with N;
Gamma/reconstruction are already small;
the remaining floor is mostly I/O + small setup/postprocess;
no further single-trace/per_type numerical optimization is justified;
proceed to G2.
```

But do not force that conclusion if the measurements disagree.

---

# 18. Checklist

- [ ] current G1.3 optimized GPU route reproduced
- [ ] r4 profiler closure PASS
- [ ] r6 profiler closure PASS
- [ ] r8 profiler closure PASS
- [ ] residual fixed floor fully decomposed
- [ ] P_misc <5% or decomposed
- [ ] full host Gamma allocation checked
- [ ] full host Gamma fill checked
- [ ] obsolete CPU calculate_gamma_nm checked
- [ ] CPU BLAS reconstruction on CUDA route checked
- [ ] scalar reconstruction on CUDA route checked
- [ ] full moment D2H remains zero
- [ ] unused full host moment allocation checked
- [ ] reduced-result D2H recorded
- [ ] CUDA/context startup excluded from steady transport timing
- [ ] energy/Chebyshev invariant tables audited
- [ ] invariant reuse implemented only where beneficial
- [ ] output preparation separated from output I/O
- [ ] benchmark-only no-write diagnostic implemented
- [ ] normal production output unchanged
- [ ] I/O floor quantified
- [ ] I/O mechanics optimized only if measurably inefficient
- [ ] remaining >5% host stages explicitly timed
- [ ] OMP 1/2/4/8 tested for any material host numerical stage
- [ ] r4 spin correctness passes
- [ ] r6/r8 correctness passes
- [ ] charge/orbital r4 checks pass where supported
- [ ] before/after exclusive profile table produced
- [ ] CPU/GPU scaling context table produced
- [ ] I/O attribution table produced
- [ ] next action explicitly stated
- [ ] no G2 implementation added
- [ ] no physics changed
- [ ] no tolerances loosened
- [ ] no output silently suppressed in production

---

# 19. Commit message

If code optimization is retained:

```text
Reduce fixed KPM transport overhead
```

If profiling shows no worthwhile code optimization:

```text
Characterize optimized KPM transport floor
```

---

# 20. Final mindset

Do not ask:

```text
"What else can we put on the GPU?"
```

Ask:

```text
"Is the remaining wall time still computational waste, or is it now mainly
the unavoidable host/output floor?"
```

G1.1 and G1.3 already solved the major numerical bottlenecks.

G1.4 should be the final cleanup/attribution pass for the single-trace
`per_type` production route.

If the remaining floor is mostly required I/O and small host work, stop and
move to KPM-G2 rather than continuing to optimize diminishing returns.
