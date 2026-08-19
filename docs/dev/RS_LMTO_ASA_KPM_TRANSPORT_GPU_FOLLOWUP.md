# RS-LMTO-ASA Phase III-A Follow-up
## KPM / Kubo–Bastin transport GPU optimization and fair benchmark campaign

**Target branch:** `fable_v3`

**Purpose:** establish the true GPU performance of the real-space KPM/Kubo–Bastin transport implementation under realistic production workloads, improve the stochastic-moment GPU path where profiling justifies it, and produce like-for-like CPU/GPU evidence suitable for scientific and performance discussions.

---

# 0. Why this follow-up is needed

ACC-13 established that the existing CUDA stochastic-moment machinery can support charge, spin, and orbital Kubo–Bastin transport without duplicating the downstream physics. Its correctness evidence was useful, but its published performance probe was too small and too mixed in precision to answer the production-performance question.

The previous small benchmark used values such as:

```text
cond_ll = 20 / 40
```

which are not representative of serious transport calculations.

A mandatory production-like reference for the new campaign is:

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

The new work must also distinguish:

1. **the stochastic Chebyshev/Kubo moment kernel speedup**;
2. **the complete Kubo–Bastin transport speedup**;
3. **the whole executable speedup**.

Do not use one number for all three.

The final benchmark must use real physical Hamiltonians and realistic system sizes, and CPU/GPU performance claims must be **like-for-like in numerical precision**.

---

# 1. Precision policy for all work packages

This campaign must explicitly distinguish the following routes.

## Scientific FP64 comparison

```text
CPU FP64
vs
GPU FP64
```

Both must use:
- identical physical H;
- identical Chebyshev scaling/bounds;
- identical current operators;
- identical moment order;
- identical energy grid;
- identical estimator/projector/random vectors;
- identical physical prefactors.

This is the principal equal-precision scientific performance comparison.

## Optimized FP32 comparison

```text
CPU FP32
vs
GPU FP32
```

Again, all physics and stochastic realizations must be identical.

This is the principal equal-precision performance comparison.

## Mixed production route

If the existing normal production route remains:

```text
CPU scientific arrays FP64
GPU recurrence/moments FP32
host output widened/accumulated to FP64
```

benchmark it separately and label it **mixed precision**.

Do not compare:

```text
CPU legacy FP64
vs
GPU FP32
```

and call that the FP32 GPU speedup.

It may remain a useful practical production comparison, but it is not a like-for-like precision benchmark.

---

# 2. Real-system benchmark philosophy

Synthetic sparse matrices may be used for low-level kernel tests only.

All performance conclusions must include real production RS-LMTO Hamiltonians.

Mandatory physics classes:

1. **SOC Pt transport**
   - use the existing validated Pt transport fixture/state;
   - spin conductivity is mandatory;
   - charge and orbital conductivity should also be covered where already supported.

2. **Magnetic metallic transport**
   - use an existing validated real magnetic material in the repository, preferably the current bcc-Fe/spd reference if it supports the required transport route cleanly;
   - retain the production magnetic/SOC conventions.

3. **Larger real-space workload**
   - use validated periodic replication/supercell construction of a real material;
   - prefer the existing Pt transport state when a clean replicated SOC workload can be constructed;
   - use several increasing real-space sizes rather than one hand-selected favorable case.

If the repository already contains a larger multi-type SOC/interface transport example, it may be included as an additional production workload.

Do not invent a toy material merely to obtain favorable sparsity.

---

# 3. Mandatory transport parameters

The new campaign must include at least one row with:

```fortran
cond_ll = 500
lld = 150
recur = 'chebyshev'
cond_type = 'spin'
cond_calctype = 'per_type'
```

using a real SOC material.

For scaling analysis, use physically meaningful moment orders such as:

```text
cond_ll = 250, 500, 750
```

and add `1000` only if memory/runtime is practical.

The main report must not use `cond_ll=20` or `40` as representative production performance.

Small orders may remain microbenchmarks.

Keep `lld=150` fixed for the mandatory production series unless the code semantics or a current physics fixture require another value.

If `lld` has a different meaning in current HEAD than assumed here, inspect and document that before changing the campaign. Do not reinterpret an input keyword from memory.

---

# 4. Stochastic-trace coverage

`cond_calctype='per_type'` is a mandatory user-relevant workflow and must remain fully benchmarked.

In addition, inspect the current code for the supported stochastic/random-vector trace mode.

If a random-vector estimator exists, benchmark its real input keyword and semantics; do not invent a new keyword merely for this campaign.

For stochastic CPU/GPU correctness comparisons:

- use the same random seed;
- use the exact same generated vectors/projectors;
- preserve their distribution and normalization;
- compare the same stochastic realization.

This separates backend numerical error from stochastic sampling error.

Separately, benchmark statistical convergence over multiple seeds when needed for physical validation.

---

# COMMON LUNA HEADER

Prepend this to each implementation/benchmark task:

```text
You are working on CURRENT HEAD of:

https://github.com/rslmtoasa/rslmtoasa/tree/fable_v3

This task is a focused follow-up to the completed ACC accelerator campaign.

The purpose is NOT to prove that GPU is faster.

The purpose is to determine honestly where the KPM/Kubo-Bastin GPU path is
faster, why, at what precision, and for which real production workloads.

============================================================
NON-NEGOTIABLE RULES
============================================================

1. INSPECT CURRENT HEAD FIRST

Read:
- recursion_transport.f90 and its submodules;
- conductivity.f90 and related Kubo-Bastin reconstruction code;
- current CUDA RS-recursion / stochastic-moment implementation;
- ACC-13 report and tests;
- current Chebyshev SCF benchmark harness;
- current CPU fast/legacy Chebyshev implementations;
- precision definitions;
- current transport examples/validation fixtures.

Do not implement from old task wording if HEAD differs.

2. REAL MATERIALS ARE MANDATORY

Synthetic matrices are allowed only for focused unit tests.

Performance conclusions require real production RS-LMTO Hamiltonians.

3. LIKE-FOR-LIKE PRECISION IS MANDATORY

Headline performance comparisons are:

    CPU FP64 vs GPU FP64
    CPU FP32 vs GPU FP32

Mixed-precision routes are labelled separately.

Never silently compare CPU FP64 against GPU FP32 and describe it as an
equal-precision GPU speedup.

4. SAME STOCHASTIC REALIZATION

For a CPU/GPU comparison, use the same stochastic vectors/projectors and the
same random seed.

Do not allow stochastic noise to masquerade as numerical backend error.

5. SAME PHYSICS

Do not modify:
- Hamiltonian;
- current operators;
- spin-current convention;
- orbital-current convention;
- Chebyshev scaling;
- spectral bounds;
- Jackson/kernel factors;
- Kubo-Bastin prefactors;
- energy grid;
- Fermi level;
- broadening;
- stochastic normalization

to improve GPU timing.

6. PRESERVE CANONICAL OUTPUTS

Do not create a second GPU-specific conductivity physics path.

GPU kernels may accelerate numerical stages, but the output must retain the
current canonical moments/conductivity semantics and existing downstream API.

7. NO SILENT FALLBACK

Explicit GPU means GPU.

Unsupported GPU precision/estimator combinations must be reported explicitly.

8. SEPARATE TIMING STAGES

At minimum distinguish:

    setup/operator construction
    stochastic Chebyshev moment stage
    host/device transfers
    Gamma/reconstruction
    energy integration
    total transport phase
    whole executable

Do not call whole-executable speedup "kernel speedup".

9. CPU THREADING MUST BE FAIR

For CPU algorithms benchmark:

    OMP_NUM_THREADS = 1, 2, 4, 8

Use controlled affinity and report it.

Avoid nested BLAS/OpenMP oversubscription.

Report which CPU thread count wins each workload.

GPU speedup should be quoted against the BEST reasonable CPU result for the
same precision and workload, not automatically against one CPU thread.

10. PERFORMANCE RUNS ARE PERSISTENT AND REPEATED

Warm up before timing.

Use multiple measured repetitions in one process where the application
architecture permits.

Report median and spread.

Do not base the final conclusion on one lucky wall-clock sample.

11. DO NOT OVER-ENGINEER

Prefer:
- persistent workspaces;
- block/vector batching;
- existing CUDA libraries;
- narrow kernels;
- current class-like Fortran/CUDA ownership style.

Do not introduce:
- a universal GPU memory manager;
- a generic solver/provider hierarchy;
- a second transport architecture.

12. PROTECTED PHYSICS

Do not modify legacy atomic LMTO routines or unrelated SCF physics.

13. NEGATIVE RESULTS ARE VALID

If a proposed optimization does not improve end-to-end performance:
- revert or do not promote it;
- record the result.

14. ONE FOCUSED COMMIT PER TASK

At completion:
- tick completed boxes only;
- state exact hardware/build;
- list correctness evidence;
- list performance evidence;
- list rejected optimizations;
- make one focused commit using the supplied message.
```

---

# KPM-G0 — Establish realistic transport baselines and a stage-level profiler

## Objective

Before optimizing the transport CUDA implementation, establish a trustworthy
baseline using realistic moment orders and real materials.

This task primarily adds instrumentation and benchmark fixtures.

Do not perform major kernel rewrites.

## A. Map the current transport algorithm

Document the production call chain for:

```text
Hamiltonian construction
current/operator construction
Chebyshev scaling
stochastic/projector state construction
left/right recurrence
mu_nm accumulation
transfer
Gamma_nm(E)
Gamma*mu contraction
energy integration
output
```

For each stage state:
- CPU or GPU;
- precision;
- OpenMP/CUDA library use;
- asymptotic scaling in N, M=cond_ll, NE/lld, and Ntrace where relevant.

Use the current implementation as truth.

## B. Add stage-level timing

Add interval-consistent timers for:

```text
T_operator
T_trace_setup
T_cheb_moments
T_H2D
T_D2H
T_gamma
T_gamma_mu
T_energy_integral
T_transport_total
```

If some stages are naturally inseparable, document that rather than inventing
fake precision.

Also record:
- matrix dimension / orbital dimension;
- sparse nnz;
- cond_ll;
- lld;
- trace/projector count;
- precision;
- bytes transferred.

## C. Mandatory realistic Pt anchor

Use the current validated SOC Pt transport state and run:

```fortran
calctype = 'B'
nsp = 2
cond_ll = 500
lld = 150
recur = 'chebyshev'
cond_type = 'spin'
cond_calctype = 'per_type'
```

Do not reduce `cond_ll` to make the benchmark convenient.

## D. Size scaling

Create or reuse real replicated Pt workloads at several sizes.

Prefer a geometric sequence that spans:
- clearly CPU-small;
- intermediate;
- genuinely large RS workload.

For example, if existing production replication semantics support it cleanly:

```text
replication ~ 4
replication ~ 6
replication ~ 8
```

and one larger case if memory allows.

Record the actual:
- atom count;
- basis/orbital dimension;
- nnz;
- memory footprint.

Do not report only replication factors.

## E. Additional real material

Include at least one real magnetic metal, preferably the existing validated
bcc-Fe reference when compatible with the transport route.

Do not create a new physical validation burden solely for benchmarking.

## F. Baseline precision routes

Measure, where currently available:

```text
CPU FP64 legacy/reference
CPU FP32 optimized
GPU FP32 production
GPU FP64 reference/production-capable route
```

If the GPU FP64 route is only a low-level validator and cannot yet run the
real transport path, record that explicitly; KPM-G1 will address it.

## G. CPU threading sweep

For each CPU route:

```text
OMP_NUM_THREADS=1
OMP_NUM_THREADS=2
OMP_NUM_THREADS=4
OMP_NUM_THREADS=8
```

Control:
- OMP_PROC_BIND;
- OMP_PLACES;
- BLAS/MKL thread count;
- nested parallelism.

Do not oversubscribe.

Report the best CPU result per precision.

## Deliverable

### Current implementation map

The current production call chain is now instrumented by the `KPM_PROFILE`
record emitted by `kpm_profile_mod`. The ownership and precision at CURRENT
HEAD are:

| Stage | Current owner and precision | Scaling/measurement note |
|---|---|---|
| Hamiltonian construction | `prepare_post_processing_stack` / host FP64 | Existing RS-LMTO preprocessing; included in `T_transport_total`, not silently folded into the moment kernel. |
| Current/operator construction | `compute_moments_stochastic` / host FP64 | `build_realspace_velocity_operators` plus charge/spin/orbital specialization; measured as `T_operator`. |
| Chebyshev scaling | `resolve_chebyshev_window` / host FP64 | The same affine `(H-b)/a` window is passed to CPU and CUDA routes. |
| Trace/projector setup | `compute_moments_stochastic` / host FP64 | `per_type` uses one identity block on the selected type; `random_vec` uses the existing normalized random phase block; measured as `T_trace_setup`. |
| Left/right recurrence and `mu_nm` | CPU legacy FP64, CPU `fast` FP32 with FP64 host output, CPU `fast_dp` FP64, or CUDA block-ELL FP32 with FP64 host output by default | Measured as `T_cheb_moments`; current CUDA API reports its arithmetic interval separately from its synchronous transfers. |
| H2D / D2H | CUDA plugin | Operator upload and per-trace `psiref`/`mu_nm` copies are measured separately; byte counters are recorded as `bytes_h2d` and `bytes_d2h`. |
| `Gamma_nm(E)` | `conductivity%calculate_gamma_nm` / host FP64 | `O(NE*M^2)` construction; measured as `T_gamma`. |
| `Gamma*mu` contraction | `conductivity%calculate_conductivity_tensor` / host FP64 + OpenMP | Current `ntype * NE * M^2 * nb` diagonal contraction; measured as `T_gamma_mu`. |
| Energy integration/output | `conductivity%calculate_conductivity_tensor` / host FP64 + host I/O | Simpson integrations and per-type output; measured as `T_energy_integral`. |
| Complete transport phase | `post_processing_conductivity` / mixed according to the selected route | Includes stack construction, moments, reconstruction, and output; measured as `T_transport_total`. |

Here `N = nb * lattice%kk`, and profile `nnz` is the current block-scalar
workload proxy `count(nn(:,2:) > 0) * nb^2`; the first neighbor-list column
stores the per-site neighbor count and is excluded. It is not a claim that
every element of every stored block is numerically nonzero. The existing top-level `Calculation`
timer remains the whole-executable measurement.

The CUDA moment engine's arithmetic and final moment copy are naturally inside
one synchronous API call. The CUDA boundary therefore reports the three
non-overlapping intervals `T_H2D`, `T_cheb_moments`, and `T_D2H` directly; CPU
`T_H2D`/`T_D2H` are zero by construction. This avoids presenting a host-side
wall interval that double-counts a transfer as arithmetic.

The production baseline command is:

```bash
python3 tests/benchmarks/kpm_g0_transport.py \
  --binary build-acc01-cuda/bin/rslmto.x \
  --scratch-root /tmp/rslmto-kpm-g0 \
  --output results/benchmarks/kpm_g0.json \
  --cond-ll 500 --lld 150 --channels 2500 \
  --replications 4 6 8 --omp-threads 1
```

Run the CPU reference route with `--skip-gpu --cheb-backend legacy` for each
controlled `--omp-threads 1`, `2`, `4`, and `8` setting. Repeat with
`--cheb-backend fast` for the optimized CPU FP32 route (or `fast_dp` for the
optimized FP64 route). The driver stages the existing SOC Pt and magnetic
bcc-Fe fixtures and records actual dimensions, block workload, stage
intervals, and transfer bytes. Its smoke mode may use a smaller `M`, but such
output is not production evidence.

A stage-level baseline table such as:

```text
material
size
N
nnz
M
lld
estimator
Ntrace
precision
backend
OMP threads
operator
moments
D2H
gamma
gamma_mu
integration
transport total
whole run
```

## Checklist

- [x] production call chain documented
- [x] precision of every numerical stage documented
- [x] stage-level timers added
- [x] transfer bytes recorded
- [ ] realistic Pt M=500/lld=150 anchor run
- [ ] several real Pt sizes run
- [ ] one magnetic real-material case run
- [ ] CPU FP64 baseline measured
- [ ] CPU FP32 baseline measured if supported
- [ ] GPU FP32 baseline measured
- [ ] GPU FP64 status established
- [ ] CPU OMP 1/2/4/8 sweep complete
- [ ] best CPU comparator recorded
- [x] no major GPU optimization introduced
- [ ] no small-M result presented as production representative

**Commit message:** `Profile realistic KPM transport workloads`

---

# KPM-G1 — Make precision-matched GPU transport and persistent workspaces production-capable

## Objective

Create clean, benchmarkable FP64 and FP32 GPU stochastic-moment paths and
remove repeated allocation/setup overhead from the hot path.

Do not yet redesign the stochastic estimator.

## A. Explicit transport precision strategy

Provide narrow internal strategies such as:

```text
gpu_fp64
gpu_fp32
```

or the closest naming consistent with current backend style.

Do not introduce a global application-wide precision abstraction.

Fortran scientific arrays remain in the current application precision.

The GPU strategy controls the numerical precision of the stochastic
Chebyshev/moment stage.

## B. FP64 GPU path

If current FP64 CUDA support exists only in low-level validation, wire it
through the real production transport call chain with identical physics.

Use:
- complex FP64 state vectors;
- FP64 moment accumulation;
- current FP64 output semantics.

Validate against CPU FP64 before performance claims.

## C. FP32 GPU path

Retain the production FP32 recurrence/moment route.

Document exactly:
- which arrays are FP32;
- which accumulations are FP32;
- which outputs are widened to FP64;
- where conversion occurs.

Do not call widened output "FP64 computation".

## D. Persistent GPU workspaces

Audit the stochastic moment request for repeated:

```text
cudaMalloc/cudaFree
cuSPARSE descriptor creation/destruction
cuBLAS handle creation/destruction
temporary recurrence-vector allocation
left-state allocation
mu_nm allocation
operator staging allocation
timing-event creation/destruction
```

Move repeatedly reusable resources to context/workspace lifetime.

Use capacity-based growth:

```text
if capacity >= required:
    reuse
else:
    grow
```

Do not repeatedly shrink.

## E. Persistent sparse operators

Where H and current operators are unchanged over repeated transport calls,
avoid re-uploading/rebuilding identical sparse metadata within the measured
transport phase.

Preserve explicit invalidation when:
- H changes;
- scaling changes;
- operator changes;
- system dimension/nnz changes.

Do not create stale device-state hazards.

## F. Correctness

For Pt M=500/lld=150 and at least one larger Pt size:

FP64 GPU vs CPU FP64:
- moment tensor;
- charge/spin/orbital conductivity as applicable;
- selected energy-resolved output;
- sum/integral invariants already used in validation.

FP32 GPU vs CPU FP32:
same comparisons.

Also compare FP32 against FP64 scientifically, but label this as a precision
study rather than backend correctness.

## G. Performance

Show before/after for:
- stochastic moment stage;
- transport total.

Use the same real workloads and repeated persistent timing from G0.

## Checklist

- [ ] explicit GPU FP64 transport route available
- [ ] explicit GPU FP32 transport route available
- [ ] FP32 conversion/accumulation semantics documented
- [ ] persistent recurrence workspace implemented
- [ ] persistent moment workspace implemented
- [ ] repeated CUDA allocations removed where material
- [ ] sparse operator upload/reuse audited
- [ ] safe invalidation implemented where reused
- [ ] Pt M=500 FP64 correctness passes
- [ ] Pt M=500 FP32 same-precision correctness passes
- [ ] larger Pt correctness passes
- [ ] pre/post optimization timing recorded
- [ ] no global precision framework introduced

**Commit message:** `Add persistent precision-matched KPM CUDA paths`

---

# KPM-G2 — Block stochastic traces and expose vector-level GPU parallelism

## Objective

Increase GPU utilization for stochastic KPM by processing several independent
trace/projector vectors together where the production estimator permits it.

This task is particularly important for genuine stochastic/random-vector
transport calculations.

Do not alter the estimator mathematically.

## A. Inspect current estimator semantics

Identify all supported `cond_calctype` modes.

For each state:
- what the starting vector/projector represents;
- how many vectors are used;
- whether vectors are independent;
- how averaging/normalization works.

Do not assume `per_type` is a random-vector estimator.

## B. Preserve the Chebyshev recurrence exactly

For scaled Hamiltonian Htilde:

```text
T0 = R
T1 = Htilde R
T_{m+1} = 2 Htilde T_m - T_{m-1}
```

where R may contain B independent starting vectors as columns.

The block implementation must be mathematically identical to running B scalar
recurrences separately.

No changed normalization or random-vector distribution.

## C. Block sparse multiply

Where beneficial, evaluate:

```text
H * [r1 r2 ... rB]
```

using the appropriate sparse-matrix/multiple-vector GPU primitive rather than
B independent host-dispatched SpMV calls.

Prefer existing cuSPARSE functionality.

Benchmark block widths, for example:

```text
B = 1, 2, 4, 8, 16
```

subject to memory.

Do not hard-code a universal B.

## D. Block operator application

If current velocity/spin/orbital-current operators are sparse and applied
independently to the same vector block, reuse their descriptors and evaluate
them as block operations where profitable.

Preserve operator ordering exactly.

## E. Moment contractions

For the current two-index Kubo moments:

```text
mu_nm
```

use matrix-matrix contractions/cuBLAS where the block representation exposes
them naturally.

Do not change the current definition, complex conjugation, transposition, spin
or orbital-current convention.

## F. Memory model

Before implementing large B, record a memory model for:

```text
recurrence states
stored left Chebyshev states
right recurrence states
operator-applied states
mu_nm
workspace
```

as functions of:

```text
N
M
B
precision
```

For M=500, do not allow block-vector storage to consume the device
unpredictably.

Implement a capacity-aware block size or explicit benchmark control.

Do not silently reduce M.

## G. `per_type` path

Do not force block stochastic machinery onto `per_type` if there is only one
independent projector for elemental Pt.

However, if several type/projector vectors are mathematically independent,
the same block engine may process them together.

Keep the simple B=1 path efficient.

## H. Random-vector path

For the actual stochastic trace estimator, benchmark increasing Ntrace and
block width.

Mandatory:
- identical random vectors CPU/GPU;
- deterministic seed recording;
- same normalization;
- same final averaging.

## I. Correctness

For B>1 compare the entire batched result with the concatenation/average of
B scalar reference runs.

Test:
- FP64;
- FP32;
- several M;
- at least one real Pt case.

## J. Performance

Report:
- time per trace vector;
- total moment time;
- traces/second;
- end-to-end transport speedup.

The central question is whether GPU efficiency improves as Ntrace grows.

## Checklist

- [ ] all estimator modes documented
- [ ] block recurrence mathematically verified
- [ ] block sparse multiply implemented/evaluated
- [ ] block operator application evaluated
- [ ] moment contraction uses appropriate batched/GEMM path where profitable
- [ ] memory model documented
- [ ] B=1 remains supported
- [ ] per_type semantics unchanged
- [ ] stochastic vectors/seeds identical CPU/GPU
- [ ] B=1/2/4/8 tested
- [ ] B=16 tested where memory permits
- [ ] FP64 block correctness passes
- [ ] FP32 block correctness passes
- [ ] speedup vs Ntrace reported
- [ ] no estimator physics changed

**Commit message:** `Batch stochastic KPM trace vectors on CUDA`

---

# KPM-G3 — Evidence-gated GPU Kubo–Bastin reconstruction

## Objective

Only if G0--G2 show that CPU post-processing materially limits the accelerated
transport calculation, move the expensive Kubo–Bastin reconstruction/
contraction to the GPU.

This task is conditional.

A valid outcome is: **no implementation because profiling does not justify it.**

## Entry criterion

Proceed with production code only if the stage profile shows a material
fraction in one or more of:

```text
Gamma_nm(E)
Gamma_nm(E) * mu_nm contraction
energy integration
host transfer of full mu_nm
```

Use a documented threshold such as >15--20% of transport wall time in at least
one realistic M=500+ workload, or otherwise justify the threshold.

## A. Keep mu_nm resident

If reconstruction is GPU-accelerated, avoid:

```text
GPU mu_nm
 -> D2H full mu_nm
 -> CPU reconstruction
```

when the full moment tensor is not otherwise required on host.

Provide an explicit path to download moments for diagnostics/output.

Do not remove canonical host output capability.

## B. Preserve the production formula exactly

Port the current production Gamma/Kubo-Bastin algebra verbatim at the
mathematical level.

Do not:
- change kernel damping;
- change energy quadrature;
- change prefactors;
- change complex convention;
- change spin/orbital-current normalization.

## C. Numerical strategy

Implement both FP64 and FP32/mixed variants only if justified by G1.

For mixed precision, consider FP32 moments with FP64 accumulation/reconstruction
only if it has a clear numerical/performance rationale and is explicitly
labelled.

Do not silently promote FP32 inputs to "FP64 transport".

## D. Algorithmic implementation

Prefer:
- GPU matrix contractions;
- batched energy blocks;
- reusable Gamma/work buffers.

Avoid one CUDA kernel launch for every `(E,n,m)` scalar.

Use tiling over energy and/or moment indices if full Gamma storage is excessive.

## E. Correctness

Compare GPU reconstruction with the existing CPU reconstruction for identical
moments.

This isolates reconstruction error from stochastic/recurrence error.

Require:
- full energy-resolved conductivity comparison;
- tensor component comparison;
- integrated values;
- charge/spin/orbital modes currently validated.

## F. Performance

Report:
- D2H bytes avoided;
- reconstruction time CPU/GPU;
- total transport improvement.

Do not retain the GPU reconstruction if it increases total wall time.

## Checklist

- [ ] G0--G2 profile satisfies entry criterion OR task closed as unjustified
- [ ] CPU postprocessing bottleneck quantified
- [ ] resident moment handoff designed narrowly
- [ ] full moment download remains available
- [ ] production formula preserved
- [ ] FP64 reconstruction correctness passes
- [ ] FP32/mixed reconstruction labelled correctly
- [ ] energy-resolved outputs agree
- [ ] tensor outputs agree
- [ ] total speedup measured
- [ ] implementation retained only if beneficial

**Commit message if implemented:** `Accelerate Kubo-Bastin reconstruction on CUDA`

**Commit message if evidence says no code change:** `Document Kubo-Bastin reconstruction profile`

---

# KPM-B0 — Build a fair, reusable CPU/GPU benchmark harness

## Objective

Create the benchmark harness that will generate the numbers used in papers,
presentations, documentation, and discussions with collaborators.

This harness must make unfair comparisons difficult.

Do not optimize kernels in this task.

## A. Benchmark dimensions

Every row must uniquely specify:

```text
material
physical input/reference state
replication/supercell
Natom
Norb / Hamiltonian dimension
nnz
nsp
SOC state
cond_type
cond_calctype
Ntrace/projector count
cond_ll
lld
precision
backend
OMP threads
GPU strategy
random seed
```

## B. CPU configurations

For every CPU algorithm/precision combination run:

```text
OMP_NUM_THREADS=1
OMP_NUM_THREADS=2
OMP_NUM_THREADS=4
OMP_NUM_THREADS=8
```

Control and record:

```text
OMP_PROC_BIND
OMP_PLACES
MKL_NUM_THREADS / BLAS threads
compiler
optimization flags
CPU frequency/governor if available
NUMA placement
```

Use one MPI rank unless the benchmark explicitly studies MPI.

Avoid oversubscription.

For the headline speedup select the fastest reasonable CPU configuration for
the same precision.

Retain all thread-scaling rows.

## C. GPU configurations

For GPU record:

```text
GPU model
VRAM
driver
CUDA
compute capability
selected device
precision
block stochastic width
resident/nonresident moments
```

Use one GPU for the primary campaign.

Do not mix multi-GPU scaling into the headline.

## D. Warm-up and repetitions

For each configuration:

1. create/init backend;
2. warm up at least twice where practical;
3. reset timing counters;
4. run repeated measured transport work;
5. report median;
6. report min/max or IQR/MAD.

Use at least 5 repetitions for moderate workloads.

For very expensive realistic rows, at least 3 measured repetitions is
acceptable if documented consistently.

Do not spawn a fresh CUDA process for every measured kernel sample if a
persistent benchmark route exists.

## E. Precision-matched comparisons

Generate explicit paired rows:

```text
CPU_FP64 vs GPU_FP64
CPU_FP32 vs GPU_FP32
```

The harness should reject a request to compute a headline equal-precision
speedup from mismatched precision rows.

Mixed-precision rows get a separate label:

```text
MIXED
```

## F. Same random states

For stochastic trace:
- pre-generate/store or deterministically regenerate vectors;
- CPU and GPU consume the same vectors;
- record seed;
- record Ntrace.

The benchmark should fail if paired rows used different stochastic seeds.

## G. Stage timings

Collect:

```text
operator/setup
moment kernel
H2D
D2H
Gamma
Gamma*mu
integration
transport total
whole process
```

If G3 moves reconstruction to GPU, preserve compatible stage names.

## H. Correctness attached to every performance row

For every CPU/GPU paired benchmark, record at least:

```text
max/rms mu_nm difference where available
max/rms conductivity-spectrum difference
selected integrated/tensor differences
```

The performance table must indicate PASS/FAIL against the established
same-precision numerical tolerance.

No speedup claim for a failed numerical row.

## I. Output formats

Write:
- CSV;
- JSON with full metadata;
- Markdown summary.

Keep raw logs.

Include git commit hash and dirty-state check.

## J. Speedup definitions

Define explicitly:

```text
S_moment = T_bestCPU_moment / T_GPU_moment

S_transport = T_bestCPU_transport / T_GPU_transport

S_total = T_bestCPU_whole / T_GPU_whole
```

Never mix these.

## Checklist

- [ ] all benchmark dimensions recorded
- [ ] OMP 1/2/4/8 sweep automated
- [ ] BLAS/OpenMP oversubscription controlled
- [ ] same-precision pairing enforced
- [ ] mixed precision separately labelled
- [ ] stochastic seeds paired
- [ ] warmups automated
- [ ] repeated medians/spread reported
- [ ] stage timing captured
- [ ] correctness result attached to each paired row
- [ ] CPU-best comparator selected transparently
- [ ] CSV generated
- [ ] JSON generated
- [ ] Markdown generated
- [ ] git/hardware/build metadata generated

**Commit message:** `Add fair KPM CPU GPU benchmark harness`

---

# KPM-B1 — Run the systematic real-material benchmark campaign and publish the report

## Objective

Generate a performance dataset credible enough to show collaborators exactly
where GPU acceleration matters.

Do not modify implementation while collecting the final campaign.

Freeze the code first.

## A. Freeze environment

Record:

```text
git commit
dirty state
compiler/version
CUDA/version
driver
GPU
CPU
RAM
MPI
OpenMP settings
BLAS/MKL
build type/flags
```

Run correctness gates before timing.

## B. Mandatory Pt production anchor

At minimum:

```fortran
calctype = 'B'
nsp = 2
cond_ll = 500
lld = 150
recur = 'chebyshev'
cond_type = 'spin'
cond_calctype = 'per_type'
```

Use the real validated SOC Pt state.

Run several real-space sizes.

## C. Moment-order scaling

For at least one medium and one large real Pt system:

```text
M = cond_ll = 250, 500, 750
```

Add `1000` where practical.

Keep the physical system and all other parameters fixed.

This isolates M scaling.

## D. Real-space size scaling

At fixed realistic:

```text
M = 500
lld = 150
```

run several increasing real Pt supercell/replication sizes.

Report:
- Natom;
- Hamiltonian dimension;
- nnz.

This is the primary system-size scaling plot/table.

## E. Estimator scaling

For `per_type`, report the production case.

For the supported stochastic random-vector estimator:

```text
Ntrace = 1, 2, 4, 8, 16
```

and larger only where scientifically meaningful.

Use fixed seeds/vectors across backends.

Report:
- total time;
- time/trace;
- traces/s;
- statistical uncertainty where relevant.

## F. Conductivity modes

Mandatory:
- spin.

Where supported by the same physical fixture:
- charge;
- orbital.

Do not treat separate conductivity modes as different hardware systems; use
the same state where possible.

## G. Magnetic real-material campaign

Run at least one validated magnetic metallic system.

Use realistic M, not a tiny smoke-test order.

This checks whether spin structure/operator complexity changes the performance
picture.

## H. Precision pairs

For every headline workload run:

```text
CPU FP64:
    OMP 1,2,4,8

GPU FP64

CPU FP32:
    OMP 1,2,4,8

GPU FP32
```

Where a route is scientifically unsupported, report `unsupported`.

Do not silently substitute another precision.

## I. Headline tables

Produce separate headline tables for:

### 1. Equal-precision FP64

```text
system | N | nnz | M | estimator | best CPU threads |
CPU FP64 moment | GPU FP64 moment | S_moment |
CPU FP64 transport | GPU FP64 transport | S_transport |
correctness
```

### 2. Equal-precision FP32

Same schema.

### 3. Mixed/current-production route

Separate table.

## J. Required plots

Generate at least:

1. **moment-kernel speedup vs system size**;
2. **transport-total speedup vs system size**;
3. **speedup vs cond_ll**;
4. **CPU thread scaling for representative systems**;
5. **stochastic trace throughput vs Ntrace** if supported;
6. stage-fraction stacked data/table showing where time remains.

Plots must show data points and identify precision.

Do not combine FP32 and FP64 into one unlabeled speedup curve.

## K. Interpretation

The report must explicitly answer:

1. At what system size does GPU cross over?
2. At what M does GPU cross over?
3. Is the stochastic moment kernel accelerated much more than end-to-end
   transport?
4. What limits end-to-end speedup after G1/G2?
5. Does block stochastic processing increase GPU advantage?
6. Does CPU thread scaling materially change the crossover?
7. Is FP32 scientifically acceptable for the scoped transport workflows?
8. Is FP64 GPU still worthwhile?
9. Is GPU reconstruction G3 justified/beneficial?
10. Which workload should users actually select GPU for?

## L. No cherry-picking

Publish unfavorable rows.

If:
- small systems are CPU-faster;
- OMP=8 beats GPU for a regime;
- FP64 GPU has weak crossover;
- FP32 fails a scientific tolerance;

state this explicitly.

## Checklist

- [ ] environment frozen
- [ ] correctness gates passed
- [ ] realistic Pt M=500 anchor complete
- [ ] Pt size scaling complete
- [ ] M=250/500/750 scaling complete
- [ ] M=1000 measured or explicitly skipped
- [ ] per_type campaign complete
- [ ] stochastic-trace scaling complete if supported
- [ ] spin conductivity complete
- [ ] charge conductivity complete where supported
- [ ] orbital conductivity complete where supported
- [ ] magnetic real-material case complete
- [ ] CPU FP64 OMP 1/2/4/8 complete
- [ ] CPU FP32 OMP 1/2/4/8 complete
- [ ] GPU FP64 complete
- [ ] GPU FP32 complete
- [ ] equal-precision FP64 table published
- [ ] equal-precision FP32 table published
- [ ] mixed-precision table separately published
- [ ] kernel and end-to-end speedups distinguished
- [ ] plots generated
- [ ] unfavorable results retained
- [ ] user-facing backend recommendation written

**Commit message:** `Publish realistic KPM GPU performance campaign`

---

# 5. Recommended execution order

Run the work packages in this order:

```text
KPM-G0
    realistic baseline + profiler
        |
        v
KPM-G1
    precision parity + persistent GPU workspaces
        |
        v
KPM-G2
    block stochastic trace/vector acceleration
        |
        v
KPM-G3
    ONLY if CPU reconstruction is now material
        |
        v
KPM-B0
    fair reusable benchmark harness
        |
        v
KPM-B1
    final real-material benchmark campaign/report
```

If G0 already shows that CPU Gamma/reconstruction dominates, do not skip G1/G2:
first establish whether the GPU moment path itself is efficient and persistent,
then reconsider the post-processing boundary.

If G2 is irrelevant to the user's principal `per_type` workload, it may still
be completed for the true stochastic estimator, but its performance claim must
be scoped accordingly.

---

# 6. Acceptance criteria for the overall follow-up

This follow-up is successful if it produces all of the following, even if the
final GPU speedup is modest in some regimes:

1. a realistic `cond_ll=500`, `lld=150` production transport benchmark;
2. real-system size scaling;
3. equal-precision FP64 CPU/GPU performance;
4. equal-precision FP32 CPU/GPU performance;
5. fair CPU OpenMP scaling through 1/2/4/8 threads;
6. stage-level explanation of total runtime;
7. quantified stochastic-moment kernel speedup;
8. quantified end-to-end transport speedup;
9. a measured crossover map;
10. a clear user recommendation for CPU/GPU selection.

The campaign must make it impossible to confuse:

```text
fast GPU kernel
```

with

```text
fast complete transport calculation.
```

That distinction is the main scientific/performance question left open by ACC-13.
