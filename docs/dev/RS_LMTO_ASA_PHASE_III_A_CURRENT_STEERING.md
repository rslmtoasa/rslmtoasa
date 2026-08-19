# RS-LMTO-ASA Phase III-A — Current GPU Performance Steering
## Post ACC-P0 / ACC-P1 / ACC-P2 / ACC-P3 / ACC-P4 update — ACC-10 / ACC-11 / ACC-12 / ACC-13 complete; active entry point: ACC-14

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

4. **The post-P2 SCF-shaped campaign confirms that the useful regime survives
   eigenvector transfer.**  With persistent `vectors=yes` requests, CUDA is
   still clearly preferred for `nmat=1152` and `2250`; `nmat=486` remains a
   near-boundary case.  This is a measured RTX A4000 policy input, not a
   portable global default.

Therefore the current Phase III-A goal is no longer:

> "Can reciprocal CUDA work at all?"

It is now:

> **Determine the best precision/solver/transfer strategy as a function of real LMTO matrix size and workload, optimize the already-useful medium/large-matrix CUDA path, and then decide which downstream GPU work packages are justified.**

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
| **ACC-P4** | **Completed** | Production CPU/GPU/precision policy frozen from post-P2 `vectors=yes` evidence | ACC-10/11/12/13 closed; proceed to ACC-14 |
| **ACC-10** | **Completed on current HEAD** | Backend-owned CUDA Lehmann contraction over the existing host-eigenpair contract | Reuse; do not reopen unless a new residency or mixed-routing question is justified |
| **ACC-11** | **Completed on current HEAD** | Narrow backend-owned CUDA FP64 resident eigensystem handoff for Lehmann | Reuse only for the explicit normal-mesh CUDA Lehmann path; no general residency framework |
| **ACC-12** | **Completed — not justified** | Post-P2/P4 profiling leaves CPU H(k) assembly below the materiality gate; no GPU assembler added | Keep the canonical CPU Fourier assembler; reopen only with new evidence |
| **ACC-13** | **Completed on current HEAD** | Existing RS CUDA KPM moment path validated for production charge/spin/orbital transport; no new kernel justified by the measured audit | Reuse; keep the focused production probe as the hardware evidence |
| ACC-14 | Final gate | Accelerator support/performance matrix | Active next task |

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

These were inputs to ACC-P4, not reasons to reopen ACC-P2.

# 10. ACC-P3 — completed evidence corpus

ACC-P3 is complete.

Its real-material crossover corpus is an **evidence source**, not a pending work package.

After ACC-P2, only rerun the specific established P3 rows needed to quantify the effect of P2 on current HEAD. Do not recreate the full P3 development task.

Use the same real fixtures and benchmark conventions so pre-P2 and post-P2 results remain comparable.

ACC-P4 used the latest valid P3/P2 combination as its performance dataset.

---

# 11. ACC-P4 — completed production policy

ACC-P4 is complete.  It consolidated the existing P1/P1b/P2/P3 evidence and
ran the current post-P2 persistent real-material campaign with `vectors=yes`.
No new solver or hot-path implementation was needed for this decision.

## Post-P2 campaign used for the policy

The decisive run used the existing `accp2_real_material.py` driver with one
warm-up, three measured repetitions, `Nk=1`, tile size 1, FP64
`zheevd_serial`, and the production Si/Fe fixtures.  The CUDA host was an
RTX A4000 with CUDA 13.3.  Times below are persistent end-to-end reciprocal
request times and include host H(k) assembly, transfers, eigensolution, and
host synchronization.

| Fixture | nmat | CPU LAPACK | CUDA FP64 Zheevd | CPU/CUDA | P4 classification |
|---|---:|---:|---:|---:|---|
| diamond Si primitive, vectors=yes | 16 | 3.03e-5 s | 6.35e-4 s | 0.048x | CPU-preferred |
| bcc Fe 1^3, vectors=yes | 18 | 4.62e-5 s | 7.01e-4 s | 0.066x | CPU-preferred |
| bcc Fe 3^3, vectors=yes | 486 | 5.29e-2 s | 3.61e-2 s | 1.46x | boundary; retain CPU default |
| bcc Fe 4^3, vectors=yes | 1152 | 6.49e-1 s | 1.63e-1 s | 3.98x | FP64 CUDA-preferred |
| bcc Fe 5^3, vectors=yes | 2250 | 3.66e0 s | 8.58e-1 s | 4.27x | FP64 CUDA-preferred |

The Fe `3^3` row is a real SCF-shaped CUDA win, but it is close enough to
the promotion boundary that it is not a reason to introduce automatic
dispatch.  The earlier values-only corpus independently keeps `nmat=486`
CPU-preferred and confirms the same large-matrix direction: values-only
CUDA is preferred at `nmat=1152` and `2250`.

The measured scope is `Nk=1`, tile 1, persistent steady state.  It does not
claim a universal threshold for cold startup, other GPUs, or different
mesh/tile shapes.  Cold startup remains CPU-preferred for small matrices;
larger CUDA claims require a persistent workload.

## Frozen reciprocal policy

| Workload | Production policy |
|---|---|
| `nmat <= 144`, values or vectors | CPU LAPACK |
| `nmat=486`, values-only | CPU LAPACK |
| `nmat=486`, vectors=yes | CUDA is supported and competitive; CPU remains the conservative default |
| `nmat >= 1152`, FP64 standard-Hermitian, persistent | CUDA `zheevd_serial` is the preferred measured route on the campaign GPU |
| generalized/overlap eigensolution | CPU LAPACK; CUDA remains unsupported |
| FP32 reciprocal/SCF | Not established; experimental probes only, no production selection |

LAPACK remains the global default.  Backend selection remains explicit:
`reciprocal_backend='lapack'` is the safe default and
`reciprocal_backend='cuda'` requests CUDA without a CPU fallback.  A
conservative `auto` mode is **not justified yet**: the evidence is from one
GPU, one mesh/tile shape, and does not provide enough hardware- or
request-shape calibration to make a portable automatic threshold safe.

This preserves the no-silent-fallback rule: an explicit CUDA request that is
unsupported or has no usable device fails at the typed backend boundary; it
does not silently execute LAPACK.

## Downstream decisions

- **FP32:** P1b's completed local physical probe had no usable CUDA device.
  Its status is `unsupported`, not pass or fail; FP32 remains unestablished.
  Application scientific arrays stay FP64.
- **Pinned staging:** transfer components improve in isolated rows, but total
  time changes were `-3.7%`, `+1.2%`, and `-1.5%` for Fe `3^3`, `4^3`, and
  `5^3`, respectively.  No general pinned-memory promotion is made.
- **ACC-11:** the prior defer decision is now closed by focused combined
  evidence.  The implementation is deliberately narrow: CUDA FP64 normal-mesh
  Lehmann requests reuse the eigensystem still owned by the CUDA context, while
  the public host eigenpair contract remains materialized for compatibility.
  The measured benefit is removal of the redundant ACC-10 eigenpair H2D copy;
  the contraction still dominates total production time, so no broad end-to-end
  speedup claim is made.
- **ACC-12:** close GPU H(R)->H(k) assembly as not justified by the measured
  workload.  CPU H(k) is about `11.9%`, `10.0%`, and `6.3%` of the post-P2
  Fe `3^3`, `4^3`, and `5^3` requests while the GPU eigensolver dominates;
  the canonical CPU Fourier assembler remains appropriate.
- **ACC-10:** the existing host-eigenpair CUDA contraction is now validated and
  closed.  Its measured production route uses CUDA eigenpairs in the large-
  matrix regime; CPU-preferred reciprocal cases remain on the CPU route.
- **ACC-13:** is now closed for the existing RS CUDA KPM transport path; the
  focused production evidence is recorded in `ACC-13_KPM_TRANSPORT_CUDA.md`.

ACC-P4 checklist is therefore complete: policy frozen, FP32 not promoted,
LAPACK default retained, no automatic dispatch added, ACC-11 completed only in
its evidence-backed narrow scope, ACC-12 closed without a GPU assembler,
ACC-10 closed, ACC-13 completed for the existing RS KPM transport path, and
ACC-14 is the next implementation/documentation track.

# 12. Revised relevance of original ACC-10...ACC-14

This is the main steering section for the work that follows P4.

---

## ACC-10 — GPU Lehmann Green-function contractions

### Status
**Completed on current HEAD for the initial host-eigenpair contract.**

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
The initial implementation is now closed and uses the P4-compatible host boundary:
CPU Fourier assembly and the selected reciprocal eigensolver produce ordinary host
eigenpairs, one CUDA request evaluates all directed pair blocks for the normal and
eta contours, and the result is copied into the existing canonical Green arrays.
The established Pauli decomposition still owns the torque-resolved arrays.

The production selection remains explicit.  The P4 large-matrix route may use
CUDA-produced host eigenpairs, while CPU-preferred reciprocal cases remain on the
CPU route.  No automatic mixed CPU-eigensolver/CUDA-contraction dispatch or device
residency framework is added by ACC-10.

Use P4 to select whether ACC-10 starts from:
- CPU eigenpairs;
- GPU eigenpairs;
- both, depending on regime.

For the current implementation, the selected measured route is the CUDA
eigensolver plus host-eigenpair handoff in the large-matrix regime.  A separate
CPU-eigensolver-to-CUDA-contraction mode remains outside this focused completion;
it must be justified by a new workload measurement before expanding the API.

### Completion evidence on current HEAD

The focused CUDA build and tests completed successfully on the NVIDIA host:

```text
Acc10LehmannSource       PASS
UnitLehmannChain         PASS
UnitDysonEquivalence     PASS
ReciprocalCudaLehmann   PASS
```

The isolated complete pair/energy/block request reports:

```text
max_error=5.64785e-17
```

The current production bcc-Fe VAL-05 CUDA campaign passed all k-mesh,
broadening, energy-window, selected onsite/intersite Green, and Sigma=0
Lehmann/Dyson checks.  The three CUDA route triads also passed with the
established values:

| Consumer | CUDA Lehmann | Dyson(Sigma=0) | Recursion |
|---|---:|---:|---:|
| bcc-Fe Jij, pair 1-335 | 0.2547380660120 | 0.2547380660129 | 0.5078779010703 |
| bcc-Fe conductivity `sigma_xx` | 3.235054 | 3.235054 | 3.239956 |
| bcc-Fe damping `alpha` | 0.002527619 | 0.002527619 | 0.001341155 |

Full implementation scope and earlier end-to-end timing evidence remain in
[`ACC-10_LEHMANN_CUDA.md`](ACC-10_LEHMANN_CUDA.md).  ACC-11 is now complete in
the narrow scope recorded below; the host-eigenpair ACC-10 contract remains
available and unchanged for other routes.

---

## ACC-11 — device-resident eigensystem handoff

### Status
**Completed on current HEAD, narrowly for the explicit normal-mesh CUDA FP64
Lehmann route.**

### Decision and measured result

The ACC-10 CUDA contraction was worthwhile in the established large-matrix
route, and the combined eigensolver/Lehmann path now reuses the CUDA context's
latest FP64 eigensystem.  The isolated resident contraction passed with
`resident_error=3.14018e-16`.  The full CUDA VAL-05 campaign passed all seven
mesh/broadening/window cases, with every production timing line reporting
`resident=1`.

The representative bcc-Fe `k12`, `eta=0.02` request reports:

```text
resident=1 total=2.798 h2d=1.70432001e-4 contraction=2.79478687 d2h=6.54079974e-4
```

For the larger `k16` case, the old host-eigenpair route measured
`h2d=3.58076811e-3`, while the resident route measured
`h2d=2.25664e-4` (about a 16x reduction).  Total time remains contraction
dominated, so this is a focused transfer reduction rather than a claim that
the whole production request is materially faster.

### Implementation boundary

- The CUDA context owns the device eigensystem and exposes only a
  token/shape-validity query; no device pointer crosses into Fortran physics
  code.
- A generation token invalidates residency after a new solve, a values-only
  solve, FP32 solve, or backend release.  The resident request also checks
  matrix size, mesh size, batch shape, and token before launching.
- `fill_green_lehmann` requests the resident route explicitly.  For that route,
  the normal mesh is solved as one CUDA tile so the context owns the complete
  eigensystem consumed by the contraction.  Ordinary CUDA tiling is unchanged.
- Host eigenvalues/eigenvectors remain available and the canonical Green and
  Pauli arrays are populated through the existing interfaces.
- CPU execution, arbitrary-k execution, values-only execution, FP32, and
  generalized/overlap paths remain on their existing routes.  No generic
  multi-consumer residency abstraction or automatic CPU/GPU dispatch was added.

### Completion checks

- CUDA C++ resident contraction unit: pass; stale-token rejection: pass.
- CUDA source-contract tests and resident source-contract test: pass.
- Focused CUDA CTest selection: 3/3 pass.
- CPU `rslmto.x` rebuild after the interface changes: pass.
- CUDA VAL-05 production campaign: pass, with all seven cases reporting
  `resident=1` and the established direct-Green/Dyson convergence checks
  unchanged.

---

## ACC-12 — GPU H(R)->H(k) assembly

### Status
**Completed on current HEAD — no production-code change justified.**

### Decision and measured result

The post-P2 `vectors=yes` real-material campaign measured the CPU Fourier
assembly fraction of the persistent CUDA reciprocal requests as:

| Fixture | nmat | CPU H(k) assembly fraction | Decision |
|---|---:|---:|---|
| bcc Fe `3^3` | 486 | 11.9% | retain CPU assembly |
| bcc Fe `4^3` | 1152 | 10.0% | retain CPU assembly |
| bcc Fe `5^3` | 2250 | 6.3% | retain CPU assembly |

The CUDA eigensolver remains the dominant component in all three measured
requests.  Earlier ACC-07 evidence also established that the public H(k)
compatibility cache is filled from the already assembled host tile, so there
is no redundant device-to-host H(k) transfer for a GPU assembler to remove.
H(k) remains a real compatibility product for Dyson, BSF, and legacy/bands
consumers.

The quantitative action gate therefore closes ACC-12 without code: H(k) has
not become a material end-to-end bottleneck after solver optimization.  The
canonical CPU `reciprocal_assembler` is retained, no GPU Fourier kernel or
new residency/cache flag is added, and the established phase conventions and
physics are unchanged.  ACC-12 may be reopened only if a new validated
workload makes H(k) material or creates a repeated-residency consumer that
changes this cost balance.

### Completion checklist

- [x] post-P2/P4 H(k) assembly fraction measured
- [x] representative real Fe `3^3`, `4^3`, and `5^3` workloads measured
- [x] quantitative no-port decision recorded
- [x] no GPU assembler added when the evidence did not justify one
- [x] H(k) phase conventions and downstream compatibility retained

---

## ACC-13 — RS KPM/Kubo-Bastin GPU completeness and performance

### Status
**Completed on current HEAD for the existing CUDA moment path.**

### Completion decision

The audit found no missing high-value transport kernel.  The existing CUDA
plugin already covers the production dataflow:

```text
CPU Hamiltonian/current-operator construction
    -> CUDA stochastic Chebyshev moments
    -> existing CPU Gamma/Kubo-Bastin reconstruction and output
```

The same backend-owned stochastic moment entry point is used for charge, spin,
and orbital conductivity.  The separate orbital-moment route is also covered
by the existing CUDA orbital-moment entry point.  CPU postprocessing remains
the canonical consumer; no transport formula or operator convention was
changed.

### Current evidence

The low-level CUDA ABI validator passed all 15 existing recursion,
stochastic/orbital-moment, and DOS/GF comparisons at a maximum reported
relative error of approximately `2.2e-15` in the FP64 validation mode.

The focused production probe
[`acc13_kpm_cuda.py`](../../tests/validation/acc13_kpm_cuda.py) then ran the
real SOC fcc-Pt conductivity fixture through the same CUDA-enabled binary in
CPU and CUDA modes.  For `cond_ll=20`, replication 4, and `rc=20`, charge,
spin, and orbital conductivity all passed the `5e-3` complex-observable
envelope; the measured CPU/GPU wall ratios were `0.61x`, `0.66x`, and `0.62x`
respectively, so this small workload remains CPU-preferred.

On the larger `cond_ll=40`, replication-6 probe, the same three observables
passed with complex relative errors below `1.0e-6`; CPU/GPU wall ratios were
`1.26x`, `1.33x`, and `1.34x` for charge, spin, and orbital respectively.
This establishes a useful measured regime without claiming a portable
automatic threshold or a universal transport speedup.

Full scope, limitations, and reproduction commands are recorded in
[`ACC-13_KPM_TRANSPORT_CUDA.md`](ACC-13_KPM_TRANSPORT_CUDA.md).

### Post-run verification on the current CUDA build

On 2026-08-19, the complete requested validation was rerun after restoring
the explicit `execution_backend%synchronize()` boundary required by the
existing arbitrary-k source contract.  The low-level validator again passed
all 15 CUDA ABI routes.  The CUDA-enabled unit/source gate passed `58/58`,
including the CUDA lifecycle, device-mapping, reciprocal, KPM transport, and
ACC source-contract tests.

The full VAL-09 campaign passed all `59/59` CPU cases and all `59/59` CUDA
cases using `--gpu-plugin --gpu-backend csr`.  The campaigns covered the
charge, spin, and orbital tensor/sign checks, Fermi/order/replication/window
sweeps, and recursion-versus-Lehmann consistency.  The configured GPU matrix
also passed all `10/10` CPU backend regression cases and `8/8` CPU-vs-GPU
consistency checks; its two MKL-specific cases were explicitly skipped because
`ENABLE_MKL_KERNELS=OFF`.

### ACC-13 closeout rerun on current hardware

The closeout was independently rerun on 2026-08-19 with the RTX A4000.  The
low-level validator again passed all 15 CUDA ABI routes, with a maximum
displayed FP64 relative error of `2.20e-15`.  The focused production probe
passed charge, spin, and orbital conductivity at both workload sizes:

| Probe | Charge CPU/GPU | Spin CPU/GPU | Orbital CPU/GPU | max complex error |
|---|---:|---:|---:|---:|
| `cond_ll=20`, replication 4 | 0.622x | 0.655x | 0.631x | `7.3e-7` |
| `cond_ll=40`, replication 6 | 1.382x | 1.341x | 1.348x | `1.0e-6` |

These are single current-host process-wall samples and are consistent with
the existing policy: the small case remains CPU-preferred, while the larger
case is CUDA-beneficial on this RTX A4000.  The full current CUDA VAL-09
campaign also passed all 59 cases, including charge/spin/orbital tensor and
sign checks, convergence sweeps, and the recursion-versus-Lehmann check.

### Retained boundary and negative results

- Existing CPU Hamiltonian/current-operator construction and CPU
  Gamma/Kubo-Bastin postprocessing are retained.
- No redundant charge/spin/orbital-specific CUDA kernel was added: the
  shared stochastic moment contract is the correct reuse point.
- The default CUDA Chebyshev transport arithmetic remains FP32 for moments
  with FP64 host outputs; the low-level FP64 run is a reference check, not a
  production precision-policy promotion.
- Small production transport cases are CPU-preferred; larger cases can be
  CUDA-beneficial on the measured RTX A4000.  No automatic dispatch was added.
- Structured FFT/conv orbital moments, `ccor_2c+hoh`, and local-axis routes
  remain explicitly unsupported as recorded by the existing guards.

### Completion checklist

- [x] transport CPU/GPU dataflow mapped
- [x] existing GPU moment support validated
- [x] charge conductivity CPU/GPU compared
- [x] spin conductivity CPU/GPU compared
- [x] orbital conductivity CPU/GPU compared
- [x] symmetry/operator conventions retained through the existing VAL-09 path
- [x] performance hotspots measured on small and larger real workloads
- [x] no redundant kernel added
- [x] end-to-end transport speedup recorded, including the CPU-preferred result
- [x] support matrix input prepared for ACC-14

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
ACC-P4  <-- COMPLETE
    |
    +-----------------------------+
    |                             |
    v                             v
ACC-10  <-- COMPLETE           ACC-13  <-- COMPLETE
host-eigenpair CUDA             RS KPM/transport
    |
    +--> ACC-11  <-- COMPLETE (narrow CUDA resident Lehmann handoff)

ACC-12  <-- COMPLETE (CPU H(k) assembly retained; GPU port not justified)

    |
    v
ACC-14 final accelerator support/performance gate  <-- NEXT
```

ACC-P1b, ACC-P2, and ACC-P3 are **completed evidence inputs to P4**, not preliminary gates that the next assistant should execute again.

P4 recorded FP32 as unestablished and proceeded with the FP64 policy. It does not block the rest of Phase III-A.

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

# 17. CURRENT HANDOFF — instruction for the next assistant

**Do not start at the top of this document.**

The completed history is:

```text
ACC-00...09  completed
ACC-P0       completed
ACC-P1       completed
ACC-P2       completed
ACC-P3       completed
ACC-P4       completed
ACC-10       completed
ACC-11       completed (narrow CUDA resident Lehmann handoff)
ACC-12       completed (CPU H(k) assembly retained; GPU port not justified)
ACC-13       completed (existing CUDA KPM transport path)
```

ACC-P1b was an evidence input.  Its physical CUDA rows were unsupported on
the completed local probe, so P4 left FP32 unestablished.

## The next active work package is ACC-14.

ACC-10 is closed with the existing host-eigenpair contract, ACC-11 is complete
in its evidence-backed narrow CUDA scope, ACC-12 is closed without a GPU
assembler, and ACC-13 is complete for the existing RS KPM transport path.
ACC-14 should proceed next.  Do not reopen ACC-P1b/P2/P3 or expand ACC-11/12
beyond the evidence gates stated above.

Do not repeat:
- benchmark-harness construction;
- supercell construction;
- true-batching implementation;
- P2 hot-path optimization;
- P3 crossover campaign.

Only rerun established benchmark rows when needed to obtain current post-P2 numbers on the same corpus.

The remaining sequence is:

```text
ACC-10 is complete; ACC-11 is complete in its narrow resident-Lehmann scope;
ACC-12 is complete with the measured no-port decision;
ACC-13 is complete for the existing CUDA KPM transport path.
ACC-14 is the next implementation/documentation track.
ACC-14 closes the accelerator campaign.
```

The current steering principle is:

> **The benchmarking/solver-rescue phase is complete. The measured execution policy is frozen; spend engineering effort on the downstream workloads that can deliver additional real application speedup.**
