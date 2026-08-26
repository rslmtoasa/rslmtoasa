# SCF-B1R2 — lean desktop performance closeout

## 1. Build and provenance

The completed artifact is the lean tier at
`results/benchmarks/scf_b1r2/lean/`. Its manifest is `rslmto.scf-b1r2.manifest.v2`
and records **18/18 completed cases**, all with row status `PASS`.

The Release executable was built in
`build-b1r2-release-cuda` with `/usr/bin/gfortran`, effective `-O3`, CUDA
enabled, and tracked source commit
`5967fbba82d9ccaf1de59614230c28c22de6284c`. Tracked source was clean when the
measurements were taken; unrelated untracked build/result artifacts were
present. No expensive runs were repeated during this closeout.

The lean tier intentionally uses one measured process per case, no warmup,
`nstep=1`, and is not a replacement for the full HPC tier. The canonical
records are [campaign.json](../../results/benchmarks/scf_b1r2/lean/campaign.json),
[campaign.md](../../results/benchmarks/scf_b1r2/lean/campaign.md), and the
[manifest](../../results/benchmarks/scf_b1r2/lean/manifest.json).

## 2. RS block recursion context

The preserved B1R Fe3 comparison remains the reference context: the best
8-thread FP64 CPU route measured about 69.7 s kernel / 78.3 s per iteration,
while the mixed CUDA route measured about 73.4 s kernel / 78.0 s per iteration.
That is approximately parity at this measured size, not an equal-precision GPU
speedup.

## 3. RS Chebyshev

The lean production rows use fast Chebyshev with `M=200`, Gaussian smearing
`sigma=0.01`, the production `-1.5..1.0` energy window and 2000-point DOS
grid. The CPU rows are FP64; the CUDA rows use an FP32 working recurrence with
FP64 canonical state and are therefore labelled `numeric_mode=mixed`.

| case | best CPU | GPU | CPU/GPU iteration ratio | kernel ratio | phase ratio | result |
|---|---:|---:|---:|---:|---:|---|
| Si1, `nmat=16` | 1.286 s (OMP8) | 0.709 s | 1.81x | 2.14x | 1.87x | GPU-beneficial, mixed |
| Si2, `nmat=128` | 88.042 s (OMP8) | 5.942 s | 14.82x | 19.54x | 15.74x | GPU-beneficial, mixed |

Both ratios are eligible production-route comparisons because the paired
common-state correctness gates pass. The measured trend is strongly more
favorable than the preserved Fe3 block-recursion parity result, although the
materials and size ladders are not identical enough to make a strict
cross-algorithm scaling claim.

The lean CPU sweep is Si1 OMP1/4/8 and Si2 OMP8; OMP2 and Si3/Si4 are reserved
for the full HPC tier. Thus this evidence answers the desktop landscape
question but does not establish the full size-scaling curve.

## 4. Reciprocal SCF iteration

The new reciprocal rows are full production SCF-iteration measurements from
the same common-state CPU/GPU pairings. The executable reported `nmat=144`
for Fe2 and `nmat=486` for Fe3, with `Nk_unique=65` in these fixtures.

| case | CPU FP64 | GPU FP64 | CPU/GPU iteration ratio | classification | correctness |
|---|---:|---:|---:|---|---|
| Fe2, `nmat=144` | 0.992 s | 1.168 s | 0.850x | CPU-preferred | PASS |
| Fe3, `nmat=486` | 7.333 s | 4.820 s | 1.521x | GPU-beneficial | PASS |

These are equal-precision FP64 comparisons. The ratios are promoted only after
the corresponding common-state checks passed. The Fe3 result shows the
application-level crossover in the measured lean workload; it is consistent
with, but independent of, the preserved K2 solver evidence.

## 5. Large reciprocal FP64 eigensolver

The current-build K2 frozen-potential evidence remains authoritative and was
not rerun: GPU/CPU speedups are approximately 1.49x at `nmat=486`, 3.24x at
`nmat=1152`, and 4.28x at `nmat=2250`. These are equal-precision FP64 solver
measurements with eigenvectors enabled, not full-SCF convergence claims.

## 6. CPU OpenMP scaling

The lean tier is intentionally not a complete OMP sweep. Its directly measured
Chebyshev CPU rows are:

| system | OMP1 | OMP4 | OMP8 |
|---|---:|---:|---:|
| Si1 iteration | 5.815 s | 1.855 s | 1.286 s |

Si2 has the selected desktop comparison at OMP8 (88.042 s). The full tier
adds OMP2 and the Si3/Si4 ladder with warmups and repeated measurements.

## 7. Precision and correctness qualifications

The Si1 and Si2 Chebyshev CPU/GPU common-state comparisons both pass the
observable and DOS gates. DOS uses the same energy grid and production
broadening; the relative L2 differences are `7.44e-7` and `5.07e-7`,
respectively. The reciprocal Fe2 and Fe3 common-state comparisons also pass.

Chebyshev GPU ratios are mixed-vs-FP64 production ratios, not equal-precision
speedups. Reciprocal Fe2/Fe3 ratios are equal-precision FP64. SCF cycle count
is retained as a correctness/stability diagnostic and is not used as the
performance metric.

## 8. Final workload map

| route | measured conclusion | precision caveat |
|---|---|---|
| RS block recursion, Fe3 | parity between mixed CUDA and best OMP8 CPU | mixed CUDA vs FP64 CPU |
| RS Chebyshev, Si1/Si2 | GPU-beneficial; ratio grows from 1.81x to 14.82x | mixed CUDA vs FP64 CPU |
| reciprocal Fe2 | CPU-preferred, ratio 0.850x | equal FP64 |
| reciprocal Fe3 | GPU-beneficial, ratio 1.521x | equal FP64 |
| large reciprocal K2 | GPU-beneficial from `nmat=486` through `2250` | equal FP64, frozen potential |

Within the measured lean desktop scope, Chebyshev is a credible production-SCF
GPU candidate and has the strongest observed full-iteration scaling. The full
HPC tier is still needed before claiming complete Si1–Si4 or OMP1/2/4/8
scaling.

## 9. K1 timer and deferred validation

The invalid CPU `P_eigensolver` component timer remains explicitly suppressed
as `solver_component_timer_status=INVALID_FOR_K1`. Complete reciprocal
electronic-structure/SCF iteration timing is authoritative for K1 application
performance, and K2 remains authoritative for dense-solver performance.

The RS-vs-k-space fixed-potential accuracy issue remains
`INCONCLUSIVE` and deferred. No formulation-level timing conclusion is drawn
from it.

## Closeout decision

The **lean desktop SCF performance scope can close**: all 18 planned cases are
present, correctness gates pass, eligible ratios are recorded, and plots were
regenerated under `results/benchmarks/scf_b1r2/lean/plots/`.

The **full HPC B1R2 scope does not close yet**, because it was deliberately not
run on the desktop. Run `bash tests/benchmarks/run_scf_b1r2_all.sh --tier full`
on an appropriate HPC system if complete Si1–Si4 scaling and repeated timing
statistics are required.
