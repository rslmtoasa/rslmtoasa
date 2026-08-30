# KPM-G1.4 report — optimized CUDA transport fixed-floor audit

Status: complete.  The optimized CUDA route has no obsolete full-host Gamma or
moment work in normal execution, the residual profile is decomposed, and the
remaining floor is explained well enough to stop G1.4 and proceed to G2.

Measurements were run on 2026-08-21 from HEAD `5a1ace1a1db51ba9ccb13b336ba1159aad190473`
using an NVIDIA RTX A4000 (driver 610.57.04, 16376 MiB).  The production fixture
is real SOC Pt, `M=500`, `NE=2510`, `lld=150`, spin, `per_type`, with `OMP_NUM_THREADS=1`
and `BLAS_NUM_THREADS=1`.  CUDA rows use FP32 moments and FP32 CUDA/cuBLAS
reconstruction unless stated otherwise.  Each r4/r6/r8 row has two warmups and
three fresh-process repetitions; values below are medians.

## Exclusive CUDA budget

The required exclusive phases all close exactly (`profile_closure_error=0` for
every row used below).  `P_stack_setup` and `P_moment_finalize` are additional
named phases split out of the former residual so that `P_misc` is not a generic
catch-all.

| Size | P_operator_setup | P_trace_setup | P_moments_total | P_gamma_basis_setup | P_gamma_generation | P_reconstruction_total | P_result_unpack | P_energy_integration | P_tensor_postprocess | P_output_prepare | P_output_io | P_stack_setup | P_moment_finalize | P_misc | T_transport_total | whole_wall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| r4 | 0.2946 | 0.0416 | 0.4435 | 0.5877 | 0.1029 | 0.0243 | 0.0009 | 0.0026 | 0.0007 | 0.0646 | 0.0046 | 0.5020 | 0.0720 | 0.0146 | 2.1580 | 2.3731 |
| r6 | 0.2972 | 0.1387 | 1.0524 | 0.5879 | 0.1029 | 0.0244 | 0.0009 | 0.0022 | 0.0006 | 0.0666 | 0.0052 | 0.5947 | 0.0714 | 0.0189 | 2.9713 | 3.1857 |
| r8 | 0.3117 | 0.3306 | 2.3217 | 0.5872 | 0.1031 | 0.0233 | 0.0009 | 0.0022 | 0.0007 | 0.0638 | 0.0052 | 0.9488 | 0.0701 | 0.0437 | 4.8112 | 5.0255 |

`P_misc` is 0.67%, 0.64%, and 0.91% of transport time for r4/r6/r8.  The
previously hidden roughly 3 s Simpson scan is now `P_energy_integration`; the
batched prefix implementation reduces it to about 2.6, 2.2, and 2.2 ms while
preserving the legacy Simpson domain and output bytes.

## Fixed-floor interpretation

The post-moment fixed work below is:

```text
P_gamma_basis_setup + P_gamma_generation + P_reconstruction_total
+ P_result_unpack + P_energy_integration + P_tensor_postprocess
+ P_output_prepare + P_output_io + P_moment_finalize + P_misc
```

`P_stack_setup` is shown separately because it is pre-moment host stack and
Hamiltonian setup, not post-moment transport work.

| Size | GPU moments | GPU fixed post-moment work | fraction moments of total | fraction fixed floor of total |
| --- | ---: | ---: | ---: | ---: |
| r4 | 0.4435 s | 0.8749 s | 20.5% | 40.5% |
| r6 | 1.0524 s | 0.8810 s | 35.4% | 29.7% |
| r8 | 2.3217 s | 0.9001 s | 48.3% | 18.7% |

The post-moment floor is therefore nearly independent of real-space size.  The
remaining non-moment, non-floor setup is explicit: `P_operator_setup +
P_trace_setup + P_stack_setup` is 0.84, 1.03, and 1.59 s.  It grows with the
real-space fixture because the production call constructs the lattice,
Hamiltonian, recursion stack, and trace state inside the measured transport
region; it is not obsolete Gamma or reconstruction work.

For context, the best CPU fast-FP64 rows from the G1.2 campaign were:

| Size | CPU fast FP64 moments | CPU fast FP64 total | GPU FP32 moments | GPU FP32 total |
| --- | ---: | ---: | ---: | ---: |
| r4 | 10.3861 s | 23.9821 s | 0.4435 s | 2.1580 s |
| r6 | 34.3592 s | 48.2105 s | 1.0524 s | 2.9713 s |
| r8 | 82.0803 s | 96.5290 s | 2.3217 s | 4.8112 s |

This is scaling context, not a precision-matched speedup claim: the CPU table
is FP64 and the headline CUDA table is FP32/FP32.

## Legacy-work and resident-memory audit

| Item on optimized CUDA route | Executed? | Required? | Evidence/action |
| --- | --- | --- | --- |
| Full host `gamma_nm(NE,M,M)` allocation | No | No | `calculate_gamma_nm` returns on resident CUDA moments; tiled Gamma basis/fill runs in the CUDA reconstruction backend. |
| Full host Gamma fill | No | No | No host Gamma allocation/fill in the normal resident path; CPU/reference path is retained. |
| CPU `calculate_gamma_nm` for an unused tensor | No | No | The CUDA resident check precedes the legacy CPU constructor. |
| CPU mu packing | No | No | Normal CUDA reconstruction consumes resident packed diagonal moments; CPU packing remains for CPU/reference routes. |
| CPU ZGEMM reconstruction | No | No | CUDA path uses the CUDA reconstruction backend and cuBLAS; CPU BLAS remains as reference functionality. |
| Full moment download | No | No | Normal rows report `full_moment_d2h_bytes=0`; the explicit diagnostic download path remains available. |
| Legacy scalar Gamma*mu loop | No | No | CUDA path receives reduced reconstruction results; legacy scalar/CPU consumers remain on reference routes. |
| Legacy full integrand reconstruction pass | No | No | Only reduced-result unpack and required tensor postprocessing remain. |
| Full host moment allocation/zero/fill | No | No | `host_full_moment_bytes=0` on all normal CUDA rows; the allocation is retained only for CPU or explicit diagnostic download. |

Normal FP32 memory counters were identical for r4/r6/r8:

| resident device moment bytes | host full-moment bytes | full-moment D2H bytes | reduced-result D2H bytes |
| ---: | ---: | ---: | ---: |
| 36,000,000 | 0 | 0 | 361,440 |

The explicit diagnostic mode can still request the full host moment tensor and
download; it is not used by the production or no-write campaigns.

## Output and no-write diagnostic

Output preparation is now separated from filesystem I/O.  Normal production
continues to write the existing files with the existing names, ordering, and
formats.  The benchmark-only `RSLMTO_KPM_BENCHMARK_NO_WRITE=1` mode executes
the same numerical work and formatting, but suppresses file opens/writes.

| Size | normal production wall | benchmark no-write wall | normal - no-write wall delta | exclusive `P_output_io` | direct I/O fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| r4 | 2.3731 s | 2.3695 s | +0.0036 s | 0.0046 s | 0.20% |
| r6 | 3.1857 s | 3.1873 s | -0.0016 s | 0.0052 s | 0.16% |
| r8 | 5.0255 s | 5.0542 s | -0.0287 s | 0.0052 s | 0.10% |

The signed whole-process deltas are below the run-to-run noise, especially at
r8; the exclusive phase is the authoritative attribution.  Normal output file
sizes were unchanged at r4/r6/r8: `cond_total.out` and `fort.123` 122,990 B,
and each orbital output 765,550 B.  The no-write scratch directories contained
none of those output files.  Thus the old 2.8--2.9 s “output” number was not
genuine filesystem I/O; it was the repeated energy integration that had been
timed inside the old output phase.

## Correctness and invariant checks

- The r4, r6, and r8 FP32 production outputs are byte-identical to the
  pre-G1.4 CUDA outputs for `fort.123`, total real/imaginary conductivity,
  orbital real/imaginary conductivity, and the Pt per-type files.
- CUDA r4 charge and orbital smoke rows passed profile validation with FP32
  moments/reconstruction and zero closure error.
- The CUDA r4 FP64 reference row passed with explicit FP64 moments and FP64
  reconstruction (`P_moments_total=8.6829 s`, `T_transport_total=11.1115 s`).
- The G1.2 CPU campaign already covered `OMP_NUM_THREADS=1,2,4,8` with BLAS
  held at one thread for the material CPU moment stage.  The newly exposed
  stack phase is setup/object construction rather than a material host
  numerical kernel, so no additional OpenMP optimization was justified.
- Energy scaling, Gamma basis, kernel weights, normalization, and conductivity
  conventions were not changed.  The Simpson batching reuses fixed weights
  within a call and reproduces the legacy `NPTS+8` endpoint domain exactly.
  No persistent cache was added: the current production process performs one
  transport call and the measured basis setup is already below the dominant
  moment stage at r6/r8.

## Validation and stop decision

The focused CUDA CTest set passed: `UnitKpmTransport`, `CudaPluginSurface`,
and `BenchmarkHarnessSchema` (3/3).  Python byte-compilation, the benchmark
harness contract test, and `git diff --check` also pass.  Both the CPU and
CUDA production targets build successfully.

The G1.4 questions resolve as follows:

1. The former 4--5 s CUDA floor is explained by a nearly constant 0.875--0.900
   s post-moment floor, explicit pre-moment stack/operator/trace setup, and
   GPU moments that scale with real-space size.
2. No obsolete CPU Gamma/reconstruction work executes on the optimized CUDA
   route.
3. No full host Gamma tensor is allocated or filled on that route.
4. Full moment D2H remains zero in normal optimized execution.
5. Genuine measured filesystem I/O is about 4.6--5.2 ms, not seconds.
6. Output preparation/formatting is about 64.6, 66.6, and 63.8 ms.
7. The largest remaining named post-moment numerical stage is the tiled Gamma
   basis setup at about 0.587--0.588 s; it is already a CUDA stage and is
   nearly size-independent.  The largest host setup stage is `P_stack_setup`,
   which is explicit, size-dependent construction work rather than a redundant
   postprocessing kernel.
8. No remaining single-trace/per-type numerical host stage justifies another
   pre-G2 optimization work package.
9. A factorized-Gamma task is not justified by this profile: Gamma generation
   and reconstruction are already small, resident, and size-independent here.
10. G1.4 is complete; proceed to G2.  No G2 implementation was added in this
    work package.
