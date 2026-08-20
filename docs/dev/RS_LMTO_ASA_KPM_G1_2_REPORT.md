# KPM-G1.2 benchmark methodology

## Status

The KPM transport profile now has an exclusive phase contract, nested
implementation timers, precision/backend metadata, controlled thread
settings, parser validation, and a Pt campaign driver. The complete CPU
r4/r6/r8 campaign was run on 2026-08-20; the GPU portion was not run because
this host has no active NVIDIA driver. The measured CPU results and that
limitation are recorded below.

## G1.1 timer audit

The old fields were useful for locating work, but they were not an additive
partition of transport time:

| Old field | Actual scope after the audit | Why it was not additive |
| --- | --- | --- |
| `T_operator` | one-time operator construction/upload | setup work was inside the outer transport timer |
| `T_trace_setup` | per-trace random vectors and trace initialization | repeated inside `T_transport_total` |
| `T_cheb_moments` | CPU wall time, or a GPU event-derived moment time | GPU H2D/kernel/D2H values were mixed into related fields and CPU/GPU scopes differed |
| `T_H2D` | GPU host-to-device transfer detail | nested in moment work; not a phase |
| `T_D2H` | GPU device-to-host transfer detail | nested in moment work; not a phase |
| `T_gamma` | Gamma basis construction and fill | nested in total transport work |
| `T_mu_pack` | diagonal mu packing | nested in reconstruction |
| `T_gamma_mu` | reconstruction contraction | nested in reconstruction and total transport |
| `T_energy_integral` | numerical integration plus output-file work | I/O was included in an alleged integration timer |
| `T_transport_total` | inclusive end-to-end KPM transport scope | it enclosed every other old timer, so summing the fields double-counted work |

The old names remain as compatibility aliases in `KPM_PROFILE` records, but
new analysis uses only the exclusive `P_*` fields and the nested `D_*` fields.

## Exclusive phase contract

Every transport record partitions the inclusive total as

```text
T_transport_total = P_operator + P_trace_setup + P_moments_total
                    + P_gamma + P_reconstruction_total
                    + P_energy_integration + P_output_io + P_other
```

`P_other` is an explicit residual for work not assigned to a named phase. The
emitter reports `profile_closure_error` and `profile_child_error`, and sets
`PROFILE_STATUS=FAIL` if the residual is materially negative, the exclusive
sum does not close, or a detail timer exceeds its parent. The Python parser
applies the same checks and rejects invalid rows in the G1.2 campaign driver.

The detail timers are:

| Parent phase | Detail timers |
| --- | --- |
| `P_moments_total` | `D_moment_H2D`, `D_moment_GPU_kernel`, `D_moment_D2H`, `D_conversion` |
| `P_reconstruction_total` | `D_mu_pack`, `D_reconstruction_BLAS` |
| `P_gamma` | `D_gamma_basis`, `D_gamma_fill` |

CPU and GPU moments use the same parent boundary: the trace setup has finished
and moment generation has not yet returned to the transport driver. For GPU,
that parent includes H2D, the CUDA kernel, D2H, and host precision conversion.
The CPU routes use the same semantic boundary even though their detail timers
are zero.

Energy integration stops before any file is opened. `P_output_io` covers the
total and per-type conductivity output writes. Gamma and packed-mu byte
counts are emitted separately as `bytes_gamma` and `bytes_mu_pack`; transfer
bytes remain `bytes_h2d` and `bytes_d2h`.

## Precision and environment metadata

Each row records:

- `moment_backend`: `cpu_legacy`, `cpu_fast`, or `cuda`;
- `moment_precision`: `fp32` or `fp64`;
- `reconstruction_backend`: currently `cpu_blas`;
- `reconstruction_precision`: currently `fp64` in the production route;
- `OMP_NUM_THREADS`, controlled BLAS thread count, MPI rank context, compiler,
  BLAS/LAPACK, CUDA, GPU, CPU, `OMP_PROC_BIND`, and `OMP_PLACES` where the
  harness can observe them.

The CPU fast route is therefore recorded honestly as FP32 moments with FP64
reconstruction. A minimal `CGEMM` analogue is available as
`gamma_mu_cblas` for a future validated FP32 reconstruction route; the current
production reconstruction remains FP64 and is not labelled a full FP32 run.

The campaign driver controls `OMP_NUM_THREADS`, `BLAS_NUM_THREADS`,
`MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS`, and defaults to an OpenMP sweep
of 1, 2, 4, and 8 with BLAS held at 1. Each repetition launches a fresh
process, with two warm-ups and five recorded repetitions by default.

## Campaign definition

The reproducible driver is:

```bash
python3 tests/benchmarks/kpm_g12_transport.py \
  --binary build-ci-reference-serial/bin/rslmto.x \
  --build-dir build-ci-reference-serial \
  --scratch-root /tmp/rslmto-kpm-g12 \
  --output results/benchmarks/kpm_g12_pt.json \
  --warmups 2 --repetitions 5
```

Its default scope is real Pt SOC, r4/r6/r8, `M=500`, `NE=2510`, `lld=150`,
per-type spin transport, and the legacy/fast/fast-double CPU routes. Charge
and orbital r4 rows can be selected with `--replications 4 --cond-types
charge orbital`. GPU fp32/fp64 rows are opt-in with `--gpu`; the driver only
computes speedups inside equal moment/reconstruction precision groups.

For publication-quality claims, the paired rows must additionally pass the
existing moment and conductivity correctness comparisons at the same Pt
dimensions and estimator. A profile closure pass is a timing-integrity gate,
not a physics-correctness result.

## Initial post-change check

The serial fast-double Pt r4 smoke run (`M=500`, `NE=2510`, `lld=150`, one
OpenMP/BLAS thread, four repetitions) emitted `PROFILE_STATUS=PASS` with zero
reported closure error. Its representative profile showed:

| Field | Seconds |
| --- | ---: |
| `P_moments_total` | 10.1805 |
| `P_gamma` | 8.2344 |
| `P_reconstruction_total` | 1.9708 |
| `P_energy_integration` | 0.00125 |
| `P_output_io` | 2.8243 |
| `P_other` | 0.5757 |
| `T_transport_total` | 23.8301 |

This check demonstrates the methodology change: numerical integration and
output I/O are no longer presented as one timer. It is a local smoke result,
superseded as a benchmark result by the full campaign below.

## Full CPU campaign results

The default campaign was run with the committed driver and serial reference
binary:

```bash
python3 tests/benchmarks/kpm_g12_transport.py \
  --binary build-ci-reference-serial/bin/rslmto.x \
  --build-dir build-ci-reference-serial \
  --scratch-root /tmp/rslmto-kpm-g12-full \
  --output results/benchmarks/kpm_g12_pt_full.json \
  --warmups 2 --repetitions 5
```

The run produced 36 valid rows and 180 measured samples: 12 rows each for
r4, r6, and r8, covering legacy FP64 moments, fast FP32 moments, fast FP64
moments, and OpenMP settings 1/2/4/8 with BLAS held at one thread. All 180
profiles passed the exclusive-phase closure and nested-child checks. The
workloads were real Pt SOC, `M=500`, `NE=2510`, `lld=150`, spin transport,
and the `per_type` estimator.

The table reports the best transport median over the four OpenMP settings for
each route; `wall` is the corresponding whole-process median. `fast FP32` is
mixed precision because reconstruction remains FP64.

| Size | N | nnz | Legacy FP64 transport / wall (s) | Fast FP32 transport / wall (s) | Fast FP64 transport / wall (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| r4 | 1,152 | 215,784 | 34.082665 / 34.399029 | 23.930090 / 24.232400 | 23.982069 / 24.286550 |
| r6 | 3,888 | 884,520 | 82.527604 / 82.854602 | 48.250907 / 48.573595 | 48.210506 / 48.529115 |
| r8 | 9,216 | 2,299,752 | 177.596385 / 177.925224 | 96.586468 / 96.943705 | 96.529018 / 96.860246 |

For the best fast FP64 rows, the exclusive phase medians were:

| Size | OMP | P_moments_total | P_gamma | P_reconstruction_total | P_energy_integration | P_output_io | P_other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| r4 | 4 | 10.386106 | 8.197732 | 1.938904 | 0.001352 | 2.870264 | 0.584480 |
| r6 | 1 | 34.359155 | 8.202176 | 1.944229 | 0.001330 | 2.852123 | 0.715698 |
| r8 | 2 | 82.080269 | 8.210503 | 1.943107 | 0.001305 | 2.858237 | 1.110225 |

The post-G1.1 bottleneck is therefore still Chebyshev moment generation on
the optimized CPU route. It grows from about 10.4 s at r4 to 82.1 s at r8,
while Gamma generation stays near 8.2 s and BLAS reconstruction near 1.94 s.
Integration is negligible at these dimensions; output I/O is a separate,
roughly 2.8--2.9 s phase. The equal FP64 legacy-to-fast transport ratios are
about 1.42x, 1.71x, and 1.84x for r4/r6/r8 respectively, but these are timing
observations only until the required precision-matched moment and conductivity
comparisons are attached.

No GPU rows are included. `nvidia-smi` is present on the host, but reports
that it cannot communicate with the NVIDIA driver, so the driver correctly
ran no CUDA route and did not silently fall back to CPU. The JSON result is
machine-local evidence at `results/benchmarks/kpm_g12_pt_full.json`; it is
not a source-controlled timing artifact. A CUDA-host rerun with `--gpu` is
required before answering the GPU FP32/FP64 and equal-precision speedup
questions.

The campaign's `correctness_status` remains
`attach validated moment/conductivity comparison before speedup claim` for
all rows. A profile `PASS` establishes timing integrity only; it does not
establish physical equivalence between the legacy and fast routes.

## Verification

The implementation was checked with:

```text
cmake --build build-openblas-mpi-diagnose --target rslmto.x -j2
cmake --build build-acc01-cuda --target rslmto.x -j2
cmake --build build-acc01-cuda --target UnitKpmTransport -j2
ctest --test-dir build-acc01-cuda -R 'UnitKpmTransport|BenchmarkHarness' --output-on-failure
python3 tests/benchmarks/test_benchmark_harness.py
python3 -m py_compile tests/benchmarks/kpm_g12_transport.py
```

The CUDA build validates the expanded stochastic profiling ABI, while the
serial Pt run validates a real production profile and the parser's closure
checks. The full CPU campaign above validates all 180 measured profiles with
the same parser and closure contract.
