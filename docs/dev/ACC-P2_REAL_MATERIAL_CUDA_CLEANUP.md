# ACC-P2 real-material CUDA cleanup

ACC-P1b was classified before the P2 edits. The completed local probe found no
usable CUDA device, so its FP32 SCF rows are `unsupported`, not `acceptable`
or `unacceptable`. No FP32 production claim is inferred from that run. P2
therefore keeps FP64 as the production precision and uses FP64 Zheevd for the
large-matrix optimization target. The FP32 benchmark/eigensystem evidence is
retained in the ACC-P1b driver for a CUDA-host rerun.

The primary P2 workload is reciprocal SCF with eigenvectors requested. The
values-only route remains available through ACC-P0 for bands and solver
characterization, but it is not used to claim SCF acceleration.

## Instrumentation and ownership

`source/cuda/reciprocal_cuda.cpp` now owns persistent solver resources and
reports interval-local metrics for:

| component | reported field |
|---|---|
| H(k) assembly | `T_Hk_CPU_s` |
| host staging/conversion | `T_host_staging_s`, `H64_to_H32_s` |
| H2D | `H2D_s`, `H2D_bytes` |
| eigensolver | `solver_s` |
| D2H values/vectors | `T_D2H_values_s`, `T_D2H_vectors_s`, byte counters |
| host synchronization | `T_sync_s` |
| remaining backend/call overhead | `T_other_backend_s` |
| end-to-end request | `T_total_s` |

CUDA events are created and destroyed once at backend-context lifetime. Device
Hamiltonian, eigenvalue, info, and solver-work buffers grow to the required
capacity and are reused for smaller/equal requests. The benchmark records
allocation, workspace-query, event, and pinned-buffer counters both before
the measured interval and after it.

Pinned staging is opt-in with `RSLMTO_CUDA_PINNED_HOST=1` and is gated to
`n >= 486`. It is persistent, backend-owned, and does not register application
arrays. The P2 driver compares pageable and pinned staging for Fe 3^3, 4^3,
and 5^3 with vectors requested; it should be retained only when the
`total_improvement` is meaningful and positive.

## Synchronization dependency map

The CUDA solve enqueues H2D, cuSOLVER, and D2H operations on its one backend
stream, then performs `cudaStreamSynchronize` before returning host
eigenpairs. `cuda_backend_execute_batch` returns only after that solve has
completed, and the Fortran Fourier consumers copy `result%eigenvalues` and
`result%eigenvectors` immediately afterward. The benchmark uses the same
contract. The redundant Fortran `execution_backend%synchronize()` calls at
those result-consumption boundaries were removed; no asynchronous pipeline or
overlap assumption was introduced. The explicit synchronization API remains
available for callers that submit other backend work without a host-result
return boundary.

## Reproducible campaign

On a CUDA host:

```bash
python3 tests/benchmarks/accp2_real_material.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda \
  --output-dir results/benchmarks/accp2 \
  --warmups 2 --repetitions 5
```

For a CPU-only evidence run, add `--skip-cuda`; it produces the CPU control
rows and leaves CUDA budget/speedup fields unavailable rather than inventing
GPU timings. An optional `--before-results` file supplies the prior persistent
campaign to the mandatory before/after table. Without it, `before_total` is
null by policy.

The generated `accp2_results.json` contains the full timing budget for Fe
primitive, Fe 3^3 (`n=486`), Fe 4^3 (`n=1152`), and Fe 5^3 (`n=2250`), plus
the before/after table, pinned comparison, actual transferred bytes, and the
allocation audit. For complex FP64 vectors, the expected D2H payload is
`16*n^2` bytes per matrix, which is the evidence passed forward to ACC-11.
The H(k) fraction is recorded for ACC-12; ACC-P2 introduces no H(k) GPU
kernel.
The P2 wrapper suppresses the separate, expensive ACC-P1 eigenpair validation
inside each timing campaign (`--skip-cuda-validation`); the focused CUDA CTest
suite and the ACC-P1/P1b validation drivers remain the correctness gates. The
pinned comparison reuses the pageable run's LAPACK controls and runs only the
CUDA rows in its second pass.

## Current result status

The source/build checks pass on the available host. CUDA correctness and
timing rows remain unavailable because the host reports zero usable CUDA
devices; a CUDA-host run is required before accepting a performance checkbox.
`compute-sanitizer` is installed, but its memcheck run also stopped before the
first instrumented API call for the same reason. No GPU before/after speedup,
pinned-memory win, or FP32 scientific pass is claimed from this host; those
fields remain unavailable until the mandatory vector campaign runs on a
CUDA-equipped host.
