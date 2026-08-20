# KPM-G1.3 GPU reconstruction

## Implementation status

The CUDA transport path now follows the G1.1 flattened reconstruction:

```text
G(NE,M*M) * U(M*M,nb) -> C(NE,nb)
```

The production formula, `q = n + (m-1)*M` layout, factor, and per-type/trace
semantics are unchanged. The optimized CUDA route:

- retains only diagonal-packed `U(M*M,nb)` matrices in the backend context;
- generates Gamma in energy tiles on the device from the production basis;
- dispatches FP64 through cuBLAS `ZGEMM` and FP32 through `CGEMM3M`;
- downloads only `C_block( BE, nb )` by default;
- widens FP32 `C_block` to the existing FP64 integration/output layer;
- leaves energy integration and file output on the CPU;
- keeps the old full-moment API and an explicit `gpu_moment_download` diagnostic
  option.

Residency is invalidated by Hamiltonian, velocity, precision, or explicit
clear operations. No device pointer is exposed to Fortran. `RSLMTO_KPM_GPU_BE`
selects the requested energy tile; the default is 64 and the runtime clamps it
to available device memory and the final partial block.

The profile now reports `D_gpu_mu_pack`, `D_reconstruction_D2H`,
`gpu_energy_block`, reduced-result D2H bytes, and the CUDA reconstruction
backend/precision metadata.

## Validation evidence

Completed:

```text
CPU build: build-openblas-mpi-diagnose
CUDA build: build-acc01-cuda
CUDA plugin surface test: PASS
benchmark harness contract test: PASS
UnitKpmTransport: PASS
CudaPluginSurface: PASS
```

The CUDA build compiles and links the resident-moment ABI, tiled Gamma kernel,
and both cuBLAS precision paths. The production-source checks also verify that
the optimized route contains no scalar GPU Gamma-mu loop and does not upload a
full Gamma tensor.

## Host-GPU correctness

The low-level CUDA harness was run on an RTX A4000 with CUDA 13.3. All prior
CUDA routes passed. The new G1.3 checks passed with:

| Check | Result |
| --- | ---: |
| FP64 resident moments versus full diagnostic moments | relative error 0 |
| FP64 reconstruction, `NE=5`, `M=6`, `BE=3` | relative error `6.66e-16` |
| FP32 reconstruction, same partial tile | relative error `3.77e-7` |
| optimized resident moment full-tensor D2H | 0 bytes |
| FP64 reduced-result D2H | 480 bytes |

The first host run exposed and fixed a real Gamma layout defect: the CUDA fill
kernel had stored energy-major data while cuBLAS consumed column-major data.
The corrected test now covers the final partial block and no longer reports the
previous M>1 mismatch.

## Pt performance evidence

The production Pt fixture was run with `M=500`, `NE=2510`, `lld=150`, spin,
and `per_type` on the RTX A4000 host. All rows below had `PROFILE_STATUS=PASS`
and `reconstruction_backend=cuda_blas`; the r4 tile sweep used one warmup and
two measured repetitions.

| Route | BE | `P_gamma` (s) | `P_reconstruction_total` (s) | `T_transport_total` (s) | result D2H |
| --- | ---: | ---: | ---: | ---: | ---: |
| r4 FP64 | 16 | 0.702 | 0.707 | 13.993 | 722,880 B |
| r4 FP64 | 32 | 0.718 | 0.525 | 16.397 | 722,880 B |
| r4 FP64 | 64 | 0.765 | 0.689 | 16.155 | 722,880 B |
| r4 FP64 | 128 | 0.725 | 0.696 | 14.165 | 722,880 B |
| r4 FP32 | 128 | 0.727 | 0.022 | 5.237 | 361,440 B |

The BE rows were collected across the two available A4000 devices, so the
small differences should not be treated as a universal optimum. BE=16 and
BE=128 were the most stable r4 FP64 points in this sample; the default remains
64 pending a larger repeated sweep.

The FP32 BE=128 scale check remained nearly N-independent in reconstruction:

| Size | `P_moments_total` (s) | `P_reconstruction_total` (s) | `T_transport_total` (s) |
| --- | ---: | ---: | ---: |
| r4 | 0.439 | 0.021 | 5.225 |
| r6 | 1.022 | 0.020 | 6.190 |
| r8 | 2.248 | 0.021 | 7.751 |

Single-measurement FP64 scale checks at BE=128 gave `T_transport_total` of
37.023 s for r6 and 73.386 s for r8, with reconstruction times of 0.699 s
and 0.699 s respectively. These are performance evidence, not matched
CPU/GPU speedup claims; the runs used fresh stochastic processes and were not
paired with a same-process CPU reference.
