# ACC-01 — Existing RS CUDA coverage

**Scope:** audit of the existing RS CUDA plugin, not a redesign of its
architecture. **Audit revision:** working tree after `7da224b9`.

## Result

The plugin has one Fortran wrapper, one public C ABI, one C++ dispatcher, and
one CUDA implementation. The API is broader than the previous production
test matrix: all recursion and reconstruction entry points are implemented,
but only part of that surface had a numerical CPU/GPU contract test.

The audit added:

- `tests/cuda/test_plugin_surface.py`, a no-GPU source/API contract test;
- a guard preventing the unsupported structured FFT/conv orbital route from
  entering CUDA;
- CPU fallback for Green reconstruction when `gpu_plugin=.true.` is used in
  an executable built without the CUDA plugin;
- this support and evidence matrix.

The source/API test passes locally. The low-level validator is now a typed ABI
contract and includes direct NumPy comparisons for orbital moments and all
four DOS/GF reconstruction entry points. The hardware campaign subsequently
passed all low-level routes (15 reported comparisons) at approximately 1e-15 relative error, the
CPU regression matrix at 10/10, and the production CPU/GPU Chebyshev matrix
at 8/8. The two optional MKL cases were correctly skipped because the build
reported `ENABLE_MKL_KERNELS=OFF`.

## Public route inventory

| Route | CPU production counterpart | CUDA entry point | Inputs/limitations | Existing evidence |
|---|---|---|---|---|
| Scalar Lanczos | `recursion_haydock.f90:recur` | `rsrec_scalar_lanczos` | `nsp=1`; scalar orbital chains | `tests/cuda/rsrec_validate.py`; production GPU matrix does not exercise it directly |
| Block Lanczos | `recursion_haydock.f90:recur_b` | `rsrec_block_lanczos` | block coefficients; HOH supported; `ccor_2c+hoh` rejected | `tests/cuda/rsrec_validate.py`; CPU regression build matrix |
| Block intersite Lanczos | `recursion_transport.f90:recur_b_ij` | `rsrec_block_lanczos` | four signed pair starts; same HOH/CCOR guards | pair initialization is present in `rsrec_validate.py`; no production GPU pair-output comparison |
| Chebyshev moments | `recursion_chebyshev.f90:chebyshev_recur` | `rsrec_chebyshev_moments` | one starting block per call; structured backend is FP64 and `nrhs=nb` | `rsrec_validate.py`; `run_gpu_matrix.py` through SCF observables |
| Chebyshev intersite moments | `recursion_transport.f90:chebyshev_recur_ij` | `rsrec_chebyshev_moments` | four signed pair starts; structured route supports only compatible periodic ee-only input | pair-start unit in `rsrec_validate.py`; no committed production GPU pair contract |
| Stochastic moments | `recursion_transport.f90:compute_moments_stochastic` | `rsrec_stochastic_moments` | requires Hamiltonian and velocity; HOH requires `vo_a/vo_b`; device memory scales with `cond_ll` | `rsrec_validate.py`; no end-to-end GPU conductivity comparison |
| Orbital moments | `recursion_transport.f90:chebyshev_orbital_mod` | `rsrec_orbital_moments` | host builds the orbital-current left state; structured FFT/conv is unsupported; HOH is supported by block-ELL | direct NumPy contract in `rsrec_validate.py`; hardware pending |
| Chebyshev DOS | `green_chebyshev.f90:chebyshev_green_core` | `rsrec_chebyshev_dos` | Jackson reconstruction; FP32 default, FP64 validation; real energy mesh | direct array contract plus production SCF/DOS matrix; hardware pending |
| Chebyshev eta GF | `green_chebyshev.f90:chebyshev_green_core` / intersite eta consumer | `rsrec_chebyshev_gf_eta` | host builds complex transfer factors; eta contour and pair batching | direct NumPy contract in `rsrec_validate.py`; hardware pending |
| Block DOS | `green_block.f90:block_green_core` | `rsrec_block_dos` | continued fraction and diagonal terminator inputs; FP32 default, FP64 validation | direct NumPy contract in `rsrec_validate.py`; hardware pending |
| Block eta GF | `green_block.f90:calculate_intersite_gf_eta_gpu` | `rsrec_block_gf_eta` | one Fermi energy, multiple complex eta values, four combinations per pair | direct NumPy contract in `rsrec_validate.py`; hardware pending |

The low-level C API also exposes lifecycle/configuration operations
(`create`, `destroy`, backend selection, lattice/Hamiltonian/velocity upload,
and precision selection). These are covered by the ABI surface test and by
the production callers; they are not physical-output contracts by themselves.

## Backend and combination policy

| Combination | Status | Reason/evidence |
|---|---|---|
| CSR/BSR + ordinary Hamiltonian | supported path | Main RS CUDA route; BSR uploads through the existing sparse export seam |
| CSR/BSR + HOH | supported for block-ELL recursion/reconstruction | Two-sweep CUDA implementation; sensitive block normalization remains FP64 |
| `ccor_2c` without HOH | supported through merged operator work arrays | `ensure_ccor_operator_blocks` supplies the effective operators |
| `ccor_2c` with HOH | unsupported | Explicit `gpu_plugin_ready` rejection; no fallback is mislabeled as GPU support |
| local-axis rotation | unsupported | Explicit `gpu_plugin_ready` rejection |
| scalar Lanczos with `nsp != 1` | unsupported | Explicit `require_nsp1` guard |
| FFT/conv with non-periodic, partial-periodic, impurity, or non-ee input | unsupported | Explicit periodic/full-direction/`nmax=0` checks |
| FFT/conv orbital moments | unsupported | CUDA implementation rejects structured orbital moments; ACC-01 now rejects before dispatch |
| stochastic HOH without velocity-orthogonalization operators | unsupported | CUDA ABI returns an explicit error |
| GPU requested in a CPU-only executable for recursion or Green reconstruction | CPU fallback | Recursion and Green dispatch now both check compiled-plugin availability |

No CPU fallback is counted as GPU coverage for a route whose numerical CUDA
entry point is not selected. Unsupported combinations remain visible through
warnings or explicit errors.

## Test and benchmark evidence

| Evidence | What it checks | Hardware required | Status |
|---|---|---:|---|
| `CudaPluginSurface` | Fortran binding ↔ C header ↔ C++ dispatcher ↔ CUDA symbol and production route markers | no | pass |
| `tests/cuda/rsrec_validate.py` | typed public ABI; matvec, Chebyshev site/pair moments, Block Lanczos, scalar Lanczos, stochastic/orbital moments, Chebyshev DOS/GF, and Block DOS/GF against NumPy mirrors | yes for CUDA build/library | all 15 reported comparisons pass at ~1e-15 |
| `tests/cuda/build_and_validate.sh` | reproducibly builds the standalone CUDA library and runs the complete low-level validator | yes | pass |
| `tests/run_gpu_matrix.sh` | CPU regression backend matrix plus Chebyshev CPU/GPU `etot`, `ws_r`, `vmad` comparison | yes | CPU 10/10; GPU consistency 8/8 |
| ACC-00 harness | repeatable timing schema and comparison tooling | no for harness; yes for GPU records | GPU Block/Chebyshev records captured |

The correctness and GPU performance campaigns are complete. The first
manifest campaign did not set `control%gpu_plugin=true`; those five records
are CPU timings and are retained as baselines. The explicit GPU campaign
captured the following medians on an NVIDIA RTX A4000:

| Route | CPU median (s) | GPU median (s) | CPU/GPU speedup |
|---|---:|---:|---:|
| RS Block Fe | 3.267535 | 2.872801 | 1.14x |
| RS Chebyshev Si | 11.020424 | 1.922788 | 5.73x |

These are end-to-end fixture timings, including setup and output checks, from
the same build and three repetitions after one warm-up. The build links Intel
oneMKL for base BLAS/LAPACK, but `ENABLE_MKL_KERNELS=OFF`; the optional MKL
backend cases remain correctly excluded.

## Reproduction command

On a CUDA workstation, build with `-DENABLE_CUDA_PLUGIN=ON`, run the source
contract test and existing manual matrix, then record route timings with the
ACC-00 harness. Do not update CPU references from GPU differences.

```bash
python3 tests/cuda/test_plugin_surface.py
tests/cuda/build_and_validate.sh /tmp/rslmto-acc01-lowlevel
tests/run_gpu_matrix.sh build-cuda
python3 tests/benchmarks/benchmark_harness.py run-manifest \
  --manifest tests/benchmarks/manifest.json \
  --binary build-cuda/bin/rslmto.x \
  --build-dir build-cuda --warmups 1 --repetitions 3 \
  --gpu-plugin \
  --scratch-root /tmp/rslmto-acc01-bench \
  --output-dir results/benchmarks/acc01-gpu
```

The production build must be configured with `-DENABLE_CUDA_PLUGIN=ON`.
For meaningful optional-MKL timings, use the project’s MKL-enabled toolchain
and confirm `ENABLE_MKL_KERNELS=ON` in `CMakeCache.txt`; a CUDA build with that
option off must not be treated as an MKL benchmark.
