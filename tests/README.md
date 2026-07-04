# Tests

Top-level index for the test suites. Each suite has its own README with
the case-file format and how to add cases:

- `regression/cases.json` — backend regression matrix (10 cases: 5 CPU
  Chebyshev backends × hoh × ccor_2c + block sp/dp), compared against
  references from the trusted pre-refactor build. See
  `tests/regression/run_matrix.py` and `run_gpu_matrix.py` (below) for the
  runners.
- [`scf/README.md`](scf/README.md) — example SCF suite (bulk/surface/impurity, all
  recursion modes).
- [`postproc/README.md`](postproc/README.md) — example post-processing suite (exchange,
  conductivity, bands, DOS, orbital moments, PAOFLOW import).

## GPU coverage strategy

GitHub-hosted CI runners have no GPU, so GPU coverage splits into three tiers:

1. **MKL (runs in CI):** the `mkl` job in `.github/workflows/tests.yml` installs
   Intel oneAPI MKL, configures with `-DENABLE_MKL_KERNELS=ON`, and runs the
   `mkl_batch`/`mkl_sparse` regression cases plus one MKL-backend example case.
2. **CUDA compile-only (runs in CI):** the `cuda_compile` job installs
   `nvidia-cuda-toolkit` (no driver needed to build), configures with
   `-DENABLE_CUDA_PLUGIN=ON`, and runs the CPU regression subset against that
   binary. This catches interface drift between `rsrec_cuda_plugin.f90`,
   `source/cuda/rsrec.h`, and `source/cuda/rsrec_gpu.cu` on every push without
   needing GPU hardware — the most common GPU-side breakage mode.
3. **GPU (manual, off-CI):** `tests/run_gpu_matrix.sh [build_dir]` — run by hand
   on a workstation with an NVIDIA GPU, after building with
   `-DENABLE_CUDA_PLUGIN=ON -DRUN_REG_TESTS=ON`. It runs the CPU regression
   backend matrix against the CUDA-enabled binary, then
   `tests/regression/run_gpu_matrix.py`, which re-runs every
   `recur='chebyshev'` regression case with `gpu_plugin=true` and diffs
   `etot`/`ws_r`/`vmad` against the same case run with `gpu_plugin=false`
   (`gpu_plugin` bypasses `cheb_backend` entirely and dispatches straight to
   CUDA, so this is a CPU-vs-GPU consistency check, not a comparison against
   the committed CPU references). Treat this as the manual pre-release step
   before shipping a build meant to run on GPU hardware; it is not wired into
   any CI-required check.
