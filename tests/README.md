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

## Regenerating references

Each suite's `references/<TestName>/ref.json` (or `.nml` for the legacy
regression baselines) is a **committed value**, not a cache — a test failure
means the code changed, not that the reference is stale. Only regenerate a
reference when you *intend* the physics/output to change, e.g. a deliberate
numerical improvement, a namelist default change, or adding a brand-new case
that has no reference yet. The workflow:

1. Build a **trusted** binary — one you have independently verified is
   correct for the change you're making (not just "it compiles").
2. Regenerate: `python3 tests/generate_references.py --binary <path> --cases-json <suite>/cases.json --references-dir <suite>/references [--case <TestName> ...]`.
   Omit `--case` to regenerate every case in that suite.
3. **Review the diff of the reference files like any other code change.**
   A `ref.json` diff is a physics/behavior diff — read it with the same
   scrutiny as a source diff, not as machine-generated noise to skip past.
   If a value moved and you can't explain why, that's a bug to find, not a
   reference to accept.
4. Commit the reference update in the same commit as (or immediately
   following) the change that caused it, with a commit message explaining
   *why* the values moved.

Never regenerate references to make an unexplained failure go away — that
converts a real regression into a silently-accepted one.

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

## CI trigger and matrix hygiene

- **Triggers:** `binaries.yml` builds on `push`/`pull_request`/`release`.
  `tests.yml` runs the example/regression matrix three ways: via
  `workflow_run` after a successful `binaries.yml` push build (downloads the
  built binary), directly on `pull_request` (builds from source itself,
  since there's no upstream artifact yet for an unmerged PR), and nightly on
  `schedule` (also self-built). Both workflows have `concurrency` groups with
  `cancel-in-progress: true` so superseded pushes/PR updates don't queue up
  redundant runs.
- **Quick vs. full:** `pull_request` runs only cases tagged `"quick": true`
  in `cases.json` (~10 representative cases across bulk/surface/impurity/
  k-space/exchange/paoflow/hubbard, picked to keep total runtime under a
  couple of minutes) via `ctest --label-regex quick`. Every other trigger
  (push, nightly schedule, manual dispatch) runs the full matrix
  (`--label-regex example`). This keeps PR turnaround fast while nothing
  rots — the nightly run still exercises every case every day.
- **Windows:** `binaries.yml` builds a Windows binary but never tested it.
  The `windows_quick` job in `tests.yml` runs the same quick subset under
  MSYS2 (mirroring `binaries.yml`'s Windows build steps, with
  `mingw-w64-x86_64-python`/`-python-pip` added for `f90nml`), with
  `ENABLE_MPI=OFF` to avoid MSYS2 MPI toolchain coupling issues (same as the
  Windows build job). This job carries `continue-on-error: true`: it was
  authored and validated on Linux only (no Windows runner was available to
  verify MSYS2 path/toolchain behavior end-to-end), so a first-run failure
  there should not block the rest of CI. If it turns out to need real
  iteration to get green, the fallback is to drop it and document Windows as
  build-only instead.
