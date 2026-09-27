# B4 — GPU k-space: batched eigensolver + batched inversion

**Effort:** M. **Lead:** SONNET (the plugin pattern is established; kernels
are vendor-library calls), OPUS only if numerical policy questions arise.
**Depends on:** B2 defines the consumers. **Accelerates:** B2 both
backends, B3, B5, B8, B10, B11 (χ₀ sums) — one investment, many payoffs.

## 1. Design (follow the existing plugin pattern exactly)

One plugin `rsksp_gpu` exposing **three batched primitives** behind the
same C-API + `iso_c_binding` + CPU-fallback pattern already used by the
CUDA/MKL kernels (template: `rsrec_cuda_plugin.f90` and the MKL kernel
gating fixed in Phase 2 P3). Device memory reuse follows the
fingerprint/context lifecycle from Phase 1 (one shared context,
`ensure_context`); the Fortran side decides batching over k (and energy
for backend D).

1. **`batched_zheevd`** — diagonalize N_batch Hermitian N×N matrices
   (H(k) over a k-tile). CUDA: cuSOLVER `cusolverDnZheevdBatched` (or
   syevjBatched for small N — benchmark both; Jacobi often wins at
   N ≤ 32); HIP: hipSOLVER equivalents.
2. **`batched_zgetrf_zgetrs`** — batched LU + solve for backend D's
   per-(k,z) resolvents and later CPA. cuBLAS
   `cublasZgetrfBatched`/`ZgetrsBatched`; HIP equivalents. (Solve, not
   full inversion, unless a consumer genuinely needs G's full matrix —
   backend D accumulating pair blocks usually does; provide `zgetri` too
   but prefer `getrs` where the RHS set is small.)
3. **`batched_zgemm`** — eigenvector-product accumulation for Lehmann
   fills and moment/GF assembly (B2.2/B5.1 inner loops).

**Explicitly out of scope:** GPU-porting tetrahedron bookkeeping or
symmetry reduction — host-side, cheap.

**GPU strux (from the original wishlist): gated on a profile.**
Structure-constant setup normally amortizes against recursion time; port
only if profiling a large interface/surface workload shows it as a
bottleneck, and then via the same plugin pattern at the
`lattice_strux.f90` ↔ `include_codes/strux_lib` boundary. Do not open
this task without the profile.

Rules:
- **CPU fallback mandatory** and bit-comparable in tests at fp64 (same
  LAPACK routines, looped). Build-time options mirror the existing
  `ENABLE_*_KERNELS` CMake pattern with the fatal-with-hint behavior
  validated in Phase 2 P10 — copy that gating verbatim.
- **Precision (Gate G-B4-1):** default fp64. A single build-time knob may
  later allow fp32 for backend-D energy scans (spectra tolerate it;
  eigen-decompositions for charge self-consistency do not). Do NOT
  implement mixed precision until Anders signs the policy — expose the
  knob's plumbing, wire only fp64.
- **Batching layout:** k-tiles sized to fill the device (thousands of
  18–72×N matrices); pinned host staging buffers; one H(k) build stream on
  CPU overlapped with device compute (double-buffer, as in `librsrec`'s
  lessons — avoid tiny per-matrix transfers, the Apple-Accelerate
  call-overhead lesson generalizes: batch or die).
- No physics lives here. If a physics question surfaces, it belongs in
  B2's spec — escalate.

## 2. Tests
1. Fallback equivalence: GPU vs CPU eigenvalues ≤ 1e-12 relative,
   eigenvectors up to phase (compare |U†U'| ≈ 1), on random Hermitian
   batches and on real bcc-Fe H(k) tiles.
2. B2's backend-equivalence and chain-oracle tests re-run with GPU enabled
   — zero tolerance changes allowed.
3. CI: compile-only GPU job (pattern from Phase 2 P3); the run-matrix
   script `tests/run_gpu_matrix.sh` gains the two new kernels for real-
   hardware validation.
4. Benchmark artifact: k-points/second vs CPU for N = 18, 36, 72;
   recorded in `docs/dev/gpu_kspace_bench.md`.

## 3. Tasks
- **B4.1 [SONNET]** C API header for all three primitives + CPU fallback +
  CMake gating + tests (runs everywhere, no GPU needed). *Kit:* this file;
  the existing plugin triplet as template (`rsrec_cuda_plugin.f90` +
  header + CMake stanza under `source/cuda/`).
- **B4.2 [SONNET]** CUDA implementations of the three primitives + stream/
  staging layer + `ensure_context` lifecycle. *Kit:* B4.1 files; the
  librsrec lesson summary in §1.
- **B4.3 [SONNET]** HIP port (mechanical; hipify then hand-verify API
  deltas). *Kit:* B4.2 files.
- **B4.4 [SONNET]** Wire into B2 backends behind the existing kernel-
  dispatch pattern; re-run B2 test suite; bench doc. *Kit:* B2 module;
  B4.1–B4.3.

## 4. Checklist
- [ ] B4.1 API + fallback + gating
- [ ] B4.2 CUDA
- [ ] B4.3 HIP
- [ ] B4.4 B2 wiring + bench (G-B4-1 signed if fp32 knob activated)
