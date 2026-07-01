/* ===========================================================================
 * hambuild_gpu.cu -- CUDA (with HIP hooks) device kernels for Hamiltonian
 *                    assembly.  Phase-0 scaffold: no numerics yet.
 *
 * Real kernels are added one blueprint phase at a time:
 *   Phase 1   pack_spinor, hcpx_dev  (device primitives, reused everywhere)
 *             + on-site blocks: build_lsham / build_obarm / build_enim
 *   Phase 1.5 geometry maps (clusba filter + CUB compaction, neigh_map builder)
 *   Phase 2   build_bulkham / build_locham hot loop (ham0m_nc batched,
 *             cuBLAS ZgemmStridedBatched for the hoh matmuls)
 *   Phase 3   CCOR (eecc / hallcc)
 *   Phase 4   LDA+U / Hubbard add-in
 *
 * --- HIP portability hooks ---------------------------------------------------
 * This file is authored in CUDA. To target AMD via HIP, hipify-perl converts
 * the cuda* symbols to hip*; the shim below lets the same source compile under
 * either toolchain. Guard device-specific intrinsics behind these aliases
 * rather than raw cuda* / __CUDA* names so the HIP port stays mechanical.
 * =========================================================================== */

#if defined(__HIP_PLATFORM_AMD__) || defined(USE_HIP)
  #include <hip/hip_runtime.h>
  #define HB_MALLOC       hipMalloc
  #define HB_FREE         hipFree
  #define HB_MEMCPY       hipMemcpy
  #define HB_MEMCPY_H2D   hipMemcpyHostToDevice
  #define HB_MEMCPY_D2H   hipMemcpyDeviceToHost
  #define HB_DEVICE_SYNC  hipDeviceSynchronize
  #define HB_SUCCESS      hipSuccess
#else
  #include <cuda_runtime.h>
  #define HB_MALLOC       cudaMalloc
  #define HB_FREE         cudaFree
  #define HB_MEMCPY       cudaMemcpy
  #define HB_MEMCPY_H2D   cudaMemcpyHostToDevice
  #define HB_MEMCPY_D2H   cudaMemcpyDeviceToHost
  #define HB_DEVICE_SYNC  cudaDeviceSynchronize
  #define HB_SUCCESS      cudaSuccess
#endif

/* Device primitives declared here for Phase 1. Kept as forward stubs so the
 * translation unit compiles cleanly and the ABI is stable. */

/* pack_spinor: assemble the 2x2-spin nb x nb block from the cartesian
 * (H0, Hx, Hy, Hz) 4-component per-(orbital,orbital) block. See §2a of the
 * blueprint. Implemented in Phase 1. */

/* hcpx_dev: on-device cart<->sph similarity transform per l-block (math.f90
 * hcpx). Constant transform matrices uploaded once to device constant memory.
 * See §2b. Implemented in Phase 1. */

/* No host-callable kernels are exported yet; hambuild_cuda.cpp forwards nothing
 * to this unit during Phase 0. This keeps the .cu in the build (validating the
 * CUDA toolchain wiring) without introducing numerics. */
