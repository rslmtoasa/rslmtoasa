/* ===========================================================================
 * hambuild.h -- GPU Hamiltonian-assembly backend ABI for RS-LMTO-ASA.
 *
 * Scope (see docs/HAMILTONIAN_GPU_PORT_BLUEPRINT.md):
 *   Move the Hamiltonian *assembly* (build_lsham / build_obarm / build_enim /
 *   build_bulkham / build_locham / CCOR) onto the GPU so the ee/eeo/hall/...
 *   arrays are built on the device and feed straight into the rsrec_cuda
 *   recursion backend without a host round-trip.
 *
 * This is Phase-0 scaffolding: the context lifecycle and error plumbing are in
 * place; the numerical set_* / onsite / bulk / local / ccor entry points are
 * declared so the Fortran shim and the CPU/GPU dispatchers can be filled in
 * phase by phase. The CPU path in hamiltonian.f90 is the regression oracle and
 * bit-identical when control%gpu_hambuild = .false.
 *
 * Layering mirrors the recursion plugin exactly:
 *   hambuild_cuda_plugin.f90  (iso_c_binding shims)
 *     -> hambuild_cuda.cpp     (backend dispatcher: cpu vs gpu)
 *          -> hambuild_gpu.cu   (CUDA kernels; HIP hooks via macros)
 *          -> hambuild_cpu.cpp  (thin CPU fallback -- optional)
 *   hambuild_stub.c            (no-op symbols for CPU-only builds)
 * =========================================================================== */
#ifndef HAMBUILD_H
#define HAMBUILD_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hambuild_ctx hambuild_ctx;

/* Backend selection (mirrors rsrec_cuda's enum style). */
enum {
    HAMBUILD_BACKEND_CPU = 0, /* reference / fallback           */
    HAMBUILD_BACKEND_GPU = 1  /* CUDA (or HIP) device kernels   */
};

/* --- Context lifecycle ------------------------------------------------------
 * Dimensions:
 *   nb     spin-orbital block size (2*norb)
 *   norb   orbitals per spin channel
 *   nnmax  max neighbours per site (build_bulkham/build_locham loop bound)
 *   ntype  number of atom types
 *   nmax   number of local (impurity interaction-zone) sites; 0 for bulk-only
 *   device CUDA device ordinal
 */
hambuild_ctx *hambuild_cuda_create(int nb, int norb, int nnmax, int ntype,
                                   int nmax, int device);
void hambuild_cuda_destroy(hambuild_ctx *ctx);
const char *hambuild_cuda_last_error(void);

int hambuild_cuda_set_backend(hambuild_ctx *ctx, int backend);

#ifdef __cplusplus
}
#endif

#endif /* HAMBUILD_H */
