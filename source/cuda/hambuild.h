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

/* --- Phase 1: on-site blocks (build_lsham / build_obarm / build_enim) --------
 *
 * All complex arrays are passed as interleaved double (re, im) pairs, matching
 * Fortran complex(rp) with rp=real64 and column-major storage. Sizes below use
 * norb (orbitals per spin channel) and nb = 2*norb (spinor block).
 *
 * Constant, geometry-independent inputs uploaded once (they only depend on the
 * basis, not the potential):
 *   v, vc : norb*norb complex, the hcpx cart<->sph transform matrices (the exact
 *           matrices hcpx() builds for this norb; identity when norb is neither
 *           4 nor 9). Column-major.
 *   lx, ly, lz : norb*norb complex, orbital angular-momentum operators in the
 *           CARTESIAN real-harmonic basis (as L_x/L_y/L_z before hcpx). The
 *           lsham kernel applies hcpx_dev to these internally, mirroring the CPU.
 */
int hambuild_cuda_set_constants(hambuild_ctx *ctx, const void *v, const void *vc,
                                const void *lx, const void *ly, const void *lz);

/* Per-type potential inputs for the on-site blocks. Each pointer covers all
 * ntype types, contiguous per type:
 *   xi_p, xi_d, rac : real64, 2 per type (spin up/down)  [rac used iff orb_pol]
 *   obx0, obx1      : complex, norb per type (diagonal seeds for obarm)
 *   cx, cex         : complex, norb*2 per type (col-major (norb,2): spin up/down)
 *   mom             : real64, 3 per type (normalized moment direction)
 *   lmom            : real64, 3 per type (orbital moment; z used iff orb_pol)
 *   orb_pol         : int flag (0/1); enables the rac/lmom polarization term
 */
int hambuild_cuda_set_potential_onsite(hambuild_ctx *ctx, const double *xi_p,
                                       const double *xi_d, const double *rac,
                                       const void *obx0, const void *obx1,
                                       const void *cx, const void *cex,
                                       const double *mom, const double *lmom,
                                       int orb_pol);

/* Build the three on-site blocks. Outputs are nb*nb*ntype complex, col-major,
 * written to the host buffers provided (device->host copy inside). Any of the
 * output pointers may be NULL to skip that block. */
int hambuild_cuda_onsite(hambuild_ctx *ctx, void *lsham, void *obarm,
                         void *enim);

#ifdef __cplusplus
}
#endif

#endif /* HAMBUILD_H */
