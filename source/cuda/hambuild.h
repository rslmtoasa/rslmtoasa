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

/* --- Phase 1.5: geometry maps (once per run; geometry is fixed across SCF) ----
 *
 * The per-SCF hot loop gathers from resident tables built here, so the neighbour
 * search (clusba distance filter + hmfind match) runs once at setup rather than
 * every SCF. All arrays are Fortran column-major.
 *
 * Static geometry upload (integer arrays are 1-based Fortran indices, kept as-is;
 * the kernels convert to 0-based internally):
 *   cr     : real64 3*kk   cluster positions (in lattice units; *alat on device)
 *   num    : int    kk     bravais/structure-constant type per cluster atom
 *   iz     : int    kk     atom type per cluster atom
 *   nn     : int    ndi*nn_max  neighbour table (col-major (ndi,nn_max)); nn(ia,1)
 *                          holds the neighbour count nr for atom ia
 *   atlist : int    ntype  cluster atom index (ia) for each atom type
 * Scalars: kk, nn_max, ndi, alat, r2. pbc!=0 selects wrapped coordinate diffs
 * (currently the non-pbc straight difference is implemented; pbc falls back to
 * host-provided vet -- see set_geometry_vet).
 */
int hambuild_cuda_set_geometry(hambuild_ctx *ctx, const double *cr,
                               const int *num, const int *iz, const int *nn,
                               const int *atlist, int kk, int nn_max, int ndi,
                               double alat, double r2, int pbc);

/* For pbc (or any case where the host already has the exact vet per (type,m)),
 * upload the precomputed displacement vectors: vet is real64 3*nn_max*ntype,
 * col-major (3, nn_max, ntype). When provided, the geometry-map builder uses
 * these directly instead of recomputing the straight coordinate difference. */
int hambuild_cuda_set_geometry_vet(hambuild_ctx *ctx, const double *vet);

/* Build the resident neigh_map + vet tables from the uploaded geometry (runs the
 * clusba filter + hmfind match on-device). Must be called after set_geometry
 * (and set_geometry_vet if pbc) and before the Phase-2 hot loop. */
int hambuild_cuda_build_geometry_maps(hambuild_ctx *ctx);

/* Validation getter: copy the built maps back to the host.
 *   valid : int    nn_max*ntype   1 if neighbour found (hmfind ni), else 0
 *   shell : int    nn_max*ntype   sbar shell index m (1-based, as CPU uses)
 *   ino   : int    nn_max*ntype   sbar type index num(ia) (1-based)
 *   vet   : real64 3*nn_max*ntype the displacement vectors used
 * Any pointer may be NULL. Lets Phase 1.5 be diffed against CPU hmfind exactly. */
int hambuild_cuda_get_geometry_maps(hambuild_ctx *ctx, int *valid, int *shell,
                                    int *ino, double *vet);

/* --- Phase 2: bulk/local hot loop (chbar_nc: ham0m_nc + hcpx + pack) ----------
 *
 * Per-type potential inputs for the neighbour Hamiltonian (ham0m_nc). All
 * complex arrays are norb per type, contiguous per type, interleaved re/im:
 *   wx0, wx1     : overlap-weight combos (up+dn, up-dn)
 *   cx0, cx1     : band-centre combos    (on-site, non-hoh)
 *   cex0, cex1   : band-centre combos    (on-site, hoh)
 *   mom          : real64 3 per type (moment direction; overridden by spiral)
 * Spin spiral: q_ss (3) and theta_ss (1) scalar; when |q_ss|>1e-5 or
 * |sin(theta_ss)|>1e-8 the kernel recomputes mom_ia/mom_ja from cr and q_ss.
 */
int hambuild_cuda_set_potential_bulk(hambuild_ctx *ctx, const void *wx0,
                                     const void *wx1, const void *cx0,
                                     const void *cx1, const void *cex0,
                                     const void *cex1, const double *mom,
                                     const double *q_ss, double theta_ss);

/* Upload the (real part of the) structure constants sbar. Refreshed per SCF.
 * sbar is complex norb*norb*nm_store*num_types, col-major; the kernel uses
 * real(sbar) matching the CPU real() cast. nm_store = 3rd dim, ntot = 4th dim. */
int hambuild_cuda_set_sbar(hambuild_ctx *ctx, const void *sbar, int nm_store,
                           int ntot);

/* Build the bulk neighbour Hamiltonian on-device (batched over type*neighbour),
 * then pack into ee/hxc and, if hoh, the two hoh matmuls into eeo/eeoee (cuBLAS
 * ZgemmStridedBatched). obarm must already be resident (from a prior onsite
 * call with want_obarm). Outputs copied back to the host buffers:
 *   ee, hxc          : nb*nb*nn_max*ntype complex
 *   eeo, eeoee       : nb*nb*nn_max*ntype complex (only written when hoh)
 * eeo/eeoee may be NULL when hoh==0. */
int hambuild_cuda_bulk(hambuild_ctx *ctx, int hoh, void *ee, void *hxc,
                       void *eeo, void *eeoee);

#ifdef __cplusplus
}
#endif

#endif /* HAMBUILD_H */
