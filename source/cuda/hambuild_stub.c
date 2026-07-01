/* ===========================================================================
 * hambuild_stub.c -- no-op fallback for the hambuild_* GPU assembly backend.
 *
 * hambuild_cuda_plugin.f90 is compiled unconditionally and binds the
 * hambuild_cuda_* C interface (see hambuild.h). The real implementations live
 * in cuda/hambuild_cuda.cpp + cuda/hambuild_gpu.cu and are only built with
 * -DENABLE_CUDA_PLUGIN=ON.
 *
 * For ordinary CPU builds we still need those symbols to link, but we do NOT
 * want to pull in any compute kernels. This stub supplies them:
 * hambuild_cuda_create returns NULL, which the Fortran ensure_context() detects
 * (c_associated test) and turns into a clean fatal error. Every other entry
 * point returns a nonzero error code. None of these are ever reached unless the
 * user requests the GPU path at runtime (control%gpu_hambuild = .true.).
 * =========================================================================== */
#include "hambuild.h"

static const char *g_msg =
    "hambuild: GPU Hamiltonian-assembly backend not built in this binary. "
    "Reconfigure with -DENABLE_CUDA_PLUGIN=ON to enable the gpu_hambuild path.";

const char *hambuild_cuda_last_error(void) { return g_msg; }

hambuild_ctx *hambuild_cuda_create(int nb, int norb, int nnmax, int ntype,
                                   int nmax, int device)
{
    (void)nb; (void)norb; (void)nnmax; (void)ntype; (void)nmax; (void)device;
    return 0; /* NULL ctx -> ensure_context() aborts cleanly */
}

void hambuild_cuda_destroy(hambuild_ctx *ctx) { (void)ctx; }

int hambuild_cuda_set_backend(hambuild_ctx *ctx, int backend)
{
    (void)ctx; (void)backend;
    return 1;
}

int hambuild_cuda_set_constants(hambuild_ctx *ctx, const void *v, const void *vc,
                                const void *lx, const void *ly, const void *lz)
{
    (void)ctx; (void)v; (void)vc; (void)lx; (void)ly; (void)lz;
    return 1;
}

int hambuild_cuda_set_potential_onsite(hambuild_ctx *ctx, const double *xi_p,
                                       const double *xi_d, const double *rac,
                                       const void *obx0, const void *obx1,
                                       const void *cx, const void *cex,
                                       const double *mom, const double *lmom,
                                       int orb_pol)
{
    (void)ctx; (void)xi_p; (void)xi_d; (void)rac; (void)obx0; (void)obx1;
    (void)cx; (void)cex; (void)mom; (void)lmom; (void)orb_pol;
    return 1;
}

int hambuild_cuda_onsite(hambuild_ctx *ctx, void *lsham, void *obarm,
                         void *enim)
{
    (void)ctx; (void)lsham; (void)obarm; (void)enim;
    return 1;
}

int hambuild_cuda_set_geometry(hambuild_ctx *ctx, const double *cr,
                               const int *num, const int *iz, const int *nn,
                               const int *atlist, int kk, int nn_max, int ndi,
                               double alat, double r2, int pbc)
{
    (void)ctx; (void)cr; (void)num; (void)iz; (void)nn; (void)atlist;
    (void)kk; (void)nn_max; (void)ndi; (void)alat; (void)r2; (void)pbc;
    return 1;
}

int hambuild_cuda_set_geometry_vet(hambuild_ctx *ctx, const double *vet)
{
    (void)ctx; (void)vet;
    return 1;
}

int hambuild_cuda_build_geometry_maps(hambuild_ctx *ctx)
{
    (void)ctx;
    return 1;
}

int hambuild_cuda_get_geometry_maps(hambuild_ctx *ctx, int *valid, int *shell,
                                    int *ino, double *vet)
{
    (void)ctx; (void)valid; (void)shell; (void)ino; (void)vet;
    return 1;
}

int hambuild_cuda_set_potential_bulk(hambuild_ctx *ctx, const void *wx0,
                                     const void *wx1, const void *cx0,
                                     const void *cx1, const void *cex0,
                                     const void *cex1, const double *mom,
                                     const double *q_ss, double theta_ss)
{
    (void)ctx; (void)wx0; (void)wx1; (void)cx0; (void)cx1; (void)cex0;
    (void)cex1; (void)mom; (void)q_ss; (void)theta_ss;
    return 1;
}

int hambuild_cuda_set_sbar(hambuild_ctx *ctx, const void *sbar, int nm_store,
                           int ntot)
{
    (void)ctx; (void)sbar; (void)nm_store; (void)ntot;
    return 1;
}

int hambuild_cuda_bulk(hambuild_ctx *ctx, int hoh, void *ee, void *hxc,
                       void *eeo, void *eeoee)
{
    (void)ctx; (void)hoh; (void)ee; (void)hxc; (void)eeo; (void)eeoee;
    return 1;
}

int hambuild_cuda_set_local_sites(hambuild_ctx *ctx, const int *site_list,
                                  int nmax)
{
    (void)ctx; (void)site_list; (void)nmax;
    return 1;
}

int hambuild_cuda_build_local_geometry_maps(hambuild_ctx *ctx)
{
    (void)ctx;
    return 1;
}

int hambuild_cuda_local(hambuild_ctx *ctx, int hoh, void *hall, void *hallo)
{
    (void)ctx; (void)hoh; (void)hall; (void)hallo;
    return 1;
}
