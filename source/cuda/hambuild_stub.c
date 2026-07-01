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
