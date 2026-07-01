/* ===========================================================================
 * hambuild_cuda.cpp -- backend dispatcher for GPU Hamiltonian assembly.
 *
 * Owns the opaque hambuild_ctx (dimensions + device pointers + backend
 * selector) and dispatches the numerical entry points to either the CPU
 * reference path or the CUDA/HIP kernels in hambuild_gpu.cu. Phase 0 only wires
 * the lifecycle and backend selection; the set_* / onsite / bulk / local / ccor
 * are added phase by phase.
 *
 * Mirrors the rsrec_cuda.cpp pattern (error string, create/destroy, backend
 * validation) so the two plugins stay structurally identical.
 * =========================================================================== */
#include "hambuild.h"

#include <cuda_runtime_api.h>

#include <string>

struct hambuild_ctx {
    int nb = 0;
    int norb = 0;
    int nnmax = 0;
    int ntype = 0;
    int nmax = 0;
    int device = 0;
    int backend = HAMBUILD_BACKEND_GPU;
    /* Device pointers for resident Hamiltonian / potential / geometry arrays
     * are added in Phase 1+. Keeping them here (not in Fortran) is what lets
     * the assembled arrays stay resident and feed recursion without a host
     * round-trip. */
};

static std::string g_hambuild_err;
static void set_error(const std::string &msg) { g_hambuild_err = msg; }

extern "C" const char *hambuild_cuda_last_error(void) {
    return g_hambuild_err.c_str();
}

extern "C" hambuild_ctx *hambuild_cuda_create(int nb, int norb, int nnmax,
                                              int ntype, int nmax, int device) {
    if (nb <= 0 || norb <= 0 || nnmax <= 0 || ntype <= 0 || nmax < 0) {
        set_error("hambuild_cuda_create: bad dimensions");
        return nullptr;
    }

    if (cudaSetDevice(device) != cudaSuccess) {
        set_error("hambuild_cuda_create: cudaSetDevice failed");
        return nullptr;
    }

    auto *ctx = new hambuild_ctx();
    ctx->nb = nb;
    ctx->norb = norb;
    ctx->nnmax = nnmax;
    ctx->ntype = ntype;
    ctx->nmax = nmax;
    ctx->device = device;
    ctx->backend = HAMBUILD_BACKEND_GPU;
    return ctx;
}

extern "C" void hambuild_cuda_destroy(hambuild_ctx *ctx) {
    if (!ctx) return;
    delete ctx;
}

extern "C" int hambuild_cuda_set_backend(hambuild_ctx *ctx, int backend) {
    if (!ctx) {
        set_error("hambuild_cuda_set_backend: null ctx");
        return 1;
    }
    if (backend != HAMBUILD_BACKEND_CPU && backend != HAMBUILD_BACKEND_GPU) {
        set_error("hambuild_cuda_set_backend: invalid backend");
        return 1;
    }
    ctx->backend = backend;
    return 0;
}
