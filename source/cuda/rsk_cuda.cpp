#include "rsk_cuda.h"

#include <string>

struct rsk_cuda_ctx {
    int nb = 0;
    int nsite = 0;
    int ntype = 0;
    int nnmax = 0;
    int device = 0;
};

static std::string g_rsk_cuda_err;

static void set_error(const char *msg) {
    g_rsk_cuda_err = msg;
}

extern "C" const char *rsk_cuda_last_error(void) {
    return g_rsk_cuda_err.c_str();
}

extern "C" rsk_cuda_ctx *rsk_cuda_create(int nb, int nsite, int ntype,
                                         int nnmax, int device) {
    if (nb <= 0 || nsite <= 0 || ntype <= 0 || nnmax <= 0) {
        set_error("rsk_cuda_create: bad dimensions");
        return nullptr;
    }

    auto *ctx = new rsk_cuda_ctx();
    ctx->nb = nb;
    ctx->nsite = nsite;
    ctx->ntype = ntype;
    ctx->nnmax = nnmax;
    ctx->device = device;
    set_error("rsk_cuda: k-space backend scaffold only; kernels not implemented");
    return ctx;
}

extern "C" void rsk_cuda_destroy(rsk_cuda_ctx *ctx) {
    delete ctx;
}
