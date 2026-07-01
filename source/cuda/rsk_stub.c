#include "rsk.h"

struct rsk_ctx {
    int unused;
};

static const char *g_msg =
    "rsk: GPU k-space backend not built in this binary. "
    "Reconfigure with -DENABLE_CUDA_PLUGIN=ON when the rsk backend is implemented.";

const char *rsk_last_error(void) { return g_msg; }

rsk_ctx *rsk_create(int nb, int nsite, int ntype, int nnmax, int device) {
    (void)nb;
    (void)nsite;
    (void)ntype;
    (void)nnmax;
    (void)device;
    return 0;
}

void rsk_destroy(rsk_ctx *ctx) { (void)ctx; }

int rsk_set_geometry(rsk_ctx *ctx, const int *site_type, const int *atlist,
                     const int *nn, const int *iz, const double *r_direct) {
    (void)ctx;
    (void)site_type;
    (void)atlist;
    (void)nn;
    (void)iz;
    (void)r_direct;
    return 1;
}

int rsk_set_operators(rsk_ctx *ctx, const void *ee, const void *eeo,
                      const void *eecc, const void *enim,
                      const void *lsham) {
    (void)ctx;
    (void)ee;
    (void)eeo;
    (void)eecc;
    (void)enim;
    (void)lsham;
    return 1;
}

int rsk_build_hk_batch(rsk_ctx *ctx, const double *kpoints, int nk,
                       int order, void *hk_out) {
    (void)ctx;
    (void)kpoints;
    (void)nk;
    (void)order;
    (void)hk_out;
    return 1;
}

int rsk_diagonalize_hk_batch(rsk_ctx *ctx, int jobz, int generalized,
                             void *evals_out, void *evecs_out) {
    (void)ctx;
    (void)jobz;
    (void)generalized;
    (void)evals_out;
    (void)evecs_out;
    return 1;
}
