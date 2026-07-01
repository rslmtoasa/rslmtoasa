/* ===========================================================================
 * hambuild_cuda.cpp -- backend dispatcher for GPU Hamiltonian assembly.
 *
 * Owns the opaque hambuild_ctx (dimensions + device pointers + backend
 * selector) and dispatches the numerical entry points to the CUDA/HIP kernels
 * in hambuild_gpu.cu. Phase 1 implements the on-site blocks (lsham/obarm/enim);
 * later phases add the neighbour hot loop, CCOR and Hubbard.
 *
 * Device buffers for the resident constant / potential / output arrays live in
 * the context so that (from Phase 2 on) assembled arrays can stay on the device
 * and feed recursion without a host round-trip. In Phase 1 the outputs are
 * copied back for element-wise validation against the CPU path.
 *
 * Mirrors the rsrec_cuda.cpp pattern (error string, create/destroy, backend
 * validation) so the two plugins stay structurally identical.
 * =========================================================================== */
#include "hambuild.h"

#include <cuComplex.h>
#include <cuda_runtime.h>

#include <string>

typedef cuDoubleComplex hbc;

/* Kernel launcher implemented in hambuild_gpu.cu. */
extern "C" int hambuild_gpu_launch_onsite(
    const hbc *d_v, const hbc *d_vc, const hbc *d_lx, const hbc *d_ly,
    const hbc *d_lz, const double *d_xi_p, const double *d_xi_d,
    const double *d_rac, const hbc *d_obx0, const hbc *d_obx1, const hbc *d_cx,
    const hbc *d_cex, const double *d_mom, const double *d_lmom, int orb_pol,
    hbc *d_lsham, hbc *d_obarm, hbc *d_enim, int norb, int nb, int ntype,
    int want_lsham);

struct hambuild_ctx {
    int nb = 0;
    int norb = 0;
    int nnmax = 0;
    int ntype = 0;
    int nmax = 0;
    int device = 0;
    int backend = HAMBUILD_BACKEND_GPU;

    /* --- Phase 1 resident device buffers --- */
    /* constants (basis-dependent, uploaded once) */
    hbc *d_v = nullptr, *d_vc = nullptr;
    hbc *d_lx = nullptr, *d_ly = nullptr, *d_lz = nullptr;
    bool have_constants = false;
    /* per-type potential inputs */
    double *d_xi_p = nullptr, *d_xi_d = nullptr, *d_rac = nullptr;
    double *d_mom = nullptr, *d_lmom = nullptr;
    hbc *d_obx0 = nullptr, *d_obx1 = nullptr, *d_cx = nullptr, *d_cex = nullptr;
    int orb_pol = 0;
    bool have_potential = false;
    /* outputs */
    hbc *d_lsham = nullptr, *d_obarm = nullptr, *d_enim = nullptr;
};

static std::string g_hambuild_err;
static void set_error(const std::string &msg) { g_hambuild_err = msg; }

extern "C" const char *hambuild_cuda_last_error(void) {
    return g_hambuild_err.c_str();
}

static bool dev_alloc(void **p, size_t bytes) {
    return cudaMalloc(p, bytes) == cudaSuccess;
}
static bool up(void *dst, const void *src, size_t bytes) {
    return cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice) == cudaSuccess;
}
static bool down(void *dst, const void *src, size_t bytes) {
    return cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
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

    const size_t nn = (size_t)norb * norb;
    const size_t blk = (size_t)nb * nb * ntype;
    bool ok = true;
    /* constants */
    ok = ok && dev_alloc((void **)&ctx->d_v, nn * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_vc, nn * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_lx, nn * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_ly, nn * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_lz, nn * sizeof(hbc));
    /* potential */
    ok = ok && dev_alloc((void **)&ctx->d_xi_p, 2 * (size_t)ntype * sizeof(double));
    ok = ok && dev_alloc((void **)&ctx->d_xi_d, 2 * (size_t)ntype * sizeof(double));
    ok = ok && dev_alloc((void **)&ctx->d_rac, 2 * (size_t)ntype * sizeof(double));
    ok = ok && dev_alloc((void **)&ctx->d_mom, 3 * (size_t)ntype * sizeof(double));
    ok = ok && dev_alloc((void **)&ctx->d_lmom, 3 * (size_t)ntype * sizeof(double));
    ok = ok && dev_alloc((void **)&ctx->d_obx0, (size_t)norb * ntype * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_obx1, (size_t)norb * ntype * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_cx, (size_t)norb * 2 * ntype * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_cex, (size_t)norb * 2 * ntype * sizeof(hbc));
    /* outputs */
    ok = ok && dev_alloc((void **)&ctx->d_lsham, blk * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_obarm, blk * sizeof(hbc));
    ok = ok && dev_alloc((void **)&ctx->d_enim, blk * sizeof(hbc));

    if (!ok) {
        set_error("hambuild_cuda_create: device allocation failed");
        hambuild_cuda_destroy(ctx);
        return nullptr;
    }
    return ctx;
}

extern "C" void hambuild_cuda_destroy(hambuild_ctx *ctx) {
    if (!ctx) return;
    cudaFree(ctx->d_v); cudaFree(ctx->d_vc);
    cudaFree(ctx->d_lx); cudaFree(ctx->d_ly); cudaFree(ctx->d_lz);
    cudaFree(ctx->d_xi_p); cudaFree(ctx->d_xi_d); cudaFree(ctx->d_rac);
    cudaFree(ctx->d_mom); cudaFree(ctx->d_lmom);
    cudaFree(ctx->d_obx0); cudaFree(ctx->d_obx1);
    cudaFree(ctx->d_cx); cudaFree(ctx->d_cex);
    cudaFree(ctx->d_lsham); cudaFree(ctx->d_obarm); cudaFree(ctx->d_enim);
    delete ctx;
}

extern "C" int hambuild_cuda_set_backend(hambuild_ctx *ctx, int backend) {
    if (!ctx) { set_error("hambuild_cuda_set_backend: null ctx"); return 1; }
    if (backend != HAMBUILD_BACKEND_CPU && backend != HAMBUILD_BACKEND_GPU) {
        set_error("hambuild_cuda_set_backend: invalid backend");
        return 1;
    }
    ctx->backend = backend;
    return 0;
}

extern "C" int hambuild_cuda_set_constants(hambuild_ctx *ctx, const void *v,
                                           const void *vc, const void *lx,
                                           const void *ly, const void *lz) {
    if (!ctx) { set_error("hambuild_cuda_set_constants: null ctx"); return 1; }
    const size_t nn = (size_t)ctx->norb * ctx->norb * sizeof(hbc);
    bool ok = up(ctx->d_v, v, nn) && up(ctx->d_vc, vc, nn) &&
              up(ctx->d_lx, lx, nn) && up(ctx->d_ly, ly, nn) &&
              up(ctx->d_lz, lz, nn);
    if (!ok) { set_error("hambuild_cuda_set_constants: upload failed"); return 1; }
    ctx->have_constants = true;
    return 0;
}

extern "C" int hambuild_cuda_set_potential_onsite(
    hambuild_ctx *ctx, const double *xi_p, const double *xi_d,
    const double *rac, const void *obx0, const void *obx1, const void *cx,
    const void *cex, const double *mom, const double *lmom, int orb_pol) {
    if (!ctx) { set_error("hambuild_cuda_set_potential_onsite: null ctx"); return 1; }
    const int nt = ctx->ntype, norb = ctx->norb;
    bool ok = true;
    ok = ok && up(ctx->d_xi_p, xi_p, 2 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_xi_d, xi_d, 2 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_rac, rac, 2 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_mom, mom, 3 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_lmom, lmom, 3 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_obx0, obx0, (size_t)norb * nt * sizeof(hbc));
    ok = ok && up(ctx->d_obx1, obx1, (size_t)norb * nt * sizeof(hbc));
    ok = ok && up(ctx->d_cx, cx, (size_t)norb * 2 * nt * sizeof(hbc));
    ok = ok && up(ctx->d_cex, cex, (size_t)norb * 2 * nt * sizeof(hbc));
    if (!ok) { set_error("hambuild_cuda_set_potential_onsite: upload failed"); return 1; }
    ctx->orb_pol = orb_pol;
    ctx->have_potential = true;
    return 0;
}

extern "C" int hambuild_cuda_onsite(hambuild_ctx *ctx, void *lsham, void *obarm,
                                    void *enim) {
    if (!ctx) { set_error("hambuild_cuda_onsite: null ctx"); return 1; }
    if (!ctx->have_constants || !ctx->have_potential) {
        set_error("hambuild_cuda_onsite: constants/potential not set");
        return 1;
    }
    const int want_lsham = (lsham != nullptr) ? 1 : 0;
    const int rc = hambuild_gpu_launch_onsite(
        ctx->d_v, ctx->d_vc, ctx->d_lx, ctx->d_ly, ctx->d_lz,
        ctx->d_xi_p, ctx->d_xi_d, ctx->d_rac, ctx->d_obx0, ctx->d_obx1,
        ctx->d_cx, ctx->d_cex, ctx->d_mom, ctx->d_lmom, ctx->orb_pol,
        lsham ? ctx->d_lsham : nullptr, obarm ? ctx->d_obarm : nullptr,
        enim ? ctx->d_enim : nullptr, ctx->norb, ctx->nb, ctx->ntype,
        want_lsham);
    if (rc != 0) {
        set_error("hambuild_cuda_onsite: kernel launch failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }
    const size_t blk = (size_t)ctx->nb * ctx->nb * ctx->ntype * sizeof(hbc);
    bool ok = true;
    if (lsham) ok = ok && down(lsham, ctx->d_lsham, blk);
    if (obarm) ok = ok && down(obarm, ctx->d_obarm, blk);
    if (enim)  ok = ok && down(enim, ctx->d_enim, blk);
    if (!ok) { set_error("hambuild_cuda_onsite: device->host copy failed"); return 1; }
    return 0;
}
