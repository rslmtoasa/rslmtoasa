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

/* Kernel launchers implemented in hambuild_gpu.cu. */
extern "C" int hambuild_gpu_launch_onsite(
    const hbc *d_v, const hbc *d_vc, const hbc *d_lx, const hbc *d_ly,
    const hbc *d_lz, const double *d_xi_p, const double *d_xi_d,
    const double *d_rac, const hbc *d_obx0, const hbc *d_obx1, const hbc *d_cx,
    const hbc *d_cex, const double *d_mom, const double *d_lmom, int orb_pol,
    hbc *d_lsham, hbc *d_obarm, hbc *d_enim, int norb, int nb, int ntype,
    int want_lsham);

extern "C" int hambuild_gpu_launch_geometry_maps(
    const double *d_cr, const int *d_num, const int *d_iz, const int *d_nn,
    const int *d_atlist, const double *d_vet_in, int use_vet_in, int kk,
    int nn_max, int ndi, double alat, double r2, int ntype, double *d_hv_scratch,
    int *d_valid, int *d_shell, int *d_ino, double *d_vet);

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

    /* --- Phase 1.5 geometry maps --- */
    int kk = 0, nn_max = 0, ndi = 0, pbc = 0;
    double alat = 0.0, r2 = 0.0;
    double *d_cr = nullptr;             /* 3*kk */
    int *d_num = nullptr, *d_iz = nullptr; /* kk */
    int *d_nn = nullptr;                /* ndi*nn_max */
    int *d_atlist = nullptr;            /* ntype */
    double *d_vet_in = nullptr;         /* 3*nn_max*ntype (pbc/precomputed) */
    bool have_vet_in = false;
    double *d_hv_scratch = nullptr;     /* 3*kk*ntype clusba scratch */
    /* resident maps: (nn_max, ntype) col-major */
    int *d_valid = nullptr, *d_shell = nullptr, *d_ino = nullptr;
    double *d_vet = nullptr;            /* 3*nn_max*ntype */
    bool have_geometry = false;
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
    cudaFree(ctx->d_cr); cudaFree(ctx->d_num); cudaFree(ctx->d_iz);
    cudaFree(ctx->d_nn); cudaFree(ctx->d_atlist); cudaFree(ctx->d_vet_in);
    cudaFree(ctx->d_hv_scratch);
    cudaFree(ctx->d_valid); cudaFree(ctx->d_shell); cudaFree(ctx->d_ino);
    cudaFree(ctx->d_vet);
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

/* --- Phase 1.5: geometry maps -------------------------------------------- */

extern "C" int hambuild_cuda_set_geometry(hambuild_ctx *ctx, const double *cr,
                                          const int *num, const int *iz,
                                          const int *nn, const int *atlist,
                                          int kk, int nn_max, int ndi,
                                          double alat, double r2, int pbc) {
    if (!ctx) { set_error("hambuild_cuda_set_geometry: null ctx"); return 1; }
    if (kk <= 0 || nn_max <= 0 || ndi <= 0) {
        set_error("hambuild_cuda_set_geometry: bad geometry dimensions");
        return 1;
    }
    const int nt = ctx->ntype;
    /* (Re)allocate geometry buffers if the dimensions changed. */
    bool realloc = (kk != ctx->kk || nn_max != ctx->nn_max || ndi != ctx->ndi ||
                    ctx->d_cr == nullptr);
    ctx->kk = kk; ctx->nn_max = nn_max; ctx->ndi = ndi;
    ctx->alat = alat; ctx->r2 = r2; ctx->pbc = pbc;
    bool ok = true;
    if (realloc) {
        cudaFree(ctx->d_cr); cudaFree(ctx->d_num); cudaFree(ctx->d_iz);
        cudaFree(ctx->d_nn); cudaFree(ctx->d_atlist); cudaFree(ctx->d_vet_in);
        cudaFree(ctx->d_hv_scratch);
        cudaFree(ctx->d_valid); cudaFree(ctx->d_shell); cudaFree(ctx->d_ino);
        cudaFree(ctx->d_vet);
        ctx->d_cr = ctx->d_vet_in = ctx->d_hv_scratch = ctx->d_vet = nullptr;
        ctx->d_num = ctx->d_iz = ctx->d_nn = ctx->d_atlist = nullptr;
        ctx->d_valid = ctx->d_shell = ctx->d_ino = nullptr;
        ok = ok && dev_alloc((void **)&ctx->d_cr, 3 * (size_t)kk * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_num, (size_t)kk * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_iz, (size_t)kk * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_nn, (size_t)ndi * nn_max * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_atlist, (size_t)nt * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_vet_in, 3 * (size_t)nn_max * nt * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_hv_scratch, 3 * (size_t)kk * nt * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_valid, (size_t)nn_max * nt * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_shell, (size_t)nn_max * nt * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_ino, (size_t)nn_max * nt * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_vet, 3 * (size_t)nn_max * nt * sizeof(double));
        if (!ok) { set_error("hambuild_cuda_set_geometry: allocation failed"); return 1; }
    }
    ok = ok && up(ctx->d_cr, cr, 3 * (size_t)kk * sizeof(double));
    ok = ok && up(ctx->d_num, num, (size_t)kk * sizeof(int));
    ok = ok && up(ctx->d_iz, iz, (size_t)kk * sizeof(int));
    ok = ok && up(ctx->d_nn, nn, (size_t)ndi * nn_max * sizeof(int));
    ok = ok && up(ctx->d_atlist, atlist, (size_t)nt * sizeof(int));
    if (!ok) { set_error("hambuild_cuda_set_geometry: upload failed"); return 1; }
    ctx->have_vet_in = false;
    ctx->have_geometry = false; /* maps must be rebuilt */
    return 0;
}

extern "C" int hambuild_cuda_set_geometry_vet(hambuild_ctx *ctx,
                                              const double *vet) {
    if (!ctx) { set_error("hambuild_cuda_set_geometry_vet: null ctx"); return 1; }
    if (!ctx->d_vet_in) {
        set_error("hambuild_cuda_set_geometry_vet: call set_geometry first");
        return 1;
    }
    const size_t bytes = 3 * (size_t)ctx->nn_max * ctx->ntype * sizeof(double);
    if (!up(ctx->d_vet_in, vet, bytes)) {
        set_error("hambuild_cuda_set_geometry_vet: upload failed");
        return 1;
    }
    ctx->have_vet_in = true;
    return 0;
}

extern "C" int hambuild_cuda_build_geometry_maps(hambuild_ctx *ctx) {
    if (!ctx) { set_error("hambuild_cuda_build_geometry_maps: null ctx"); return 1; }
    if (!ctx->d_cr) {
        set_error("hambuild_cuda_build_geometry_maps: geometry not set");
        return 1;
    }
    const int rc = hambuild_gpu_launch_geometry_maps(
        ctx->d_cr, ctx->d_num, ctx->d_iz, ctx->d_nn, ctx->d_atlist,
        ctx->d_vet_in, ctx->have_vet_in ? 1 : 0, ctx->kk, ctx->nn_max, ctx->ndi,
        ctx->alat, ctx->r2, ctx->ntype, ctx->d_hv_scratch, ctx->d_valid,
        ctx->d_shell, ctx->d_ino, ctx->d_vet);
    if (rc != 0) {
        set_error("hambuild_cuda_build_geometry_maps: kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }
    ctx->have_geometry = true;
    return 0;
}

extern "C" int hambuild_cuda_get_geometry_maps(hambuild_ctx *ctx, int *valid,
                                               int *shell, int *ino,
                                               double *vet) {
    if (!ctx) { set_error("hambuild_cuda_get_geometry_maps: null ctx"); return 1; }
    if (!ctx->have_geometry) {
        set_error("hambuild_cuda_get_geometry_maps: maps not built");
        return 1;
    }
    const size_t ni = (size_t)ctx->nn_max * ctx->ntype;
    bool ok = true;
    if (valid) ok = ok && down(valid, ctx->d_valid, ni * sizeof(int));
    if (shell) ok = ok && down(shell, ctx->d_shell, ni * sizeof(int));
    if (ino)   ok = ok && down(ino, ctx->d_ino, ni * sizeof(int));
    if (vet)   ok = ok && down(vet, ctx->d_vet, 3 * ni * sizeof(double));
    if (!ok) { set_error("hambuild_cuda_get_geometry_maps: copy failed"); return 1; }
    return 0;
}
