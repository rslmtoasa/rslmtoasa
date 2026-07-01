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
#include <cublas_v2.h>

#include <string>
#include <vector>

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

extern "C" int hambuild_gpu_launch_bulk_ham(
    const double *d_sbar_re, int norb_s, int nm_store, int ntot,
    const hbc *d_wx0, const hbc *d_wx1, const hbc *d_cx0, const hbc *d_cx1,
    const hbc *d_cex0, const hbc *d_cex1, const double *d_mom,
    const double *d_q_ss, double theta_ss, const hbc *d_v, const hbc *d_vc,
    const int *d_iz, const int *d_nn, const int *d_atlist, const double *d_cr,
    const int *d_valid, const int *d_shell, const int *d_ino, const double *d_vet,
    int ndi, int kk, int hoh, hbc *d_ee, hbc *d_hxc, int norb, int nb,
    int nn_max, int ntype);

extern "C" int hambuild_gpu_launch_ccor(
    const double *d_sbar_re, const double *d_sdot_re, int norb_s, int nm_store,
    int ntot, const hbc *d_wx0, const hbc *d_wx1, const double *d_mom,
    const double *d_q_ss, double theta_ss, const double *d_ccd,
    const double *d_lambda, double avw, const hbc *d_v, const hbc *d_vc,
    const int *d_iz, const int *d_nn, const int *d_atlist, const double *d_cr,
    const int *d_shell, const int *d_ino, const double *d_vet, int ndi, int kk,
    hbc *d_out, int norb, int nb, int nn_max, int nsite, int lambda_stride);

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

    /* --- Phase 2 bulk/local hot loop --- */
    cublasHandle_t cublas = nullptr;
    hbc *d_wx0 = nullptr, *d_wx1 = nullptr, *d_cx0 = nullptr, *d_cx1 = nullptr;
    hbc *d_cex0 = nullptr, *d_cex1 = nullptr;
    double *d_mom_bulk = nullptr;       /* 3*ntype */
    double *d_q_ss = nullptr;           /* 3 */
    double theta_ss = 0.0;
    bool have_potential_bulk = false;
    double *d_sbar_re = nullptr;        /* norb*norb*nm_store*ntot */
    int nm_store = 0, ntot = 0;
    bool have_sbar = false;
    hbc *d_ee = nullptr, *d_hxc = nullptr, *d_eeo = nullptr, *d_eeoee = nullptr;
    /* host-side geometry mirror to compute ji for hoh batching */
    std::vector<int> h_iz, h_nn, h_atlist, h_valid;
    int h_ndi = 0;

    /* --- Phase 2b local (impurity) path (nmax reused from create) --- */
    int *d_local_atlist = nullptr;      /* nmax (1-based cluster indices) */
    int *d_local_valid = nullptr, *d_local_shell = nullptr, *d_local_ino = nullptr;
    double *d_local_vet = nullptr;      /* 3*nn_max*nmax */
    double *d_local_hv = nullptr;       /* 3*kk*nmax clusba scratch */
    hbc *d_hall = nullptr, *d_hallo = nullptr, *d_hall_hxc = nullptr; /* hxc scratch */
    std::vector<int> h_local_atlist, h_local_valid;
    bool have_local_geometry = false;

    /* --- Phase 3 CCOR --- */
    double *d_ccd = nullptr;            /* norb*3*ntype */
    double *d_lambda = nullptr;         /* ntype*ntype */
    double *d_sdot_re = nullptr;        /* norb*norb*nm_store*ntot */
    double ccor_avw = 0.0;
    bool have_ccor = false;
    hbc *d_eecc = nullptr;              /* nb*nb*nn_max*ntype */
    hbc *d_hallcc = nullptr;           /* nb*nb*nn_max*nmax */
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
    if (cublasCreate(&ctx->cublas) != CUBLAS_STATUS_SUCCESS) {
        set_error("hambuild_cuda_create: cublasCreate failed");
        delete ctx;
        return nullptr;
    }
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
    cudaFree(ctx->d_wx0); cudaFree(ctx->d_wx1); cudaFree(ctx->d_cx0);
    cudaFree(ctx->d_cx1); cudaFree(ctx->d_cex0); cudaFree(ctx->d_cex1);
    cudaFree(ctx->d_mom_bulk); cudaFree(ctx->d_q_ss); cudaFree(ctx->d_sbar_re);
    cudaFree(ctx->d_ee); cudaFree(ctx->d_hxc); cudaFree(ctx->d_eeo);
    cudaFree(ctx->d_eeoee);
    cudaFree(ctx->d_local_atlist); cudaFree(ctx->d_local_valid);
    cudaFree(ctx->d_local_shell); cudaFree(ctx->d_local_ino);
    cudaFree(ctx->d_local_vet); cudaFree(ctx->d_local_hv);
    cudaFree(ctx->d_hall); cudaFree(ctx->d_hallo); cudaFree(ctx->d_hall_hxc);
    cudaFree(ctx->d_ccd); cudaFree(ctx->d_lambda); cudaFree(ctx->d_sdot_re);
    cudaFree(ctx->d_eecc); cudaFree(ctx->d_hallcc);
    if (ctx->cublas) cublasDestroy(ctx->cublas);
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
    /* Host mirror of iz/nn/atlist for hoh ji computation (see hambuild_cuda_bulk). */
    ctx->h_iz.assign(iz, iz + kk);
    ctx->h_nn.assign(nn, nn + (size_t)ndi * nn_max);
    ctx->h_atlist.assign(atlist, atlist + nt);
    ctx->h_ndi = ndi;
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
    /* Mirror valid to host for hoh ji computation. */
    ctx->h_valid.resize((size_t)ctx->nn_max * ctx->ntype);
    if (!down(ctx->h_valid.data(), ctx->d_valid,
              ctx->h_valid.size() * sizeof(int))) {
        set_error("hambuild_cuda_build_geometry_maps: valid mirror copy failed");
        return 1;
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

/* --- Phase 2: bulk/local hot loop ---------------------------------------- */

extern "C" int hambuild_cuda_set_potential_bulk(
    hambuild_ctx *ctx, const void *wx0, const void *wx1, const void *cx0,
    const void *cx1, const void *cex0, const void *cex1, const double *mom,
    const double *q_ss, double theta_ss) {
    if (!ctx) { set_error("hambuild_cuda_set_potential_bulk: null ctx"); return 1; }
    const int nt = ctx->ntype, norb = ctx->norb;
    const size_t vbytes = (size_t)norb * nt * sizeof(hbc);
    bool ok = true;
    if (!ctx->d_wx0) {
        ok = ok && dev_alloc((void **)&ctx->d_wx0, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_wx1, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_cx0, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_cx1, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_cex0, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_cex1, vbytes);
        ok = ok && dev_alloc((void **)&ctx->d_mom_bulk, 3 * (size_t)nt * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_q_ss, 3 * sizeof(double));
        if (!ok) { set_error("hambuild_cuda_set_potential_bulk: alloc failed"); return 1; }
    }
    ok = ok && up(ctx->d_wx0, wx0, vbytes);
    ok = ok && up(ctx->d_wx1, wx1, vbytes);
    ok = ok && up(ctx->d_cx0, cx0, vbytes);
    ok = ok && up(ctx->d_cx1, cx1, vbytes);
    ok = ok && up(ctx->d_cex0, cex0, vbytes);
    ok = ok && up(ctx->d_cex1, cex1, vbytes);
    ok = ok && up(ctx->d_mom_bulk, mom, 3 * (size_t)nt * sizeof(double));
    ok = ok && up(ctx->d_q_ss, q_ss, 3 * sizeof(double));
    if (!ok) { set_error("hambuild_cuda_set_potential_bulk: upload failed"); return 1; }
    ctx->theta_ss = theta_ss;
    ctx->have_potential_bulk = true;
    return 0;
}

extern "C" int hambuild_cuda_set_sbar(hambuild_ctx *ctx, const void *sbar,
                                      int nm_store, int ntot) {
    if (!ctx) { set_error("hambuild_cuda_set_sbar: null ctx"); return 1; }
    const int norb = ctx->norb;
    const size_t n = (size_t)norb * norb * nm_store * ntot;
    /* Extract the real part on the host (CPU uses real(sbar)); upload as double. */
    const double *src = static_cast<const double *>(sbar); /* interleaved re,im */
    std::vector<double> re(n);
    for (size_t i = 0; i < n; ++i) re[i] = src[2 * i];
    if (nm_store != ctx->nm_store || ntot != ctx->ntot || !ctx->d_sbar_re) {
        cudaFree(ctx->d_sbar_re);
        ctx->d_sbar_re = nullptr;
        if (!dev_alloc((void **)&ctx->d_sbar_re, n * sizeof(double))) {
            set_error("hambuild_cuda_set_sbar: alloc failed");
            return 1;
        }
        ctx->nm_store = nm_store;
        ctx->ntot = ntot;
    }
    if (!up(ctx->d_sbar_re, re.data(), n * sizeof(double))) {
        set_error("hambuild_cuda_set_sbar: upload failed");
        return 1;
    }
    ctx->have_sbar = true;
    return 0;
}

extern "C" int hambuild_cuda_bulk(hambuild_ctx *ctx, int hoh, void *ee,
                                  void *hxc, void *eeo, void *eeoee) {
    if (!ctx) { set_error("hambuild_cuda_bulk: null ctx"); return 1; }
    if (!ctx->have_geometry || !ctx->have_potential_bulk || !ctx->have_sbar) {
        set_error("hambuild_cuda_bulk: geometry/potential/sbar not set");
        return 1;
    }
    const int nb = ctx->nb, norb = ctx->norb, nt = ctx->ntype, nm = ctx->nn_max;
    const size_t blk = (size_t)nb * nb * nm * nt;
    bool ok = true;
    if (!ctx->d_ee) {
        ok = ok && dev_alloc((void **)&ctx->d_ee, blk * sizeof(hbc));
        ok = ok && dev_alloc((void **)&ctx->d_hxc, blk * sizeof(hbc));
        ok = ok && dev_alloc((void **)&ctx->d_eeo, blk * sizeof(hbc));
        ok = ok && dev_alloc((void **)&ctx->d_eeoee, blk * sizeof(hbc));
        if (!ok) { set_error("hambuild_cuda_bulk: output alloc failed"); return 1; }
    }

    int rc = hambuild_gpu_launch_bulk_ham(
        ctx->d_sbar_re, norb, ctx->nm_store, ctx->ntot, ctx->d_wx0, ctx->d_wx1,
        ctx->d_cx0, ctx->d_cx1, ctx->d_cex0, ctx->d_cex1, ctx->d_mom_bulk,
        ctx->d_q_ss, ctx->theta_ss, ctx->d_v, ctx->d_vc, ctx->d_iz, ctx->d_nn,
        ctx->d_atlist, ctx->d_cr, ctx->d_valid, ctx->d_shell, ctx->d_ino,
        ctx->d_vet, ctx->h_ndi, ctx->kk, hoh, ctx->d_ee, ctx->d_hxc, norb, nb,
        nm, nt);
    if (rc != 0) {
        set_error("hambuild_cuda_bulk: ham kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }

    if (hoh) {
        /* eeo(:,:,m,t) = ee(:,:,m,t) * obarm(:,:,ji);
         * eeoee(:,:,m,t) = eeo(:,:,m,t) * ee(:,:,m,t)^H.
         * ji = iz(jj) with jj = nn(ia,m) (m>0) or ia (m==0); absent -> zero eeo.
         * ji varies per (t,m), so use cublasZgemmBatched with pointer arrays. */
        cudaMemset(ctx->d_eeo, 0, blk * sizeof(hbc));
        cudaMemset(ctx->d_eeoee, 0, blk * sizeof(hbc));

        const int ndi = ctx->h_ndi;
        std::vector<const hbc *> Aptr, Bptr; std::vector<hbc *> Cptr;
        std::vector<const hbc *> A2, B2; std::vector<hbc *> C2;
        const size_t slab = (size_t)nb * nb;
        for (int t = 0; t < nt; ++t) {
            int ia = ctx->h_atlist[t] - 1;
            int nr = ctx->h_nn[ia + (size_t)ndi * 0];
            for (int m = 0; m < nr; ++m) {
                int jj = (m == 0) ? (ia + 1) : ctx->h_nn[ia + (size_t)ndi * m];
                int vld = ctx->h_valid[m + (size_t)nm * t];
                int ji = 0;
                if (jj != 0 && vld != 0) ji = ctx->h_iz[jj - 1];
                if (ji <= 0) continue; /* leave eeo/eeoee zero for this (t,m) */
                size_t off = slab * (m + (size_t)nm * t);
                Aptr.push_back(ctx->d_ee + off);
                Bptr.push_back(ctx->d_obarm + slab * (ji - 1));
                Cptr.push_back(ctx->d_eeo + off);
            }
        }
        int batch = (int)Aptr.size();
        if (batch > 0) {
            const hbc **dA, **dB; hbc **dC;
            cudaMalloc((void **)&dA, batch * sizeof(hbc *));
            cudaMalloc((void **)&dB, batch * sizeof(hbc *));
            cudaMalloc((void **)&dC, batch * sizeof(hbc *));
            cudaMemcpy(dA, Aptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            cudaMemcpy(dB, Bptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            cudaMemcpy(dC, Cptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            hbc one = make_cuDoubleComplex(1.0, 0.0);
            hbc zero = make_cuDoubleComplex(0.0, 0.0);
            /* eeo = ee * obarm  (n,n) */
            cublasZgemmBatched(ctx->cublas, CUBLAS_OP_N, CUBLAS_OP_N, nb, nb, nb,
                               &one, dA, nb, dB, nb, &zero, dC, nb, batch);
            /* eeoee = eeo * ee^H  ('n','c'); reuse eeo(dC) as A2, ee(dA) as B2 */
            std::vector<hbc *> C2p;
            for (int t = 0; t < nt; ++t) {
                int ia = ctx->h_atlist[t] - 1;
                int nr = ctx->h_nn[ia + (size_t)ndi * 0];
                for (int m = 0; m < nr; ++m) {
                    int jj = (m == 0) ? (ia + 1) : ctx->h_nn[ia + (size_t)ndi * m];
                    int vld = ctx->h_valid[m + (size_t)nm * t];
                    int ji = 0;
                    if (jj != 0 && vld != 0) ji = ctx->h_iz[jj - 1];
                    if (ji <= 0) continue;
                    size_t off = slab * (m + (size_t)nm * t);
                    C2p.push_back(ctx->d_eeoee + off);
                }
            }
            hbc **dC2;
            cudaMalloc((void **)&dC2, batch * sizeof(hbc *));
            cudaMemcpy(dC2, C2p.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            cublasZgemmBatched(ctx->cublas, CUBLAS_OP_N, CUBLAS_OP_C, nb, nb, nb,
                               &one, (const hbc **)dC, nb, dA, nb, &zero, dC2, nb,
                               batch);
            cudaDeviceSynchronize();
            cudaFree(dA); cudaFree(dB); cudaFree(dC); cudaFree(dC2);
        }
    }

    /* Copy outputs back to host. */
    ok = true;
    if (ee)  ok = ok && down(ee, ctx->d_ee, blk * sizeof(hbc));
    if (hxc) ok = ok && down(hxc, ctx->d_hxc, blk * sizeof(hbc));
    if (hoh && eeo)   ok = ok && down(eeo, ctx->d_eeo, blk * sizeof(hbc));
    if (hoh && eeoee) ok = ok && down(eeoee, ctx->d_eeoee, blk * sizeof(hbc));
    if (!ok) { set_error("hambuild_cuda_bulk: copy back failed"); return 1; }
    return 0;
}

/* --- Phase 2b: local (impurity interaction-zone) Hamiltonian ------------- */

extern "C" int hambuild_cuda_set_local_sites(hambuild_ctx *ctx,
                                             const int *site_list, int nmax) {
    if (!ctx) { set_error("hambuild_cuda_set_local_sites: null ctx"); return 1; }
    if (nmax <= 0) { set_error("hambuild_cuda_set_local_sites: nmax<=0"); return 1; }
    if (!ctx->d_cr) {
        set_error("hambuild_cuda_set_local_sites: call set_geometry first");
        return 1;
    }
    const int nm = ctx->nn_max, kk = ctx->kk;
    bool ok = true;
    if (nmax != ctx->nmax || !ctx->d_local_atlist) {
        cudaFree(ctx->d_local_atlist); cudaFree(ctx->d_local_valid);
        cudaFree(ctx->d_local_shell); cudaFree(ctx->d_local_ino);
        cudaFree(ctx->d_local_vet); cudaFree(ctx->d_local_hv);
        cudaFree(ctx->d_hall); cudaFree(ctx->d_hallo); cudaFree(ctx->d_hall_hxc);
        ctx->d_local_atlist = ctx->d_local_valid = ctx->d_local_shell = nullptr;
        ctx->d_local_ino = nullptr;
        ctx->d_local_vet = ctx->d_local_hv = nullptr;
        ctx->d_hall = ctx->d_hallo = ctx->d_hall_hxc = nullptr;
        ctx->nmax = nmax;
        ok = ok && dev_alloc((void **)&ctx->d_local_atlist, (size_t)nmax * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_local_valid, (size_t)nm * nmax * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_local_shell, (size_t)nm * nmax * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_local_ino, (size_t)nm * nmax * sizeof(int));
        ok = ok && dev_alloc((void **)&ctx->d_local_vet, 3 * (size_t)nm * nmax * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_local_hv, 3 * (size_t)kk * nmax * sizeof(double));
        const size_t lblk = (size_t)ctx->nb * ctx->nb * nm * nmax;
        ok = ok && dev_alloc((void **)&ctx->d_hall, lblk * sizeof(hbc));
        ok = ok && dev_alloc((void **)&ctx->d_hallo, lblk * sizeof(hbc));
        ok = ok && dev_alloc((void **)&ctx->d_hall_hxc, lblk * sizeof(hbc));
        if (!ok) { set_error("hambuild_cuda_set_local_sites: alloc failed"); return 1; }
    }
    if (!up(ctx->d_local_atlist, site_list, (size_t)nmax * sizeof(int))) {
        set_error("hambuild_cuda_set_local_sites: upload failed");
        return 1;
    }
    ctx->h_local_atlist.assign(site_list, site_list + nmax);
    ctx->have_local_geometry = false;
    return 0;
}

extern "C" int hambuild_cuda_build_local_geometry_maps(hambuild_ctx *ctx) {
    if (!ctx) { set_error("hambuild_cuda_build_local_geometry_maps: null ctx"); return 1; }
    if (!ctx->d_local_atlist) {
        set_error("hambuild_cuda_build_local_geometry_maps: local sites not set");
        return 1;
    }
    /* Reuse the geometry-map kernel with the local site list as the "atlist" and
     * nmax as the "type count". use_vet_in=0 (non-pbc); pbc impurity would need a
     * local vet upload -- deferred (impurity examples are non-pbc). */
    const int rc = hambuild_gpu_launch_geometry_maps(
        ctx->d_cr, ctx->d_num, ctx->d_iz, ctx->d_nn, ctx->d_local_atlist,
        nullptr, 0, ctx->kk, ctx->nn_max, ctx->h_ndi, ctx->alat, ctx->r2,
        ctx->nmax, ctx->d_local_hv, ctx->d_local_valid, ctx->d_local_shell,
        ctx->d_local_ino, ctx->d_local_vet);
    if (rc != 0) {
        set_error("hambuild_cuda_build_local_geometry_maps: kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }
    ctx->h_local_valid.resize((size_t)ctx->nn_max * ctx->nmax);
    if (!down(ctx->h_local_valid.data(), ctx->d_local_valid,
              ctx->h_local_valid.size() * sizeof(int))) {
        set_error("hambuild_cuda_build_local_geometry_maps: valid mirror failed");
        return 1;
    }
    ctx->have_local_geometry = true;
    return 0;
}

extern "C" int hambuild_cuda_local(hambuild_ctx *ctx, int hoh, void *hall,
                                   void *hallo) {
    if (!ctx) { set_error("hambuild_cuda_local: null ctx"); return 1; }
    if (!ctx->have_local_geometry || !ctx->have_potential_bulk || !ctx->have_sbar) {
        set_error("hambuild_cuda_local: local geometry/potential/sbar not set");
        return 1;
    }
    const int nb = ctx->nb, norb = ctx->norb, nm = ctx->nn_max, nmax = ctx->nmax;
    const size_t lblk = (size_t)nb * nb * nm * nmax;

    /* Reuse the bulk ham kernel over the local site list. hall <- ee slot;
     * hxc goes to scratch (local path has no hxc output). hoh flag drives the
     * on-site cx/cex selection identically. */
    int rc = hambuild_gpu_launch_bulk_ham(
        ctx->d_sbar_re, norb, ctx->nm_store, ctx->ntot, ctx->d_wx0, ctx->d_wx1,
        ctx->d_cx0, ctx->d_cx1, ctx->d_cex0, ctx->d_cex1, ctx->d_mom_bulk,
        ctx->d_q_ss, ctx->theta_ss, ctx->d_v, ctx->d_vc, ctx->d_iz, ctx->d_nn,
        ctx->d_local_atlist, ctx->d_cr, ctx->d_local_valid, ctx->d_local_shell,
        ctx->d_local_ino, ctx->d_local_vet, ctx->h_ndi, ctx->kk, hoh,
        ctx->d_hall, ctx->d_hall_hxc, norb, nb, nm, nmax);
    if (rc != 0) {
        set_error("hambuild_cuda_local: ham kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }

    if (hoh) {
        /* hallo(:,:,m,nlim) = hall(:,:,m,nlim) * obarm(:,:,ji);  single gemm.
         * ji = iz(jj), jj = nn(nlim,m) (m>0) or nlim (m==0); absent -> zero. */
        cudaMemset(ctx->d_hallo, 0, lblk * sizeof(hbc));
        const int ndi = ctx->h_ndi;
        std::vector<const hbc *> Aptr, Bptr; std::vector<hbc *> Cptr;
        const size_t slab = (size_t)nb * nb;
        for (int s = 0; s < nmax; ++s) {
            int ia = ctx->h_local_atlist[s] - 1;
            int nr = ctx->h_nn[ia + (size_t)ndi * 0];
            for (int m = 0; m < nr; ++m) {
                int jj = (m == 0) ? (ia + 1) : ctx->h_nn[ia + (size_t)ndi * m];
                int vld = ctx->h_local_valid[m + (size_t)nm * s];
                int ji = 0;
                if (jj != 0 && vld != 0) ji = ctx->h_iz[jj - 1];
                if (ji <= 0) continue;
                size_t off = slab * (m + (size_t)nm * s);
                Aptr.push_back(ctx->d_hall + off);
                Bptr.push_back(ctx->d_obarm + slab * (ji - 1));
                Cptr.push_back(ctx->d_hallo + off);
            }
        }
        int batch = (int)Aptr.size();
        if (batch > 0) {
            const hbc **dA, **dB; hbc **dC;
            cudaMalloc((void **)&dA, batch * sizeof(hbc *));
            cudaMalloc((void **)&dB, batch * sizeof(hbc *));
            cudaMalloc((void **)&dC, batch * sizeof(hbc *));
            cudaMemcpy(dA, Aptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            cudaMemcpy(dB, Bptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            cudaMemcpy(dC, Cptr.data(), batch * sizeof(hbc *), cudaMemcpyHostToDevice);
            hbc one = make_cuDoubleComplex(1.0, 0.0);
            hbc zero = make_cuDoubleComplex(0.0, 0.0);
            cublasZgemmBatched(ctx->cublas, CUBLAS_OP_N, CUBLAS_OP_N, nb, nb, nb,
                               &one, dA, nb, dB, nb, &zero, dC, nb, batch);
            cudaDeviceSynchronize();
            cudaFree(dA); cudaFree(dB); cudaFree(dC);
        }
    }

    bool ok = true;
    if (hall)  ok = ok && down(hall, ctx->d_hall, lblk * sizeof(hbc));
    if (hoh && hallo) ok = ok && down(hallo, ctx->d_hallo, lblk * sizeof(hbc));
    if (!ok) { set_error("hambuild_cuda_local: copy back failed"); return 1; }
    return 0;
}

/* --- Phase 3: CCOR (H_cc) ------------------------------------------------- */

extern "C" int hambuild_cuda_set_ccor(hambuild_ctx *ctx, const double *ccd,
                                      const double *lambda, const void *sdot,
                                      double avw) {
    if (!ctx) { set_error("hambuild_cuda_set_ccor: null ctx"); return 1; }
    if (!ctx->have_sbar) {
        set_error("hambuild_cuda_set_ccor: set_sbar first (nm_store/ntot needed)");
        return 1;
    }
    const int norb = ctx->norb, nt = ctx->ntype;
    const size_t nccd = (size_t)norb * 3 * ctx->kk; /* per cluster site */
    const size_t nlam = (size_t)nt * nt;
    const size_t nsd = (size_t)norb * norb * ctx->nm_store * ctx->ntot;
    bool ok = true;
    if (!ctx->d_ccd) {
        ok = ok && dev_alloc((void **)&ctx->d_ccd, nccd * sizeof(double));
        ok = ok && dev_alloc((void **)&ctx->d_lambda, nlam * sizeof(double));
        if (!ok) { set_error("hambuild_cuda_set_ccor: alloc failed"); return 1; }
    }
    if (!ctx->d_sdot_re) {
        if (!dev_alloc((void **)&ctx->d_sdot_re, nsd * sizeof(double))) {
            set_error("hambuild_cuda_set_ccor: sdot alloc failed"); return 1;
        }
    }
    /* Extract real part of sdot (CPU uses real(sdot)). */
    const double *src = static_cast<const double *>(sdot); /* interleaved re,im */
    std::vector<double> re(nsd);
    for (size_t i = 0; i < nsd; ++i) re[i] = src[2 * i];
    ok = ok && up(ctx->d_ccd, ccd, nccd * sizeof(double));
    ok = ok && up(ctx->d_lambda, lambda, nlam * sizeof(double));
    ok = ok && up(ctx->d_sdot_re, re.data(), nsd * sizeof(double));
    if (!ok) { set_error("hambuild_cuda_set_ccor: upload failed"); return 1; }
    ctx->ccor_avw = avw;
    ctx->have_ccor = true;
    return 0;
}

extern "C" int hambuild_cuda_ccor_bulk(hambuild_ctx *ctx, void *eecc) {
    if (!ctx) { set_error("hambuild_cuda_ccor_bulk: null ctx"); return 1; }
    if (!ctx->have_geometry || !ctx->have_potential_bulk || !ctx->have_sbar ||
        !ctx->have_ccor) {
        set_error("hambuild_cuda_ccor_bulk: geometry/potential/sbar/ccor not set");
        return 1;
    }
    const int nb = ctx->nb, norb = ctx->norb, nt = ctx->ntype, nm = ctx->nn_max;
    const size_t blk = (size_t)nb * nb * nm * nt;
    if (!ctx->d_eecc) {
        if (!dev_alloc((void **)&ctx->d_eecc, blk * sizeof(hbc))) {
            set_error("hambuild_cuda_ccor_bulk: alloc failed"); return 1;
        }
    }
    int rc = hambuild_gpu_launch_ccor(
        ctx->d_sbar_re, ctx->d_sdot_re, norb, ctx->nm_store, ctx->ntot,
        ctx->d_wx0, ctx->d_wx1, ctx->d_mom_bulk, ctx->d_q_ss, ctx->theta_ss,
        ctx->d_ccd, ctx->d_lambda, ctx->ccor_avw, ctx->d_v, ctx->d_vc, ctx->d_iz,
        ctx->d_nn, ctx->d_atlist, ctx->d_cr, ctx->d_shell, ctx->d_ino,
        ctx->d_vet, ctx->h_ndi, ctx->kk, ctx->d_eecc, norb, nb, nm, nt, nt);
    if (rc != 0) {
        set_error("hambuild_cuda_ccor_bulk: kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }
    if (eecc && !down(eecc, ctx->d_eecc, blk * sizeof(hbc))) {
        set_error("hambuild_cuda_ccor_bulk: copy back failed"); return 1;
    }
    return 0;
}

extern "C" int hambuild_cuda_ccor_local(hambuild_ctx *ctx, void *hallcc) {
    if (!ctx) { set_error("hambuild_cuda_ccor_local: null ctx"); return 1; }
    if (!ctx->have_local_geometry || !ctx->have_potential_bulk ||
        !ctx->have_sbar || !ctx->have_ccor) {
        set_error("hambuild_cuda_ccor_local: local geometry/potential/sbar/ccor not set");
        return 1;
    }
    const int nb = ctx->nb, norb = ctx->norb, nm = ctx->nn_max, nmax = ctx->nmax;
    const size_t lblk = (size_t)nb * nb * nm * nmax;
    if (!ctx->d_hallcc) {
        if (!dev_alloc((void **)&ctx->d_hallcc, lblk * sizeof(hbc))) {
            set_error("hambuild_cuda_ccor_local: alloc failed"); return 1;
        }
    }
    int rc = hambuild_gpu_launch_ccor(
        ctx->d_sbar_re, ctx->d_sdot_re, norb, ctx->nm_store, ctx->ntot,
        ctx->d_wx0, ctx->d_wx1, ctx->d_mom_bulk, ctx->d_q_ss, ctx->theta_ss,
        ctx->d_ccd, ctx->d_lambda, ctx->ccor_avw, ctx->d_v, ctx->d_vc, ctx->d_iz,
        ctx->d_nn, ctx->d_local_atlist, ctx->d_cr, ctx->d_local_shell,
        ctx->d_local_ino, ctx->d_local_vet, ctx->h_ndi, ctx->kk, ctx->d_hallcc,
        norb, nb, nm, nmax, ctx->ntype);
    if (rc != 0) {
        set_error("hambuild_cuda_ccor_local: kernel failed (code " +
                  std::to_string(rc) + ")");
        return rc;
    }
    if (hallcc && !down(hallcc, ctx->d_hallcc, lblk * sizeof(hbc))) {
        set_error("hambuild_cuda_ccor_local: copy back failed"); return 1;
    }
    return 0;
}
