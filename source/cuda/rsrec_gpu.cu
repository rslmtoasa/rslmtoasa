/* ===========================================================================
 * rsrec_gpu.cu -- CUDA implementation of rsrec.h (v3)
 *
 * Strategy (general GPUs; HPC parts A100/H100 run FP64 at 1/2 FP32,
 * GeForce/RTX-A at 1/64 -- both are served):
 *
 *   matvec     -> ONE fused templated kernel for BOTH precisions:
 *                 template<CT, NB, TC>, NB in {2, 8, 18, 32}
 *                 (s/sp/spd/spdf x 2 spinors), dynamic fallback otherwise.
 *                 Padded shared memory (LDS = NB+1 for even NB: conflict-
 *                 free column-strided reads), fully unrolled inner product,
 *                 register tiling TC = 2 columns/thread for NB >= 16,
 *                 column slabs keep dynamic shared <= 40 KB at any nrhs.
 *                 Chebyshev recurrence and (y - b x)/a shift fused into the
 *                 epilogue: zero extra field passes.
 *                 The per-shell cublasZgemmBatched pointer-array matvec is
 *                 REMOVED: 10-30% tiny-GEMM efficiency + per-shell output
 *                 re-read made it dominated on every architecture.
 *   layout     -> ALL drivers use site-major psi: (ld = nb*kk) x nrhs
 *                 column-major. Every block reduction sum_k A_k^H B_k is
 *                 ONE tall-skinny GEMM (Zgemm / Cgemm3m) and every right-
 *                 multiply psi*M is ONE (ld x nb x nb) GEMM with pointer
 *                 swap (no copy-back, no batched scratch, no gemv, no
 *                 atomics in hot paths).
 *   batching   -> rsrec_chebyshev_moments_batch runs several starting
 *                 states as extra RHS columns of the same recurrence:
 *                 H-block loads amortize over nstates (cheapest remaining
 *                 speedup; directly serves the J_ij exchange driver).
 *   precision  -> fp32 Chebyshev engine by default (KPM-safe: RMP 78, 275;
 *                 KITE/GPUQT practice), fp64 via rsrec_set_precision(1).
 *                 Lanczos drivers always fp64; the nb x nb eig-sqrt stays
 *                 on host LAPACK zheev, bit-comparable with crecal_b.
 *   L2         -> per-type bulk blocks pinned in the L2 persistence window
 *                 (CUDA >= 11) so psi streaming cannot evict them.
 *   structured -> cuFFT stencil + correction backend retained (fp64,
 *                 nrhs == nb), adapted to the site-major layout.
 *
 * Build:  nvcc -O3 -arch=native -Xcompiler -fPIC -shared rsrec_gpu.cu \
 *              -o librsrec_gpu.so -lcublas -lcufft -llapack
 * =========================================================================== */
#include "rsrec.h"

#include <cublas_v2.h>
#include <cufft.h>
#include <cuComplex.h>
#include <cuda_runtime.h>

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>

typedef cuDoubleComplex zc;
typedef cuComplex fc;
typedef std::complex<double> cplx;

extern "C" void zheev_(const char *, const char *, const int *, cplx *,
                       const int *, double *, cplx *, const int *, double *,
                       int *);

static std::string g_err;
extern "C" const char *rsrec_last_error(void) { return g_err.c_str(); }
#define FAIL(msg) do { g_err = (msg); return 1; } while (0)
#define CUCHK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    g_err = std::string("CUDA: ") + cudaGetErrorString(e); return 1; } } while (0)
#define CBCHK(x) do { if ((x) != CUBLAS_STATUS_SUCCESS) \
    { g_err = "cuBLAS error"; return 1; } } while (0)

/* --------------------------------------------------------------------------
 * Context
 * ------------------------------------------------------------------------ */
struct rsrec_ctx {
    int kk, nb, nnmax, ntype, nmax, device;

    zc *d_ee = nullptr, *d_hall = nullptr, *d_va = nullptr, *d_vb = nullptr;
    int *d_nn = nullptr, *d_iz = nullptr;
    bool have_h = false, have_v = false;
    int ham_ver = 0;

    /* Orthogonalisation (hoh) operands for the two-sweep apply. Uploaded only
     * when set_hamiltonian receives eeo/hallo/enim. d_ee_bare/d_hall_bare are
     * the BARE h (no lsham fold), d_eeo/d_hallo the eeo/hallo factor, and
     * d_hons the on-site (enim + lsham) per type (block-diagonal). The combine
     *   y = ( h*x - eeo*(h*x) + hons*x - b*x ) / a
     * mirrors ham_hoh_vec_matmul. Stored UNSCALED (raw); the 1/a and -b*I
     * scaling is applied in the combine kernel. f_* are fp32 mirrors used by
     * the default fp32 Chebyshev path (raw down-conversion, no scaling).     */
    zc *d_ee_bare = nullptr, *d_hall_bare = nullptr;
    zc *d_eeo = nullptr, *d_hallo = nullptr, *d_hons = nullptr;
    fc *f_ee_bare = nullptr, *f_hall_bare = nullptr;
    fc *f_eeo = nullptr, *f_hallo = nullptr, *f_hons = nullptr;
    int f_hoh_ver = -1;
    bool have_hoh = false;

    int cheb_prec = 0;                       /* 0 = fp32 (default), 1 = fp64 */
    fc *f_ee = nullptr, *f_hall = nullptr;   /* (H - b)/a fp32 copies        */
    double f_a = 0.0, f_b = 0.0;
    int f_ver = -1;

    /* site-major scratch fields, each ld*nb cplx (single-state size)        */
    zc *d_s0 = nullptr, *d_s1 = nullptr, *d_s2 = nullptr, *d_s3 = nullptr;
    zc *d_blk = nullptr;
    double *d_col = nullptr;
    cublasHandle_t blas = nullptr;

    /* structured (cuFFT stencil + correction) path                          */
    bool use_struct = false;
    int t_ref = 0, nshell_ref = 0;
    int N[3] = {0, 0, 0}, Np[3] = {0, 0, 0};
    size_t Npc = 0;
    std::vector<int> roff;
    int *d_gcell = nullptr;
    zc *d_gin = nullptr, *d_gw = nullptr;
    zc *d_Hq = nullptr, *d_VAq = nullptr, *d_VBq = nullptr;
    zc *d_negW = nullptr;
    cufftHandle fft_plan = 0;
    bool have_plan = false;
    std::vector<int> cg_off;
    int *d_c_atom = nullptr, *d_c_nbr = nullptr;
    const zc **d_cA_H = nullptr, **d_cA_VA = nullptr, **d_cA_VB = nullptr;
    const zc **d_cB = nullptr;
    zc **d_cC = nullptr;
    size_t c_nnz = 0;
};

static inline size_t fieldsz(const rsrec_ctx *c) {
    return (size_t)c->nb * c->nb * c->kk;
}

extern "C" void rsrec_destroy(rsrec_ctx *c) {
    if (!c) return;
    for (void *p : {(void *)c->d_ee, (void *)c->d_hall, (void *)c->d_va,
                    (void *)c->d_vb, (void *)c->d_s0, (void *)c->d_s1,
                    (void *)c->d_s2, (void *)c->d_s3, (void *)c->d_blk,
                    (void *)c->d_nn, (void *)c->d_iz, (void *)c->f_ee,
                    (void *)c->f_hall, (void *)c->d_gcell, (void *)c->d_gin,
                    (void *)c->d_gw, (void *)c->d_Hq, (void *)c->d_VAq,
                    (void *)c->d_VBq, (void *)c->d_negW, (void *)c->d_c_atom,
                    (void *)c->d_c_nbr, (void *)c->d_cA_H,
                    (void *)c->d_cA_VA, (void *)c->d_cA_VB, (void *)c->d_cB,
                    (void *)c->d_cC, (void *)c->d_col,
                    (void *)c->d_ee_bare, (void *)c->d_hall_bare,
                    (void *)c->d_eeo, (void *)c->d_hallo, (void *)c->d_hons,
                    (void *)c->f_ee_bare, (void *)c->f_hall_bare,
                    (void *)c->f_eeo, (void *)c->f_hallo, (void *)c->f_hons})
        if (p) cudaFree(p);
    if (c->have_plan) cufftDestroy(c->fft_plan);
    if (c->blas) cublasDestroy(c->blas);
    delete c;
}

extern "C" rsrec_ctx *rsrec_create(int kk, int nb, int nnmax, int ntype,
                                   int nmax, int device) {
    if (kk <= 0 || nb <= 0 || nb > 32 || nnmax <= 0 || ntype <= 0 ||
        nmax < 0 || nmax > kk) {
        g_err = "rsrec_create: bad dimensions (need 0 < nb <= 32)";
        return nullptr;
    }
    rsrec_ctx *c = new rsrec_ctx();
    c->kk = kk; c->nb = nb; c->nnmax = nnmax; c->ntype = ntype;
    c->nmax = nmax; c->device = device;
    if (cudaSetDevice(device) != cudaSuccess ||
        cublasCreate(&c->blas) != CUBLAS_STATUS_SUCCESS) {
        g_err = "rsrec_create: CUDA/cuBLAS init failed";
        delete c; return nullptr;
    }
    size_t nf = fieldsz(c) * sizeof(zc);
    if (cudaMalloc(&c->d_s0, nf) || cudaMalloc(&c->d_s1, nf) ||
        cudaMalloc(&c->d_s2, nf) || cudaMalloc(&c->d_s3, nf) ||
        cudaMalloc(&c->d_blk, 8 * (size_t)nb * nb * sizeof(zc)) ||
        cudaMalloc(&c->d_col, 2 * (size_t)nb * sizeof(double))) {
        g_err = "rsrec_create: out of device memory for scratch fields";
        rsrec_destroy(c); return nullptr;
    }
    return c;
}

/* L2 persistence: pin per-type operator blocks (small, reused by every atom
 * every step). No-op below CUDA 11 / unsupported hardware.                  */
static void l2_pin(const void *ptr, size_t bytes) {
#if defined(CUDART_VERSION) && CUDART_VERSION >= 11000
    int dev = 0; cudaGetDevice(&dev);
    cudaDeviceProp prop; cudaGetDeviceProperties(&prop, dev);
    if (prop.persistingL2CacheMaxSize <= 0) return;
    size_t win = bytes < (size_t)prop.persistingL2CacheMaxSize
                     ? bytes : (size_t)prop.persistingL2CacheMaxSize;
    cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize, win);
    cudaStreamAttrValue attr;
    attr.accessPolicyWindow.base_ptr = const_cast<void *>(ptr);
    attr.accessPolicyWindow.num_bytes = win;
    attr.accessPolicyWindow.hitRatio = 1.0f;
    attr.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;
    attr.accessPolicyWindow.missProp = cudaAccessPropertyStreaming;
    cudaStreamSetAttribute(cudaStreamDefault,
                           cudaStreamAttributeAccessPolicyWindow, &attr);
#else
    (void)ptr; (void)bytes;
#endif
}

/* --------------------------------------------------------------------------
 * Upload (lsham folded into shell-1 once, like the CPU reference)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_set_hamiltonian(rsrec_ctx *c, const void *ee_,
                                     const void *hall_, const void *lsham_,
                                     const int *nn, const int *iz,
                                     const void *eeo_, const void *hallo_,
                                     const void *enim_) {
    const cplx *ee = (const cplx *)ee_;
    const cplx *hall = (const cplx *)hall_;
    const cplx *ls = (const cplx *)lsham_;
    const cplx *eeo = (const cplx *)eeo_;
    const cplx *hallo = (const cplx *)hallo_;
    const cplx *enim = (const cplx *)enim_;
    if (!ee || !nn || !iz) FAIL("set_hamiltonian: null input");
    if (c->nmax > 0 && !hall) FAIL("set_hamiltonian: nmax>0 but hall is null");
    const int nb = c->nb, nnmax = c->nnmax;
    const size_t bb = (size_t)nb * nb;
    const bool hoh = (eeo != nullptr);
    if (hoh && c->nmax > 0 && !hallo)
        FAIL("set_hamiltonian: hoh with nmax>0 but hallo is null");

    std::vector<cplx> hee(ee, ee + bb * nnmax * c->ntype);
    std::vector<cplx> hha;
    if (c->nmax > 0) hha.assign(hall, hall + bb * nnmax * c->nmax);
    if (ls) {
        for (int t = 0; t < c->ntype; ++t)
            for (size_t i = 0; i < bb; ++i)
                hee[i + bb * (size_t)nnmax * t] += ls[i + bb * t];
        for (int k = 0; k < c->nmax; ++k) {
            int t = iz[k] - 1;
            for (size_t i = 0; i < bb; ++i)
                hha[i + bb * (size_t)nnmax * k] += ls[i + bb * t];
        }
    }
    if (!c->d_ee) CUCHK(cudaMalloc(&c->d_ee, hee.size() * sizeof(zc)));
    CUCHK(cudaMemcpy(c->d_ee, hee.data(), hee.size() * sizeof(zc),
                     cudaMemcpyHostToDevice));
    if (c->nmax > 0) {
        if (!c->d_hall) CUCHK(cudaMalloc(&c->d_hall, hha.size() * sizeof(zc)));
        CUCHK(cudaMemcpy(c->d_hall, hha.data(), hha.size() * sizeof(zc),
                         cudaMemcpyHostToDevice));
    }
    if (!c->d_nn)
        CUCHK(cudaMalloc(&c->d_nn, (size_t)c->kk * nnmax * sizeof(int)));
    CUCHK(cudaMemcpy(c->d_nn, nn, (size_t)c->kk * nnmax * sizeof(int),
                     cudaMemcpyHostToDevice));
    if (!c->d_iz) CUCHK(cudaMalloc(&c->d_iz, c->kk * sizeof(int)));
    CUCHK(cudaMemcpy(c->d_iz, iz, c->kk * sizeof(int),
                     cudaMemcpyHostToDevice));
    l2_pin(c->d_ee, bb * nnmax * c->ntype * sizeof(zc));

    /* --- two-sweep hoh operands (uploaded raw; combine applies 1/a) ------ */
    c->have_hoh = false;
    if (hoh) {
        /* bare h: ee/hall WITHOUT the lsham fold (the eeo sweep must see the
         * bare h*x; lsham + enim are carried on-site by d_hons).            */
        if (!c->d_ee_bare) CUCHK(cudaMalloc(&c->d_ee_bare, bb * nnmax * c->ntype * sizeof(zc)));
        CUCHK(cudaMemcpy(c->d_ee_bare, ee, bb * nnmax * c->ntype * sizeof(zc),
                         cudaMemcpyHostToDevice));
        if (!c->d_eeo) CUCHK(cudaMalloc(&c->d_eeo, bb * nnmax * c->ntype * sizeof(zc)));
        CUCHK(cudaMemcpy(c->d_eeo, eeo, bb * nnmax * c->ntype * sizeof(zc),
                         cudaMemcpyHostToDevice));
        if (c->nmax > 0) {
            if (!c->d_hall_bare) CUCHK(cudaMalloc(&c->d_hall_bare, bb * nnmax * c->nmax * sizeof(zc)));
            CUCHK(cudaMemcpy(c->d_hall_bare, hall, bb * nnmax * c->nmax * sizeof(zc),
                             cudaMemcpyHostToDevice));
            if (!c->d_hallo) CUCHK(cudaMalloc(&c->d_hallo, bb * nnmax * c->nmax * sizeof(zc)));
            CUCHK(cudaMemcpy(c->d_hallo, hallo, bb * nnmax * c->nmax * sizeof(zc),
                             cudaMemcpyHostToDevice));
        }
        /* on-site (enim + lsham) per type, block-diagonal (bb per type). */
        std::vector<cplx> hons(bb * c->ntype, cplx(0.0, 0.0));
        for (int t = 0; t < c->ntype; ++t)
            for (size_t i = 0; i < bb; ++i) {
                cplx v(0.0, 0.0);
                if (enim) v += enim[i + bb * t];
                if (ls) v += ls[i + bb * t];
                hons[i + bb * t] = v;
            }
        if (!c->d_hons) CUCHK(cudaMalloc(&c->d_hons, bb * c->ntype * sizeof(zc)));
        CUCHK(cudaMemcpy(c->d_hons, hons.data(), bb * c->ntype * sizeof(zc),
                         cudaMemcpyHostToDevice));
        c->have_hoh = true;
    }
    c->have_h = true;
    c->ham_ver++;
    return 0;
}

extern "C" int rsrec_set_velocity(rsrec_ctx *c, const void *va,
                                  const void *vb) {
    if (!va || !vb) FAIL("set_velocity: null input");
    size_t n = (size_t)c->nb * c->nb * c->nnmax * c->ntype * sizeof(zc);
    if (!c->d_va) CUCHK(cudaMalloc(&c->d_va, n));
    if (!c->d_vb) CUCHK(cudaMalloc(&c->d_vb, n));
    CUCHK(cudaMemcpy(c->d_va, va, n, cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(c->d_vb, vb, n, cudaMemcpyHostToDevice));
    c->have_v = true;
    return 0;
}

extern "C" int rsrec_set_precision(rsrec_ctx *c, int prec) {
    if (prec < 0 || prec > 1) FAIL("set_precision: 0 = fp32, 1 = fp64");
    c->cheb_prec = prec;
    return 0;
}

/* --------------------------------------------------------------------------
 * Device complex helpers (overloaded for zc / fc)
 * ------------------------------------------------------------------------ */
__device__ __forceinline__ zc cfma(zc a, zc b, zc s) {
    s.x = fma(a.x, b.x, fma(-a.y, b.y, s.x));
    s.y = fma(a.x, b.y, fma(a.y, b.x, s.y));
    return s;
}
__device__ __forceinline__ fc cfma(fc a, fc b, fc s) {
    s.x = fmaf(a.x, b.x, fmaf(-a.y, b.y, s.x));
    s.y = fmaf(a.x, b.y, fmaf(a.y, b.x, s.y));
    return s;
}
template <class CT> __device__ __forceinline__ CT czero_v();
template <> __device__ __forceinline__ zc czero_v<zc>() {
    return make_cuDoubleComplex(0.0, 0.0);
}
template <> __device__ __forceinline__ fc czero_v<fc>() {
    return make_cuFloatComplex(0.0f, 0.0f);
}

/* --------------------------------------------------------------------------
 * Fused block-ELL step kernel, templated:
 *   y(:,c) = alpha * ((Op x1)(:,c) - bsc*x1(:,c)) * inva + beta * x0(:,c)
 * over columns [col0, col0+ncols). Site-major: element (l, col, atom k)
 * lives at l + NB*k + ld*col. One thread block per atom, threads
 * (NB, ncols/TC), each thread owns TC adjacent columns; shared staging
 * padded with LDS = NB+1 (even NB) for conflict-free column reads; inner
 * product fully unrolled into FMA registers.
 * ------------------------------------------------------------------------ */
template <class CT, int NB, int TC, class RT>
__global__ void k_step_t(int kk, int nnmax, int nmax,
                         const int *__restrict__ nn,
                         const int *__restrict__ iz,
                         const CT *__restrict__ btype,
                         const CT *__restrict__ bimp, int use_imp,
                         const CT *__restrict__ x1,
                         const CT *__restrict__ x0, CT *__restrict__ y,
                         RT alpha, RT beta, RT inva, RT bsc,
                         size_t ld, int col0, int ncols) {
    constexpr int LDS = (NB % 2 == 0) ? NB + 1 : NB;
    extern __shared__ unsigned char smem_raw[];
    CT *sH = (CT *)smem_raw;
    CT *sX = sH + (size_t)LDS * NB;
    const int k = blockIdx.x;
    if (k >= kk) return;
    const int l = threadIdx.x, ct = threadIdx.y;
    const int tid = l + NB * ct, nthr = NB * blockDim.y;
    const size_t bb = (size_t)NB * NB;
    const bool imp = use_imp && (k < nmax);
    const CT *base = imp ? bimp + bb * (size_t)nnmax * k
                         : btype + bb * (size_t)nnmax * (iz[k] - 1);
    CT acc[TC];
#pragma unroll
    for (int j = 0; j < TC; ++j) acc[j] = czero_v<CT>();

    const int nr = nn[k];
    for (int s = 0; s < nr; ++s) {
        const int nbr = (s == 0) ? k : nn[k + (size_t)kk * s] - 1;
        if (nbr < 0) continue;
        for (int i = tid; i < NB * NB; i += nthr)
            sH[(i % NB) + LDS * (i / NB)] = base[bb * s + i];
        for (int i = tid; i < NB * ncols; i += nthr) {
            int r = i % NB, cc = i / NB;
            sX[r + LDS * cc] =
                x1[(size_t)r + (size_t)NB * nbr + ld * (col0 + cc)];
        }
        __syncthreads();
#pragma unroll
        for (int m = 0; m < NB; ++m) {
            const CT h = sH[l + LDS * m];
#pragma unroll
            for (int j = 0; j < TC; ++j)
                acc[j] = cfma(h, sX[m + LDS * (ct * TC + j)], acc[j]);
        }
        __syncthreads();
    }
#pragma unroll
    for (int j = 0; j < TC; ++j) {
        const int cc = ct * TC + j;
        if (cc >= ncols) continue;
        const size_t idx = (size_t)l + (size_t)NB * k + ld * (col0 + cc);
        CT v = acc[j];
        if (bsc != (RT)0) {
            CT xs = x1[idx];
            v.x -= bsc * xs.x;
            v.y -= bsc * xs.y;
        }
        v.x *= inva; v.y *= inva;
        if (beta != (RT)0) {
            CT p = x0[idx];
            v.x = alpha * v.x + beta * p.x;
            v.y = alpha * v.y + beta * p.y;
        } else {
            v.x *= alpha; v.y *= alpha;
        }
        y[idx] = v;
    }
}

/* Dynamic-NB fallback (TC = 1).                                            */
template <class CT, class RT>
__global__ void k_step_dyn(int kk, int nb, int nnmax, int nmax,
                           const int *__restrict__ nn,
                           const int *__restrict__ iz,
                           const CT *__restrict__ btype,
                           const CT *__restrict__ bimp, int use_imp,
                           const CT *__restrict__ x1,
                           const CT *__restrict__ x0, CT *__restrict__ y,
                           RT alpha, RT beta, RT inva, RT bsc,
                           size_t ld, int col0, int ncols) {
    extern __shared__ unsigned char smem_raw[];
    CT *sH = (CT *)smem_raw;
    const int LDS = nb + 1;
    CT *sX = sH + (size_t)LDS * nb;
    const int k = blockIdx.x;
    if (k >= kk) return;
    const int l = threadIdx.x, ct = threadIdx.y;
    const int tid = l + nb * ct, nthr = nb * blockDim.y;
    const size_t bb = (size_t)nb * nb;
    const bool imp = use_imp && (k < nmax);
    const CT *base = imp ? bimp + bb * (size_t)nnmax * k
                         : btype + bb * (size_t)nnmax * (iz[k] - 1);
    CT acc = czero_v<CT>();
    const int nr = nn[k];
    for (int s = 0; s < nr; ++s) {
        const int nbr = (s == 0) ? k : nn[k + (size_t)kk * s] - 1;
        if (nbr < 0) continue;
        for (int i = tid; i < nb * nb; i += nthr)
            sH[(i % nb) + LDS * (i / nb)] = base[bb * s + i];
        for (int i = tid; i < nb * ncols; i += nthr) {
            int r = i % nb, cc = i / nb;
            sX[r + LDS * cc] =
                x1[(size_t)r + (size_t)nb * nbr + ld * (col0 + cc)];
        }
        __syncthreads();
        for (int m = 0; m < nb; ++m)
            acc = cfma(sH[l + LDS * m], sX[m + LDS * ct], acc);
        __syncthreads();
    }
    if (ct < ncols) {
        const size_t idx = (size_t)l + (size_t)nb * k + ld * (col0 + ct);
        CT v = acc;
        if (bsc != (RT)0) {
            CT xs = x1[idx];
            v.x -= bsc * xs.x;
            v.y -= bsc * xs.y;
        }
        v.x *= inva; v.y *= inva;
        if (beta != (RT)0) {
            CT p = x0[idx];
            v.x = alpha * v.x + beta * p.x;
            v.y = alpha * v.y + beta * p.y;
        } else {
            v.x *= alpha; v.y *= alpha;
        }
        y[idx] = v;
    }
}

/* Host dispatcher: template instantiation by NB, column slabs sized so the
 * dynamic shared stays <= 40 KB and blockDim <= 1024.                       */
template <class CT, class RT>
static int step_apply(rsrec_ctx *c, const CT *btype, const CT *bimp,
                      int use_imp, const CT *x1, const CT *x0, CT *y,
                      int nrhs, RT alpha, RT beta, RT inva, RT bsc) {
    const int nb = c->nb;
    const int LDS = (nb % 2 == 0) ? nb + 1 : nb;
    const size_t ld = (size_t)nb * c->kk;
    int cs_max = (40 * 1024 / (int)sizeof(CT) - LDS * nb) / LDS;
    cs_max = std::min(cs_max, 1024 / nb);
    cs_max = std::max(1, cs_max - (cs_max % 2));     /* keep even for TC=2  */
    for (int c0 = 0; c0 < nrhs; c0 += cs_max) {
        int nc = std::min(cs_max, nrhs - c0);
        size_t shmem = (size_t)LDS * (nb + nc) * sizeof(CT);
        bool done = false;
#define RSREC_LAUNCH(NBv, TCv)                                               \
        if (!done && nb == NBv && nc % TCv == 0) {                           \
            dim3 thr(NBv, nc / TCv);                                         \
            k_step_t<CT, NBv, TCv, RT><<<c->kk, thr, shmem>>>(               \
                c->kk, c->nnmax, c->nmax, c->d_nn, c->d_iz, btype, bimp,     \
                use_imp, x1, x0, y, alpha, beta, inva, bsc, ld, c0, nc);     \
            done = true;                                                     \
        }
        RSREC_LAUNCH(32, 2)
        RSREC_LAUNCH(18, 2)
        RSREC_LAUNCH(8, 1)
        RSREC_LAUNCH(2, 1)
#undef RSREC_LAUNCH
        if (!done) {
            dim3 thr(nb, nc);
            k_step_dyn<CT, RT><<<c->kk, thr, shmem>>>(
                c->kk, nb, c->nnmax, c->nmax, c->d_nn, c->d_iz, btype, bimp,
                use_imp, x1, x0, y, alpha, beta, inva, bsc, ld, c0, nc);
        }
        CUCHK(cudaGetLastError());
    }
    return 0;
}

/* --------------------------------------------------------------------------
 * Small generic kernels (site-major)
 * ------------------------------------------------------------------------ */
template <class CT, class RT>
__global__ void k_axpby(size_t n, RT alpha, const CT *__restrict__ x,
                        RT beta, CT *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        y[i].x = alpha * x[i].x + beta * y[i].x;
        y[i].y = alpha * x[i].y + beta * y[i].y;
    }
}

/* pack: atom-major fp64 (nb, ncols, kk) -> site-major CT (ld x ncols)      */
template <class CT>
/* pack: atom-major fp64 input -> site-major CT (ld x nb*nstates).
 * Input is the Fortran array (nb, nb_cols, kk, nstates), i.e. per state a
 * (nb, nb_cols, kk) block with states OUTERMOST (slowest). Flat input index
 *   i = l + nb*(col + nb_cols*(k + kk*s)).
 * Output column for (state s, block-col col) is s*nb_cols + col, so the
 * site-major destination is l + nb*k + ld*(s*nb_cols + col).               */
__global__ void k_pack(size_t n, int nb, int nb_cols, int kk, size_t ld,
                       const zc *__restrict__ in, CT *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int l = (int)(i % nb);
    int col = (int)((i / nb) % nb_cols);
    size_t k = (i / ((size_t)nb * nb_cols)) % (size_t)kk;
    size_t s = i / ((size_t)nb * nb_cols * kk);
    CT v; v.x = (decltype(v.x))in[i].x; v.y = (decltype(v.y))in[i].y;
    out[(size_t)l + (size_t)nb * k + ld * ((long)s * nb_cols + col)] = v;
}

/* unpack: site-major zc -> atom-major zc (for rsrec_op_apply only)         */
__global__ void k_unpack64(size_t n, int nb, int ncols, size_t ld,
                           const zc *__restrict__ in, zc *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int l = (int)(i % nb);
    int col = (int)((i / nb) % ncols);
    size_t k = i / ((size_t)nb * ncols);
    out[i] = in[(size_t)l + (size_t)nb * k + ld * col];
}

/* Scalar-Lanczos per-column helpers (site-major: column stride ld).
 * The per-column dot products diag(A^H B) come from a single gram gemm
 * (gram64) followed by this diagonal pick -- no atomics in the hot path.
 * G is the nb x nb gram matrix (column-major, leading dim nb); its (j,j)
 * real part is the column-j dot.                                            */
__global__ void k_diag_real(int nb, const zc *__restrict__ G,
                            double *__restrict__ out) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < nb) out[j] = G[(size_t)j + (size_t)nb * j].x;
}

__global__ void k_col_axpy(size_t ld, int ncols,
                           const double *__restrict__ alpha,
                           const zc *__restrict__ x, zc *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < ld * (size_t)ncols) {
        int col = (int)(i / ld);
        y[i].x += alpha[col] * x[i].x;
        y[i].y += alpha[col] * x[i].y;
    }
}

__global__ void k_col_shift(size_t ld, int ncols,
                            const double *__restrict__ s2,
                            zc *__restrict__ psi, zc *__restrict__ pmn) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < ld * (size_t)ncols) {
        int col = (int)(i / ld);
        double sf = sqrt(s2[col]), si = (sf > 0.0) ? 1.0 / sf : 0.0;
        zc newpsi = make_cuDoubleComplex(pmn[i].x * si, pmn[i].y * si);
        pmn[i] = make_cuDoubleComplex(-psi[i].x * sf, -psi[i].y * sf);
        psi[i] = newpsi;
    }
}

/* cuBLAS scalar constants and thin precision-overloaded wrappers           */
static const zc Z_ONE = {1.0, 0.0};
static const zc Z_ZERO = {0.0, 0.0};
static const zc Z_MONE = {-1.0, 0.0};
static const fc C_ONE = {1.0f, 0.0f};
static const fc C_ZERO = {0.0f, 0.0f};

static cublasStatus_t blas_gemm(cublasHandle_t h, cublasOperation_t ta,
                                cublasOperation_t tb, int m, int n, int k,
                                const zc *al, const zc *A, int lda,
                                const zc *B, int ldb, const zc *be, zc *C,
                                int ldc) {
    return cublasZgemm(h, ta, tb, m, n, k, al, A, lda, B, ldb, be, C, ldc);
}
static cublasStatus_t blas_gemm(cublasHandle_t h, cublasOperation_t ta,
                                cublasOperation_t tb, int m, int n, int k,
                                const fc *al, const fc *A, int lda,
                                const fc *B, int ldb, const fc *be, fc *C,
                                int ldc) {
    return cublasCgemm3m(h, ta, tb, m, n, k, al, A, lda, B, ldb, be, C, ldc);
}
static cublasStatus_t blas_gemm_sb(cublasHandle_t h, cublasOperation_t ta,
                                   cublasOperation_t tb, int m, int n, int k,
                                   const zc *al, const zc *A, int lda,
                                   long long sA, const zc *B, int ldb,
                                   long long sB, const zc *be, zc *C, int ldc,
                                   long long sC, int batch) {
    return cublasZgemmStridedBatched(h, ta, tb, m, n, k, al, A, lda, sA, B,
                                     ldb, sB, be, C, ldc, sC, batch);
}
static cublasStatus_t blas_gemm_sb(cublasHandle_t h, cublasOperation_t ta,
                                   cublasOperation_t tb, int m, int n, int k,
                                   const fc *al, const fc *A, int lda,
                                   long long sA, const fc *B, int ldb,
                                   long long sB, const fc *be, fc *C, int ldc,
                                   long long sC, int batch) {
    return cublasCgemm3mStridedBatched(h, ta, tb, m, n, k, al, A, lda, sA, B,
                                       ldb, sB, be, C, ldc, sC, batch);
}

/* B = U sqrt(ev) U^H, Bi = U ev^{-1/2} U^H on the host (zheev), like
 * crecal_b. dS/dB/dBi are device pointers.                                  */
static int host_eig_sqrt(rsrec_ctx *c, const zc *dS, zc *dB, zc *dBi) {
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb;
    std::vector<cplx> S(bb), Bm(bb), Bi(bb);
    CUCHK(cudaMemcpy(S.data(), dS, bb * sizeof(zc), cudaMemcpyDeviceToHost));
    std::vector<cplx> U = S;
    std::vector<double> ev(nb), rwork(3 * nb - 2);
    int lwork = nb * nb, info = 0;
    std::vector<cplx> work(lwork);
    zheev_("V", "U", &nb, U.data(), &nb, ev.data(), work.data(), &lwork,
           rwork.data(), &info);
    if (info != 0) FAIL("block_lanczos: zheev failed");
    for (int j = 0; j < nb; ++j)
        for (int i = 0; i < nb; ++i) {
            cplx sb(0, 0), sbi(0, 0);
            for (int l = 0; l < nb; ++l) {
                double lam = sqrt(ev[l] > 0.0 ? ev[l] : 0.0);
                cplx uu = U[i + (size_t)nb * l] *
                          std::conj(U[j + (size_t)nb * l]);
                sb += uu * lam;
                if (lam > 0.0) sbi += uu / lam;
            }
            Bm[i + (size_t)nb * j] = sb;
            Bi[i + (size_t)nb * j] = sbi;
        }
    CUCHK(cudaMemcpy(dB, Bm.data(), bb * sizeof(zc), cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(dBi, Bi.data(), bb * sizeof(zc),
                     cudaMemcpyHostToDevice));
    return 0;
}

/* ==========================================================================
 * Structured (cuFFT stencil + correction) path -- fp64, nrhs == nb,
 * site-major layout (atom block = nb rows at offset nb*k, column stride ld).
 * Channel grids stay channel-major: channel ch occupies a contiguous padded
 * volume of Npc cells.
 * ======================================================================== */
__global__ void k_scatter(size_t nfield, int nb, int kk, size_t ld,
                          size_t Npc, const int *__restrict__ gcell,
                          const zc *__restrict__ x, zc *__restrict__ g) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nfield) return;
    int l = (int)(i % nb);
    int k = (int)((i / nb) % kk);
    int col = (int)(i / ((size_t)nb * kk));
    g[(size_t)(l + nb * col) * Npc + gcell[k]] =
        x[(size_t)l + (size_t)nb * k + ld * col];
}

__global__ void k_gather(size_t nfield, int nb, int kk, size_t ld,
                         size_t Npc, const int *__restrict__ gcell,
                         const zc *__restrict__ g, zc *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nfield) return;
    int l = (int)(i % nb);
    int k = (int)((i / nb) % kk);
    int col = (int)(i / ((size_t)nb * kk));
    y[(size_t)l + (size_t)nb * k + ld * col] =
        g[(size_t)(l + nb * col) * Npc + gcell[k]];
}

__global__ void k_stencil_scatter(int nshell, int nbsq, size_t Npc,
                                  int Npx, int Npy, int Npz,
                                  const int *__restrict__ roff,
                                  const zc *__restrict__ W, double scale,
                                  zc *__restrict__ g) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nshell * nbsq) return;
    int s = idx / nbsq, ch = idx % nbsq;
    int rx = ((-roff[3 * s]) % Npx + Npx) % Npx;
    int ry = ((-roff[3 * s + 1]) % Npy + Npy) % Npy;
    int rz = ((-roff[3 * s + 2]) % Npz + Npz) % Npz;
    size_t cell = ((size_t)rx * Npy + ry) * Npz + rz;
    zc w = W[(size_t)s * nbsq + ch];
    g[(size_t)ch * Npc + cell] = make_cuDoubleComplex(w.x * scale,
                                                      w.y * scale);
}

__global__ void k_pointwise_blockmul(size_t Npc, int nb,
                                     const zc *__restrict__ Hq,
                                     zc *__restrict__ g) {
    extern __shared__ unsigned char smem_raw[];
    zc *sH = (zc *)smem_raw, *sx = sH + nb * nb;
    for (size_t cell = blockIdx.x; cell < Npc; cell += gridDim.x) {
        const int l = threadIdx.x, c2 = threadIdx.y;
        const int ch = l + nb * c2;
        sH[ch] = Hq[(size_t)ch * Npc + cell];
        sx[ch] = g[(size_t)ch * Npc + cell];
        __syncthreads();
        zc acc = make_cuDoubleComplex(0.0, 0.0);
        for (int m = 0; m < nb; ++m)
            acc = cfma(sH[l + nb * m], sx[m + nb * c2], acc);
        __syncthreads();
        g[(size_t)ch * Npc + cell] = acc;
    }
}

/* Correction x/y block pointers over the site-major layout (block stride
 * nb columns; ldb = ldc = ld).                                             */
__global__ void k_make_ptrs(size_t nnz, int nb,
                            const int *__restrict__ atom,
                            const int *__restrict__ nbr,
                            const zc *__restrict__ x, zc *__restrict__ y,
                            const zc **__restrict__ Bp,
                            zc **__restrict__ Cp) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nnz) {
        Bp[i] = x + (size_t)nb * nbr[i];
        Cp[i] = y + (size_t)nb * atom[i];
    }
}

extern "C" int rsrec_set_grid(rsrec_ctx *c, const int *coords,
                              int use_structured) {
    if (!c->have_h) FAIL("set_grid: call set_hamiltonian first");
    if (!use_structured) { c->use_struct = false; return 0; }
    const int kk = c->kk, nb = c->nb;
    const size_t bb = (size_t)nb * nb;

    /* ---- host copies of nn/iz (kept on device only so far) -------------- */
    std::vector<int> nn((size_t)kk * c->nnmax), iz(kk);
    CUCHK(cudaMemcpy(nn.data(), c->d_nn, nn.size() * sizeof(int),
                     cudaMemcpyDeviceToHost));
    CUCHK(cudaMemcpy(iz.data(), c->d_iz, kk * sizeof(int),
                     cudaMemcpyDeviceToHost));
    auto NNF = [&](int k, int s) { return nn[(size_t)k + (size_t)kk * s]; };

    /* ---- bounding box, cell map ----------------------------------------- */
    int lo[3];
    for (int d = 0; d < 3; ++d) {
        int mn = coords[d], mx = coords[d];
        for (int k = 1; k < kk; ++k) {
            mn = std::min(mn, coords[3 * k + d]);
            mx = std::max(mx, coords[3 * k + d]);
        }
        lo[d] = mn; c->N[d] = mx - mn + 1;
    }
    const int Nx = c->N[0], Ny = c->N[1], Nz = c->N[2];
    std::vector<int> cellmap((size_t)Nx * Ny * Nz, -1);
    std::vector<int> gc(kk);
    for (int k = 0; k < kk; ++k) {
        int cx = coords[3 * k] - lo[0], cy = coords[3 * k + 1] - lo[1];
        int cz = coords[3 * k + 2] - lo[2];
        size_t cell = (cx * (size_t)Ny + cy) * Nz + cz;
        if (cellmap[cell] >= 0)
            FAIL("set_grid: two atoms share one lattice cell (multi-atom "
                 "bases need a basis-resolved grid; not supported yet)");
        cellmap[cell] = k;
        gc[k] = (int)cell;
    }

    /* ---- reference type / template / stencil offsets -------------------- */
    if (kk == c->nmax)
        FAIL("set_grid: hall covers all atoms; no bulk stencil to exploit "
             "-- stay on the ELL backend");
    std::vector<int> cnt(c->ntype, 0);
    for (int k = c->nmax; k < kk; ++k) cnt[iz[k] - 1]++;
    c->t_ref = (int)(std::max_element(cnt.begin(), cnt.end()) - cnt.begin());
    int k0 = -1, best = -1;
    for (int k = c->nmax; k < kk; ++k)
        if (iz[k] - 1 == c->t_ref && NNF(k, 0) > best) {
            best = NNF(k, 0); k0 = k;
        }
    if (k0 < 0) FAIL("set_grid: no bulk atom of the reference type");
    c->nshell_ref = best;
    c->roff.assign(3 * c->nshell_ref, 0);
    for (int s = 1; s < c->nshell_ref; ++s) {
        int nbr = NNF(k0, s) - 1;
        if (nbr < 0) FAIL("set_grid: template atom has a hole in its "
                          "neighbour list; pick a cluster interior template");
        for (int d = 0; d < 3; ++d)
            c->roff[3 * s + d] = coords[3 * nbr + d] - coords[3 * k0 + d];
    }

    /* ---- flag rows the stencil cannot represent -------------------------- */
    auto shells_match = [&](int k) -> bool {
        const int nr = NNF(k, 0);
        if (nr > c->nshell_ref) return false;
        int cell = gc[k];
        int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
        for (int s = 1; s < c->nshell_ref; ++s) {
            int tx = cx + c->roff[3 * s], ty = cy + c->roff[3 * s + 1];
            int tz = cz + c->roff[3 * s + 2];
            bool inbox = (tx >= 0 && tx < Nx && ty >= 0 && ty < Ny &&
                          tz >= 0 && tz < Nz);
            int occ = inbox ? cellmap[(tx * (size_t)Ny + ty) * Nz + tz] : -1;
            int nbr = (s < nr) ? NNF(k, s) - 1 : -1;
            if (occ != nbr) return false;
        }
        return true;
    };
    std::vector<char> flagH(kk);
    size_t nflagH = 0;
    for (int k = 0; k < kk; ++k) {
        bool irr = !shells_match(k) || (iz[k] - 1 != c->t_ref);
        flagH[k] = irr || (k < c->nmax);
        nflagH += flagH[k];
    }

    /* ---- correction entries, grouped by per-atom slot (race-free) --------
     * One shared (atom, nbr) list built from the H flag set (a superset of
     * the velocity flag set: flagH = flagV or k < nmax). For the velocity
     * operators, rows flagged only because of hall are stencil-regular, so
     * their (+v_type, -v_ref) entries cancel exactly -- a few wasted GEMMs,
     * no error. Entries store (shell s | negW slot) + atom, from which the
     * per-operator A-pointers are derived.                                  */
    std::vector<int> ca, cn, cs;          /* atom, nbr, shell (>=0) or
                                             -(negW slot)-1                  */
    {
        std::vector<std::vector<int>> pa, pn, ps;
        for (int k = 0; k < kk; ++k) {
            if (!flagH[k]) continue;
            std::vector<int> ea, en, es;
            const int nr = NNF(k, 0);
            for (int s = 0; s < nr; ++s) {
                int nb_at = (s == 0) ? k : NNF(k, s) - 1;
                if (nb_at < 0) continue;
                ea.push_back(k); en.push_back(nb_at); es.push_back(s);
            }
            int cell = gc[k];
            int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
            for (int s = 0; s < c->nshell_ref; ++s) {
                int tx = cx + c->roff[3 * s], ty = cy + c->roff[3 * s + 1];
                int tz = cz + c->roff[3 * s + 2];
                if (tx < 0 || tx >= Nx || ty < 0 || ty >= Ny ||
                    tz < 0 || tz >= Nz) continue;
                int occ = cellmap[(tx * (size_t)Ny + ty) * Nz + tz];
                if (occ < 0) continue;
                ea.push_back(k); en.push_back(occ); es.push_back(-s - 1);
            }
            pa.push_back(std::move(ea)); pn.push_back(std::move(en));
            ps.push_back(std::move(es));
        }
        c->cg_off.assign(1, 0);
        size_t maxlen = 0;
        for (auto &v : pa) maxlen = std::max(maxlen, v.size());
        for (size_t g = 0; g < maxlen; ++g) {
            for (size_t i = 0; i < pa.size(); ++i)
                if (g < pa[i].size()) {
                    ca.push_back(pa[i][g]); cn.push_back(pn[i][g]);
                    cs.push_back(ps[i][g]);
                }
            c->cg_off.push_back((int)ca.size());
        }
    }
    c->c_nnz = ca.size();

    /* materialize -W_ref for all three ops in one buffer                    */
    if (c->d_negW) cudaFree(c->d_negW);
    CUCHK(cudaMalloc(&c->d_negW, 3 * bb * c->nshell_ref * sizeof(zc)));
    {
        std::vector<cplx> wref(bb * c->nshell_ref);
        std::vector<cplx> neg(3 * bb * c->nshell_ref, cplx(0, 0));
        const zc *bases[3] = {c->d_ee, c->d_va, c->d_vb};
        for (int op = 0; op < 3; ++op) {
            if (!bases[op]) continue;
            CUCHK(cudaMemcpy(wref.data(),
                             bases[op] + bb * (size_t)c->nnmax * c->t_ref,
                             bb * c->nshell_ref * sizeof(zc),
                             cudaMemcpyDeviceToHost));
            for (size_t i = 0; i < wref.size(); ++i)
                neg[(size_t)op * bb * c->nshell_ref + i] = -wref[i];
        }
        CUCHK(cudaMemcpy(c->d_negW, neg.data(), neg.size() * sizeof(cplx),
                         cudaMemcpyHostToDevice));
    }

    /* per-operator A-pointer arrays over the shared entry list              */
    if (c->c_nnz > 0) {
        auto upload_A = [&](int op, const zc *type_base, const zc *imp_base,
                            const zc ***dst) -> int {
            std::vector<const zc *> A(c->c_nnz);
            for (size_t i = 0; i < c->c_nnz; ++i) {
                int k = ca[i], s = cs[i];
                if (s < 0)                                   /* -W_ref      */
                    A[i] = c->d_negW + (size_t)op * bb * c->nshell_ref +
                           bb * (size_t)(-s - 1);
                else if (op == 0 && k < c->nmax)             /* hall row    */
                    A[i] = imp_base +
                           bb * ((size_t)s + (size_t)c->nnmax * k);
                else                                         /* per type    */
                    A[i] = type_base +
                           bb * ((size_t)s + (size_t)c->nnmax * (iz[k] - 1));
            }
            if (*dst) cudaFree((void *)*dst);
            CUCHK(cudaMalloc(dst, c->c_nnz * sizeof(zc *)));
            CUCHK(cudaMemcpy((void *)*dst, A.data(),
                             c->c_nnz * sizeof(zc *),
                             cudaMemcpyHostToDevice));
            return 0;
        };
        if (upload_A(0, c->d_ee, c->d_hall, &c->d_cA_H)) return 1;
        if (c->have_v) {
            if (upload_A(1, c->d_va, nullptr, &c->d_cA_VA)) return 1;
            if (upload_A(2, c->d_vb, nullptr, &c->d_cA_VB)) return 1;
        }
        if (c->d_c_atom) cudaFree(c->d_c_atom);
        if (c->d_c_nbr) cudaFree(c->d_c_nbr);
        CUCHK(cudaMalloc(&c->d_c_atom, c->c_nnz * sizeof(int)));
        CUCHK(cudaMalloc(&c->d_c_nbr, c->c_nnz * sizeof(int)));
        CUCHK(cudaMemcpy(c->d_c_atom, ca.data(), c->c_nnz * sizeof(int),
                         cudaMemcpyHostToDevice));
        CUCHK(cudaMemcpy(c->d_c_nbr, cn.data(), c->c_nnz * sizeof(int),
                         cudaMemcpyHostToDevice));
        if (c->d_cB) cudaFree((void *)c->d_cB);
        if (c->d_cC) cudaFree((void *)c->d_cC);
        CUCHK(cudaMalloc(&c->d_cB, c->c_nnz * sizeof(zc *)));
        CUCHK(cudaMalloc(&c->d_cC, c->c_nnz * sizeof(zc *)));
    }

    /* ---- padded box, grids, H(q) tables, cuFFT plan ---------------------- */
    int rmax[3] = {0, 0, 0};
    for (int s = 0; s < c->nshell_ref; ++s)
        for (int d = 0; d < 3; ++d)
            rmax[d] = std::max(rmax[d], std::abs(c->roff[3 * s + d]));
    for (int d = 0; d < 3; ++d) c->Np[d] = c->N[d] + 2 * rmax[d];
    c->Npc = (size_t)c->Np[0] * c->Np[1] * c->Np[2];

    std::vector<int> gpad(kk);
    for (int k = 0; k < kk; ++k) {
        int cell = gc[k];
        int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
        gpad[k] = (int)((cx * (size_t)c->Np[1] + cy) * c->Np[2] + cz);
    }
    if (c->d_gcell) cudaFree(c->d_gcell);
    CUCHK(cudaMalloc(&c->d_gcell, kk * sizeof(int)));
    CUCHK(cudaMemcpy(c->d_gcell, gpad.data(), kk * sizeof(int),
                     cudaMemcpyHostToDevice));

    size_t gridbytes = c->Npc * bb * sizeof(zc);
    for (zc **p : {&c->d_gin, &c->d_gw, &c->d_Hq})
        { if (*p) cudaFree(*p); CUCHK(cudaMalloc(p, gridbytes)); }
    if (c->have_v)
        for (zc **p : {&c->d_VAq, &c->d_VBq})
            { if (*p) cudaFree(*p); CUCHK(cudaMalloc(p, gridbytes)); }

    if (c->have_plan) cufftDestroy(c->fft_plan);
    if (cufftPlanMany(&c->fft_plan, 3, c->Np, nullptr, 1, (int)c->Npc,
                      nullptr, 1, (int)c->Npc, CUFFT_Z2Z, nb * nb)
        != CUFFT_SUCCESS)
        FAIL("set_grid: cufftPlanMany failed");
    c->have_plan = true;

    /* build H(q) per operator: scatter W(-R)/Npc then forward FFT           */
    int *d_roff;
    CUCHK(cudaMalloc(&d_roff, c->roff.size() * sizeof(int)));
    CUCHK(cudaMemcpy(d_roff, c->roff.data(), c->roff.size() * sizeof(int),
                     cudaMemcpyHostToDevice));
    const zc *opbase[3] = {c->d_ee, c->d_va, c->d_vb};
    zc *optab[3] = {c->d_Hq, c->d_VAq, c->d_VBq};
    for (int op = 0; op < 3; ++op) {
        if (!opbase[op] || !optab[op]) continue;
        CUCHK(cudaMemset(c->d_gw, 0, gridbytes));
        int nthread = c->nshell_ref * (int)bb;
        k_stencil_scatter<<<(nthread + 255) / 256, 256>>>(
            c->nshell_ref, (int)bb, c->Npc, c->Np[0], c->Np[1], c->Np[2],
            d_roff, opbase[op] + bb * (size_t)c->nnmax * c->t_ref,
            1.0 / (double)c->Npc, c->d_gw);
        CUCHK(cudaGetLastError());
        if (cufftExecZ2Z(c->fft_plan, (cufftDoubleComplex *)c->d_gw,
                         (cufftDoubleComplex *)optab[op], CUFFT_FORWARD)
            != CUFFT_SUCCESS)
            FAIL("set_grid: cufftExecZ2Z failed building H(q)");
    }
    cudaFree(d_roff);
    CUCHK(cudaMemset(c->d_gin, 0, gridbytes));   /* pristine zero background */

    fprintf(stderr,
            "rsrec_set_grid: box %dx%dx%d (padded %dx%dx%d), ref type %d, "
            "%d shells, corrected rows: %zu/%d, device tables %.0f MB\n",
            Nx, Ny, Nz, c->Np[0], c->Np[1], c->Np[2], c->t_ref + 1,
            c->nshell_ref, nflagH, kk,
            (double)gridbytes * (2 + (c->have_v ? 3 : 1)) / 1e6);
    c->use_struct = true;
    return 0;
}

/* Structured y = Op x (UNSCALED, site-major, nrhs == nb). Callers apply the
 * shift/recurrence via k_combine.                                          */
static int struct_apply(rsrec_ctx *c, int which, const zc *x, zc *y) {
    const int nb = c->nb;
    const size_t nf = fieldsz(c), bb = (size_t)nb * nb;
    const size_t ld = (size_t)nb * c->kk;
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    const zc *Hq = (which == 1) ? c->d_VAq : (which == 2) ? c->d_VBq
                                                          : c->d_Hq;
    if (!Hq) FAIL("struct_apply: operator table missing (call set_velocity "
                  "before set_grid)");
    k_scatter<<<nbl, tpb>>>(nf, nb, c->kk, ld, c->Npc, c->d_gcell, x,
                            c->d_gin);
    CUCHK(cudaGetLastError());
    if (cufftExecZ2Z(c->fft_plan, (cufftDoubleComplex *)c->d_gin,
                     (cufftDoubleComplex *)c->d_gw, CUFFT_FORWARD)
        != CUFFT_SUCCESS) FAIL("struct_apply: forward FFT failed");
    {
        dim3 thr(nb, nb);
        int nblocks = (int)std::min<size_t>(c->Npc, 65535);
        k_pointwise_blockmul<<<nblocks, thr, 2 * bb * sizeof(zc)>>>(
            c->Npc, nb, Hq, c->d_gw);
        CUCHK(cudaGetLastError());
    }
    if (cufftExecZ2Z(c->fft_plan, (cufftDoubleComplex *)c->d_gw,
                     (cufftDoubleComplex *)c->d_gw, CUFFT_INVERSE)
        != CUFFT_SUCCESS) FAIL("struct_apply: inverse FFT failed");
    k_gather<<<nbl, tpb>>>(nf, nb, c->kk, ld, c->Npc, c->d_gcell, c->d_gw, y);
    CUCHK(cudaGetLastError());

    if (c->c_nnz > 0) {
        const zc **Ap = (which == 1) ? c->d_cA_VA :
                        (which == 2) ? c->d_cA_VB : c->d_cA_H;
        k_make_ptrs<<<(int)((c->c_nnz + tpb - 1) / tpb), tpb>>>(
            c->c_nnz, nb, c->d_c_atom, c->d_c_nbr, x, y, c->d_cB, c->d_cC);
        CUCHK(cudaGetLastError());
        for (size_t g = 0; g + 1 < c->cg_off.size(); ++g) {
            int off = c->cg_off[g], cnt2 = c->cg_off[g + 1] - off;
            if (cnt2 <= 0) continue;
            CBCHK(cublasZgemmBatched(c->blas, CUBLAS_OP_N, CUBLAS_OP_N,
                                     nb, nb, nb, &Z_ONE,
                                     Ap + off, nb, c->d_cB + off, (int)ld,
                                     &Z_ONE, c->d_cC + off, (int)ld, cnt2));
        }
    }
    return 0;
}

/* y = alpha*((t - bsc*x1)*inva) + beta*x0, applied after struct_apply.      */
__global__ void k_combine(size_t n, const zc *__restrict__ t,
                          const zc *__restrict__ x1,
                          const zc *__restrict__ x0, zc *__restrict__ y,
                          double alpha, double beta, double inva,
                          double bsc) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double vx = (t[i].x - bsc * x1[i].x) * inva;
    double vy = (t[i].y - bsc * x1[i].y) * inva;
    if (beta != 0.0)
        y[i] = make_cuDoubleComplex(alpha * vx + beta * x0[i].x,
                                    alpha * vy + beta * x0[i].y);
    else
        y[i] = make_cuDoubleComplex(alpha * vx, alpha * vy);
}

/* On-site block-diagonal multiply o = hons[type(k)] * x1, site-major.
 * One block per atom; hons is (bb per type), x1/o are (ld * ncols).         */
template <class CT>
__global__ void k_onsite_blockmul(int kk, int nb, size_t ld, int ncols,
                                  const int *__restrict__ iz,
                                  const CT *__restrict__ hons,
                                  const CT *__restrict__ x1,
                                  CT *__restrict__ o) {
    const int k = blockIdx.x;
    if (k >= kk) return;
    const int l = threadIdx.x, cc = threadIdx.y;
    if (l >= nb || cc >= ncols) return;
    const size_t bb = (size_t)nb * nb;
    const CT *H = hons + bb * (size_t)(iz[k] - 1);
    CT acc = czero_v<CT>();
    for (int m = 0; m < nb; ++m) {
        CT h = H[l + nb * m];
        CT xv = x1[(size_t)m + (size_t)nb * k + ld * cc];
        acc.x += h.x * xv.x - h.y * xv.y;
        acc.y += h.x * xv.y + h.y * xv.x;
    }
    o[(size_t)l + (size_t)nb * k + ld * cc] = acc;
}

/* hoh combine: y = alpha*((t - e + o - bsc*x1)*inva) + beta*x0, elementwise.
 * t = h*x1, e = eeo*(h*x1), o = hons*x1 (all unscaled).                      */
template <class CT, class RT>
__global__ void k_combine_hoh(size_t n, const CT *__restrict__ t,
                              const CT *__restrict__ e,
                              const CT *__restrict__ o,
                              const CT *__restrict__ x1,
                              const CT *__restrict__ x0, CT *__restrict__ y,
                              RT alpha, RT beta, RT inva, RT bsc) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    RT vx = (t[i].x - e[i].x + o[i].x - bsc * x1[i].x) * inva;
    RT vy = (t[i].y - e[i].y + o[i].y - bsc * x1[i].y) * inva;
    CT r;
    if (beta != (RT)0) {
        r.x = alpha * vx + beta * x0[i].x;
        r.y = alpha * vy + beta * x0[i].y;
    } else {
        r.x = alpha * vx;
        r.y = alpha * vy;
    }
    y[i] = r;
}

/* Two-sweep hoh step (block-ELL), mirrors ham_hoh_vec_matmul:
 *   sweep A: t = h * x1            (bare ee/hall, unscaled)
 *   sweep B: e = eeo * t           (eeo/hallo, unscaled)
 *   on-site: o = (enim+lsham)*x1
 *   combine: y = alpha*((t - e + o - bsc*x1)*inva) + beta*x0
 * Operator pointers (bare/eeo/hons) and scratch (t/e/o) are caller-provided
 * in the engine's precision CT, so the whole apply runs at that precision.  */
template <class CT, class RT>
static int step_apply_hoh(rsrec_ctx *c, const CT *bare, const CT *bare_imp,
                          const CT *eeo, const CT *eeo_imp, const CT *hons,
                          CT *t, CT *e, CT *o, const CT *x1, const CT *x0,
                          CT *y, int nrhs, RT alpha, RT beta, RT inva,
                          RT bsc) {
    const size_t ld = (size_t)c->nb * c->kk;
    const size_t n = ld * (size_t)nrhs;
    /* sweep A: t = h*x1 (alpha=1,beta=0,inva=1,bsc=0 -> bare apply) */
    if (step_apply<CT, RT>(c, bare, bare_imp, c->nmax > 0, x1, nullptr, t, nrhs,
                           (RT)1, (RT)0, (RT)1, (RT)0)) return 1;
    /* sweep B: e = eeo*t */
    if (step_apply<CT, RT>(c, eeo, eeo_imp, c->nmax > 0, t, nullptr, e, nrhs,
                           (RT)1, (RT)0, (RT)1, (RT)0)) return 1;
    /* on-site: o = hons*x1 (loop column slabs so blockDim.y stays <= 1024) */
    {
        int cs = std::max(1, 1024 / c->nb);
        for (int c0 = 0; c0 < nrhs; c0 += cs) {
            int nc = std::min(cs, nrhs - c0);
            dim3 t2(c->nb, nc);
            k_onsite_blockmul<CT><<<c->kk, t2>>>(c->kk, c->nb, ld, nc, c->d_iz,
                                                 hons, x1 + ld * (size_t)c0,
                                                 o + ld * (size_t)c0);
            CUCHK(cudaGetLastError());
        }
    }
    /* combine */
    {
        const int tpb = 256;
        k_combine_hoh<CT, RT><<<(int)((n + tpb - 1) / tpb), tpb>>>(
            n, t, e, o, x1, x0, y, alpha, beta, inva, bsc);
        CUCHK(cudaGetLastError());
    }
    return 0;
}

/* fp64 generalized step (site-major): y = alpha*(Op x1 - bsc x1)*inva +
 * beta*x0. Routes through the structured backend when active (needs tmp).   */
static int step64(rsrec_ctx *c, int which, const zc *x1, const zc *x0,
                  zc *y, zc *tmp, int nrhs, double alpha, double beta,
                  double inva, double bsc) {
    if (c->use_struct) {
        if (nrhs != c->nb)
            FAIL("structured backend supports nrhs == nb only");
        if (struct_apply(c, which, x1, tmp)) return 1;
        size_t n = (size_t)c->nb * c->kk * nrhs;
        const int tpb = 256;
        k_combine<<<(int)((n + tpb - 1) / tpb), tpb>>>(n, tmp, x1, x0, y,
                                                       alpha, beta, inva,
                                                       bsc);
        CUCHK(cudaGetLastError());
        return 0;
    }
    const zc *bt = (which == 1) ? c->d_va : (which == 2) ? c->d_vb : c->d_ee;
    return step_apply<zc, double>(c, bt, c->d_hall, which == 0 ? 1 : 0,
                                  x1, x0, y, nrhs, alpha, beta, inva, bsc);
}

/* C = A^H B over site-major fields: ONE tall-skinny gemm (na x nbc x ld).   */
static int gram64(rsrec_ctx *c, const zc *A, int na, const zc *B, int nbc,
                  zc *C) {
    const int ld = c->nb * c->kk;
    CBCHK(blas_gemm(c->blas, CUBLAS_OP_C, CUBLAS_OP_N, na, nbc, ld, &Z_ONE,
                    A, ld, B, ld, &Z_ZERO, C, na));
    return 0;
}

/* --------------------------------------------------------------------------
 * API: matvec (host arrays atom-major in/out; pack/unpack at the edges)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_op_apply(rsrec_ctx *c, int which, const void *x_,
                              void *y_, int nrhs, double a, double b) {
    if (!c->have_h) FAIL("op_apply: Hamiltonian not set");
    if (which != 0 && !c->have_v) FAIL("op_apply: velocity operators not set");
    if (nrhs != c->nb) FAIL("op_apply: this build expects nrhs == nb");
    const size_t nf = fieldsz(c), ld = (size_t)c->nb * c->kk;
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    CUCHK(cudaMemcpy(c->d_s3, x_, nf * sizeof(zc), cudaMemcpyHostToDevice));
    k_pack<zc><<<nbl, tpb>>>(nf, c->nb, nrhs, c->kk, ld, c->d_s3, c->d_s0);
    CUCHK(cudaGetLastError());
    if (step64(c, which, c->d_s0, nullptr, c->d_s1, c->d_s2, nrhs,
               1.0, 0.0, 1.0 / a, b)) return 1;
    k_unpack64<<<nbl, tpb>>>(nf, c->nb, nrhs, ld, c->d_s1, c->d_s3);
    CUCHK(cudaGetLastError());
    CUCHK(cudaMemcpy(y_, c->d_s3, nf * sizeof(zc), cudaMemcpyDeviceToHost));
    return 0;
}

/* ==========================================================================
 * Chebyshev block moments: templated engine over CT in {fc, zc}, multi-state
 * RHS batching (nrhs = nb * nstates). States are column groups
 * [s*nb, (s+1)*nb); per-state moments via gemm strided-batched over states.
 * ======================================================================== */
__global__ void k_ham_to_f32(size_t n, int nb, int nnmax, double inva,
                             double b, const zc *__restrict__ in,
                             fc *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    size_t bb = (size_t)nb * nb;
    int l = (int)(i % nb), m = (int)((i / nb) % nb);
    int s = (int)((i / bb) % nnmax);
    double re = in[i].x, im = in[i].y;
    if (l == m && s == 0) re -= b;
    out[i] = make_cuFloatComplex((float)(re * inva), (float)(im * inva));
}

static int ensure_f32_ham(rsrec_ctx *c, double a, double b) {
    if (c->f_ee && c->f_ver == c->ham_ver && c->f_a == a && c->f_b == b)
        return 0;
    const size_t bb = (size_t)c->nb * c->nb;
    const size_t nee = bb * c->nnmax * c->ntype;
    const size_t nha = bb * c->nnmax * (size_t)c->nmax;
    if (!c->f_ee) CUCHK(cudaMalloc(&c->f_ee, nee * sizeof(fc)));
    if (c->nmax > 0 && !c->f_hall)
        CUCHK(cudaMalloc(&c->f_hall, nha * sizeof(fc)));
    const int tpb = 256;
    k_ham_to_f32<<<(int)((nee + tpb - 1) / tpb), tpb>>>(
        nee, c->nb, c->nnmax, 1.0 / a, b, c->d_ee, c->f_ee);
    CUCHK(cudaGetLastError());
    if (c->nmax > 0) {
        k_ham_to_f32<<<(int)((nha + tpb - 1) / tpb), tpb>>>(
            nha, c->nb, c->nnmax, 1.0 / a, b, c->d_hall, c->f_hall);
        CUCHK(cudaGetLastError());
    }
    l2_pin(c->f_ee, nee * sizeof(fc));
    c->f_a = a; c->f_b = b; c->f_ver = c->ham_ver;
    return 0;
}

/* Plain fp64 -> fp32 down-conversion (no scaling): out[i] = (fc)in[i]. */
__global__ void k_zc_to_fc(size_t n, const zc *__restrict__ in,
                           fc *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = make_cuFloatComplex((float)in[i].x, (float)in[i].y);
}

/* fp32 mirrors of the (unscaled) hoh operands. The 1/a, -b*I scaling is
 * applied in k_combine_hoh, so these are pure precision down-conversions and
 * depend only on ham_ver (not on a/b). */
static int ensure_f32_hoh(rsrec_ctx *c) {
    if (!c->have_hoh) return 0;
    if (c->f_ee_bare && c->f_hoh_ver == c->ham_ver) return 0;
    const size_t bb = (size_t)c->nb * c->nb;
    const size_t nee = bb * c->nnmax * c->ntype;
    const size_t nha = bb * c->nnmax * (size_t)c->nmax;
    const size_t nons = bb * c->ntype;
    if (!c->f_ee_bare) CUCHK(cudaMalloc(&c->f_ee_bare, nee * sizeof(fc)));
    if (!c->f_eeo) CUCHK(cudaMalloc(&c->f_eeo, nee * sizeof(fc)));
    if (!c->f_hons) CUCHK(cudaMalloc(&c->f_hons, nons * sizeof(fc)));
    if (c->nmax > 0) {
        if (!c->f_hall_bare) CUCHK(cudaMalloc(&c->f_hall_bare, nha * sizeof(fc)));
        if (!c->f_hallo) CUCHK(cudaMalloc(&c->f_hallo, nha * sizeof(fc)));
    }
    const int tpb = 256;
#define RSREC_DOWNCONV(N, SRC, DST)                                          \
    k_zc_to_fc<<<(int)(((N) + tpb - 1) / tpb), tpb>>>((N), (SRC), (DST));    \
    CUCHK(cudaGetLastError());
    RSREC_DOWNCONV(nee, c->d_ee_bare, c->f_ee_bare)
    RSREC_DOWNCONV(nee, c->d_eeo, c->f_eeo)
    RSREC_DOWNCONV(nons, c->d_hons, c->f_hons)
    if (c->nmax > 0) {
        RSREC_DOWNCONV(nha, c->d_hall_bare, c->f_hall_bare)
        RSREC_DOWNCONV(nha, c->d_hallo, c->f_hallo)
    }
#undef RSREC_DOWNCONV
    c->f_hoh_ver = c->ham_ver;
    return 0;
}

/* Templated multi-state engine (block-ELL path only). Moment buffer layout:
 * dmu has 2*lld+2 slots, each slot a (bb x ns) block (state-major); the host
 * combine writes mu as (bb, nmom, ns) state-major.                          */
template <class CT, class RT>
static int cheb_engine(rsrec_ctx *c, const void *psi0_h, int ns, int lld,
                       double a, double b, void *mu_h,
                       const CT *bt, const CT *bh, RT inva, RT bsc,
                       const CT *h_bare = nullptr, const CT *h_bare_imp = nullptr,
                       const CT *h_eeo = nullptr, const CT *h_eeo_imp = nullptr,
                       const CT *h_hons = nullptr) {
    const int nb = c->nb, kk = c->kk;
    const size_t bb = (size_t)nb * nb;
    const size_t ld = (size_t)nb * kk;
    const int nrhs = nb * ns;
    const size_t nfield = ld * (size_t)nrhs;
    const size_t nmom = 2 * (size_t)lld + 2;
    const bool do_hoh = c->have_hoh;

    size_t freeb, totalb;
    cudaMemGetInfo(&freeb, &totalb);
    /* p0/p1/p2 (+ 3 hoh scratch when active) + stage(zc) + dmu */
    if ((do_hoh ? 6 : 3) * nfield * sizeof(CT) + nfield * sizeof(zc) +
        nmom * bb * ns * sizeof(CT) > freeb * 9 / 10)
        FAIL("chebyshev_moments: state batch does not fit; reduce nstates");

    CT *p0, *p1, *p2, *dmu;
    CT *t_hoh = nullptr, *e_hoh = nullptr, *o_hoh = nullptr;
    zc *stage;
    CUCHK(cudaMalloc(&p0, nfield * sizeof(CT)));
    CUCHK(cudaMalloc(&p1, nfield * sizeof(CT)));
    CUCHK(cudaMalloc(&p2, nfield * sizeof(CT)));
    CUCHK(cudaMalloc(&stage, nfield * sizeof(zc)));
    CUCHK(cudaMalloc(&dmu, nmom * bb * ns * sizeof(CT)));
    if (do_hoh) {
        CUCHK(cudaMalloc(&t_hoh, nfield * sizeof(CT)));
        CUCHK(cudaMalloc(&e_hoh, nfield * sizeof(CT)));
        CUCHK(cudaMalloc(&o_hoh, nfield * sizeof(CT)));
    }

    const int tpb = 256;
    const int nbl = (int)((nfield + tpb - 1) / tpb);
    CUCHK(cudaMemcpy(stage, psi0_h, nfield * sizeof(zc),
                     cudaMemcpyHostToDevice));
    k_pack<CT><<<nbl, tpb>>>(nfield, nb, nb, kk, ld, stage, p0);
    CUCHK(cudaGetLastError());

    const CT one = {(RT)1, (RT)0};
    const CT zero = {(RT)0, (RT)0};
    const long long scol = (long long)ld * nb;   /* per-state column stride  */
    const long long smu = (long long)bb;
    auto moments = [&](const CT *A, const CT *B, size_t slot) -> int {
        CBCHK(blas_gemm_sb(c->blas, CUBLAS_OP_C, CUBLAS_OP_N, nb, nb,
                           (int)ld, &one, A, (int)ld, scol, B, (int)ld,
                           scol, &zero, dmu + slot * bb * ns, nb, smu, ns));
        return 0;
    };
    /* Matvec dispatch: the two-sweep hoh apply when active, otherwise the
     * single-operator step. Both run at the engine precision CT/RT.         */
    auto matvec = [&](const CT *x1, const CT *x0, CT *yv, RT al, RT be) -> int {
        if (do_hoh)
            return step_apply_hoh<CT, RT>(c, h_bare, h_bare_imp, h_eeo,
                                          h_eeo_imp, h_hons, t_hoh, e_hoh,
                                          o_hoh, x1, x0, yv, nrhs, al, be,
                                          inva, bsc);
        return step_apply<CT, RT>(c, bt, bh, 1, x1, x0, yv, nrhs, al, be,
                                  inva, bsc);
    };

    if (moments(p0, p0, 0)) return 1;
    if (matvec(p0, nullptr, p1, (RT)1, (RT)0)) return 1;
    if (moments(p0, p1, 1)) return 1;
    for (int ll = 1; ll <= lld; ++ll) {
        if (matvec(p1, p0, p2, (RT)2, (RT)-1)) return 1;
        if (moments(p1, p1, (size_t)(2 * ll))) return 1;
        if (moments(p2, p1, (size_t)(2 * ll + 1))) return 1;
        CT *t = p0; p0 = p1; p1 = p2; p2 = t;
    }

    std::vector<CT> h(nmom * bb * ns);
    CUCHK(cudaMemcpy(h.data(), dmu, h.size() * sizeof(CT),
                     cudaMemcpyDeviceToHost));
    cudaFree(p0); cudaFree(p1); cudaFree(p2); cudaFree(stage); cudaFree(dmu);
    if (t_hoh) cudaFree(t_hoh);
    if (e_hoh) cudaFree(e_hoh);
    if (o_hoh) cudaFree(o_hoh);
    cplx *mu = (cplx *)mu_h;                       /* (bb, nmom, ns)         */
    for (int s = 0; s < ns; ++s) {
        cplx *mus = mu + (size_t)s * bb * nmom;
        const CT *h1 = h.data() + 0 * bb * ns + (size_t)s * bb;
        const CT *h2 = h.data() + 1 * bb * ns + (size_t)s * bb;
        for (size_t i = 0; i < bb; ++i) {
            mus[i] = cplx(h1[i].x, h1[i].y);
            mus[bb + i] = cplx(h2[i].x, h2[i].y);
        }
        for (int ll = 1; ll <= lld; ++ll) {
            const CT *d1 = h.data() + (size_t)(2 * ll) * bb * ns +
                           (size_t)s * bb;
            const CT *d2 = h.data() + (size_t)(2 * ll + 1) * bb * ns +
                           (size_t)s * bb;
            for (size_t i = 0; i < bb; ++i) {
                mus[(size_t)(2 * ll) * bb + i] =
                    2.0 * cplx(d1[i].x, d1[i].y) - mus[i];
                mus[(size_t)(2 * ll + 1) * bb + i] =
                    2.0 * cplx(d2[i].x, d2[i].y) - mus[bb + i];
            }
        }
    }
    return 0;
}

/* Single-state structured-backend Chebyshev moments (nrhs == nb): explicit
 * recurrence using step64 (FFT stencil + correction) and gram64 reductions.
 * psi0 is one (nb, nb, kk) state; mu one (bb, 2*lld+2) slot block.          */
static int cheb_struct_one(rsrec_ctx *c, const void *psi0, int lld, double a,
                           double b, void *mu) {
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);
    const size_t ld = (size_t)nb * c->kk;
    const size_t nmom = 2 * (size_t)lld + 2;
    zc *p0 = c->d_s0, *p1 = c->d_s1, *p2 = c->d_s2, *tmp = c->d_s3;
    zc *stage, *dmu;
    CUCHK(cudaMalloc(&stage, nf * sizeof(zc)));
    CUCHK(cudaMalloc(&dmu, nmom * bb * sizeof(zc)));
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    CUCHK(cudaMemcpy(stage, psi0, nf * sizeof(zc), cudaMemcpyHostToDevice));
    k_pack<zc><<<nbl, tpb>>>(nf, nb, nb, c->kk, ld, stage, p0);
    CUCHK(cudaGetLastError());
    if (gram64(c, p0, nb, p0, nb, dmu)) { cudaFree(stage); cudaFree(dmu);
        return 1; }
    if (step64(c, 0, p0, nullptr, p1, tmp, nb, 1.0, 0.0, 1.0 / a, b)) {
        cudaFree(stage); cudaFree(dmu); return 1; }
    if (gram64(c, p0, nb, p1, nb, dmu + bb)) { cudaFree(stage); cudaFree(dmu);
        return 1; }
    for (int ll = 1; ll <= lld; ++ll) {
        if (step64(c, 0, p1, p0, p2, tmp, nb, 2.0, -1.0, 1.0 / a, b)) {
            cudaFree(stage); cudaFree(dmu); return 1; }
        if (gram64(c, p1, nb, p1, nb, dmu + (size_t)(2 * ll) * bb)) {
            cudaFree(stage); cudaFree(dmu); return 1; }
        if (gram64(c, p2, nb, p1, nb, dmu + (size_t)(2 * ll + 1) * bb)) {
            cudaFree(stage); cudaFree(dmu); return 1; }
        zc *t = p0; p0 = p1; p1 = p2; p2 = t;
    }
    std::vector<cplx> h(nmom * bb);
    CUCHK(cudaMemcpy(h.data(), dmu, h.size() * sizeof(zc),
                     cudaMemcpyDeviceToHost));
    cudaFree(stage); cudaFree(dmu);
    cplx *mu_ = (cplx *)mu;
    for (size_t i = 0; i < 2 * bb; ++i) mu_[i] = h[i];
    for (int ll = 1; ll <= lld; ++ll)
        for (size_t i = 0; i < bb; ++i) {
            mu_[(size_t)(2 * ll) * bb + i] =
                2.0 * h[(size_t)(2 * ll) * bb + i] - h[i];
            mu_[(size_t)(2 * ll + 1) * bb + i] =
                2.0 * h[(size_t)(2 * ll + 1) * bb + i] - h[bb + i];
        }
    return 0;
}

extern "C" int rsrec_chebyshev_moments_batch(rsrec_ctx *c, const void *psi0,
                                             int nstates, int lld, double a,
                                             double b, void *mu) {
    if (!c->have_h) FAIL("chebyshev_moments: Hamiltonian not set");
    if (nstates < 1) FAIL("chebyshev_moments: nstates < 1");
    /* hoh: two-sweep block-ELL path. The structured/FFT backend assumes an
     * ee-only Hamiltonian, so hoh always uses the block-ELL engine. The
     * operators are stored unscaled; the engine's combine applies (1/a, -b),
     * hence inva = 1/a and bsc = b are passed (NOT 1, 0 as in the prescaled
     * non-hoh fp32 path). fp32 by default for speed; fp64 if requested.     */
    if (c->have_hoh) {
        if (c->cheb_prec == 0) {
            if (ensure_f32_hoh(c)) return 1;
            return cheb_engine<fc, float>(
                c, psi0, nstates, lld, a, b, mu, c->f_ee_bare, c->f_hall_bare,
                (float)(1.0 / a), (float)b, c->f_ee_bare, c->f_hall_bare,
                c->f_eeo, c->f_hallo, c->f_hons);
        }
        return cheb_engine<zc, double>(
            c, psi0, nstates, lld, a, b, mu, c->d_ee_bare, c->d_hall_bare,
            1.0 / a, b, c->d_ee_bare, c->d_hall_bare, c->d_eeo, c->d_hallo,
            c->d_hons);
    }
    if (c->cheb_prec == 0 && !c->use_struct) {
        if (ensure_f32_ham(c, a, b)) return 1;
        return cheb_engine<fc, float>(c, psi0, nstates, lld, a, b, mu,
                                      c->f_ee, c->f_hall, 1.0f, 0.0f);
    }
    if (!c->use_struct)
        return cheb_engine<zc, double>(c, psi0, nstates, lld, a, b, mu,
                                       c->d_ee, c->d_hall, 1.0 / a, b);
    /* structured backend (fp64): one FFT-stencil sweep per state. The
     * stencil apply is fixed at nrhs == nb, so states are looped rather than
     * stacked as RHS columns; the moments are identical state-by-state.     */
    const size_t nf = fieldsz(c);
    const size_t per_mu = (size_t)c->nb * c->nb * (2 * (size_t)lld + 2);
    const cplx *psi = (const cplx *)psi0;
    cplx *mu_ = (cplx *)mu;
    for (int s = 0; s < nstates; ++s)
        if (cheb_struct_one(c, psi + (size_t)s * nf, lld, a, b,
                            mu_ + (size_t)s * per_mu)) return 1;
    return 0;
}

extern "C" int rsrec_chebyshev_moments(rsrec_ctx *c, const void *psi0,
                                       int lld, double a, double b,
                                       void *mu) {
    return rsrec_chebyshev_moments_batch(c, psi0, 1, lld, a, b, mu);
}

/* --------------------------------------------------------------------------
 * API: block Lanczos (fp64, site-major: each reduction/right-multiply is one
 * big Zgemm; host zheev per step like crecal_b)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_block_lanczos(rsrec_ctx *c, const void *psi0_, int lld,
                                   void *a_b_, void *b2_b_) {
    if (!c->have_h) FAIL("block_lanczos: Hamiltonian not set");
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);
    const size_t ld = (size_t)nb * c->kk;
    cplx *a_b = (cplx *)a_b_, *b2_b = (cplx *)b2_b_;

    zc *psi = c->d_s0, *pmn = c->d_s1, *hpsi = c->d_s2, *tmp = c->d_s3;
    zc *gt;
    CUCHK(cudaMalloc(&gt, nf * sizeof(zc)));
    zc *dAn = c->d_blk, *dSum = c->d_blk + bb;
    zc *dB = c->d_blk + 2 * bb, *dBi = c->d_blk + 3 * bb;

    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    CUCHK(cudaMemcpy(gt, psi0_, nf * sizeof(zc), cudaMemcpyHostToDevice));
    k_pack<zc><<<nbl, tpb>>>(nf, nb, nb, c->kk, ld, gt, psi);
    CUCHK(cudaGetLastError());
    CUCHK(cudaMemset(pmn, 0, nf * sizeof(zc)));

    std::vector<cplx> sum_b(bb, cplx(0, 0));
    for (int i = 0; i < nb; ++i) sum_b[i + (size_t)nb * i] = 1.0;
    std::memset(a_b, 0, bb * lld * sizeof(cplx));
    std::memset(b2_b, 0, bb * lld * sizeof(cplx));

    for (int ll = 0; ll < lld - 1; ++ll) {
        if (step64(c, 0, psi, nullptr, hpsi, tmp, nb, 1.0, 0.0, 1.0, 0.0))
            return 1;
        if (gram64(c, psi, nb, hpsi, nb, dAn)) return 1;          /* A_n    */
        CUCHK(cudaMemcpy(a_b + (size_t)ll * bb, dAn, bb * sizeof(zc),
                         cudaMemcpyDeviceToHost));
        k_axpby<zc, double><<<nbl, tpb>>>(nf, 1.0, hpsi, -1.0, pmn);
        CUCHK(cudaGetLastError());                       /* pmn = hpsi - pmn */
        std::memcpy(b2_b + (size_t)ll * bb, sum_b.data(), bb * sizeof(cplx));

        /* pmn -= psi * A_n: ONE (ld x nb x nb) gemm, beta = 1               */
        CBCHK(blas_gemm(c->blas, CUBLAS_OP_N, CUBLAS_OP_N, (int)ld, nb, nb,
                        &Z_MONE, psi, (int)ld, dAn, nb, &Z_ONE, pmn,
                        (int)ld));
        if (gram64(c, pmn, nb, pmn, nb, dSum)) return 1;          /* B^2    */
        CUCHK(cudaMemcpy(sum_b.data(), dSum, bb * sizeof(zc),
                         cudaMemcpyDeviceToHost));
        if (host_eig_sqrt(c, dSum, dB, dBi)) return 1;

        /* psi_{n+1} = pmn * Bi (gemm into gt); pmn = psi_n * B; swap        */
        CBCHK(blas_gemm(c->blas, CUBLAS_OP_N, CUBLAS_OP_N, (int)ld, nb, nb,
                        &Z_ONE, pmn, (int)ld, dBi, nb, &Z_ZERO, gt,
                        (int)ld));
        CBCHK(blas_gemm(c->blas, CUBLAS_OP_N, CUBLAS_OP_N, (int)ld, nb, nb,
                        &Z_ONE, psi, (int)ld, dB, nb, &Z_ZERO, pmn,
                        (int)ld));
        zc *t = psi; psi = gt; gt = t;        /* psi <- new; gt <- old psi   */
    }
    std::memcpy(b2_b + (size_t)(lld - 1) * bb, sum_b.data(),
                bb * sizeof(cplx));
    /* free whichever rotated buffer is the extra allocation (not s0..s3)    */
    zc *extra = (psi != c->d_s0 && psi != c->d_s1 && psi != c->d_s2 &&
                 psi != c->d_s3) ? psi : gt;
    cudaFree(extra);
    return 0;
}

/* --------------------------------------------------------------------------
 * API: scalar Lanczos (fp64, site-major, all nb chains batched as columns)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_scalar_lanczos(rsrec_ctx *c, int site_j, int lld,
                                    double *a_out, double *b2_out) {
    if (!c->have_h) FAIL("scalar_lanczos: Hamiltonian not set");
    if (site_j < 1 || site_j > c->kk) FAIL("scalar_lanczos: bad site");
    const int nb = c->nb;
    const size_t nf = fieldsz(c), ld = (size_t)nb * c->kk;

    zc *psi = c->d_s0, *pmn = c->d_s1, *hpsi = c->d_s2, *tmp = c->d_s3;
    CUCHK(cudaMemset(psi, 0, nf * sizeof(zc)));
    CUCHK(cudaMemset(pmn, 0, nf * sizeof(zc)));
    {   /* psi(row l, col l, atom site_j) = 1 (site-major)                   */
        std::vector<cplx> col(nb, cplx(0, 0));
        for (int l = 0; l < nb; ++l) {
            col.assign(nb, cplx(0, 0));
            col[l] = 1.0;
            CUCHK(cudaMemcpy(psi + (size_t)nb * (site_j - 1) + ld * l,
                             col.data(), nb * sizeof(zc),
                             cudaMemcpyHostToDevice));
        }
    }
    double *dA = c->d_col, *dS = c->d_col + nb;
    zc *dG = c->d_blk;                       /* nb x nb gram scratch          */
    std::vector<double> acol(nb), s2(nb), summ(nb, 1.0);
    std::memset(a_out, 0, sizeof(double) * lld * nb);
    std::memset(b2_out, 0, sizeof(double) * lld * nb);

    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    const int ndg = (nb + 255) / 256;

    for (int ll = 0; ll < lld - 1; ++ll) {
        if (step64(c, 0, psi, nullptr, hpsi, tmp, nb, 1.0, 0.0, 1.0, 0.0))
            return 1;
        /* a_col = Re diag(psi^H H psi): one gram gemm + diagonal pick        */
        if (gram64(c, psi, nb, hpsi, nb, dG)) return 1;
        k_diag_real<<<ndg, tpb>>>(nb, dG, dA);
        CUCHK(cudaGetLastError());
        k_axpby<zc, double><<<nbl, tpb>>>(nf, 1.0, hpsi, 1.0, pmn);
        CUCHK(cudaMemcpy(acol.data(), dA, nb * sizeof(double),
                         cudaMemcpyDeviceToHost));
        for (int col = 0; col < nb; ++col) {
            a_out[ll + (size_t)lld * col] = acol[col];
            b2_out[ll + (size_t)lld * col] = summ[col];
            acol[col] = -acol[col];
        }
        CUCHK(cudaMemcpy(dA, acol.data(), nb * sizeof(double),
                         cudaMemcpyHostToDevice));
        k_col_axpy<<<nbl, tpb>>>(ld, nb, dA, psi, pmn);  /* pmn -= a psi    */
        /* B^2 = Re diag(pmn^H pmn): one gram gemm + diagonal pick           */
        if (gram64(c, pmn, nb, pmn, nb, dG)) return 1;
        k_diag_real<<<ndg, tpb>>>(nb, dG, dS);
        CUCHK(cudaGetLastError());
        CUCHK(cudaMemcpy(s2.data(), dS, nb * sizeof(double),
                         cudaMemcpyDeviceToHost));
        k_col_shift<<<nbl, tpb>>>(ld, nb, dS, psi, pmn);
        CUCHK(cudaGetLastError());
        summ = s2;
    }
    for (int col = 0; col < nb; ++col)
        b2_out[(lld - 1) + (size_t)lld * col] = summ[col];
    return 0;
}

/* --------------------------------------------------------------------------
 * API: stochastic conductivity moments (fp64, site-major; the (n, m)
 * contraction is one gram gemm per stored left state)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_stochastic_moments(rsrec_ctx *c, const void *psiref_,
                                        int lld, double a, double b,
                                        void *mu_) {
    if (!c->have_h) FAIL("stochastic_moments: Hamiltonian not set");
    if (!c->have_v) FAIL("stochastic_moments: velocity operators not set");
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);
    const size_t ld = (size_t)nb * c->kk;

    size_t freeb, totalb;
    cudaMemGetInfo(&freeb, &totalb);
    size_t need = ((size_t)lld * nf + nf + bb * (size_t)lld * lld)
                  * sizeof(zc);
    if (need > freeb * 9 / 10)
        FAIL("stochastic_moments: left states do not fit on the device; "
             "reduce lld or split over reference vectors / devices");

    zc *left, *dmu, *R;
    CUCHK(cudaMalloc(&left, (size_t)lld * nf * sizeof(zc)));
    CUCHK(cudaMalloc(&dmu, bb * (size_t)lld * lld * sizeof(zc)));
    CUCHK(cudaMalloc(&R, nf * sizeof(zc)));
    CUCHK(cudaMemset(dmu, 0, bb * (size_t)lld * lld * sizeof(zc)));

    zc *w0 = c->d_s0, *w1 = c->d_s1, *w2 = c->d_s2, *tmp = c->d_s3;
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    CUCHK(cudaMemcpy(R, psiref_, nf * sizeof(zc), cudaMemcpyHostToDevice));
    k_pack<zc><<<nbl, tpb>>>(nf, nb, nb, c->kk, ld, R, w1);    /* w1 = psiref       */
    CUCHK(cudaGetLastError());

    /* left states L_m = T_{m-1}(H~)|psiref>, stored device-resident         */
    CUCHK(cudaMemcpy(left, w1, nf * sizeof(zc), cudaMemcpyDeviceToDevice));
    for (int m = 2; m <= lld; ++m) {
        if (m == 2) {
            zc *t = w0; w0 = w1; w1 = t;
            if (step64(c, 0, w0, nullptr, w1, tmp, nb, 1.0, 0.0, 1.0 / a, b))
                return 1;
        } else {
            if (step64(c, 0, w1, w0, w2, tmp, nb, 2.0, -1.0, 1.0 / a, b))
                return 1;
            zc *t = w0; w0 = w1; w1 = w2; w2 = t;
        }
        CUCHK(cudaMemcpy(left + (size_t)(m - 1) * nf, w1, nf * sizeof(zc),
                         cudaMemcpyDeviceToDevice));
    }

    /* right recursion v_n = T_{n-1}(H~) V_b|psiref>; left slot 0 = psiref   */
    if (step64(c, 2, left, nullptr, w0, tmp, nb, 1.0, 0.0, 1.0, 0.0))
        return 1;                                /* w0 = V_b psiref          */
    zc *v0 = w0, *v1 = w1, *v2 = w2;
    for (int n = 1; n <= lld; ++n) {
        if (n == 1) {
            CUCHK(cudaMemcpy(v1, v0, nf * sizeof(zc),
                             cudaMemcpyDeviceToDevice));
        } else if (n == 2) {
            CUCHK(cudaMemcpy(v0, v1, nf * sizeof(zc),
                             cudaMemcpyDeviceToDevice));
            if (step64(c, 0, v0, nullptr, v1, tmp, nb, 1.0, 0.0, 1.0 / a, b))
                return 1;
        } else {
            if (step64(c, 0, v1, v0, v2, tmp, nb, 2.0, -1.0, 1.0 / a, b))
                return 1;
            zc *t = v0; v0 = v1; v1 = v2; v2 = t;
        }
        if (step64(c, 1, v1, nullptr, R, tmp, nb, 1.0, 0.0, 1.0, 0.0))
            return 1;                            /* R = V_a v_n              */
        /* (n, m) contraction for ALL m at once: one strided-batched gemm.
         * batch m -> C_m = left_m^H R (nb x nb). left states are contiguous
         * (stride nf); R is shared (strideB = 0). The dmu slot (n-1, m) sits
         * at base bb*(n-1) with per-m stride bb*lld, so the batched output
         * lands directly in place -- lld launches collapse to one.          */
        CBCHK(blas_gemm_sb(c->blas, CUBLAS_OP_C, CUBLAS_OP_N, nb, nb, (int)ld,
                           &Z_ONE, left, (int)ld, (long long)nf,
                           R, (int)ld, 0,
                           &Z_ZERO, dmu + bb * (size_t)(n - 1), nb,
                           (long long)bb * lld, lld));
    }
    CUCHK(cudaMemcpy(mu_, dmu, bb * (size_t)lld * lld * sizeof(zc),
                     cudaMemcpyDeviceToHost));
    cudaFree(left); cudaFree(dmu); cudaFree(R);
    return 0;
}

/* --------------------------------------------------------------------------
 * Chebyshev Green-function / DOS reconstruction (GPU port of
 * chebyshev_green in bands.f90; identical math and Jackson convention).
 *
 *   g0(:,:,ie,n) = sum_i mu(:,:,i,n) F(i,ie),
 *   F(i,ie) = g_J(i) c_i (-i) e^{-i(i-1)acos(x_ie)} / sqrt(a^2-(E_ie-b)^2)
 *
 * F is built by one kernel; the contraction over moments is ONE
 * cublas[Z|C]gemmStridedBatched over the local atoms (strideB = 0 -- all
 * atoms share F). Precision follows rsrec_set_precision: fp32 (default,
 * cgemm3m, consistent with the fp32 moment engine) or fp64
 * (bit-comparable with the Fortran loop). Atom chunking keeps device
 * memory bounded for large nv/lld.
 * ------------------------------------------------------------------------ */
__global__ void k_build_F64(int N, int nv, const double *__restrict__ ene,
                            double a, double b, zc *__restrict__ F) {
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= (size_t)N * nv) return;
    const double PI = 3.14159265358979323846;
    int i = (int)(idx % N), ie = (int)(idx / N);
    double x = (ene[ie] - b) / a;
    double th = acos(x);
    double pref = rsqrt(a * a - (ene[ie] - b) * (ene[ie] - b));
    double thl = PI * i / (N + 1.0);
    double cot = cos(PI / (N + 1.0)) / sin(PI / (N + 1.0));
    double gj = ((N - i + 1) * cos(thl) + sin(thl) * cot) / (N + 1.0);
    double cc = (i == 0) ? 1.0 : 2.0;
    double ang = i * th;                       /* (-i) e^{-i ang}           */
    F[idx] = make_cuDoubleComplex(-gj * cc * pref * sin(ang),
                                  -gj * cc * pref * cos(ang));
}

__global__ void k_z2c(size_t n, const zc *__restrict__ in,
                      fc *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = make_cuFloatComplex((float)in[i].x, (float)in[i].y);
}

__global__ void k_c2z(size_t n, const fc *__restrict__ in,
                      zc *__restrict__ out) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = make_cuDoubleComplex(in[i].x, in[i].y);
}

extern "C" int rsrec_chebyshev_dos(rsrec_ctx *c, const void *mu_, int n_mom,
                                   int natoms, const double *ene, int nv,
                                   double a, double b, void *g0_) {
    if (n_mom < 1 || natoms < 1 || nv < 1) FAIL("chebyshev_dos: bad sizes");
    const int bb = c->nb * c->nb;
    const cplx *mu = (const cplx *)mu_;
    cplx *g0 = (cplx *)g0_;
    const bool f32 = (c->cheb_prec == 0);
    const int tpb = 256;

    /* F (always built in fp64; converted once if needed)                   */
    double *d_ene; zc *d_F64 = nullptr; fc *d_F32 = nullptr;
    CUCHK(cudaMalloc(&d_ene, nv * sizeof(double)));
    CUCHK(cudaMemcpy(d_ene, ene, nv * sizeof(double),
                     cudaMemcpyHostToDevice));
    size_t nF = (size_t)n_mom * nv;
    CUCHK(cudaMalloc(&d_F64, nF * sizeof(zc)));
    k_build_F64<<<(int)((nF + tpb - 1) / tpb), tpb>>>(n_mom, nv, d_ene,
                                                      a, b, d_F64);
    CUCHK(cudaGetLastError());
    if (f32) {
        CUCHK(cudaMalloc(&d_F32, nF * sizeof(fc)));
        k_z2c<<<(int)((nF + tpb - 1) / tpb), tpb>>>(nF, d_F64, d_F32);
        CUCHK(cudaGetLastError());
    }

    /* atom chunking against free device memory                             */
    size_t per_atom = (size_t)bb * n_mom * sizeof(zc)        /* mu          */
                    + (size_t)bb * nv * sizeof(zc)           /* g0          */
                    + (f32 ? ((size_t)bb * n_mom + (size_t)bb * nv)
                                 * sizeof(fc) : 0);
    size_t freeb, totalb;
    cudaMemGetInfo(&freeb, &totalb);
    int chunk = (int)std::min<size_t>((size_t)natoms,
                                      std::max<size_t>(1, (freeb * 7 / 10)
                                                          / per_atom));
    zc *d_mu, *d_g0; fc *d_muf = nullptr, *d_g0f = nullptr;
    CUCHK(cudaMalloc(&d_mu, (size_t)bb * n_mom * chunk * sizeof(zc)));
    CUCHK(cudaMalloc(&d_g0, (size_t)bb * nv * chunk * sizeof(zc)));
    if (f32) {
        CUCHK(cudaMalloc(&d_muf, (size_t)bb * n_mom * chunk * sizeof(fc)));
        CUCHK(cudaMalloc(&d_g0f, (size_t)bb * nv * chunk * sizeof(fc)));
    }
    for (int n0 = 0; n0 < natoms; n0 += chunk) {
        int nc = std::min(chunk, natoms - n0);
        size_t nmu = (size_t)bb * n_mom * nc, ng = (size_t)bb * nv * nc;
        CUCHK(cudaMemcpy(d_mu, mu + (size_t)bb * n_mom * n0,
                         nmu * sizeof(zc), cudaMemcpyHostToDevice));
        if (f32) {
            k_z2c<<<(int)((nmu + tpb - 1) / tpb), tpb>>>(nmu, d_mu, d_muf);
            CUCHK(cudaGetLastError());
            CBCHK(blas_gemm_sb(c->blas, CUBLAS_OP_N, CUBLAS_OP_N, bb, nv,
                               n_mom, &C_ONE, d_muf, bb,
                               (long long)bb * n_mom, d_F32, n_mom, 0,
                               &C_ZERO, d_g0f, bb, (long long)bb * nv, nc));
            k_c2z<<<(int)((ng + tpb - 1) / tpb), tpb>>>(ng, d_g0f, d_g0);
            CUCHK(cudaGetLastError());
        } else {
            CBCHK(blas_gemm_sb(c->blas, CUBLAS_OP_N, CUBLAS_OP_N, bb, nv,
                               n_mom, &Z_ONE, d_mu, bb,
                               (long long)bb * n_mom, d_F64, n_mom, 0,
                               &Z_ZERO, d_g0, bb, (long long)bb * nv, nc));
        }
        CUCHK(cudaMemcpy(g0 + (size_t)bb * nv * n0, d_g0, ng * sizeof(zc),
                         cudaMemcpyDeviceToHost));
    }
    for (void *p : {(void *)d_ene, (void *)d_F64, (void *)d_F32,
                    (void *)d_mu, (void *)d_g0, (void *)d_muf,
                    (void *)d_g0f})
        if (p) cudaFree(p);
    return 0;
}
