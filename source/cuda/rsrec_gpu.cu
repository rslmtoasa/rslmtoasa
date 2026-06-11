/* ===========================================================================
 * rsrec_gpu.cu -- CUDA implementation of rsrec.h
 *
 * Mirrors rsrec_cpu.cpp routine-for-routine (the CPU file is validated to
 * machine precision against a NumPy transcription of recursion.f90). All
 * hot paths are cuBLAS:
 *
 *   op_apply   -> cublasZgemmBatched, one batch per neighbour shell over
 *                 precomputed pointer arrays. The (atom, neighbour) pair
 *                 lists are static, so A-pointers (Hamiltonian blocks) are
 *                 built once at upload; x/y pointers are refreshed by a
 *                 trivial kernel per call. Within a shell every atom owns
 *                 its output block exactly once -> race-free beta=1
 *                 accumulation. The (y - b x)/a shift is one fused kernel.
 *   block_dot  -> sum_k A_k^H B_k as cublasZgemmStridedBatched (per-atom
 *                 C_k = A_k^H B_k) + cublasZgemv against a ones vector
 *                 (reduction over atoms). Two library calls, no atomics.
 *   right-muls -> cublasZgemmStridedBatched with strideB = 0 (same nb x nb
 *                 matrix for all atoms): psi*B, psi*B^-1, pmn -= psi*A_n.
 *   eig sqrt   -> the nb x nb B^2 is copied to the HOST and solved with
 *                 LAPACK zheev, exactly like crecal_b (bit-comparable with
 *                 the CPU path; ~5 KB of traffic per recursion step).
 *
 * A fused custom matvec kernel is kept as a fallback for memory-constrained
 * cards (export RSREC_MATVEC=fused); it avoids the pointer arrays at the
 * cost of not using cuBLAS. Default is the cuBLAS path.
 *
 * Pointer-array memory: 8 B per (atom,shell) block per operator + 16 B for
 * the shared x/y arrays, i.e. ~24 B per stored Hamiltonian block -- small
 * compared to the wavefunction fields (16 B * nb^2 per atom each).
 *
 * Everything is double complex; recurrences stay device-resident, only
 * nb x nb moment/coefficient blocks cross PCIe per step.
 *
 * Build:  nvcc -O3 -arch=native -Xcompiler -fPIC -shared rsrec_gpu.cu \
 *              -o librsrec_gpu.so -lcublas -llapack
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

/* Custom atomicAdd for double (CUDA 12.8 compatibility) */
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ < 600
__device__ double atomicAddDouble(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(
                            __longlong_as_double(assumed) + val));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#else
__device__ inline double atomicAddDouble(double *address, double val) {
    return atomicAdd(address, val);
}
#endif

/* --------------------------------------------------------------------------
 * Context
 * ------------------------------------------------------------------------ */
struct rsrec_ctx {
    int kk, nb, nnmax, ntype, nmax, device;
    bool fused_matvec = false;            /* RSREC_MATVEC=fused             */

    /* Device-resident operators (lsham pre-folded into shell-1 blocks)     */
    zc *d_ee = nullptr, *d_hall = nullptr, *d_va = nullptr, *d_vb = nullptr;
    int *d_nn = nullptr, *d_iz = nullptr;
    bool have_h = false, have_v = false;

    /* Static (atom, neighbour) pair lists, compacted per shell, and the
     * per-operator A-pointer arrays over the same ordering.                 */
    int n_shells = 0;
    std::vector<int> sh_off;              /* host: shell start offsets       */
    int *d_pair_atom = nullptr;           /* [nnz] 0-based atom k            */
    int *d_pair_nbr = nullptr;            /* [nnz] 0-based neighbour         */
    const zc **d_Aptr_H = nullptr;        /* [nnz] Hamiltonian block ptrs    */
    const zc **d_Aptr_VA = nullptr;       /* [nnz] velocity blocks (lazy)    */
    const zc **d_Aptr_VB = nullptr;
    const zc **d_Bptr = nullptr;          /* [nnz] x blocks (per call)       */
    zc **d_Cptr = nullptr;                /* [nnz] y blocks (per call)       */
    size_t nnz = 0;

    /* Device scratch fields, each (nb, nb, kk)                              */
    zc *d_s0 = nullptr, *d_s1 = nullptr, *d_s2 = nullptr, *d_s3 = nullptr;
    zc *d_bd = nullptr;                   /* (nb*nb, kk) block-dot scratch   */
    zc *d_ones = nullptr;                 /* (kk) ones for the reduction     */
    zc *d_blk = nullptr;                  /* a few nb x nb accumulators      */
    double *d_col = nullptr;              /* per-column scalars              */
    cublasHandle_t blas = nullptr;

    /* ---- structured (cuFFT stencil + batched-GEMM correction) path ----- */
    bool use_struct = false;
    int t_ref = 0, nshell_ref = 0;
    int N[3] = {0, 0, 0}, Np[3] = {0, 0, 0};
    size_t Npc = 0;                       /* padded cell count               */
    std::vector<int> roff;                /* (3, nshell_ref) host            */
    int *d_gcell = nullptr;               /* atom -> padded cell index       */
    zc *d_gin = nullptr, *d_gw = nullptr; /* channel-major grids (nb^2 ch)   */
    zc *d_Hq = nullptr, *d_VAq = nullptr, *d_VBq = nullptr;  /* H(q) per op */
    zc *d_negW = nullptr;                 /* -W_ref blocks, 3 ops packed     */
    cufftHandle fft_plan = 0;
    bool have_plan = false;
    /* correction lists, grouped by per-atom slot for race-free beta=1      */
    std::vector<int> cg_off;              /* host group offsets              */
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
    for (zc *p : {c->d_ee, c->d_hall, c->d_va, c->d_vb, c->d_s0, c->d_s1,
                  c->d_s2, c->d_s3, c->d_bd, c->d_ones, c->d_blk})
        if (p) cudaFree(p);
    for (void *p : {(void *)c->d_nn, (void *)c->d_iz, (void *)c->d_pair_atom,
                    (void *)c->d_pair_nbr, (void *)c->d_Aptr_H,
                    (void *)c->d_Aptr_VA, (void *)c->d_Aptr_VB,
                    (void *)c->d_Bptr, (void *)c->d_Cptr,
                    (void *)c->d_gcell, (void *)c->d_gin, (void *)c->d_gw,
                    (void *)c->d_Hq, (void *)c->d_VAq, (void *)c->d_VBq,
                    (void *)c->d_negW, (void *)c->d_c_atom,
                    (void *)c->d_c_nbr, (void *)c->d_cA_H,
                    (void *)c->d_cA_VA, (void *)c->d_cA_VB,
                    (void *)c->d_cB, (void *)c->d_cC})
        if (p) cudaFree(p);
    if (c->have_plan) cufftDestroy(c->fft_plan);
    if (c->d_col) cudaFree(c->d_col);
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
    const char *mv = getenv("RSREC_MATVEC");
    c->fused_matvec = (mv && std::string(mv) == "fused");
    if (cudaSetDevice(device) != cudaSuccess ||
        cublasCreate(&c->blas) != CUBLAS_STATUS_SUCCESS) {
        g_err = "rsrec_create: CUDA/cuBLAS init failed";
        delete c; return nullptr;
    }
    size_t nf = fieldsz(c) * sizeof(zc);
    std::vector<cplx> ones(kk, cplx(1.0, 0.0));
    if (cudaMalloc(&c->d_s0, nf) || cudaMalloc(&c->d_s1, nf) ||
        cudaMalloc(&c->d_s2, nf) || cudaMalloc(&c->d_s3, nf) ||
        cudaMalloc(&c->d_bd, nf) ||
        cudaMalloc(&c->d_ones, kk * sizeof(zc)) ||
        cudaMalloc(&c->d_blk, 8 * (size_t)nb * nb * sizeof(zc)) ||
        cudaMalloc(&c->d_col, 2 * (size_t)nb * sizeof(double)) ||
        cudaMemcpy(c->d_ones, ones.data(), kk * sizeof(zc),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
        g_err = "rsrec_create: out of device memory for scratch fields";
        rsrec_destroy(c); return nullptr;
    }
    return c;
}

/* --------------------------------------------------------------------------
 * Upload + static pointer-array construction
 * ------------------------------------------------------------------------ */
static int build_pairs(rsrec_ctx *c, const int *nn, const int *iz) {
    /* Compact (atom, neighbour) pairs per shell; shell 0 is on-site.       */
    const int kk = c->kk, nnmax = c->nnmax;
    std::vector<int> atom, nbr;
    c->sh_off.assign(1, 0);
    int n_shells = 0;
    for (int s = 0; s < nnmax; ++s) {
        bool any = false;
        for (int k = 0; k < kk; ++k) {
            int nr = nn[k];                       /* nn(k,1) = shell count   */
            if (s >= nr) continue;
            int nb_at = (s == 0) ? k : nn[k + (size_t)kk * s] - 1;
            if (nb_at < 0) continue;
            atom.push_back(k); nbr.push_back(nb_at);
            any = true;
        }
        c->sh_off.push_back((int)atom.size());
        if (any) n_shells = s + 1;
    }
    c->n_shells = n_shells;
    c->sh_off.resize(n_shells + 1);
    c->nnz = atom.size();
    if (c->nnz == 0) FAIL("set_hamiltonian: empty neighbour list");

    CUCHK(cudaMalloc(&c->d_pair_atom, c->nnz * sizeof(int)));
    CUCHK(cudaMalloc(&c->d_pair_nbr, c->nnz * sizeof(int)));
    CUCHK(cudaMemcpy(c->d_pair_atom, atom.data(), c->nnz * sizeof(int),
                     cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(c->d_pair_nbr, nbr.data(), c->nnz * sizeof(int),
                     cudaMemcpyHostToDevice));
    CUCHK(cudaMalloc(&c->d_Bptr, c->nnz * sizeof(zc *)));
    CUCHK(cudaMalloc(&c->d_Cptr, c->nnz * sizeof(zc *)));

    /* Hamiltonian A-pointers over the same ordering (device addresses are
     * known on the host after upload, so build there).                      */
    const size_t bb = (size_t)c->nb * c->nb;
    std::vector<const zc *> Ah(c->nnz);
    size_t idx = 0;
    for (int s = 0; s < n_shells; ++s)
        for (int k = 0; k < kk; ++k) {
            int nr = nn[k];
            if (s >= nr) continue;
            if (s > 0 && nn[k + (size_t)kk * s] - 1 < 0) continue;
            bool imp = (k < c->nmax);
            int slot = imp ? k : (iz[k] - 1);
            const zc *base = imp ? c->d_hall : c->d_ee;
            Ah[idx++] = base + bb * ((size_t)s + (size_t)c->nnmax * slot);
        }
    CUCHK(cudaMalloc(&c->d_Aptr_H, c->nnz * sizeof(zc *)));
    CUCHK(cudaMemcpy(c->d_Aptr_H, Ah.data(), c->nnz * sizeof(zc *),
                     cudaMemcpyHostToDevice));
    return 0;
}

extern "C" int rsrec_set_hamiltonian(rsrec_ctx *c, const void *ee_,
                                     const void *hall_, const void *lsham_,
                                     const int *nn, const int *iz) {
    const cplx *ee = (const cplx *)ee_;
    const cplx *hall = (const cplx *)hall_;
    const cplx *ls = (const cplx *)lsham_;
    if (!ee || !nn || !iz) FAIL("set_hamiltonian: null input");
    if (c->nmax > 0 && !hall) FAIL("set_hamiltonian: nmax>0 but hall is null");
    const int nb = c->nb, nnmax = c->nnmax;
    const size_t bb = (size_t)nb * nb;

    /* Fold lsham into the on-site (shell-1) blocks once, like the CPU.     */
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

    if (!c->d_pair_atom)                  /* static topology: build once    */
        if (build_pairs(c, nn, iz)) return 1;
    c->have_h = true;
    return 0;
}

extern "C" int rsrec_set_velocity(rsrec_ctx *c, const void *va, const void *vb) {
    if (!va || !vb) FAIL("set_velocity: null input");
    if (!c->d_pair_atom) FAIL("set_velocity: call set_hamiltonian first");
    const size_t bb = (size_t)c->nb * c->nb;
    size_t n = bb * c->nnmax * c->ntype * sizeof(zc);
    if (!c->d_va) CUCHK(cudaMalloc(&c->d_va, n));
    if (!c->d_vb) CUCHK(cudaMalloc(&c->d_vb, n));
    CUCHK(cudaMemcpy(c->d_va, va, n, cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(c->d_vb, vb, n, cudaMemcpyHostToDevice));

    /* Velocity A-pointers: per-type blocks for every atom (no hall part).  */
    std::vector<int> h_atom(c->nnz);
    CUCHK(cudaMemcpy(h_atom.data(), c->d_pair_atom, c->nnz * sizeof(int),
                     cudaMemcpyDeviceToHost));
    std::vector<int> h_iz(c->kk);
    CUCHK(cudaMemcpy(h_iz.data(), c->d_iz, c->kk * sizeof(int),
                     cudaMemcpyDeviceToHost));
    std::vector<const zc *> Aa(c->nnz), Ab(c->nnz);
    for (int s = 0, idx = 0; s < c->n_shells; ++s)
        for (int i = c->sh_off[s]; i < c->sh_off[s + 1]; ++i, ++idx) {
            int t = h_iz[h_atom[idx]] - 1;
            size_t off = bb * ((size_t)s + (size_t)c->nnmax * t);
            Aa[idx] = c->d_va + off;
            Ab[idx] = c->d_vb + off;
        }
    if (!c->d_Aptr_VA) CUCHK(cudaMalloc(&c->d_Aptr_VA, c->nnz * sizeof(zc *)));
    if (!c->d_Aptr_VB) CUCHK(cudaMalloc(&c->d_Aptr_VB, c->nnz * sizeof(zc *)));
    CUCHK(cudaMemcpy(c->d_Aptr_VA, Aa.data(), c->nnz * sizeof(zc *),
                     cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(c->d_Aptr_VB, Ab.data(), c->nnz * sizeof(zc *),
                     cudaMemcpyHostToDevice));
    c->have_v = true;
    return 0;
}

/* --------------------------------------------------------------------------
 * Kernels
 * ------------------------------------------------------------------------ */
__device__ __forceinline__ zc zfma(zc a, zc b, zc s) {
    s.x = fma(a.x, b.x, fma(-a.y, b.y, s.x));
    s.y = fma(a.x, b.y, fma(a.y, b.x, s.y));
    return s;
}

/* Refresh the x/y block pointers over the static pair list.                */
__global__ void k_make_ptrs(size_t nnz, size_t bb,
                            const int *__restrict__ atom,
                            const int *__restrict__ nbr,
                            const zc *__restrict__ x, zc *__restrict__ y,
                            const zc **__restrict__ Bp, zc **__restrict__ Cp) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nnz) {
        Bp[i] = x + bb * (size_t)nbr[i];
        Cp[i] = y + bb * (size_t)atom[i];
    }
}

/* y = alpha*x + beta*y (elementwise over a full field)                     */
__global__ void k_axpby(size_t n, double alpha, const zc *__restrict__ x,
                        double beta, zc *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        y[i].x = alpha * x[i].x + beta * y[i].x;
        y[i].y = alpha * x[i].y + beta * y[i].y;
    }
}

/* Fallback fused block-ELL matvec (RSREC_MATVEC=fused): thread block per
 * atom, (nb x nrhs) threads, blocks staged through shared memory.          */
__global__ void k_op_apply(int kk, int nb, int nnmax, int nmax, int nrhs,
                           int which,
                           const zc *__restrict__ ee,
                           const zc *__restrict__ hall,
                           const zc *__restrict__ va,
                           const zc *__restrict__ vb,
                           const int *__restrict__ nn,
                           const int *__restrict__ iz,
                           const zc *__restrict__ x, zc *__restrict__ y,
                           double a, double b) {
    extern __shared__ zc sh[];
    zc *shH = sh, *shX = sh + nb * nb;
    const int k = blockIdx.x;
    if (k >= kk) return;
    const int l = threadIdx.x, col = threadIdx.y;
    const size_t bb = (size_t)nb * nb;
    const bool imp = (which == 0) && (k < nmax);
    const zc *blocks = (which == 1) ? va : (which == 2) ? vb
                       : (imp ? hall : ee);
    const int slot = imp ? k : (iz[k] - 1);
    const zc *base = blocks + bb * (size_t)nnmax * slot;

    zc acc = make_cuDoubleComplex(0.0, 0.0);
    const int nr = nn[k];
    for (int s = 0; s < nr; ++s) {
        int nbr = (s == 0) ? k : nn[k + (size_t)kk * s] - 1;
        if (nbr >= 0) {
            for (int i = l + nb * col; i < nb * nb; i += nb * nrhs)
                shH[i] = base[bb * s + i];
            __syncthreads();
            shX[l + nb * col] = x[(size_t)l + nb * ((size_t)col + nrhs * nbr)];
            __syncthreads();
            for (int m = 0; m < nb; ++m)
                acc = zfma(shH[l + nb * m], shX[m + nb * col], acc);
            __syncthreads();
        }
    }
    const size_t idx = (size_t)l + nb * ((size_t)col + nrhs * k);
    if (b != 0.0 || a != 1.0) {
        zc xk = x[idx];
        acc.x = (acc.x - b * xk.x) / a;
        acc.y = (acc.y - b * xk.y) / a;
    }
    y[idx] = acc;
}

/* Scalar-Lanczos per-column kernels (columns are independent chains)       */
__global__ void k_col_dot(int kk, int nb, const zc *__restrict__ A,
                          const zc *__restrict__ B, double *__restrict__ out,
                          int atoms_per_block) {
    const int col = threadIdx.y, l = threadIdx.x;
    const int k0 = blockIdx.x * atoms_per_block;
    const int k1 = min(kk, k0 + atoms_per_block);
    double s = 0.0;
    for (int k = k0; k < k1; ++k) {
        size_t id = (size_t)l + nb * ((size_t)col + nb * k);
        s += A[id].x * B[id].x + A[id].y * B[id].y;
    }
    atomicAddDouble(&out[col], s);
}

__global__ void k_col_axpy(size_t kk, int nb, const double *__restrict__ alpha,
                           const zc *__restrict__ x, zc *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t n = (size_t)nb * nb * kk;
    if (i < n) {
        int col = (int)((i / nb) % nb);
        y[i].x += alpha[col] * x[i].x;
        y[i].y += alpha[col] * x[i].y;
    }
}

__global__ void k_col_shift(size_t kk, int nb, const double *__restrict__ s2,
                            zc *__restrict__ psi, zc *__restrict__ pmn) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t n = (size_t)nb * nb * kk;
    if (i < n) {
        int col = (int)((i / nb) % nb);
        double sf = sqrt(s2[col]), si = 1.0 / sf;
        zc newpsi = make_cuDoubleComplex(pmn[i].x * si, pmn[i].y * si);
        pmn[i] = make_cuDoubleComplex(-psi[i].x * sf, -psi[i].y * sf);
        psi[i] = newpsi;
    }
}

/* cuBLAS host-pointer-mode scalar constants                                */
static const zc Z_ONE = {1.0, 0.0};
static const zc Z_ZERO = {0.0, 0.0};
static const zc Z_MONE = {-1.0, 0.0};

/* --------------------------------------------------------------------------
 * Structured path kernels (cuFFT stencil + batched-GEMM correction).
 * Grids are channel-major: channel ch = l + nb*col (or l + nb*m for H(q)),
 * each channel a contiguous padded volume of Npc cells.
 * ------------------------------------------------------------------------ */
__global__ void k_scatter(size_t nfield, int nbsq, size_t Npc,
                          const int *__restrict__ gcell,
                          const zc *__restrict__ x, zc *__restrict__ g) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nfield) {
        int ch = (int)(i % nbsq), k = (int)(i / nbsq);
        g[(size_t)ch * Npc + gcell[k]] = x[i];
    }
}

__global__ void k_gather(size_t nfield, int nbsq, size_t Npc,
                         const int *__restrict__ gcell,
                         const zc *__restrict__ g, zc *__restrict__ y) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nfield) {
        int ch = (int)(i % nbsq), k = (int)(i / nbsq);
        y[i] = g[(size_t)ch * Npc + gcell[k]];
    }
}

/* Place W(R_s)/Npc at cell (-R_s mod Np) per channel: forward FFT of this
 * gives H(q) = sum_R W(R) e^{+i q R} pre-scaled for the unnormalized IFFT. */
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
    g[(size_t)ch * Npc + cell] = make_cuDoubleComplex(w.x * scale, w.y * scale);
}

/* In-place per-cell block multiply: g(:,:,cell) <- H(q=cell) g(:,:,cell).
 * One thread block per cell, (nb x nb) threads, shared staging.            */
__global__ void k_pointwise_blockmul(size_t Npc, int nb,
                                     const zc *__restrict__ Hq,
                                     zc *__restrict__ g) {
    extern __shared__ zc sh[];
    zc *sH = sh, *sx = sh + nb * nb;
    for (size_t cell = blockIdx.x; cell < Npc; cell += gridDim.x) {
        const int l = threadIdx.x, c2 = threadIdx.y;
        const int ch = l + nb * c2;
        sH[ch] = Hq[(size_t)ch * Npc + cell];
        sx[ch] = g[(size_t)ch * Npc + cell];
        __syncthreads();
        zc acc = make_cuDoubleComplex(0.0, 0.0);
        for (int m = 0; m < nb; ++m)
            acc = zfma(sH[l + nb * m], sx[m + nb * c2], acc);
        __syncthreads();
        g[(size_t)ch * Npc + cell] = acc;
    }
}

/* --------------------------------------------------------------------------
 * rsrec_set_grid: host-side builder (identical logic to rsrec_cpu.cpp,
 * which validates it), then device uploads, H(q) tables and cuFFT plan.
 * coords: (3, kk) integer lattice coordinates. use_structured = 0 reverts
 * to the block-ELL backend.
 *
 * Device memory: 2 grid buffers + one H(q) table per operator, each
 * Npc * nb^2 * 16 B (e.g. 36^3 cells, nb = 18: ~242 MB each).
 * ------------------------------------------------------------------------ */
struct CorrEntryG { int atom, nbr, op_blk; };  /* op_blk: encoded A index   */

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
    std::vector<char> flagH(kk), flagV(kk);
    size_t nflagH = 0, nflagV = 0;
    for (int k = 0; k < kk; ++k) {
        bool irr = !shells_match(k) || (iz[k] - 1 != c->t_ref);
        flagV[k] = irr;
        flagH[k] = irr || (k < c->nmax);
        nflagH += flagH[k]; nflagV += flagV[k];
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
            "%d shells, corrected rows: H %zu/%d, V %zu/%d, "
            "device tables %.0f MB\n",
            Nx, Ny, Nz, c->Np[0], c->Np[1], c->Np[2], c->t_ref + 1,
            c->nshell_ref, nflagH, kk, c->have_v ? nflagV : 0, kk,
            (double)gridbytes * (2 + (c->have_v ? 3 : 1)) / 1e6);
    c->use_struct = true;
    return 0;
}

/* Structured y = (Op x - b x)/a: scatter -> FFT -> H(q) -> IFFT -> gather,
 * then batched-GEMM corrections, then the fused shift.                      */
static int dev_op_apply_struct(rsrec_ctx *c, int which, const zc *x, zc *y,
                               double a, double b) {
    const int nb = c->nb;
    const size_t nf = fieldsz(c), bb = (size_t)nb * nb;
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    const zc *Hq = (which == 1) ? c->d_VAq : (which == 2) ? c->d_VBq : c->d_Hq;

    k_scatter<<<nbl, tpb>>>(nf, (int)bb, c->Npc, c->d_gcell, x, c->d_gin);
    CUCHK(cudaGetLastError());
    if (cufftExecZ2Z(c->fft_plan, (cufftDoubleComplex *)c->d_gin,
                     (cufftDoubleComplex *)c->d_gw, CUFFT_FORWARD)
        != CUFFT_SUCCESS) FAIL("op_apply: forward FFT failed");
    {
        dim3 thr(nb, nb);
        int nblocks = (int)std::min<size_t>(c->Npc, 65535);
        k_pointwise_blockmul<<<nblocks, thr, 2 * bb * sizeof(zc)>>>(
            c->Npc, nb, Hq, c->d_gw);
        CUCHK(cudaGetLastError());
    }
    if (cufftExecZ2Z(c->fft_plan, (cufftDoubleComplex *)c->d_gw,
                     (cufftDoubleComplex *)c->d_gw, CUFFT_INVERSE)
        != CUFFT_SUCCESS) FAIL("op_apply: inverse FFT failed");
    k_gather<<<nbl, tpb>>>(nf, (int)bb, c->Npc, c->d_gcell, c->d_gw, y);
    CUCHK(cudaGetLastError());

    /* corrections */
    if (c->c_nnz > 0) {
        const zc **Ap = (which == 1) ? c->d_cA_VA :
                        (which == 2) ? c->d_cA_VB : c->d_cA_H;
        k_make_ptrs<<<(int)((c->c_nnz + tpb - 1) / tpb), tpb>>>(
            c->c_nnz, bb, c->d_c_atom, c->d_c_nbr, x, y, c->d_cB, c->d_cC);
        CUCHK(cudaGetLastError());
        for (size_t g = 0; g + 1 < c->cg_off.size(); ++g) {
            int off = c->cg_off[g], cnt2 = c->cg_off[g + 1] - off;
            if (cnt2 <= 0) continue;
            CBCHK(cublasZgemmBatched(c->blas, CUBLAS_OP_N, CUBLAS_OP_N,
                                     nb, nb, nb, &Z_ONE,
                                     Ap + off, nb, c->d_cB + off, nb,
                                     &Z_ONE, c->d_cC + off, nb, cnt2));
        }
    }
    if (b != 0.0 || a != 1.0) {
        k_axpby<<<nbl, tpb>>>(nf, -b / a, x, 1.0 / a, y);
        CUCHK(cudaGetLastError());
    }
    return 0;
}



/* y = (Op x - b x)/a, Op in {H, VA, VB}; x, y device, nrhs == nb.          */
static int dev_op_apply(rsrec_ctx *c, int which, const zc *x, zc *y,
                        double a, double b) {
    const int nb = c->nb;
    const size_t nf = fieldsz(c), bb = (size_t)nb * nb;

    /* Structured (cuFFT stencil + correction) path, when enabled and the
     * operator's H(q) table exists (velocity tables require set_velocity
     * to have been called BEFORE rsrec_set_grid).                           */
    if (c->use_struct) {
        const zc *tab = (which == 1) ? c->d_VAq :
                        (which == 2) ? c->d_VBq : c->d_Hq;
        if (tab) return dev_op_apply_struct(c, which, x, y, a, b);
    }

    if (c->fused_matvec) {
        dim3 thr(nb, nb);
        size_t shmem = 2 * bb * sizeof(zc);
        k_op_apply<<<c->kk, thr, shmem>>>(c->kk, nb, c->nnmax, c->nmax, nb,
                                          which, c->d_ee, c->d_hall, c->d_va,
                                          c->d_vb, c->d_nn, c->d_iz, x, y,
                                          a, b);
        CUCHK(cudaGetLastError());
        return 0;
    }

    const zc **Ap = (which == 1) ? c->d_Aptr_VA :
                    (which == 2) ? c->d_Aptr_VB : c->d_Aptr_H;
    const int tpb = 256;
    k_make_ptrs<<<(int)((c->nnz + tpb - 1) / tpb), tpb>>>(
        c->nnz, bb, c->d_pair_atom, c->d_pair_nbr, x, y, c->d_Bptr, c->d_Cptr);
    CUCHK(cudaGetLastError());
    CUCHK(cudaMemsetAsync(y, 0, nf * sizeof(zc)));

    for (int s = 0; s < c->n_shells; ++s) {
        int off = c->sh_off[s], cnt = c->sh_off[s + 1] - off;
        if (cnt == 0) continue;
        CBCHK(cublasZgemmBatched(c->blas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nb, nb, nb, &Z_ONE,
                                 Ap + off, nb, c->d_Bptr + off, nb,
                                 &Z_ONE, c->d_Cptr + off, nb, cnt));
    }
    if (b != 0.0 || a != 1.0) {
        const int nbl = (int)((nf + tpb - 1) / tpb);
        k_axpby<<<nbl, tpb>>>(nf, -b / a, x, 1.0 / a, y);
        CUCHK(cudaGetLastError());
    }
    return 0;
}

/* C (+)= sum_k A_k^H B_k: strided-batched GEMM into d_bd, then one Zgemv
 * against the ones vector reduces over atoms. C is an nb x nb device blk.  */
static int dev_block_dot(rsrec_ctx *c, const zc *A, const zc *B, zc *C,
                         bool zero_first) {
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb;
    CBCHK(cublasZgemmStridedBatched(c->blas, CUBLAS_OP_C, CUBLAS_OP_N,
                                    nb, nb, nb, &Z_ONE,
                                    A, nb, (long long)bb,
                                    B, nb, (long long)bb, &Z_ZERO,
                                    c->d_bd, nb, (long long)bb, c->kk));
    CBCHK(cublasZgemv(c->blas, CUBLAS_OP_N, (int)bb, c->kk, &Z_ONE,
                      c->d_bd, (int)bb, c->d_ones, 1,
                      zero_first ? &Z_ZERO : &Z_ONE, C, 1));
    return 0;
}

/* psi(:,:,k) <- psi(:,:,k) * M for all k (strideB = 0)                     */
static int dev_right_mul(rsrec_ctx *c, zc *psi, const zc *M, zc *tmp) {
    const int nb = c->nb;
    CBCHK(cublasZgemmStridedBatched(c->blas, CUBLAS_OP_N, CUBLAS_OP_N,
                                    nb, nb, nb, &Z_ONE,
                                    psi, nb, (long long)nb * nb, M, nb, 0,
                                    &Z_ZERO, tmp, nb, (long long)nb * nb,
                                    c->kk));
    CUCHK(cudaMemcpy(psi, tmp, fieldsz(c) * sizeof(zc),
                     cudaMemcpyDeviceToDevice));
    return 0;
}

/* B = U sqrt(ev) U^H, Bi = U ev^{-1/2} U^H on the HOST (LAPACK zheev),
 * identical to crecal_b. S, B, Bi are device pointers.                     */
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
                cplx uu = U[i + (size_t)nb * l] * std::conj(U[j + (size_t)nb * l]);
                sb += uu * lam;
                if (lam > 0.0) sbi += uu / lam;
            }
            Bm[i + (size_t)nb * j] = sb;
            Bi[i + (size_t)nb * j] = sbi;
        }
    CUCHK(cudaMemcpy(dB, Bm.data(), bb * sizeof(zc), cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(dBi, Bi.data(), bb * sizeof(zc), cudaMemcpyHostToDevice));
    return 0;
}

/* --------------------------------------------------------------------------
 * API: matvec (host arrays in/out; drivers below stay device-resident)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_op_apply(rsrec_ctx *c, int which, const void *x_,
                              void *y_, int nrhs, double a, double b) {
    if (!c->have_h) FAIL("op_apply: Hamiltonian not set");
    if (which != 0 && !c->have_v) FAIL("op_apply: velocity operators not set");
    if (nrhs != c->nb) FAIL("op_apply: this build expects nrhs == nb");
    size_t n = fieldsz(c) * sizeof(zc);
    CUCHK(cudaMemcpy(c->d_s0, x_, n, cudaMemcpyHostToDevice));
    if (dev_op_apply(c, which, c->d_s0, c->d_s1, a, b)) return 1;
    CUCHK(cudaMemcpy(y_, c->d_s1, n, cudaMemcpyDeviceToHost));
    return 0;
}

/* --------------------------------------------------------------------------
 * API: Chebyshev block moments (device-resident loop)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_chebyshev_moments(rsrec_ctx *c, const void *psi0_,
                                       int lld, double a, double b,
                                       void *mu_) {
    if (!c->have_h) FAIL("chebyshev_moments: Hamiltonian not set");
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);
    const size_t nmom = 2 * (size_t)lld + 2;

    zc *p0 = c->d_s0, *p1 = c->d_s1, *p2 = c->d_s2, *pref = c->d_s3;
    zc *dmu;
    CUCHK(cudaMalloc(&dmu, nmom * bb * sizeof(zc)));
    CUCHK(cudaMemset(dmu, 0, nmom * bb * sizeof(zc)));
    CUCHK(cudaMemcpy(pref, psi0_, nf * sizeof(zc), cudaMemcpyHostToDevice));
    CUCHK(cudaMemcpy(p0, pref, nf * sizeof(zc), cudaMemcpyDeviceToDevice));

    dev_block_dot(c, pref, p0, dmu + 0 * bb, true);       /* mu_1 */
    if (dev_op_apply(c, 0, p0, p1, a, b)) return 1;
    dev_block_dot(c, pref, p1, dmu + 1 * bb, true);       /* mu_2 */

    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    for (int ll = 1; ll <= lld; ++ll) {
        if (dev_op_apply(c, 0, p1, p2, a, b)) return 1;
        k_axpby<<<nbl, tpb>>>(nf, -1.0, p0, 2.0, p2);     /* p2 = 2p2 - p0 */
        CUCHK(cudaGetLastError());
        zc *m1 = dmu + (size_t)(2 * ll) * bb;
        zc *m2 = dmu + (size_t)(2 * ll + 1) * bb;
        dev_block_dot(c, p1, p1, m1, true);               /* dum1          */
        dev_block_dot(c, p2, p1, m2, true);               /* dum2          */
        k_axpby<<<1, (int)bb>>>(bb, -1.0, dmu + 0 * bb, 2.0, m1);
        k_axpby<<<1, (int)bb>>>(bb, -1.0, dmu + 1 * bb, 2.0, m2);
        CUCHK(cudaGetLastError());
        zc *t = p0; p0 = p1; p1 = p2; p2 = t;
    }
    CUCHK(cudaMemcpy(mu_, dmu, nmom * bb * sizeof(zc),
                     cudaMemcpyDeviceToHost));
    cudaFree(dmu);
    return 0;
}

/* --------------------------------------------------------------------------
 * API: block Lanczos (device loop, host zheev per step like crecal_b)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_block_lanczos(rsrec_ctx *c, const void *psi0_, int lld,
                                   void *a_b_, void *b2_b_) {
    if (!c->have_h) FAIL("block_lanczos: Hamiltonian not set");
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);
    cplx *a_b = (cplx *)a_b_, *b2_b = (cplx *)b2_b_;

    zc *psi = c->d_s0, *pmn = c->d_s1, *hpsi = c->d_s2, *tmp = c->d_s3;
    zc *dAn = c->d_blk, *dSum = c->d_blk + bb;
    zc *dB = c->d_blk + 2 * bb, *dBi = c->d_blk + 3 * bb;

    CUCHK(cudaMemcpy(psi, psi0_, nf * sizeof(zc), cudaMemcpyHostToDevice));
    CUCHK(cudaMemset(pmn, 0, nf * sizeof(zc)));

    std::vector<cplx> sum_b(bb, cplx(0, 0));
    for (int i = 0; i < nb; ++i) sum_b[i + (size_t)nb * i] = 1.0;
    std::memset(a_b, 0, bb * lld * sizeof(cplx));
    std::memset(b2_b, 0, bb * lld * sizeof(cplx));

    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    for (int ll = 0; ll < lld - 1; ++ll) {
        if (dev_op_apply(c, 0, psi, hpsi, 1.0, 0.0)) return 1;
        dev_block_dot(c, psi, hpsi, dAn, true);           /* A_n            */
        CUCHK(cudaMemcpy(a_b + (size_t)ll * bb, dAn, bb * sizeof(zc),
                         cudaMemcpyDeviceToHost));
        k_axpby<<<nbl, tpb>>>(nf, 1.0, hpsi, -1.0, pmn);  /* pmn = hpsi-pmn */
        CUCHK(cudaGetLastError());
        std::memcpy(b2_b + (size_t)ll * bb, sum_b.data(), bb * sizeof(cplx));

        /* pmn -= psi * A_n (strided batched, strideB = 0)                  */
        CBCHK(cublasZgemmStridedBatched(c->blas, CUBLAS_OP_N, CUBLAS_OP_N,
                                        nb, nb, nb, &Z_MONE,
                                        psi, nb, (long long)bb, dAn, nb, 0,
                                        &Z_ONE, pmn, nb, (long long)bb,
                                        c->kk));
        dev_block_dot(c, pmn, pmn, dSum, true);           /* B^2            */
        CUCHK(cudaMemcpy(sum_b.data(), dSum, bb * sizeof(zc),
                         cudaMemcpyDeviceToHost));
        if (host_eig_sqrt(c, dSum, dB, dBi)) return 1;

        CUCHK(cudaMemcpy(tmp, psi, nf * sizeof(zc),
                         cudaMemcpyDeviceToDevice));      /* psi_t          */
        if (dev_right_mul(c, pmn, dBi, hpsi)) return 1;   /* pmn = pmn*Bi   */
        zc *t = psi; psi = pmn; pmn = t;                  /* psi = new      */
        CUCHK(cudaMemcpy(pmn, tmp, nf * sizeof(zc),
                         cudaMemcpyDeviceToDevice));
        if (dev_right_mul(c, pmn, dB, hpsi)) return 1;    /* pmn = psi_t*B  */
    }
    std::memcpy(b2_b + (size_t)(lld - 1) * bb, sum_b.data(),
                bb * sizeof(cplx));
    return 0;
}

/* --------------------------------------------------------------------------
 * API: scalar Lanczos, all nb orbital chains batched as columns
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_scalar_lanczos(rsrec_ctx *c, int site_j, int lld,
                                    double *a_out, double *b2_out) {
    if (!c->have_h) FAIL("scalar_lanczos: Hamiltonian not set");
    if (site_j < 1 || site_j > c->kk) FAIL("scalar_lanczos: bad site");
    const int nb = c->nb;
    const size_t nf = fieldsz(c);

    zc *psi = c->d_s0, *pmn = c->d_s1, *hpsi = c->d_s2;
    CUCHK(cudaMemset(psi, 0, nf * sizeof(zc)));
    CUCHK(cudaMemset(pmn, 0, nf * sizeof(zc)));
    {
        std::vector<cplx> blk((size_t)nb * nb, cplx(0, 0));
        for (int l = 0; l < nb; ++l) blk[l + (size_t)nb * l] = 1.0;
        CUCHK(cudaMemcpy(psi + (size_t)nb * nb * (site_j - 1), blk.data(),
                         (size_t)nb * nb * sizeof(zc),
                         cudaMemcpyHostToDevice));
    }
    double *dA = c->d_col, *dS = c->d_col + nb;
    std::vector<double> acol(nb), s2(nb), summ(nb, 1.0);
    std::memset(a_out, 0, sizeof(double) * lld * nb);
    std::memset(b2_out, 0, sizeof(double) * lld * nb);

    const int apb = 64;
    const int nblk = (c->kk + apb - 1) / apb;
    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);
    dim3 thr(nb, nb);

    for (int ll = 0; ll < lld - 1; ++ll) {
        if (dev_op_apply(c, 0, psi, hpsi, 1.0, 0.0)) return 1;
        CUCHK(cudaMemset(dA, 0, nb * sizeof(double)));
        k_col_dot<<<nblk, thr>>>(c->kk, nb, psi, hpsi, dA, apb);
        k_axpby<<<nbl, tpb>>>(nf, 1.0, hpsi, 1.0, pmn);   /* pmn += hpsi    */
        CUCHK(cudaMemcpy(acol.data(), dA, nb * sizeof(double),
                         cudaMemcpyDeviceToHost));
        for (int col = 0; col < nb; ++col) {
            a_out[ll + (size_t)lld * col] = acol[col];
            b2_out[ll + (size_t)lld * col] = summ[col];
            acol[col] = -acol[col];
        }
        CUCHK(cudaMemcpy(dA, acol.data(), nb * sizeof(double),
                         cudaMemcpyHostToDevice));
        k_col_axpy<<<nbl, tpb>>>(c->kk, nb, dA, psi, pmn); /* pmn -= a psi  */
        CUCHK(cudaMemset(dS, 0, nb * sizeof(double)));
        k_col_dot<<<nblk, thr>>>(c->kk, nb, pmn, pmn, dS, apb);
        CUCHK(cudaMemcpy(s2.data(), dS, nb * sizeof(double),
                         cudaMemcpyDeviceToHost));
        k_col_shift<<<nbl, tpb>>>(c->kk, nb, dS, psi, pmn);
        CUCHK(cudaGetLastError());
        summ = s2;
    }
    for (int col = 0; col < nb; ++col)
        b2_out[(lld - 1) + (size_t)lld * col] = summ[col];
    return 0;
}

/* --------------------------------------------------------------------------
 * API: stochastic conductivity moments (left states device-resident)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_stochastic_moments(rsrec_ctx *c, const void *psiref_,
                                        int lld, double a, double b,
                                        void *mu_) {
    if (!c->have_h) FAIL("stochastic_moments: Hamiltonian not set");
    if (!c->have_v) FAIL("stochastic_moments: velocity operators not set");
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb, nf = fieldsz(c);

    size_t freeb, totalb;
    cudaMemGetInfo(&freeb, &totalb);
    size_t need = ((size_t)lld * nf + nf + bb * (size_t)lld * lld) * sizeof(zc);
    if (need > freeb * 9 / 10)
        FAIL("stochastic_moments: left states do not fit on the device; "
             "reduce lld or split over reference vectors / devices");

    zc *left;
    CUCHK(cudaMalloc(&left, (size_t)lld * nf * sizeof(zc)));
    zc *dmu;
    CUCHK(cudaMalloc(&dmu, bb * (size_t)lld * lld * sizeof(zc)));
    CUCHK(cudaMemset(dmu, 0, bb * (size_t)lld * lld * sizeof(zc)));
    zc *R;
    CUCHK(cudaMalloc(&R, nf * sizeof(zc)));

    zc *w0 = c->d_s0, *w1 = c->d_s1, *w2 = c->d_s2, *pref = c->d_s3;
    CUCHK(cudaMemcpy(pref, psiref_, nf * sizeof(zc), cudaMemcpyHostToDevice));

    const int tpb = 256;
    const int nbl = (int)((nf + tpb - 1) / tpb);

    /* Left states: L_m = T_{m-1}(H~)|psiref>                                */
    CUCHK(cudaMemcpy(w1, pref, nf * sizeof(zc), cudaMemcpyDeviceToDevice));
    CUCHK(cudaMemcpy(left, w1, nf * sizeof(zc), cudaMemcpyDeviceToDevice));
    for (int m = 2; m <= lld; ++m) {
        if (m == 2) {
            zc *t = w0; w0 = w1; w1 = t;
            if (dev_op_apply(c, 0, w0, w1, a, b)) return 1;
        } else {
            if (dev_op_apply(c, 0, w1, w2, a, b)) return 1;
            k_axpby<<<nbl, tpb>>>(nf, -1.0, w0, 2.0, w2);
            CUCHK(cudaGetLastError());
            zc *t = w0; w0 = w1; w1 = w2; w2 = t;
        }
        CUCHK(cudaMemcpy(left + (size_t)(m - 1) * nf, w1, nf * sizeof(zc),
                         cudaMemcpyDeviceToDevice));
    }

    /* Right recursion: v_n = T_{n-1}(H~) V_b|psiref>, R = V_a v_n           */
    if (dev_op_apply(c, 2, pref, w0, 1.0, 0.0)) return 1;
    zc *v0 = w0, *v1 = w1, *v2 = w2;
    for (int n = 1; n <= lld; ++n) {
        if (n == 1) {
            CUCHK(cudaMemcpy(v1, v0, nf * sizeof(zc),
                             cudaMemcpyDeviceToDevice));
        } else if (n == 2) {
            CUCHK(cudaMemcpy(v0, v1, nf * sizeof(zc),
                             cudaMemcpyDeviceToDevice));
            if (dev_op_apply(c, 0, v0, v1, a, b)) return 1;
        } else {
            if (dev_op_apply(c, 0, v1, v2, a, b)) return 1;
            k_axpby<<<nbl, tpb>>>(nf, -1.0, v0, 2.0, v2);
            CUCHK(cudaGetLastError());
            zc *t = v0; v0 = v1; v1 = v2; v2 = t;
        }
        if (dev_op_apply(c, 1, v1, R, 1.0, 0.0)) return 1;
        for (int m = 0; m < lld; ++m)
            dev_block_dot(c, left + (size_t)m * nf, R,
                          dmu + bb * ((size_t)(n - 1) + (size_t)lld * m),
                          true);
    }
    CUCHK(cudaMemcpy(mu_, dmu, bb * (size_t)lld * lld * sizeof(zc),
                     cudaMemcpyDeviceToHost));
    cudaFree(left); cudaFree(dmu); cudaFree(R);
    return 0;
}
