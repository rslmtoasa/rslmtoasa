/* ===========================================================================
 * rsrec_cpu_fast.cpp -- optimized CPU backend, applying the GPU lessons.
 *
 * Lessons transferred from the CUDA port:
 *   1. SITE-MAJOR layout: psi packed as (ld = nb*kk) x nrhs column-major.
 *      Then every block reduction  C = sum_k A_k^H B_k  is ONE tall-skinny
 *      zgemm('C','N', nb, nb, ld), and every right-multiply psi*M is ONE
 *      zgemm('N','N', ld, nb, nb). The hand-rolled triple loops + omp-critical
 *      reductions of the reference are replaced by BLAS-3 (cache-blocked,
 *      threaded, vectorized by OpenBLAS/MKL).
 *   2. The matvec STAYS a custom block-ELL kernel (not tiny per-shell zgemm):
 *      on CPU a 18x18x18 zgemm is dominated by call/packing overhead, the
 *      opposite of the GPU where batching them won. We instead keep a fused,
 *      branch-free, site-major manual kernel with a local accumulator -- the
 *      structural analogue of the GPU fused step kernel.
 *   3. Constants folded once (lsham into shell 0; scaling done in the matvec
 *      epilogue), no separate scale/shift passes, minimal temporaries.
 *
 * fp64 throughout (correctness reference + portable backend). An fp32 KPM
 * path mirrors the GPU one for ~2x SIMD/throughput; see the note at the end.
 *
 * Build (versioned libs needed in minimal containers):
 *   g++ -O3 -march=native -fopenmp -fPIC -shared rsrec_cpu_fast.cpp \
 *       -o librsrec_cpu.so -lopenblas        # or: liblapack.so.3 libblas.so.3
 * =========================================================================== */
#include "rsrec.h"

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <type_traits>
#include <string>

typedef std::complex<double> cplx;
typedef std::complex<float> fcplx;

extern "C" {
void zheev_(const char *, const char *, const int *, cplx *, const int *,
            double *, cplx *, const int *, double *, int *);
void zgemm_(const char *, const char *, const int *, const int *, const int *,
            const cplx *, const cplx *, const int *, const cplx *, const int *,
            const cplx *, cplx *, const int *);
void cgemm_(const char *, const char *, const int *, const int *, const int *,
            const fcplx *, const fcplx *, const int *, const fcplx *,
            const int *, const fcplx *, fcplx *, const int *);
}

static std::string g_err;
extern "C" const char *rsrec_last_error(void) { return g_err.c_str(); }
#define FAIL(msg) do { g_err = (msg); return 1; } while (0)

static const cplx Z1(1, 0), Z0(0, 0), ZM1(-1, 0);

/* --------------------------------------------------------------------------
 * Context. Operators stay in the upload (Fortran) block layout; wavefunction
 * fields are SITE-MAJOR (ld = nb*kk rows, columns = right-hand sides).
 * ------------------------------------------------------------------------ */
struct rsrec_ctx {
    int kk = 0, nb = 0, nnmax = 0, ntype = 0, nmax = 0;
    std::vector<cplx> ee, hall, va, vb;     /* (nb,nb,nnmax,*) Fortran      */
    std::vector<int> nn, iz;
    bool have_h = false, have_v = false;
    std::vector<cplx> s0, s1, s2, s3;        /* site-major scratch fields    */

    /* fp32 KPM: scaled (H-b)/a copies, baked once per (a,b)/upload.         */
    int prec = 0;                            /* 0 = fp32 (default), 1 = fp64 */
    std::vector<fcplx> f_ee, f_hall;
    double f_a = 0.0, f_b = 0.0;
    bool f_ok = false;
};

static inline size_t BLK(const rsrec_ctx *c, int s, int slot) {
    return (size_t)c->nb * c->nb * ((size_t)s + (size_t)c->nnmax * slot);
}
static inline int NNF(const rsrec_ctx *c, int k, int s) {
    return c->nn[(size_t)k + (size_t)c->kk * s];
}
static inline size_t LD(const rsrec_ctx *c) { return (size_t)c->nb * c->kk; }

extern "C" rsrec_ctx *rsrec_create(int kk, int nb, int nnmax, int ntype,
                                   int nmax, int) {
    if (kk <= 0 || nb <= 0 || nnmax <= 0 || ntype <= 0 || nmax < 0 ||
        nmax > kk) { g_err = "rsrec_create: bad dimensions"; return nullptr; }
    rsrec_ctx *c = new rsrec_ctx();
    c->kk = kk; c->nb = nb; c->nnmax = nnmax; c->ntype = ntype; c->nmax = nmax;
    size_t nf = (size_t)nb * nb * kk;
    c->s0.resize(nf); c->s1.resize(nf); c->s2.resize(nf); c->s3.resize(nf);
    return c;
}
extern "C" void rsrec_destroy(rsrec_ctx *c) { delete c; }

extern "C" int rsrec_set_hamiltonian(rsrec_ctx *c, const void *ee_,
                                     const void *hall_, const void *lsham_,
                                     const int *nn, const int *iz,
                                     const void *eeo_, const void *hallo_,
                                     const void *enim_) {
    (void)eeo_; (void)hallo_; (void)enim_;
    const cplx *ee = (const cplx *)ee_, *hall = (const cplx *)hall_;
    const cplx *ls = (const cplx *)lsham_;
    if (!ee || !nn || !iz) FAIL("set_hamiltonian: null input");
    if (c->nmax > 0 && !hall) FAIL("set_hamiltonian: nmax>0 but hall is null");
    const int nb = c->nb;
    size_t nee = (size_t)nb * nb * c->nnmax * c->ntype;
    c->ee.assign(ee, ee + nee);
    if (c->nmax > 0)
        c->hall.assign(hall, hall + (size_t)nb * nb * c->nnmax * c->nmax);
    c->nn.assign(nn, nn + (size_t)c->kk * c->nnmax);
    c->iz.assign(iz, iz + c->kk);
    if (ls) {
        for (int t = 0; t < c->ntype; ++t)
            for (size_t i = 0; i < (size_t)nb * nb; ++i)
                c->ee[BLK(c, 0, t) + i] += ls[(size_t)nb * nb * t + i];
        for (int k = 0; k < c->nmax; ++k) {
            int t = c->iz[k] - 1;
            for (size_t i = 0; i < (size_t)nb * nb; ++i)
                c->hall[BLK(c, 0, k) + i] += ls[(size_t)nb * nb * t + i];
        }
    }
    c->have_h = true;
    c->f_ok = false;
    return 0;
}

extern "C" int rsrec_set_hamiltonian_additive(rsrec_ctx *c,
                                               const void *ee_add_,
                                               const void *hall_add_) {
    (void)c; (void)ee_add_; (void)hall_add_;
    return 0;
}

extern "C" int rsrec_set_velocity(rsrec_ctx *c, const void *va_,
                                  const void *vb_, const void *voa_,
                                  const void *vob_) {
    (void)voa_; (void)vob_;
    const cplx *va = (const cplx *)va_, *vb = (const cplx *)vb_;
    if (!va || !vb) FAIL("set_velocity: null input");
    size_t n = (size_t)c->nb * c->nb * c->nnmax * c->ntype;
    c->va.assign(va, va + n); c->vb.assign(vb, vb + n);
    c->have_v = true;
    return 0;
}

extern "C" int rsrec_set_precision(rsrec_ctx *c, int prec) {
    if (prec < 0 || prec > 1) FAIL("set_precision: 0 = fp32, 1 = fp64");
    c->prec = prec;
    return 0;
}

/* Bake the scaled operator Ht = (H - b I)/a into fp32 copies (cached by
 * (a,b)). The on-site -b applies only to shell-0 diagonal, like the GPU. */
static void ensure_f32(rsrec_ctx *c, double a, double b) {
    if (c->f_ok && c->f_a == a && c->f_b == b) return;
    const int nb = c->nb;
    const double inva = 1.0 / a;
    c->f_ee.resize(c->ee.size());
    for (int t = 0; t < c->ntype; ++t)
        for (int s = 0; s < c->nnmax; ++s)
            for (int m = 0; m < nb; ++m)
                for (int l = 0; l < nb; ++l) {
                    size_t idx = BLK(c, s, t) + (size_t)l + (size_t)nb * m;
                    cplx v = c->ee[idx];
                    if (s == 0 && l == m) v -= b;
                    c->f_ee[idx] = fcplx((float)(v.real() * inva),
                                         (float)(v.imag() * inva));
                }
    if (c->nmax > 0) {
        c->f_hall.resize(c->hall.size());
        for (int k = 0; k < c->nmax; ++k)
            for (int s = 0; s < c->nnmax; ++s)
                for (int m = 0; m < nb; ++m)
                    for (int l = 0; l < nb; ++l) {
                        size_t idx = BLK(c, s, k) + (size_t)l + (size_t)nb * m;
                        cplx v = c->hall[idx];
                        if (s == 0 && l == m) v -= b;
                        c->f_hall[idx] = fcplx((float)(v.real() * inva),
                                               (float)(v.imag() * inva));
                    }
    }
    c->f_a = a; c->f_b = b; c->f_ok = true;
}

/* --------------------------------------------------------------------------
 * Layout helpers: atom-major (nb, ncols, kk) fp64 -> site-major (ld, ncols)
 * in either precision (CT = cplx or fcplx).
 * ------------------------------------------------------------------------ */
template <class CT>
static void pack(const rsrec_ctx *c, const cplx *am, CT *sm, int ncols) {
    const int nb = c->nb; const size_t ld = LD(c);
#pragma omp parallel for schedule(static)
    for (int k = 0; k < c->kk; ++k)
        for (int col = 0; col < ncols; ++col)
            for (int l = 0; l < nb; ++l) {
                const cplx v = am[(size_t)l + nb * (col + (size_t)ncols * k)];
                sm[(size_t)nb * k + l + ld * col] =
                    CT((typename CT::value_type)v.real(),
                       (typename CT::value_type)v.imag());
            }
}
static void unpack(const rsrec_ctx *c, const cplx *sm, cplx *am, int ncols) {
    const int nb = c->nb; const size_t ld = LD(c);
#pragma omp parallel for schedule(static)
    for (int k = 0; k < c->kk; ++k)
        for (int col = 0; col < ncols; ++col)
            for (int l = 0; l < nb; ++l)
                am[(size_t)l + nb * (col + (size_t)ncols * k)] =
                    sm[(size_t)nb * k + l + ld * col];
}

/* C = A^H B over site-major fields: ONE tall-skinny gemm (z or c).         */
static void gram(const rsrec_ctx *c, const cplx *A, int na, const cplx *B,
                 int nbc, cplx *C) {
    const int ld = (int)LD(c);
    zgemm_("C", "N", &na, &nbc, &ld, &Z1, A, &ld, B, &ld, &Z0, C, &na);
}
static void gram(const rsrec_ctx *c, const fcplx *A, int na, const fcplx *B,
                 int nbc, fcplx *C) {
    const int ld = (int)LD(c);
    static const fcplx F1(1, 0), F0(0, 0);
    cgemm_("C", "N", &na, &nbc, &ld, &F1, A, &ld, B, &ld, &F0, C, &na);
}

/* psi <- psi * M  (ld x nb x nb), into out (caller swaps).                 */
static void rmul(const rsrec_ctx *c, const cplx *psi, const cplx *M,
                 cplx *out) {
    const int ld = (int)LD(c), nb = c->nb;
    zgemm_("N", "N", &ld, &nb, &nb, &Z1, psi, &ld, M, &nb, &Z0, out, &ld);
}

/* --------------------------------------------------------------------------
 * Fused site-major block-Chebyshev operator apply:  y = (Op x - bsc x)*inva.
 *
 * This is a BLOCK recursion: x carries nb right-hand sides, so the per-shell
 * operation is a small GEMM B(nb x nb) * x_nbr(nb x nb), NOT a matvec. The
 * loop order below is register-blocked GEMM order: the k-index (m) is
 * outermost, so each operator column B(:,m) is loaded ONCE and reused as a
 * rank-1 update across all nb psi columns -> arithmetic intensity ~nb, the
 * compute-bound regime (the reason the GPU liked this kernel). The output
 * tile acc(nb x nb) stays resident; the nb-length complex FMA vectorizes.
 *
 * Operator blocks are passed explicitly (btype per type, bimp per atom) so
 * the identical code serves fp64 (ee/hall) and fp32 (baked (H-b)/a copies).
 * For the fp32 path the scaling is already baked, so inva=1, bsc=0.
 * ------------------------------------------------------------------------ */
template <class CT, int NB>
static void opcore_t(const rsrec_ctx *c, const CT *btype, const CT *bimp,
                     int use_imp, const CT *x, CT *y,
                     typename CT::value_type inva,
                     typename CT::value_type bsc) {
    const size_t ld = LD(c);
    const bool scale = (bsc != (typename CT::value_type)0 ||
                        inva != (typename CT::value_type)1);
#pragma omp parallel
    {
        CT acc[NB * NB];
#pragma omp for schedule(dynamic, 64)
        for (int k = 0; k < c->kk; ++k) {
            const bool imp = use_imp && (k < c->nmax);
            const CT *base = imp
                ? bimp + (size_t)NB * NB * ((size_t)c->nnmax * k)
                : btype + (size_t)NB * NB * ((size_t)c->nnmax * (c->iz[k] - 1));
            for (int i = 0; i < NB * NB; ++i) acc[i] = CT(0);
            const int nr = NNF(c, k, 0);
            for (int s = 0; s < nr; ++s) {
                int nbr = (s == 0) ? k : NNF(c, k, s) - 1;
                if (nbr < 0) continue;
                const CT *B = base + (size_t)NB * NB * s;
                const CT *xb = x + (size_t)NB * nbr;
                for (int m = 0; m < NB; ++m) {       /* GEMM k-index outer   */
                    const CT *Bm = B + (size_t)NB * m;
                    for (int col = 0; col < NB; ++col) {
                        const CT xm = xb[m + ld * col];
                        CT *ac = acc + NB * col;
#pragma omp simd
                        for (int l = 0; l < NB; ++l) ac[l] += Bm[l] * xm;
                    }
                }
            }
            CT *yk = y + (size_t)NB * k;
            const CT *xk = x + (size_t)NB * k;
            for (int col = 0; col < NB; ++col)
                for (int l = 0; l < NB; ++l) {
                    CT v = acc[NB * col + l];
                    if (scale) v = (v - bsc * xk[l + ld * col]) * inva;
                    yk[l + ld * col] = v;
                }
        }
    }
}

template <class CT>
static void opcore_dyn(const rsrec_ctx *c, const CT *btype, const CT *bimp,
                       int use_imp, const CT *x, CT *y, int nrhs,
                       typename CT::value_type inva,
                       typename CT::value_type bsc) {
    const int nb = c->nb; const size_t ld = LD(c);
    const bool scale = (bsc != (typename CT::value_type)0 ||
                        inva != (typename CT::value_type)1);
#pragma omp parallel
    {
        std::vector<CT> acc((size_t)nb * nrhs);
#pragma omp for schedule(dynamic, 64)
        for (int k = 0; k < c->kk; ++k) {
            const bool imp = use_imp && (k < c->nmax);
            const CT *base = imp
                ? bimp + (size_t)nb * nb * ((size_t)c->nnmax * k)
                : btype + (size_t)nb * nb * ((size_t)c->nnmax * (c->iz[k] - 1));
            std::fill(acc.begin(), acc.end(), CT(0));
            const int nr = NNF(c, k, 0);
            for (int s = 0; s < nr; ++s) {
                int nbr = (s == 0) ? k : NNF(c, k, s) - 1;
                if (nbr < 0) continue;
                const CT *B = base + (size_t)nb * nb * s;
                const CT *xb = x + (size_t)nb * nbr;
                for (int m = 0; m < nb; ++m) {
                    const CT *Bm = B + (size_t)nb * m;
                    for (int col = 0; col < nrhs; ++col) {
                        const CT xm = xb[m + ld * col];
                        CT *ac = acc.data() + (size_t)nb * col;
                        for (int l = 0; l < nb; ++l) ac[l] += Bm[l] * xm;
                    }
                }
            }
            CT *yk = y + (size_t)nb * k;
            const CT *xk = x + (size_t)nb * k;
            for (int col = 0; col < nrhs; ++col)
                for (int l = 0; l < nb; ++l) {
                    CT v = acc[(size_t)nb * col + l];
                    if (scale) v = (v - bsc * xk[l + ld * col]) * inva;
                    yk[l + ld * col] = v;
                }
        }
    }
}

/* nb-dispatch (s/sp/spd/spdf x 2 spinors), generic over precision.         */
template <class CT>
static void opapply(const rsrec_ctx *c, const CT *bt, const CT *bi,
                    int use_imp, const CT *x, CT *y, int nrhs,
                    typename CT::value_type inva,
                    typename CT::value_type bsc) {
    if (nrhs == c->nb) {
        switch (c->nb) {
            case 2:  opcore_t<CT, 2>(c, bt, bi, use_imp, x, y, inva, bsc); return;
            case 8:  opcore_t<CT, 8>(c, bt, bi, use_imp, x, y, inva, bsc); return;
            case 16: opcore_t<CT,16>(c, bt, bi, use_imp, x, y, inva, bsc); return;
            case 18: opcore_t<CT,18>(c, bt, bi, use_imp, x, y, inva, bsc); return;
            case 32: opcore_t<CT,32>(c, bt, bi, use_imp, x, y, inva, bsc); return;
            default: break;
        }
    }
    opcore_dyn<CT>(c, bt, bi, use_imp, x, y, nrhs, inva, bsc);
}

/* fp64 wrapper used by the Lanczos/stochastic drivers: y = (Op x - b x)/a. */
static void op_apply_sm(const rsrec_ctx *c, int which, const cplx *x,
                        cplx *y, int nrhs, double a, double b) {
    const cplx *bt = (which == 1) ? c->va.data()
                   : (which == 2) ? c->vb.data() : c->ee.data();
    opapply<cplx>(c, bt, c->hall.data(), which == 0 ? 1 : 0,
                  x, y, nrhs, 1.0 / a, b);
}

extern "C" int rsrec_op_apply(rsrec_ctx *c, int which, const void *x,
                              void *y, int nrhs, double a, double b) {
    if (!c->have_h) FAIL("op_apply: Hamiltonian not set");
    if (which != 0 && !c->have_v) FAIL("op_apply: velocity operators not set");
    cplx *xs = c->s0.data(), *ys = c->s1.data();
    pack<cplx>(c, (const cplx *)x, xs, nrhs);
    op_apply_sm(c, which, xs, ys, nrhs, a, b);
    unpack(c, ys, (cplx *)y, nrhs);
    return 0;
}

/* B = U sqrt(ev) U^H, Bi = U ev^{-1/2} U^H from Hermitian S (= crecal_b).   */
static int eig_sqrt(int nb, const cplx *S, cplx *B, cplx *Bi) {
    std::vector<cplx> U(S, S + (size_t)nb * nb);
    std::vector<double> ev(nb), rwork(3 * nb - 2);
    int lwork = nb * nb, info = 0;
    std::vector<cplx> work(lwork);
    zheev_("V", "U", &nb, U.data(), &nb, ev.data(), work.data(), &lwork,
           rwork.data(), &info);
    if (info != 0) return info;
    for (int j = 0; j < nb; ++j)
        for (int i = 0; i < nb; ++i) {
            cplx sb(0, 0), sbi(0, 0);
            for (int l = 0; l < nb; ++l) {
                double lam = std::sqrt(std::max(0.0, ev[l]));
                cplx uu = U[i + (size_t)nb * l] *
                          std::conj(U[j + (size_t)nb * l]);
                sb += uu * lam;
                if (lam > 0.0) sbi += uu / lam;
            }
            B[i + (size_t)nb * j] = sb;
            Bi[i + (size_t)nb * j] = sbi;
        }
    return 0;
}

/* --------------------------------------------------------------------------
 * Chebyshev block moments (double-moment trick). gram() does every reduction.
 *
 * Precision (rsrec_set_precision): 0 = fp32 (default, KPM-safe ~5e-7, ~2x
 * from doubled SIMD lanes + halved memory traffic on the bandwidth-sensitive
 * block recursion), 1 = fp64 (bit-reference). The fp32 path bakes the scaled
 * operator Ht = (H-b)/a into fp32 copies, so its apply runs unscaled; moment
 * blocks are combined in fp64 from the fp32 grams. Mirrors the GPU engine.
 * ------------------------------------------------------------------------ */
template <class CT>
static void cheb_engine(rsrec_ctx *c, const CT *bt, const CT *bi,
                        const cplx *psi0, int lld, double a, double b,
                        cplx *mu) {
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb, nf = blk * c->kk;
    const bool f32 = !std::is_same<CT, cplx>::value;
    const typename CT::value_type inva =
        f32 ? (typename CT::value_type)1 : (typename CT::value_type)(1.0 / a);
    const typename CT::value_type bsc =
        f32 ? (typename CT::value_type)0 : (typename CT::value_type)b;

    std::memset(mu, 0, sizeof(cplx) * blk * (2 * (size_t)lld + 2));
    std::vector<CT> P0(nf), P1(nf), P2(nf), D(blk);
    CT *p0 = P0.data(), *p1 = P1.data(), *p2 = P2.data();
    pack<CT>(c, psi0, p0, nb);

    auto store = [&](const CT *src, cplx *dst) {
        for (size_t i = 0; i < blk; ++i)
            dst[i] = cplx((double)src[i].real(), (double)src[i].imag());
    };
    gram(c, p0, nb, p0, nb, D.data()); store(D.data(), mu + 0 * blk);  /* mu1 */
    opapply<CT>(c, bt, bi, 1, p0, p1, nb, inva, bsc);                  /* Ht p0 */
    gram(c, p0, nb, p1, nb, D.data()); store(D.data(), mu + 1 * blk);  /* mu2 */
    for (int ll = 1; ll <= lld; ++ll) {
        opapply<CT>(c, bt, bi, 1, p1, p2, nb, inva, bsc);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)nf; ++i)
            p2[i] = (typename CT::value_type)2 * p2[i] - p0[i];
        cplx *m1 = mu + (size_t)(2 * ll) * blk;
        cplx *m2 = mu + (size_t)(2 * ll + 1) * blk;
        gram(c, p1, nb, p1, nb, D.data()); store(D.data(), m1);   /* dum1 */
        gram(c, p2, nb, p1, nb, D.data()); store(D.data(), m2);   /* dum2 */
        for (size_t i = 0; i < blk; ++i) {
            m1[i] = 2.0 * m1[i] - mu[i];
            m2[i] = 2.0 * m2[i] - mu[blk + i];
        }
        CT *t = p0; p0 = p1; p1 = p2; p2 = t;
    }
}

extern "C" int rsrec_chebyshev_moments(rsrec_ctx *c, const void *psi0_,
                                       int lld, double a, double b, void *mu_) {
    if (!c->have_h) FAIL("chebyshev_moments: Hamiltonian not set");
    const cplx *psi0 = (const cplx *)psi0_;
    cplx *mu = (cplx *)mu_;
    if (c->prec == 0) {
        ensure_f32(c, a, b);
        cheb_engine<fcplx>(c, c->f_ee.data(),
                           c->nmax > 0 ? c->f_hall.data() : nullptr,
                           psi0, lld, a, b, mu);
    } else {
        cheb_engine<cplx>(c, c->ee.data(),
                          c->nmax > 0 ? c->hall.data() : nullptr,
                          psi0, lld, a, b, mu);
    }
    return 0;
}

/* --------------------------------------------------------------------------
 * Block Lanczos (crecal_b / hop_b), site-major + BLAS.
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_block_lanczos(rsrec_ctx *c, const void *psi0_, int lld,
                                   void *a_b_, void *b2_b_, int prec) {
    (void)prec;  /* this CPU engine is ee-only fp64; precision flag unused */
    if (!c->have_h) FAIL("block_lanczos: Hamiltonian not set");
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb, nf = blk * c->kk;
    cplx *a_b = (cplx *)a_b_, *b2_b = (cplx *)b2_b_;

    std::vector<cplx> PSI(nf), PMN(nf), HPS(nf), GT(nf);
    cplx *psi = PSI.data(), *pmn = PMN.data(), *hpsi = HPS.data(),
         *gt = GT.data();
    pack<cplx>(c, (const cplx *)psi0_, psi, nb);
    std::memset(pmn, 0, sizeof(cplx) * nf);

    std::vector<cplx> An(blk), sum_b(blk, Z0), B(blk), Bi(blk);
    for (int i = 0; i < nb; ++i) sum_b[i + (size_t)nb * i] = 1.0;
    std::memset(a_b, 0, sizeof(cplx) * blk * lld);
    std::memset(b2_b, 0, sizeof(cplx) * blk * lld);
    const int ld = (int)LD(c);

    for (int ll = 0; ll < lld - 1; ++ll) {
        op_apply_sm(c, 0, psi, hpsi, nb, 1.0, 0.0);
        gram(c, psi, nb, hpsi, nb, An.data());                /* A_n        */
        std::memcpy(a_b + (size_t)ll * blk, An.data(), sizeof(cplx) * blk);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)nf; ++i) pmn[i] = hpsi[i] - pmn[i];
        std::memcpy(b2_b + (size_t)ll * blk, sum_b.data(), sizeof(cplx) * blk);
        /* pmn -= psi * A_n : one zgemm, beta = 1                            */
        zgemm_("N", "N", &ld, &nb, &nb, &ZM1, psi, &ld, An.data(), &nb,
               &Z1, pmn, &ld);
        gram(c, pmn, nb, pmn, nb, sum_b.data());              /* B^2        */
        if (eig_sqrt(nb, sum_b.data(), B.data(), Bi.data()) != 0)
            FAIL("block_lanczos: zheev failed");
        rmul(c, pmn, Bi.data(), gt);            /* gt = pmn*Bi = psi_{n+1}  */
        rmul(c, psi, B.data(), pmn);            /* pmn = psi_n * B          */
        cplx *t = psi; psi = gt; gt = t;        /* psi <- new               */
    }
    std::memcpy(b2_b + (size_t)(lld - 1) * blk, sum_b.data(),
                sizeof(cplx) * blk);
    return 0;
}

/* --------------------------------------------------------------------------
 * Scalar Lanczos (nb independent chains as columns). The per-column dots
 * a_col and |pmn_col|^2 are the DIAGONAL of a gram() -- one zgemm replaces
 * the hand reduction (mirrors the GPU recommendation; no atomics).
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_scalar_lanczos(rsrec_ctx *c, int site_j, int lld,
                                    double *a_out, double *b2_out) {
    if (!c->have_h) FAIL("scalar_lanczos: Hamiltonian not set");
    if (site_j < 1 || site_j > c->kk) FAIL("scalar_lanczos: bad site");
    const int nb = c->nb, j0 = site_j - 1;
    const size_t blk = (size_t)nb * nb, nf = blk * c->kk, ld = LD(c);

    std::vector<cplx> PSI(nf, Z0), PMN(nf, Z0), HPS(nf), G(blk);
    cplx *psi = PSI.data(), *pmn = PMN.data(), *hpsi = HPS.data();
    for (int l = 0; l < nb; ++l)                  /* column l = e_l on site j */
        psi[(size_t)nb * j0 + l + ld * l] = 1.0;

    std::vector<double> summ(nb, 1.0), acol(nb), s2(nb);
    for (int i = 0; i < lld * nb; ++i) { a_out[i] = 0.0; b2_out[i] = 0.0; }

    for (int ll = 0; ll < lld - 1; ++ll) {
        op_apply_sm(c, 0, psi, hpsi, nb, 1.0, 0.0);
        gram(c, psi, nb, hpsi, nb, G.data());          /* diag = a_col     */
        for (int col = 0; col < nb; ++col) {
            acol[col] = std::real(G[col + (size_t)nb * col]);
            a_out[ll + (size_t)lld * col] = acol[col];
            b2_out[ll + (size_t)lld * col] = summ[col];
        }
#pragma omp parallel for schedule(static)             /* pmn += hpsi - a psi */
        for (int col = 0; col < nb; ++col)
            for (size_t r = 0; r < ld; ++r) {
                size_t id = r + ld * col;
                pmn[id] += hpsi[id] - acol[col] * psi[id];
            }
        gram(c, pmn, nb, pmn, nb, G.data());            /* diag = |pmn|^2   */
        for (int col = 0; col < nb; ++col)
            s2[col] = std::real(G[col + (size_t)nb * col]);
#pragma omp parallel for schedule(static)
        for (int col = 0; col < nb; ++col) {
            double sf = std::sqrt(s2[col]), si = sf > 0 ? 1.0 / sf : 0.0;
            for (size_t r = 0; r < ld; ++r) {
                size_t id = r + ld * col;
                cplx np = pmn[id] * si;
                pmn[id] = -psi[id] * sf;
                psi[id] = np;
            }
        }
        summ = s2;
    }
    for (int col = 0; col < nb; ++col)
        b2_out[(lld - 1) + (size_t)lld * col] = summ[col];
    return 0;
}

/* --------------------------------------------------------------------------
 * Stochastic conductivity moments. Left states stored site-major and stacked
 * as ONE (ld, nb*lld) matrix, so the (n, all-m) contraction L^H R is ONE
 * zgemm (lld^2 hand reductions -> lld zgemms; mirrors the GPU recommendation).
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_stochastic_moments(rsrec_ctx *c, const void *psiref_,
                                        int lld, double a, double b,
                                        void *mu_nm_) {
    if (!c->have_h) FAIL("stochastic_moments: Hamiltonian not set");
    if (!c->have_v) FAIL("stochastic_moments: velocity operators not set");
    cplx *mu = (cplx *)mu_nm_;
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb, nf = blk * c->kk, ld = LD(c);
    std::memset(mu, 0, sizeof(cplx) * blk * (size_t)lld * lld);

    std::vector<cplx> left((size_t)lld * nf);          /* (ld, nb*lld)      */
    std::vector<cplx> W0(nf), W1(nf), W2(nf), R(nf), Cbig((size_t)nb * lld * nb);
    cplx *w0 = W0.data(), *w1 = W1.data(), *w2 = W2.data();

    pack<cplx>(c, (const cplx *)psiref_, w1, nb);
    std::memcpy(left.data(), w1, sizeof(cplx) * nf);
    for (int m = 2; m <= lld; ++m) {
        if (m == 2) { std::memcpy(w0, w1, sizeof(cplx) * nf);
                      op_apply_sm(c, 0, w0, w1, nb, a, b); }
        else {
            op_apply_sm(c, 0, w1, w2, nb, a, b);
#pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)nf; ++i) w2[i] = 2.0 * w2[i] - w0[i];
            cplx *t = w0; w0 = w1; w1 = w2; w2 = t;
        }
        std::memcpy(left.data() + (size_t)(m - 1) * nf, w1, sizeof(cplx) * nf);
    }

    /* right recursion; left slot 0 = psiref (site-major already)            */
    cplx *v0 = w0, *v1 = w1, *v2 = w2;
    op_apply_sm(c, 2, left.data(), v0, nb, 1.0, 0.0);  /* v0 = V_b psiref   */
    const int ldi = (int)ld, M = nb * lld;
    for (int n = 1; n <= lld; ++n) {
        if (n == 1) std::memcpy(v1, v0, sizeof(cplx) * nf);
        else if (n == 2) { std::memcpy(v0, v1, sizeof(cplx) * nf);
                           op_apply_sm(c, 0, v0, v1, nb, a, b); }
        else {
            op_apply_sm(c, 0, v1, v2, nb, a, b);
#pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)nf; ++i) v2[i] = 2.0 * v2[i] - v0[i];
            cplx *t = v0; v0 = v1; v1 = v2; v2 = t;
        }
        op_apply_sm(c, 1, v1, R.data(), nb, 1.0, 0.0); /* R = V_a v_n       */
        /* Cbig (nb*lld x nb) = left^H R : block m at rows [nb*m, nb*m+nb)   */
        zgemm_("C", "N", &M, &nb, &ldi, &Z1, left.data(), &ldi, R.data(),
               &ldi, &Z0, Cbig.data(), &M);
        for (int m = 0; m < lld; ++m) {
            cplx *mum = mu + blk * ((size_t)(n - 1) + (size_t)lld * m);
            for (int j = 0; j < nb; ++j)
                for (int i = 0; i < nb; ++i)
                    mum[i + (size_t)nb * j] = Cbig[(size_t)nb * m + i +
                                                   (size_t)M * j];
        }
    }
    return 0;
}

/* Structured (cuFFT-stencil) backend is GPU-only. On CPU we always run the
 * block-ELL path, so set_grid is a no-op: succeed for use_structured == 0
 * (explicit ELL request), and report that structured mode is unavailable
 * otherwise so callers don't silently assume a backend they aren't getting. */
extern "C" int rsrec_set_grid(rsrec_ctx *c, const int *coords,
                              int use_structured) {
    (void)c; (void)coords;
    if (use_structured)
        FAIL("set_grid: structured/FFT backend is GPU-only; the CPU library "
             "always uses the block-ELL path");
    return 0;
}

/* Transfer matrix F(i,ie) = g_J(i) c_i (-i) e^{-i(i-1)acos(x)} /
 *                           sqrt(a^2 - (E-b)^2),  x = (E-b)/a.
 * Both i-dependent angles (Jackson thl = pi*i/(N+1) and the phase ang = i*th)
 * advance by a constant step, so they are produced by complex-rotation
 * recurrences -- one acos per energy point instead of 2 trig calls per
 * element. Rounding over the recurrence is ~N*eps (~1e-13 at N~1e3), far
 * inside KPM tolerance. */
static void build_F(int N, int nv, const double *ene, double a, double b,
                    cplx *F) {
    const double PI = 3.14159265358979323846;
    const double sd = std::sin(PI / (N + 1.0));
    const double cot = sd > 0.0 ? std::cos(PI / (N + 1.0)) / sd : 0.0;
    const double jc = std::cos(PI / (N + 1.0)), js = sd;   /* Jackson step    */
#pragma omp parallel for schedule(static)
    for (int ie = 0; ie < nv; ++ie) {
        const double d = ene[ie] - b;
        const double th = std::acos(d / a);
        const double pref = 1.0 / std::sqrt(a * a - d * d);
        const double rc = std::cos(th), rs = -std::sin(th);  /* rot=e^{-i th} */
        /* Jackson angle (cos,sin) and phase (-i)e^{-i*i*th}, both i=0 init    */
        double jcos = 1.0, jsin = 0.0;                       /* thl_0 = 0      */
        double pr = 0.0, pi = -1.0;                          /* phase_0 = -i   */
        cplx *Fie = F + (size_t)N * ie;
        for (int i = 0; i < N; ++i) {
            const double gj = ((N - i + 1) * jcos + jsin * cot) / (N + 1.0);
            const double w = gj * (i == 0 ? 1.0 : 2.0) * pref;
            Fie[i] = cplx(w * pr, w * pi);
            /* advance Jackson angle by pi/(N+1) and phase by th             */
            double nc = jcos * jc - jsin * js;
            jsin = jsin * jc + jcos * js;  jcos = nc;
            double npr = pr * rc - pi * rs;
            pi = pr * rs + pi * rc;        pr = npr;
        }
    }
}

extern "C" int rsrec_chebyshev_dos(rsrec_ctx *c, const void *mu_, int n_mom,
                                   int natoms, const double *ene, int nv,
                                   double a, double b, void *g0_) {
    if (n_mom < 1 || natoms < 1 || nv < 1) FAIL("chebyshev_dos: bad sizes");
    const int bb = c->nb * c->nb;
    const cplx *mu = (const cplx *)mu_;
    cplx *g0 = (cplx *)g0_;
    std::vector<cplx> F((size_t)n_mom * nv);
    build_F(n_mom, nv, ene, a, b, F.data());

    /* g0_n = M_n * F (bb x nv x n_mom), one BLAS-3 gemm per atom; F is the
     * shared right operand. Each gemm is sizable (bb=nb^2, nv, n_mom in the
     * hundreds), so the threaded BLAS handles parallelism -- keep the atom
     * loop serial to avoid oversubscribing it. Precision follows
     * rsrec_set_precision (prec==1 fp64 bit-reference, else fp32 ~1e-6,
     * consistent with the moment engine). */
    if (c->prec == 1) {
        const cplx one(1.0, 0.0), zero(0.0, 0.0);
        for (int n = 0; n < natoms; ++n)
            zgemm_("N", "N", &bb, &nv, &n_mom, &one,
                   mu + (size_t)bb * n_mom * n, &bb, F.data(), &n_mom,
                   &zero, g0 + (size_t)bb * nv * n, &bb);
        return 0;
    }
    typedef std::complex<float> cf;
    std::vector<cf> Ff(F.size()), Mf((size_t)bb * n_mom), Gf((size_t)bb * nv);
    for (size_t i = 0; i < F.size(); ++i)
        Ff[i] = cf((float)F[i].real(), (float)F[i].imag());
    const cf onef(1.0f, 0.0f), zerof(0.0f, 0.0f);
    for (int n = 0; n < natoms; ++n) {
        const cplx *Mn = mu + (size_t)bb * n_mom * n;
        for (size_t i = 0; i < Mf.size(); ++i)
            Mf[i] = cf((float)Mn[i].real(), (float)Mn[i].imag());
        cgemm_("N", "N", &bb, &nv, &n_mom, &onef, Mf.data(), &bb,
               Ff.data(), &n_mom, &zerof, Gf.data(), &bb);
        cplx *Gn = g0 + (size_t)bb * nv * n;
        for (size_t i = 0; i < Gf.size(); ++i)
            Gn[i] = cplx(Gf[i].real(), Gf[i].imag());
    }
    return 0;
}
