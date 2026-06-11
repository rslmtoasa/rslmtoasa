/* ===========================================================================
 * rsrec_cpu.cpp -- CPU reference implementation of rsrec.h
 *
 * Purpose: (1) correctness reference for the CUDA backend (identical code
 * structure, routine-for-routine), (2) portable fallback when no GPU is
 * present. OpenMP-parallel; uses LAPACK zheev for the 18x18 eig-sqrt in the
 * block Lanczos, exactly like crecal_b.
 * =========================================================================== */
#include "rsrec.h"

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

typedef std::complex<double> cplx;

extern "C" void zheev_(const char *jobz, const char *uplo, const int *n,
                       cplx *a, const int *lda, double *w, cplx *work,
                       const int *lwork, double *rwork, int *info);
extern "C" void zgemm_(const char *, const char *, const int *, const int *,
                       const int *, const std::complex<double> *,
                       const std::complex<double> *, const int *,
                       const std::complex<double> *, const int *,
                       const std::complex<double> *, std::complex<double> *,
                       const int *);
extern "C" void cgemm_(const char *, const char *, const int *, const int *,
                       const int *, const std::complex<float> *,
                       const std::complex<float> *, const int *,
                       const std::complex<float> *, const int *,
                       const std::complex<float> *, std::complex<float> *,
                       const int *);

static std::string g_err;
extern "C" const char *rsrec_last_error(void) { return g_err.c_str(); }
#define FAIL(msg) do { g_err = (msg); return 1; } while (0)

/* --------------------------------------------------------------------------
 * Context
 * ------------------------------------------------------------------------ */
struct rsrec_ctx {
    int kk = 0, nb = 0, nnmax = 0, ntype = 0, nmax = 0;
    /* Hamiltonian data, Fortran layout, with lsham pre-added to the on-site
     * blocks so the matvec inner loop has one block per shell. */
    std::vector<cplx> ee;     /* (nb,nb,nnmax,ntype): shell 1 = ee+lsham   */
    std::vector<cplx> hall;   /* (nb,nb,nnmax,nmax) : shell 1 = hall+lsham */
    std::vector<cplx> va, vb; /* (nb,nb,nnmax,ntype) velocity operators    */
    std::vector<int> nn;      /* (kk,nnmax) 1-based, nn(k,1)=count          */
    std::vector<int> iz;      /* (kk)       1-based type                    */
    bool have_h = false, have_v = false;
    int cheb_prec = 0;                   /* 0 = fp32 (default), 1 = fp64     */

    /* ---- structured (stencil + correction) path, see rsrec_set_grid ---- */
    bool use_struct = false;
    int t_ref = 0;                       /* 0-based reference type           */
    int N[3] = {0, 0, 0}, lo[3] = {0, 0, 0};
    std::vector<int> cellmap;            /* cell -> 0-based atom or -1       */
    std::vector<int> gcell;              /* atom -> cell index               */
    std::vector<int> roff;               /* (3, nshell_ref) stencil offsets  */
    int nshell_ref = 0;
    std::vector<cplx> negW_h, negW_va, negW_vb;  /* -W_ref per operator     */
    struct CorrEntry { int atom, nbr; const cplx *blk; };
    std::vector<CorrEntry> corr_h, corr_va, corr_vb;  /* per-op corrections */

    /* scratch */
    std::vector<cplx> s0, s1, s2, s3;
};

static inline size_t EE(const rsrec_ctx *c, int l, int m, int s, int t) {
    return (size_t)l + (size_t)c->nb * (m + (size_t)c->nb * (s + (size_t)c->nnmax * t));
}
static inline size_t PSI(const rsrec_ctx *c, int l, int col, int nrhs, int k) {
    return (size_t)l + (size_t)c->nb * (col + (size_t)nrhs * k);
}
static inline int NNF(const rsrec_ctx *c, int k, int s) { /* 0-based k,s */
    return c->nn[(size_t)k + (size_t)c->kk * s];
}

extern "C" rsrec_ctx *rsrec_create(int kk, int nb, int nnmax, int ntype,
                                   int nmax, int /*device*/) {
    if (kk <= 0 || nb <= 0 || nnmax <= 0 || ntype <= 0 || nmax < 0 ||
        nmax > kk) { g_err = "rsrec_create: bad dimensions"; return nullptr; }
    rsrec_ctx *c = new rsrec_ctx();
    c->kk = kk; c->nb = nb; c->nnmax = nnmax; c->ntype = ntype; c->nmax = nmax;
    size_t nfield = (size_t)nb * nb * kk;
    c->s0.resize(nfield); c->s1.resize(nfield);
    c->s2.resize(nfield); c->s3.resize(nfield);
    return c;
}

extern "C" void rsrec_destroy(rsrec_ctx *ctx) { delete ctx; }

extern "C" int rsrec_set_hamiltonian(rsrec_ctx *c, const void *ee_,
                                     const void *hall_, const void *lsham_,
                                     const int *nn, const int *iz) {
    const cplx *ee = (const cplx *)ee_;
    const cplx *hall = (const cplx *)hall_;
    const cplx *lsham = (const cplx *)lsham_;
    if (!ee || !nn || !iz) FAIL("set_hamiltonian: null input");
    if (c->nmax > 0 && !hall) FAIL("set_hamiltonian: nmax>0 but hall is null");

    size_t nee = (size_t)c->nb * c->nb * c->nnmax * c->ntype;
    c->ee.assign(ee, ee + nee);
    if (c->nmax > 0) {
        size_t nh = (size_t)c->nb * c->nb * c->nnmax * c->nmax;
        c->hall.assign(hall, hall + nh);
    }
    c->nn.assign(nn, nn + (size_t)c->kk * c->nnmax);
    c->iz.assign(iz, iz + c->kk);

    /* Fold lsham into the on-site (shell-1) blocks once. */
    if (lsham) {
        for (int t = 0; t < c->ntype; ++t)
            for (int m = 0; m < c->nb; ++m)
                for (int l = 0; l < c->nb; ++l)
                    c->ee[EE(c, l, m, 0, t)] +=
                        lsham[(size_t)l + (size_t)c->nb * (m + (size_t)c->nb * t)];
        for (int k = 0; k < c->nmax; ++k) {
            int t = c->iz[k] - 1;
            for (int m = 0; m < c->nb; ++m)
                for (int l = 0; l < c->nb; ++l)
                    c->hall[EE(c, l, m, 0, k)] +=
                        lsham[(size_t)l + (size_t)c->nb * (m + (size_t)c->nb * t)];
        }
    }
    c->have_h = true;
    return 0;
}

extern "C" int rsrec_set_velocity(rsrec_ctx *c, const void *va_, const void *vb_) {
    const cplx *va = (const cplx *)va_, *vb = (const cplx *)vb_;
    if (!va || !vb) FAIL("set_velocity: null input");
    size_t n = (size_t)c->nb * c->nb * c->nnmax * c->ntype;
    c->va.assign(va, va + n);
    c->vb.assign(vb, vb + n);
    c->have_v = true;
    return 0;
}

/* --------------------------------------------------------------------------
 * Core block-ELL matvec:  y_k = sum_shells B(k,s) x_{nn(k,s)}  (s=1 onsite)
 * then y = (y - b x)/a.  `which`: 0=H, 1=v_a, 2=v_b.
 * The CUDA backend implements this same loop as one fused kernel.
 * ------------------------------------------------------------------------ */
static void op_apply_ell(rsrec_ctx *c, int which, const cplx *x, cplx *y,
                         int nrhs, double a, double b) {
    const int nb = c->nb;
#pragma omp parallel for schedule(dynamic, 64)
    for (int k = 0; k < c->kk; ++k) {
        const bool imp = (which == 0) && (k < c->nmax);
        const int t = c->iz[k] - 1;
        const cplx *blocks =
            (which == 1) ? c->va.data() :
            (which == 2) ? c->vb.data() :
            (imp ? c->hall.data() : c->ee.data());
        const int slot = imp ? k : t;   /* 4th index of the block array */

        cplx *yk = y + PSI(c, 0, 0, nrhs, k);
        for (int i = 0; i < nb * nrhs; ++i) yk[i] = cplx(0.0, 0.0);

        const int nr = NNF(c, k, 0);
        for (int s = 0; s < nr; ++s) {           /* s=0 is the on-site shell */
            int nbr;
            if (s == 0) nbr = k;
            else {
                nbr = NNF(c, k, s) - 1;          /* to 0-based; -1 => absent */
                if (nbr < 0) continue;
            }
            const cplx *B = blocks + EE(c, 0, 0, s, slot);
            const cplx *xn = x + PSI(c, 0, 0, nrhs, nbr);
            /* yk(:,col) += B * xn(:,col) */
            for (int col = 0; col < nrhs; ++col)
                for (int m = 0; m < nb; ++m) {
                    const cplx xm = xn[m + nb * col];
                    if (xm == cplx(0.0, 0.0)) continue;
                    const cplx *Bm = B + (size_t)nb * m;
                    cplx *yc = yk + nb * col;
                    for (int l = 0; l < nb; ++l) yc[l] += Bm[l] * xm;
                }
        }
        if (b != 0.0 || a != 1.0) {
            const cplx *xk = x + PSI(c, 0, 0, nrhs, k);
            const double inva = 1.0 / a;
            for (int i = 0; i < nb * nrhs; ++i)
                yk[i] = (yk[i] - b * xk[i]) * inva;
        }
    }
}

/* --------------------------------------------------------------------------
 * Structured path: y = Stencil_ref(x) + Correction(x), then (y - b x)/a.
 *
 * This is the KPM conv/FFT backend from the Python prototypes mapped onto
 * the RS-LMTO data model. The reference stencil is the per-type bulk
 * Hamiltonian ee(:,:,:,t_ref) applied as a LINEAR (zero-padded) convolution
 * over the lattice grid: vacuum/outside cells carry psi = 0, so open
 * boundaries are automatic. Every row that deviates from the stencil --
 * the hall impurity region, atoms of type /= t_ref, irregular neighbour
 * patterns -- is repaired by a sparse correction  dH = H_true - H_stencil
 * built once at rsrec_set_grid. On the CPU the convolution is evaluated
 * directly (identical algebra to the cuFFT version in rsrec_gpu.cu, which
 * this file validates).
 * ------------------------------------------------------------------------ */
static void op_apply_struct(rsrec_ctx *c, int which, const cplx *x, cplx *y,
                            int nrhs, double a, double b) {
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb;
    const cplx *Wref =
        (which == 1) ? c->va.data() + EE(c, 0, 0, 0, c->t_ref) :
        (which == 2) ? c->vb.data() + EE(c, 0, 0, 0, c->t_ref) :
                       c->ee.data() + EE(c, 0, 0, 0, c->t_ref);
    const int Nx = c->N[0], Ny = c->N[1], Nz = c->N[2];

    /* Stencil sweep (== zero-padded convolution gathered at atom cells)    */
#pragma omp parallel for schedule(static)
    for (int k = 0; k < c->kk; ++k) {
        cplx *yk = y + PSI(c, 0, 0, nrhs, k);
        for (int i = 0; i < nb * nrhs; ++i) yk[i] = cplx(0.0, 0.0);
        int cell = c->gcell[k];
        int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
        for (int s = 0; s < c->nshell_ref; ++s) {
            int tx = cx + c->roff[3 * s], ty = cy + c->roff[3 * s + 1];
            int tz = cz + c->roff[3 * s + 2];
            if (tx < 0 || tx >= Nx || ty < 0 || ty >= Ny ||
                tz < 0 || tz >= Nz) continue;
            int nbr = c->cellmap[(tx * Ny + ty) * (size_t)Nz + tz];
            if (nbr < 0) continue;                    /* vacuum cell        */
            const cplx *B = Wref + bb * s;
            const cplx *xn = x + PSI(c, 0, 0, nrhs, nbr);
            for (int col = 0; col < nrhs; ++col)
                for (int m = 0; m < nb; ++m) {
                    const cplx xm = xn[m + nb * col];
                    if (xm == cplx(0.0, 0.0)) continue;
                    const cplx *Bm = B + (size_t)nb * m;
                    cplx *yc = yk + nb * col;
                    for (int l = 0; l < nb; ++l) yc[l] += Bm[l] * xm;
                }
        }
    }
    /* Sparse correction dH = H_true - H_stencil on flagged rows            */
    const auto &corr = (which == 1) ? c->corr_va :
                       (which == 2) ? c->corr_vb : c->corr_h;
    for (const auto &e : corr) {                /* serial: few entries, and
                                                   atoms repeat -> races    */
        cplx *yk = y + PSI(c, 0, 0, nrhs, e.atom);
        const cplx *xn = x + PSI(c, 0, 0, nrhs, e.nbr);
        for (int col = 0; col < nrhs; ++col)
            for (int m = 0; m < nb; ++m) {
                const cplx xm = xn[m + nb * col];
                if (xm == cplx(0.0, 0.0)) continue;
                const cplx *Bm = e.blk + (size_t)nb * m;
                cplx *yc = yk + nb * col;
                for (int l = 0; l < nb; ++l) yc[l] += Bm[l] * xm;
            }
    }
    if (b != 0.0 || a != 1.0) {
        const double inva = 1.0 / a;
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)((size_t)nb * nrhs * c->kk); ++i)
            y[i] = (y[i] - b * x[i]) * inva;
    }
}

static void op_apply(rsrec_ctx *c, int which, const cplx *x, cplx *y,
                     int nrhs, double a, double b) {
    if (c->use_struct)
        op_apply_struct(c, which, x, y, nrhs, a, b);
    else
        op_apply_ell(c, which, x, y, nrhs, a, b);
}

/* --------------------------------------------------------------------------
 * rsrec_set_grid: map the clust onto a regular lattice grid, derive the
 * reference stencil and build the correction lists. coords: (3, kk) integer
 * lattice coordinates (Fortran order). use_structured = 0 reverts to ELL.
 * ------------------------------------------------------------------------ */
static void grid_corr_for_op(rsrec_ctx *c, const std::vector<char> &flag,
                             const cplx *blocks_type, const cplx *blocks_imp,
                             bool imp_rows, std::vector<cplx> &negW,
                             std::vector<rsrec_ctx::CorrEntry> &corr) {
    const int nb = c->nb;
    const size_t bb = (size_t)nb * nb;
    const int Nx = c->N[0], Ny = c->N[1], Nz = c->N[2];
    /* materialize -W_ref once */
    negW.resize(bb * c->nshell_ref);
    const cplx *Wref = blocks_type + EE(c, 0, 0, 0, c->t_ref);
    for (size_t i = 0; i < negW.size(); ++i) negW[i] = -Wref[i];
    corr.clear();
    for (int k = 0; k < c->kk; ++k) {
        if (!flag[k]) continue;
        /* + true row */
        const bool imp = imp_rows && (k < c->nmax);
        const cplx *base = imp ? blocks_imp : blocks_type;
        const int slot = imp ? k : (c->iz[k] - 1);
        const int nr = NNF(c, k, 0);
        for (int s = 0; s < nr; ++s) {
            int nbr = (s == 0) ? k : NNF(c, k, s) - 1;
            if (nbr < 0) continue;
            corr.push_back({k, nbr, base + EE(c, 0, 0, s, slot)});
        }
        /* - stencil row */
        int cell = c->gcell[k];
        int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
        for (int s = 0; s < c->nshell_ref; ++s) {
            int tx = cx + c->roff[3 * s], ty = cy + c->roff[3 * s + 1];
            int tz = cz + c->roff[3 * s + 2];
            if (tx < 0 || tx >= Nx || ty < 0 || ty >= Ny ||
                tz < 0 || tz >= Nz) continue;
            int nbr = c->cellmap[(tx * Ny + ty) * (size_t)Nz + tz];
            if (nbr < 0) continue;
            corr.push_back({k, nbr, negW.data() + bb * s});
        }
    }
}

extern "C" int rsrec_set_grid(rsrec_ctx *c, const int *coords,
                              int use_structured) {
    if (!c->have_h) FAIL("set_grid: call set_hamiltonian first");
    if (!use_structured) { c->use_struct = false; return 0; }
    const int kk = c->kk;

    /* bounding box and cell map */
    for (int d = 0; d < 3; ++d) {
        int mn = coords[d], mx = coords[d];
        for (int k = 1; k < kk; ++k) {
            mn = std::min(mn, coords[3 * k + d]);
            mx = std::max(mx, coords[3 * k + d]);
        }
        c->lo[d] = mn; c->N[d] = mx - mn + 1;
    }
    const int Nx = c->N[0], Ny = c->N[1], Nz = c->N[2];
    c->cellmap.assign((size_t)Nx * Ny * Nz, -1);
    c->gcell.resize(kk);
    for (int k = 0; k < kk; ++k) {
        int cx = coords[3 * k] - c->lo[0];
        int cy = coords[3 * k + 1] - c->lo[1];
        int cz = coords[3 * k + 2] - c->lo[2];
        size_t cell = (cx * (size_t)Ny + cy) * Nz + cz;
        if (c->cellmap[cell] >= 0)
            FAIL("set_grid: two atoms share one lattice cell (multi-atom "
                 "bases need a basis-resolved grid; not supported yet)");
        c->cellmap[cell] = k;
        c->gcell[k] = (int)cell;
    }

    /* reference type = majority among bulk atoms; template = bulk atom of
     * that type with the largest shell count                                */
    std::vector<int> cnt(c->ntype, 0);
    for (int k = c->nmax; k < kk; ++k) cnt[c->iz[k] - 1]++;
    if (kk == c->nmax)
        FAIL("set_grid: hall covers all atoms; the structured path has no "
             "bulk stencil to exploit -- stay on the ELL backend");
    c->t_ref = (int)(std::max_element(cnt.begin(), cnt.end()) - cnt.begin());
    int k0 = -1, best = -1;
    for (int k = c->nmax; k < kk; ++k)
        if (c->iz[k] - 1 == c->t_ref && NNF(c, k, 0) > best) {
            best = NNF(c, k, 0); k0 = k;
        }
    if (k0 < 0) FAIL("set_grid: no bulk atom of the reference type");
    c->nshell_ref = best;
    c->roff.assign(3 * c->nshell_ref, 0);
    for (int s = 1; s < c->nshell_ref; ++s) {
        int nbr = NNF(c, k0, s) - 1;
        if (nbr < 0) FAIL("set_grid: template atom has a hole in its "
                          "neighbour list; pick a cluster interior template");
        for (int d = 0; d < 3; ++d)
            c->roff[3 * s + d] = coords[3 * nbr + d] - coords[3 * k0 + d];
    }

    /* flag every row the stencil cannot represent                          */
    auto shells_match = [&](int k) -> bool {
        const int nr = NNF(c, k, 0);
        if (nr > c->nshell_ref) return false;
        int cell = c->gcell[k];
        int cz = cell % Nz, cy = (cell / Nz) % Ny, cx = cell / (Nz * Ny);
        for (int s = 1; s < c->nshell_ref; ++s) {
            int tx = cx + c->roff[3 * s], ty = cy + c->roff[3 * s + 1];
            int tz = cz + c->roff[3 * s + 2];
            bool inbox = (tx >= 0 && tx < Nx && ty >= 0 && ty < Ny &&
                          tz >= 0 && tz < Nz);
            int occ = inbox ? c->cellmap[(tx * (size_t)Ny + ty) * Nz + tz] : -1;
            int nbr = (s < nr) ? NNF(c, k, s) - 1 : -1;
            if (occ != nbr) return false;   /* stencil and nn disagree       */
        }
        return true;
    };
    std::vector<char> flagH(kk), flagV(kk);
    size_t nflagH = 0, nflagV = 0;
    for (int k = 0; k < kk; ++k) {
        bool irregular = !shells_match(k) || (c->iz[k] - 1 != c->t_ref);
        flagV[k] = irregular;
        flagH[k] = irregular || (k < c->nmax);
        nflagH += flagH[k]; nflagV += flagV[k];
    }
    grid_corr_for_op(c, flagH, c->ee.data(), c->hall.data(), true,
                     c->negW_h, c->corr_h);
    if (c->have_v) {
        grid_corr_for_op(c, flagV, c->va.data(), nullptr, false,
                         c->negW_va, c->corr_va);
        grid_corr_for_op(c, flagV, c->vb.data(), nullptr, false,
                         c->negW_vb, c->corr_vb);
    }
    fprintf(stderr,
            "rsrec_set_grid: box %dx%dx%d, ref type %d, %d stencil shells, "
            "corrected rows: H %zu/%d, V %zu/%d\n",
            Nx, Ny, Nz, c->t_ref + 1, c->nshell_ref, nflagH, kk,
            c->have_v ? nflagV : 0, kk);
    c->use_struct = true;
    return 0;
}

extern "C" int rsrec_op_apply(rsrec_ctx *c, int which, const void *x, void *y,
                              int nrhs, double a, double b) {
    if (!c->have_h) FAIL("op_apply: Hamiltonian not set");
    if (which != 0 && !c->have_v) FAIL("op_apply: velocity operators not set");
    op_apply(c, which, (const cplx *)x, (cplx *)y, nrhs, a, b);
    return 0;
}

/* --------------------------------------------------------------------------
 * Block reductions:  C += sum_k A_k^H B_k   and  C = sum_k A_k^H B_k
 * (zgemm('c','n') per atom in the Fortran code)
 * ------------------------------------------------------------------------ */
static void block_dot(const rsrec_ctx *c, const cplx *A, const cplx *B,
                      cplx *C /* nb*nb, zeroed by caller */) {
    const int nb = c->nb;
#pragma omp parallel
    {
        std::vector<cplx> loc((size_t)nb * nb, cplx(0.0, 0.0));
#pragma omp for schedule(static) nowait
        for (int k = 0; k < c->kk; ++k) {
            const cplx *Ak = A + PSI(c, 0, 0, nb, k);
            const cplx *Bk = B + PSI(c, 0, 0, nb, k);
            for (int j = 0; j < nb; ++j)
                for (int i = 0; i < nb; ++i) {
                    cplx s(0.0, 0.0);
                    const cplx *ai = Ak + nb * i, *bj = Bk + nb * j;
                    for (int l = 0; l < nb; ++l) s += std::conj(ai[l]) * bj[l];
                    loc[i + nb * j] += s;
                }
        }
#pragma omp critical
        for (size_t i = 0; i < loc.size(); ++i) C[i] += loc[i];
    }
}

/* psi(:,:,k) <- psi(:,:,k) * M  for all atoms (right block multiply) */
static void right_mul(const rsrec_ctx *c, cplx *psi, const cplx *M, cplx *tmp) {
    const int nb = c->nb;
#pragma omp parallel for schedule(static)
    for (int k = 0; k < c->kk; ++k) {
        cplx *pk = psi + PSI(c, 0, 0, nb, k);
        cplx *tk = tmp + PSI(c, 0, 0, nb, k);
        for (int j = 0; j < nb; ++j)
            for (int l = 0; l < nb; ++l) {
                cplx s(0.0, 0.0);
                for (int m = 0; m < nb; ++m) s += pk[l + nb * m] * M[m + nb * j];
                tk[l + nb * j] = s;
            }
        std::memcpy(pk, tk, sizeof(cplx) * nb * nb);
    }
}

/* B = U sqrt(ev) U^H and Bi = U ev^{-1/2} U^H from Hermitian S (= crecal_b) */
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
            cplx sb(0.0, 0.0), sbi(0.0, 0.0);
            for (int l = 0; l < nb; ++l) {
                double lam = std::sqrt(std::max(0.0, ev[l]));
                cplx uu = U[i + (size_t)nb * l] * std::conj(U[j + (size_t)nb * l]);
                sb += uu * lam;
                if (lam > 0.0) sbi += uu / lam;
            }
            B[i + (size_t)nb * j] = sb;
            Bi[i + (size_t)nb * j] = sbi;
        }
    return 0;
}

extern "C" int rsrec_set_precision(rsrec_ctx *c, int prec) {
    if (prec < 0 || prec > 1) FAIL("set_precision: 0 = fp32, 1 = fp64");
    c->cheb_prec = prec;
    return 0;
}

/* --------------------------------------------------------------------------
 * FP32 block-Chebyshev path: exact algorithmic mirror of the CUDA engine
 * (site-major packed psi, scaled fp32 Hamiltonian copy, fused recurrence,
 * raw moments combined in double on the host). Validates the GPU math and
 * doubles as a faster CPU option.
 * ------------------------------------------------------------------------ */
typedef std::complex<float> cflt;

static int cheb_moments_f32_cpu(rsrec_ctx *c, const void *psi0_, int lld,
                                double a, double b, void *mu_) {
    const int nb = c->nb, kk = c->kk, nnmax = c->nnmax;
    const size_t bb = (size_t)nb * nb;
    const size_t ld = (size_t)nb * kk;
    const double inva = 1.0 / a;

    /* scaled fp32 Hamiltonian copies (lsham already folded in ee/hall)     */
    auto cvt = [&](const std::vector<cplx> &in, std::vector<cflt> &out) {
        out.resize(in.size());
        for (size_t i = 0; i < in.size(); ++i) {
            int l = (int)(i % nb), m = (int)((i / nb) % nb);
            int s = (int)((i / bb) % nnmax);
            double re = in[i].real(), im = in[i].imag();
            if (l == m && s == 0) re -= b;
            out[i] = cflt((float)(re * inva), (float)(im * inva));
        }
    };
    std::vector<cflt> fee, fhall;
    cvt(c->ee, fee);
    if (c->nmax > 0) cvt(c->hall, fhall);

    /* pack psi0 site-major fp32                                            */
    const cplx *p0d = (const cplx *)psi0_;
    std::vector<cflt> p0(ld * nb), p1(ld * nb), p2(ld * nb);
    for (int k = 0; k < kk; ++k)
        for (int col = 0; col < nb; ++col)
            for (int l = 0; l < nb; ++l)
                p0[(size_t)l + nb * k + ld * col] =
                    cflt((float)p0d[l + nb * ((size_t)col + nb * k)].real(),
                         (float)p0d[l + nb * ((size_t)col + nb * k)].imag());

    auto step = [&](const std::vector<cflt> &x1, const std::vector<cflt> &x0,
                    std::vector<cflt> &y, float alpha, float beta) {
#pragma omp parallel for schedule(static)
        for (int k = 0; k < kk; ++k) {
            const cflt *base = (k < c->nmax)
                ? fhall.data() + bb * (size_t)nnmax * k
                : fee.data() + bb * (size_t)nnmax * (c->iz[k] - 1);
            const int nr = NNF(c, k, 0);
            for (int col = 0; col < nb; ++col) {
                for (int l = 0; l < nb; ++l) {
                    cflt acc(0.0f, 0.0f);
                    for (int s = 0; s < nr; ++s) {
                        const int nbr = (s == 0) ? k : NNF(c, k, s) - 1;
                        if (nbr < 0) continue;
                        const cflt *B = base + bb * s;
                        const cflt *xv = x1.data() + (size_t)nb * nbr +
                                         ld * col;
                        for (int m = 0; m < nb; ++m)
                            acc += B[l + nb * m] * xv[m];
                    }
                    size_t idx = (size_t)l + nb * k + ld * col;
                    y[idx] = (beta != 0.0f)
                                 ? alpha * acc + beta * x0[idx]
                                 : alpha * acc;
                }
            }
        }
    };
    /* fp32 block moments C = A^H B over the packed layout                  */
    auto fdot = [&](const std::vector<cflt> &A, const std::vector<cflt> &B,
                    cflt *C) {
        for (int j = 0; j < nb; ++j)
            for (int i = 0; i < nb; ++i) {
                cflt s(0.0f, 0.0f);
                const cflt *ai = A.data() + ld * i, *bj = B.data() + ld * j;
                for (size_t r = 0; r < ld; ++r) s += std::conj(ai[r]) * bj[r];
                C[i + (size_t)nb * j] = s;
            }
    };

    const size_t nmom = 2 * (size_t)lld + 2;
    std::vector<cflt> raw(nmom * bb);
    fdot(p0, p0, raw.data());
    step(p0, p0, p1, 1.0f, 0.0f);
    fdot(p0, p1, raw.data() + bb);
    for (int ll = 1; ll <= lld; ++ll) {
        step(p1, p0, p2, 2.0f, -1.0f);
        fdot(p1, p1, raw.data() + (size_t)(2 * ll) * bb);
        fdot(p2, p1, raw.data() + (size_t)(2 * ll + 1) * bb);
        std::swap(p0, p1); std::swap(p1, p2);
    }
    cplx *mu = (cplx *)mu_;
    for (size_t i = 0; i < bb; ++i) {
        mu[i] = cplx(raw[i].real(), raw[i].imag());
        mu[bb + i] = cplx(raw[bb + i].real(), raw[bb + i].imag());
    }
    for (int ll = 1; ll <= lld; ++ll)
        for (size_t i = 0; i < bb; ++i) {
            cplx d1(raw[(size_t)(2 * ll) * bb + i].real(),
                    raw[(size_t)(2 * ll) * bb + i].imag());
            cplx d2(raw[(size_t)(2 * ll + 1) * bb + i].real(),
                    raw[(size_t)(2 * ll + 1) * bb + i].imag());
            mu[(size_t)(2 * ll) * bb + i] = 2.0 * d1 - mu[i];
            mu[(size_t)(2 * ll + 1) * bb + i] = 2.0 * d2 - mu[bb + i];
        }
    return 0;
}

/* --------------------------------------------------------------------------
 * Chebyshev block moments (mirrors cheb_0th_mom / cheb_1st_mom /
 * chebyshev_recur_ll with the double-moment trick)
 * ------------------------------------------------------------------------ */
static int cheb_moments_f64_cpu(rsrec_ctx *c, const void *psi0_, int lld,
                                double a, double b, void *mu_) {
    if (!c->have_h) FAIL("chebyshev_moments: Hamiltonian not set");
    const cplx *psiref = (const cplx *)psi0_;
    cplx *mu = (cplx *)mu_;
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb;
    const size_t nfield = blk * c->kk;

    cplx *p0 = c->s0.data(), *p1 = c->s1.data(), *p2 = c->s2.data();
    std::memcpy(p0, psiref, sizeof(cplx) * nfield);
    std::memset(mu, 0, sizeof(cplx) * blk * (2 * (size_t)lld + 2));

    /* mu_1 = sum_k psiref^H psi0 */
    block_dot(c, psiref, p0, mu + 0 * blk);
    /* psi1 = (H psi0 - b psi0)/a ; mu_2 = sum psiref^H psi1 */
    op_apply(c, 0, p0, p1, nb, a, b);
    block_dot(c, psiref, p1, mu + 1 * blk);

    for (int ll = 1; ll <= lld; ++ll) {
        /* psi2 = 2*(H~ psi1) - psi0 */
        op_apply(c, 0, p1, p2, nb, a, b);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)nfield; ++i)
            p2[i] = 2.0 * p2[i] - p0[i];

        cplx *m1 = mu + (size_t)(2 * ll) * blk;       /* mu(2ll+1), 0-based */
        cplx *m2 = mu + (size_t)(2 * ll + 1) * blk;   /* mu(2ll+2)          */
        block_dot(c, p1, p1, m1);                     /* dum1 = psi1^H psi1 */
        block_dot(c, p2, p1, m2);                     /* dum2 = psi2^H psi1 */
        for (size_t i = 0; i < blk; ++i) {
            m1[i] = 2.0 * m1[i] - mu[i];              /* 2 dum1 - mu_1 */
            m2[i] = 2.0 * m2[i] - mu[blk + i];        /* 2 dum2 - mu_2 */
        }
        /* shift registers */
        cplx *t = p0; p0 = p1; p1 = p2; p2 = t;
    }
    return 0;
}

extern "C" int rsrec_chebyshev_moments(rsrec_ctx *c, const void *psi0_,
                                       int lld, double a, double b,
                                       void *mu_) {
    if (!c->have_h) FAIL("chebyshev_moments: Hamiltonian not set");
    if (c->cheb_prec == 0 && !c->use_struct)
        return cheb_moments_f32_cpu(c, psi0_, lld, a, b, mu_);
    return cheb_moments_f64_cpu(c, psi0_, lld, a, b, mu_);
}

/* --------------------------------------------------------------------------
 * Block Lanczos (mirrors crecal_b / hop_b)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_block_lanczos(rsrec_ctx *c, const void *psi0_, int lld,
                                   void *a_b_, void *b2_b_) {
    if (!c->have_h) FAIL("block_lanczos: Hamiltonian not set");
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb;
    const size_t nfield = blk * c->kk;
    cplx *a_b = (cplx *)a_b_, *b2_b = (cplx *)b2_b_;

    cplx *psi = c->s0.data(), *pmn = c->s1.data();
    cplx *hpsi = c->s2.data(), *tmp = c->s3.data();
    std::memcpy(psi, psi0_, sizeof(cplx) * nfield);
    std::memset(pmn, 0, sizeof(cplx) * nfield);

    std::vector<cplx> sum_b(blk, cplx(0.0, 0.0));
    for (int i = 0; i < nb; ++i) sum_b[i + (size_t)nb * i] = 1.0;  /* b2(1)=I */
    std::vector<cplx> B(blk), Bi(blk), An(blk);

    std::memset(a_b, 0, sizeof(cplx) * blk * lld);
    std::memset(b2_b, 0, sizeof(cplx) * blk * lld);

    for (int ll = 0; ll < lld - 1; ++ll) {
        /* hop_b: hpsi = H psi ; A_n = sum psi^H hpsi ; pmn = hpsi - pmn   */
        op_apply(c, 0, psi, hpsi, nb, 1.0, 0.0);
        std::fill(An.begin(), An.end(), cplx(0.0, 0.0));
        block_dot(c, psi, hpsi, An.data());
        std::memcpy(a_b + (size_t)ll * blk, An.data(), sizeof(cplx) * blk);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)nfield; ++i) pmn[i] = hpsi[i] - pmn[i];

        /* crecal_b: b2(ll) = sum_b ; pmn -= psi*A_n ; sum_b = pmn^H pmn   */
        std::memcpy(b2_b + (size_t)ll * blk, sum_b.data(), sizeof(cplx) * blk);
#pragma omp parallel for schedule(static)
        for (int k = 0; k < c->kk; ++k) {
            cplx *pk = pmn + PSI(c, 0, 0, nb, k);
            const cplx *sk = psi + PSI(c, 0, 0, nb, k);
            for (int j = 0; j < nb; ++j)
                for (int l = 0; l < nb; ++l) {
                    cplx s(0.0, 0.0);
                    for (int m = 0; m < nb; ++m)
                        s += sk[l + nb * m] * An[m + (size_t)nb * j];
                    pk[l + nb * j] -= s;
                }
        }
        std::fill(sum_b.begin(), sum_b.end(), cplx(0.0, 0.0));
        block_dot(c, pmn, pmn, sum_b.data());

        /* B_{n+1} = sqrt(sum_b) by zheev, exactly like crecal_b           */
        if (eig_sqrt(nb, sum_b.data(), B.data(), Bi.data()) != 0)
            FAIL("block_lanczos: zheev failed");

        /* psi_{n+1} = pmn * B^{-1};  pmn <- psi_n * B (for next subtract) */
        std::memcpy(tmp, psi, sizeof(cplx) * nfield);   /* psi_t = psi_n  */
        right_mul(c, pmn, Bi.data(), hpsi);             /* pmn  = pmn*Bi  */
        cplx *t = psi; psi = pmn; pmn = t;              /* swap            */
        /* now pmn holds old psi_t? no: pmn holds old psi (swapped). We
         * stored psi_t in tmp; recover: pmn = tmp * B                     */
        std::memcpy(pmn, tmp, sizeof(cplx) * nfield);
        right_mul(c, pmn, B.data(), hpsi);
    }
    std::memcpy(b2_b + (size_t)(lld - 1) * blk, sum_b.data(),
                sizeof(cplx) * blk);
    return 0;
}

/* --------------------------------------------------------------------------
 * Scalar Lanczos batched over the nb orbital chains of one site
 * (mirrors recur()/crecal()/hop(); each column is an independent chain)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_scalar_lanczos(rsrec_ctx *c, int site_j, int lld,
                                    double *a_out, double *b2_out) {
    if (!c->have_h) FAIL("scalar_lanczos: Hamiltonian not set");
    if (site_j < 1 || site_j > c->kk) FAIL("scalar_lanczos: bad site");
    const int nb = c->nb, j0 = site_j - 1;
    const size_t nfield = (size_t)nb * nb * c->kk;

    cplx *psi = c->s0.data(), *pmn = c->s1.data(), *hpsi = c->s2.data();
    std::memset(psi, 0, sizeof(cplx) * nfield);
    std::memset(pmn, 0, sizeof(cplx) * nfield);
    for (int l = 0; l < nb; ++l)               /* column l = e_l on site j  */
        psi[PSI(c, l, l, nb, j0)] = 1.0;

    std::vector<double> summ(nb, 1.0);          /* b2temp(1) = 1            */
    for (int i = 0; i < lld * nb; ++i) { a_out[i] = 0.0; b2_out[i] = 0.0; }

    for (int ll = 0; ll < lld - 1; ++ll) {
        op_apply(c, 0, psi, hpsi, nb, 1.0, 0.0);
        /* per-column: a(ll,col) = Re <psi_col | hpsi_col>; pmn += hpsi    */
        std::vector<double> acol(nb, 0.0);
#pragma omp parallel
        {
            std::vector<double> loc(nb, 0.0);
#pragma omp for schedule(static) nowait
            for (int k = 0; k < c->kk; ++k)
                for (int col = 0; col < nb; ++col)
                    for (int l = 0; l < nb; ++l) {
                        size_t id = PSI(c, l, col, nb, k);
                        loc[col] += std::real(std::conj(psi[id]) * hpsi[id]);
                        /* pmn updated below to keep this loop reduction-only */
                    }
#pragma omp critical
            for (int col = 0; col < nb; ++col) acol[col] += loc[col];
        }
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)nfield; ++i) pmn[i] += hpsi[i];

        for (int col = 0; col < nb; ++col) {
            a_out[ll + (size_t)lld * col] = acol[col];
            b2_out[ll + (size_t)lld * col] = summ[col];
        }
        /* pmn(:,col,:) -= a_col * psi(:,col,:); summ_col = |pmn_col|^2     */
        std::vector<double> s2(nb, 0.0);
#pragma omp parallel
        {
            std::vector<double> loc(nb, 0.0);
#pragma omp for schedule(static) nowait
            for (int k = 0; k < c->kk; ++k)
                for (int col = 0; col < nb; ++col)
                    for (int l = 0; l < nb; ++l) {
                        size_t id = PSI(c, l, col, nb, k);
                        pmn[id] -= acol[col] * psi[id];
                        loc[col] += std::norm(pmn[id]);
                    }
#pragma omp critical
            for (int col = 0; col < nb; ++col) s2[col] += loc[col];
        }
        /* psi <- pmn/sqrt(summ); pmn <- -psi_old*sqrt(summ)               */
#pragma omp parallel for schedule(static)
        for (int k = 0; k < c->kk; ++k)
            for (int col = 0; col < nb; ++col) {
                double si = 1.0 / std::sqrt(s2[col]), sf = std::sqrt(s2[col]);
                for (int l = 0; l < nb; ++l) {
                    size_t id = PSI(c, l, col, nb, k);
                    cplx newpsi = pmn[id] * si;
                    pmn[id] = -psi[id] * sf;
                    psi[id] = newpsi;
                }
            }
        for (int col = 0; col < nb; ++col) summ[col] = s2[col];
    }
    for (int col = 0; col < nb; ++col)
        b2_out[(lld - 1) + (size_t)lld * col] = summ[col];
    return 0;
}

/* --------------------------------------------------------------------------
 * Stochastic conductivity moments (non-hoh path of compute_moments_stochastic)
 * ------------------------------------------------------------------------ */
extern "C" int rsrec_stochastic_moments(rsrec_ctx *c, const void *psiref_,
                                        int lld, double a, double b,
                                        void *mu_nm_) {
    if (!c->have_h) FAIL("stochastic_moments: Hamiltonian not set");
    if (!c->have_v) FAIL("stochastic_moments: velocity operators not set");
    const cplx *psiref = (const cplx *)psiref_;
    cplx *mu = (cplx *)mu_nm_;
    const int nb = c->nb;
    const size_t blk = (size_t)nb * nb, nfield = blk * c->kk;

    /* Left states L_m = T_{m-1}(H~)|psiref>, all kept resident.            */
    std::vector<cplx> left((size_t)lld * nfield);
    {
        cplx *w0 = c->s0.data(), *w1 = c->s1.data(), *w2 = c->s2.data();
        std::memcpy(w1, psiref, sizeof(cplx) * nfield);
        std::memcpy(left.data(), w1, sizeof(cplx) * nfield);
        for (int m = 2; m <= lld; ++m) {
            if (m == 2) {
                std::memcpy(w0, w1, sizeof(cplx) * nfield);
                op_apply(c, 0, w0, w1, nb, a, b);
            } else {
                op_apply(c, 0, w1, w2, nb, a, b);
#pragma omp parallel for schedule(static)
                for (long i = 0; i < (long)nfield; ++i)
                    w2[i] = 2.0 * w2[i] - w0[i];
                cplx *t = w0; w0 = w1; w1 = w2; w2 = t;
            }
            std::memcpy(left.data() + (size_t)(m - 1) * nfield, w1,
                        sizeof(cplx) * nfield);
        }
    }

    /* Right recursion: v0 = V_b psiref; v_n = T_{n-1}(H~) v0; R = V_a v_n. */
    std::vector<cplx> rightbuf(nfield);
    cplx *v0 = c->s0.data(), *v1 = c->s1.data(), *v2 = c->s2.data();
    cplx *R = c->s3.data();
    op_apply(c, 2, psiref, v0, nb, 1.0, 0.0);     /* V_b |psiref> */
    std::memset(mu, 0, sizeof(cplx) * blk * (size_t)lld * lld);

    for (int n = 1; n <= lld; ++n) {
        if (n == 1) {
            std::memcpy(v1, v0, sizeof(cplx) * nfield);
        } else if (n == 2) {
            std::memcpy(v0, v1, sizeof(cplx) * nfield);
            op_apply(c, 0, v0, v1, nb, a, b);
        } else {
            op_apply(c, 0, v1, v2, nb, a, b);
#pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)nfield; ++i)
                v2[i] = 2.0 * v2[i] - v0[i];
            cplx *t = v0; v0 = v1; v1 = v2; v2 = t;
        }
        op_apply(c, 1, v1, R, nb, 1.0, 0.0);      /* V_a T_{n-1} v0 */

        /* mu(:,:,n,m) = sum_k L_m(k)^H R(k)  for all m                     */
#pragma omp parallel for schedule(dynamic)
        for (int m = 0; m < lld; ++m) {
            cplx *mum = mu + blk * ((size_t)(n - 1) + (size_t)lld * m);
            const cplx *Lm = left.data() + (size_t)m * nfield;
            for (int k = 0; k < c->kk; ++k) {
                const cplx *Ak = Lm + PSI(c, 0, 0, nb, k);
                const cplx *Bk = R + PSI(c, 0, 0, nb, k);
                for (int j = 0; j < nb; ++j)
                    for (int i = 0; i < nb; ++i) {
                        cplx s(0.0, 0.0);
                        for (int l = 0; l < nb; ++l)
                            s += std::conj(Ak[l + nb * i]) * Bk[l + nb * j];
                        mum[i + nb * j] += s;
                    }
            }
        }
    }
    return 0;
}

/* --------------------------------------------------------------------------
 * Chebyshev Green-function / DOS reconstruction (mirrors chebyshev_green in
 * bands.f90, including its Jackson kernel convention):
 *
 *   g0(:,:,ie,n) = sum_i mu(:,:,i,n) * F(i,ie)
 *   F(i,ie) = g_J(i) * c_i * (-i) e^{-i (i-1) arccos(x_ie)}
 *                                       / sqrt(a^2 - (E_ie - b)^2)
 *   x = (E - b)/a,  c_1 = 1, c_i>=2 = 2,
 *   g_J(i) = [(N-(i-1)+1) cos(th_i) + sin(th_i)/tan(pi/(N+1))]/(N+1),
 *   th_i = pi (i-1)/(N+1),  N = n_moments.
 *
 * Evaluated as ONE GEMM per atom (G_n = M_n * F): the scalar exp/acos
 * triple loop of the Fortran routine becomes a level-3 BLAS call; the
 * unused t_polynomial evaluation is dropped. DOS = -Im g0 / pi as usual.
 * Precision follows rsrec_set_precision (fp32 default, matching the
 * moment engine; fp64 bit-comparable with the Fortran loop).
 * Energies outside the (a, b) window produce NaN exactly like the CPU
 * code (sqrt of a negative argument) -- keep the grid inside the window.
 * ------------------------------------------------------------------------ */
static void build_F(int N, int nv, const double *ene, double a, double b,
                    cplx *F) {
    const double PI = 3.14159265358979323846;
    const double cot = sin(PI / (N + 1.0)) > 0.0
                           ? cos(PI / (N + 1.0)) / sin(PI / (N + 1.0)) : 0.0;
#pragma omp parallel for schedule(static)
    for (int ie = 0; ie < nv; ++ie) {
        double x = (ene[ie] - b) / a;
        double th = acos(x);
        double pref = 1.0 / sqrt(a * a - (ene[ie] - b) * (ene[ie] - b));
        for (int i = 0; i < N; ++i) {
            double thl = PI * i / (N + 1.0);
            double gj = ((N - i + 1) * cos(thl) + sin(thl) * cot) / (N + 1.0);
            double cc = (i == 0) ? 1.0 : 2.0;
            double ang = i * th;             /* e^{-i ang} * (-i)           */
            double fr = -sin(ang), fi = -cos(ang);
            F[(size_t)i + (size_t)N * ie] = gj * cc * pref * cplx(fr, fi);
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

    if (c->cheb_prec == 1) {                 /* fp64: one zgemm per atom    */
        const cplx one(1.0, 0.0), zero(0.0, 0.0);
        for (int n = 0; n < natoms; ++n)
            zgemm_("N", "N", &bb, &nv, &n_mom, &one,
                   mu + (size_t)bb * n_mom * n, &bb, F.data(), &n_mom,
                   &zero, g0 + (size_t)bb * nv * n, &bb);
        return 0;
    }
    /* fp32 (default): mirrors the CUDA cgemm3m path                        */
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
