/* ===========================================================================
 * rsrec_cuda.h -- CUDA-oriented recursion plugin ABI
 *
 * Scope:
 *   - chebyshev_recur()
 *   - chebyshev_recur_ij()
 *   - compute_moments_stochastic()
 *   - recur_b()
 *   - recur_b_ij()
 *   - recur()  [nsp == 1 / spin-diagonal path]
 *
 * Backend selection:
 *   0 = csr
 *   1 = bsr
 *   2 = fft   (periodic ee-only path)
 *   3 = conv  (periodic ee-only path)
 *
 * The current implementation keeps the ABI/backend plumbing separate from the
 * numerical kernels. Backends that are not yet specialized still route through
 * the reference rsrec core after validating that the requested configuration is
 * supported.
 * =========================================================================== */
#ifndef RSREC_CUDA_H
#define RSREC_CUDA_H

#include "rsrec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rsrec_cuda_ctx rsrec_cuda_ctx;

enum {
    RSREC_CUDA_BACKEND_CSR = 0,
    RSREC_CUDA_BACKEND_BSR = 1,
    RSREC_CUDA_BACKEND_FFT = 2,
    RSREC_CUDA_BACKEND_CONV = 3
};

int rsrec_cuda_device_count(int *count);
rsrec_cuda_ctx *rsrec_cuda_create(int kk, int nb, int nnmax, int ntype,
                                  int nmax, int device);
void rsrec_cuda_destroy(rsrec_cuda_ctx *ctx);
const char *rsrec_cuda_last_error(void);

int rsrec_cuda_set_backend(rsrec_cuda_ctx *ctx, int backend);
int rsrec_cuda_set_periodic_lattice(rsrec_cuda_ctx *ctx, int pbc, int n1,
                                    int n2, int n3, const double *a,
                                    const double *crd, int nbas);
int rsrec_cuda_set_hamiltonian(rsrec_cuda_ctx *ctx, const void *ee,
                               const void *hall, const void *lsham,
                               const int *nn, const int *iz, const void *eeo,
                               const void *hallo, const void *enim);
int rsrec_cuda_set_velocity(rsrec_cuda_ctx *ctx, const void *v_a,
                            const void *v_b, const void *vo_a,
                            const void *vo_b);
int rsrec_cuda_orbital_moments(rsrec_cuda_ctx *ctx, const void *left,
                               const void *psiref, int lld, double a, double b,
                               void *mu);

int rsrec_cuda_chebyshev_moments(rsrec_cuda_ctx *ctx, const void *psi0,
                                 int lld, double a, double b, void *mu_out);
int rsrec_cuda_block_lanczos(rsrec_cuda_ctx *ctx, const void *psi0, int lld,
                             void *a_b, void *b2_b, int prec);
int rsrec_cuda_scalar_lanczos(rsrec_cuda_ctx *ctx, int site_j, int lld,
                              double *a_out, double *b2_out);
int rsrec_cuda_stochastic_moments(rsrec_cuda_ctx *ctx, const void *psiref,
                                  int lld, double a, double b, void *mu_nm);
int rsrec_cuda_stochastic_profile(rsrec_cuda_ctx *ctx, double *h2d_seconds,
                                  double *cheb_seconds, double *d2h_seconds,
                                  long long *h2d_bytes, long long *d2h_bytes);
int rsrec_cuda_set_precision(rsrec_cuda_ctx *ctx, int prec);
int rsrec_cuda_chebyshev_dos(rsrec_cuda_ctx *ctx, const void *mu, int n_moments,
                             int natoms, const double *ene, int nv,
                             double a, double b, void *g0);
int rsrec_cuda_chebyshev_gf_eta(rsrec_cuda_ctx *ctx, const void *mu,
                                int n_moments, int natoms, const void *F,
                                int n_eta, void *g0);
int rsrec_cuda_block_dos(rsrec_cuda_ctx *ctx, const void *a_b,
                         const void *b2_b, const double *a_inf,
                         const double *b_inf, const double *ene, int nv,
                         double eta_re, double eta_im, int natoms, int lld,
                         int sym, void *g0);
int rsrec_cuda_block_gf_eta(rsrec_cuda_ctx *ctx, const void *a_b,
                            const void *b2_b, const double *a_inf,
                            const double *b_inf, double ef,
                            const double *eta_re, const double *eta_im,
                            int n_eta, int natoms, int lld, int sym,
                            void *g0);

#ifdef __cplusplus
}
#endif

#endif /* RSREC_CUDA_H */
