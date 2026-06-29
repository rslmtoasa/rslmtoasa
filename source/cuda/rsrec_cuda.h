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
int rsrec_cuda_set_hamiltonian_additive(rsrec_cuda_ctx *ctx,
                                        const void *ee_add,
                                        const void *hall_add);
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

#ifdef __cplusplus
}
#endif

#endif /* RSREC_CUDA_H */
