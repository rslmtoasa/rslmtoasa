/* ===========================================================================
 * rsk.h -- k-space GPU backend scaffold for RS-LMTO-ASA.
 *
 * The first implementation target is batched construction of TB-LMTO H(k)
 * from separate real-space ingredients:
 *   h(k)    = Bloch(ee)
 *   o(k)    = Bloch(eeo)
 *   Hcc(k)  = Bloch(eecc)
 *   H(k)    = h(k) - o(k)*h(k) + Hcc(k) + onsite
 *
 * Keeping Hcc separate from ee preserves the CCOR+HOH invariant validated in
 * the real-space recursion backend.
 * =========================================================================== */
#ifndef RSK_H
#define RSK_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rsk_ctx rsk_ctx;

rsk_ctx *rsk_create(int nb, int nsite, int ntype, int nnmax, int device);
void rsk_destroy(rsk_ctx *ctx);
const char *rsk_last_error(void);

int rsk_set_geometry(rsk_ctx *ctx, const int *site_type, const int *atlist,
                     const int *nn, const int *iz, const double *r_direct);

int rsk_set_operators(rsk_ctx *ctx, const void *ee, const void *eeo,
                      const void *eecc, const void *enim,
                      const void *lsham);

int rsk_build_hk_batch(rsk_ctx *ctx, const double *kpoints, int nk,
                       int order, void *hk_out);

int rsk_diagonalize_hk_batch(rsk_ctx *ctx, int jobz, int generalized,
                             void *evals_out, void *evecs_out);

#ifdef __cplusplus
}
#endif

#endif /* RSK_H */
