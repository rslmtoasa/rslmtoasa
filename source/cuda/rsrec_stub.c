/* ===========================================================================
 * rsrec_stub.c -- no-op fallback for the rsrec_* recursion backend.
 *
 * cuda/recursion_gpu_mod.f90 is compiled unconditionally and binds to the
 * rsrec_* C interface (see rsrec.h). The real implementations live in
 * cuda/rsrec_gpu.cu (CUDA) and are only built with -DENABLE_CUDA_PLUGIN=ON.
 *
 * For ordinary CPU builds we still need those symbols to link, but we do NOT
 * want to pull in any compute kernels. This stub supplies them: rsrec_create
 * returns NULL, which rsgpu_init() detects (c_associated test) and turns into
 * a clean `error stop` via die()/rsrec_last_error(). Every other entry point
 * returns a nonzero error code. None of these are ever reached unless the
 * user requests the GPU path at runtime (control%gpu_plugin = .true.).
 * =========================================================================== */
#include "rsrec.h"

static const char *g_msg =
    "rsrec: GPU recursion backend not built in this binary. "
    "Reconfigure with -DENABLE_CUDA_PLUGIN=ON to enable the gpu_plugin path.";

const char *rsrec_last_error(void) { return g_msg; }

rsrec_ctx *rsrec_create(int kk, int nb, int nnmax, int ntype, int nmax,
                        int device)
{
    (void)kk; (void)nb; (void)nnmax; (void)ntype; (void)nmax; (void)device;
    return 0; /* NULL ctx -> rsgpu_init() aborts cleanly via die() */
}

void rsrec_destroy(rsrec_ctx *ctx) { (void)ctx; }

int rsrec_set_hamiltonian(rsrec_ctx *ctx, const void *ee, const void *hall,
                          const void *lsham, const int *nn, const int *iz,
                          const void *eeo, const void *hallo, const void *enim)
{
    (void)ctx; (void)ee; (void)hall; (void)lsham; (void)nn; (void)iz;
    (void)eeo; (void)hallo; (void)enim;
    return 1;
}

int rsrec_set_velocity(rsrec_ctx *ctx, const void *v_a, const void *v_b,
                       const void *vo_a, const void *vo_b)
{
    (void)ctx; (void)v_a; (void)v_b; (void)vo_a; (void)vo_b;
    return 1;
}

int rsrec_orbital_moments(rsrec_ctx *ctx, const void *left, const void *psiref,
                          int lld, double a, double b, void *mu)
{
    (void)ctx; (void)left; (void)psiref; (void)lld; (void)a; (void)b; (void)mu;
    return 1;
}

int rsrec_set_grid(rsrec_ctx *ctx, const int *coords, int use_structured)
{
    (void)ctx; (void)coords; (void)use_structured;
    return 1;
}

int rsrec_set_precision(rsrec_ctx *ctx, int prec)
{
    (void)ctx; (void)prec;
    return 1;
}

int rsrec_op_apply(rsrec_ctx *ctx, int which, const void *x, void *y,
                   int nrhs, double a, double b)
{
    (void)ctx; (void)which; (void)x; (void)y; (void)nrhs; (void)a; (void)b;
    return 1;
}

int rsrec_chebyshev_dos(rsrec_ctx *ctx, const void *mu, int n_moments,
                        int natoms, const double *ene, int nv,
                        double a, double b, void *g0)
{
    (void)ctx; (void)mu; (void)n_moments; (void)natoms; (void)ene; (void)nv;
    (void)a; (void)b; (void)g0;
    return 1;
}

int rsrec_block_dos(rsrec_ctx *ctx, const void *a_b, const void *b2_b,
                    const double *a_inf, const double *b_inf,
                    const double *ene, int nv, double eta_re, double eta_im,
                    int natoms, int lld, int sym, void *g0)
{
    (void)ctx; (void)a_b; (void)b2_b; (void)a_inf; (void)b_inf; (void)ene;
    (void)nv; (void)eta_re; (void)eta_im; (void)natoms; (void)lld; (void)sym;
    (void)g0;
    return 1;
}

int rsrec_block_gf_eta(rsrec_ctx *ctx, const void *a_b, const void *b2_b,
                       const double *a_inf, const double *b_inf, double ef,
                       const double *eta_re, const double *eta_im, int n_eta,
                       int natoms, int lld, int sym, void *g0)
{
    (void)ctx; (void)a_b; (void)b2_b; (void)a_inf; (void)b_inf; (void)ef;
    (void)eta_re; (void)eta_im; (void)n_eta; (void)natoms; (void)lld; (void)sym;
    (void)g0;
    return 1;
}

int rsrec_chebyshev_moments(rsrec_ctx *ctx, const void *psi0, int lld,
                            double a, double b, void *mu_out)
{
    (void)ctx; (void)psi0; (void)lld; (void)a; (void)b; (void)mu_out;
    return 1;
}

int rsrec_block_lanczos(rsrec_ctx *ctx, const void *psi0, int lld,
                        void *a_b, void *b2_b, int prec)
{
    (void)ctx; (void)psi0; (void)lld; (void)a_b; (void)b2_b; (void)prec;
    return 1;
}

int rsrec_scalar_lanczos(rsrec_ctx *ctx, int site_j, int lld,
                         double *a_out, double *b2_out)
{
    (void)ctx; (void)site_j; (void)lld; (void)a_out; (void)b2_out;
    return 1;
}

int rsrec_stochastic_moments(rsrec_ctx *ctx, const void *psiref, int lld,
                             double a, double b, void *mu_nm)
{
    (void)ctx; (void)psiref; (void)lld; (void)a; (void)b; (void)mu_nm;
    return 1;
}
