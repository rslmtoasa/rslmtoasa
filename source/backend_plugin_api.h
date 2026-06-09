#ifndef RSLMTO_BACKEND_PLUGIN_API_H
#define RSLMTO_BACKEND_PLUGIN_API_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

void *rslmto_backend_create(const char *precision_mode, int block_dim, int n_sites, const char *library_name);
int rslmto_backend_destroy(void *backend);
int rslmto_backend_upload_operator(
    void *backend,
    const char *op_name,
    int block_dim,
    int n_sites,
    int nnzb,
    const int *row_ptr,
    const int *col_ind,
    const int *site_types,
    const double _Complex *blocks);
int rslmto_backend_apply_operator(
    void *backend,
    const char *op_name,
    char trans_mode,
    int block_dim,
    int n_sites,
    const double _Complex *psi_in,
    const int *active_in,
    double _Complex *psi_out,
    int *active_out);
int rslmto_backend_run_scalar_lanczos(
    void *backend,
    int llmax,
    int n_targets,
    const int *target_sites,
    int hoh_enabled,
    double *a_out,
    double *b2_out);
int rslmto_backend_run_block_lanczos(
    void *backend,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    double _Complex *a_b_out,
    double _Complex *b2_b_out,
    double *a_diag_out,
    double *b2_diag_out);
int rslmto_backend_run_block_chebyshev(
    void *backend,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    double a_scale,
    double b_shift,
    double _Complex *mu_out);
int rslmto_backend_run_transport_stochastic(
    void *backend,
    int llmax,
    int loop_over,
    int calc_mode,
    int hoh_enabled,
    const int *type_sites,
    double a_scale,
    double b_shift,
    double _Complex *mu_nm_out);
const char *rslmto_backend_name(void);
const char *rslmto_backend_capabilities(void);

#ifdef __cplusplus
}
#endif

#endif
