/*
 * Narrow C ABI for the reciprocal CUDA execution backend.
 *
 * ACC-03 adds a deliberately narrow standard-Hermitian solve entry point.
 * Complex arrays use the native Fortran/CUDA double-complex representation:
 * interleaved (real, imaginary) values in column-major [n,n,batch] order.
 */
#ifndef RSLMTO_RECIPROCAL_CUDA_H
#define RSLMTO_RECIPROCAL_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rslmto_reciprocal_cuda_context rslmto_reciprocal_cuda_context;

int rslmto_reciprocal_cuda_device_count(int *count);
rslmto_reciprocal_cuda_context *rslmto_reciprocal_cuda_create(int device);
const char *rslmto_reciprocal_cuda_last_error(void);

int rslmto_reciprocal_cuda_prepare_operator(
    rslmto_reciprocal_cuda_context *context, int operator_generation);
int rslmto_reciprocal_cuda_solve_zheevd_batch(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    const void *host_hamiltonians,
    double *host_eigenvalues,
    void *host_eigenvectors,
    int request_eigenvectors);
int rslmto_reciprocal_cuda_contract_lehmann(
    rslmto_reciprocal_cuda_context *context,
    int nmat,
    int nk,
    int ne,
    int npair,
    int nblk,
    const double *host_eigenvalues,
    const void *host_eigenvectors,
    const double *host_k_points,
    const void *host_z_contour,
    const double *host_dr,
    const int *host_ioffset,
    const int *host_joffset,
    void *host_blocks,
    double *h2d_seconds,
    double *contraction_seconds,
    double *d2h_seconds);
void rslmto_reciprocal_cuda_get_timings(
    rslmto_reciprocal_cuda_context *context,
    double *h2d_seconds,
    double *solve_seconds,
    double *d2h_seconds,
    int *calls);
void rslmto_reciprocal_cuda_reset_timings(
    rslmto_reciprocal_cuda_context *context);
int rslmto_reciprocal_cuda_get_memory(
    rslmto_reciprocal_cuda_context *context,
    long long *free_bytes,
    long long *total_bytes);
int rslmto_reciprocal_cuda_synchronize(
    rslmto_reciprocal_cuda_context *context);
void rslmto_reciprocal_cuda_destroy(
    rslmto_reciprocal_cuda_context *context);

#ifdef __cplusplus
}
#endif

#endif /* RSLMTO_RECIPROCAL_CUDA_H */
