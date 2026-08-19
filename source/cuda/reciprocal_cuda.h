/*
 * Narrow C ABI for the reciprocal CUDA execution backend.
 *
 * ACC-03 adds a deliberately narrow standard-Hermitian solve entry point.
 * Complex arrays use the native Fortran/CUDA complex representations:
 * interleaved (real, imaginary) values in column-major [n,n,batch] order.
 */
#ifndef RSLMTO_RECIPROCAL_CUDA_H
#define RSLMTO_RECIPROCAL_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rslmto_reciprocal_cuda_context rslmto_reciprocal_cuda_context;

/* Keep the strategy surface deliberately explicit.  FP32 routes are an
 * experimental precision study; they are never selected implicitly. */
enum rslmto_reciprocal_cuda_solver_strategy {
    RSLMTO_RECIPROCAL_CUDA_ZHEEVD_SERIAL = 0,       /* fp64_zheevd */
    RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED = 1,      /* fp64_zheevj_batched */
    RSLMTO_RECIPROCAL_CUDA_CHEEVD_SERIAL = 2,       /* fp32_cheevd */
    RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED = 3       /* fp32_cheevj_batched */
};

int rslmto_reciprocal_cuda_device_count(int *count);
rslmto_reciprocal_cuda_context *rslmto_reciprocal_cuda_create(int device);
const char *rslmto_reciprocal_cuda_last_error(void);

int rslmto_reciprocal_cuda_prepare_operator(
    rslmto_reciprocal_cuda_context *context, int operator_generation);
int rslmto_reciprocal_cuda_set_solver_strategy(
    rslmto_reciprocal_cuda_context *context, int strategy);
/* Returns 0 when the selected strategy can accept this request, 1 when the
 * request is explicitly unsupported, and -1 for invalid/uninitialized state.
 * In particular, cuSOLVER ZheevjBatched requires n <= 32. */
int rslmto_reciprocal_cuda_solver_strategy_supported(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    int request_eigenvectors);
int rslmto_reciprocal_cuda_solve_zheevd_batch(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    const void *host_hamiltonians,
    double *host_eigenvalues,
    void *host_eigenvectors,
    int request_eigenvectors);
/* The FP32 entry point deliberately accepts the original host H64 and returns
 * widened FP64 eigenpairs.  Conversion is owned by the persistent CUDA
 * context, so callers do not allocate conversion buffers per matrix. */
int rslmto_reciprocal_cuda_solve_cheevd_batch(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    const void *host_hamiltonians,
    double *host_eigenvalues,
    void *host_eigenvectors,
    int request_eigenvectors);
/* ACC-11: reuse the eigenpairs left in the persistent solver context.  The
 * context validates the shape and generation token before this entry point
 * can consume its device buffers; no device pointer is exposed to Fortran or
 * physics code. */
int rslmto_reciprocal_cuda_get_resident_token(
    rslmto_reciprocal_cuda_context *context);
int rslmto_reciprocal_cuda_resident_eigensystem_matches(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    int token);
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
int rslmto_reciprocal_cuda_contract_lehmann_resident(
    rslmto_reciprocal_cuda_context *context,
    int nmat,
    int nk,
    int ne,
    int npair,
    int nblk,
    int resident_token,
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
void rslmto_reciprocal_cuda_get_detailed_timings(
    rslmto_reciprocal_cuda_context *context,
    double *host_conversion_seconds,
    double *host_staging_seconds,
    double *h2d_seconds,
    double *solve_seconds,
    double *d2h_seconds,
    double *d2h_values_seconds,
    double *d2h_vectors_seconds,
    double *sync_seconds,
    double *host_widen_seconds,
    double *total_seconds,
    long long *h2d_bytes,
    long long *d2h_values_bytes,
    long long *d2h_vectors_bytes,
    int *pinned_host_active,
    int *calls);
void rslmto_reciprocal_cuda_get_resource_counters(
    rslmto_reciprocal_cuda_context *context,
    long long *cuda_malloc_count,
    long long *cuda_free_count,
    long long *workspace_query_count,
    long long *workspace_reuse_count,
    long long *event_create_count,
    long long *event_destroy_count,
    long long *pinned_alloc_count,
    long long *pinned_free_count);
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
