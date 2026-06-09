#include "backend_plugin_api.h"

#include <complex.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <cmath>
#include <cstdlib>
#include <cstring>

extern "C" void zheev_(const char *jobz, const char *uplo, const int *n, cuDoubleComplex *a, const int *lda, double *w,
                       cuDoubleComplex *work, const int *lwork, double *rwork, int *info);

typedef struct {
    char *name;
    int block_dim;
    int n_sites;
    int nnzb;
    int *row_ptr;
    int *col_ind;
    int *site_types;
    cuDoubleComplex *blocks;
    int scalar_dim;
    int scalar_nnz;
    int *csr_row_ptr;
    int *csr_col_ind;
    cuDoubleComplex *csr_values;
    int *d_csr_row_ptr;
    int *d_csr_col_ind;
    cuDoubleComplex *d_csr_values;
    cusparseSpMatDescr_t spmat;
} plugin_operator;

typedef struct {
    char *precision_mode;
    char *library_name;
    int block_dim;
    int n_sites;
    int scalar_dim;
    int capacity;
    int count;
    plugin_operator *ops;
    cusparseHandle_t cusparse_handle;
    cublasHandle_t cublas_handle;
    cuDoubleComplex *d_state_in;
    cuDoubleComplex *d_state_out;
    int *d_active_in;
    int *d_active_out;
    cuDoubleComplex *d_x;
    cuDoubleComplex *d_y;
    cusparseDnVecDescr_t vec_x;
    cusparseDnVecDescr_t vec_y;
    void *d_spmv_buffer;
    size_t spmv_buffer_size;
    size_t state_len;
} plugin_backend;

extern "C" int rslmto_backend_destroy(void *backend_handle);

static __host__ __device__ cuDoubleComplex c_zero(void) {
    return make_cuDoubleComplex(0.0, 0.0);
}

static __host__ __device__ cuDoubleComplex c_one(void) {
    return make_cuDoubleComplex(1.0, 0.0);
}

static __host__ __device__ cuDoubleComplex c_real(double x) {
    return make_cuDoubleComplex(x, 0.0);
}

static __host__ __device__ cuDoubleComplex c_scale(cuDoubleComplex z, double s) {
    return make_cuDoubleComplex(cuCreal(z) * s, cuCimag(z) * s);
}

static __host__ __device__ cuDoubleComplex c_sub(cuDoubleComplex a, cuDoubleComplex b) {
    return make_cuDoubleComplex(cuCreal(a) - cuCreal(b), cuCimag(a) - cuCimag(b));
}

static char *dup_string(const char *input) {
    size_t len = 0U;
    char *copy = NULL;
    if (input == NULL) {
        return NULL;
    }
    len = strlen(input);
    copy = (char *) malloc(len + 1U);
    if (copy == NULL) {
        return NULL;
    }
    memcpy(copy, input, len + 1U);
    return copy;
}

static int cuda_ok(cudaError_t status) {
    return status == cudaSuccess ? 0 : 1;
}

static int cusparse_ok(cusparseStatus_t status) {
    return status == CUSPARSE_STATUS_SUCCESS ? 0 : 1;
}

static int cublas_ok(cublasStatus_t status) {
    return status == CUBLAS_STATUS_SUCCESS ? 0 : 1;
}

static size_t state_len_of(const plugin_backend *backend) {
    return (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) backend->n_sites;
}

static size_t site_block_offset(const plugin_backend *backend, int site) {
    return (size_t) site * (size_t) backend->block_dim * (size_t) backend->block_dim;
}

static size_t matrix_offset(int n, int i, int j) {
    return (size_t) i + (size_t) n * (size_t) j;
}

static size_t scalar_coeff_offset(int llmax, int nb, int ll, int orb, int wf) {
    return (size_t) ll + (size_t) llmax * (size_t) orb + (size_t) llmax * (size_t) nb * (size_t) wf;
}

static size_t a_b_offset(int nb, int llmax, int n_workflow, int i, int j, int ll, int wf) {
    return (size_t) i + (size_t) nb * (size_t) j + (size_t) nb * (size_t) nb * (size_t) ll +
           (size_t) nb * (size_t) nb * (size_t) llmax * (size_t) wf;
}

static size_t mu_offset(int nb, int nmom, int i, int j, int n, int wf) {
    return (size_t) i + (size_t) nb * (size_t) j + (size_t) nb * (size_t) nb * (size_t) n +
           (size_t) nb * (size_t) nb * (size_t) nmom * (size_t) wf;
}

static size_t mu_nm_offset(int nb, int llmax, int i, int j, int n, int m, int wf) {
    return (size_t) i + (size_t) nb * (size_t) j + (size_t) nb * (size_t) nb * (size_t) n +
           (size_t) nb * (size_t) nb * (size_t) llmax * (size_t) m +
           (size_t) nb * (size_t) nb * (size_t) llmax * (size_t) llmax * (size_t) wf;
}

static size_t moment_state_offset(size_t state_len, int idx) {
    return (size_t) idx * state_len;
}

static void state_zero(const plugin_backend *backend, cuDoubleComplex *state) {
    memset(state, 0, state_len_of(backend) * sizeof(*state));
}

static void mask_zero(const plugin_backend *backend, int *mask) {
    memset(mask, 0, (size_t) backend->n_sites * sizeof(*mask));
}

static void matrix_zero(cuDoubleComplex *mat, int n) {
    memset(mat, 0, (size_t) n * (size_t) n * sizeof(*mat));
}

static void matrix_identity(cuDoubleComplex *mat, int n) {
    int i = 0;
    matrix_zero(mat, n);
    for (i = 0; i < n; ++i) {
        mat[matrix_offset(n, i, i)] = c_one();
    }
}

static void matrix_copy(cuDoubleComplex *dst, const cuDoubleComplex *src, int n) {
    memcpy(dst, src, (size_t) n * (size_t) n * sizeof(*dst));
}

static void state_copy(const plugin_backend *backend, cuDoubleComplex *dst, const cuDoubleComplex *src) {
    memcpy(dst, src, state_len_of(backend) * sizeof(*dst));
}

static void matrix_right_multiply(cuDoubleComplex *out, const cuDoubleComplex *left, const cuDoubleComplex *right, int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            cuDoubleComplex accum = c_zero();
            for (k = 0; k < n; ++k) {
                accum = cuCadd(accum, cuCmul(left[matrix_offset(n, i, k)], right[matrix_offset(n, k, j)]));
            }
            out[matrix_offset(n, i, j)] = accum;
        }
    }
}

static void accumulate_left_dagger_right(cuDoubleComplex *accum, const cuDoubleComplex *left, const cuDoubleComplex *right, int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            cuDoubleComplex entry = c_zero();
            for (k = 0; k < n; ++k) {
                entry = cuCadd(entry, cuCmul(cuConj(left[matrix_offset(n, k, i)]), right[matrix_offset(n, k, j)]));
            }
            accum[matrix_offset(n, i, j)] = cuCadd(accum[matrix_offset(n, i, j)], entry);
        }
    }
}

static void scale_shift_state(const plugin_backend *backend, cuDoubleComplex *out, const cuDoubleComplex *in, double a_scale, double b_shift) {
    size_t idx = 0U;
    size_t total = state_len_of(backend);
    for (idx = 0U; idx < total; ++idx) {
        out[idx] = c_scale(c_sub(out[idx], c_scale(in[idx], b_shift)), 1.0 / a_scale);
    }
}

static void active_copy(const plugin_backend *backend, int *dst, const int *src) {
    memcpy(dst, src, (size_t) backend->n_sites * sizeof(*dst));
}

static void active_union(const plugin_backend *backend, int *dst, const int *a, const int *b) {
    int i = 0;
    for (i = 0; i < backend->n_sites; ++i) {
        dst[i] = (a[i] != 0 || b[i] != 0) ? 1 : 0;
    }
}

static double state_norm2(const plugin_backend *backend, const cuDoubleComplex *state) {
    size_t idx = 0U;
    size_t total = state_len_of(backend);
    double accum = 0.0;
    for (idx = 0U; idx < total; ++idx) {
        accum += cuCreal(cuCmul(cuConj(state[idx]), state[idx]));
    }
    return accum;
}

static cuDoubleComplex state_inner_product(const plugin_backend *backend, const cuDoubleComplex *left, const cuDoubleComplex *right) {
    size_t idx = 0U;
    size_t total = state_len_of(backend);
    cuDoubleComplex accum = c_zero();
    for (idx = 0U; idx < total; ++idx) {
        accum = cuCadd(accum, cuCmul(cuConj(left[idx]), right[idx]));
    }
    return accum;
}

static void free_operator(plugin_operator *op) {
    if (op == NULL) {
        return;
    }
    free(op->name);
    free(op->row_ptr);
    free(op->col_ind);
    free(op->site_types);
    free(op->blocks);
    free(op->csr_row_ptr);
    free(op->csr_col_ind);
    free(op->csr_values);
    if (op->spmat != NULL) {
        cusparseDestroySpMat(op->spmat);
    }
    if (op->d_csr_row_ptr != NULL) {
        cudaFree(op->d_csr_row_ptr);
    }
    if (op->d_csr_col_ind != NULL) {
        cudaFree(op->d_csr_col_ind);
    }
    if (op->d_csr_values != NULL) {
        cudaFree(op->d_csr_values);
    }
    memset(op, 0, sizeof(*op));
}

static plugin_operator *find_operator(plugin_backend *backend, const char *name) {
    int i = 0;
    if (backend == NULL || name == NULL) {
        return NULL;
    }
    for (i = 0; i < backend->count; ++i) {
        if (backend->ops[i].name != NULL && strcmp(backend->ops[i].name, name) == 0) {
            return &backend->ops[i];
        }
    }
    return NULL;
}

static plugin_operator *upsert_operator(plugin_backend *backend, const char *name) {
    plugin_operator *op = NULL;
    plugin_operator *grown = NULL;

    op = find_operator(backend, name);
    if (op != NULL) {
        free_operator(op);
        op->name = dup_string(name);
        return op;
    }

    if (backend->count == backend->capacity) {
        int next_capacity = backend->capacity == 0 ? 8 : 2 * backend->capacity;
        grown = (plugin_operator *) realloc(backend->ops, (size_t) next_capacity * sizeof(*grown));
        if (grown == NULL) {
            return NULL;
        }
        memset(grown + backend->capacity, 0, (size_t) (next_capacity - backend->capacity) * sizeof(*grown));
        backend->ops = grown;
        backend->capacity = next_capacity;
    }

    op = &backend->ops[backend->count++];
    memset(op, 0, sizeof(*op));
    op->name = dup_string(name);
    return op;
}

static int ensure_spmv_buffer(plugin_backend *backend, plugin_operator *op) {
    size_t buffer_size_n = 0U;
    size_t buffer_size_c = 0U;
    size_t needed = 0U;
    cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);

    if (backend == NULL || op == NULL || op->spmat == NULL) {
        return 1;
    }

    if (cusparse_ok(cusparseSpMV_bufferSize(
            backend->cusparse_handle,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha,
            op->spmat,
            backend->vec_x,
            &beta,
            backend->vec_y,
            CUDA_C_64F,
            CUSPARSE_SPMV_CSR_ALG2,
            &buffer_size_n)) != 0) {
        return 1;
    }
    if (cusparse_ok(cusparseSpMV_bufferSize(
            backend->cusparse_handle,
            CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE,
            &alpha,
            op->spmat,
            backend->vec_x,
            &beta,
            backend->vec_y,
            CUDA_C_64F,
            CUSPARSE_SPMV_CSR_ALG2,
            &buffer_size_c)) != 0) {
        return 1;
    }

    needed = buffer_size_n > buffer_size_c ? buffer_size_n : buffer_size_c;
    if (needed <= backend->spmv_buffer_size && backend->d_spmv_buffer != NULL) {
        return 0;
    }

    if (backend->d_spmv_buffer != NULL) {
        cudaFree(backend->d_spmv_buffer);
        backend->d_spmv_buffer = NULL;
        backend->spmv_buffer_size = 0U;
    }
    if (needed == 0U) {
        return 0;
    }
    if (cuda_ok(cudaMalloc(&backend->d_spmv_buffer, needed)) != 0) {
        return 1;
    }
    backend->spmv_buffer_size = needed;
    return 0;
}

static int ensure_operator_spmat(plugin_operator *op) {
    if (op == NULL) {
        return 1;
    }
    if (op->spmat != NULL) {
        return 0;
    }
    if (cusparse_ok(cusparseCreateCsr(
            &op->spmat,
            op->scalar_dim,
            op->scalar_dim,
            op->scalar_nnz,
            op->d_csr_row_ptr,
            op->d_csr_col_ind,
            op->d_csr_values,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO,
            CUDA_C_64F)) != 0) {
        return 1;
    }
    return 0;
}

__global__ static void pack_block_column_kernel(
    int block_dim,
    int n_sites,
    int column,
    const cuDoubleComplex *state_in,
    const int *active_in,
    cuDoubleComplex *x_out) {
    int tid = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int scalar_dim = block_dim * n_sites;

    if (tid >= scalar_dim) {
        return;
    }

    {
        int site = tid / block_dim;
        int row = tid - site * block_dim;
        size_t state_idx = (size_t) row +
                           (size_t) column * (size_t) block_dim +
                           (size_t) site * (size_t) block_dim * (size_t) block_dim;
        if (active_in[site] != 0) {
            x_out[tid] = state_in[state_idx];
        } else {
            x_out[tid] = c_zero();
        }
    }
}

__global__ static void unpack_block_column_kernel(
    int block_dim,
    int n_sites,
    int column,
    const cuDoubleComplex *y_in,
    cuDoubleComplex *state_out) {
    int tid = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int scalar_dim = block_dim * n_sites;

    if (tid >= scalar_dim) {
        return;
    }

    {
        int site = tid / block_dim;
        int row = tid - site * block_dim;
        size_t state_idx = (size_t) row +
                           (size_t) column * (size_t) block_dim +
                           (size_t) site * (size_t) block_dim * (size_t) block_dim;
        state_out[state_idx] = y_in[tid];
    }
}

__global__ static void detect_active_sites_kernel(
    int block_dim,
    int n_sites,
    const cuDoubleComplex *state_out,
    int *active_out) {
    int site = (int) (blockIdx.x * blockDim.x + threadIdx.x);

    if (site >= n_sites) {
        return;
    }

    {
        int touched = 0;
        int idx = 0;
        size_t base = (size_t) site * (size_t) block_dim * (size_t) block_dim;
        for (idx = 0; idx < block_dim * block_dim; ++idx) {
            cuDoubleComplex value = state_out[base + (size_t) idx];
            if (cuCreal(value) != 0.0 || cuCimag(value) != 0.0) {
                touched = 1;
                break;
            }
        }
        active_out[site] = touched;
    }
}

static int apply_named_operator(
    plugin_backend *backend,
    const char *name,
    char trans_mode,
    const cuDoubleComplex *psi_in,
    const int *active_in,
    cuDoubleComplex *psi_out,
    int *active_out) {
    plugin_operator *op = NULL;
    cusparseOperation_t op_mode = CUSPARSE_OPERATION_NON_TRANSPOSE;
    cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
    int column = 0;
    int threads = 256;
    int vector_blocks = 0;
    int site_blocks = 0;

    if (backend == NULL || name == NULL || psi_in == NULL || active_in == NULL || psi_out == NULL || active_out == NULL) {
        return 1;
    }

    op = find_operator(backend, name);
    if (op == NULL) {
        return 1;
    }
    if (ensure_operator_spmat(op) != 0) {
        return 1;
    }
    if (ensure_spmv_buffer(backend, op) != 0) {
        return 1;
    }

    if (trans_mode == 'c' || trans_mode == 'C') {
        op_mode = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    } else if (trans_mode == 'n' || trans_mode == 'N') {
        op_mode = CUSPARSE_OPERATION_NON_TRANSPOSE;
    } else {
        return 1;
    }

    vector_blocks = (backend->scalar_dim + threads - 1) / threads;
    site_blocks = (backend->n_sites + threads - 1) / threads;

    if (cuda_ok(cudaMemcpy(backend->d_state_in, psi_in, backend->state_len * sizeof(*backend->d_state_in), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(backend->d_active_in, active_in, (size_t) backend->n_sites * sizeof(*backend->d_active_in), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemset(backend->d_state_out, 0, backend->state_len * sizeof(*backend->d_state_out))) != 0 ||
        cuda_ok(cudaMemset(backend->d_active_out, 0, (size_t) backend->n_sites * sizeof(*backend->d_active_out))) != 0) {
        return 1;
    }

    for (column = 0; column < backend->block_dim; ++column) {
        pack_block_column_kernel<<<vector_blocks, threads>>>(
            backend->block_dim,
            backend->n_sites,
            column,
            backend->d_state_in,
            backend->d_active_in,
            backend->d_x);
        if (cuda_ok(cudaGetLastError()) != 0) {
            return 1;
        }

        if (cusparse_ok(cusparseSpMV(
                backend->cusparse_handle,
                op_mode,
                &alpha,
                op->spmat,
                backend->vec_x,
                &beta,
                backend->vec_y,
                CUDA_C_64F,
                CUSPARSE_SPMV_CSR_ALG2,
                backend->d_spmv_buffer)) != 0) {
            return 1;
        }

        unpack_block_column_kernel<<<vector_blocks, threads>>>(
            backend->block_dim,
            backend->n_sites,
            column,
            backend->d_y,
            backend->d_state_out);
        if (cuda_ok(cudaGetLastError()) != 0) {
            return 1;
        }
    }

    detect_active_sites_kernel<<<site_blocks, threads>>>(
        backend->block_dim,
        backend->n_sites,
        backend->d_state_out,
        backend->d_active_out);
    if (cuda_ok(cudaGetLastError()) != 0 || cuda_ok(cudaDeviceSynchronize()) != 0) {
        return 1;
    }

    if (cuda_ok(cudaMemcpy(psi_out, backend->d_state_out, backend->state_len * sizeof(*backend->d_state_out), cudaMemcpyDeviceToHost)) != 0 ||
        cuda_ok(cudaMemcpy(active_out, backend->d_active_out, (size_t) backend->n_sites * sizeof(*backend->d_active_out), cudaMemcpyDeviceToHost)) != 0) {
        return 1;
    }

    return 0;
}

static int apply_hamiltonian(
    plugin_backend *backend,
    int hoh_enabled,
    const cuDoubleComplex *psi_in,
    const int *active_in,
    double a_scale,
    double b_shift,
    cuDoubleComplex *psi_out,
    int *active_out) {
    const char *name = hoh_enabled != 0 ? "hamiltonian_hoh" : "hamiltonian";
    if (apply_named_operator(backend, name, 'n', psi_in, active_in, psi_out, active_out) != 0) {
        return 1;
    }
    scale_shift_state(backend, psi_out, psi_in, a_scale, b_shift);
    return 0;
}

static int apply_velocity_hoh_named(
    plugin_backend *backend,
    const char *velocity_name,
    const char *overlap_name,
    const cuDoubleComplex *psi_in,
    const int *active_in,
    cuDoubleComplex *psi_out,
    int *active_out) {
    cuDoubleComplex *tmp_h = NULL;
    cuDoubleComplex *tmp_v = NULL;
    cuDoubleComplex *tmp_vo = NULL;
    int *tmp_active_v = NULL;
    int *tmp_active_h = NULL;
    int *tmp_active_vo = NULL;
    size_t idx = 0U;

    tmp_h = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*tmp_h));
    tmp_v = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*tmp_v));
    tmp_vo = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*tmp_vo));
    tmp_active_v = (int *) malloc((size_t) backend->n_sites * sizeof(*tmp_active_v));
    tmp_active_h = (int *) malloc((size_t) backend->n_sites * sizeof(*tmp_active_h));
    tmp_active_vo = (int *) malloc((size_t) backend->n_sites * sizeof(*tmp_active_vo));
    if (tmp_h == NULL || tmp_v == NULL || tmp_vo == NULL || tmp_active_v == NULL || tmp_active_h == NULL || tmp_active_vo == NULL) {
        free(tmp_h);
        free(tmp_v);
        free(tmp_vo);
        free(tmp_active_v);
        free(tmp_active_h);
        free(tmp_active_vo);
        return 1;
    }

    if (apply_named_operator(backend, velocity_name, 'n', psi_in, active_in, tmp_v, tmp_active_v) != 0 ||
        apply_named_operator(backend, "hamiltonian_core", 'n', psi_in, active_in, tmp_h, tmp_active_h) != 0 ||
        apply_named_operator(backend, overlap_name, 'n', tmp_h, tmp_active_h, tmp_vo, tmp_active_vo) != 0) {
        free(tmp_h);
        free(tmp_v);
        free(tmp_vo);
        free(tmp_active_v);
        free(tmp_active_h);
        free(tmp_active_vo);
        return 1;
    }

    for (idx = 0U; idx < state_len_of(backend); ++idx) {
        psi_out[idx] = c_sub(tmp_v[idx], tmp_vo[idx]);
    }
    active_union(backend, active_out, tmp_active_v, tmp_active_vo);

    free(tmp_h);
    free(tmp_v);
    free(tmp_vo);
    free(tmp_active_v);
    free(tmp_active_h);
    free(tmp_active_vo);
    return 0;
}

static int hermitian_sqrt_and_inverse(int n, const cuDoubleComplex *b2, cuDoubleComplex *b, cuDoubleComplex *b_inv) {
    cuDoubleComplex *u = NULL;
    cuDoubleComplex *work = NULL;
    cuDoubleComplex work_query;
    double *evals = NULL;
    double *rwork = NULL;
    cuDoubleComplex *tmp = NULL;
    cuDoubleComplex *tmp2 = NULL;
    cuDoubleComplex *lam = NULL;
    cuDoubleComplex *lam_inv = NULL;
    int info = 0;
    int lwork = -1;
    int i = 0;
    char jobz = 'V';
    char uplo = 'U';

    u = (cuDoubleComplex *) malloc((size_t) n * (size_t) n * sizeof(*u));
    evals = (double *) malloc((size_t) n * sizeof(*evals));
    rwork = (double *) malloc((size_t) (3 * n - 2) * sizeof(*rwork));
    tmp = (cuDoubleComplex *) malloc((size_t) n * (size_t) n * sizeof(*tmp));
    tmp2 = (cuDoubleComplex *) malloc((size_t) n * (size_t) n * sizeof(*tmp2));
    lam = (cuDoubleComplex *) malloc((size_t) n * (size_t) n * sizeof(*lam));
    lam_inv = (cuDoubleComplex *) malloc((size_t) n * (size_t) n * sizeof(*lam_inv));
    if (u == NULL || evals == NULL || rwork == NULL || tmp == NULL || tmp2 == NULL || lam == NULL || lam_inv == NULL) {
        free(u);
        free(evals);
        free(rwork);
        free(tmp);
        free(tmp2);
        free(lam);
        free(lam_inv);
        return 1;
    }

    matrix_copy(u, b2, n);
    zheev_(&jobz, &uplo, &n, u, &n, evals, &work_query, &lwork, rwork, &info);
    if (info != 0) {
        free(u);
        free(evals);
        free(rwork);
        free(tmp);
        free(tmp2);
        free(lam);
        free(lam_inv);
        return 1;
    }

    lwork = (int) cuCreal(work_query);
    if (lwork < 2 * n) {
        lwork = 2 * n;
    }
    work = (cuDoubleComplex *) malloc((size_t) lwork * sizeof(*work));
    if (work == NULL) {
        free(u);
        free(evals);
        free(rwork);
        free(tmp);
        free(tmp2);
        free(lam);
        free(lam_inv);
        return 1;
    }

    matrix_copy(u, b2, n);
    zheev_(&jobz, &uplo, &n, u, &n, evals, work, &lwork, rwork, &info);
    if (info != 0) {
        free(u);
        free(evals);
        free(rwork);
        free(tmp);
        free(tmp2);
        free(lam);
        free(lam_inv);
        free(work);
        return 1;
    }

    matrix_zero(lam, n);
    matrix_zero(lam_inv, n);
    for (i = 0; i < n; ++i) {
        double ev = evals[i] > 0.0 ? evals[i] : 0.0;
        double root = sqrt(ev);
        lam[matrix_offset(n, i, i)] = c_real(root);
        lam_inv[matrix_offset(n, i, i)] = root > 0.0 ? c_real(1.0 / root) : c_zero();
    }

    matrix_right_multiply(tmp, u, lam, n);
    {
        int col = 0;
        int row = 0;
        int k = 0;
        for (col = 0; col < n; ++col) {
            for (row = 0; row < n; ++row) {
                cuDoubleComplex accum = c_zero();
                for (k = 0; k < n; ++k) {
                    accum = cuCadd(accum, cuCmul(tmp[matrix_offset(n, row, k)], cuConj(u[matrix_offset(n, col, k)])));
                }
                b[matrix_offset(n, row, col)] = accum;
            }
        }
    }

    matrix_right_multiply(tmp2, u, lam_inv, n);
    {
        int col = 0;
        int row = 0;
        int k = 0;
        for (col = 0; col < n; ++col) {
            for (row = 0; row < n; ++row) {
                cuDoubleComplex accum = c_zero();
                for (k = 0; k < n; ++k) {
                    accum = cuCadd(accum, cuCmul(tmp2[matrix_offset(n, row, k)], cuConj(u[matrix_offset(n, col, k)])));
                }
                b_inv[matrix_offset(n, row, col)] = accum;
            }
        }
    }

    free(u);
    free(evals);
    free(rwork);
    free(tmp);
    free(tmp2);
    free(lam);
    free(lam_inv);
    free(work);
    return 0;
}

static void set_workflow_state(plugin_backend *backend, int site_i, int site_j, int variant_id, cuDoubleComplex *state, int *active) {
    double inv_sqrt2 = 0.70710678118654752440;
    cuDoubleComplex a_sign = c_zero();
    cuDoubleComplex b_sign = c_zero();
    int l = 0;

    state_zero(backend, state);
    mask_zero(backend, active);

    switch (variant_id) {
        case 1:
            a_sign = c_real(inv_sqrt2);
            b_sign = c_real(inv_sqrt2);
            break;
        case 2:
            a_sign = c_real(inv_sqrt2);
            b_sign = c_real(-inv_sqrt2);
            break;
        case 3:
            a_sign = c_real(inv_sqrt2);
            b_sign = make_cuDoubleComplex(0.0, inv_sqrt2);
            break;
        default:
            a_sign = c_real(inv_sqrt2);
            b_sign = make_cuDoubleComplex(0.0, -inv_sqrt2);
            break;
    }

    if (site_i == site_j) {
        a_sign = c_one();
        b_sign = c_one();
    }

    active[site_i - 1] = 1;
    active[site_j - 1] = 1;
    for (l = 0; l < backend->block_dim; ++l) {
        state[site_block_offset(backend, site_i - 1) + matrix_offset(backend->block_dim, l, l)] = a_sign;
        state[site_block_offset(backend, site_j - 1) + matrix_offset(backend->block_dim, l, l)] = b_sign;
    }
}

extern "C" void *rslmto_backend_create(const char *precision_mode, int block_dim, int n_sites, const char *library_name) {
    plugin_backend *backend = NULL;
    int scalar_dim = 0;
    size_t state_len = 0U;

    if (precision_mode == NULL || strcmp(precision_mode, "complex_fp64") != 0) {
        return NULL;
    }
    if (block_dim <= 0 || n_sites <= 0) {
        return NULL;
    }

    backend = (plugin_backend *) calloc(1U, sizeof(*backend));
    if (backend == NULL) {
        return NULL;
    }

    backend->precision_mode = dup_string(precision_mode);
    backend->library_name = dup_string(library_name != NULL ? library_name : "cuda");
    backend->block_dim = block_dim;
    backend->n_sites = n_sites;
    backend->scalar_dim = block_dim * n_sites;
    scalar_dim = backend->scalar_dim;
    state_len = (size_t) block_dim * (size_t) block_dim * (size_t) n_sites;
    backend->state_len = state_len;

    if (backend->precision_mode == NULL || backend->library_name == NULL) {
        rslmto_backend_destroy(backend);
        return NULL;
    }

    if (cusparse_ok(cusparseCreate(&backend->cusparse_handle)) != 0 ||
        cublas_ok(cublasCreate(&backend->cublas_handle)) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_state_in, state_len * sizeof(*backend->d_state_in))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_state_out, state_len * sizeof(*backend->d_state_out))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_active_in, (size_t) n_sites * sizeof(*backend->d_active_in))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_active_out, (size_t) n_sites * sizeof(*backend->d_active_out))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_x, (size_t) scalar_dim * sizeof(*backend->d_x))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_y, (size_t) scalar_dim * sizeof(*backend->d_y))) != 0 ||
        cusparse_ok(cusparseCreateDnVec(&backend->vec_x, scalar_dim, backend->d_x, CUDA_C_64F)) != 0 ||
        cusparse_ok(cusparseCreateDnVec(&backend->vec_y, scalar_dim, backend->d_y, CUDA_C_64F)) != 0) {
        rslmto_backend_destroy(backend);
        return NULL;
    }

    return backend;
}

extern "C" int rslmto_backend_destroy(void *backend_handle) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    int i = 0;

    if (backend == NULL) {
        return 0;
    }

    for (i = 0; i < backend->count; ++i) {
        free_operator(&backend->ops[i]);
    }
    free(backend->ops);
    free(backend->precision_mode);
    free(backend->library_name);
    if (backend->vec_x != NULL) {
        cusparseDestroyDnVec(backend->vec_x);
    }
    if (backend->vec_y != NULL) {
        cusparseDestroyDnVec(backend->vec_y);
    }
    if (backend->d_spmv_buffer != NULL) {
        cudaFree(backend->d_spmv_buffer);
    }
    if (backend->d_state_in != NULL) {
        cudaFree(backend->d_state_in);
    }
    if (backend->d_state_out != NULL) {
        cudaFree(backend->d_state_out);
    }
    if (backend->d_active_in != NULL) {
        cudaFree(backend->d_active_in);
    }
    if (backend->d_active_out != NULL) {
        cudaFree(backend->d_active_out);
    }
    if (backend->d_x != NULL) {
        cudaFree(backend->d_x);
    }
    if (backend->d_y != NULL) {
        cudaFree(backend->d_y);
    }
    if (backend->cusparse_handle != NULL) {
        cusparseDestroy(backend->cusparse_handle);
    }
    if (backend->cublas_handle != NULL) {
        cublasDestroy(backend->cublas_handle);
    }
    free(backend);
    return 0;
}

extern "C" int rslmto_backend_upload_operator(
    void *backend_handle,
    const char *op_name,
    int block_dim,
    int n_sites,
    int nnzb,
    const int *row_ptr,
    const int *col_ind,
    const int *site_types,
    const double _Complex *blocks) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    plugin_operator *op = NULL;
    size_t row_ptr_len = 0U;
    size_t block_len = 0U;
    int scalar_dim = 0;
    int scalar_nnz = 0;
    int site = 0;
    int local_row = 0;
    int entry = 0;
    int cursor = 0;

    if (backend == NULL || op_name == NULL || row_ptr == NULL || col_ind == NULL || site_types == NULL || blocks == NULL) {
        return 1;
    }
    if (block_dim != backend->block_dim || n_sites != backend->n_sites || nnzb < 0) {
        return 1;
    }

    op = upsert_operator(backend, op_name);
    if (op == NULL || op->name == NULL) {
        return 1;
    }

    op->block_dim = block_dim;
    op->n_sites = n_sites;
    op->nnzb = nnzb;
    op->scalar_dim = block_dim * n_sites;
    op->scalar_nnz = nnzb * block_dim * block_dim;
    scalar_dim = op->scalar_dim;
    scalar_nnz = op->scalar_nnz;
    row_ptr_len = (size_t) (n_sites + 1);
    block_len = (size_t) block_dim * (size_t) block_dim * (size_t) nnzb;

    op->row_ptr = (int *) malloc(row_ptr_len * sizeof(*op->row_ptr));
    op->col_ind = (int *) malloc((size_t) nnzb * sizeof(*op->col_ind));
    op->site_types = (int *) malloc((size_t) n_sites * sizeof(*op->site_types));
    op->blocks = (cuDoubleComplex *) malloc(block_len * sizeof(*op->blocks));
    op->csr_row_ptr = (int *) malloc((size_t) (scalar_dim + 1) * sizeof(*op->csr_row_ptr));
    op->csr_col_ind = (int *) malloc((size_t) scalar_nnz * sizeof(*op->csr_col_ind));
    op->csr_values = (cuDoubleComplex *) malloc((size_t) scalar_nnz * sizeof(*op->csr_values));
    if (op->row_ptr == NULL || op->col_ind == NULL || op->site_types == NULL || op->blocks == NULL ||
        op->csr_row_ptr == NULL || op->csr_col_ind == NULL || op->csr_values == NULL) {
        free_operator(op);
        return 1;
    }

    memcpy(op->row_ptr, row_ptr, row_ptr_len * sizeof(*op->row_ptr));
    memcpy(op->col_ind, col_ind, (size_t) nnzb * sizeof(*op->col_ind));
    memcpy(op->site_types, site_types, (size_t) n_sites * sizeof(*op->site_types));
    memcpy(op->blocks, blocks, block_len * sizeof(*op->blocks));

    cursor = 0;
    op->csr_row_ptr[0] = 0;
    for (site = 0; site < n_sites; ++site) {
        int block_start = row_ptr[site] - 1;
        int block_stop = row_ptr[site + 1] - 1;
        int block_count = block_stop - block_start;
        if (block_start < 0 || block_stop < block_start || block_stop > nnzb) {
            free_operator(op);
            return 1;
        }
        for (local_row = 0; local_row < block_dim; ++local_row) {
            int scalar_row = site * block_dim + local_row;
            for (entry = block_start; entry < block_stop; ++entry) {
                int block_col = col_ind[entry] - 1;
                int local_col = 0;
                if (block_col < 0 || block_col >= n_sites) {
                    free_operator(op);
                    return 1;
                }
                for (local_col = 0; local_col < block_dim; ++local_col) {
                    size_t block_idx = (size_t) local_row +
                                       (size_t) local_col * (size_t) block_dim +
                                       (size_t) entry * (size_t) block_dim * (size_t) block_dim;
                    if (cursor >= scalar_nnz) {
                        free_operator(op);
                        return 1;
                    }
                    op->csr_col_ind[cursor] = block_col * block_dim + local_col;
                    op->csr_values[cursor] = op->blocks[block_idx];
                    ++cursor;
                }
            }
            op->csr_row_ptr[scalar_row + 1] = op->csr_row_ptr[scalar_row] + block_count * block_dim;
        }
    }

    if (cursor != scalar_nnz) {
        free_operator(op);
        return 1;
    }

    if (cuda_ok(cudaMalloc((void **) &op->d_csr_row_ptr, (size_t) (scalar_dim + 1) * sizeof(*op->d_csr_row_ptr))) != 0 ||
        cuda_ok(cudaMalloc((void **) &op->d_csr_col_ind, (size_t) scalar_nnz * sizeof(*op->d_csr_col_ind))) != 0 ||
        cuda_ok(cudaMalloc((void **) &op->d_csr_values, (size_t) scalar_nnz * sizeof(*op->d_csr_values))) != 0) {
        free_operator(op);
        return 1;
    }

    if (cuda_ok(cudaMemcpy(op->d_csr_row_ptr, op->csr_row_ptr, (size_t) (scalar_dim + 1) * sizeof(*op->d_csr_row_ptr), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(op->d_csr_col_ind, op->csr_col_ind, (size_t) scalar_nnz * sizeof(*op->d_csr_col_ind), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(op->d_csr_values, op->csr_values, (size_t) scalar_nnz * sizeof(*op->d_csr_values), cudaMemcpyHostToDevice)) != 0) {
        free_operator(op);
        return 1;
    }

    return 0;
}

extern "C" int rslmto_backend_apply_operator(
    void *backend_handle,
    const char *op_name,
    char trans_mode,
    int block_dim,
    int n_sites,
    const double _Complex *psi_in,
    const int *active_in,
    double _Complex *psi_out,
    int *active_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;

    if (backend == NULL || block_dim != backend->block_dim || n_sites != backend->n_sites) {
        return 1;
    }

    return apply_named_operator(
        backend,
        op_name,
        trans_mode,
        (const cuDoubleComplex *) psi_in,
        active_in,
        (cuDoubleComplex *) psi_out,
        active_out);
}

extern "C" int rslmto_backend_run_scalar_lanczos(
    void *backend_handle,
    int llmax,
    int n_targets,
    const int *target_sites,
    int hoh_enabled,
    double *a_out,
    double *b2_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psi = NULL;
    cuDoubleComplex *pmn = NULL;
    cuDoubleComplex *hpsi = NULL;
    int *active = NULL;
    int *active_out = NULL;
    int wf = 0;
    int orb = 0;
    int ll = 0;

    if (backend == NULL || llmax <= 0 || n_targets <= 0 || target_sites == NULL || a_out == NULL || b2_out == NULL) {
        return 1;
    }

    psi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi));
    pmn = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*pmn));
    hpsi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*hpsi));
    active = (int *) malloc((size_t) backend->n_sites * sizeof(*active));
    active_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_out));
    if (psi == NULL || pmn == NULL || hpsi == NULL || active == NULL || active_out == NULL) {
        free(psi);
        free(pmn);
        free(hpsi);
        free(active);
        free(active_out);
        return 1;
    }

    memset(a_out, 0, (size_t) llmax * (size_t) backend->block_dim * (size_t) n_targets * sizeof(*a_out));
    memset(b2_out, 0, (size_t) llmax * (size_t) backend->block_dim * (size_t) n_targets * sizeof(*b2_out));

    for (wf = 0; wf < n_targets; ++wf) {
        int site = target_sites[wf] - 1;
        if (site < 0 || site >= backend->n_sites) {
            free(psi);
            free(pmn);
            free(hpsi);
            free(active);
            free(active_out);
            return 1;
        }
        for (orb = 0; orb < backend->block_dim; ++orb) {
            double summ = 1.0;
            state_zero(backend, psi);
            state_zero(backend, pmn);
            mask_zero(backend, active);
            psi[site_block_offset(backend, site) + matrix_offset(backend->block_dim, orb, orb)] = c_one();
            active[site] = 1;
            for (ll = 0; ll < llmax - 1; ++ll) {
                cuDoubleComplex aj;
                size_t idx = 0U;
                if (apply_hamiltonian(backend, hoh_enabled, psi, active, 1.0, 0.0, hpsi, active_out) != 0) {
                    free(psi);
                    free(pmn);
                    free(hpsi);
                    free(active);
                    free(active_out);
                    return 1;
                }
                for (idx = 0U; idx < state_len_of(backend); ++idx) {
                    pmn[idx] = c_sub(hpsi[idx], pmn[idx]);
                }
                aj = state_inner_product(backend, psi, hpsi);
                a_out[scalar_coeff_offset(llmax, backend->block_dim, ll, orb, wf)] = cuCreal(aj);
                b2_out[scalar_coeff_offset(llmax, backend->block_dim, ll, orb, wf)] = summ;
                for (idx = 0U; idx < state_len_of(backend); ++idx) {
                    pmn[idx] = c_sub(pmn[idx], c_scale(psi[idx], cuCreal(aj)));
                }
                summ = state_norm2(backend, pmn);
                if (summ > 0.0) {
                    double inv_norm = 1.0 / sqrt(summ);
                    for (idx = 0U; idx < state_len_of(backend); ++idx) {
                        cuDoubleComplex prev = psi[idx];
                        psi[idx] = c_scale(pmn[idx], inv_norm);
                        pmn[idx] = c_scale(prev, -sqrt(summ));
                    }
                    active_copy(backend, active, active_out);
                } else {
                    state_zero(backend, psi);
                    state_zero(backend, pmn);
                    mask_zero(backend, active);
                }
            }
            b2_out[scalar_coeff_offset(llmax, backend->block_dim, llmax - 1, orb, wf)] = summ;
        }
    }

    free(psi);
    free(pmn);
    free(hpsi);
    free(active);
    free(active_out);
    return 0;
}

extern "C" int rslmto_backend_run_block_lanczos(
    void *backend_handle,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    double _Complex *a_b_out,
    double _Complex *b2_b_out,
    double *a_diag_out,
    double *b2_diag_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psi = NULL;
    cuDoubleComplex *pmn = NULL;
    cuDoubleComplex *hpsi = NULL;
    cuDoubleComplex *psi_t = NULL;
    cuDoubleComplex *sum_b = NULL;
    cuDoubleComplex *sum_a = NULL;
    cuDoubleComplex *b = NULL;
    cuDoubleComplex *b_inv = NULL;
    cuDoubleComplex *tmp = NULL;
    int *active = NULL;
    int *active_out = NULL;
    int wf = 0;
    int ll = 0;
    int site = 0;
    int i = 0;
    int j = 0;

    if (backend == NULL || llmax <= 0 || n_workflows <= 0) {
        return 1;
    }

    psi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi));
    pmn = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*pmn));
    hpsi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*hpsi));
    psi_t = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi_t));
    sum_b = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*sum_b));
    sum_a = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*sum_a));
    b = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*b));
    b_inv = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*b_inv));
    tmp = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
    active = (int *) malloc((size_t) backend->n_sites * sizeof(*active));
    active_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_out));
    if (psi == NULL || pmn == NULL || hpsi == NULL || psi_t == NULL || sum_b == NULL || sum_a == NULL ||
        b == NULL || b_inv == NULL || tmp == NULL || active == NULL || active_out == NULL) {
        free(psi);
        free(pmn);
        free(hpsi);
        free(psi_t);
        free(sum_b);
        free(sum_a);
        free(b);
        free(b_inv);
        free(tmp);
        free(active);
        free(active_out);
        return 1;
    }

    memset(a_b_out, 0, (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) llmax * (size_t) n_workflows * sizeof(*a_b_out));
    memset(b2_b_out, 0, (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) llmax * (size_t) n_workflows * sizeof(*b2_b_out));
    memset(a_diag_out, 0, (size_t) llmax * (size_t) backend->block_dim * (size_t) n_workflows * sizeof(*a_diag_out));
    memset(b2_diag_out, 0, (size_t) llmax * (size_t) backend->block_dim * (size_t) n_workflows * sizeof(*b2_diag_out));

    for (wf = 0; wf < n_workflows; ++wf) {
        if (site_i[wf] == site_j[wf] && variant_id[wf] != 1) {
            continue;
        }

        set_workflow_state(backend, site_i[wf], site_j[wf], variant_id[wf], psi, active);
        state_zero(backend, pmn);
        matrix_identity(sum_b, backend->block_dim);

        for (ll = 0; ll < llmax - 1; ++ll) {
            if (apply_hamiltonian(backend, hoh_enabled, psi, active, 1.0, 0.0, hpsi, active_out) != 0) {
                free(psi);
                free(pmn);
                free(hpsi);
                free(psi_t);
                free(sum_b);
                free(sum_a);
                free(b);
                free(b_inv);
                free(tmp);
                free(active);
                free(active_out);
                return 1;
            }

            matrix_zero(sum_a, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                size_t off = 0U;
                int idx = 0;
                if (active_out[site] == 0) {
                    continue;
                }
                off = site_block_offset(backend, site);
                for (idx = 0; idx < backend->block_dim * backend->block_dim; ++idx) {
                    pmn[off + (size_t) idx] = c_sub(hpsi[off + (size_t) idx], pmn[off + (size_t) idx]);
                }
                accumulate_left_dagger_right(sum_a, psi + off, hpsi + off, backend->block_dim);
            }

            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    ((cuDoubleComplex *) a_b_out)[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, ll, wf)] =
                        sum_a[matrix_offset(backend->block_dim, i, j)];
                }
                a_diag_out[scalar_coeff_offset(llmax, backend->block_dim, ll, j, wf)] =
                    cuCreal(sum_a[matrix_offset(backend->block_dim, j, j)]);
            }

            state_copy(backend, psi_t, psi);
            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    ((cuDoubleComplex *) b2_b_out)[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, ll, wf)] =
                        sum_b[matrix_offset(backend->block_dim, i, j)];
                }
                b2_diag_out[scalar_coeff_offset(llmax, backend->block_dim, ll, j, wf)] =
                    cuCreal(sum_b[matrix_offset(backend->block_dim, j, j)]);
            }

            matrix_zero(sum_b, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                size_t off = 0U;
                int idx = 0;
                if (active_out[site] == 0) {
                    continue;
                }
                off = site_block_offset(backend, site);
                matrix_right_multiply(tmp, psi + off, sum_a, backend->block_dim);
                for (idx = 0; idx < backend->block_dim * backend->block_dim; ++idx) {
                    pmn[off + (size_t) idx] = c_sub(pmn[off + (size_t) idx], tmp[idx]);
                }
                accumulate_left_dagger_right(sum_b, pmn + off, pmn + off, backend->block_dim);
            }

            if (hermitian_sqrt_and_inverse(backend->block_dim, sum_b, b, b_inv) != 0) {
                free(psi);
                free(pmn);
                free(hpsi);
                free(psi_t);
                free(sum_b);
                free(sum_a);
                free(b);
                free(b_inv);
                free(tmp);
                free(active);
                free(active_out);
                return 1;
            }

            for (site = 0; site < backend->n_sites; ++site) {
                size_t off = 0U;
                if (active_out[site] == 0) {
                    continue;
                }
                off = site_block_offset(backend, site);
                matrix_right_multiply(tmp, pmn + off, b_inv, backend->block_dim);
                memcpy(psi + off, tmp, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
                matrix_right_multiply(tmp, psi_t + off, b, backend->block_dim);
                memcpy(pmn + off, tmp, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
            }

            active_copy(backend, active, active_out);
        }

        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                ((cuDoubleComplex *) b2_b_out)[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, llmax - 1, wf)] =
                    sum_b[matrix_offset(backend->block_dim, i, j)];
            }
            b2_diag_out[scalar_coeff_offset(llmax, backend->block_dim, llmax - 1, j, wf)] =
                cuCreal(sum_b[matrix_offset(backend->block_dim, j, j)]);
        }
    }

    free(psi);
    free(pmn);
    free(hpsi);
    free(psi_t);
    free(sum_b);
    free(sum_a);
    free(b);
    free(b_inv);
    free(tmp);
    free(active);
    free(active_out);
    return 0;
}

extern "C" int rslmto_backend_run_block_chebyshev(
    void *backend_handle,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    double a_scale,
    double b_shift,
    double _Complex *mu_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psi0 = NULL;
    cuDoubleComplex *psi1 = NULL;
    cuDoubleComplex *psi2 = NULL;
    cuDoubleComplex *psiref = NULL;
    cuDoubleComplex *dum1 = NULL;
    cuDoubleComplex *dum2 = NULL;
    int *active0 = NULL;
    int *active = NULL;
    int *active_out = NULL;
    int wf = 0;
    int ll = 0;
    int site = 0;
    int i = 0;
    int j = 0;
    int nmom = 2 * llmax + 2;

    if (backend == NULL || llmax <= 0 || n_workflows <= 0) {
        return 1;
    }

    psi0 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi0));
    psi1 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi1));
    psi2 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi2));
    psiref = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psiref));
    dum1 = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*dum1));
    dum2 = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*dum2));
    active0 = (int *) malloc((size_t) backend->n_sites * sizeof(*active0));
    active = (int *) malloc((size_t) backend->n_sites * sizeof(*active));
    active_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_out));
    if (psi0 == NULL || psi1 == NULL || psi2 == NULL || psiref == NULL || dum1 == NULL || dum2 == NULL ||
        active0 == NULL || active == NULL || active_out == NULL) {
        free(psi0);
        free(psi1);
        free(psi2);
        free(psiref);
        free(dum1);
        free(dum2);
        free(active0);
        free(active);
        free(active_out);
        return 1;
    }

    memset(mu_out, 0, (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) nmom * (size_t) n_workflows * sizeof(*mu_out));

    for (wf = 0; wf < n_workflows; ++wf) {
        if (site_i[wf] == site_j[wf] && variant_id[wf] != 1) {
            continue;
        }

        set_workflow_state(backend, site_i[wf], site_j[wf], variant_id[wf], psi0, active0);
        state_copy(backend, psiref, psi0);
        mask_zero(backend, active);
        mask_zero(backend, active_out);

        matrix_zero(dum1, backend->block_dim);
        for (site = 0; site < backend->n_sites; ++site) {
            if (active0[site] == 0) {
                continue;
            }
            accumulate_left_dagger_right(dum1, psiref + site_block_offset(backend, site), psi0 + site_block_offset(backend, site), backend->block_dim);
        }
        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 0, wf)] = dum1[matrix_offset(backend->block_dim, i, j)];
            }
        }

        if (apply_hamiltonian(backend, hoh_enabled, psi0, active0, a_scale, b_shift, psi1, active_out) != 0) {
            free(psi0);
            free(psi1);
            free(psi2);
            free(psiref);
            free(dum1);
            free(dum2);
            free(active0);
            free(active);
            free(active_out);
            return 1;
        }

        matrix_zero(dum1, backend->block_dim);
        for (site = 0; site < backend->n_sites; ++site) {
            if (active0[site] == 0) {
                continue;
            }
            accumulate_left_dagger_right(dum1, psiref + site_block_offset(backend, site), psi1 + site_block_offset(backend, site), backend->block_dim);
        }
        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 1, wf)] = dum1[matrix_offset(backend->block_dim, i, j)];
            }
        }

        active_copy(backend, active, active_out);
        for (ll = 0; ll < llmax; ++ll) {
            size_t idx = 0U;
            if (apply_hamiltonian(backend, hoh_enabled, psi1, active, a_scale, b_shift, psi2, active_out) != 0) {
                free(psi0);
                free(psi1);
                free(psi2);
                free(psiref);
                free(dum1);
                free(dum2);
                free(active0);
                free(active);
                free(active_out);
                return 1;
            }
            for (idx = 0U; idx < state_len_of(backend); ++idx) {
                psi2[idx] = c_scale(psi2[idx], 2.0);
            }

            matrix_zero(dum1, backend->block_dim);
            matrix_zero(dum2, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                size_t off = 0U;
                int idm = 0;
                if (active_out[site] == 0) {
                    continue;
                }
                off = site_block_offset(backend, site);
                for (idm = 0; idm < backend->block_dim * backend->block_dim; ++idm) {
                    psi2[off + (size_t) idm] = c_sub(psi2[off + (size_t) idm], psi0[off + (size_t) idm]);
                }
                accumulate_left_dagger_right(dum1, psi1 + off, psi1 + off, backend->block_dim);
                accumulate_left_dagger_right(dum2, psi2 + off, psi1 + off, backend->block_dim);
                memcpy(psi0 + off, psi1 + off, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*psi0));
                memcpy(psi1 + off, psi2 + off, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*psi1));
            }

            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 2 * ll + 2, wf)] =
                        c_sub(c_scale(dum1[matrix_offset(backend->block_dim, i, j)], 2.0),
                              ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 0, wf)]);
                    ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 2 * ll + 3, wf)] =
                        c_sub(c_scale(dum2[matrix_offset(backend->block_dim, i, j)], 2.0),
                              ((cuDoubleComplex *) mu_out)[mu_offset(backend->block_dim, nmom, i, j, 1, wf)]);
                }
            }

            active_copy(backend, active, active_out);
        }
    }

    free(psi0);
    free(psi1);
    free(psi2);
    free(psiref);
    free(dum1);
    free(dum2);
    free(active0);
    free(active);
    free(active_out);
    return 0;
}

extern "C" int rslmto_backend_run_transport_stochastic(
    void *backend_handle,
    int llmax,
    int loop_over,
    int calc_mode,
    int hoh_enabled,
    const int *type_sites,
    double a_scale,
    double b_shift,
    double _Complex *mu_nm_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psiref = NULL;
    cuDoubleComplex *w0 = NULL;
    cuDoubleComplex *w1 = NULL;
    cuDoubleComplex *w2 = NULL;
    cuDoubleComplex *right_vec = NULL;
    cuDoubleComplex *v0 = NULL;
    cuDoubleComplex *v1 = NULL;
    cuDoubleComplex *v2 = NULL;
    cuDoubleComplex *left_vec = NULL;
    cuDoubleComplex *dum = NULL;
    int *active_ref = NULL;
    int *active_w = NULL;
    int *active_w_out = NULL;
    int *active_v = NULL;
    int *active_v_out = NULL;
    int loop_idx = 0;
    int m = 0;
    int n = 0;
    int i = 0;
    int j = 0;
    int site = 0;
    const double twopi = 6.28318530717958647692;

    if (backend == NULL || llmax <= 0 || loop_over <= 0 || mu_nm_out == NULL) {
        return 1;
    }

    psiref = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psiref));
    w0 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*w0));
    w1 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*w1));
    w2 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*w2));
    right_vec = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*right_vec));
    v0 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*v0));
    v1 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*v1));
    v2 = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*v2));
    left_vec = (cuDoubleComplex *) malloc((size_t) llmax * state_len_of(backend) * sizeof(*left_vec));
    dum = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*dum));
    active_ref = (int *) malloc((size_t) backend->n_sites * sizeof(*active_ref));
    active_w = (int *) malloc((size_t) backend->n_sites * sizeof(*active_w));
    active_w_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_w_out));
    active_v = (int *) malloc((size_t) backend->n_sites * sizeof(*active_v));
    active_v_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_v_out));
    if (psiref == NULL || w0 == NULL || w1 == NULL || w2 == NULL || right_vec == NULL || v0 == NULL || v1 == NULL || v2 == NULL ||
        left_vec == NULL || dum == NULL || active_ref == NULL || active_w == NULL || active_w_out == NULL || active_v == NULL || active_v_out == NULL) {
        free(psiref);
        free(w0);
        free(w1);
        free(w2);
        free(right_vec);
        free(v0);
        free(v1);
        free(v2);
        free(left_vec);
        free(dum);
        free(active_ref);
        free(active_w);
        free(active_w_out);
        free(active_v);
        free(active_v_out);
        return 1;
    }

    memset(mu_nm_out, 0, (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) llmax * (size_t) llmax * (size_t) loop_over * sizeof(*mu_nm_out));
    srand(12345U);

    for (loop_idx = 0; loop_idx < loop_over; ++loop_idx) {
        state_zero(backend, psiref);
        state_zero(backend, w0);
        state_zero(backend, w1);
        state_zero(backend, w2);
        state_zero(backend, right_vec);
        state_zero(backend, v0);
        state_zero(backend, v1);
        state_zero(backend, v2);
        memset(left_vec, 0, (size_t) llmax * state_len_of(backend) * sizeof(*left_vec));
        mask_zero(backend, active_ref);
        mask_zero(backend, active_w);
        mask_zero(backend, active_w_out);
        mask_zero(backend, active_v);
        mask_zero(backend, active_v_out);

        if (calc_mode == 1) {
            int target_site = type_sites != NULL ? type_sites[loop_idx] - 1 : -1;
            if (target_site < 0 || target_site >= backend->n_sites) {
                free(psiref);
                free(w0);
                free(w1);
                free(w2);
                free(right_vec);
                free(v0);
                free(v1);
                free(v2);
                free(left_vec);
                free(dum);
                free(active_ref);
                free(active_w);
                free(active_w_out);
                free(active_v);
                free(active_v_out);
                return 1;
            }
            active_ref[target_site] = 1;
            for (i = 0; i < backend->block_dim; ++i) {
                psiref[site_block_offset(backend, target_site) + matrix_offset(backend->block_dim, i, i)] = c_one();
            }
        } else {
            for (site = 0; site < backend->n_sites; ++site) {
                double rng = ((double) rand()) / ((double) RAND_MAX);
                cuDoubleComplex phase = make_cuDoubleComplex(cos(twopi * rng), sin(twopi * rng));
                active_ref[site] = 1;
                for (i = 0; i < backend->block_dim; ++i) {
                    psiref[site_block_offset(backend, site) + matrix_offset(backend->block_dim, i, i)] = phase;
                }
            }
            {
                double norm = sqrt((double) backend->n_sites);
                size_t idx = 0U;
                for (idx = 0U; idx < state_len_of(backend); ++idx) {
                    psiref[idx] = c_scale(psiref[idx], 1.0 / norm);
                }
            }
        }

        for (m = 0; m < llmax; ++m) {
            if (m == 0) {
                state_copy(backend, w1, psiref);
                active_copy(backend, active_w, active_ref);
            } else if (m == 1) {
                state_copy(backend, w0, w1);
                if (apply_hamiltonian(backend, hoh_enabled, w0, active_w, a_scale, b_shift, w1, active_w_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
                active_copy(backend, active_w, active_w_out);
            } else {
                size_t idx = 0U;
                if (apply_hamiltonian(backend, hoh_enabled, w1, active_w, a_scale, b_shift, w2, active_w_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
                for (idx = 0U; idx < state_len_of(backend); ++idx) {
                    w2[idx] = c_sub(c_scale(w2[idx], 2.0), w0[idx]);
                }
                state_copy(backend, w0, w1);
                state_copy(backend, w1, w2);
                state_zero(backend, w2);
                active_copy(backend, active_w, active_w_out);
            }
            memcpy(left_vec + moment_state_offset(state_len_of(backend), m), w1, state_len_of(backend) * sizeof(*left_vec));
        }

        if (hoh_enabled != 0) {
            if (apply_velocity_hoh_named(backend, "transport_velocity_b", "transport_velocity_overlap_b", psiref, active_ref, v0, active_v) != 0) {
                free(psiref);
                free(w0);
                free(w1);
                free(w2);
                free(right_vec);
                free(v0);
                free(v1);
                free(v2);
                free(left_vec);
                free(dum);
                free(active_ref);
                free(active_w);
                free(active_w_out);
                free(active_v);
                free(active_v_out);
                return 1;
            }
        } else {
            if (apply_named_operator(backend, "transport_velocity_b", 'n', psiref, active_ref, v0, active_v) != 0) {
                free(psiref);
                free(w0);
                free(w1);
                free(w2);
                free(right_vec);
                free(v0);
                free(v1);
                free(v2);
                free(left_vec);
                free(dum);
                free(active_ref);
                free(active_w);
                free(active_w_out);
                free(active_v);
                free(active_v_out);
                return 1;
            }
        }

        for (n = 0; n < llmax; ++n) {
            if (n == 0) {
                state_copy(backend, v1, v0);
            } else if (n == 1) {
                state_copy(backend, v0, v1);
                if (apply_hamiltonian(backend, hoh_enabled, v0, active_v, a_scale, b_shift, v1, active_v_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
                active_copy(backend, active_v, active_v_out);
            } else {
                size_t idx = 0U;
                if (apply_hamiltonian(backend, hoh_enabled, v1, active_v, a_scale, b_shift, v2, active_v_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
                for (idx = 0U; idx < state_len_of(backend); ++idx) {
                    v2[idx] = c_sub(c_scale(v2[idx], 2.0), v0[idx]);
                }
                state_copy(backend, v0, v1);
                state_copy(backend, v1, v2);
                state_zero(backend, v2);
                active_copy(backend, active_v, active_v_out);
            }

            if (hoh_enabled != 0) {
                if (apply_velocity_hoh_named(backend, "transport_velocity_a", "transport_velocity_overlap_a", v1, active_v, right_vec, active_v_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
            } else {
                if (apply_named_operator(backend, "transport_velocity_a", 'n', v1, active_v, right_vec, active_v_out) != 0) {
                    free(psiref);
                    free(w0);
                    free(w1);
                    free(w2);
                    free(right_vec);
                    free(v0);
                    free(v1);
                    free(v2);
                    free(left_vec);
                    free(dum);
                    free(active_ref);
                    free(active_w);
                    free(active_w_out);
                    free(active_v);
                    free(active_v_out);
                    return 1;
                }
            }

            for (m = 0; m < llmax; ++m) {
                matrix_zero(dum, backend->block_dim);
                for (site = 0; site < backend->n_sites; ++site) {
                    accumulate_left_dagger_right(
                        dum,
                        left_vec + moment_state_offset(state_len_of(backend), m) + site_block_offset(backend, site),
                        right_vec + site_block_offset(backend, site),
                        backend->block_dim);
                }
                for (j = 0; j < backend->block_dim; ++j) {
                    for (i = 0; i < backend->block_dim; ++i) {
                        ((cuDoubleComplex *) mu_nm_out)[mu_nm_offset(backend->block_dim, llmax, i, j, n, m, loop_idx)] =
                            dum[matrix_offset(backend->block_dim, i, j)];
                    }
                }
            }
        }
    }

    free(psiref);
    free(w0);
    free(w1);
    free(w2);
    free(right_vec);
    free(v0);
    free(v1);
    free(v2);
    free(left_vec);
    free(dum);
    free(active_ref);
    free(active_w);
    free(active_w_out);
    free(active_v);
    free(active_v_out);
    return 0;
}

extern "C" const char *rslmto_backend_name(void) {
    return "rslmto_cuda_sparse_plugin";
}

extern "C" const char *rslmto_backend_capabilities(void) {
    return "execution=nvidia_gpu;library=cuda;precision=complex_fp64;operators=block_sparse,velocity,hoh;recurrences=scalar_lanczos,block_lanczos,block_chebyshev,transport;storage=scalar_csr_lowering;spmv=cusparse_spmv;dense=cublas_ready";
}
