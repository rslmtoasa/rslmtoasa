#include <cmath>
#include <cstdlib>
#include <cstring>

#include <cuda_runtime.h>
#include <cuComplex.h>

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
    int *d_row_ptr;
    int *d_col_ind;
    int *d_site_types;
    cuDoubleComplex *d_blocks;
} plugin_operator;

typedef struct {
    char *precision_mode;
    char *library_name;
    int block_dim;
    int n_sites;
    int capacity;
    int count;
    plugin_operator *ops;
    cuDoubleComplex *d_psi_in;
    cuDoubleComplex *d_psi_out;
    int *d_active_in;
    int *d_active_out;
    size_t state_capacity;
    size_t site_capacity;
} plugin_backend;

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
    size_t len = 0;
    char *copy = 0;
    if (input == 0) {
        return 0;
    }
    len = strlen(input);
    copy = (char *) malloc(len + 1U);
    if (copy == 0) {
        return 0;
    }
    memcpy(copy, input, len + 1U);
    return copy;
}

static int cuda_ok(cudaError_t status) {
    return status == cudaSuccess ? 0 : 1;
}

static void free_operator(plugin_operator *op) {
    if (op == 0) {
        return;
    }
    free(op->name);
    free(op->row_ptr);
    free(op->col_ind);
    free(op->site_types);
    free(op->blocks);
    if (op->d_row_ptr != 0) cudaFree(op->d_row_ptr);
    if (op->d_col_ind != 0) cudaFree(op->d_col_ind);
    if (op->d_site_types != 0) cudaFree(op->d_site_types);
    if (op->d_blocks != 0) cudaFree(op->d_blocks);
    memset(op, 0, sizeof(*op));
}

static plugin_operator *find_operator(plugin_backend *backend, const char *name) {
    int i = 0;
    if (backend == 0 || name == 0) {
        return 0;
    }
    for (i = 0; i < backend->count; ++i) {
        if (backend->ops[i].name != 0 && strcmp(backend->ops[i].name, name) == 0) {
            return &backend->ops[i];
        }
    }
    return 0;
}

static plugin_operator *upsert_operator(plugin_backend *backend, const char *name) {
    plugin_operator *op = 0;
    plugin_operator *grown = 0;

    op = find_operator(backend, name);
    if (op != 0) {
        free_operator(op);
        op->name = dup_string(name);
        return op;
    }

    if (backend->count == backend->capacity) {
        int next_capacity = backend->capacity == 0 ? 8 : 2 * backend->capacity;
        grown = (plugin_operator *) realloc(backend->ops, (size_t) next_capacity * sizeof(*grown));
        if (grown == 0) {
            return 0;
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

static int ensure_state_buffers(plugin_backend *backend, int block_dim, int n_sites) {
    size_t state_len = (size_t) block_dim * (size_t) block_dim * (size_t) n_sites;
    size_t site_len = (size_t) n_sites;

    if (backend == 0) {
        return 1;
    }
    if (backend->state_capacity == state_len && backend->site_capacity == site_len &&
        backend->d_psi_in != 0 && backend->d_psi_out != 0 &&
        backend->d_active_in != 0 && backend->d_active_out != 0) {
        return 0;
    }

    if (backend->d_psi_in != 0) cudaFree(backend->d_psi_in);
    if (backend->d_psi_out != 0) cudaFree(backend->d_psi_out);
    if (backend->d_active_in != 0) cudaFree(backend->d_active_in);
    if (backend->d_active_out != 0) cudaFree(backend->d_active_out);
    backend->d_psi_in = 0;
    backend->d_psi_out = 0;
    backend->d_active_in = 0;
    backend->d_active_out = 0;
    backend->state_capacity = 0U;
    backend->site_capacity = 0U;

    if (cuda_ok(cudaMalloc((void **) &backend->d_psi_in, state_len * sizeof(*backend->d_psi_in))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_psi_out, state_len * sizeof(*backend->d_psi_out))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_active_in, site_len * sizeof(*backend->d_active_in))) != 0 ||
        cuda_ok(cudaMalloc((void **) &backend->d_active_out, site_len * sizeof(*backend->d_active_out))) != 0) {
        return 1;
    }

    backend->state_capacity = state_len;
    backend->site_capacity = site_len;
    return 0;
}

static size_t state_len_of(const plugin_backend *backend) {
    return (size_t) backend->block_dim * (size_t) backend->block_dim * (size_t) backend->n_sites;
}

static size_t site_block_offset(const plugin_backend *backend, int site) {
    return (size_t) site * (size_t) backend->block_dim * (size_t) backend->block_dim;
}

static size_t matrix_offset(int nrow, int ncol, int i, int j) {
    return (size_t) i + (size_t) nrow * (size_t) j;
}

static size_t a_b_offset(int nb, int llmax, int n_workflow, int i, int j, int ll, int wf) {
    return (size_t) i + (size_t) nb * (size_t) j + (size_t) nb * (size_t) nb * (size_t) ll +
           (size_t) nb * (size_t) nb * (size_t) llmax * (size_t) wf;
}

static size_t mu_offset(int nb, int nmom, int i, int j, int n, int wf) {
    return (size_t) i + (size_t) nb * (size_t) j + (size_t) nb * (size_t) nb * (size_t) n +
           (size_t) nb * (size_t) nb * (size_t) nmom * (size_t) wf;
}

static void state_zero(plugin_backend *backend, cuDoubleComplex *state) {
    memset(state, 0, state_len_of(backend) * sizeof(*state));
}

static void mask_zero(plugin_backend *backend, int *mask) {
    memset(mask, 0, (size_t) backend->n_sites * sizeof(*mask));
}

static void matrix_zero(cuDoubleComplex *mat, int n) {
    memset(mat, 0, (size_t) n * (size_t) n * sizeof(*mat));
}

static void matrix_identity(cuDoubleComplex *mat, int n) {
    int i = 0;
    matrix_zero(mat, n);
    for (i = 0; i < n; ++i) {
        mat[matrix_offset(n, n, i, i)] = c_one();
    }
}

static void matrix_copy(cuDoubleComplex *dst, const cuDoubleComplex *src, int n) {
    memcpy(dst, src, (size_t) n * (size_t) n * sizeof(*dst));
}

static void state_copy(plugin_backend *backend, cuDoubleComplex *dst, const cuDoubleComplex *src) {
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
                accum = cuCadd(accum, cuCmul(left[matrix_offset(n, n, i, k)], right[matrix_offset(n, n, k, j)]));
            }
            out[matrix_offset(n, n, i, j)] = accum;
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
                entry = cuCadd(entry, cuCmul(cuConj(left[matrix_offset(n, n, k, i)]), right[matrix_offset(n, n, k, j)]));
            }
            accum[matrix_offset(n, n, i, j)] = cuCadd(accum[matrix_offset(n, n, i, j)], entry);
        }
    }
}

static void scale_shift_state(plugin_backend *backend, cuDoubleComplex *out, const cuDoubleComplex *in, double a_scale, double b_shift) {
    size_t idx = 0;
    size_t total = state_len_of(backend);
    for (idx = 0; idx < total; ++idx) {
        out[idx] = c_scale(c_sub(out[idx], c_scale(in[idx], b_shift)), 1.0 / a_scale);
    }
}

static int apply_named_operator(plugin_backend *backend, const char *name, char trans_mode, const cuDoubleComplex *psi_in, const int *active_in, cuDoubleComplex *psi_out, int *active_out);
__global__ static void apply_block_operator_kernel(
    int block_dim,
    int n_sites,
    const int *row_ptr,
    const int *col_ind,
    const cuDoubleComplex *blocks,
    char trans_mode,
    const int *active_in,
    const cuDoubleComplex *psi_in,
    cuDoubleComplex *psi_out,
    int *active_out);

static int apply_hamiltonian(plugin_backend *backend, int hoh_enabled, const cuDoubleComplex *psi_in, const int *active_in, double a_scale, double b_shift, cuDoubleComplex *psi_out, int *active_out) {
    const char *name = hoh_enabled != 0 ? "hamiltonian_hoh" : "hamiltonian";
    if (apply_named_operator(backend, name, 'n', psi_in, active_in, psi_out, active_out) != 0) {
        return 1;
    }
    scale_shift_state(backend, psi_out, psi_in, a_scale, b_shift);
    return 0;
}

static int apply_named_operator(plugin_backend *backend, const char *name, char trans_mode, const cuDoubleComplex *psi_in, const int *active_in, cuDoubleComplex *psi_out, int *active_out) {
    plugin_operator *op = 0;
    size_t state_len = state_len_of(backend);
    int threads_per_block = 256;
    int blocks_per_grid = 0;

    if (backend == 0 || psi_in == 0 || active_in == 0 || psi_out == 0 || active_out == 0) {
        return 1;
    }
    if (ensure_state_buffers(backend, backend->block_dim, backend->n_sites) != 0) {
        return 1;
    }

    op = find_operator(backend, name);
    if (op == 0 || op->d_row_ptr == 0 || op->d_col_ind == 0 || op->d_blocks == 0) {
        return 1;
    }

    if (cuda_ok(cudaMemcpy(backend->d_psi_in, psi_in, state_len * sizeof(*backend->d_psi_in), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(backend->d_active_in, active_in, (size_t) backend->n_sites * sizeof(*backend->d_active_in), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemset(backend->d_psi_out, 0, state_len * sizeof(*backend->d_psi_out))) != 0 ||
        cuda_ok(cudaMemset(backend->d_active_out, 0, (size_t) backend->n_sites * sizeof(*backend->d_active_out))) != 0) {
        return 1;
    }

    blocks_per_grid = (int) ((state_len + (size_t) threads_per_block - 1U) / (size_t) threads_per_block);
    apply_block_operator_kernel<<<blocks_per_grid, threads_per_block>>>(
        backend->block_dim,
        backend->n_sites,
        op->d_row_ptr,
        op->d_col_ind,
        op->d_blocks,
        trans_mode,
        backend->d_active_in,
        backend->d_psi_in,
        backend->d_psi_out,
        backend->d_active_out);
    if (cuda_ok(cudaGetLastError()) != 0 || cuda_ok(cudaDeviceSynchronize()) != 0) {
        return 1;
    }

    if (cuda_ok(cudaMemcpy(psi_out, backend->d_psi_out, state_len * sizeof(*backend->d_psi_out), cudaMemcpyDeviceToHost)) != 0 ||
        cuda_ok(cudaMemcpy(active_out, backend->d_active_out, (size_t) backend->n_sites * sizeof(*backend->d_active_out), cudaMemcpyDeviceToHost)) != 0) {
        return 1;
    }

    return 0;
}

static int hermitian_sqrt_and_inverse(int n, const cuDoubleComplex *b2, cuDoubleComplex *b, cuDoubleComplex *b_inv) {
    cuDoubleComplex *u = 0;
    cuDoubleComplex *work = 0;
    cuDoubleComplex work_query;
    double *evals = 0;
    double *rwork = 0;
    cuDoubleComplex *tmp = 0;
    cuDoubleComplex *tmp2 = 0;
    cuDoubleComplex *lam = 0;
    cuDoubleComplex *lam_inv = 0;
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
    if (u == 0 || evals == 0 || rwork == 0 || tmp == 0 || tmp2 == 0 || lam == 0 || lam_inv == 0) {
        free(u); free(evals); free(rwork); free(tmp); free(tmp2); free(lam); free(lam_inv);
        return 1;
    }

    matrix_copy(u, b2, n);
    zheev_(&jobz, &uplo, &n, u, &n, evals, &work_query, &lwork, rwork, &info);
    if (info != 0) {
        free(u); free(evals); free(rwork); free(tmp); free(tmp2); free(lam); free(lam_inv);
        return 1;
    }

    lwork = (int) cuCreal(work_query);
    if (lwork < 2 * n) lwork = 2 * n;
    work = (cuDoubleComplex *) malloc((size_t) lwork * sizeof(*work));
    if (work == 0) {
        free(u); free(evals); free(rwork); free(tmp); free(tmp2); free(lam); free(lam_inv);
        return 1;
    }

    matrix_copy(u, b2, n);
    zheev_(&jobz, &uplo, &n, u, &n, evals, work, &lwork, rwork, &info);
    if (info != 0) {
        free(u); free(evals); free(rwork); free(tmp); free(tmp2); free(lam); free(lam_inv); free(work);
        return 1;
    }

    matrix_zero(lam, n);
    matrix_zero(lam_inv, n);
    for (i = 0; i < n; ++i) {
      double ev = evals[i] > 0.0 ? evals[i] : 0.0;
      double root = sqrt(ev);
      lam[matrix_offset(n, n, i, i)] = c_real(root);
      lam_inv[matrix_offset(n, n, i, i)] = root > 0.0 ? c_real(1.0 / root) : c_zero();
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
                    accum = cuCadd(accum, cuCmul(tmp[matrix_offset(n, n, row, k)], cuConj(u[matrix_offset(n, n, col, k)])));
                }
                b[matrix_offset(n, n, row, col)] = accum;
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
                    accum = cuCadd(accum, cuCmul(tmp2[matrix_offset(n, n, row, k)], cuConj(u[matrix_offset(n, n, col, k)])));
                }
                b_inv[matrix_offset(n, n, row, col)] = accum;
            }
        }
    }

    free(u); free(evals); free(rwork); free(tmp); free(tmp2); free(lam); free(lam_inv); free(work);
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
        state[site_block_offset(backend, site_i - 1) + matrix_offset(backend->block_dim, backend->block_dim, l, l)] = a_sign;
        state[site_block_offset(backend, site_j - 1) + matrix_offset(backend->block_dim, backend->block_dim, l, l)] = b_sign;
    }
}

extern "C" void *rslmto_backend_create(const char *precision_mode, int block_dim, int n_sites, const char *library_name) {
    plugin_backend *backend = 0;
    if (precision_mode == 0 || strcmp(precision_mode, "complex_fp64") != 0) {
        return 0;
    }
    backend = (plugin_backend *) calloc(1U, sizeof(*backend));
    if (backend == 0) {
        return 0;
    }
    backend->precision_mode = dup_string(precision_mode);
    backend->library_name = dup_string(library_name != 0 ? library_name : "cuda");
    backend->block_dim = block_dim;
    backend->n_sites = n_sites;
    if (ensure_state_buffers(backend, block_dim, n_sites) != 0) {
        free(backend->precision_mode);
        free(backend->library_name);
        free(backend);
        return 0;
    }
    return backend;
}

extern "C" int rslmto_backend_destroy(void *backend_handle) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    int i = 0;
    if (backend == 0) {
        return 0;
    }
    for (i = 0; i < backend->count; ++i) {
        free_operator(&backend->ops[i]);
    }
    free(backend->ops);
    free(backend->precision_mode);
    free(backend->library_name);
    if (backend->d_psi_in != 0) cudaFree(backend->d_psi_in);
    if (backend->d_psi_out != 0) cudaFree(backend->d_psi_out);
    if (backend->d_active_in != 0) cudaFree(backend->d_active_in);
    if (backend->d_active_out != 0) cudaFree(backend->d_active_out);
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
    const cuDoubleComplex *blocks) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    plugin_operator *op = 0;
    size_t row_ptr_len = 0U;
    size_t block_len = 0U;

    if (backend == 0 || op_name == 0 || row_ptr == 0 || col_ind == 0 || site_types == 0 || blocks == 0) {
        return 1;
    }
    if (block_dim != backend->block_dim || n_sites != backend->n_sites) {
        return 1;
    }

    op = upsert_operator(backend, op_name);
    if (op == 0 || op->name == 0) {
        return 1;
    }

    op->block_dim = block_dim;
    op->n_sites = n_sites;
    op->nnzb = nnzb;
    row_ptr_len = (size_t) (n_sites + 1);
    block_len = (size_t) block_dim * (size_t) block_dim * (size_t) nnzb;

    op->row_ptr = (int *) malloc(row_ptr_len * sizeof(*op->row_ptr));
    op->col_ind = (int *) malloc((size_t) nnzb * sizeof(*op->col_ind));
    op->site_types = (int *) malloc((size_t) n_sites * sizeof(*op->site_types));
    op->blocks = (cuDoubleComplex *) malloc(block_len * sizeof(*op->blocks));
    if (op->row_ptr == 0 || op->col_ind == 0 || op->site_types == 0 || op->blocks == 0) {
        free_operator(op);
        return 1;
    }

    memcpy(op->row_ptr, row_ptr, row_ptr_len * sizeof(*op->row_ptr));
    memcpy(op->col_ind, col_ind, (size_t) nnzb * sizeof(*op->col_ind));
    memcpy(op->site_types, site_types, (size_t) n_sites * sizeof(*op->site_types));
    memcpy(op->blocks, blocks, block_len * sizeof(*op->blocks));

    if (cuda_ok(cudaMalloc((void **) &op->d_row_ptr, row_ptr_len * sizeof(*op->d_row_ptr))) != 0 ||
        cuda_ok(cudaMalloc((void **) &op->d_col_ind, (size_t) nnzb * sizeof(*op->d_col_ind))) != 0 ||
        cuda_ok(cudaMalloc((void **) &op->d_site_types, (size_t) n_sites * sizeof(*op->d_site_types))) != 0 ||
        cuda_ok(cudaMalloc((void **) &op->d_blocks, block_len * sizeof(*op->d_blocks))) != 0) {
        free_operator(op);
        return 1;
    }

    if (cuda_ok(cudaMemcpy(op->d_row_ptr, op->row_ptr, row_ptr_len * sizeof(*op->d_row_ptr), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(op->d_col_ind, op->col_ind, (size_t) nnzb * sizeof(*op->d_col_ind), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(op->d_site_types, op->site_types, (size_t) n_sites * sizeof(*op->d_site_types), cudaMemcpyHostToDevice)) != 0 ||
        cuda_ok(cudaMemcpy(op->d_blocks, op->blocks, block_len * sizeof(*op->d_blocks), cudaMemcpyHostToDevice)) != 0) {
        free_operator(op);
        return 1;
    }

    return 0;
}

__global__ static void apply_block_operator_kernel(
    int block_dim,
    int n_sites,
    const int *row_ptr,
    const int *col_ind,
    const cuDoubleComplex *blocks,
    char trans_mode,
    const int *active_in,
    const cuDoubleComplex *psi_in,
    cuDoubleComplex *psi_out,
    int *active_out) {
    unsigned long long tid = (unsigned long long) blockIdx.x * (unsigned long long) blockDim.x + (unsigned long long) threadIdx.x;
    unsigned long long block_size = (unsigned long long) block_dim * (unsigned long long) block_dim;
    unsigned long long state_len = block_size * (unsigned long long) n_sites;
    unsigned long long row_site = 0U;
    unsigned long long local = 0U;
    int row_i = 0;
    int row_j = 0;
    int entry = 0;
    int touched = 0;
    cuDoubleComplex accum = c_zero();

    if (tid >= state_len) {
        return;
    }

    row_site = tid / block_size;
    local = tid % block_size;
    row_i = (int) (local % (unsigned long long) block_dim);
    row_j = (int) (local / (unsigned long long) block_dim);

    for (entry = row_ptr[row_site] - 1; entry < row_ptr[row_site + 1] - 1; ++entry) {
        int col = col_ind[entry] - 1;
        int l = 0;
        if (col < 0 || col >= n_sites || active_in[col] == 0) {
            continue;
        }
        touched = 1;
        for (l = 0; l < block_dim; ++l) {
            unsigned long long psi_idx = (unsigned long long) l +
                                         (unsigned long long) row_j * (unsigned long long) block_dim +
                                         (unsigned long long) col * block_size;
            cuDoubleComplex opv;
            if (trans_mode == 'c' || trans_mode == 'C') {
                unsigned long long op_idx = (unsigned long long) l +
                                            (unsigned long long) row_i * (unsigned long long) block_dim +
                                            (unsigned long long) entry * block_size;
                opv = cuConj(blocks[op_idx]);
            } else {
                unsigned long long op_idx = (unsigned long long) row_i +
                                            (unsigned long long) l * (unsigned long long) block_dim +
                                            (unsigned long long) entry * block_size;
                opv = blocks[op_idx];
            }
            accum = cuCadd(accum, cuCmul(opv, psi_in[psi_idx]));
        }
    }

    psi_out[tid] = accum;
    if (touched != 0) {
        atomicExch(active_out + row_site, 1);
    }
}

extern "C" int rslmto_backend_apply_operator(
    void *backend_handle,
    const char *op_name,
    char trans_mode,
    int block_dim,
    int n_sites,
    const cuDoubleComplex *psi_in,
    const int *active_in,
    cuDoubleComplex *psi_out,
    int *active_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    if (backend == 0 || block_dim != backend->block_dim || n_sites != backend->n_sites) {
        return 1;
    }
    return apply_named_operator(backend, op_name, trans_mode, psi_in, active_in, psi_out, active_out);
}

extern "C" int rslmto_backend_run_scalar_lanczos(
    void *backend,
    int llmax,
    int n_targets,
    const int *target_sites,
    int hoh_enabled,
    double *a_out,
    double *b2_out) {
    (void) backend;
    (void) llmax;
    (void) n_targets;
    (void) target_sites;
    (void) hoh_enabled;
    (void) a_out;
    (void) b2_out;
    return 1;
}

extern "C" int rslmto_backend_run_block_lanczos(
    void *backend_handle,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    cuDoubleComplex *a_b_out,
    cuDoubleComplex *b2_b_out,
    double *a_diag_out,
    double *b2_diag_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psi = 0;
    cuDoubleComplex *psi_prev_b = 0;
    cuDoubleComplex *hpsi = 0;
    cuDoubleComplex *pmn = 0;
    cuDoubleComplex *psi_t = 0;
    cuDoubleComplex *sum_b = 0;
    cuDoubleComplex *sum_a = 0;
    cuDoubleComplex *b = 0;
    cuDoubleComplex *b_inv = 0;
    cuDoubleComplex *tmp = 0;
    int *active = 0;
    int *active_out = 0;
    int wf = 0;
    int ll = 0;
    int site = 0;
    int i = 0;
    int j = 0;

    if (backend == 0 || llmax <= 0 || n_workflows <= 0) {
        return 1;
    }

    psi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi));
    psi_prev_b = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi_prev_b));
    hpsi = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*hpsi));
    pmn = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*pmn));
    psi_t = (cuDoubleComplex *) malloc(state_len_of(backend) * sizeof(*psi_t));
    sum_b = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*sum_b));
    sum_a = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*sum_a));
    b = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*b));
    b_inv = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*b_inv));
    tmp = (cuDoubleComplex *) malloc((size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
    active = (int *) malloc((size_t) backend->n_sites * sizeof(*active));
    active_out = (int *) malloc((size_t) backend->n_sites * sizeof(*active_out));
    if (psi == 0 || psi_prev_b == 0 || hpsi == 0 || pmn == 0 || psi_t == 0 || sum_b == 0 || sum_a == 0 ||
        b == 0 || b_inv == 0 || tmp == 0 || active == 0 || active_out == 0) {
        free(psi); free(psi_prev_b); free(hpsi); free(pmn); free(psi_t); free(sum_b); free(sum_a); free(b); free(b_inv); free(tmp); free(active); free(active_out);
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
                free(psi); free(psi_prev_b); free(hpsi); free(pmn); free(psi_t); free(sum_b); free(sum_a); free(b); free(b_inv); free(tmp); free(active); free(active_out);
                return 1;
            }

            matrix_zero(sum_a, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                if (active_out[site] == 0) continue;
                {
                    size_t off = site_block_offset(backend, site);
                    int idx = 0;
                    for (idx = 0; idx < backend->block_dim * backend->block_dim; ++idx) {
                        pmn[off + idx] = c_sub(hpsi[off + idx], pmn[off + idx]);
                    }
                    accumulate_left_dagger_right(sum_a, psi + off, hpsi + off, backend->block_dim);
                }
            }

            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    a_b_out[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, ll, wf)] = sum_a[matrix_offset(backend->block_dim, backend->block_dim, i, j)];
                }
                a_diag_out[(size_t) ll + (size_t) llmax * (size_t) j + (size_t) llmax * (size_t) backend->block_dim * (size_t) wf] =
                    cuCreal(sum_a[matrix_offset(backend->block_dim, backend->block_dim, j, j)]);
            }

            state_copy(backend, psi_t, psi);
            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    b2_b_out[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, ll, wf)] = sum_b[matrix_offset(backend->block_dim, backend->block_dim, i, j)];
                }
                b2_diag_out[(size_t) ll + (size_t) llmax * (size_t) j + (size_t) llmax * (size_t) backend->block_dim * (size_t) wf] =
                    cuCreal(sum_b[matrix_offset(backend->block_dim, backend->block_dim, j, j)]);
            }

            matrix_zero(sum_b, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                if (active_out[site] == 0) continue;
                {
                    size_t off = site_block_offset(backend, site);
                    matrix_right_multiply(tmp, psi + off, sum_a, backend->block_dim);
                    for (i = 0; i < backend->block_dim * backend->block_dim; ++i) {
                        pmn[off + i] = c_sub(pmn[off + i], tmp[i]);
                    }
                    accumulate_left_dagger_right(sum_b, pmn + off, pmn + off, backend->block_dim);
                }
            }

            if (hermitian_sqrt_and_inverse(backend->block_dim, sum_b, b, b_inv) != 0) {
                free(psi); free(psi_prev_b); free(hpsi); free(pmn); free(psi_t); free(sum_b); free(sum_a); free(b); free(b_inv); free(tmp); free(active); free(active_out);
                return 1;
            }

            for (site = 0; site < backend->n_sites; ++site) {
                if (active_out[site] == 0) continue;
                {
                    size_t off = site_block_offset(backend, site);
                    matrix_right_multiply(tmp, pmn + off, b_inv, backend->block_dim);
                    memcpy(psi + off, tmp, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
                    matrix_right_multiply(tmp, psi_t + off, b, backend->block_dim);
                    memcpy(pmn + off, tmp, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*tmp));
                }
            }

            memcpy(active, active_out, (size_t) backend->n_sites * sizeof(*active));
        }

        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                b2_b_out[a_b_offset(backend->block_dim, llmax, n_workflows, i, j, llmax - 1, wf)] = sum_b[matrix_offset(backend->block_dim, backend->block_dim, i, j)];
            }
            b2_diag_out[(size_t) (llmax - 1) + (size_t) llmax * (size_t) j + (size_t) llmax * (size_t) backend->block_dim * (size_t) wf] =
                cuCreal(sum_b[matrix_offset(backend->block_dim, backend->block_dim, j, j)]);
        }
    }

    free(psi); free(psi_prev_b); free(hpsi); free(pmn); free(psi_t); free(sum_b); free(sum_a); free(b); free(b_inv); free(tmp); free(active); free(active_out);
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
    cuDoubleComplex *mu_out) {
    plugin_backend *backend = (plugin_backend *) backend_handle;
    cuDoubleComplex *psi0 = 0;
    cuDoubleComplex *psi1 = 0;
    cuDoubleComplex *psi2 = 0;
    cuDoubleComplex *psiref = 0;
    cuDoubleComplex *dum1 = 0;
    cuDoubleComplex *dum2 = 0;
    int *active0 = 0;
    int *active = 0;
    int *active_out = 0;
    int wf = 0;
    int ll = 0;
    int site = 0;
    int i = 0;
    int j = 0;
    int nmom = 2 * llmax + 2;

    if (backend == 0 || llmax <= 0 || n_workflows <= 0) {
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
    if (psi0 == 0 || psi1 == 0 || psi2 == 0 || psiref == 0 || dum1 == 0 || dum2 == 0 || active0 == 0 || active == 0 || active_out == 0) {
        free(psi0); free(psi1); free(psi2); free(psiref); free(dum1); free(dum2); free(active0); free(active); free(active_out);
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
            if (active0[site] == 0) continue;
            accumulate_left_dagger_right(dum1, psiref + site_block_offset(backend, site), psi0 + site_block_offset(backend, site), backend->block_dim);
        }
        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                mu_out[mu_offset(backend->block_dim, nmom, i, j, 0, wf)] = dum1[matrix_offset(backend->block_dim, backend->block_dim, i, j)];
            }
        }

        if (apply_hamiltonian(backend, hoh_enabled, psi0, active0, a_scale, b_shift, psi1, active_out) != 0) {
            free(psi0); free(psi1); free(psi2); free(psiref); free(dum1); free(dum2); free(active0); free(active); free(active_out);
            return 1;
        }

        matrix_zero(dum1, backend->block_dim);
        for (site = 0; site < backend->n_sites; ++site) {
            if (active0[site] == 0) continue;
            accumulate_left_dagger_right(dum1, psiref + site_block_offset(backend, site), psi1 + site_block_offset(backend, site), backend->block_dim);
        }
        for (j = 0; j < backend->block_dim; ++j) {
            for (i = 0; i < backend->block_dim; ++i) {
                mu_out[mu_offset(backend->block_dim, nmom, i, j, 1, wf)] = dum1[matrix_offset(backend->block_dim, backend->block_dim, i, j)];
            }
        }

        memcpy(active, active_out, (size_t) backend->n_sites * sizeof(*active));
        for (ll = 0; ll < llmax; ++ll) {
            if (apply_hamiltonian(backend, hoh_enabled, psi1, active, a_scale, b_shift, psi2, active_out) != 0) {
                free(psi0); free(psi1); free(psi2); free(psiref); free(dum1); free(dum2); free(active0); free(active); free(active_out);
                return 1;
            }
            {
                size_t idx = 0;
                size_t total = state_len_of(backend);
                for (idx = 0; idx < total; ++idx) {
                    psi2[idx] = c_scale(psi2[idx], 2.0);
                }
            }

            matrix_zero(dum1, backend->block_dim);
            matrix_zero(dum2, backend->block_dim);
            for (site = 0; site < backend->n_sites; ++site) {
                if (active_out[site] == 0) continue;
                {
                    size_t off = site_block_offset(backend, site);
                    int idx = 0;
                    for (idx = 0; idx < backend->block_dim * backend->block_dim; ++idx) {
                        psi2[off + idx] = c_sub(psi2[off + idx], psi0[off + idx]);
                    }
                    accumulate_left_dagger_right(dum1, psi1 + off, psi1 + off, backend->block_dim);
                    accumulate_left_dagger_right(dum2, psi2 + off, psi1 + off, backend->block_dim);
                    memcpy(psi0 + off, psi1 + off, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*psi0));
                    memcpy(psi1 + off, psi2 + off, (size_t) backend->block_dim * (size_t) backend->block_dim * sizeof(*psi1));
                }
            }

            for (j = 0; j < backend->block_dim; ++j) {
                for (i = 0; i < backend->block_dim; ++i) {
                    mu_out[mu_offset(backend->block_dim, nmom, i, j, 2 * ll + 2, wf)] =
                        c_sub(c_scale(dum1[matrix_offset(backend->block_dim, backend->block_dim, i, j)], 2.0),
                              mu_out[mu_offset(backend->block_dim, nmom, i, j, 0, wf)]);
                    mu_out[mu_offset(backend->block_dim, nmom, i, j, 2 * ll + 3, wf)] =
                        c_sub(c_scale(dum2[matrix_offset(backend->block_dim, backend->block_dim, i, j)], 2.0),
                              mu_out[mu_offset(backend->block_dim, nmom, i, j, 1, wf)]);
                }
            }

            memcpy(active, active_out, (size_t) backend->n_sites * sizeof(*active));
        }
    }

    free(psi0); free(psi1); free(psi2); free(psiref); free(dum1); free(dum2); free(active0); free(active); free(active_out);
    return 0;
}

extern "C" int rslmto_backend_run_transport_stochastic(
    void *backend,
    int llmax,
    int loop_over,
    int calc_mode,
    int hoh_enabled,
    const int *type_sites,
    double a_scale,
    double b_shift,
    cuDoubleComplex *mu_nm_out) {
    (void) backend;
    (void) llmax;
    (void) loop_over;
    (void) calc_mode;
    (void) hoh_enabled;
    (void) type_sites;
    (void) a_scale;
    (void) b_shift;
    (void) mu_nm_out;
    return 1;
}

extern "C" const char *rslmto_backend_name(void) {
    return "rslmto_cuda_sparse_plugin";
}

extern "C" const char *rslmto_backend_capabilities(void) {
    return "execution=nvidia_gpu;library=cuda;precision=complex_fp64;operators=block_sparse,velocity,hoh;recurrences=block_chebyshev,block_lanczos,apply_only_transport;device_buffers=operator,state;kernel=custom_block_csr_apply";
}
