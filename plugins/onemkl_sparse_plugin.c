#include "backend_plugin_api.h"

#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <mkl_spblas.h>

typedef struct {
    char *name;
    int block_dim;
    int n_sites;
    int nnzb;
    int *row_ptr;
    int *col_ind;
    int *site_types;
    double _Complex *blocks;
    MKL_INT scalar_dim;
    MKL_INT scalar_nnz;
    MKL_INT *csr_row_start;
    MKL_INT *csr_row_end;
    MKL_INT *csr_col_ind;
    MKL_Complex16 *csr_values;
    sparse_matrix_t csr_handle;
    struct matrix_descr descr;
} plugin_operator;

typedef struct {
    char *precision_mode;
    char *library_name;
    int block_dim;
    int n_sites;
    int capacity;
    int count;
    plugin_operator *ops;
    MKL_INT packed_rows;
    MKL_INT packed_cols;
    MKL_Complex16 *packed_in;
    MKL_Complex16 *packed_out;
} plugin_backend;

static MKL_Complex16 to_mkl_complex(double real_part, double imag_part) {
    MKL_Complex16 value;
    memset(&value, 0, sizeof(value));
    ((double *) &value)[0] = real_part;
    ((double *) &value)[1] = imag_part;
    return value;
}

static char *dup_string(const char *input) {
    size_t len = 0;
    char *copy = NULL;
    if (input == NULL) {
        return NULL;
    }
    len = strlen(input);
    copy = (char *) malloc(len + 1);
    if (copy == NULL) {
        return NULL;
    }
    memcpy(copy, input, len + 1);
    return copy;
}

static void free_operator(plugin_operator *op) {
    if (op == NULL) {
        return;
    }
    if (op->csr_handle != NULL) {
        mkl_sparse_destroy(op->csr_handle);
    }
    free(op->name);
    free(op->row_ptr);
    free(op->col_ind);
    free(op->site_types);
    free(op->blocks);
    free(op->csr_row_start);
    free(op->csr_row_end);
    free(op->csr_col_ind);
    free(op->csr_values);
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

static int ensure_packed_buffers(plugin_backend *backend, MKL_INT rows, MKL_INT cols) {
    size_t count = 0;
    if (backend == NULL) {
        return 1;
    }
    if (backend->packed_rows == rows && backend->packed_cols == cols &&
        backend->packed_in != NULL && backend->packed_out != NULL) {
        return 0;
    }

    free(backend->packed_in);
    free(backend->packed_out);
    backend->packed_in = NULL;
    backend->packed_out = NULL;
    backend->packed_rows = 0;
    backend->packed_cols = 0;

    count = (size_t) rows * (size_t) cols;
    backend->packed_in = (MKL_Complex16 *) calloc(count, sizeof(*backend->packed_in));
    backend->packed_out = (MKL_Complex16 *) calloc(count, sizeof(*backend->packed_out));
    if (backend->packed_in == NULL || backend->packed_out == NULL) {
        free(backend->packed_in);
        free(backend->packed_out);
        backend->packed_in = NULL;
        backend->packed_out = NULL;
        return 1;
    }

    backend->packed_rows = rows;
    backend->packed_cols = cols;
    return 0;
}

static void pack_block_state(
    MKL_Complex16 *packed,
    MKL_INT packed_rows,
    int block_dim,
    int n_sites,
    const double _Complex *psi_in,
    const int *active_in) {
    int site = 0;
    int col = 0;
    int row = 0;
    memset(packed, 0, (size_t) packed_rows * (size_t) block_dim * sizeof(*packed));
    for (site = 0; site < n_sites; ++site) {
        if (active_in[site] == 0) {
            continue;
        }
        for (col = 0; col < block_dim; ++col) {
            for (row = 0; row < block_dim; ++row) {
                size_t packed_idx = (size_t) (site * block_dim + row) + (size_t) packed_rows * (size_t) col;
                size_t psi_idx = (size_t) row + (size_t) block_dim * (size_t) col +
                                 (size_t) block_dim * (size_t) block_dim * (size_t) site;
                memcpy(&packed[packed_idx], &psi_in[psi_idx], sizeof(packed[packed_idx]));
            }
        }
    }
}

static void unpack_block_state(
    double _Complex *psi_out,
    int block_dim,
    int n_sites,
    const MKL_Complex16 *packed,
    MKL_INT packed_rows) {
    int site = 0;
    int col = 0;
    int row = 0;
    for (site = 0; site < n_sites; ++site) {
        for (col = 0; col < block_dim; ++col) {
            for (row = 0; row < block_dim; ++row) {
                size_t packed_idx = (size_t) (site * block_dim + row) + (size_t) packed_rows * (size_t) col;
                size_t psi_idx = (size_t) row + (size_t) block_dim * (size_t) col +
                                 (size_t) block_dim * (size_t) block_dim * (size_t) site;
                memcpy(&psi_out[psi_idx], &packed[packed_idx], sizeof(packed[packed_idx]));
            }
        }
    }
}

static void compute_active_out(
    const plugin_operator *op,
    char trans_mode,
    const int *active_in,
    int *active_out) {
    int row = 0;
    memset(active_out, 0, (size_t) op->n_sites * sizeof(*active_out));
    if (trans_mode == 'c' || trans_mode == 'C') {
        for (row = 0; row < op->n_sites; ++row) {
            int entry = 0;
            if (active_in[row] == 0) {
                continue;
            }
            for (entry = op->row_ptr[row] - 1; entry < op->row_ptr[row + 1] - 1; ++entry) {
                int col = op->col_ind[entry] - 1;
                if (col >= 0 && col < op->n_sites) {
                    active_out[col] = 1;
                }
            }
        }
    } else {
        for (row = 0; row < op->n_sites; ++row) {
            int entry = 0;
            for (entry = op->row_ptr[row] - 1; entry < op->row_ptr[row + 1] - 1; ++entry) {
                int col = op->col_ind[entry] - 1;
                if (col >= 0 && col < op->n_sites && active_in[col] != 0) {
                    active_out[row] = 1;
                    break;
                }
            }
        }
    }
}

static sparse_operation_t to_mkl_operation(char trans_mode) {
    if (trans_mode == 'c' || trans_mode == 'C') {
        return SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    }
    return SPARSE_OPERATION_NON_TRANSPOSE;
}

void *rslmto_backend_create(const char *precision_mode, int block_dim, int n_sites, const char *library_name) {
    plugin_backend *backend = NULL;
    if (precision_mode == NULL || strcmp(precision_mode, "complex_fp64") != 0) {
        return NULL;
    }
    backend = (plugin_backend *) calloc(1, sizeof(*backend));
    if (backend == NULL) {
        return NULL;
    }
    backend->precision_mode = dup_string(precision_mode);
    backend->library_name = dup_string(library_name != NULL ? library_name : "onemkl");
    backend->block_dim = block_dim;
    backend->n_sites = n_sites;
    return backend;
}

int rslmto_backend_destroy(void *backend_handle) {
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
    free(backend->packed_in);
    free(backend->packed_out);
    free(backend);
    return 0;
}

int rslmto_backend_upload_operator(
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
    MKL_INT scalar_dim = 0;
    MKL_INT scalar_nnz = 0;
    MKL_INT value_cursor = 0;
    int row_site = 0;
    sparse_status_t status = SPARSE_STATUS_SUCCESS;

    if (backend == NULL || op_name == NULL || row_ptr == NULL || col_ind == NULL || site_types == NULL || blocks == NULL) {
        return 1;
    }
    if (block_dim != backend->block_dim || n_sites != backend->n_sites) {
        return 1;
    }

    op = upsert_operator(backend, op_name);
    if (op == NULL || op->name == NULL) {
        return 1;
    }

    op->block_dim = block_dim;
    op->n_sites = n_sites;
    op->nnzb = nnzb;
    scalar_dim = (MKL_INT) block_dim * (MKL_INT) n_sites;
    scalar_nnz = (MKL_INT) block_dim * (MKL_INT) block_dim * (MKL_INT) nnzb;
    op->scalar_dim = scalar_dim;
    op->scalar_nnz = scalar_nnz;

    op->row_ptr = (int *) malloc((size_t) (n_sites + 1) * sizeof(*op->row_ptr));
    op->col_ind = (int *) malloc((size_t) nnzb * sizeof(*op->col_ind));
    op->site_types = (int *) malloc((size_t) n_sites * sizeof(*op->site_types));
    op->blocks = (double _Complex *) malloc((size_t) block_dim * (size_t) block_dim * (size_t) nnzb * sizeof(*op->blocks));
    op->csr_row_start = (MKL_INT *) malloc((size_t) scalar_dim * sizeof(*op->csr_row_start));
    op->csr_row_end = (MKL_INT *) malloc((size_t) scalar_dim * sizeof(*op->csr_row_end));
    op->csr_col_ind = (MKL_INT *) malloc((size_t) scalar_nnz * sizeof(*op->csr_col_ind));
    op->csr_values = (MKL_Complex16 *) malloc((size_t) scalar_nnz * sizeof(*op->csr_values));
    if (op->row_ptr == NULL || op->col_ind == NULL || op->site_types == NULL || op->blocks == NULL ||
        op->csr_row_start == NULL || op->csr_row_end == NULL || op->csr_col_ind == NULL || op->csr_values == NULL) {
        free_operator(op);
        return 1;
    }

    memcpy(op->row_ptr, row_ptr, (size_t) (n_sites + 1) * sizeof(*op->row_ptr));
    memcpy(op->col_ind, col_ind, (size_t) nnzb * sizeof(*op->col_ind));
    memcpy(op->site_types, site_types, (size_t) n_sites * sizeof(*op->site_types));
    memcpy(op->blocks, blocks, (size_t) block_dim * (size_t) block_dim * (size_t) nnzb * sizeof(*op->blocks));

    for (row_site = 0; row_site < n_sites; ++row_site) {
        int local_row = 0;
        for (local_row = 0; local_row < block_dim; ++local_row) {
            MKL_INT scalar_row = (MKL_INT) row_site * (MKL_INT) block_dim + (MKL_INT) local_row;
            int entry = 0;
            op->csr_row_start[scalar_row] = value_cursor + 1;
            for (entry = row_ptr[row_site] - 1; entry < row_ptr[row_site + 1] - 1; ++entry) {
                int col_site = col_ind[entry] - 1;
                int local_col = 0;
                if (col_site < 0 || col_site >= n_sites) {
                    free_operator(op);
                    return 1;
                }
                for (local_col = 0; local_col < block_dim; ++local_col) {
                    size_t block_idx = (size_t) local_row + (size_t) block_dim * (size_t) local_col +
                                       (size_t) block_dim * (size_t) block_dim * (size_t) entry;
                    MKL_INT scalar_col = (MKL_INT) col_site * (MKL_INT) block_dim + (MKL_INT) local_col + 1;
                    op->csr_col_ind[value_cursor] = scalar_col;
                    memcpy(&op->csr_values[value_cursor], &blocks[block_idx], sizeof(op->csr_values[value_cursor]));
                    value_cursor += 1;
                }
            }
            op->csr_row_end[scalar_row] = value_cursor + 1;
        }
    }
    if (value_cursor != scalar_nnz) {
        free_operator(op);
        return 1;
    }

    op->descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    op->descr.mode = SPARSE_FILL_MODE_FULL;
    op->descr.diag = SPARSE_DIAG_NON_UNIT;

    status = mkl_sparse_z_create_csr(
        &op->csr_handle,
        SPARSE_INDEX_BASE_ONE,
        scalar_dim,
        scalar_dim,
        op->csr_row_start,
        op->csr_row_end,
        op->csr_col_ind,
        op->csr_values
    );
    if (status != SPARSE_STATUS_SUCCESS) {
        free_operator(op);
        return 1;
    }

    status = mkl_sparse_set_mm_hint(
        op->csr_handle,
        SPARSE_OPERATION_NON_TRANSPOSE,
        op->descr,
        SPARSE_LAYOUT_COLUMN_MAJOR,
        (MKL_INT) block_dim,
        64
    );
    if (status != SPARSE_STATUS_SUCCESS) {
        free_operator(op);
        return 1;
    }
    status = mkl_sparse_set_mm_hint(
        op->csr_handle,
        SPARSE_OPERATION_CONJUGATE_TRANSPOSE,
        op->descr,
        SPARSE_LAYOUT_COLUMN_MAJOR,
        (MKL_INT) block_dim,
        64
    );
    if (status != SPARSE_STATUS_SUCCESS) {
        free_operator(op);
        return 1;
    }
    status = mkl_sparse_optimize(op->csr_handle);
    if (status != SPARSE_STATUS_SUCCESS) {
        free_operator(op);
        return 1;
    }
    return 0;
}

int rslmto_backend_apply_operator(
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
    plugin_operator *op = NULL;
    MKL_INT scalar_dim = 0;
    MKL_Complex16 alpha;
    MKL_Complex16 beta;
    sparse_status_t status = SPARSE_STATUS_SUCCESS;

    if (backend == NULL || psi_in == NULL || active_in == NULL || psi_out == NULL || active_out == NULL) {
        return 1;
    }

    op = find_operator(backend, op_name);
    if (op == NULL || op->csr_handle == NULL) {
        return 1;
    }
    if (block_dim != backend->block_dim || n_sites != backend->n_sites) {
        return 1;
    }

    scalar_dim = op->scalar_dim;
    if (ensure_packed_buffers(backend, scalar_dim, (MKL_INT) block_dim) != 0) {
        return 1;
    }

    pack_block_state(backend->packed_in, scalar_dim, block_dim, n_sites, psi_in, active_in);
    memset(backend->packed_out, 0, (size_t) scalar_dim * (size_t) block_dim * sizeof(*backend->packed_out));

    alpha = to_mkl_complex(1.0, 0.0);
    beta = to_mkl_complex(0.0, 0.0);

    status = mkl_sparse_z_mm(
        to_mkl_operation(trans_mode),
        alpha,
        op->csr_handle,
        op->descr,
        SPARSE_LAYOUT_COLUMN_MAJOR,
        backend->packed_in,
        (MKL_INT) block_dim,
        scalar_dim,
        beta,
        backend->packed_out,
        scalar_dim
    );
    if (status != SPARSE_STATUS_SUCCESS) {
        return 1;
    }

    unpack_block_state(psi_out, block_dim, n_sites, backend->packed_out, scalar_dim);
    compute_active_out(op, trans_mode, active_in, active_out);
    return 0;
}

const char *rslmto_backend_name(void) {
    return "rslmto_onemkl_sparse_plugin";
}

const char *rslmto_backend_capabilities(void) {
    return "execution=vendor_cpu;library=oneMKL;precision=complex_fp64;operators=block_sparse,velocity,hoh;recurrences=apply_only;lowering=scalar_csr;kernel=mkl_sparse_z_mm;inspector_executor=yes";
}
