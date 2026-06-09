#include "backend_plugin_api.h"

#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char *name;
    int block_dim;
    int n_sites;
    int nnzb;
    int *row_ptr;
    int *col_ind;
    int *site_types;
    double _Complex *blocks;
} plugin_operator;

typedef struct {
    char *precision_mode;
    char *library_name;
    int block_dim;
    int n_sites;
    int capacity;
    int count;
    plugin_operator *ops;
} plugin_backend;

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
    free(op->name);
    free(op->row_ptr);
    free(op->col_ind);
    free(op->site_types);
    free(op->blocks);
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

void *rslmto_backend_create(const char *precision_mode, int block_dim, int n_sites, const char *library_name) {
    plugin_backend *backend = (plugin_backend *) calloc(1, sizeof(*backend));
    if (backend == NULL) {
        return NULL;
    }
    backend->precision_mode = dup_string(precision_mode != NULL ? precision_mode : "complex_fp64");
    backend->library_name = dup_string(library_name != NULL ? library_name : "cpu_sparse_plugin");
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
    size_t row_ptr_len = 0;
    size_t block_len = 0;

    if (backend == NULL || op_name == NULL || row_ptr == NULL || col_ind == NULL || site_types == NULL || blocks == NULL) {
        return 1;
    }

    op = upsert_operator(backend, op_name);
    if (op == NULL || op->name == NULL) {
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
    op->blocks = (double _Complex *) malloc(block_len * sizeof(*op->blocks));
    if (op->row_ptr == NULL || op->col_ind == NULL || op->site_types == NULL || op->blocks == NULL) {
        free_operator(op);
        return 1;
    }

    memcpy(op->row_ptr, row_ptr, row_ptr_len * sizeof(*op->row_ptr));
    memcpy(op->col_ind, col_ind, (size_t) nnzb * sizeof(*op->col_ind));
    memcpy(op->site_types, site_types, (size_t) n_sites * sizeof(*op->site_types));
    memcpy(op->blocks, blocks, block_len * sizeof(*op->blocks));
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
    int row = 0;
    int entry = 0;
    int col = 0;
    int i = 0;
    int j = 0;
    int l = 0;
    size_t state_len = (size_t) block_dim * (size_t) block_dim * (size_t) n_sites;

    if (backend == NULL || psi_in == NULL || active_in == NULL || psi_out == NULL || active_out == NULL) {
        return 1;
    }

    op = find_operator(backend, op_name);
    if (op == NULL) {
        return 1;
    }

    memset(psi_out, 0, state_len * sizeof(*psi_out));
    memset(active_out, 0, (size_t) n_sites * sizeof(*active_out));

    for (row = 0; row < op->n_sites; ++row) {
        int touched = 0;
        for (entry = op->row_ptr[row] - 1; entry < op->row_ptr[row + 1] - 1; ++entry) {
            col = op->col_ind[entry] - 1;
            if (col < 0 || col >= n_sites) {
                continue;
            }
            if (active_in[col] == 0) {
                continue;
            }
            touched = 1;
            for (j = 0; j < block_dim; ++j) {
                for (i = 0; i < block_dim; ++i) {
                    double _Complex accum = 0.0;
                    for (l = 0; l < block_dim; ++l) {
                        size_t op_idx = 0;
                        size_t psi_idx = 0;
                        op_idx = (size_t) i + (size_t) l * (size_t) block_dim + (size_t) entry * (size_t) block_dim * (size_t) block_dim;
                        psi_idx = (size_t) l + (size_t) j * (size_t) block_dim + (size_t) col * (size_t) block_dim * (size_t) block_dim;
                        if (trans_mode == 'c' || trans_mode == 'C') {
                            size_t op_conj_idx = (size_t) l + (size_t) i * (size_t) block_dim + (size_t) entry * (size_t) block_dim * (size_t) block_dim;
                            accum += conj(op->blocks[op_conj_idx]) * psi_in[psi_idx];
                        } else {
                            accum += op->blocks[op_idx] * psi_in[psi_idx];
                        }
                    }
                    psi_out[(size_t) i + (size_t) j * (size_t) block_dim + (size_t) row * (size_t) block_dim * (size_t) block_dim] += accum;
                }
            }
        }
        if (touched) {
            active_out[row] = 1;
        }
    }

    return 0;
}

const char *rslmto_backend_name(void) {
    return "rslmto_cpu_sparse_plugin";
}

const char *rslmto_backend_capabilities(void) {
    return "execution=serial_cpu;precision=complex_fp64;operators=block_sparse,velocity,hoh;recurrences=apply_only;acceleration=none";
}
