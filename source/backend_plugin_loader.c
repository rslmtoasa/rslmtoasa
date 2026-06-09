#include "backend_plugin_api.h"

#include <dlfcn.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef void *(*backend_create_fn)(const char *, int, int, const char *);
typedef int (*backend_destroy_fn)(void *);
typedef int (*backend_upload_operator_fn)(void *, const char *, int, int, int, const int *, const int *, const int *, const double _Complex *);
typedef int (*backend_apply_operator_fn)(void *, const char *, char, int, int, const double _Complex *, const int *, double _Complex *, int *);
typedef int (*backend_run_scalar_lanczos_fn)(void *, int, int, const int *, int, double *, double *);
typedef int (*backend_run_block_lanczos_fn)(void *, int, int, const int *, const int *, const int *, int, double _Complex *, double _Complex *, double *, double *);
typedef int (*backend_run_block_chebyshev_fn)(void *, int, int, const int *, const int *, const int *, int, double, double, double _Complex *);
typedef int (*backend_run_transport_stochastic_fn)(void *, int, int, int, int, const int *, double, double, double _Complex *);
typedef const char *(*backend_name_fn)(void);
typedef const char *(*backend_capabilities_fn)(void);

typedef struct {
    void *dl_handle;
    backend_create_fn create;
    backend_destroy_fn destroy;
    backend_upload_operator_fn upload_operator;
    backend_apply_operator_fn apply_operator;
    backend_run_scalar_lanczos_fn run_scalar_lanczos;
    backend_run_block_lanczos_fn run_block_lanczos;
    backend_run_block_chebyshev_fn run_block_chebyshev;
    backend_run_transport_stochastic_fn run_transport_stochastic;
    backend_name_fn name;
    backend_capabilities_fn capabilities;
} rslmto_backend_plugin_runtime;

static void set_error(char *errbuf, int errbuf_len, const char *message) {
    if (errbuf == NULL || errbuf_len <= 0) {
        return;
    }
    if (message == NULL) {
        errbuf[0] = '\0';
        return;
    }
    snprintf(errbuf, (size_t) errbuf_len, "%s", message);
}

static int load_symbol(void *handle, const char *symbol, void **out, char *errbuf, int errbuf_len) {
    dlerror();
    *out = dlsym(handle, symbol);
    if (*out == NULL) {
        set_error(errbuf, errbuf_len, dlerror());
        return 1;
    }
    return 0;
}

void *rslmto_backend_plugin_open(const char *path, char *errbuf, int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = NULL;
    void *dl_handle = NULL;

    set_error(errbuf, errbuf_len, "");
    if (path == NULL || path[0] == '\0') {
        set_error(errbuf, errbuf_len, "empty plugin path");
        return NULL;
    }

    dl_handle = dlopen(path, RTLD_NOW | RTLD_LOCAL);
    if (dl_handle == NULL) {
        set_error(errbuf, errbuf_len, dlerror());
        return NULL;
    }

    runtime = (rslmto_backend_plugin_runtime *) calloc(1, sizeof(*runtime));
    if (runtime == NULL) {
        set_error(errbuf, errbuf_len, "failed to allocate plugin runtime");
        dlclose(dl_handle);
        return NULL;
    }

    runtime->dl_handle = dl_handle;
    if (load_symbol(dl_handle, "rslmto_backend_create", (void **) &runtime->create, errbuf, errbuf_len) != 0 ||
        load_symbol(dl_handle, "rslmto_backend_destroy", (void **) &runtime->destroy, errbuf, errbuf_len) != 0 ||
        load_symbol(dl_handle, "rslmto_backend_upload_operator", (void **) &runtime->upload_operator, errbuf, errbuf_len) != 0 ||
        load_symbol(dl_handle, "rslmto_backend_apply_operator", (void **) &runtime->apply_operator, errbuf, errbuf_len) != 0 ||
        load_symbol(dl_handle, "rslmto_backend_name", (void **) &runtime->name, errbuf, errbuf_len) != 0) {
        dlclose(dl_handle);
        free(runtime);
        return NULL;
    }
    dlerror();
    runtime->capabilities = (backend_capabilities_fn) dlsym(dl_handle, "rslmto_backend_capabilities");
    runtime->run_scalar_lanczos = (backend_run_scalar_lanczos_fn) dlsym(dl_handle, "rslmto_backend_run_scalar_lanczos");
    runtime->run_block_lanczos = (backend_run_block_lanczos_fn) dlsym(dl_handle, "rslmto_backend_run_block_lanczos");
    runtime->run_block_chebyshev = (backend_run_block_chebyshev_fn) dlsym(dl_handle, "rslmto_backend_run_block_chebyshev");
    runtime->run_transport_stochastic = (backend_run_transport_stochastic_fn) dlsym(dl_handle, "rslmto_backend_run_transport_stochastic");

    return runtime;
}

void rslmto_backend_plugin_close(void *runtime_handle) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    if (runtime == NULL) {
        return;
    }
    /* Keep the image loaded until process exit. Final Fortran teardown can still
       reference plugin-owned code paths, and eager dlclose caused shutdown traps. */
    free(runtime);
}

void *rslmto_backend_plugin_create_backend(
    void *runtime_handle,
    const char *precision_mode,
    int block_dim,
    int n_sites,
    const char *library_name,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    void *backend = NULL;

    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL || runtime->create == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return NULL;
    }

    backend = runtime->create(precision_mode, block_dim, n_sites, library_name);
    if (backend == NULL) {
        set_error(errbuf, errbuf_len, "plugin backend creation failed");
    }
    return backend;
}

int rslmto_backend_plugin_destroy_backend(void *runtime_handle, void *backend_handle, char *errbuf, int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL || runtime->destroy == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    return runtime->destroy(backend_handle);
}

int rslmto_backend_plugin_upload_operator(
    void *runtime_handle,
    void *backend_handle,
    const char *op_name,
    int block_dim,
    int n_sites,
    int nnzb,
    const int *row_ptr,
    const int *col_ind,
    const int *site_types,
    const double _Complex *blocks,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL || runtime->upload_operator == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    return runtime->upload_operator(backend_handle, op_name, block_dim, n_sites, nnzb, row_ptr, col_ind, site_types, blocks);
}

int rslmto_backend_plugin_apply_operator(
    void *runtime_handle,
    void *backend_handle,
    const char *op_name,
    char trans_mode,
    int block_dim,
    int n_sites,
    const double _Complex *psi_in,
    const int *active_in,
    double _Complex *psi_out,
    int *active_out,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL || runtime->apply_operator == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    return runtime->apply_operator(backend_handle, op_name, trans_mode, block_dim, n_sites, psi_in, active_in, psi_out, active_out);
}

int rslmto_backend_plugin_run_scalar_lanczos(
    void *runtime_handle,
    void *backend_handle,
    int llmax,
    int n_targets,
    const int *target_sites,
    int hoh_enabled,
    double *a_out,
    double *b2_out,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    if (runtime->run_scalar_lanczos == NULL) {
        set_error(errbuf, errbuf_len, "plugin does not export rslmto_backend_run_scalar_lanczos");
        return 1;
    }
    return runtime->run_scalar_lanczos(backend_handle, llmax, n_targets, target_sites, hoh_enabled, a_out, b2_out);
}

int rslmto_backend_plugin_run_block_lanczos(
    void *runtime_handle,
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
    double *b2_diag_out,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    if (runtime->run_block_lanczos == NULL) {
        set_error(errbuf, errbuf_len, "plugin does not export rslmto_backend_run_block_lanczos");
        return 1;
    }
    return runtime->run_block_lanczos(backend_handle, llmax, n_workflows, site_i, site_j, variant_id, hoh_enabled, a_b_out, b2_b_out, a_diag_out, b2_diag_out);
}

int rslmto_backend_plugin_run_block_chebyshev(
    void *runtime_handle,
    void *backend_handle,
    int llmax,
    int n_workflows,
    const int *site_i,
    const int *site_j,
    const int *variant_id,
    int hoh_enabled,
    double a_scale,
    double b_shift,
    double _Complex *mu_out,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    if (runtime->run_block_chebyshev == NULL) {
        set_error(errbuf, errbuf_len, "plugin does not export rslmto_backend_run_block_chebyshev");
        return 1;
    }
    return runtime->run_block_chebyshev(backend_handle, llmax, n_workflows, site_i, site_j, variant_id, hoh_enabled, a_scale, b_shift, mu_out);
}

int rslmto_backend_plugin_run_transport_stochastic(
    void *runtime_handle,
    void *backend_handle,
    int llmax,
    int loop_over,
    int calc_mode,
    int hoh_enabled,
    const int *type_sites,
    double a_scale,
    double b_shift,
    double _Complex *mu_nm_out,
    char *errbuf,
    int errbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    set_error(errbuf, errbuf_len, "");
    if (runtime == NULL) {
        set_error(errbuf, errbuf_len, "plugin runtime is not initialized");
        return 1;
    }
    if (runtime->run_transport_stochastic == NULL) {
        set_error(errbuf, errbuf_len, "plugin does not export rslmto_backend_run_transport_stochastic");
        return 1;
    }
    return runtime->run_transport_stochastic(backend_handle, llmax, loop_over, calc_mode, hoh_enabled, type_sites, a_scale, b_shift, mu_nm_out);
}

void rslmto_backend_plugin_name(void *runtime_handle, char *namebuf, int namebuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    const char *name = "unknown";
    if (runtime != NULL && runtime->name != NULL) {
        name = runtime->name();
    }
    set_error(namebuf, namebuf_len, name);
}

void rslmto_backend_plugin_capabilities(void *runtime_handle, char *capbuf, int capbuf_len) {
    rslmto_backend_plugin_runtime *runtime = (rslmto_backend_plugin_runtime *) runtime_handle;
    const char *capabilities = "capabilities_unavailable";
    if (runtime != NULL && runtime->capabilities != NULL) {
        capabilities = runtime->capabilities();
    }
    set_error(capbuf, capbuf_len, capabilities);
}
