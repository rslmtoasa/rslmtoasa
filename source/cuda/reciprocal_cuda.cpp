#include "reciprocal_cuda.h"

#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

extern "C" int rslmto_reciprocal_cuda_launch_lehmann(
    cudaStream_t stream,
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
extern "C" int rslmto_reciprocal_cuda_launch_lehmann_device(
    cudaStream_t stream,
    int nmat,
    int nk,
    int ne,
    int npair,
    int nblk,
    const double *device_eigenvalues,
    const void *device_eigenvectors,
    const double *host_k_points,
    const void *host_z_contour,
    const double *host_dr,
    const int *host_ioffset,
    const int *host_joffset,
    void *host_blocks,
    double *h2d_seconds,
    double *contraction_seconds,
    double *d2h_seconds);

struct rslmto_reciprocal_cuda_context {
    int device = -1;
    int prepared_operator_generation = -1;
    cudaStream_t stream = nullptr;
    cusolverDnHandle_t solver = nullptr;
    syevjInfo_t jacobi_params = nullptr;
    int solver_strategy = RSLMTO_RECIPROCAL_CUDA_ZHEEVD_SERIAL;
    int workspace_n = 0;
    int workspace_batch = 0;
    int workspace_lwork = 0;
    int workspace_jobz = -1;
    int workspace_strategy = -1;
    int workspace_query_n = 0;
    int workspace_query_batch = 0;
    std::size_t workspace_matrix_capacity_bytes = 0;
    std::size_t workspace_eigenvalue_capacity_bytes = 0;
    std::size_t workspace_info_capacity_bytes = 0;
    cuDoubleComplex *device_hamiltonians = nullptr;
    double *device_eigenvalues = nullptr;
    cuDoubleComplex *device_work = nullptr;
    int *device_info = nullptr;
    cuComplex *device_hamiltonians_fp32 = nullptr;
    float *device_eigenvalues_fp32 = nullptr;
    cuComplex *device_work_fp32 = nullptr;
    int *device_info_fp32 = nullptr;
    /* Reusable host staging for the FP32 route.  These vectors grow only when
     * the requested shape grows or changes; they are never allocated per
     * matrix. */
    std::vector<cuComplex> host_hamiltonians_fp32;
    std::vector<float> host_eigenvalues_fp32;
    std::vector<int> host_info;
    /* Optional, backend-owned pinned staging.  It is enabled explicitly by
     * RSLMTO_CUDA_PINNED_HOST=1 and only used for n >= 486; application
     * arrays are never registered or pinned. */
    bool pinned_host_enabled = false;
    int pinned_host_active = 0;
    static constexpr int pinned_host_min_n = 486;
    cuDoubleComplex *pinned_hamiltonians_fp64 = nullptr;
    double *pinned_eigenvalues_fp64 = nullptr;
    cuComplex *pinned_hamiltonians_fp32 = nullptr;
    float *pinned_eigenvalues_fp32 = nullptr;
    std::size_t pinned_matrix_capacity = 0;
    std::size_t pinned_eigenvalue_capacity = 0;
    std::size_t pinned_fp32_matrix_capacity = 0;
    std::size_t pinned_fp32_eigenvalue_capacity = 0;
    /* Events are context-lifetime resources.  They are recorded repeatedly
     * but never created/destroyed in the solve hot path. */
    cudaEvent_t event_h2d_start = nullptr;
    cudaEvent_t event_h2d_stop = nullptr;
    cudaEvent_t event_solver_start = nullptr;
    cudaEvent_t event_solver_stop = nullptr;
    cudaEvent_t event_d2h_values_start = nullptr;
    cudaEvent_t event_d2h_values_stop = nullptr;
    cudaEvent_t event_d2h_vectors_start = nullptr;
    cudaEvent_t event_d2h_vectors_stop = nullptr;
    cudaEvent_t event_d2h_stop = nullptr;
    bool timing_events_ready = false;
    double h2d_seconds = 0.0;
    double solve_seconds = 0.0;
    double d2h_seconds = 0.0;
    double d2h_values_seconds = 0.0;
    double d2h_vectors_seconds = 0.0;
    double sync_seconds = 0.0;
    double host_staging_seconds = 0.0;
    double backend_wall_seconds = 0.0;
    double host_conversion_seconds = 0.0;
    double host_widen_seconds = 0.0;
    long long h2d_bytes = 0;
    long long d2h_values_bytes = 0;
    long long d2h_vectors_bytes = 0;
    int timing_calls = 0;
    long long cuda_malloc_count = 0;
    long long cuda_free_count = 0;
    long long workspace_query_count = 0;
    long long workspace_reuse_count = 0;
    long long event_create_count = 0;
    long long event_destroy_count = 0;
    long long pinned_alloc_count = 0;
    long long pinned_free_count = 0;
    int resident_token = 0;
    int resident_n = 0;
    int resident_batch = 0;
    bool resident_valid = false;
};

static std::string g_reciprocal_cuda_error;

static void set_error(const char *message) {
    g_reciprocal_cuda_error = message;
}

static void set_cuda_error(const char *operation, cudaError_t status) {
    g_reciprocal_cuda_error = std::string(operation) + ": " +
                              cudaGetErrorString(status);
}

static void set_solver_error(const char *operation, cusolverStatus_t status) {
    g_reciprocal_cuda_error = std::string(operation) + ": cuSOLVER status " +
                              std::to_string(static_cast<int>(status));
}

template <typename T>
static void release_device_buffer(T *&pointer, long long *free_count) {
    if (!pointer) return;
    const cudaError_t status = cudaFree(pointer);
    if (status == cudaSuccess && free_count) ++(*free_count);
    pointer = nullptr;
}

template <typename T>
static bool ensure_pinned_buffer(T *&pointer, std::size_t &capacity, std::size_t required,
                                 long long *alloc_count, long long *free_count) {
    if (capacity >= required) return true;
    if (pointer) {
        const cudaError_t free_status = cudaFreeHost(pointer);
        if (free_status == cudaSuccess && free_count) ++(*free_count);
        pointer = nullptr;
        capacity = 0;
    }
    if (required == 0) return true;
    const cudaError_t status = cudaHostAlloc(reinterpret_cast<void **>(&pointer),
                                              required * sizeof(T), cudaHostAllocDefault);
    if (status != cudaSuccess) return false;
    capacity = required;
    if (alloc_count) ++(*alloc_count);
    return true;
}

static void destroy_timing_events(rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    cudaEvent_t *events[] = {
        &context->event_h2d_start, &context->event_h2d_stop,
        &context->event_solver_start, &context->event_solver_stop,
        &context->event_d2h_values_start, &context->event_d2h_values_stop,
        &context->event_d2h_vectors_start, &context->event_d2h_vectors_stop,
        &context->event_d2h_stop};
    for (cudaEvent_t *event : events) {
        if (!*event) continue;
        if (cudaEventDestroy(*event) == cudaSuccess) ++context->event_destroy_count;
        *event = nullptr;
    }
    context->timing_events_ready = false;
}

static bool create_timing_events(rslmto_reciprocal_cuda_context *context) {
    if (!context) return false;
    cudaEvent_t *events[] = {
        &context->event_h2d_start, &context->event_h2d_stop,
        &context->event_solver_start, &context->event_solver_stop,
        &context->event_d2h_values_start, &context->event_d2h_values_stop,
        &context->event_d2h_vectors_start, &context->event_d2h_vectors_stop,
        &context->event_d2h_stop};
    for (cudaEvent_t *event : events) {
        if (cudaEventCreate(event) != cudaSuccess) {
            destroy_timing_events(context);
            return false;
        }
        ++context->event_create_count;
    }
    context->timing_events_ready = true;
    return true;
}

static bool record_timing_event(cudaEvent_t event, cudaStream_t stream, const char *operation) {
    const cudaError_t status = cudaEventRecord(event, stream);
    if (status != cudaSuccess) {
        set_cuda_error(operation, status);
        return false;
    }
    return true;
}

static bool elapsed_timing_event(float *milliseconds, cudaEvent_t start, cudaEvent_t stop,
                                 const char *operation) {
    const cudaError_t status = cudaEventElapsedTime(milliseconds, start, stop);
    if (status != cudaSuccess) {
        set_cuda_error(operation, status);
        return false;
    }
    return true;
}

static void release_pinned_buffers(rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    auto free_host = [&](auto *&pointer, std::size_t &capacity) {
        if (!pointer) return;
        if (cudaFreeHost(pointer) == cudaSuccess) ++context->pinned_free_count;
        pointer = nullptr;
        capacity = 0;
    };
    free_host(context->pinned_hamiltonians_fp64, context->pinned_matrix_capacity);
    free_host(context->pinned_eigenvalues_fp64, context->pinned_eigenvalue_capacity);
    free_host(context->pinned_hamiltonians_fp32, context->pinned_fp32_matrix_capacity);
    free_host(context->pinned_eigenvalues_fp32, context->pinned_fp32_eigenvalue_capacity);
}

static void release_solver_workspace(rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    release_device_buffer(context->device_hamiltonians, &context->cuda_free_count);
    release_device_buffer(context->device_eigenvalues, &context->cuda_free_count);
    release_device_buffer(context->device_work, &context->cuda_free_count);
    release_device_buffer(context->device_info, &context->cuda_free_count);
    release_device_buffer(context->device_hamiltonians_fp32, &context->cuda_free_count);
    release_device_buffer(context->device_eigenvalues_fp32, &context->cuda_free_count);
    release_device_buffer(context->device_work_fp32, &context->cuda_free_count);
    release_device_buffer(context->device_info_fp32, &context->cuda_free_count);
    context->workspace_n = 0;
    context->workspace_batch = 0;
    context->workspace_lwork = 0;
    context->workspace_jobz = -1;
    context->workspace_strategy = -1;
    context->workspace_query_n = 0;
    context->workspace_query_batch = 0;
    context->workspace_matrix_capacity_bytes = 0;
    context->workspace_eigenvalue_capacity_bytes = 0;
    context->workspace_info_capacity_bytes = 0;
}

extern "C" const char *rslmto_reciprocal_cuda_last_error(void) {
    return g_reciprocal_cuda_error.c_str();
}

extern "C" int rslmto_reciprocal_cuda_device_count(int *count) {
    if (!count) {
        set_error("rslmto_reciprocal_cuda_device_count: null output");
        return 1;
    }

    const cudaError_t status = cudaGetDeviceCount(count);
    if (status != cudaSuccess) {
        *count = 0;
        set_cuda_error("cudaGetDeviceCount", status);
        return 1;
    }
    return 0;
}

extern "C" rslmto_reciprocal_cuda_context *
rslmto_reciprocal_cuda_create(int device) {
    int device_count = 0;
    if (rslmto_reciprocal_cuda_device_count(&device_count) != 0) {
        return nullptr;
    }
    if (device < 0 || device >= device_count) {
        set_error("rslmto_reciprocal_cuda_create: invalid device index");
        return nullptr;
    }

    const cudaError_t set_status = cudaSetDevice(device);
    if (set_status != cudaSuccess) {
        set_cuda_error("cudaSetDevice", set_status);
        return nullptr;
    }

    auto *context = new rslmto_reciprocal_cuda_context();
    context->device = device;
    const char *pinned_setting = std::getenv("RSLMTO_CUDA_PINNED_HOST");
    context->pinned_host_enabled = pinned_setting &&
        (std::strcmp(pinned_setting, "1") == 0 ||
         std::strcmp(pinned_setting, "true") == 0 ||
         std::strcmp(pinned_setting, "TRUE") == 0);
    const cudaError_t stream_status = cudaStreamCreate(&context->stream);
    if (stream_status != cudaSuccess) {
        set_cuda_error("cudaStreamCreate", stream_status);
        delete context;
        return nullptr;
    }
    cusolverStatus_t solver_status = cusolverDnCreate(&context->solver);
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error("cusolverDnCreate", solver_status);
        cudaStreamDestroy(context->stream);
        delete context;
        return nullptr;
    }
    solver_status = cusolverDnSetStream(context->solver, context->stream);
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error("cusolverDnSetStream", solver_status);
        cusolverDnDestroy(context->solver);
        cudaStreamDestroy(context->stream);
        delete context;
        return nullptr;
    }
    solver_status = cusolverDnCreateSyevjInfo(&context->jacobi_params);
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error("cusolverDnCreateSyevjInfo", solver_status);
        cusolverDnDestroy(context->solver);
        cudaStreamDestroy(context->stream);
        delete context;
        return nullptr;
    }
    /* Keep the documented/default Jacobi tolerance and sweep limit.  The
     * benchmark needs sorted eigenvalues, which is also cuSOLVER's default;
     * set it explicitly so the output contract is not toolkit-version
     * dependent. */
    solver_status = cusolverDnXsyevjSetSortEig(context->jacobi_params, 1);
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error("cusolverDnXsyevjSetSortEig", solver_status);
        cusolverDnDestroySyevjInfo(context->jacobi_params);
        cusolverDnDestroy(context->solver);
        cudaStreamDestroy(context->stream);
        delete context;
        return nullptr;
    }
    /* Timing events are created once per backend context.  A toolkit that
     * cannot create them must not make the numerical route fail; the host
     * wall timer and transfer counters remain available. */
    create_timing_events(context);
    return context;
}

extern "C" int rslmto_reciprocal_cuda_prepare_operator(
    rslmto_reciprocal_cuda_context *context, int operator_generation) {
    if (!context) {
        set_error("rslmto_reciprocal_cuda_prepare_operator: null context");
        return -1;
    }
    if (context->prepared_operator_generation == operator_generation) {
        return 0;
    }

    /* The standard eigensolver currently has no generation-specific operator
     * data. Keep the marker in the backend-owned context so later reciprocal
     * stages can reuse the same lifecycle seam. */
    context->prepared_operator_generation = operator_generation;
    return 1;
}

extern "C" int rslmto_reciprocal_cuda_set_solver_strategy(
    rslmto_reciprocal_cuda_context *context, int strategy) {
    if (!context || !context->solver) {
        set_error("rslmto_reciprocal_cuda_set_solver_strategy: null context");
        return 1;
    }
    if (strategy != RSLMTO_RECIPROCAL_CUDA_ZHEEVD_SERIAL &&
        strategy != RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED &&
        strategy != RSLMTO_RECIPROCAL_CUDA_CHEEVD_SERIAL &&
        strategy != RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) {
        set_error("rslmto_reciprocal_cuda_set_solver_strategy: unknown strategy");
        return 1;
    }
    if (context->solver_strategy != strategy) {
        release_solver_workspace(context);
        context->solver_strategy = strategy;
    }
    return 0;
}

extern "C" int rslmto_reciprocal_cuda_solver_strategy_supported(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    int request_eigenvectors) {
    if (!context || !context->solver || !context->jacobi_params) {
        set_error("rslmto_reciprocal_cuda_solver_strategy_supported: null context");
        return -1;
    }
    if (n < 1 || batch_size < 1 || (request_eigenvectors != 0 && request_eigenvectors != 1)) {
        set_error("rslmto_reciprocal_cuda_solver_strategy_supported: invalid dimensions");
        return -1;
    }
    if ((context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED ||
         context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) && n > 32) {
        set_error("cheevj/zheevj_batched unsupported: cuSOLVER requires n <= 32");
        return 1;
    }
    return 0;
}

extern "C" int rslmto_reciprocal_cuda_solve_zheevd_batch(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    const void *host_hamiltonians,
    double *host_eigenvalues,
    void *host_eigenvectors,
    int request_eigenvectors) {
    if (!context || !context->solver) {
        set_error("rslmto_reciprocal_cuda_solve_zheevd_batch: null context");
        return 1;
    }
    if (context->solver_strategy != RSLMTO_RECIPROCAL_CUDA_ZHEEVD_SERIAL &&
        context->solver_strategy != RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED) {
        set_error("rslmto_reciprocal_cuda_solve_zheevd_batch: selected strategy is FP32");
        return 2;
    }
    if (n < 1 || batch_size < 1 || !host_hamiltonians || !host_eigenvalues ||
        (request_eigenvectors && !host_eigenvectors)) {
        set_error("rslmto_reciprocal_cuda_solve_zheevd_batch: invalid arguments");
        return 1;
    }
    const int strategy = context->solver_strategy;
    const int support = rslmto_reciprocal_cuda_solver_strategy_supported(
        context, n, batch_size, request_eigenvectors);
    if (support != 0) return support > 0 ? 2 : 1;

    /* Any new solve invalidates the previous handoff immediately.  The
     * monotonically increasing token prevents an old Fortran reciprocal
     * object from accidentally matching a later solve with the same shape. */
    context->resident_valid = false;
    context->resident_n = 0;
    context->resident_batch = 0;

    const std::size_t matrix_elements = static_cast<std::size_t>(n) * n;
    const std::size_t matrix_bytes = matrix_elements * batch_size * sizeof(cuDoubleComplex);
    const std::size_t eigenvalue_bytes = static_cast<std::size_t>(n) * batch_size * sizeof(double);
    const int jobz = request_eigenvectors ? 1 : 0;
    const auto request_start = std::chrono::steady_clock::now();
    cudaError_t cuda_status = cudaSuccess;
    if (matrix_bytes > context->workspace_matrix_capacity_bytes ||
        eigenvalue_bytes > context->workspace_eigenvalue_capacity_bytes ||
        static_cast<std::size_t>(batch_size) * sizeof(int) > context->workspace_info_capacity_bytes ||
        !context->device_hamiltonians || !context->device_eigenvalues ||
        !context->device_info || strategy != context->workspace_strategy) {
        release_solver_workspace(context);
        context->workspace_n = n;
        context->workspace_batch = batch_size;

        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_hamiltonians), matrix_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_hamiltonians)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_matrix_capacity_bytes = matrix_bytes;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_eigenvalues), eigenvalue_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_eigenvalues)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_eigenvalue_capacity_bytes = eigenvalue_bytes;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_info),
                                  static_cast<std::size_t>(batch_size) * sizeof(int));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_info)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_info_capacity_bytes = static_cast<std::size_t>(batch_size) * sizeof(int);

        int lwork = 0;
        cusolverStatus_t solver_status;
        if (strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED) {
            solver_status = cusolverDnZheevjBatched_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians, n,
                context->device_eigenvalues, &lwork, context->jacobi_params, batch_size);
        } else {
            solver_status = cusolverDnZheevd_bufferSize(
                context->solver, request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians, n,
                context->device_eigenvalues, &lwork);
        }
        if (solver_status != CUSOLVER_STATUS_SUCCESS || lwork < 1) {
            set_solver_error(strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED
                                 ? "cusolverDnZheevjBatched_bufferSize"
                                 : "cusolverDnZheevd_bufferSize", solver_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->workspace_query_count;
        context->workspace_query_n = n;
        context->workspace_query_batch = batch_size;
        context->workspace_lwork = lwork;
        context->workspace_jobz = jobz;
        context->workspace_strategy = strategy;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_work),
                                 static_cast<std::size_t>(lwork) * sizeof(cuDoubleComplex));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_work)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
    } else {
        ++context->workspace_reuse_count;
    }

    /* A vector/value mode switch may require a larger cuSOLVER workspace,
     * but it must not discard the matrix/eigenvalue/info buffers.  Re-query
     * only at the mode boundary and grow the work buffer only when required. */
    if (n > context->workspace_query_n || batch_size > context->workspace_query_batch ||
        jobz != context->workspace_jobz) {
        int lwork = 0;
        cusolverStatus_t solver_status;
        if (strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED) {
            solver_status = cusolverDnZheevjBatched_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians, n,
                context->device_eigenvalues, &lwork, context->jacobi_params, batch_size);
        } else {
            solver_status = cusolverDnZheevd_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians, n,
                context->device_eigenvalues, &lwork);
        }
        ++context->workspace_query_count;
        if (solver_status != CUSOLVER_STATUS_SUCCESS || lwork < 1) {
            set_solver_error(strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED
                                 ? "cusolverDnZheevjBatched_bufferSize"
                                 : "cusolverDnZheevd_bufferSize", solver_status);
            return 1;
        }
        context->workspace_query_n = std::max(context->workspace_query_n, n);
        context->workspace_query_batch = std::max(context->workspace_query_batch, batch_size);
        if (lwork > context->workspace_lwork) {
            release_device_buffer(context->device_work, &context->cuda_free_count);
            cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_work),
                                     static_cast<std::size_t>(lwork) * sizeof(cuDoubleComplex));
            if (cuda_status != cudaSuccess) {
                set_cuda_error("cudaMalloc(device_work)", cuda_status);
                return 1;
            }
            ++context->cuda_malloc_count;
            context->workspace_lwork = lwork;
        }
        context->workspace_jobz = jobz;
    }

    const bool use_pinned = context->pinned_host_enabled && n >= rslmto_reciprocal_cuda_context::pinned_host_min_n;
    context->pinned_host_active = use_pinned ? 1 : 0;
    const void *h2d_source = host_hamiltonians;
    if (use_pinned) {
        if (!ensure_pinned_buffer(context->pinned_hamiltonians_fp64, context->pinned_matrix_capacity,
                                  matrix_elements * batch_size, &context->pinned_alloc_count,
                                  &context->pinned_free_count) ||
            !ensure_pinned_buffer(context->pinned_eigenvalues_fp64, context->pinned_eigenvalue_capacity,
                                  static_cast<std::size_t>(n) * batch_size, &context->pinned_alloc_count,
                                  &context->pinned_free_count)) {
            set_error("pinned host staging allocation failed; use pageable staging");
            return 1;
        }
        const auto staging_start = std::chrono::steady_clock::now();
        std::memcpy(context->pinned_hamiltonians_fp64, host_hamiltonians, matrix_bytes);
        context->host_staging_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - staging_start).count();
        h2d_source = context->pinned_hamiltonians_fp64;
    }
    if (context->timing_events_ready &&
        !record_timing_event(context->event_h2d_start, context->stream, "cudaEventRecord(H2D start)")) return 1;

    cuda_status = cudaMemcpyAsync(
        context->device_hamiltonians, h2d_source, matrix_bytes,
        cudaMemcpyHostToDevice, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(H2D Hamiltonians)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready) {
        if (!record_timing_event(context->event_h2d_stop, context->stream, "cudaEventRecord(H2D stop)")) return 1;
        if (!record_timing_event(context->event_solver_start, context->stream, "cudaEventRecord(solver start)")) return 1;
    }

    cusolverStatus_t solver_status = CUSOLVER_STATUS_SUCCESS;
    auto *device_matrices = context->device_hamiltonians;
    if (strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED) {
        /* This is the one true batched call.  The device buffers are laid out
         * as [n,n,batch] column-major matrices, exactly the layout described
         * by cuSOLVER's ZheevjBatched contract. */
        solver_status = cusolverDnZheevjBatched(
            context->solver,
            request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER, n, device_matrices, n,
            context->device_eigenvalues, context->device_work, context->workspace_lwork,
            context->device_info, context->jacobi_params, batch_size);
    } else {
        for (int ibatch = 0; ibatch < batch_size; ++ibatch) {
            solver_status = cusolverDnZheevd(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, device_matrices + static_cast<std::size_t>(ibatch) * matrix_elements, n,
                context->device_eigenvalues + static_cast<std::size_t>(ibatch) * n,
                context->device_work, context->workspace_lwork, context->device_info + ibatch);
            if (solver_status != CUSOLVER_STATUS_SUCCESS) break;
        }
    }
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error(strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED
                             ? "cusolverDnZheevjBatched"
                             : "cusolverDnZheevd", solver_status);
        return 1;
    }
    if (context->timing_events_ready) {
        if (!record_timing_event(context->event_solver_stop, context->stream, "cudaEventRecord(solver stop)")) return 1;
        if (!record_timing_event(context->event_d2h_values_start, context->stream, "cudaEventRecord(D2H values start)")) return 1;
    }

    void *d2h_values_target = host_eigenvalues;
    if (use_pinned) d2h_values_target = context->pinned_eigenvalues_fp64;
    cuda_status = cudaMemcpyAsync(d2h_values_target, context->device_eigenvalues,
                                  eigenvalue_bytes, cudaMemcpyDeviceToHost, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(D2H eigenvalues)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready &&
        !record_timing_event(context->event_d2h_values_stop, context->stream, "cudaEventRecord(D2H values stop)")) return 1;
    if (request_eigenvectors) {
        if (context->timing_events_ready &&
            !record_timing_event(context->event_d2h_vectors_start, context->stream, "cudaEventRecord(D2H vectors start)")) return 1;
        void *d2h_vectors_target = host_eigenvectors;
        if (use_pinned) d2h_vectors_target = context->pinned_hamiltonians_fp64;
        cuda_status = cudaMemcpyAsync(d2h_vectors_target, context->device_hamiltonians,
                                      matrix_bytes, cudaMemcpyDeviceToHost, context->stream);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMemcpyAsync(D2H eigenvectors)", cuda_status);
            return 1;
        }
        if (context->timing_events_ready &&
            !record_timing_event(context->event_d2h_vectors_stop, context->stream, "cudaEventRecord(D2H vectors stop)")) return 1;
    }
    context->host_info.assign(static_cast<std::size_t>(batch_size), 0);
    cuda_status = cudaMemcpyAsync(context->host_info.data(), context->device_info,
                                  static_cast<std::size_t>(batch_size) * sizeof(int),
                                  cudaMemcpyDeviceToHost, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(D2H solver info)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready &&
        !record_timing_event(context->event_d2h_stop, context->stream, "cudaEventRecord(D2H stop)")) return 1;
    const auto sync_start = std::chrono::steady_clock::now();
    cuda_status = cudaStreamSynchronize(context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaStreamSynchronize(eigensolution)", cuda_status);
        return 1;
    }
    context->sync_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - sync_start).count();
    if (use_pinned) {
        const auto staging_start = std::chrono::steady_clock::now();
        std::memcpy(host_eigenvalues, context->pinned_eigenvalues_fp64, eigenvalue_bytes);
        if (request_eigenvectors) std::memcpy(host_eigenvectors, context->pinned_hamiltonians_fp64, matrix_bytes);
        context->host_staging_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - staging_start).count();
    }

    for (int ibatch = 0; ibatch < batch_size; ++ibatch) {
        if (context->host_info[static_cast<std::size_t>(ibatch)] != 0) {
            g_reciprocal_cuda_error =
                std::string(strategy == RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED
                                ? "cusolverDnZheevjBatched"
                                : "cusolverDnZheevd") +
                ": solver info[" + std::to_string(ibatch) + "]=" +
                std::to_string(context->host_info[static_cast<std::size_t>(ibatch)]);
            return 1;
        }
    }
    if (context->timing_events_ready) {
        float h2d_ms = 0.0f, solve_ms = 0.0f, d2h_ms = 0.0f;
        float d2h_values_ms = 0.0f, d2h_vectors_ms = 0.0f;
        if (!elapsed_timing_event(&h2d_ms, context->event_h2d_start, context->event_h2d_stop, "cudaEventElapsedTime(H2D)")) return 1;
        if (!elapsed_timing_event(&solve_ms, context->event_solver_start, context->event_solver_stop, "cudaEventElapsedTime(solver)")) return 1;
        if (!elapsed_timing_event(&d2h_values_ms, context->event_d2h_values_start, context->event_d2h_values_stop, "cudaEventElapsedTime(D2H values)")) return 1;
        if (!elapsed_timing_event(&d2h_ms, context->event_d2h_values_start, context->event_d2h_stop, "cudaEventElapsedTime(D2H)")) return 1;
        if (request_eigenvectors &&
            !elapsed_timing_event(&d2h_vectors_ms, context->event_d2h_vectors_start, context->event_d2h_vectors_stop, "cudaEventElapsedTime(D2H vectors)")) return 1;
        context->h2d_seconds += static_cast<double>(h2d_ms) * 1.0e-3;
        context->solve_seconds += static_cast<double>(solve_ms) * 1.0e-3;
        context->d2h_seconds += static_cast<double>(d2h_ms) * 1.0e-3;
        context->d2h_values_seconds += static_cast<double>(d2h_values_ms) * 1.0e-3;
        context->d2h_vectors_seconds += static_cast<double>(d2h_vectors_ms) * 1.0e-3;
        ++context->timing_calls;
    }
    ++context->resident_token;
    if (request_eigenvectors) {
        context->resident_valid = true;
        context->resident_n = n;
        context->resident_batch = batch_size;
    }
    context->h2d_bytes += static_cast<long long>(matrix_bytes);
    context->d2h_values_bytes += static_cast<long long>(eigenvalue_bytes);
    if (request_eigenvectors) context->d2h_vectors_bytes += static_cast<long long>(matrix_bytes);
    context->backend_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - request_start).count();
    return 0;
}

extern "C" int rslmto_reciprocal_cuda_solve_cheevd_batch(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    const void *host_hamiltonians,
    double *host_eigenvalues,
    void *host_eigenvectors,
    int request_eigenvectors) {
    if (!context || !context->solver) {
        set_error("rslmto_reciprocal_cuda_solve_cheevd_batch: null context");
        return 1;
    }
    if (context->solver_strategy != RSLMTO_RECIPROCAL_CUDA_CHEEVD_SERIAL &&
        context->solver_strategy != RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) {
        set_error("rslmto_reciprocal_cuda_solve_cheevd_batch: selected strategy is FP64");
        return 2;
    }
    if (n < 1 || batch_size < 1 || !host_hamiltonians || !host_eigenvalues ||
        (request_eigenvectors && !host_eigenvectors)) {
        set_error("rslmto_reciprocal_cuda_solve_cheevd_batch: invalid arguments");
        return 1;
    }
    const int support = rslmto_reciprocal_cuda_solver_strategy_supported(
        context, n, batch_size, request_eigenvectors);
    if (support != 0) return support > 0 ? 2 : 1;

    context->resident_valid = false;
    context->resident_n = 0;
    context->resident_batch = 0;

    const std::size_t matrix_elements = static_cast<std::size_t>(n) * n;
    const std::size_t matrix_count = matrix_elements * batch_size;
    const std::size_t matrix_bytes = matrix_count * sizeof(cuComplex);
    const std::size_t eigenvalue_bytes = static_cast<std::size_t>(n) * batch_size * sizeof(float);
    const int jobz = request_eigenvectors ? 1 : 0;
    const auto request_start = std::chrono::steady_clock::now();
    const bool use_pinned = context->pinned_host_enabled && n >= rslmto_reciprocal_cuda_context::pinned_host_min_n;
    context->pinned_host_active = use_pinned ? 1 : 0;

    /* Convert H64 on the host before H2D.  This is separately timed and uses
     * persistent staging vectors so a matrix solve does not allocate buffers. */
    const auto conversion_start = std::chrono::steady_clock::now();
    cuComplex *host_hamiltonians_fp32 = nullptr;
    float *host_eigenvalues_fp32 = nullptr;
    if (use_pinned) {
        if (!ensure_pinned_buffer(context->pinned_hamiltonians_fp32, context->pinned_fp32_matrix_capacity,
                                  matrix_count, &context->pinned_alloc_count, &context->pinned_free_count) ||
            !ensure_pinned_buffer(context->pinned_eigenvalues_fp32, context->pinned_fp32_eigenvalue_capacity,
                                  static_cast<std::size_t>(n) * batch_size, &context->pinned_alloc_count,
                                  &context->pinned_free_count)) {
            set_error("pinned host staging allocation failed; use pageable staging");
            return 1;
        }
        host_hamiltonians_fp32 = context->pinned_hamiltonians_fp32;
        host_eigenvalues_fp32 = context->pinned_eigenvalues_fp32;
    } else {
        context->host_hamiltonians_fp32.resize(matrix_count);
        context->host_eigenvalues_fp32.resize(static_cast<std::size_t>(n) * batch_size);
        host_hamiltonians_fp32 = context->host_hamiltonians_fp32.data();
        host_eigenvalues_fp32 = context->host_eigenvalues_fp32.data();
    }
    const auto *source = static_cast<const cuDoubleComplex *>(host_hamiltonians);
    for (std::size_t i = 0; i < matrix_count; ++i) {
        host_hamiltonians_fp32[i] = make_cuComplex(
            static_cast<float>(source[i].x), static_cast<float>(source[i].y));
    }
    context->host_info.assign(static_cast<std::size_t>(batch_size), 0);
    context->host_conversion_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - conversion_start).count();

    cudaError_t cuda_status = cudaSuccess;
    if (matrix_bytes > context->workspace_matrix_capacity_bytes ||
        eigenvalue_bytes > context->workspace_eigenvalue_capacity_bytes ||
        static_cast<std::size_t>(batch_size) * sizeof(int) > context->workspace_info_capacity_bytes ||
        !context->device_hamiltonians_fp32 || !context->device_eigenvalues_fp32 ||
        !context->device_info_fp32 || context->solver_strategy != context->workspace_strategy) {
        release_solver_workspace(context);
        context->workspace_n = n;
        context->workspace_batch = batch_size;

        cuda_status = cudaMalloc(
            reinterpret_cast<void **>(&context->device_hamiltonians_fp32), matrix_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_hamiltonians_fp32)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_matrix_capacity_bytes = matrix_bytes;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_eigenvalues_fp32), eigenvalue_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_eigenvalues_fp32)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_eigenvalue_capacity_bytes = eigenvalue_bytes;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_info_fp32),
                                 static_cast<std::size_t>(batch_size) * sizeof(int));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_info_fp32)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
        context->workspace_info_capacity_bytes = static_cast<std::size_t>(batch_size) * sizeof(int);

        int lwork = 0;
        cusolverStatus_t solver_status;
        if (context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) {
            solver_status = cusolverDnCheevjBatched_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians_fp32, n,
                context->device_eigenvalues_fp32, &lwork, context->jacobi_params, batch_size);
        } else {
            solver_status = cusolverDnCheevd_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians_fp32, n,
                context->device_eigenvalues_fp32, &lwork);
        }
        if (solver_status != CUSOLVER_STATUS_SUCCESS || lwork < 1) {
            set_solver_error(context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED
                                 ? "cusolverDnCheevjBatched_bufferSize"
                                 : "cusolverDnCheevd_bufferSize", solver_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->workspace_query_count;
        context->workspace_query_n = n;
        context->workspace_query_batch = batch_size;
        context->workspace_lwork = lwork;
        context->workspace_jobz = jobz;
        context->workspace_strategy = context->solver_strategy;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_work_fp32),
                                 static_cast<std::size_t>(lwork) * sizeof(cuComplex));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_work_fp32)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        ++context->cuda_malloc_count;
    } else {
        ++context->workspace_reuse_count;
    }

    if (n > context->workspace_query_n || batch_size > context->workspace_query_batch ||
        jobz != context->workspace_jobz) {
        int lwork = 0;
        cusolverStatus_t solver_status;
        if (context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) {
            solver_status = cusolverDnCheevjBatched_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians_fp32, n,
                context->device_eigenvalues_fp32, &lwork, context->jacobi_params, batch_size);
        } else {
            solver_status = cusolverDnCheevd_bufferSize(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians_fp32, n,
                context->device_eigenvalues_fp32, &lwork);
        }
        ++context->workspace_query_count;
        if (solver_status != CUSOLVER_STATUS_SUCCESS || lwork < 1) {
            set_solver_error(context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED
                                 ? "cusolverDnCheevjBatched_bufferSize"
                                 : "cusolverDnCheevd_bufferSize", solver_status);
            return 1;
        }
        context->workspace_query_n = std::max(context->workspace_query_n, n);
        context->workspace_query_batch = std::max(context->workspace_query_batch, batch_size);
        if (lwork > context->workspace_lwork) {
            release_device_buffer(context->device_work_fp32, &context->cuda_free_count);
            cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_work_fp32),
                                     static_cast<std::size_t>(lwork) * sizeof(cuComplex));
            if (cuda_status != cudaSuccess) {
                set_cuda_error("cudaMalloc(device_work_fp32)", cuda_status);
                return 1;
            }
            ++context->cuda_malloc_count;
            context->workspace_lwork = lwork;
        }
        context->workspace_jobz = jobz;
    }

    if (context->timing_events_ready &&
        !record_timing_event(context->event_h2d_start, context->stream, "cudaEventRecord(H32 H2D start)")) return 1;

    cuda_status = cudaMemcpyAsync(
        context->device_hamiltonians_fp32, host_hamiltonians_fp32, matrix_bytes,
        cudaMemcpyHostToDevice, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(H32 H2D Hamiltonians)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready) {
        if (!record_timing_event(context->event_h2d_stop, context->stream, "cudaEventRecord(H32 H2D stop)")) return 1;
        if (!record_timing_event(context->event_solver_start, context->stream, "cudaEventRecord(H32 solver start)")) return 1;
    }

    cusolverStatus_t solver_status = CUSOLVER_STATUS_SUCCESS;
    if (context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) {
        solver_status = cusolverDnCheevjBatched(
            context->solver,
            request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians_fp32, n,
            context->device_eigenvalues_fp32, context->device_work_fp32,
            context->workspace_lwork, context->device_info_fp32,
            context->jacobi_params, batch_size);
    } else {
        for (int ibatch = 0; ibatch < batch_size; ++ibatch) {
            solver_status = cusolverDnCheevd(
                context->solver,
                request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
                CUBLAS_FILL_MODE_UPPER, n,
                context->device_hamiltonians_fp32 + static_cast<std::size_t>(ibatch) * matrix_elements,
                n, context->device_eigenvalues_fp32 + static_cast<std::size_t>(ibatch) * n,
                context->device_work_fp32, context->workspace_lwork,
                context->device_info_fp32 + ibatch);
            if (solver_status != CUSOLVER_STATUS_SUCCESS) break;
        }
    }
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        set_solver_error(context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED
                             ? "cusolverDnCheevjBatched" : "cusolverDnCheevd", solver_status);
        return 1;
    }
    if (context->timing_events_ready) {
        if (!record_timing_event(context->event_solver_stop, context->stream, "cudaEventRecord(H32 solver stop)")) return 1;
        if (!record_timing_event(context->event_d2h_values_start, context->stream, "cudaEventRecord(H32 D2H values start)")) return 1;
    }
    cuda_status = cudaMemcpyAsync(host_eigenvalues_fp32,
                                  context->device_eigenvalues_fp32, eigenvalue_bytes,
                                  cudaMemcpyDeviceToHost, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(D2H FP32 eigenvalues)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready &&
        !record_timing_event(context->event_d2h_values_stop, context->stream, "cudaEventRecord(H32 D2H values stop)")) return 1;
    if (request_eigenvectors) {
        if (context->timing_events_ready &&
            !record_timing_event(context->event_d2h_vectors_start, context->stream, "cudaEventRecord(H32 D2H vectors start)")) return 1;
        cuda_status = cudaMemcpyAsync(host_hamiltonians_fp32,
                                      context->device_hamiltonians_fp32, matrix_bytes,
                                      cudaMemcpyDeviceToHost, context->stream);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMemcpyAsync(D2H FP32 eigenvectors)", cuda_status);
            return 1;
        }
        if (context->timing_events_ready &&
            !record_timing_event(context->event_d2h_vectors_stop, context->stream, "cudaEventRecord(H32 D2H vectors stop)")) return 1;
    }
    cuda_status = cudaMemcpyAsync(context->host_info.data(), context->device_info_fp32,
                                  static_cast<std::size_t>(batch_size) * sizeof(int),
                                  cudaMemcpyDeviceToHost, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(D2H FP32 solver info)", cuda_status);
        return 1;
    }
    if (context->timing_events_ready &&
        !record_timing_event(context->event_d2h_stop, context->stream, "cudaEventRecord(H32 D2H stop)")) return 1;
    const auto sync_start = std::chrono::steady_clock::now();
    cuda_status = cudaStreamSynchronize(context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaStreamSynchronize(FP32 eigensolution)", cuda_status);
        return 1;
    }
    context->sync_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - sync_start).count();
    for (int ibatch = 0; ibatch < batch_size; ++ibatch) {
        if (context->host_info[static_cast<std::size_t>(ibatch)] != 0) {
            g_reciprocal_cuda_error =
                std::string(context->solver_strategy == RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED
                                ? "cusolverDnCheevjBatched" : "cusolverDnCheevd") +
                ": solver info[" + std::to_string(ibatch) + "]=" +
                std::to_string(context->host_info[static_cast<std::size_t>(ibatch)]);
            return 1;
        }
    }

    const auto widen_start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < static_cast<std::size_t>(n) * batch_size; ++i) {
        host_eigenvalues[i] = static_cast<double>(host_eigenvalues_fp32[i]);
    }
    if (request_eigenvectors) {
        auto *target = static_cast<cuDoubleComplex *>(host_eigenvectors);
        for (std::size_t i = 0; i < matrix_count; ++i) {
            target[i].x = static_cast<double>(host_hamiltonians_fp32[i].x);
            target[i].y = static_cast<double>(host_hamiltonians_fp32[i].y);
        }
    }
    context->host_widen_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - widen_start).count();

    if (context->timing_events_ready) {
        float h2d_ms = 0.0f, solve_ms = 0.0f, d2h_ms = 0.0f;
        float d2h_values_ms = 0.0f, d2h_vectors_ms = 0.0f;
        if (!elapsed_timing_event(&h2d_ms, context->event_h2d_start, context->event_h2d_stop, "cudaEventElapsedTime(H32 H2D)")) return 1;
        if (!elapsed_timing_event(&solve_ms, context->event_solver_start, context->event_solver_stop, "cudaEventElapsedTime(H32 solver)")) return 1;
        if (!elapsed_timing_event(&d2h_values_ms, context->event_d2h_values_start, context->event_d2h_values_stop, "cudaEventElapsedTime(H32 D2H values)")) return 1;
        if (!elapsed_timing_event(&d2h_ms, context->event_d2h_values_start, context->event_d2h_stop, "cudaEventElapsedTime(H32 D2H)")) return 1;
        if (request_eigenvectors &&
            !elapsed_timing_event(&d2h_vectors_ms, context->event_d2h_vectors_start, context->event_d2h_vectors_stop, "cudaEventElapsedTime(H32 D2H vectors)")) return 1;
        context->h2d_seconds += static_cast<double>(h2d_ms) * 1.0e-3;
        context->solve_seconds += static_cast<double>(solve_ms) * 1.0e-3;
        context->d2h_seconds += static_cast<double>(d2h_ms) * 1.0e-3;
        context->d2h_values_seconds += static_cast<double>(d2h_values_ms) * 1.0e-3;
        context->d2h_vectors_seconds += static_cast<double>(d2h_vectors_ms) * 1.0e-3;
        ++context->timing_calls;
    }
    ++context->resident_token;
    /* The FP32 buffers have a different layout/type and are deliberately not
     * eligible for the FP64 Lehmann resident handoff. */
    context->h2d_bytes += static_cast<long long>(matrix_bytes);
    context->d2h_values_bytes += static_cast<long long>(eigenvalue_bytes);
    if (request_eigenvectors) context->d2h_vectors_bytes += static_cast<long long>(matrix_bytes);
    context->backend_wall_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - request_start).count();
    return 0;
}

extern "C" int rslmto_reciprocal_cuda_contract_lehmann(
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
    double *d2h_seconds) {
    if (!context || !context->stream) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann: null context");
        return 1;
    }
    if (nmat < 1 || nk < 1 || ne < 1 || npair < 1 || nblk < 1 || nblk > nmat ||
        !host_eigenvalues || !host_eigenvectors || !host_k_points || !host_z_contour ||
        !host_dr || !host_ioffset || !host_joffset || !host_blocks) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann: invalid arguments");
        return 1;
    }
    const int status = rslmto_reciprocal_cuda_launch_lehmann(
        context->stream, nmat, nk, ne, npair, nblk, host_eigenvalues,
        host_eigenvectors, host_k_points, host_z_contour, host_dr,
        host_ioffset, host_joffset, host_blocks, h2d_seconds,
        contraction_seconds, d2h_seconds);
    if (status != 0) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann: CUDA kernel failed");
    }
    return status;
}

extern "C" int rslmto_reciprocal_cuda_get_resident_token(
    rslmto_reciprocal_cuda_context *context) {
    if (!context) {
        set_error("rslmto_reciprocal_cuda_get_resident_token: null context");
        return 0;
    }
    return context->resident_valid ? context->resident_token : 0;
}

extern "C" int rslmto_reciprocal_cuda_resident_eigensystem_matches(
    rslmto_reciprocal_cuda_context *context,
    int n,
    int batch_size,
    int token) {
    if (!context || !context->stream) {
        set_error("rslmto_reciprocal_cuda_resident_eigensystem_matches: null context");
        return 1;
    }
    return context->resident_valid && context->resident_n == n &&
                   context->resident_batch == batch_size &&
                   context->resident_token == token
               ? 0
               : 1;
}

extern "C" int rslmto_reciprocal_cuda_contract_lehmann_resident(
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
    double *d2h_seconds) {
    if (!context || !context->stream) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann_resident: null context");
        return 1;
    }
    if (nmat < 1 || nk < 1 || ne < 1 || npair < 1 || nblk < 1 || nblk > nmat ||
        !host_k_points || !host_z_contour || !host_dr || !host_ioffset ||
        !host_joffset || !host_blocks) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann_resident: invalid arguments");
        return 1;
    }
    if (rslmto_reciprocal_cuda_resident_eigensystem_matches(
            context, nmat, nk, resident_token) != 0) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann_resident: resident eigensystem is stale or unavailable");
        return 2;
    }
    const int status = rslmto_reciprocal_cuda_launch_lehmann_device(
        context->stream, nmat, nk, ne, npair, nblk,
        context->device_eigenvalues, context->device_hamiltonians,
        host_k_points, host_z_contour, host_dr, host_ioffset, host_joffset,
        host_blocks, h2d_seconds, contraction_seconds, d2h_seconds);
    if (status != 0) {
        set_error("rslmto_reciprocal_cuda_contract_lehmann_resident: CUDA kernel failed");
    }
    return status;
}

extern "C" void rslmto_reciprocal_cuda_get_timings(
    rslmto_reciprocal_cuda_context *context,
    double *h2d_seconds, double *solve_seconds, double *d2h_seconds, int *calls) {
    if (!context) return;
    if (h2d_seconds) *h2d_seconds = context->h2d_seconds;
    if (solve_seconds) *solve_seconds = context->solve_seconds;
    if (d2h_seconds) *d2h_seconds = context->d2h_seconds;
    if (calls) *calls = context->timing_calls;
}

extern "C" void rslmto_reciprocal_cuda_get_detailed_timings(
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
    int *calls) {
    if (!context) return;
    if (host_conversion_seconds) *host_conversion_seconds = context->host_conversion_seconds;
    if (host_staging_seconds) *host_staging_seconds = context->host_staging_seconds;
    if (h2d_seconds) *h2d_seconds = context->h2d_seconds;
    if (solve_seconds) *solve_seconds = context->solve_seconds;
    if (d2h_seconds) *d2h_seconds = context->d2h_seconds;
    if (d2h_values_seconds) *d2h_values_seconds = context->d2h_values_seconds;
    if (d2h_vectors_seconds) *d2h_vectors_seconds = context->d2h_vectors_seconds;
    if (sync_seconds) *sync_seconds = context->sync_seconds;
    if (host_widen_seconds) *host_widen_seconds = context->host_widen_seconds;
    if (total_seconds) *total_seconds = context->backend_wall_seconds;
    if (h2d_bytes) *h2d_bytes = context->h2d_bytes;
    if (d2h_values_bytes) *d2h_values_bytes = context->d2h_values_bytes;
    if (d2h_vectors_bytes) *d2h_vectors_bytes = context->d2h_vectors_bytes;
    if (pinned_host_active) *pinned_host_active = context->pinned_host_active;
    if (calls) *calls = context->timing_calls;
}

extern "C" void rslmto_reciprocal_cuda_get_resource_counters(
    rslmto_reciprocal_cuda_context *context,
    long long *cuda_malloc_count,
    long long *cuda_free_count,
    long long *workspace_query_count,
    long long *workspace_reuse_count,
    long long *event_create_count,
    long long *event_destroy_count,
    long long *pinned_alloc_count,
    long long *pinned_free_count) {
    if (!context) return;
    if (cuda_malloc_count) *cuda_malloc_count = context->cuda_malloc_count;
    if (cuda_free_count) *cuda_free_count = context->cuda_free_count;
    if (workspace_query_count) *workspace_query_count = context->workspace_query_count;
    if (workspace_reuse_count) *workspace_reuse_count = context->workspace_reuse_count;
    if (event_create_count) *event_create_count = context->event_create_count;
    if (event_destroy_count) *event_destroy_count = context->event_destroy_count;
    if (pinned_alloc_count) *pinned_alloc_count = context->pinned_alloc_count;
    if (pinned_free_count) *pinned_free_count = context->pinned_free_count;
}

extern "C" void rslmto_reciprocal_cuda_reset_timings(
    rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    context->h2d_seconds = 0.0;
    context->solve_seconds = 0.0;
    context->d2h_seconds = 0.0;
    context->d2h_values_seconds = 0.0;
    context->d2h_vectors_seconds = 0.0;
    context->sync_seconds = 0.0;
    context->host_staging_seconds = 0.0;
    context->backend_wall_seconds = 0.0;
    context->host_conversion_seconds = 0.0;
    context->host_widen_seconds = 0.0;
    context->h2d_bytes = 0;
    context->d2h_values_bytes = 0;
    context->d2h_vectors_bytes = 0;
    context->timing_calls = 0;
}

extern "C" int rslmto_reciprocal_cuda_get_memory(
    rslmto_reciprocal_cuda_context *context,
    long long *free_bytes, long long *total_bytes) {
    if (!context || !free_bytes || !total_bytes) {
        set_error("rslmto_reciprocal_cuda_get_memory: invalid arguments");
        return 1;
    }
    cudaSetDevice(context->device);
    std::size_t free_size = 0;
    std::size_t total_size = 0;
    const cudaError_t status = cudaMemGetInfo(&free_size, &total_size);
    if (status != cudaSuccess) {
        set_cuda_error("cudaMemGetInfo", status);
        *free_bytes = 0;
        *total_bytes = 0;
        return 1;
    }
    *free_bytes = static_cast<long long>(free_size);
    *total_bytes = static_cast<long long>(total_size);
    return 0;
}

extern "C" int rslmto_reciprocal_cuda_synchronize(
    rslmto_reciprocal_cuda_context *context) {
    if (!context) {
        set_error("rslmto_reciprocal_cuda_synchronize: null context");
        return 1;
    }
    const cudaError_t status = cudaStreamSynchronize(context->stream);
    if (status != cudaSuccess) {
        set_cuda_error("cudaStreamSynchronize", status);
        return 1;
    }
    return 0;
}

extern "C" void rslmto_reciprocal_cuda_destroy(
    rslmto_reciprocal_cuda_context *context) {
    if (!context) {
        return;
    }
    cudaSetDevice(context->device);
    release_solver_workspace(context);
    release_pinned_buffers(context);
    destroy_timing_events(context);
    if (context->jacobi_params) {
        cusolverDnDestroySyevjInfo(context->jacobi_params);
    }
    if (context->solver) {
        cusolverDnDestroy(context->solver);
    }
    if (context->stream) {
        cudaStreamDestroy(context->stream);
    }
    delete context;
}
