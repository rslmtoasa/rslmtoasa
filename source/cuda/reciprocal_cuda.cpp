#include "reciprocal_cuda.h"

#include <cuda_runtime_api.h>
#include <cusolverDn.h>

#include <algorithm>
#include <cstddef>
#include <string>

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

struct rslmto_reciprocal_cuda_context {
    int device = -1;
    int prepared_operator_generation = -1;
    cudaStream_t stream = nullptr;
    cusolverDnHandle_t solver = nullptr;
    int workspace_n = 0;
    int workspace_batch = 0;
    int workspace_lwork = 0;
    int workspace_jobz = -1;
    cuDoubleComplex *device_hamiltonians = nullptr;
    double *device_eigenvalues = nullptr;
    cuDoubleComplex *device_work = nullptr;
    int *device_info = nullptr;
    double h2d_seconds = 0.0;
    double solve_seconds = 0.0;
    double d2h_seconds = 0.0;
    int timing_calls = 0;
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

static void release_solver_workspace(rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    cudaFree(context->device_hamiltonians);
    cudaFree(context->device_eigenvalues);
    cudaFree(context->device_work);
    cudaFree(context->device_info);
    context->device_hamiltonians = nullptr;
    context->device_eigenvalues = nullptr;
    context->device_work = nullptr;
    context->device_info = nullptr;
    context->workspace_n = 0;
    context->workspace_batch = 0;
    context->workspace_lwork = 0;
    context->workspace_jobz = -1;
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
    if (n < 1 || batch_size < 1 || !host_hamiltonians || !host_eigenvalues ||
        (request_eigenvectors && !host_eigenvectors)) {
        set_error("rslmto_reciprocal_cuda_solve_zheevd_batch: invalid arguments");
        return 1;
    }

    const std::size_t matrix_elements = static_cast<std::size_t>(n) * n;
    const std::size_t matrix_bytes = matrix_elements * batch_size * sizeof(cuDoubleComplex);
    const std::size_t eigenvalue_bytes = static_cast<std::size_t>(n) * batch_size * sizeof(double);
    const int jobz = request_eigenvectors ? 1 : 0;
    if (n != context->workspace_n || batch_size > context->workspace_batch ||
        jobz != context->workspace_jobz) {
        release_solver_workspace(context);
        context->workspace_n = n;
        context->workspace_batch = batch_size;

        cudaError_t cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_hamiltonians), matrix_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_hamiltonians)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_eigenvalues), eigenvalue_bytes);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_eigenvalues)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_info), sizeof(int));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_info)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }

        int lwork = 0;
        cusolverStatus_t solver_status = cusolverDnZheevd_bufferSize(
            context->solver, request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER, n, context->device_hamiltonians, n,
            context->device_eigenvalues, &lwork);
        if (solver_status != CUSOLVER_STATUS_SUCCESS || lwork < 1) {
            set_solver_error("cusolverDnZheevd_bufferSize", solver_status);
            release_solver_workspace(context);
            return 1;
        }
        context->workspace_lwork = lwork;
        context->workspace_jobz = jobz;
        cuda_status = cudaMalloc(reinterpret_cast<void **>(&context->device_work),
                                 static_cast<std::size_t>(lwork) * sizeof(cuDoubleComplex));
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMalloc(device_work)", cuda_status);
            release_solver_workspace(context);
            return 1;
        }
    }

    cudaEvent_t h2d_start = nullptr, h2d_stop = nullptr;
    cudaEvent_t solve_start = nullptr, solve_stop = nullptr;
    cudaEvent_t d2h_start = nullptr, d2h_stop = nullptr;
    const bool timing_enabled =
        cudaEventCreate(&h2d_start) == cudaSuccess &&
        cudaEventCreate(&h2d_stop) == cudaSuccess &&
        cudaEventCreate(&solve_start) == cudaSuccess &&
        cudaEventCreate(&solve_stop) == cudaSuccess &&
        cudaEventCreate(&d2h_start) == cudaSuccess &&
        cudaEventCreate(&d2h_stop) == cudaSuccess;
    if (timing_enabled) cudaEventRecord(h2d_start, context->stream);

    cudaError_t cuda_status = cudaMemcpyAsync(
        context->device_hamiltonians, host_hamiltonians, matrix_bytes,
        cudaMemcpyHostToDevice, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(H2D Hamiltonians)", cuda_status);
        return 1;
    }
    if (timing_enabled) {
        cudaEventRecord(h2d_stop, context->stream);
        cudaEventRecord(solve_start, context->stream);
    }

    auto *device_matrices = context->device_hamiltonians;
    for (int ibatch = 0; ibatch < batch_size; ++ibatch) {
        cusolverStatus_t solver_status = cusolverDnZheevd(
            context->solver,
            request_eigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER, n, device_matrices + static_cast<std::size_t>(ibatch) * matrix_elements, n,
            context->device_eigenvalues + static_cast<std::size_t>(ibatch) * n,
            context->device_work, context->workspace_lwork, context->device_info);
        if (solver_status != CUSOLVER_STATUS_SUCCESS) {
            set_solver_error("cusolverDnZheevd", solver_status);
            return 1;
        }
    }
    if (timing_enabled) {
        cudaEventRecord(solve_stop, context->stream);
        cudaEventRecord(d2h_start, context->stream);
    }

    cuda_status = cudaMemcpyAsync(host_eigenvalues, context->device_eigenvalues,
                                  eigenvalue_bytes, cudaMemcpyDeviceToHost, context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpyAsync(D2H eigenvalues)", cuda_status);
        return 1;
    }
    if (request_eigenvectors) {
        cuda_status = cudaMemcpyAsync(host_eigenvectors, context->device_hamiltonians,
                                      matrix_bytes, cudaMemcpyDeviceToHost, context->stream);
        if (cuda_status != cudaSuccess) {
            set_cuda_error("cudaMemcpyAsync(D2H eigenvectors)", cuda_status);
            return 1;
        }
    }
    if (timing_enabled) cudaEventRecord(d2h_stop, context->stream);
    cuda_status = cudaStreamSynchronize(context->stream);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaStreamSynchronize(eigensolution)", cuda_status);
        return 1;
    }

    int host_info = 0;
    cuda_status = cudaMemcpy(&host_info, context->device_info, sizeof(host_info), cudaMemcpyDeviceToHost);
    if (cuda_status != cudaSuccess) {
        set_cuda_error("cudaMemcpy(D2H solver info)", cuda_status);
        return 1;
    }
    if (host_info != 0) {
        g_reciprocal_cuda_error = "cusolverDnZheevd: solver info " + std::to_string(host_info);
        return 1;
    }
    if (timing_enabled) {
        float h2d_ms = 0.0f, solve_ms = 0.0f, d2h_ms = 0.0f;
        cudaEventElapsedTime(&h2d_ms, h2d_start, h2d_stop);
        cudaEventElapsedTime(&solve_ms, solve_start, solve_stop);
        cudaEventElapsedTime(&d2h_ms, d2h_start, d2h_stop);
        context->h2d_seconds += static_cast<double>(h2d_ms) * 1.0e-3;
        context->solve_seconds += static_cast<double>(solve_ms) * 1.0e-3;
        context->d2h_seconds += static_cast<double>(d2h_ms) * 1.0e-3;
        ++context->timing_calls;
    }
    if (h2d_start) cudaEventDestroy(h2d_start);
    if (h2d_stop) cudaEventDestroy(h2d_stop);
    if (solve_start) cudaEventDestroy(solve_start);
    if (solve_stop) cudaEventDestroy(solve_stop);
    if (d2h_start) cudaEventDestroy(d2h_start);
    if (d2h_stop) cudaEventDestroy(d2h_stop);
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

extern "C" void rslmto_reciprocal_cuda_get_timings(
    rslmto_reciprocal_cuda_context *context,
    double *h2d_seconds, double *solve_seconds, double *d2h_seconds, int *calls) {
    if (!context) return;
    if (h2d_seconds) *h2d_seconds = context->h2d_seconds;
    if (solve_seconds) *solve_seconds = context->solve_seconds;
    if (d2h_seconds) *d2h_seconds = context->d2h_seconds;
    if (calls) *calls = context->timing_calls;
}

extern "C" void rslmto_reciprocal_cuda_reset_timings(
    rslmto_reciprocal_cuda_context *context) {
    if (!context) return;
    context->h2d_seconds = 0.0;
    context->solve_seconds = 0.0;
    context->d2h_seconds = 0.0;
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
    if (context->solver) {
        cusolverDnDestroy(context->solver);
    }
    if (context->stream) {
        cudaStreamDestroy(context->stream);
    }
    delete context;
}
