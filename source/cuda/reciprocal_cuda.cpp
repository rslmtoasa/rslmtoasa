#include "reciprocal_cuda.h"

#include <cuda_runtime_api.h>

#include <string>

struct rslmto_reciprocal_cuda_context {
    int device = -1;
    int prepared_operator_generation = -1;
    cudaStream_t stream = nullptr;
};

static std::string g_reciprocal_cuda_error;

static void set_error(const char *message) {
    g_reciprocal_cuda_error = message;
}

static void set_cuda_error(const char *operation, cudaError_t status) {
    g_reciprocal_cuda_error = std::string(operation) + ": " +
                              cudaGetErrorString(status);
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

    /* ACC-02 has no operator buffers yet.  Keep the generation in the
     * backend-owned context so ACC-03/04 can add device storage without
     * changing the Fortran lifecycle contract. */
    context->prepared_operator_generation = operator_generation;
    return 1;
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
    if (context->stream) {
        cudaStreamDestroy(context->stream);
    }
    delete context;
}
