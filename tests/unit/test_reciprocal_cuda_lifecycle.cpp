#include "reciprocal_cuda.h"

#include <cstdio>

int main() {
    int device_count = 0;
    if (rslmto_reciprocal_cuda_device_count(&device_count) != 0 ||
        device_count == 0) {
        std::puts("SKIP: no CUDA device is available");
        return 77;
    }

    rslmto_reciprocal_cuda_context *context =
        rslmto_reciprocal_cuda_create(0);
    if (!context) {
        std::fprintf(stderr, "FAIL: CUDA context creation: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }

    const int first_prepare =
        rslmto_reciprocal_cuda_prepare_operator(context, 17);
    const int repeated_prepare =
        rslmto_reciprocal_cuda_prepare_operator(context, 17);
    const int changed_prepare =
        rslmto_reciprocal_cuda_prepare_operator(context, 18);
    const int sync_status = rslmto_reciprocal_cuda_synchronize(context);
    rslmto_reciprocal_cuda_destroy(context);

    if (first_prepare != 1 || repeated_prepare != 0 || changed_prepare != 1 ||
        sync_status != 0) {
        std::fprintf(stderr,
                     "FAIL: generation-aware lifecycle (%d, %d, %d, %d)\n",
                     first_prepare, repeated_prepare, changed_prepare,
                     sync_status);
        return 1;
    }

    std::puts("PASS: reciprocal CUDA lifecycle and generation reuse");
    return 0;
}
