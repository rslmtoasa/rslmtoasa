/*
 * Narrow C ABI for the reciprocal CUDA execution backend.
 *
 * ACC-02 owns context/lifecycle plumbing only.  Dense eigensolution is
 * intentionally not part of this ABI until the cuSOLVER API is selected in
 * ACC-03.
 */
#ifndef RSLMTO_RECIPROCAL_CUDA_H
#define RSLMTO_RECIPROCAL_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rslmto_reciprocal_cuda_context rslmto_reciprocal_cuda_context;

int rslmto_reciprocal_cuda_device_count(int *count);
rslmto_reciprocal_cuda_context *rslmto_reciprocal_cuda_create(int device);
const char *rslmto_reciprocal_cuda_last_error(void);

int rslmto_reciprocal_cuda_prepare_operator(
    rslmto_reciprocal_cuda_context *context, int operator_generation);
int rslmto_reciprocal_cuda_synchronize(
    rslmto_reciprocal_cuda_context *context);
void rslmto_reciprocal_cuda_destroy(
    rslmto_reciprocal_cuda_context *context);

#ifdef __cplusplus
}
#endif

#endif /* RSLMTO_RECIPROCAL_CUDA_H */
