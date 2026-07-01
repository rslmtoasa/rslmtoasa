#ifndef RSK_CUDA_H
#define RSK_CUDA_H

#include "rsk.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rsk_cuda_ctx rsk_cuda_ctx;

rsk_cuda_ctx *rsk_cuda_create(int nb, int nsite, int ntype, int nnmax,
                              int device);
void rsk_cuda_destroy(rsk_cuda_ctx *ctx);
const char *rsk_cuda_last_error(void);

#ifdef __cplusplus
}
#endif

#endif /* RSK_CUDA_H */
