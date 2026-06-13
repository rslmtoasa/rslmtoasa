#include <cstdint>
#include <vector>

#include "mkl.h"

extern "C" void rslmto_mkl_cgemm_batch_nn(
    const std::int32_t *batch_count_in,
    const std::int32_t *m_in,
    const std::int32_t *n_in,
    const std::int32_t *k_in,
    const MKL_Complex8 *alpha,
    const void **a_array,
    const std::int32_t *lda_in,
    const void **b_array,
    const std::int32_t *ldb_in,
    const MKL_Complex8 *beta,
    void **c_array,
    const std::int32_t *ldc_in,
    std::int32_t *status)
{
    if (!batch_count_in || !m_in || !n_in || !k_in || !alpha || !a_array ||
        !lda_in || !b_array || !ldb_in || !beta || !c_array || !ldc_in ||
        !status || *batch_count_in < 0) {
        if (status) {
            *status = -1;
        }
        return;
    }

    const MKL_INT group_count = 1;
    const MKL_INT group_size[1] = {static_cast<MKL_INT>(*batch_count_in)};
    const CBLAS_TRANSPOSE transa[1] = {CblasNoTrans};
    const CBLAS_TRANSPOSE transb[1] = {CblasNoTrans};
    const MKL_INT m[1] = {static_cast<MKL_INT>(*m_in)};
    const MKL_INT n[1] = {static_cast<MKL_INT>(*n_in)};
    const MKL_INT k[1] = {static_cast<MKL_INT>(*k_in)};
    const MKL_INT lda[1] = {static_cast<MKL_INT>(*lda_in)};
    const MKL_INT ldb[1] = {static_cast<MKL_INT>(*ldb_in)};
    const MKL_INT ldc[1] = {static_cast<MKL_INT>(*ldc_in)};
    const MKL_Complex8 alpha_array[1] = {*alpha};
    const MKL_Complex8 beta_array[1] = {*beta};

    cblas_cgemm_batch(CblasColMajor, transa, transb, m, n, k, alpha_array,
                      a_array, lda, b_array, ldb, beta_array,
                      c_array, ldc, group_count, group_size);
    *status = 0;
}
