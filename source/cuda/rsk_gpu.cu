/* rsk_gpu.cu -- CUDA k-space backend scaffold.
 *
 * Real kernels will be added in the dedicated implementation session:
 *   1. batched Bloch sums for ee/eeo/eecc
 *   2. batched GEMM for eeo(k)*ee(k)
 *   3. onsite addition and Hermiticity diagnostics
 *   4. batched Hermitian diagonalization
 */
