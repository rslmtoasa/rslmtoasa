#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>

namespace {

__device__ inline double2 add_complex(double2 a, double2 b) {
    return make_double2(a.x + b.x, a.y + b.y);
}

__device__ inline double2 multiply_complex(double2 a, double2 b) {
    return make_double2(a.x * b.x - a.y * b.y,
                        a.x * b.y + a.y * b.x);
}

__device__ inline double2 conjugate_complex(double2 a) {
    return make_double2(a.x, -a.y);
}

__device__ inline double2 divide_by_complex(double2 numerator, double denominator_real,
                                             double denominator_imag) {
    const double denominator = denominator_real * denominator_real +
                               denominator_imag * denominator_imag;
    return make_double2(
        (numerator.x * denominator_real + numerator.y * denominator_imag) / denominator,
        (numerator.y * denominator_real - numerator.x * denominator_imag) / denominator);
}

__global__ void lehmann_contract_kernel(
    int nmat,
    int nk,
    int ne,
    int npair,
    int nblk,
    const double *eigenvalues,
    const double2 *eigenvectors,
    const double *k_points,
    const double2 *z_contour,
    const double *dr,
    const int *ioffset,
    const int *joffset,
    double2 *blocks) {
    const long long index = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long long total = static_cast<long long>(nblk) * nblk * ne * npair;
    if (index >= total) return;

    long long remainder = index;
    const int row = static_cast<int>(remainder % nblk);
    remainder /= nblk;
    const int column = static_cast<int>(remainder % nblk);
    remainder /= nblk;
    const int energy = static_cast<int>(remainder % ne);
    const int pair = static_cast<int>(remainder / ne);
    const int irow = ioffset[pair] + row;
    const int jrow = joffset[pair] + column;

    double2 accumulated = make_double2(0.0, 0.0);
    for (int ik = 0; ik < nk; ++ik) {
        const double kdotr = 2.0 * 3.1415926535897932384626433832795 *
            (k_points[3 * ik] * dr[3 * pair] +
             k_points[3 * ik + 1] * dr[3 * pair + 1] +
             k_points[3 * ik + 2] * dr[3 * pair + 2]);
        const double2 phase = make_double2(cos(kdotr), sin(kdotr));
        for (int band = 0; band < nmat; ++band) {
            const long long eigenvector_i = irow + static_cast<long long>(nmat) *
                (band + static_cast<long long>(nmat) * ik);
            const long long eigenvector_j = jrow + static_cast<long long>(nmat) *
                (band + static_cast<long long>(nmat) * ik);
            const double2 outer = multiply_complex(
                eigenvectors[eigenvector_i], conjugate_complex(eigenvectors[eigenvector_j]));
            const double2 numerator = multiply_complex(phase, outer);
            const double2 z = z_contour[energy];
            accumulated = add_complex(
                accumulated, divide_by_complex(numerator, z.x - eigenvalues[band + nmat * ik], z.y));
        }
    }
    const double normalization = 1.0 / static_cast<double>(nk);
    blocks[index] = make_double2(accumulated.x * normalization,
                                 accumulated.y * normalization);
}

template <typename T>
int allocate_and_copy(T **device, const T *host, std::size_t count, cudaStream_t stream) {
    cudaError_t status = cudaMalloc(reinterpret_cast<void **>(device), count * sizeof(T));
    if (status != cudaSuccess) return 1;
    status = cudaMemcpyAsync(*device, host, count * sizeof(T), cudaMemcpyHostToDevice, stream);
    return status == cudaSuccess ? 0 : 1;
}

} // namespace

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
    double *d2h_seconds) {
    double *device_eigenvalues = nullptr;
    double2 *device_eigenvectors = nullptr;
    double *device_k_points = nullptr;
    double2 *device_z_contour = nullptr;
    double *device_dr = nullptr;
    int *device_ioffset = nullptr;
    int *device_joffset = nullptr;
    double2 *device_blocks = nullptr;
    cudaEvent_t h2d_start = nullptr, h2d_stop = nullptr;
    cudaEvent_t contraction_start = nullptr, contraction_stop = nullptr;
    cudaEvent_t d2h_start = nullptr, d2h_stop = nullptr;
    int status = 1;

    const std::size_t eigenvalue_count = static_cast<std::size_t>(nmat) * nk;
    const std::size_t eigenvector_count = static_cast<std::size_t>(nmat) * nmat * nk;
    const std::size_t block_count = static_cast<std::size_t>(nblk) * nblk * ne * npair;

    do {
        if (cudaEventCreate(&h2d_start) != cudaSuccess ||
            cudaEventCreate(&h2d_stop) != cudaSuccess ||
            cudaEventCreate(&contraction_start) != cudaSuccess ||
            cudaEventCreate(&contraction_stop) != cudaSuccess ||
            cudaEventCreate(&d2h_start) != cudaSuccess ||
            cudaEventCreate(&d2h_stop) != cudaSuccess) break;
        if (cudaEventRecord(h2d_start, stream) != cudaSuccess) break;
        if (allocate_and_copy(&device_eigenvalues, host_eigenvalues, eigenvalue_count, stream) != 0) break;
        if (allocate_and_copy(&device_eigenvectors,
                              static_cast<const double2 *>(host_eigenvectors), eigenvector_count, stream) != 0) break;
        if (allocate_and_copy(&device_k_points, host_k_points, static_cast<std::size_t>(3) * nk, stream) != 0) break;
        if (allocate_and_copy(&device_z_contour,
                              static_cast<const double2 *>(host_z_contour), ne, stream) != 0) break;
        if (allocate_and_copy(&device_dr, host_dr, static_cast<std::size_t>(3) * npair, stream) != 0) break;
        if (allocate_and_copy(&device_ioffset, host_ioffset, npair, stream) != 0) break;
        if (allocate_and_copy(&device_joffset, host_joffset, npair, stream) != 0) break;
        if (cudaMalloc(reinterpret_cast<void **>(&device_blocks), block_count * sizeof(double2)) != cudaSuccess) break;
        if (cudaEventRecord(h2d_stop, stream) != cudaSuccess) break;
        if (cudaEventRecord(contraction_start, stream) != cudaSuccess) break;

        const long long total = static_cast<long long>(block_count);
        const unsigned int grid = static_cast<unsigned int>((total + 255) / 256);
        lehmann_contract_kernel<<<grid, 256, 0, stream>>>(
            nmat, nk, ne, npair, nblk, device_eigenvalues, device_eigenvectors,
            device_k_points, device_z_contour, device_dr, device_ioffset,
            device_joffset, device_blocks);
        if (cudaGetLastError() != cudaSuccess) break;
        if (cudaEventRecord(contraction_stop, stream) != cudaSuccess) break;
        if (cudaEventRecord(d2h_start, stream) != cudaSuccess) break;
        if (cudaMemcpyAsync(host_blocks, device_blocks, block_count * sizeof(double2),
                            cudaMemcpyDeviceToHost, stream) != cudaSuccess) break;
        if (cudaEventRecord(d2h_stop, stream) != cudaSuccess) break;
        if (cudaStreamSynchronize(stream) != cudaSuccess) break;

        float h2d_ms = 0.0f, contraction_ms = 0.0f, d2h_ms = 0.0f;
        if (cudaEventElapsedTime(&h2d_ms, h2d_start, h2d_stop) != cudaSuccess ||
            cudaEventElapsedTime(&contraction_ms, contraction_start, contraction_stop) != cudaSuccess ||
            cudaEventElapsedTime(&d2h_ms, d2h_start, d2h_stop) != cudaSuccess) break;
        if (h2d_seconds) *h2d_seconds = static_cast<double>(h2d_ms) * 1.e-3;
        if (contraction_seconds) *contraction_seconds = static_cast<double>(contraction_ms) * 1.e-3;
        if (d2h_seconds) *d2h_seconds = static_cast<double>(d2h_ms) * 1.e-3;
        status = 0;
    } while (false);

    cudaFree(device_eigenvalues);
    cudaFree(device_eigenvectors);
    cudaFree(device_k_points);
    cudaFree(device_z_contour);
    cudaFree(device_dr);
    cudaFree(device_ioffset);
    cudaFree(device_joffset);
    cudaFree(device_blocks);
    cudaEventDestroy(h2d_start);
    cudaEventDestroy(h2d_stop);
    cudaEventDestroy(contraction_start);
    cudaEventDestroy(contraction_stop);
    cudaEventDestroy(d2h_start);
    cudaEventDestroy(d2h_stop);
    return status;
}
