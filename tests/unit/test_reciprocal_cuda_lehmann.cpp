#include "reciprocal_cuda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::size_t complex_index(int row, int column, int rows) {
    return static_cast<std::size_t>(2) * (row + rows * column);
}

std::size_t complex_index_3d(int row, int column, int slice, int rows, int columns) {
    return static_cast<std::size_t>(2) * (row + rows * (column + columns * slice));
}

std::size_t complex_index_4d(int row, int column, int energy, int pair,
                             int nblk, int ne) {
    return static_cast<std::size_t>(2) *
           (row + nblk * (column + nblk * (energy + ne * pair)));
}

void set_complex(std::vector<double> &values, std::size_t index,
                 double real, double imag) {
    values[index] = real;
    values[index + 1] = imag;
}

void add_complex(std::vector<double> &values, std::size_t index,
                 double real, double imag) {
    values[index] += real;
    values[index + 1] += imag;
}

void cpu_reference(int nmat, int nk, int ne, int npair, int nblk,
                   const std::vector<double> &eigenvalues,
                   const std::vector<double> &eigenvectors,
                   const std::vector<double> &k_points,
                   const std::vector<double> &z_contour,
                   const std::vector<double> &dr,
                   const std::vector<int> &ioffset,
                   const std::vector<int> &joffset,
                   std::vector<double> &blocks) {
    std::fill(blocks.begin(), blocks.end(), 0.0);
    for (int pair = 0; pair < npair; ++pair) {
        for (int energy = 0; energy < ne; ++energy) {
            for (int ik = 0; ik < nk; ++ik) {
                const double angle = 2.0 * 3.14159265358979323846 *
                    (k_points[3 * ik] * dr[3 * pair] +
                     k_points[3 * ik + 1] * dr[3 * pair + 1] +
                     k_points[3 * ik + 2] * dr[3 * pair + 2]);
                const double phase_real = std::cos(angle);
                const double phase_imag = std::sin(angle);
                for (int band = 0; band < nmat; ++band) {
                    const std::size_t vi = complex_index_3d(ioffset[pair], band, ik, nmat, nmat);
                    const std::size_t vj = complex_index_3d(joffset[pair], band, ik, nmat, nmat);
                    const double vjr = eigenvectors[vj];
                    const double vji = -eigenvectors[vj + 1];
                    const double vr = eigenvectors[vi] * vjr - eigenvectors[vi + 1] * vji;
                    const double vi_imag = eigenvectors[vi] * vji + eigenvectors[vi + 1] * vjr;
                    const double numerator_real = phase_real * vr - phase_imag * vi_imag;
                    const double numerator_imag = phase_real * vi_imag + phase_imag * vr;
                    const double denominator_real = z_contour[2 * energy] - eigenvalues[band + nmat * ik];
                    const double denominator_imag = z_contour[2 * energy + 1];
                    const double denominator = denominator_real * denominator_real +
                                               denominator_imag * denominator_imag;
                    add_complex(blocks, complex_index_4d(0, 0, energy, pair, nblk, ne),
                                (numerator_real * denominator_real + numerator_imag * denominator_imag) / denominator / nk,
                                (numerator_imag * denominator_real - numerator_real * denominator_imag) / denominator / nk);
                    for (int column = 0; column < nblk; ++column) {
                        for (int row = 0; row < nblk; ++row) {
                            if (row == 0 && column == 0) continue;
                            const std::size_t vi_block = complex_index_3d(ioffset[pair] + row, band, ik, nmat, nmat);
                            const std::size_t vj_block = complex_index_3d(joffset[pair] + column, band, ik, nmat, nmat);
                            const double vj_block_real = eigenvectors[vj_block];
                            const double vj_block_imag = -eigenvectors[vj_block + 1];
                            const double outer_real = eigenvectors[vi_block] * vj_block_real -
                                                      eigenvectors[vi_block + 1] * vj_block_imag;
                            const double outer_imag = eigenvectors[vi_block] * vj_block_imag +
                                                      eigenvectors[vi_block + 1] * vj_block_real;
                            const double nr = phase_real * outer_real - phase_imag * outer_imag;
                            const double ni = phase_real * outer_imag + phase_imag * outer_real;
                            add_complex(blocks, complex_index_4d(row, column, energy, pair, nblk, ne),
                                        (nr * denominator_real + ni * denominator_imag) / denominator / nk,
                                        (ni * denominator_real - nr * denominator_imag) / denominator / nk);
                        }
                    }
                }
            }
        }
    }
}

} // namespace

int main() {
    int device_count = 0;
    if (rslmto_reciprocal_cuda_device_count(&device_count) != 0 || device_count == 0) {
        std::puts("SKIP: no CUDA device is available");
        return 77;
    }

    constexpr int nmat = 4;
    constexpr int nk = 5;
    constexpr int ne = 4;
    constexpr int npair = 3;
    constexpr int nblk = 2;
    std::vector<double> eigenvalues(nmat * nk);
    std::vector<double> eigenvectors(2 * nmat * nmat * nk);
    std::vector<double> k_points(3 * nk);
    std::vector<double> z_contour(2 * ne);
    std::vector<double> dr(3 * npair);
    std::vector<int> ioffset{0, 2, 1};
    std::vector<int> joffset{0, 0, 2};
    std::vector<double> expected(2 * nblk * nblk * ne * npair);
    std::vector<double> actual(expected.size());

    for (int ik = 0; ik < nk; ++ik) {
        k_points[3 * ik] = (ik + 0.25) / nk;
        k_points[3 * ik + 1] = (2 * ik + 0.5) / nk;
        k_points[3 * ik + 2] = (3 * ik + 0.75) / nk;
        for (int band = 0; band < nmat; ++band) {
            eigenvalues[band + nmat * ik] = -0.8 + 0.37 * band + 0.11 * ik;
            for (int row = 0; row < nmat; ++row) {
                const std::size_t index = complex_index_3d(row, band, ik, nmat, nmat);
                set_complex(eigenvectors, index,
                            0.07 * (row + 1) - 0.013 * (band + 1) + 0.005 * ik,
                            -0.03 * (row + 1) + 0.011 * (band + 1) - 0.004 * ik);
            }
        }
    }
    for (int energy = 0; energy < ne; ++energy) {
        set_complex(z_contour, 2 * energy, -1.1 + 0.41 * energy, 0.23 + 0.03 * energy);
    }
    dr = {0.0, 0.0, 0.0, -1.0, 0.25, 0.0, 0.5, -0.75, 1.0};
    cpu_reference(nmat, nk, ne, npair, nblk, eigenvalues, eigenvectors,
                  k_points, z_contour, dr, ioffset, joffset, expected);

    rslmto_reciprocal_cuda_context *context = rslmto_reciprocal_cuda_create(0);
    if (!context) {
        std::fprintf(stderr, "FAIL: CUDA context creation: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }
    double h2d = 0.0, contraction = 0.0, d2h = 0.0;
    const int status = rslmto_reciprocal_cuda_contract_lehmann(
        context, nmat, nk, ne, npair, nblk, eigenvalues.data(),
        eigenvectors.data(), k_points.data(), z_contour.data(), dr.data(),
        ioffset.data(), joffset.data(), actual.data(), &h2d, &contraction, &d2h);
    if (status != 0) {
        std::fprintf(stderr, "FAIL: CUDA Lehmann contraction: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }

    /* ACC-11: solve a complete k tile, consume the context-owned FP64
     * eigenpairs without presenting them to the Lehmann ABI, and prove that a
     * later values-only solve invalidates the handoff token. */
    std::vector<double> resident_h(2 * nmat * nmat * nk, 0.0);
    std::vector<double> resident_values(nmat * nk);
    std::vector<double> resident_vectors(resident_h.size());
    for (int ik = 0; ik < nk; ++ik) {
        for (int iband = 0; iband < nmat; ++iband) {
            set_complex(resident_h,
                        complex_index_3d(iband, iband, ik, nmat, nmat),
                        -0.7 + 0.31 * iband + 0.09 * ik, 0.0);
        }
    }
    const int resident_solve_status = rslmto_reciprocal_cuda_solve_zheevd_batch(
        context, nmat, nk, resident_h.data(), resident_values.data(),
        resident_vectors.data(), 1);
    if (resident_solve_status != 0) {
        std::fprintf(stderr, "FAIL: ACC-11 resident eigensolve: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }
    const int resident_token = rslmto_reciprocal_cuda_get_resident_token(context);
    if (resident_token <= 0 ||
        rslmto_reciprocal_cuda_resident_eigensystem_matches(
            context, nmat, nk, resident_token) != 0) {
        std::fprintf(stderr, "FAIL: ACC-11 resident token was not valid\n");
        return 1;
    }
    std::vector<double> resident_expected(expected.size());
    std::vector<double> resident_actual(expected.size());
    cpu_reference(nmat, nk, ne, npair, nblk, resident_values, resident_vectors,
                  k_points, z_contour, dr, ioffset, joffset, resident_expected);
    double resident_h2d = 0.0, resident_contraction = 0.0, resident_d2h = 0.0;
    const int resident_status = rslmto_reciprocal_cuda_contract_lehmann_resident(
        context, nmat, nk, ne, npair, nblk, resident_token, k_points.data(), z_contour.data(),
        dr.data(), ioffset.data(), joffset.data(), resident_actual.data(),
        &resident_h2d, &resident_contraction, &resident_d2h);
    if (resident_status != 0) {
        std::fprintf(stderr, "FAIL: ACC-11 resident Lehmann contraction: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }
    double resident_error = 0.0;
    for (std::size_t i = 0; i < resident_actual.size(); i += 2) {
        resident_error = std::max(
            resident_error,
            std::hypot(resident_actual[i] - resident_expected[i],
                       resident_actual[i + 1] - resident_expected[i + 1]));
    }
    if (resident_error > 3.e-11 || resident_h2d < 0.0 ||
        resident_contraction < 0.0 || resident_d2h < 0.0) {
        std::fprintf(stderr, "FAIL: ACC-11 resident error=%g timings=(%g,%g,%g)\n",
                     resident_error, resident_h2d, resident_contraction, resident_d2h);
        return 1;
    }
    std::vector<double> resident_values_only(nmat * nk);
    if (rslmto_reciprocal_cuda_solve_zheevd_batch(
            context, nmat, nk, resident_h.data(), resident_values_only.data(),
            nullptr, 0) != 0 ||
        rslmto_reciprocal_cuda_get_resident_token(context) != 0 ||
        rslmto_reciprocal_cuda_resident_eigensystem_matches(
            context, nmat, nk, resident_token) == 0) {
        std::fprintf(stderr, "FAIL: ACC-11 stale resident token was accepted\n");
        return 1;
    }

    rslmto_reciprocal_cuda_destroy(context);
    double error = 0.0;
    for (std::size_t i = 0; i < actual.size(); i += 2) {
        error = std::max(error, std::hypot(actual[i] - expected[i], actual[i + 1] - expected[i + 1]));
    }
    if (error > 3.e-11 || h2d < 0.0 || contraction < 0.0 || d2h < 0.0) {
        std::fprintf(stderr, "FAIL: max error=%g timings=(%g,%g,%g)\n",
                     error, h2d, contraction, d2h);
        return 1;
    }
    std::printf("PASS: CUDA Lehmann contraction max_error=%g resident_error=%g resident_h2d=%g resident_contraction=%g resident_d2h=%g\n",
                error, resident_error, resident_h2d, resident_contraction, resident_d2h);
    return 0;
}
