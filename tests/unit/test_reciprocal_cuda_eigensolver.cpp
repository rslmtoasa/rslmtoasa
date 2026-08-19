#include "reciprocal_cuda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

double &component(std::vector<double> &values, int n, int batch, int row,
                  int column, int part) {
    const std::size_t scalar = static_cast<std::size_t>(row) +
                               static_cast<std::size_t>(n) *
                                   (static_cast<std::size_t>(column) +
                                    static_cast<std::size_t>(n) * batch);
    return values[2 * scalar + part];
}

double component(const std::vector<double> &values, int n, int batch, int row,
                 int column, int part) {
    const std::size_t scalar = static_cast<std::size_t>(row) +
                               static_cast<std::size_t>(n) *
                                   (static_cast<std::size_t>(column) +
                                    static_cast<std::size_t>(n) * batch);
    return values[2 * scalar + part];
}

void set_real(std::vector<double> &matrix, int n, int batch, int row,
              int column, double value) {
    component(matrix, n, batch, row, column, 0) = value;
    component(matrix, n, batch, row, column, 1) = 0.0;
}

bool check_batch(const std::vector<double> &h, const std::vector<double> &w,
                 const std::vector<double> &v, int n, int batch,
                 double tolerance) {
    double worst_residual = 0.0;
    double worst_orthogonality = 0.0;
    for (int ibatch = 0; ibatch < batch; ++ibatch) {
        for (int column = 0; column < n; ++column) {
            double norm = 0.0;
            for (int row = 0; row < n; ++row) {
                double real = 0.0;
                double imag = 0.0;
                for (int inner = 0; inner < n; ++inner) {
                    const double ar = component(h, n, ibatch, row, inner, 0);
                    const double ai = component(h, n, ibatch, row, inner, 1);
                    const double vr = component(v, n, ibatch, inner, column, 0);
                    const double vi = component(v, n, ibatch, inner, column, 1);
                    real += ar * vr - ai * vi;
                    imag += ar * vi + ai * vr;
                }
                real -= w[static_cast<std::size_t>(column) +
                          static_cast<std::size_t>(n) * ibatch] *
                        component(v, n, ibatch, row, column, 0);
                imag -= w[static_cast<std::size_t>(column) +
                          static_cast<std::size_t>(n) * ibatch] *
                        component(v, n, ibatch, row, column, 1);
                norm += real * real + imag * imag;
            }
            worst_residual = std::max(worst_residual, std::sqrt(norm));
        }
        for (int left = 0; left < n; ++left) {
            for (int right = 0; right < n; ++right) {
                double real = (left == right) ? -1.0 : 0.0;
                double imag = 0.0;
                for (int row = 0; row < n; ++row) {
                    const double lr = component(v, n, ibatch, row, left, 0);
                    const double li = component(v, n, ibatch, row, left, 1);
                    const double rr = component(v, n, ibatch, row, right, 0);
                    const double ri = component(v, n, ibatch, row, right, 1);
                    real += lr * rr + li * ri;
                    imag += lr * ri - li * rr;
                }
                worst_orthogonality = std::max(worst_orthogonality,
                                               std::hypot(real, imag));
            }
        }
    }
    if (worst_residual > tolerance || worst_orthogonality > tolerance) {
        std::fprintf(stderr, "FAIL: residual=%g orthogonality=%g\n",
                     worst_residual, worst_orthogonality);
        return false;
    }
    return true;
}

bool check_degenerate_projector(const std::vector<double> &w,
                                const std::vector<double> &v, int n) {
    double error = 0.0;
    int selected = 0;
    for (int column = 0; column < n; ++column) {
        if (std::abs(w[column] - 2.0) > 1.e-10) continue;
        ++selected;
    }
    if (selected != 2) {
        std::fprintf(stderr, "FAIL: degenerate projector selected=%d\n", selected);
        return false;
    }
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            const double expected = (row == column && row < 2) ? 1.0 : 0.0;
            double real = 0.0;
            double imag = 0.0;
            for (int eigenvector = 0; eigenvector < n; ++eigenvector) {
                if (std::abs(w[eigenvector] - 2.0) > 1.e-10) continue;
                const double vr = component(v, n, 0, row, eigenvector, 0);
                const double vi = component(v, n, 0, row, eigenvector, 1);
                const double cr = component(v, n, 0, column, eigenvector, 0);
                const double ci = component(v, n, 0, column, eigenvector, 1);
                real += vr * cr + vi * ci;
                imag += vi * cr - vr * ci;
            }
            error = std::max(error, std::hypot(real - expected, imag));
        }
    }
    if (error > 1.e-10) {
        std::fprintf(stderr, "FAIL: degenerate projector error=%g\n", error);
        return false;
    }
    return true;
}

bool solve_and_check(rslmto_reciprocal_cuda_context *context, int n, int batch,
                     std::vector<double> &h) {
    std::vector<double> w(static_cast<std::size_t>(n) * batch);
    std::vector<double> v(h.size());
    if (rslmto_reciprocal_cuda_solve_zheevd_batch(
            context, n, batch, h.data(), w.data(), v.data(), 1) != 0) {
        std::fprintf(stderr, "FAIL: eigensolver: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return false;
    }
    if (n == 2) {
        const double gap = std::sqrt(17.0);
        for (int ibatch = 0; ibatch < batch; ++ibatch) {
            const double expected_low = (7.0 + 2.0 * ibatch - gap) / 2.0;
            const double expected_high = (7.0 + 2.0 * ibatch + gap) / 2.0;
            const double low = w[static_cast<std::size_t>(2) * ibatch];
            const double high = w[static_cast<std::size_t>(2) * ibatch + 1];
            if (std::abs(low - expected_low) > 5.e-11 ||
                std::abs(high - expected_high) > 5.e-11) {
                std::fprintf(stderr, "FAIL: eigenvalue batch=%d got=(%.17g,%.17g)\n",
                             ibatch, low, high);
                return false;
            }
        }
    }
    if (!check_batch(h, w, v, n, batch, 5.e-11)) return false;
    if (n == 4 && !check_degenerate_projector(w, v, n)) return false;
    return true;
}

bool solve_values_only_and_check(rslmto_reciprocal_cuda_context *context,
                                 const std::vector<double> &h, int n, int batch) {
    std::vector<double> w(static_cast<std::size_t>(n) * batch);
    if (rslmto_reciprocal_cuda_solve_zheevd_batch(
            context, n, batch, h.data(), w.data(), nullptr, 0) != 0) {
        std::fprintf(stderr, "FAIL: eigenvalues-only solve: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return false;
    }
    const double gap = std::sqrt(17.0);
    for (int ibatch = 0; ibatch < batch; ++ibatch) {
        const double low = (7.0 + 2.0 * ibatch - gap) / 2.0;
        const double high = (7.0 + 2.0 * ibatch + gap) / 2.0;
        if (std::abs(w[static_cast<std::size_t>(2) * ibatch] - low) > 5.e-11 ||
            std::abs(w[static_cast<std::size_t>(2) * ibatch + 1] - high) > 5.e-11) {
            std::fprintf(stderr, "FAIL: eigenvalues-only result batch=%d\n", ibatch);
            return false;
        }
    }
    return true;
}

bool solve_fp32_and_check(rslmto_reciprocal_cuda_context *context, int n,
                          int batch, const std::vector<double> &h) {
    std::vector<double> w(static_cast<std::size_t>(n) * batch);
    std::vector<double> v(h.size());
    if (rslmto_reciprocal_cuda_solve_cheevd_batch(
            context, n, batch, h.data(), w.data(), v.data(), 1) != 0) {
        std::fprintf(stderr, "FAIL: FP32 eigensolver: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return false;
    }
    if (n == 2) {
        const double gap = std::sqrt(17.0);
        for (int ibatch = 0; ibatch < batch; ++ibatch) {
            const double expected_low = (7.0 + 2.0 * ibatch - gap) / 2.0;
            const double expected_high = (7.0 + 2.0 * ibatch + gap) / 2.0;
            if (std::abs(w[static_cast<std::size_t>(2) * ibatch] - expected_low) > 2.e-5 ||
                std::abs(w[static_cast<std::size_t>(2) * ibatch + 1] - expected_high) > 2.e-5) {
                std::fprintf(stderr, "FAIL: FP32 eigenvalue batch=%d\n", ibatch);
                return false;
            }
        }
    }
    return check_batch(h, w, v, n, batch, 2.e-5);
}

bool select_strategy(rslmto_reciprocal_cuda_context *context, int strategy) {
    if (rslmto_reciprocal_cuda_set_solver_strategy(context, strategy) != 0) {
        std::fprintf(stderr, "FAIL: selecting CUDA solver strategy %d: %s\n",
                     strategy, rslmto_reciprocal_cuda_last_error());
        return false;
    }
    return true;
}

} // namespace

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
        std::fprintf(stderr, "FAIL: context creation: %s\n",
                     rslmto_reciprocal_cuda_last_error());
        return 1;
    }

    std::vector<double> h2(2 * 2 * 2 * 3, 0.0);
    for (int ibatch = 0; ibatch < 3; ++ibatch) {
        set_real(h2, 2, ibatch, 0, 0, 2.0 + ibatch);
        set_real(h2, 2, ibatch, 1, 1, 5.0 + ibatch);
        component(h2, 2, ibatch, 0, 1, 0) = 1.0;
        component(h2, 2, ibatch, 0, 1, 1) = 1.0;
        component(h2, 2, ibatch, 1, 0, 0) = 1.0;
        component(h2, 2, ibatch, 1, 0, 1) = -1.0;
    }
    bool passed = solve_and_check(context, 2, 3, h2);
    passed = solve_values_only_and_check(context, h2, 2, 3) && passed;

    std::vector<double> h4(2 * 4 * 4, 0.0);
    set_real(h4, 4, 0, 0, 0, 2.0);
    set_real(h4, 4, 0, 1, 1, 2.0);
    set_real(h4, 4, 0, 2, 2, 5.0);
    set_real(h4, 4, 0, 3, 3, 7.0);
    passed = solve_and_check(context, 4, 1, h4) && passed;

    std::vector<double> h8(2 * 8 * 8 * 2, 0.0);
    for (int ibatch = 0; ibatch < 2; ++ibatch) {
        for (int i = 0; i < 8; ++i) set_real(h8, 8, ibatch, i, i, -2.0 + i + ibatch);
    }
    passed = solve_and_check(context, 8, 2, h8) && passed;

    std::vector<double> h18(2 * 18 * 18, 0.0);
    for (int i = 0; i < 18; ++i) set_real(h18, 18, 0, i, i, -3.0 + 0.25 * i);
    passed = solve_and_check(context, 18, 1, h18) && passed;

    std::vector<double> h36(2 * 36 * 36, 0.0);
    for (int i = 0; i < 36; ++i) set_real(h36, 36, 0, i, i, -4.0 + 0.125 * i);
    passed = solve_and_check(context, 36, 1, h36) && passed;

    /* ACC-P1: exercise one true same-size batch call for every supported
     * small real/complex layout used by the reciprocal path. */
    passed = select_strategy(context, RSLMTO_RECIPROCAL_CUDA_ZHEEVJ_BATCHED) && passed;
    passed = solve_and_check(context, 2, 3, h2) && passed;
    passed = solve_values_only_and_check(context, h2, 2, 3) && passed;
    passed = solve_and_check(context, 4, 1, h4) && passed;
    passed = solve_and_check(context, 8, 2, h8) && passed;
    passed = solve_and_check(context, 18, 1, h18) && passed;
    if (rslmto_reciprocal_cuda_solver_strategy_supported(context, 36, 1, 1) != 1) {
        std::fprintf(stderr, "FAIL: ZheevjBatched n=36 was not reported unsupported\n");
        passed = false;
    }
    std::vector<double> unsupported_w(36, 0.0);
    std::vector<double> unsupported_v(h36.size(), 0.0);
    if (rslmto_reciprocal_cuda_solve_zheevd_batch(
            context, 36, 1, h36.data(), unsupported_w.data(), unsupported_v.data(), 1) != 2) {
        std::fprintf(stderr, "FAIL: unsupported ZheevjBatched request did not fail explicitly\n");
        passed = false;
    }
    /* Restore the reference strategy explicitly; no automatic threshold or
     * hidden CPU fallback is permitted by ACC-P1. */
    passed = select_strategy(context, RSLMTO_RECIPROCAL_CUDA_ZHEEVD_SERIAL) && passed;

    /* ACC-P1b: explicit FP32 direct and small-batch routes use the same H64
     * input boundary and widen eigenpairs back to the C ABI. */
    passed = select_strategy(context, RSLMTO_RECIPROCAL_CUDA_CHEEVD_SERIAL) && passed;
    passed = solve_fp32_and_check(context, 2, 3, h2) && passed;
    passed = solve_fp32_and_check(context, 18, 1, h18) && passed;
    passed = select_strategy(context, RSLMTO_RECIPROCAL_CUDA_CHEEVJ_BATCHED) && passed;
    passed = solve_fp32_and_check(context, 2, 3, h2) && passed;
    passed = solve_fp32_and_check(context, 18, 1, h18) && passed;
    if (rslmto_reciprocal_cuda_solver_strategy_supported(context, 36, 1, 1) != 1) {
        std::fprintf(stderr, "FAIL: CheevjBatched n=36 was not reported unsupported\n");
        passed = false;
    }

    rslmto_reciprocal_cuda_destroy(context);
    if (!passed) return 1;
    std::puts("PASS: standard complex-Hermitian CUDA eigensolver");
    return 0;
}
