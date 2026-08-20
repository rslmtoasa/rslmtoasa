#include "rsrec_cuda.h"

#include <cuda_runtime_api.h>
#include <cusparse.h>
#include <cufft.h>

#include <string>
#include <vector>

struct rsrec_cuda_ctx {
    rsrec_ctx *inner = nullptr;
    int kk = 0, nb = 0, nnmax = 0, ntype = 0, nmax = 0, device = 0;
    int backend = RSREC_CUDA_BACKEND_CSR;
    bool have_h = false;
    bool have_v = false;
    bool pbc = false;
    int n1 = 0, n2 = 0, n3 = 0, nbas = 0;
    std::vector<double> a;
    std::vector<double> crd;
    cusparseHandle_t cusparse = nullptr;
    cufftHandle fft_plan = 0;
    bool fft_plan_ready = false;
};

static std::string g_cuda_err;
static void set_error(const std::string &msg) { g_cuda_err = msg; }

extern "C" const char *rsrec_cuda_last_error(void) {
    return g_cuda_err.c_str();
}

extern "C" int rsrec_cuda_device_count(int *count) {
    if (!count) {
        set_error("rsrec_cuda_device_count: null output");
        return 1;
    }
    const cudaError_t status = cudaGetDeviceCount(count);
    if (status != cudaSuccess) {
        *count = 0;
        set_error(std::string("rsrec_cuda_device_count: ") + cudaGetErrorString(status));
        return 1;
    }
    return 0;
}

static bool fft_eligible(const rsrec_cuda_ctx *ctx) {
    return ctx->pbc && ctx->nmax == 0 && ctx->nbas > 0 &&
           ctx->n1 > 0 && ctx->n2 > 0 && ctx->n3 > 0 &&
           ctx->a.size() == 9 && !ctx->crd.empty();
}

static int validate_backend(const rsrec_cuda_ctx *ctx) {
    if (ctx->backend == RSREC_CUDA_BACKEND_FFT ||
        ctx->backend == RSREC_CUDA_BACKEND_CONV) {
        if (!fft_eligible(ctx)) {
            set_error("Periodic FFT/conv backend requires lattice%pbc with "
                      "valid a/crd data and ee-only Hamiltonian (nmax=0).");
            return 1;
        }
    }
    return 0;
}

extern "C" rsrec_cuda_ctx *rsrec_cuda_create(int kk, int nb, int nnmax,
                                             int ntype, int nmax,
                                             int device) {
    if (kk <= 0 || nb <= 0 || nnmax <= 0 || ntype <= 0 || nmax < 0 ||
        nmax > kk) {
        set_error("rsrec_cuda_create: bad dimensions");
        return nullptr;
    }

    auto *ctx = new rsrec_cuda_ctx();
    ctx->kk = kk;
    ctx->nb = nb;
    ctx->nnmax = nnmax;
    ctx->ntype = ntype;
    ctx->nmax = nmax;
    ctx->device = device;

    if (cudaSetDevice(device) != cudaSuccess) {
        set_error("rsrec_cuda_create: cudaSetDevice failed");
        delete ctx;
        return nullptr;
    }

    if (cusparseCreate(&ctx->cusparse) != CUSPARSE_STATUS_SUCCESS) {
        set_error("rsrec_cuda_create: cusparseCreate failed");
        delete ctx;
        return nullptr;
    }

    ctx->inner = rsrec_create(kk, nb, nnmax, ntype, nmax, device);
    if (!ctx->inner) {
        set_error(std::string("rsrec_cuda_create: ") + rsrec_last_error());
        if (ctx->cusparse) cusparseDestroy(ctx->cusparse);
        delete ctx;
        return nullptr;
    }

    return ctx;
}

extern "C" void rsrec_cuda_destroy(rsrec_cuda_ctx *ctx) {
    if (!ctx) return;
    if (ctx->fft_plan_ready) cufftDestroy(ctx->fft_plan);
    if (ctx->inner) rsrec_destroy(ctx->inner);
    if (ctx->cusparse) cusparseDestroy(ctx->cusparse);
    delete ctx;
}

extern "C" int rsrec_cuda_set_backend(rsrec_cuda_ctx *ctx, int backend) {
    if (!ctx) {
        set_error("rsrec_cuda_set_backend: null ctx");
        return 1;
    }
    if (backend < RSREC_CUDA_BACKEND_CSR ||
        backend > RSREC_CUDA_BACKEND_CONV) {
        set_error("rsrec_cuda_set_backend: invalid backend");
        return 1;
    }
    ctx->backend = backend;
    return validate_backend(ctx);
}

extern "C" int rsrec_cuda_set_periodic_lattice(rsrec_cuda_ctx *ctx, int pbc,
                                               int n1, int n2, int n3,
                                               const double *a,
                                               const double *crd, int nbas) {
    if (!ctx) {
        set_error("rsrec_cuda_set_periodic_lattice: null ctx");
        return 1;
    }

    ctx->pbc = (pbc != 0);
    ctx->n1 = n1;
    ctx->n2 = n2;
    ctx->n3 = n3;
    ctx->nbas = nbas;

    ctx->a.assign(a, a + 9);
    if (crd && nbas > 0) {
        ctx->crd.assign(crd, crd + 3 * static_cast<size_t>(nbas));
    } else {
        ctx->crd.clear();
    }

    if (ctx->fft_plan_ready) {
        cufftDestroy(ctx->fft_plan);
        ctx->fft_plan_ready = false;
    }

    if (ctx->pbc && ctx->n1 > 0 && ctx->n2 > 0 && ctx->n3 > 0) {
        if (cufftPlan3d(&ctx->fft_plan, ctx->n1, ctx->n2, ctx->n3,
                        CUFFT_Z2Z) == CUFFT_SUCCESS) {
            ctx->fft_plan_ready = true;
        }
    }

    return validate_backend(ctx);
}

extern "C" int rsrec_cuda_set_hamiltonian(rsrec_cuda_ctx *ctx, const void *ee,
                                          const void *hall,
                                          const void *lsham, const int *nn,
                                          const int *iz, const void *eeo,
                                          const void *hallo, const void *enim) {
    if (!ctx) {
        set_error("rsrec_cuda_set_hamiltonian: null ctx");
        return 1;
    }
    const int status = rsrec_set_hamiltonian(ctx->inner, ee, hall, lsham, nn, iz,
                                             eeo, hallo, enim);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_set_hamiltonian: ") +
                  rsrec_last_error());
        return status;
    }
    ctx->have_h = true;
    return validate_backend(ctx);
}

extern "C" int rsrec_cuda_set_velocity(rsrec_cuda_ctx *ctx, const void *v_a,
                                       const void *v_b, const void *vo_a,
                                       const void *vo_b) {
    if (!ctx) {
        set_error("rsrec_cuda_set_velocity: null ctx");
        return 1;
    }
    const int status = rsrec_set_velocity(ctx->inner, v_a, v_b, vo_a, vo_b);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_set_velocity: ") + rsrec_last_error());
        return status;
    }
    ctx->have_v = true;
    return validate_backend(ctx);
}

extern "C" int rsrec_cuda_orbital_moments(rsrec_cuda_ctx *ctx, const void *left,
                                          const void *psiref, int lld, double a,
                                          double b, void *mu) {
    if (!ctx || !ctx->have_h) {
        set_error("rsrec_cuda_orbital_moments: Hamiltonian not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status =
        rsrec_orbital_moments(ctx->inner, left, psiref, lld, a, b, mu);
    if (status != 0)
        set_error(std::string("rsrec_cuda_orbital_moments: ") +
                  rsrec_last_error());
    return status;
}

extern "C" int rsrec_cuda_chebyshev_moments(rsrec_cuda_ctx *ctx,
                                            const void *psi0, int lld,
                                            double a, double b,
                                            void *mu_out) {
    if (!ctx || !ctx->have_h) {
        set_error("rsrec_cuda_chebyshev_moments: Hamiltonian not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status =
        rsrec_chebyshev_moments(ctx->inner, psi0, lld, a, b, mu_out);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_chebyshev_moments: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_block_lanczos(rsrec_cuda_ctx *ctx, const void *psi0,
                                        int lld, void *a_b, void *b2_b,
                                        int prec) {
    if (!ctx || !ctx->have_h) {
        set_error("rsrec_cuda_block_lanczos: Hamiltonian not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status =
        rsrec_block_lanczos(ctx->inner, psi0, lld, a_b, b2_b, prec);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_block_lanczos: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_scalar_lanczos(rsrec_cuda_ctx *ctx, int site_j,
                                         int lld, double *a_out,
                                         double *b2_out) {
    if (!ctx || !ctx->have_h) {
        set_error("rsrec_cuda_scalar_lanczos: Hamiltonian not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status =
        rsrec_scalar_lanczos(ctx->inner, site_j, lld, a_out, b2_out);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_scalar_lanczos: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_stochastic_moments(rsrec_cuda_ctx *ctx,
                                             const void *psiref, int lld,
                                             double a, double b,
                                             void *mu_nm) {
    if (!ctx || !ctx->have_h || !ctx->have_v) {
        set_error("rsrec_cuda_stochastic_moments: Hamiltonian/velocity not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status =
        rsrec_stochastic_moments(ctx->inner, psiref, lld, a, b, mu_nm);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_stochastic_moments: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_stochastic_moments_resident(
    rsrec_cuda_ctx *ctx, const void *psiref, int lld, double a, double b,
    int trace_index, int download_host, void *mu_nm) {
    if (!ctx || !ctx->have_h || !ctx->have_v) {
        set_error("rsrec_cuda_stochastic_moments_resident: Hamiltonian/velocity not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status = rsrec_stochastic_moments_resident(
        ctx->inner, psiref, lld, a, b, trace_index, download_host, mu_nm);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_stochastic_moments_resident: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_clear_resident_moments(rsrec_cuda_ctx *ctx) {
    if (!ctx) {
        set_error("rsrec_cuda_clear_resident_moments: null ctx");
        return 1;
    }
    const int status = rsrec_clear_resident_moments(ctx->inner);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_clear_resident_moments: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_resident_count(rsrec_cuda_ctx *ctx, int *count) {
    if (!ctx || !count) {
        set_error("rsrec_cuda_resident_count: null argument");
        return 1;
    }
    return rsrec_resident_count(ctx->inner, count);
}

extern "C" int rsrec_cuda_reconstruct_conductivity(
    rsrec_cuda_ctx *ctx, const double *energy, int n_energy, int moments,
    double a, double b, double factor, int ntrace, int energy_block,
    void *c_out, double *gamma_seconds, double *gamma_basis_seconds,
    double *gamma_fill_seconds, double *gemm_seconds,
    double *result_d2h_seconds, long long *gamma_h2d_bytes,
    long long *gamma_block_bytes, long long *result_d2h_bytes,
    int *actual_energy_block) {
    if (!ctx || !ctx->have_h) {
        set_error("rsrec_cuda_reconstruct_conductivity: Hamiltonian not set");
        return 1;
    }
    if (validate_backend(ctx) != 0) return 1;
    const int status = rsrec_reconstruct_conductivity(
        ctx->inner, energy, n_energy, moments, a, b, factor, ntrace,
        energy_block, c_out, gamma_seconds, gamma_basis_seconds,
        gamma_fill_seconds, gemm_seconds, result_d2h_seconds,
        gamma_h2d_bytes, gamma_block_bytes, result_d2h_bytes,
        actual_energy_block);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_reconstruct_conductivity: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_stochastic_profile(rsrec_cuda_ctx *ctx,
                                              double *h2d_seconds,
                                              double *cheb_seconds,
                                              double *d2h_seconds,
                                              double *conversion_seconds,
                                              double *mu_pack_seconds,
                                              long long *h2d_bytes,
                                              long long *d2h_bytes,
                                              long long *mu_pack_bytes) {
    if (!ctx) {
        set_error("rsrec_cuda_stochastic_profile: null ctx");
        return 1;
    }
    const int status = rsrec_stochastic_profile(ctx->inner, h2d_seconds,
                                                cheb_seconds, d2h_seconds,
                                                conversion_seconds,
                                                mu_pack_seconds, h2d_bytes,
                                                d2h_bytes, mu_pack_bytes);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_stochastic_profile: ") +
                  rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_set_precision(rsrec_cuda_ctx *ctx, int prec) {
    if (!ctx) {
        set_error("rsrec_cuda_set_precision: null ctx");
        return 1;
    }
    const int status = rsrec_set_precision(ctx->inner, prec);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_set_precision: ") + rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_chebyshev_dos(rsrec_cuda_ctx *ctx, const void *mu,
                                        int n_moments, int natoms,
                                        const double *ene, int nv, double a,
                                        double b, void *g0) {
    if (!ctx) {
        set_error("rsrec_cuda_chebyshev_dos: null ctx");
        return 1;
    }
    const int status =
        rsrec_chebyshev_dos(ctx->inner, mu, n_moments, natoms, ene, nv, a, b, g0);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_chebyshev_dos: ") + rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_chebyshev_gf_eta(rsrec_cuda_ctx *ctx,
                                           const void *mu, int n_moments,
                                           int natoms, const void *F,
                                           int n_eta, void *g0) {
    if (!ctx) {
        set_error("rsrec_cuda_chebyshev_gf_eta: null ctx");
        return 1;
    }
    const int status =
        rsrec_chebyshev_gf_eta(ctx->inner, mu, n_moments, natoms, F, n_eta, g0);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_chebyshev_gf_eta: ") + rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_block_dos(rsrec_cuda_ctx *ctx, const void *a_b,
                                    const void *b2_b, const double *a_inf,
                                    const double *b_inf, const double *ene,
                                    int nv, double eta_re, double eta_im,
                                    int natoms, int lld, int sym, void *g0) {
    if (!ctx) {
        set_error("rsrec_cuda_block_dos: null ctx");
        return 1;
    }
    const int status = rsrec_block_dos(ctx->inner, a_b, b2_b, a_inf, b_inf, ene,
                                       nv, eta_re, eta_im, natoms, lld, sym, g0);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_block_dos: ") + rsrec_last_error());
    }
    return status;
}

extern "C" int rsrec_cuda_block_gf_eta(rsrec_cuda_ctx *ctx, const void *a_b,
                                       const void *b2_b, const double *a_inf,
                                       const double *b_inf, double ef,
                                       const double *eta_re,
                                       const double *eta_im, int n_eta,
                                       int natoms, int lld, int sym, void *g0) {
    if (!ctx) {
        set_error("rsrec_cuda_block_gf_eta: null ctx");
        return 1;
    }
    const int status =
        rsrec_block_gf_eta(ctx->inner, a_b, b2_b, a_inf, b_inf, ef, eta_re,
                           eta_im, n_eta, natoms, lld, sym, g0);
    if (status != 0) {
        set_error(std::string("rsrec_cuda_block_gf_eta: ") + rsrec_last_error());
    }
    return status;
}
