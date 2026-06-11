/* ===========================================================================
 * rsrec.h -- GPU/CPU recursion engine for RS-LMTO-ASA (recursion_mod)
 *
 * Drop-in compute backends for:
 *   recur()                    -> rsrec_scalar_lanczos   (batched over nb orbitals)
 *   recur_b(), recur_b_ij()    -> rsrec_block_lanczos    (any psi0: single site
 *                                                         or the 4 ij sign combos)
 *   chebyshev_recur(),
 *   chebyshev_recur_ij()       -> rsrec_chebyshev_moments
 *   compute_moments_stochastic -> rsrec_stochastic_moments
 *
 * Hamiltonian model (mirrors hamiltonian_mod exactly):
 *   H_{k,k}       = hall(:,:,1,k) + lsham(:,:,iz(k))      for k <= nmax
 *                 = ee(:,:,1,iz(k)) + lsham(:,:,iz(k))    for k >  nmax
 *   H_{k,nn(k,j)} = hall(:,:,j,k)                          for k <= nmax
 *                 = ee(:,:,j,iz(k))                        for k >  nmax
 * with j = 2..nn(k,1). nmax may be 0 (pure bulk) or kk (hall covers everything).
 *
 * All complex arrays are double complex in FORTRAN COLUMN-MAJOR order with
 * the exact shapes used in recursion_mod; pass them directly from Fortran.
 * Atom indices in nn(:,:) and site arguments are 1-BASED (Fortran style).
 *
 * The izero/idum "lightcone" bookkeeping of the CPU code is intentionally
 * dropped: atoms outside the lightcone carry exactly zero wavefunction, so
 * a full sweep gives bit-identical results; on the GPU the dense sweep is
 * faster than divergent masking for typical cluster/recursion sizes.
 *
 * Thread safety: one context per thread/MPI rank. All routines synchronous.
 * Returns: 0 on success, nonzero on error (message via rsrec_last_error).
 * =========================================================================== */
#ifndef RSREC_H
#define RSREC_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rsrec_ctx rsrec_ctx;

/* --- lifecycle ----------------------------------------------------------- */
/* kk: cluster size, nb: block size (18), nnmax: 2nd dim of nn(:,:),
 * ntype: number of atom types, nmax: number of "impurity" atoms described
 * by hall (0 => bulk only), device: CUDA device id (ignored by CPU build). */
rsrec_ctx *rsrec_create(int kk, int nb, int nnmax, int ntype, int nmax,
                        int device);
void rsrec_destroy(rsrec_ctx *ctx);
const char *rsrec_last_error(void);

/* --- data upload (call once per SCF iteration / Hamiltonian change) ------ */
/* ee   : (nb, nb, nnmax, ntype)   complex(rp)
 * hall : (nb, nb, nnmax, nmax)    complex(rp)   (may be NULL if nmax == 0)
 * lsham: (nb, nb, ntype)          complex(rp)   (may be NULL => zero)
 * nn   : (kk, nnmax)              integer, nn(k,1) = #neighbour entries,
 *                                 nn(k,j) = 1-based neighbour atom (0 = none)
 * iz   : (kk)                     integer, 1-based atom type                 */
int rsrec_set_hamiltonian(rsrec_ctx *ctx, const void *ee, const void *hall,
                          const void *lsham, const int *nn, const int *iz);

/* Velocity (or current) operators for the stochastic conductivity moments.
 * Same neighbour structure as ee, per-type, no lsham term:
 * v_a, v_b: (nb, nb, nnmax, ntype) complex(rp). */
int rsrec_set_velocity(rsrec_ctx *ctx, const void *v_a, const void *v_b);

/* --- structured (FFT stencil + correction) acceleration ------------------ */
/* Map the clust onto a regular lattice grid and switch the matvec to
 *      H = Stencil(ee[:, :, :, t_ref])  +  dH,   dH = H_true - H_stencil
 * The stencil is evaluated as a zero-padded linear convolution (cuFFT on
 * the GPU): vacuum cells carry psi = 0, so OPEN BOUNDARIES ARE AUTOMATIC.
 * dH collects every deviating row -- the hall impurity region, atoms of
 * type /= t_ref, irregular neighbour patterns -- and is applied with the
 * block-sparse batched-GEMM machinery. Bit-compatible with the ELL backend
 * (validated to machine precision); profitable when most rows are bulk.
 *
 *   coords        : (3, kk) INTEGER lattice coordinates of each atom
 *                   (any consistent integer cell indexing; one atom per
 *                   cell -- multi-atom bases are not supported yet)
 *   use_structured: 1 = enable, 0 = revert to the block-ELL backend
 *
 * Call AFTER rsrec_set_hamiltonian (and after rsrec_set_velocity if the
 * stochastic driver will be used). Requires nmax < kk (a bulk must exist).
 * Re-call after each Hamiltonian upload to refresh the stencil tables.     */
int rsrec_set_grid(rsrec_ctx *ctx, const int *coords, int use_structured);

/* --- Chebyshev arithmetic precision --------------------------------------
 * prec = 0 (DEFAULT): the block-Chebyshev driver runs in COMPLEX SINGLE
 *          precision on a device-resident fused engine. KPM is robust in
 *          fp32 (iterates bounded on [-1,1], Jackson damping suppresses
 *          high-order noise; cf. Weisse et al. RMP 78, 275; KITE, GPUQT);
 *          validated here to ~1e-6 relative on the moments. On FP64-capped
 *          GPUs (GeForce/RTX A-series: FP64 = FP32/64) this is the only
 *          route to large speedups.
 * prec = 1: full double precision, bit-comparable with the CPU reference.
 *          Recommended once per new system as a cross-check.
 * Affects rsrec_chebyshev_moments only; the Lanczos drivers (whose
 * orthogonalisation is precision-sensitive) always run in fp64.            */
int rsrec_set_precision(rsrec_ctx *ctx, int prec);

/* --- Chebyshev Green-function / DOS reconstruction ------------------------
 * GPU/GEMM port of chebyshev_green (bands.f90): for each local atom n
 *   g0(:,:,ie,n) = sum_i mu(:,:,i,n) * gJ(i) c_i (-i) e^{-i(i-1)acos(x)}
 *                  / sqrt(a^2 - (E_ie - b)^2),   x = (E_ie - b)/a
 * with the Jackson kernel gJ in the exact bands.f90 convention and
 * c_1 = 1, c_i = 2 otherwise. The transfer matrix F(i,ie) is built once;
 * the moment contraction is one strided-batched GEMM over atoms.
 *
 *   mu      : (nb, nb, n_moments, natoms) complex(rp) -- the LOCAL slice
 *             mu_n(:,:,:,g2l_map(start_atom:end_atom)); n_moments = 2*lld+2
 *   ene     : (nv) real(rp) energy grid; keep it inside the (a,b) window
 *             (outside points give NaN, exactly like the CPU loop)
 *   g0      : (nb, nb, nv, natoms) complex(rp) out -- same layout as
 *             this%g0(:,:,1:nv, local atoms); DOS = -aimag(g0)/pi
 * a, b are the window coefficients, a = (emax-emin)/(2-0.3),
 * b = (emax+emin)/2, as in chebyshev_green. Precision follows
 * rsrec_set_precision (fp32 default ~1e-6, matching the moment engine;
 * fp64 bit-comparable with the Fortran triple loop).                       */
int rsrec_chebyshev_dos(rsrec_ctx *ctx, const void *mu, int n_moments,
                        int natoms, const double *ene, int nv,
                        double a, double b, void *g0);

/* --- core matvec (exposed mostly for testing / custom drivers) ----------- */
/* y = (H x - b x)/a   with x, y: (nb, nrhs, kk).  a=1,b=0 gives plain H x.
 * which: 0 = Hamiltonian, 1 = v_a, 2 = v_b (no shift/scale applied to 1,2
 * in the drivers; available here for completeness).                         */
int rsrec_op_apply(rsrec_ctx *ctx, int which, const void *x, void *y,
                   int nrhs, double a, double b);

/* --- drivers -------------------------------------------------------------- */

/* Chebyshev block moments, mirrors cheb_0th_mom/cheb_1st_mom/
 * chebyshev_recur_ll for a general starting block state.
 *   psi0  : (nb, nb, kk)  starting state (identity block on one site for
 *           chebyshev_recur; the +/-,i sign combos on sites i,j for _ij)
 *   lld   : number of recursion steps (control%lld)
 *   a, b  : Chebyshev scale and shift (resolve_chebyshev_window)
 *   mu_out: (nb, nb, 2*lld+2)  block moments, identical ordering to mu_n   */
int rsrec_chebyshev_moments(rsrec_ctx *ctx, const void *psi0, int lld,
                            double a, double b, void *mu_out);

/* Block Lanczos (Haydock) recursion, mirrors crecal_b/hop_b.
 *   psi0   : (nb, nb, kk) starting block state
 *   lld    : recursion steps (control%lld)
 *   a_b    : (nb, nb, lld) out, = atemp_b
 *   b2_b   : (nb, nb, lld) out, = b2temp_b  (B^2 BEFORE zsqr, like the CPU) */
int rsrec_block_lanczos(rsrec_ctx *ctx, const void *psi0, int lld,
                        void *a_b, void *b2_b);

/* Scalar Lanczos, mirrors recur()/crecal()/hop() but runs all nb orbital
 * chains of one site simultaneously as independent columns.
 *   site_j : 1-based atom index in the clust (lattice%irec(i))
 *   lld    : recursion steps
 *   a_out  : (lld, nb) out  -> this%a (ll, l, i, 1)
 *   b2_out : (lld, nb) out  -> this%b2(ll, l, i, 1)
 * NOTE: equals the nsp==1 CPU path exactly when the spin-off-diagonal
 * blocks of H are zero (the only regime recur() is used in).               */
int rsrec_scalar_lanczos(rsrec_ctx *ctx, int site_j, int lld,
                         double *a_out, double *b2_out);

/* Stochastic / conductivity Chebyshev moment matrix, mirrors the non-hoh
 * path of compute_moments_stochastic for one reference state:
 *   left  : T_{m-1}(H~) |psiref>                       m = 1..lld
 *   right : v_a T_{n-1}(H~) v_b |psiref>               n = 1..lld
 *   mu_nm(:,:,n,m) = sum_k  left_m(k)^H right_n(k)
 *   psiref: (nb, nb, kk),  mu_nm: (nb, nb, lld, lld)
 * Memory: keeps all lld left states resident (lld*kk*nb*nb*16 bytes); the
 * call fails with a clear message if this does not fit on the device.      */
int rsrec_stochastic_moments(rsrec_ctx *ctx, const void *psiref, int lld,
                             double a, double b, void *mu_nm);

#ifdef __cplusplus
}
#endif
#endif /* RSREC_H */
