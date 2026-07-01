/* ===========================================================================
 * hambuild_gpu.cu -- CUDA (with HIP hooks) device kernels for Hamiltonian
 *                    assembly.
 *
 * Phase 1: on-site blocks (build_lsham / build_obarm / build_enim) plus the two
 * reusable device primitives pack_spinor and hcpx_dev. Later phases add the
 * neighbour hot loop (Phase 2), CCOR (Phase 3) and Hubbard (Phase 4).
 *
 * Bit-exactness: the kernels mirror the CPU routines in hamiltonian.f90
 * (build_lsham / build_obarm / build_enim) and math.f90 (hcpx) operation for
 * operation. hcpx_dev performs the same two matmuls (vc*ham*v) hcpx does; the
 * v/vc/L matrices are uploaded from the host so no transform logic is
 * re-derived on the device.
 *
 * --- HIP portability hooks ---------------------------------------------------
 * Authored in CUDA. hipify-perl maps cuda* -> hip*; the shim below lets the same
 * source compile under either toolchain. Guard device API calls behind the HB_*
 * aliases rather than raw cuda or hip names.
 * =========================================================================== */

#if defined(__HIP_PLATFORM_AMD__) || defined(USE_HIP)
  #include <hip/hip_runtime.h>
  #define HB_MALLOC       hipMalloc
  #define HB_FREE         hipFree
  #define HB_MEMCPY       hipMemcpy
  #define HB_MEMSET       hipMemset
  #define HB_MEMCPY_H2D   hipMemcpyHostToDevice
  #define HB_MEMCPY_D2H   hipMemcpyDeviceToHost
  #define HB_DEVICE_SYNC  hipDeviceSynchronize
  #define HB_GET_LAST_ERR hipGetLastError
  #define HB_SUCCESS      hipSuccess
  typedef hipError_t hb_error_t;
#else
  #include <cuda_runtime.h>
  #define HB_MALLOC       cudaMalloc
  #define HB_FREE         cudaFree
  #define HB_MEMCPY       cudaMemcpy
  #define HB_MEMSET       cudaMemset
  #define HB_MEMCPY_H2D   cudaMemcpyHostToDevice
  #define HB_MEMCPY_D2H   cudaMemcpyDeviceToHost
  #define HB_DEVICE_SYNC  cudaDeviceSynchronize
  #define HB_GET_LAST_ERR cudaGetLastError
  #define HB_SUCCESS      cudaSuccess
  typedef cudaError_t hb_error_t;
#endif

#include <cuComplex.h>

/* Match the Fortran math_mod pi exactly (bit-for-bit spin-spiral phase). */
#define HB_PI 3.14159265358979323846

/* Interleaved complex(rp) with rp=real64 maps 1:1 to cuDoubleComplex. */
typedef cuDoubleComplex hbc;

__device__ __forceinline__ hbc hb_add(hbc a, hbc b) { return cuCadd(a, b); }
__device__ __forceinline__ hbc hb_sub(hbc a, hbc b) { return cuCsub(a, b); }
__device__ __forceinline__ hbc hb_mul(hbc a, hbc b) { return cuCmul(a, b); }
__device__ __forceinline__ hbc hb_scale(hbc a, double s) {
    return make_cuDoubleComplex(cuCreal(a) * s, cuCimag(a) * s);
}
/* i * z */
__device__ __forceinline__ hbc hb_imul(hbc z) {
    return make_cuDoubleComplex(-cuCimag(z), cuCreal(z));
}

/* Column-major (n x n) indexing. */
__device__ __forceinline__ int cmidx(int r, int c, int n) { return r + c * n; }

/* ---------------------------------------------------------------------------
 * hcpx_dev: cart<->sph transform of an (norb x norb) block, mirroring
 * math.f90:hcpx. cart2sph:  out = vc * blk * v.  Only cart2sph is needed for
 * the on-site blocks. v, vc are resident (norb x norb) constants.
 *
 * Runs as one block of norb*norb threads: each thread computes one output
 * element out(i,j) = sum_{a,b} vc(i,a) blk(a,b) v(b,j). norb <= 16 so the
 * O(norb^2) inner work per element is cheap and matches the CPU matmul order.
 * ------------------------------------------------------------------------- */
__device__ void hcpx_dev_cart2sph(const hbc *blk, const hbc *v, const hbc *vc,
                                  hbc *out, int norb) {
    /* tmp = blk * v ; out = vc * tmp   (matches htmp=matmul(ham,v); vc*htmp) */
    int tid = threadIdx.x;
    int nn = norb * norb;
    extern __shared__ hbc s_tmp[]; /* norb*norb scratch for (blk*v) */

    for (int idx = tid; idx < nn; idx += blockDim.x) {
        int i = idx % norb;
        int j = idx / norb;
        hbc acc = make_cuDoubleComplex(0.0, 0.0);
        for (int b = 0; b < norb; ++b)
            acc = hb_add(acc, hb_mul(blk[cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
        s_tmp[idx] = acc;
    }
    __syncthreads();
    for (int idx = tid; idx < nn; idx += blockDim.x) {
        int i = idx % norb;
        int j = idx / norb;
        hbc acc = make_cuDoubleComplex(0.0, 0.0);
        for (int a = 0; a < norb; ++a)
            acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], s_tmp[cmidx(a, j, norb)]));
        out[idx] = acc;
    }
    __syncthreads();
}

/* ===========================================================================
 * build_obarm  (mirror of hamiltonian.f90:build_obarm)
 *
 * Per type: obm0/obm1 diagonal from obx0/obx1; pack into 4 quadrants using
 * mom(1:3); then hcpx cart2sph on each quadrant. One thread block per type.
 * ------------------------------------------------------------------------- */
__global__ void k_build_obarm(const hbc *obx0, const hbc *obx1,
                              const double *mom, const hbc *v, const hbc *vc,
                              hbc *obarm, int norb, int nb, int ntype) {
    int t = blockIdx.x;
    if (t >= ntype) return;
    int tid = threadIdx.x;
    int nn = norb * norb;

    const hbc *o0 = obx0 + (size_t)t * norb;
    const hbc *o1 = obx1 + (size_t)t * norb;
    double m1 = mom[3 * t + 0], m2 = mom[3 * t + 1], m3 = mom[3 * t + 2];
    hbc *out = obarm + (size_t)t * nb * nb;

    /* Build the four cartesian norb x norb quadrants into shared scratch. */
    extern __shared__ hbc sh[]; /* layout: q11,q22,q12,q21 (nn each) then hcpx tmp (nn) */
    hbc *q11 = sh, *q22 = sh + nn, *q12 = sh + 2 * nn, *q21 = sh + 3 * nn;

    for (int idx = tid; idx < nn; idx += blockDim.x) {
        int m = idx % norb; /* row */
        int l = idx / norb; /* col */
        /* obm0/obm1 are diagonal: obm(m,l) nonzero only for m==l. */
        hbc o0ml = (m == l) ? o0[m] : make_cuDoubleComplex(0.0, 0.0);
        hbc o1ml = (m == l) ? o1[m] : make_cuDoubleComplex(0.0, 0.0);
        /* obarm(m,l)         = obm0 + obm1*mom(3)           (H11) */
        q11[cmidx(m, l, norb)] = hb_add(o0ml, hb_scale(o1ml, m3));
        /* obarm(m+off,l+off) = obm0 - obm1*mom(3)           (H22) */
        q22[cmidx(m, l, norb)] = hb_sub(o0ml, hb_scale(o1ml, m3));
        /* CPU writes obarm(l,m+off) = obm1(m,l)*(mom1 - i*mom2), i.e. the (l,m)
         * element of the upper-right quadrant. Store into q12 at (l,m). */
        hbc off = hb_sub(hb_scale(o1ml, m1), hb_imul(hb_scale(o1ml, m2)));
        q12[cmidx(l, m, norb)] = off;
        /* obarm(l+off,m) = obm1(m,l)*(mom1 + i*mom2) -> q21 at (l,m). */
        hbc ofc = hb_add(hb_scale(o1ml, m1), hb_imul(hb_scale(o1ml, m2)));
        q21[cmidx(l, m, norb)] = ofc;
    }
    __syncthreads();

    /* hcpx cart2sph on each quadrant, then scatter into nb x nb output. */
    hbc *tmp = sh + 4 * nn; /* hcpx scratch (nn) */
    hbc *res = tmp + nn;    /* result of one quadrant (nn) */
    hbc *quads[4] = {q11, q22, q12, q21};
    for (int qq = 0; qq < 4; ++qq) {
        /* hcpx_dev_cart2sph uses its own extern shared; call an inlined version
         * over tmp to avoid clashing. Reuse the same two-matmul logic here. */
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int b = 0; b < norb; ++b)
                acc = hb_add(acc, hb_mul(quads[qq][cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
            tmp[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int a = 0; a < norb; ++a)
                acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], tmp[cmidx(a, j, norb)]));
            res[idx] = acc;
        }
        __syncthreads();
        /* scatter res into the proper nb x nb quadrant */
        int roff = (qq == 1 || qq == 3) ? norb : 0; /* q22,q21 shift rows */
        int coff = (qq == 1 || qq == 2) ? norb : 0; /* q22,q12 shift cols */
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            out[cmidx(i + roff, j + coff, nb)] = res[idx];
        }
        __syncthreads();
    }
}

/* ===========================================================================
 * build_enim  (mirror of hamiltonian.f90:build_enim)
 *
 * Per type: eu=cx(m,1)-cex(m,1), ed=cx(m,2)-cex(m,2); ex0=(eu+ed)/2,
 * ex1=(eu-ed)/2; diagonal seeds em0/em1; pack 4 quadrants with mom; hcpx.
 * ------------------------------------------------------------------------- */
__global__ void k_build_enim(const hbc *cx, const hbc *cex, const double *mom,
                             const hbc *v, const hbc *vc, hbc *enim, int norb,
                             int nb, int ntype) {
    int t = blockIdx.x;
    if (t >= ntype) return;
    int tid = threadIdx.x;
    int nn = norb * norb;

    const hbc *cxt = cx + (size_t)t * norb * 2;   /* (norb,2) col-major */
    const hbc *cext = cex + (size_t)t * norb * 2;
    double m1 = mom[3 * t + 0], m2 = mom[3 * t + 1], m3 = mom[3 * t + 2];
    hbc *out = enim + (size_t)t * nb * nb;

    extern __shared__ hbc sh[];
    hbc *q11 = sh, *q22 = sh + nn, *q12 = sh + 2 * nn, *q21 = sh + 3 * nn;

    for (int idx = tid; idx < nn; idx += blockDim.x) {
        int m = idx % norb, l = idx / norb;
        hbc e0ml = make_cuDoubleComplex(0.0, 0.0);
        hbc e1ml = make_cuDoubleComplex(0.0, 0.0);
        if (m == l) {
            hbc eu = hb_sub(cxt[cmidx(m, 0, norb)], cext[cmidx(m, 0, norb)]);
            hbc ed = hb_sub(cxt[cmidx(m, 1, norb)], cext[cmidx(m, 1, norb)]);
            e0ml = hb_scale(hb_add(eu, ed), 0.5);
            e1ml = hb_scale(hb_sub(eu, ed), 0.5);
        }
        q11[cmidx(m, l, norb)] = hb_add(e0ml, hb_scale(e1ml, m3));
        q22[cmidx(m, l, norb)] = hb_sub(e0ml, hb_scale(e1ml, m3));
        q12[cmidx(l, m, norb)] = hb_sub(hb_scale(e1ml, m1), hb_imul(hb_scale(e1ml, m2)));
        q21[cmidx(l, m, norb)] = hb_add(hb_scale(e1ml, m1), hb_imul(hb_scale(e1ml, m2)));
    }
    __syncthreads();

    hbc *tmp = sh + 4 * nn;
    hbc *res = tmp + nn;
    hbc *quads[4] = {q11, q22, q12, q21};
    for (int qq = 0; qq < 4; ++qq) {
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int b = 0; b < norb; ++b)
                acc = hb_add(acc, hb_mul(quads[qq][cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
            tmp[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int a = 0; a < norb; ++a)
                acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], tmp[cmidx(a, j, norb)]));
            res[idx] = acc;
        }
        __syncthreads();
        int roff = (qq == 1 || qq == 3) ? norb : 0;
        int coff = (qq == 1 || qq == 2) ? norb : 0;
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            out[cmidx(i + roff, j + coff, nb)] = res[idx];
        }
        __syncthreads();
    }
}

/* ===========================================================================
 * build_lsham  (mirror of hamiltonian.f90:build_lsham)
 *
 * Lx/Ly/Lz are hcpx-transformed ONCE (cart2sph) on the host-uploaded cartesian
 * operators -- but build_lsham transforms them before the type loop, so the
 * transform is type-independent. We transform them once here into shared, then
 * every type packs the L.S block. prefac depends on orbital l-block via soc_p /
 * soc_d (p: idx 2..4, d: idx 5..9, f: 10..16 -> 0). orb_pol adds a rac/lmom term
 * to the diagonal (H11/H22) via Lz only.
 *
 * Because the transformed L operators are shared across types, we transform them
 * in block 0-style: every block re-transforms them into its own shared memory
 * (cheap, norb<=16) to avoid cross-block globals. One block per type.
 * ------------------------------------------------------------------------- */
__global__ void k_build_lsham(const hbc *lx, const hbc *ly, const hbc *lz,
                              const double *xi_p, const double *xi_d,
                              const double *rac, const double *lmom,
                              const hbc *v, const hbc *vc, int orb_pol,
                              hbc *lsham, int norb, int nb, int ntype) {
    int t = blockIdx.x;
    if (t >= ntype) return;
    int tid = threadIdx.x;
    int nn = norb * norb;

    extern __shared__ hbc sh[];
    hbc *sLx = sh, *sLy = sh + nn, *sLz = sh + 2 * nn;
    hbc *tmp = sh + 3 * nn; /* hcpx scratch */

    /* transform Lx,Ly,Lz cart2sph into sLx,sLy,sLz */
    const hbc *Lsrc[3] = {lx, ly, lz};
    hbc *Ldst[3] = {sLx, sLy, sLz};
    for (int c = 0; c < 3; ++c) {
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int b = 0; b < norb; ++b)
                acc = hb_add(acc, hb_mul(Lsrc[c][cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
            tmp[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int a = 0; a < norb; ++a)
                acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], tmp[cmidx(a, j, norb)]));
            Ldst[c][idx] = acc;
        }
        __syncthreads();
    }

    double soc_p = sqrt(xi_p[2 * t + 0] * xi_p[2 * t + 1]);
    double soc_d = sqrt(xi_d[2 * t + 0] * xi_d[2 * t + 1]);
    /* orb_pol terms: rac(1:2), lz_loc = sqrt(xi_d(1)*lmom(3)); rac already
     * sqrt(xi_d(1)*rac) on host. See note below -- host precomputes rac_eff and
     * lz_loc to keep bit-exactness with the CPU expression order. */
    double rac1 = rac[2 * t + 0];
    double rac2 = rac[2 * t + 1];
    double lz_loc = lmom[3 * t + 2]; /* host passes sqrt(xi_d(1)*lmom(3)) here */

    hbc *out = lsham + (size_t)t * nb * nb;
    for (int idx = tid; idx < nb * nb; idx += blockDim.x)
        out[idx] = make_cuDoubleComplex(0.0, 0.0);
    __syncthreads();

    /* pack: loop over (i,j) in norb x norb (CPU loops i=1..norb, j=1..norb and
     * writes lsham(j,i,...)). prefac is per (i,j) block. */
    for (int idx = tid; idx < nn; idx += blockDim.x) {
        int i = idx / norb; /* CPU outer i */
        int j = idx % norb; /* CPU inner j */
        /* CPU keeps prefac "sticky" across (i,j) (set only inside the l-block
         * conditionals, never reset per element). For cross-l elements no
         * condition fires, but there Lz/Lx/Ly are exactly zero (L is block-
         * diagonal in l), so pf*0 == stale_pf*0. Resetting pf=0 per element is
         * therefore bit-identical to the CPU. */
        double pf = 0.0;
        if (i >= 1 && i <= 3 && j >= 1 && j <= 3) pf = 0.5 * soc_p; /* p: 2..4 (0-based 1..3) */
        if (i >= 4 && i <= 8 && j >= 4 && j <= 8) pf = 0.5 * soc_d; /* d: 5..9 (0-based 4..8) */
        /* f block (0-based 9..15) -> 0 */

        hbc lzji = sLz[cmidx(j, i, norb)];
        hbc lxji = sLx[cmidx(j, i, norb)];
        hbc lyji = sLy[cmidx(j, i, norb)];

        /* H11: prefac*Lz + Lz*rac(1)*lz_loc */
        hbc h11 = hb_scale(lzji, pf);
        if (orb_pol) h11 = hb_add(h11, hb_scale(lzji, rac1 * lz_loc));
        out[cmidx(j, i, nb)] = hb_add(out[cmidx(j, i, nb)], h11);

        /* H12: prefac*(Lx - i*Ly) at (j, i+off) */
        hbc h12 = hb_scale(hb_sub(lxji, hb_imul(lyji)), pf);
        out[cmidx(j, i + norb, nb)] = hb_add(out[cmidx(j, i + norb, nb)], h12);

        /* H21: prefac*(Lx + i*Ly) at (j+off, i) */
        hbc h21 = hb_scale(hb_add(lxji, hb_imul(lyji)), pf);
        out[cmidx(j + norb, i, nb)] = hb_add(out[cmidx(j + norb, i, nb)], h21);

        /* H22: -prefac*Lz - Lz*rac(2)*lz_loc at (j+off, i+off) */
        hbc h22 = hb_scale(lzji, -pf);
        if (orb_pol) h22 = hb_sub(h22, hb_scale(lzji, rac2 * lz_loc));
        out[cmidx(j + norb, i + norb, nb)] = hb_add(out[cmidx(j + norb, i + norb, nb)], h22);
    }
}

/* ===========================================================================
 * Host-callable entry points (called by hambuild_cuda.cpp).
 * Return 0 on success, nonzero on error (message via cudaGetLastError path).
 * These allocate/keep device buffers on the ctx side; here they take raw device
 * pointers already staged by the dispatcher.
 * ========================================================================= */
extern "C" int hambuild_gpu_launch_onsite(
    /* constants (device) */ const hbc *d_v, const hbc *d_vc, const hbc *d_lx,
    const hbc *d_ly, const hbc *d_lz,
    /* potential (device) */ const double *d_xi_p, const double *d_xi_d,
    const double *d_rac, const hbc *d_obx0, const hbc *d_obx1, const hbc *d_cx,
    const hbc *d_cex, const double *d_mom, const double *d_lmom, int orb_pol,
    /* outputs (device) */ hbc *d_lsham, hbc *d_obarm, hbc *d_enim,
    int norb, int nb, int ntype, int want_lsham) {
    int nn = norb * norb;
    int threads = 256;
    dim3 grid(ntype);

    if (d_obarm) {
        size_t shmem = (size_t)(4 * nn + 2 * nn) * sizeof(hbc);
        k_build_obarm<<<grid, threads, shmem>>>(d_obx0, d_obx1, d_mom, d_v, d_vc,
                                                d_obarm, norb, nb, ntype);
    }
    if (d_enim) {
        size_t shmem = (size_t)(4 * nn + 2 * nn) * sizeof(hbc);
        k_build_enim<<<grid, threads, shmem>>>(d_cx, d_cex, d_mom, d_v, d_vc,
                                               d_enim, norb, nb, ntype);
    }
    if (want_lsham && d_lsham) {
        size_t shmem = (size_t)(4 * nn) * sizeof(hbc); /* sLx,sLy,sLz + tmp */
        k_build_lsham<<<grid, threads, shmem>>>(d_lx, d_ly, d_lz, d_xi_p, d_xi_d,
                                                d_rac, d_lmom, d_v, d_vc, orb_pol,
                                                d_lsham, norb, nb, ntype);
    }
    if (HB_DEVICE_SYNC() != HB_SUCCESS) return 1;
    if (HB_GET_LAST_ERR() != HB_SUCCESS) return 2;
    return 0;
}

/* ===========================================================================
 * Phase 1.5: geometry maps (clusba filter + hmfind match), once per run.
 *
 * Mirrors hamiltonian.f90:chbar_nc's per-(type,neighbour) preamble and
 * lattice.f90:clusba / hamiltonian.f90:hmfind. One thread block per atom type.
 *
 * All integer geometry arrays are the raw 1-based Fortran indices; we subtract 1
 * where used as a 0-based offset. Column-major throughout:
 *   cr(3,kk):   cr[i + 3*a]            (a = 0-based cluster atom)
 *   nn(ndi,nn_max): nn[ia + ndi*(m)]   (ia,m 0-based; nn(ia,1)->m=0 is the count)
 *
 * Outputs per (m, type) with column-major (nn_max, ntype):
 *   valid[m + nn_max*t]  = ni (1 found / 0 not found), 0 when jj==0
 *   shell[m + nn_max*t]  = m+1 (1-based sbar shell index; matches CPU)
 *   ino[m + nn_max*t]    = num(ia) (1-based sbar type index)
 *   vet[3*(m + nn_max*t) + i] = displacement vector components (already *alat)
 * ------------------------------------------------------------------------- */
__global__ void k_build_geometry_maps(
    const double *cr, const int *num, const int *iz, const int *nn,
    const int *atlist, const double *vet_in, int use_vet_in, int kk, int nn_max,
    int ndi, double alat, double r2, double *hv_scratch, int *valid, int *shell,
    int *ino, double *vet) {
    int t = blockIdx.x;        /* atom type (0-based) */
    int tid = threadIdx.x;

    int ia = atlist[t] - 1;    /* cluster atom index for this type (0-based) */
    int numia = num[ia];       /* 1-based sbar type index */
    int nr = nn[ia + ndi * 0]; /* nn(ia,1): neighbour count */

    /* clusba: build ham_vec for this ia -- the compacted list of displacement
     * vectors of cluster atoms within r2, in sequential order (index 1 is the
     * zero self-vector). Stored in a per-block slice of global scratch (3 x kk)
     * so it scales past the shared-memory limit; this kernel runs once per run.
     * The sequential compaction (ii=ii+1) is done by a single thread to preserve
     * the exact ordering hmfind searches. */
    double *hv = hv_scratch + (size_t)3 * kk * t; /* this block's ham_vec */
    __shared__ int s_count;

    /* CPU passes cralat = cr*alat to clusba, so ham_vec (and the r2 threshold)
     * are in alat units. We mirror that: hv holds (cr diff)*alat and r2 is the
     * same alat^2 threshold the CPU uses. */
    if (tid == 0) {
        int ii = 0;            /* 0-based; CPU ii starts at 1 with entry 1 = 0 */
        hv[3 * ii + 0] = 0.0;
        hv[3 * ii + 1] = 0.0;
        hv[3 * ii + 2] = 0.0;
        for (int a = 0; a < kk; ++a) {
            double d0 = (cr[3 * a + 0] - cr[3 * ia + 0]) * alat;
            double d1 = (cr[3 * a + 1] - cr[3 * ia + 1]) * alat;
            double d2 = (cr[3 * a + 2] - cr[3 * ia + 2]) * alat;
            double s1 = d0 * d0 + d1 * d1 + d2 * d2;
            if (s1 < r2 && s1 > 0.0001) {
                ii += 1;
                hv[3 * ii + 0] = d0;
                hv[3 * ii + 1] = d1;
                hv[3 * ii + 2] = d2;
            }
        }
        s_count = ii + 1; /* number of valid entries (CPU n = ii) */
    }
    __syncthreads();

    const double eps = 0.0001;
    /* Each thread handles a set of neighbours m. */
    for (int m = tid; m < nr; m += blockDim.x) {
        int outidx = m + nn_max * t;
        /* jj = nn(ia, m+1); for m==0 (CPU m==1) jj = ia */
        int jj;
        if (m == 0)
            jj = ia + 1; /* keep 1-based to match nn storage; converted below */
        else
            jj = nn[ia + ndi * m];

        shell[outidx] = m + 1;
        ino[outidx] = numia;

        if (jj == 0) {
            valid[outidx] = 0;
            vet[3 * outidx + 0] = 0.0;
            vet[3 * outidx + 1] = 0.0;
            vet[3 * outidx + 2] = 0.0;
            continue;
        }
        int jj0 = jj - 1; /* 0-based */

        double v0, v1, v2;
        if (use_vet_in) {
            v0 = vet_in[3 * outidx + 0];
            v1 = vet_in[3 * outidx + 1];
            v2 = vet_in[3 * outidx + 2];
        } else {
            v0 = (cr[3 * jj0 + 0] - cr[3 * ia + 0]) * alat;
            v1 = (cr[3 * jj0 + 1] - cr[3 * ia + 1]) * alat;
            v2 = (cr[3 * jj0 + 2] - cr[3 * ia + 2]) * alat;
        }
        vet[3 * outidx + 0] = v0;
        vet[3 * outidx + 1] = v1;
        vet[3 * outidx + 2] = v2;

        /* hmfind: search ham_vec (s_hv, alat units) for a match to vet (also
         * alat units); ni=1 if found else 0. */
        int ni = 0;
        for (int i = 0; i < s_count; ++i) {
            double a0 = v0 - hv[3 * i + 0];
            double a1 = v1 - hv[3 * i + 1];
            double a2 = v2 - hv[3 * i + 2];
            if (a0 * a0 + a1 * a1 + a2 * a2 < eps) { ni = 1; break; }
        }
        valid[outidx] = ni;
    }
}

extern "C" int hambuild_gpu_launch_geometry_maps(
    const double *d_cr, const int *d_num, const int *d_iz, const int *d_nn,
    const int *d_atlist, const double *d_vet_in, int use_vet_in, int kk,
    int nn_max, int ndi, double alat, double r2, int ntype, double *d_hv_scratch,
    int *d_valid, int *d_shell, int *d_ino, double *d_vet) {
    int threads = 128;
    dim3 grid(ntype);
    k_build_geometry_maps<<<grid, threads>>>(
        d_cr, d_num, d_iz, d_nn, d_atlist, d_vet_in, use_vet_in, kk, nn_max, ndi,
        alat, r2, d_hv_scratch, d_valid, d_shell, d_ino, d_vet);
    if (HB_DEVICE_SYNC() != HB_SUCCESS) return 1;
    if (HB_GET_LAST_ERR() != HB_SUCCESS) return 2;
    return 0;
}

/* ===========================================================================
 * Phase 2: bulk/local neighbour Hamiltonian (chbar_nc body).
 *
 * Mirrors hamiltonian.f90:ham0m_nc + the 4x hcpx in chbar_nc + the spinor pack
 * in build_bulkham. One thread block per (type t, neighbour m). Threads cover
 * norb*norb elements.
 *
 * hhh(ilm,jlm) = real(sbar(jlm, ilm, shell, ino))  (note transpose + real, as
 * the CPU hmfind does). shell/ino are 1-based; convert to 0-based here.
 * ------------------------------------------------------------------------- */
__global__ void k_build_bulk_ham(
    const double *sbar_re, int norb_s, int nm_store, int ntot,
    const hbc *wx0, const hbc *wx1, const hbc *cx0, const hbc *cx1,
    const hbc *cex0, const hbc *cex1, const double *mom, const double *q_ss,
    double theta_ss, const hbc *v, const hbc *vc,
    const int *iz, const int *nn, const int *atlist, const double *cr,
    const int *valid, const int *shell, const int *ino, const double *vet,
    int ndi, int kk, int hoh, hbc *ee, hbc *hxc, int norb, int nb, int nn_max,
    int ntype) {
    int t = blockIdx.y;         /* atom type (0-based) */
    int m = blockIdx.x;         /* neighbour (0-based) */
    if (t >= ntype || m >= nn_max) return;
    int tid = threadIdx.x;
    int nn2 = norb * norb;

    int ia = atlist[t] - 1;     /* 0-based cluster atom */
    int nr = nn[ia + ndi * 0];  /* neighbour count */
    int mapidx = m + nn_max * t;

    hbc *out_ee = ee + (size_t)nb * nb * (m + nn_max * t);
    hbc *out_hxc = hxc + (size_t)nb * nb * (m + nn_max * t);

    /* Zero the whole nb x nb output block first (covers absent neighbours). */
    for (int idx = tid; idx < nb * nb; idx += blockDim.x) {
        out_ee[idx] = make_cuDoubleComplex(0.0, 0.0);
        out_hxc[idx] = make_cuDoubleComplex(0.0, 0.0);
    }
    __syncthreads();

    if (m >= nr) return;        /* only nr neighbours are meaningful */

    /* jj = nn(ia,m+1); m==0 -> jj=ia (1-based bookkeeping like CPU). */
    int jj1 = (m == 0) ? (ia + 1) : nn[ia + ndi * m];
    if (jj1 == 0) return;       /* absent neighbour: leave zeros */
    int jj = jj1 - 1;
    int it = iz[ia];            /* 1-based type indices for potential lookup */
    int jt = iz[jj];

    int sh = shell[mapidx] - 1; /* 0-based sbar shell */
    int io = ino[mapidx] - 1;   /* 0-based sbar type */

    /* moments (with optional spin-spiral override) */
    double mia[3], mja[3];
    const double *momit = mom + 3 * (it - 1);
    const double *momjt = mom + 3 * (jt - 1);
    mia[0] = momit[0]; mia[1] = momit[1]; mia[2] = momit[2];
    mja[0] = momjt[0]; mja[1] = momjt[1]; mja[2] = momjt[2];
    double qn = q_ss[0]*q_ss[0] + q_ss[1]*q_ss[1] + q_ss[2]*q_ss[2];
    if (sqrt(qn) > 1.0e-5 || fabs(sin(theta_ss)) > 1.0e-8) {
        const double *ria = cr + 3 * ia;
        const double *rja = cr + 3 * jj;
        double pia = 2.0 * HB_PI * (ria[0]*q_ss[0] + ria[1]*q_ss[1] + ria[2]*q_ss[2]);
        double pja = 2.0 * HB_PI * (rja[0]*q_ss[0] + rja[1]*q_ss[1] + rja[2]*q_ss[2]);
        double st = sin(theta_ss), ct = cos(theta_ss);
        mia[0] = cos(pia)*st; mia[1] = sin(pia)*st; mia[2] = ct;
        mja[0] = cos(pja)*st; mja[1] = sin(pja)*st; mja[2] = ct;
    }
    double dotp = mia[0]*mja[0] + mia[1]*mja[1] + mia[2]*mja[2];
    double crx = mia[1]*mja[2] - mia[2]*mja[1];
    double cry = mia[2]*mja[0] - mia[0]*mja[2];
    double crz = mia[0]*mja[1] - mia[1]*mja[0];
    double crs[3] = {crx, cry, crz};

    double vx = vet[3*mapidx+0], vy = vet[3*mapidx+1], vz = vet[3*mapidx+2];
    double vv = sqrt(vx*vx + vy*vy + vz*vz);

    /* per-type potential vectors */
    const hbc *wx0i = wx0 + (size_t)norb*(it-1), *wx1i = wx1 + (size_t)norb*(it-1);
    const hbc *wx0j = wx0 + (size_t)norb*(jt-1), *wx1j = wx1 + (size_t)norb*(jt-1);
    const hbc *cx0i = cx0 + (size_t)norb*(it-1), *cx1i = cx1 + (size_t)norb*(it-1);
    const hbc *cex0i = cex0 + (size_t)norb*(it-1), *cex1i = cex1 + (size_t)norb*(it-1);

    /* shared: 4 cart matrices (nn2 each) + hcpx tmp (nn2) + hcpx res (nn2) */
    extern __shared__ hbc shp[];
    hbc *H1 = shp, *H2 = shp + nn2, *H3 = shp + 2*nn2, *H4 = shp + 3*nn2;
    hbc *tmp = shp + 4*nn2, *res = shp + 5*nn2;

    /* Build the four cart matrices H(:,:,mdir) exactly as ham0m_nc. Indices:
     * ilm = row, jlm = col (Fortran hhmag(ilm,jlm,mdir)); we store col-major. */
    for (int idx = tid; idx < nn2; idx += blockDim.x) {
        int ilm = idx % norb, jlm = idx / norb;
        double hh = sbar_re[(size_t)jlm + norb_s*((size_t)ilm + norb_s*((size_t)sh + nm_store*(size_t)io))];
        hbc hc = make_cuDoubleComplex(hh, 0.0);
        hbc a = hb_mul(hb_mul(wx0i[ilm], hc), wx0j[jlm]);
        hbc b = hb_scale(hb_mul(hb_mul(wx1i[ilm], hc), wx1j[jlm]), dotp);
        hbc H4v = hb_add(a, b);
        if (vv <= 0.01 && ilm == jlm) {
            H4v = hb_add(H4v, hoh ? cex0i[ilm] : cx0i[ilm]);
        }
        H4[idx] = H4v;
        hbc w1i_hh_w0j = hb_mul(hb_mul(wx1i[ilm], hc), wx0j[jlm]);
        hbc w0i_hh_w1j = hb_mul(hb_mul(wx0i[ilm], hc), wx1j[jlm]);
        hbc w1i_hh_w1j = hb_mul(hb_mul(wx1i[ilm], hc), wx1j[jlm]);
        hbc *Hm[3] = {H1, H2, H3};
        for (int mm = 0; mm < 3; ++mm) {
            hbc term = hb_add(hb_scale(w1i_hh_w0j, mia[mm]),
                              hb_scale(w0i_hh_w1j, mja[mm]));
            term = hb_add(term, hb_imul(hb_scale(w1i_hh_w1j, crs[mm])));
            if (vv <= 0.01 && ilm == jlm) {
                hbc add = hb_scale(hoh ? cex1i[ilm] : cx1i[ilm], mia[mm]);
                term = hb_add(term, add);
            }
            Hm[mm][idx] = term;
        }
    }
    __syncthreads();

    /* hcpx cart2sph each of the 4 in place (via tmp/res). */
    hbc *slabs[4] = {H1, H2, H3, H4};
    for (int s = 0; s < 4; ++s) {
        for (int idx = tid; idx < nn2; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int b = 0; b < norb; ++b)
                acc = hb_add(acc, hb_mul(slabs[s][cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
            tmp[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn2; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int a = 0; a < norb; ++a)
                acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], tmp[cmidx(a, j, norb)]));
            res[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn2; idx += blockDim.x) slabs[s][idx] = res[idx];
        __syncthreads();
    }

    /* pack into ee/hxc: CPU indexes hmag(j,i,m,mdir), j=row, i=col. */
    for (int idx = tid; idx < nn2; idx += blockDim.x) {
        int j = idx % norb, i = idx / norb;
        hbc h1 = H1[cmidx(j, i, norb)];
        hbc h2 = H2[cmidx(j, i, norb)];
        hbc h3 = H3[cmidx(j, i, norb)];
        hbc h4 = H4[cmidx(j, i, norb)];
        hbc h1mih2 = hb_sub(h1, hb_imul(h2));
        hbc h1pih2 = hb_add(h1, hb_imul(h2));
        out_ee[cmidx(j, i, nb)]            = hb_add(h4, h3);
        out_ee[cmidx(j+norb, i+norb, nb)]  = hb_sub(h4, h3);
        out_ee[cmidx(j, i+norb, nb)]       = h1mih2;
        out_ee[cmidx(j+norb, i, nb)]       = h1pih2;
        out_hxc[cmidx(j, i, nb)]           = h3;
        out_hxc[cmidx(j+norb, i+norb, nb)] = hb_sub(make_cuDoubleComplex(0.0,0.0), h3);
        out_hxc[cmidx(j, i+norb, nb)]      = h1mih2;
        out_hxc[cmidx(j+norb, i, nb)]      = h1pih2;
    }
}

extern "C" int hambuild_gpu_launch_bulk_ham(
    const double *d_sbar_re, int norb_s, int nm_store, int ntot,
    const hbc *d_wx0, const hbc *d_wx1, const hbc *d_cx0, const hbc *d_cx1,
    const hbc *d_cex0, const hbc *d_cex1, const double *d_mom,
    const double *d_q_ss, double theta_ss, const hbc *d_v, const hbc *d_vc,
    const int *d_iz, const int *d_nn, const int *d_atlist, const double *d_cr,
    const int *d_valid, const int *d_shell, const int *d_ino, const double *d_vet,
    int ndi, int kk, int hoh, hbc *d_ee, hbc *d_hxc, int norb, int nb,
    int nn_max, int ntype) {
    int threads = 128;
    dim3 grid(nn_max, ntype);
    size_t shmem = (size_t)6 * norb * norb * sizeof(hbc);
    k_build_bulk_ham<<<grid, threads, shmem>>>(
        d_sbar_re, norb_s, nm_store, ntot, d_wx0, d_wx1, d_cx0, d_cx1, d_cex0,
        d_cex1, d_mom, d_q_ss, theta_ss, d_v, d_vc, d_iz, d_nn, d_atlist, d_cr,
        d_valid, d_shell, d_ino, d_vet, ndi, kk, hoh, d_ee, d_hxc, norb, nb,
        nn_max, ntype);
    if (HB_DEVICE_SYNC() != HB_SUCCESS) return 1;
    if (HB_GET_LAST_ERR() != HB_SUCCESS) return 2;
    return 0;
}

/* ===========================================================================
 * Phase 3: combined correction H_cc (mirror build_ccor_pair_block_noncollinear).
 *
 * Per (site s, neighbour m): compute the 4-component D (from sbar) and Ddot
 * (from -avw^2*sdot) via the same ham0m arithmetic as Phase 2 but with the
 * on-site band-centre add disabled (onsite='none'); combine with the host-
 * computed coefficients ccd; add the on-site D^(0) diagonal for m==1; hcpx each
 * of the 4; then lambda-pack into eecc/hallcc. ccd is (norb, 0:2, ntype) col-
 * major; lambda is (ntype, ntype) col-major (it,jt). One block per (site,m).
 * ------------------------------------------------------------------------- */
__global__ void k_build_ccor(
    const double *sbar_re, const double *sdot_re, int norb_s, int nm_store,
    int ntot, const hbc *wx0, const hbc *wx1, const double *mom,
    const double *q_ss, double theta_ss, const double *ccd, const double *lambda,
    double avw, const hbc *v, const hbc *vc, const int *iz, const int *nn,
    const int *atlist, const double *cr, const int *shell, const int *ino,
    const double *vet, int ndi, int kk, hbc *out, int norb, int nb, int nn_max,
    int ntype, int lambda_stride) {
    int t = blockIdx.y;
    int m = blockIdx.x;
    if (t >= ntype || m >= nn_max) return;
    int tid = threadIdx.x;
    int nn2 = norb * norb;

    int ia = atlist[t] - 1;
    int nr = nn[ia + ndi * 0];
    int mapidx = m + nn_max * t;
    hbc *o = out + (size_t)nb * nb * (m + nn_max * t);

    for (int idx = tid; idx < nb * nb; idx += blockDim.x)
        o[idx] = make_cuDoubleComplex(0.0, 0.0);
    __syncthreads();
    if (m >= nr) return;

    int jj1 = (m == 0) ? (ia + 1) : nn[ia + ndi * m];
    if (jj1 == 0) return;
    int jj = jj1 - 1;
    int it = iz[ia];   /* 1-based */
    int jt = iz[jj];
    int sh = shell[mapidx] - 1;
    int io = ino[mapidx] - 1;

    /* moments (spin-spiral override on ia, matching build_ccor_d_components +
     * ccor_apply_spin_spiral; ja uses jt's mom unless spiral). */
    double mia[3], mja[3];
    const double *momit = mom + 3 * (it - 1);
    const double *momjt = mom + 3 * (jt - 1);
    mia[0]=momit[0]; mia[1]=momit[1]; mia[2]=momit[2];
    mja[0]=momjt[0]; mja[1]=momjt[1]; mja[2]=momjt[2];
    double qn = q_ss[0]*q_ss[0]+q_ss[1]*q_ss[1]+q_ss[2]*q_ss[2];
    bool spiral = (sqrt(qn) > 1.0e-5 || fabs(sin(theta_ss)) > 1.0e-8);
    if (spiral) {
        const double *ria = cr + 3*ia, *rja = cr + 3*jj;
        double pia = 2.0*HB_PI*(ria[0]*q_ss[0]+ria[1]*q_ss[1]+ria[2]*q_ss[2]);
        double pja = 2.0*HB_PI*(rja[0]*q_ss[0]+rja[1]*q_ss[1]+rja[2]*q_ss[2]);
        double st=sin(theta_ss), ct=cos(theta_ss);
        mia[0]=cos(pia)*st; mia[1]=sin(pia)*st; mia[2]=ct;
        mja[0]=cos(pja)*st; mja[1]=sin(pja)*st; mja[2]=ct;
    }
    /* mom_i for the on-site D^(0) term (build_ccor_d_components returns mom of it,
     * then applies spin spiral on ia). */
    double momi[3] = {mia[0], mia[1], mia[2]};

    double dotp = mia[0]*mja[0]+mia[1]*mja[1]+mia[2]*mja[2];
    double crs[3] = { mia[1]*mja[2]-mia[2]*mja[1],
                      mia[2]*mja[0]-mia[0]*mja[2],
                      mia[0]*mja[1]-mia[1]*mja[0] };

    const hbc *wx0i = wx0 + (size_t)norb*(it-1), *wx1i = wx1 + (size_t)norb*(it-1);
    const hbc *wx0j = wx0 + (size_t)norb*(jt-1), *wx1j = wx1 + (size_t)norb*(jt-1);
    /* ccd is per cluster SITE (not per type): ccd[ilm + norb*(c + 3*site)],
     * c in {0,1,2}. build_ccor_coefficients depends on num(ia)/alpha(ia), so it
     * is site-dependent; index by ia and jj. */
    const double *ccdi = ccd + (size_t)norb*3*ia;
    const double *ccdj = ccd + (size_t)norb*3*jj;
    double lam = lambda[(it-1) + lambda_stride*(jt-1)];

    /* shared: K1..K4 (nn2 each) + tmp + res */
    extern __shared__ hbc shc[];
    hbc *K1=shc, *K2=shc+nn2, *K3=shc+2*nn2, *K4=shc+3*nn2;
    hbc *tmp=shc+4*nn2, *res=shc+5*nn2;

    for (int idx = tid; idx < nn2; idx += blockDim.x) {
        int ilm = idx % norb, jlm = idx / norb;
        /* D from sbar; Ddot from -avw^2 * sdot. Both use the ham0m arithmetic
         * (onsite='none' -> no CX/CEX add). */
        double sb = sbar_re[(size_t)jlm + norb_s*((size_t)ilm + norb_s*((size_t)sh + nm_store*(size_t)io))];
        double sd = -avw*avw*sdot_re[(size_t)jlm + norb_s*((size_t)ilm + norb_s*((size_t)sh + nm_store*(size_t)io))];
        hbc dcb = make_cuDoubleComplex(sb, 0.0);
        hbc ddb = make_cuDoubleComplex(sd, 0.0);

        /* coefficient factor (ccd_i(ilm,1)+ccd_j(jlm,1)) */
        double cf = ccdi[ilm + norb*1] + ccdj[jlm + norb*1];

        /* mdir=4 */
        hbc D4 = hb_add(hb_mul(hb_mul(wx0i[ilm], dcb), wx0j[jlm]),
                        hb_scale(hb_mul(hb_mul(wx1i[ilm], dcb), wx1j[jlm]), dotp));
        hbc Dd4 = hb_add(hb_mul(hb_mul(wx0i[ilm], ddb), wx0j[jlm]),
                         hb_scale(hb_mul(hb_mul(wx1i[ilm], ddb), wx1j[jlm]), dotp));
        hbc K4v = hb_add(Dd4, hb_scale(D4, cf));
        /* mdir=1..3 */
        hbc *Kp[3] = {K1, K2, K3};
        for (int mm = 0; mm < 3; ++mm) {
            hbc Dm = hb_add(hb_scale(hb_mul(hb_mul(wx1i[ilm], dcb), wx0j[jlm]), mia[mm]),
                            hb_scale(hb_mul(hb_mul(wx0i[ilm], dcb), wx1j[jlm]), mja[mm]));
            Dm = hb_add(Dm, hb_imul(hb_scale(hb_mul(hb_mul(wx1i[ilm], dcb), wx1j[jlm]), crs[mm])));
            hbc Ddm = hb_add(hb_scale(hb_mul(hb_mul(wx1i[ilm], ddb), wx0j[jlm]), mia[mm]),
                             hb_scale(hb_mul(hb_mul(wx0i[ilm], ddb), wx1j[jlm]), mja[mm]));
            Ddm = hb_add(Ddm, hb_imul(hb_scale(hb_mul(hb_mul(wx1i[ilm], ddb), wx1j[jlm]), crs[mm])));
            Kp[mm][idx] = hb_add(Ddm, hb_scale(Dm, cf));
        }

        /* on-site D^(0) diagonal (only same-site shell m==1 -> 0-based m==0) */
        if (m == 0 && ilm == jlm) {
            double c0 = ccdi[ilm + norb*0];
            double w0 = cuCreal(wx0i[ilm]); /* wx0/wx1 are real (stored cmplx(x,0)) */
            double w1 = cuCreal(wx1i[ilm]);
            K4v = hb_add(K4v, make_cuDoubleComplex(c0*(w0*w0 + w1*w1), 0.0));
            for (int mm = 0; mm < 3; ++mm)
                Kp[mm][idx] = hb_add(Kp[mm][idx],
                    make_cuDoubleComplex(c0*2.0*w0*w1*momi[mm], 0.0));
        }
        K4[idx] = K4v;
    }
    __syncthreads();

    /* hcpx cart2sph each K in place */
    hbc *slabs[4] = {K1, K2, K3, K4};
    for (int s = 0; s < 4; ++s) {
        for (int idx = tid; idx < nn2; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int b = 0; b < norb; ++b)
                acc = hb_add(acc, hb_mul(slabs[s][cmidx(i, b, norb)], v[cmidx(b, j, norb)]));
            tmp[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn2; idx += blockDim.x) {
            int i = idx % norb, j = idx / norb;
            hbc acc = make_cuDoubleComplex(0.0, 0.0);
            for (int a = 0; a < norb; ++a)
                acc = hb_add(acc, hb_mul(vc[cmidx(i, a, norb)], tmp[cmidx(a, j, norb)]));
            res[idx] = acc;
        }
        __syncthreads();
        for (int idx = tid; idx < nn2; idx += blockDim.x) slabs[s][idx] = res[idx];
        __syncthreads();
    }

    /* lambda-pack: hcc(1:norb,1:norb)=lam*(K4+K3); (off,off)=lam*(K4-K3);
     * (1:norb,off)=lam*(K1-iK2); (off,1:norb)=lam*(K1+iK2). CPU uses (ilm,jlm)
     * row/col directly (no j,i swap here -- build_ccor packs kcomp(:,:) as-is). */
    for (int idx = tid; idx < nn2; idx += blockDim.x) {
        int ilm = idx % norb, jlm = idx / norb;
        hbc k1 = K1[cmidx(ilm, jlm, norb)];
        hbc k2 = K2[cmidx(ilm, jlm, norb)];
        hbc k3 = K3[cmidx(ilm, jlm, norb)];
        hbc k4 = K4[cmidx(ilm, jlm, norb)];
        o[cmidx(ilm, jlm, nb)]             = hb_scale(hb_add(k4, k3), lam);
        o[cmidx(ilm+norb, jlm+norb, nb)]   = hb_scale(hb_sub(k4, k3), lam);
        o[cmidx(ilm, jlm+norb, nb)]        = hb_scale(hb_sub(k1, hb_imul(k2)), lam);
        o[cmidx(ilm+norb, jlm, nb)]        = hb_scale(hb_add(k1, hb_imul(k2)), lam);
    }
}

extern "C" int hambuild_gpu_launch_ccor(
    const double *d_sbar_re, const double *d_sdot_re, int norb_s, int nm_store,
    int ntot, const hbc *d_wx0, const hbc *d_wx1, const double *d_mom,
    const double *d_q_ss, double theta_ss, const double *d_ccd,
    const double *d_lambda, double avw, const hbc *d_v, const hbc *d_vc,
    const int *d_iz, const int *d_nn, const int *d_atlist, const double *d_cr,
    const int *d_shell, const int *d_ino, const double *d_vet, int ndi, int kk,
    hbc *d_out, int norb, int nb, int nn_max, int nsite, int lambda_stride) {
    int threads = 128;
    dim3 grid(nn_max, nsite);
    size_t shmem = (size_t)6 * norb * norb * sizeof(hbc);
    k_build_ccor<<<grid, threads, shmem>>>(
        d_sbar_re, d_sdot_re, norb_s, nm_store, ntot, d_wx0, d_wx1, d_mom,
        d_q_ss, theta_ss, d_ccd, d_lambda, avw, d_v, d_vc, d_iz, d_nn, d_atlist,
        d_cr, d_shell, d_ino, d_vet, ndi, kk, d_out, norb, nb, nn_max, nsite,
        lambda_stride);
    if (HB_DEVICE_SYNC() != HB_SUCCESS) return 1;
    if (HB_GET_LAST_ERR() != HB_SUCCESS) return 2;
    return 0;
}
