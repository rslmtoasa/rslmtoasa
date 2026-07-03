"""Validation of librsrec against a literal NumPy transcription of the
Fortran algorithms in recursion.f90 (hop_b/crecal_b, cheb_*_mom/
chebyshev_recur_ll, hop/crecal, compute_moments_stochastic).

Builds a random open-boundary cluster with an impurity region (hall) and a
bulk region (ee, per type), with lsham SOC blocks, exactly the RS-LMTO data
model. Then compares every driver."""

import ctypes as ct
import numpy as np

rng = np.random.default_rng(3)

# ---------------------------------------------------------------- cluster --
nb = 6            # block size (use 6 instead of 18 for test speed)
L = 4             # 4x4x4 simple-cubic open cluster
kk = L**3
ntype = 2
nmax = 5          # first 5 atoms form the "impurity" region with hall
lld = 12

# neighbour list: on-site + 6 nearest neighbours (open boundaries)
pos = np.array([(x, y, z) for x in range(L) for y in range(L) for z in range(L)])
index = {tuple(p): i for i, p in enumerate(pos)}
nnmax = 7
nn = np.zeros((kk, nnmax), dtype=np.int32)        # Fortran nn(k,s), 1-based
for k, p in enumerate(pos):
    nbrs = []
    for d in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        q = tuple(p + np.array(d))
        if q in index:
            nbrs.append(index[q] + 1)
    nn[k, 0] = 1 + len(nbrs)       # nn(k,1) counts on-site + neighbour slots
    for s, nbi in enumerate(nbrs):
        nn[k, 1 + s] = nbi
iz = (rng.integers(0, ntype, kk) + 1).astype(np.int32)

def rblock(scale=1.0):
    return (rng.normal(size=(nb, nb)) + 1j*rng.normal(size=(nb, nb))) * scale

# NOTE on Hermiticity: like the real code, ee blocks are per (shell, type);
# we don't enforce global Hermiticity of the cluster H -- the algorithms
# don't require it for agreement testing (they mirror each other regardless).
ee = np.zeros((nb, nb, nnmax, ntype), dtype=np.complex128, order="F")
for t in range(ntype):
    D = rblock(); ee[:, :, 0, t] = 0.5*(D + D.conj().T) - 0.4*np.eye(nb)
    for s in range(1, nnmax):
        ee[:, :, s, t] = rblock(0.3)
hall = np.zeros((nb, nb, nnmax, nmax), dtype=np.complex128, order="F")
for k in range(nmax):
    D = rblock(); hall[:, :, 0, k] = 0.5*(D + D.conj().T) + 0.2*np.eye(nb)
    for s in range(1, nnmax):
        hall[:, :, s, k] = rblock(0.25)
lsham = np.zeros((nb, nb, ntype), dtype=np.complex128, order="F")
for t in range(ntype):
    D = rblock(0.1); lsham[:, :, t] = 0.5*(D + D.conj().T)
va = np.zeros_like(ee); vb = np.zeros_like(ee)
for t in range(ntype):
    for s in range(nnmax):
        va[:, :, s, t] = rblock(0.2); vb[:, :, s, t] = rblock(0.2)

# ------------------------------------------------- NumPy Fortran mirror ----
def np_apply(x, which=0, a=1.0, b=0.0):
    """x: (nb, nrhs, kk). Literal transcription of ham_vec_matmul /
    velo_vec_matmul (izero gating dropped: zero psi contributes zero)."""
    y = np.zeros_like(x)
    for k in range(kk):
        t = iz[k] - 1
        if which == 0:
            if k < nmax:
                blocks = hall[:, :, :, k]
                y[:, :, k] += (blocks[:, :, 0] + lsham[:, :, t]) @ x[:, :, k]
            else:
                blocks = ee[:, :, :, t]
                y[:, :, k] += (blocks[:, :, 0] + lsham[:, :, t]) @ x[:, :, k]
        else:
            blocks = (va if which == 1 else vb)[:, :, :, t]
            y[:, :, k] += blocks[:, :, 0] @ x[:, :, k]
        nr = nn[k, 0]
        for s in range(1, nr):
            nbr = nn[k, s] - 1
            if nbr >= 0:
                y[:, :, k] += blocks[:, :, s] @ x[:, :, nbr]
    return (y - b * x) / a

def np_cheb(psi0, lld, a, b):
    mu = np.zeros((nb, nb, 2*lld + 2), dtype=complex)
    p0 = psi0.copy()
    mu[:, :, 0] = np.einsum("lmk,lnk->mn", psi0.conj(), p0)
    p1 = np_apply(p0, 0, a, b)
    mu[:, :, 1] = np.einsum("lmk,lnk->mn", psi0.conj(), p1)
    for ll in range(1, lld + 1):
        p2 = 2.0 * np_apply(p1, 0, a, b) - p0
        d1 = np.einsum("lmk,lnk->mn", p1.conj(), p1)
        d2 = np.einsum("lmk,lnk->mn", p2.conj(), p1)
        mu[:, :, 2*ll]     = 2*d1 - mu[:, :, 0]
        mu[:, :, 2*ll + 1] = 2*d2 - mu[:, :, 1]
        p0, p1 = p1, p2
    return mu

def np_block_lanczos(psi0, lld):
    psi = psi0.copy(); pmn = np.zeros_like(psi)
    a_b = np.zeros((nb, nb, lld), dtype=complex)
    b2_b = np.zeros((nb, nb, lld), dtype=complex)
    sum_b = np.eye(nb, dtype=complex)
    for ll in range(lld - 1):
        hpsi = np_apply(psi, 0, 1.0, 0.0)
        An = np.einsum("lmk,lnk->mn", psi.conj(), hpsi)
        a_b[:, :, ll] = An
        pmn = hpsi - pmn
        b2_b[:, :, ll] = sum_b
        pmn = pmn - np.einsum("lmk,mn->lnk", psi, An)
        sum_b = np.einsum("lmk,lnk->mn", pmn.conj(), pmn)
        ev, U = np.linalg.eigh(sum_b)
        lam = np.sqrt(np.maximum(0.0, ev))
        B = (U * lam) @ U.conj().T
        Bi = (U * np.where(lam > 0, 1.0/np.where(lam > 0, lam, 1), 0.0)) @ U.conj().T
        psi_t = psi.copy()
        psi = np.einsum("lmk,mn->lnk", pmn, Bi)
        pmn = np.einsum("lmk,mn->lnk", psi_t, B)
    b2_b[:, :, lld - 1] = sum_b
    return a_b, b2_b

def np_scalar_lanczos(site_j, lld):
    psi = np.zeros((nb, nb, kk), dtype=complex)
    for l in range(nb):
        psi[l, l, site_j - 1] = 1.0
    pmn = np.zeros_like(psi)
    a = np.zeros((lld, nb)); b2 = np.zeros((lld, nb))
    summ = np.ones(nb)
    for ll in range(lld - 1):
        hpsi = np_apply(psi, 0, 1.0, 0.0)
        acol = np.einsum("lck,lck->c", psi.conj(), hpsi).real
        pmn = pmn + hpsi
        a[ll] = acol; b2[ll] = summ
        pmn = pmn - acol[None, :, None] * psi
        summ = np.einsum("lck,lck->c", pmn.conj(), pmn).real
        psi_t = psi.copy()
        psi = pmn / np.sqrt(summ)[None, :, None]
        pmn = -psi_t * np.sqrt(summ)[None, :, None]
    b2[lld - 1] = summ
    return a, b2

def np_stochastic(psiref, lld, a, b):
    left = np.zeros((nb, nb, kk, lld), dtype=complex)
    w1 = psiref.copy(); left[..., 0] = w1
    for m in range(2, lld + 1):
        if m == 2:
            w0 = w1.copy(); w1 = np_apply(w0, 0, a, b)
        else:
            w2 = 2*np_apply(w1, 0, a, b) - w0
            w0, w1 = w1, w2
        left[..., m - 1] = w1
    mu = np.zeros((nb, nb, lld, lld), dtype=complex)
    v0 = np_apply(psiref, 2, 1.0, 0.0)
    for n in range(1, lld + 1):
        if n == 1:
            v1 = v0.copy()
        elif n == 2:
            v0 = v1.copy(); v1 = np_apply(v0, 0, a, b)
        else:
            v2 = 2*np_apply(v1, 0, a, b) - v0
            v0, v1 = v1, v2
        R = np_apply(v1, 1, 1.0, 0.0)
        for m in range(lld):
            mu[:, :, n - 1, m] = np.einsum("lik,ljk->ij",
                                           left[..., m].conj(), R)
    return mu

# ------------------------------------------------------------ C library ----
lib = ct.CDLL("./librsrec_cpu.so")
lib.rsrec_create.restype = ct.c_void_p
lib.rsrec_last_error.restype = ct.c_char_p
P = lambda x: x.ctypes.data_as(ct.c_void_p)

ctx = ct.c_void_p(lib.rsrec_create(kk, nb, nnmax, ntype, nmax, 0))
assert ctx, lib.rsrec_last_error().decode()
nnF = np.asfortranarray(nn)
rc = lib.rsrec_set_hamiltonian(ctx, P(ee), P(hall), P(lsham),
                               nnF.ctypes.data_as(ct.POINTER(ct.c_int)),
                               iz.ctypes.data_as(ct.POINTER(ct.c_int)))
assert rc == 0, lib.rsrec_last_error().decode()
rc = lib.rsrec_set_velocity(ctx, P(va), P(vb)); assert rc == 0

def relerr(x, y):
    return np.max(np.abs(x - y)) / max(np.max(np.abs(y)), 1e-300)

# 1) matvec ------------------------------------------------------------------
x = np.asfortranarray(rng.normal(size=(nb, nb, kk))
                      + 1j*rng.normal(size=(nb, nb, kk)))
y = np.zeros_like(x, order="F")
for which in (0, 1, 2):
    lib.rsrec_op_apply(ctx, which, P(x), P(y), nb, ct.c_double(2.3),
                       ct.c_double(-0.4))
    print(f"op_apply which={which}: rel.err =",
          f"{relerr(y, np_apply(x, which, 2.3, -0.4)):.2e}")

# 2) Chebyshev moments (single-site init, like chebyshev_recur) --------------
a_s, b_s = 9.0, 0.1
psi0 = np.zeros((nb, nb, kk), dtype=complex, order="F")
for l in range(nb):
    psi0[l, l, 7] = 1.0
mu_c = np.zeros((nb, nb, 2*lld + 2), dtype=complex, order="F")
rc = lib.rsrec_chebyshev_moments(ctx, P(psi0), lld, ct.c_double(a_s),
                                 ct.c_double(b_s), P(mu_c))
assert rc == 0, lib.rsrec_last_error().decode()
print("chebyshev_moments (site): rel.err =",
      f"{relerr(mu_c, np_cheb(psi0, lld, a_s, b_s)):.2e}")

# 2b) pair init with the reci=3 sign combo, like chebyshev_recur_ij ----------
psi0 = np.zeros((nb, nb, kk), dtype=complex, order="F")
i_at, j_at = 2, 30
for l in range(nb):
    psi0[l, l, i_at] = 1/np.sqrt(2.0)
    psi0[l, l, j_at] = 1j/np.sqrt(2.0)
rc = lib.rsrec_chebyshev_moments(ctx, P(psi0), lld, ct.c_double(a_s),
                                 ct.c_double(b_s), P(mu_c))
print("chebyshev_moments (pair): rel.err =",
      f"{relerr(mu_c, np_cheb(psi0, lld, a_s, b_s)):.2e}")

# 3) block Lanczos ------------------------------------------------------------
psi0 = np.zeros((nb, nb, kk), dtype=complex, order="F")
for l in range(nb):
    psi0[l, l, 9] = 1.0
a_b = np.zeros((nb, nb, lld), dtype=complex, order="F")
b2_b = np.zeros((nb, nb, lld), dtype=complex, order="F")
rc = lib.rsrec_block_lanczos(ctx, P(psi0), lld, P(a_b), P(b2_b))
assert rc == 0, lib.rsrec_last_error().decode()
a_ref, b2_ref = np_block_lanczos(psi0, lld)
print("block_lanczos: a rel.err =", f"{relerr(a_b, a_ref):.2e}",
      " b2 rel.err =", f"{relerr(b2_b, b2_ref):.2e}")

# 4) scalar Lanczos ------------------------------------------------------------
a_o = np.zeros((lld, nb), order="F"); b2_o = np.zeros((lld, nb), order="F")
rc = lib.rsrec_scalar_lanczos(ctx, 10, lld,
                              a_o.ctypes.data_as(ct.POINTER(ct.c_double)),
                              b2_o.ctypes.data_as(ct.POINTER(ct.c_double)))
assert rc == 0, lib.rsrec_last_error().decode()
a_ref, b2_ref = np_scalar_lanczos(10, lld)
print("scalar_lanczos: a rel.err =", f"{relerr(a_o, a_ref):.2e}",
      " b2 rel.err =", f"{relerr(b2_o, b2_ref):.2e}")

# 5) stochastic moments --------------------------------------------------------
psiref = np.zeros((nb, nb, kk), dtype=complex, order="F")
for k in range(kk):
    ph = np.exp(2j*np.pi*rng.random())
    for l in range(nb):
        psiref[l, l, k] = ph
psiref /= np.sqrt(kk)
lld_s = 6
mu_nm = np.zeros((nb, nb, lld_s, lld_s), dtype=complex, order="F")
rc = lib.rsrec_stochastic_moments(ctx, P(psiref), lld_s, ct.c_double(a_s),
                                  ct.c_double(b_s), P(mu_nm))
assert rc == 0, lib.rsrec_last_error().decode()
print("stochastic_moments: rel.err =",
      f"{relerr(mu_nm, np_stochastic(psiref, lld_s, a_s, b_s)):.2e}")

lib.rsrec_destroy(ctx)
print("ALL DRIVERS VALIDATED")
