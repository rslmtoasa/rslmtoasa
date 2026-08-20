"""Validate the public RS recursion CUDA ABI against NumPy references.

This is an executable, hardware-side contract test.  It deliberately binds
every C argument with ``ctypes``: an ABI mismatch must fail with a useful
error instead of turning omitted arguments into garbage pointers and a
segmentation fault.  The test covers matvec, all recursion drivers, orbital
moments, and the four DOS/GF reconstruction entry points.
"""

from __future__ import annotations

import argparse
import ctypes as ct
from pathlib import Path

import numpy as np


rng = np.random.default_rng(3)
nb, L, kk, ntype, nmax, lld = 6, 4, 4**3, 2, 5, 4
nnmax = 7

pos = np.array([(x, y, z) for x in range(L) for y in range(L) for z in range(L)])
index = {tuple(p): i for i, p in enumerate(pos)}
nn = np.zeros((kk, nnmax), dtype=np.int32)
for k, p in enumerate(pos):
    nbrs = []
    for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        q = tuple(p + np.array(d))
        if q in index:
            nbrs.append(index[q] + 1)
    nn[k, 0] = 1 + len(nbrs)
    for s, nbi in enumerate(nbrs):
        nn[k, s + 1] = nbi
iz = (rng.integers(0, ntype, kk) + 1).astype(np.int32)


def rblock(scale: float = 1.0) -> np.ndarray:
    return (rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))) * scale


ee = np.zeros((nb, nb, nnmax, ntype), dtype=np.complex128, order="F")
hall = np.zeros((nb, nb, nnmax, nmax), dtype=np.complex128, order="F")
lsham = np.zeros((nb, nb, ntype), dtype=np.complex128, order="F")
va = np.zeros_like(ee)
vb = np.zeros_like(ee)
for t in range(ntype):
    d = rblock()
    ee[:, :, 0, t] = 0.5 * (d + d.conj().T) - 0.4 * np.eye(nb)
    d = rblock(0.1)
    lsham[:, :, t] = 0.5 * (d + d.conj().T)
    for s in range(1, nnmax):
        ee[:, :, s, t] = rblock(0.3)
        va[:, :, s, t] = rblock(0.2)
        vb[:, :, s, t] = rblock(0.2)
for k in range(nmax):
    d = rblock()
    hall[:, :, 0, k] = 0.5 * (d + d.conj().T) + 0.2 * np.eye(nb)
    for s in range(1, nnmax):
        hall[:, :, s, k] = rblock(0.25)


def np_apply(x: np.ndarray, which: int = 0, a: float = 1.0, b: float = 0.0) -> np.ndarray:
    y = np.zeros_like(x)
    for k in range(kk):
        t = iz[k] - 1
        if which == 0:
            blocks = hall[:, :, :, k] if k < nmax else ee[:, :, :, t]
            y[:, :, k] += (blocks[:, :, 0] + lsham[:, :, t]) @ x[:, :, k]
        else:
            blocks = (va if which == 1 else vb)[:, :, :, t]
            y[:, :, k] += blocks[:, :, 0] @ x[:, :, k]
        for s in range(1, nn[k, 0]):
            nbr = nn[k, s] - 1
            if nbr >= 0:
                y[:, :, k] += blocks[:, :, s] @ x[:, :, nbr]
    return (y - b * x) / a


def np_cheb(psi0: np.ndarray, steps: int, a: float, b: float) -> np.ndarray:
    mu = np.zeros((nb, nb, 2 * steps + 2), dtype=complex)
    p0 = psi0.copy()
    mu[:, :, 0] = np.einsum("lmk,lnk->mn", psi0.conj(), p0)
    p1 = np_apply(p0, 0, a, b)
    mu[:, :, 1] = np.einsum("lmk,lnk->mn", psi0.conj(), p1)
    for ll in range(1, steps + 1):
        p2 = 2 * np_apply(p1, 0, a, b) - p0
        mu[:, :, 2 * ll] = 2 * np.einsum("lmk,lnk->mn", p1.conj(), p1) - mu[:, :, 0]
        mu[:, :, 2 * ll + 1] = 2 * np.einsum("lmk,lnk->mn", p2.conj(), p1) - mu[:, :, 1]
        p0, p1 = p1, p2
    return mu


def np_orbital(left: np.ndarray, psi: np.ndarray, steps: int, a: float, b: float) -> np.ndarray:
    out = np.zeros((nb, nb, steps), dtype=complex)
    p0 = psi.copy()
    out[:, :, 0] = np.einsum("lmk,lnk->mn", left.conj(), p0)
    if steps == 1:
        return out
    p1 = np_apply(p0, 0, a, b)
    out[:, :, 1] = np.einsum("lmk,lnk->mn", left.conj(), p1)
    for n in range(2, steps):
        p2 = 2 * np_apply(p1, 0, a, b) - p0
        out[:, :, n] = np.einsum("lmk,lnk->mn", left.conj(), p2)
        p0, p1 = p1, p2
    return out


def np_block_lanczos(psi0: np.ndarray, steps: int) -> tuple[np.ndarray, np.ndarray]:
    psi = psi0.copy()
    pmn = np.zeros_like(psi)
    aa = np.zeros((nb, nb, steps), dtype=complex)
    bb = np.zeros((nb, nb, steps), dtype=complex)
    summ = np.eye(nb, dtype=complex)
    for ll in range(steps - 1):
        hpsi = np_apply(psi)
        an = np.einsum("lmk,lnk->mn", psi.conj(), hpsi)
        aa[:, :, ll] = an
        pmn = hpsi - pmn - np.einsum("lmk,mn->lnk", psi, an)
        bb[:, :, ll] = summ
        summ = np.einsum("lmk,lnk->mn", pmn.conj(), pmn)
        ev, u = np.linalg.eigh(summ)
        lam = np.sqrt(np.maximum(0.0, ev))
        bmat = (u * lam) @ u.conj().T
        inv = (u * np.where(lam > 0, 1.0 / np.where(lam > 0, lam, 1), 0.0)) @ u.conj().T
        psi_old = psi.copy()
        psi = np.einsum("lmk,mn->lnk", pmn, inv)
        pmn = np.einsum("lmk,mn->lnk", psi_old, bmat)
    bb[:, :, steps - 1] = summ
    return aa, bb


def np_scalar_lanczos(site: int, steps: int) -> tuple[np.ndarray, np.ndarray]:
    psi = np.zeros((nb, nb, kk), dtype=complex)
    for l in range(nb):
        psi[l, l, site - 1] = 1
    pmn = np.zeros_like(psi)
    aa = np.zeros((steps, nb))
    bb = np.zeros((steps, nb))
    summ = np.ones(nb)
    for ll in range(steps - 1):
        hpsi = np_apply(psi)
        acol = np.einsum("lck,lck->c", psi.conj(), hpsi).real
        pmn = hpsi + pmn
        aa[ll] = acol
        bb[ll] = summ
        pmn = pmn - acol[None, :, None] * psi
        summ = np.einsum("lck,lck->c", pmn.conj(), pmn).real
        psi_old = psi.copy()
        psi = pmn / np.sqrt(summ)[None, :, None]
        pmn = -psi_old * np.sqrt(summ)[None, :, None]
    bb[steps - 1] = summ
    return aa, bb


def np_stochastic(psiref: np.ndarray, steps: int, a: float, b: float) -> np.ndarray:
    left = np.zeros((nb, nb, kk, steps), dtype=complex)
    w1 = psiref.copy()
    left[..., 0] = w1
    for m in range(2, steps + 1):
        if m == 2:
            w0, w1 = w1.copy(), np_apply(w1, 0, a, b)
        else:
            w0, w1 = w1, 2 * np_apply(w1, 0, a, b) - w0
        left[..., m - 1] = w1
    mu = np.zeros((nb, nb, steps, steps), dtype=complex)
    v0 = np_apply(psiref, 2)
    for n in range(1, steps + 1):
        if n == 1:
            v1 = v0.copy()
        elif n == 2:
            v0, v1 = v1.copy(), np_apply(v1, 0, a, b)
        else:
            v0, v1 = v1, 2 * np_apply(v1, 0, a, b) - v0
        right = np_apply(v1, 1)
        for m in range(steps):
            mu[:, :, n - 1, m] = np.einsum("lik,ljk->ij", left[..., m].conj(), right)
    return mu


def transport_reconstruction_reference(energies: np.ndarray, moments: int,
                                       a: float, b: float, factor: float,
                                       mu: np.ndarray) -> np.ndarray:
    """Reference for the G1.3 flattened Gamma*U contraction."""
    x = (energies - b) / a
    root = np.sqrt(1.0 - x * x)
    theta = np.arccos(x)
    cheb = np.empty((len(energies), moments), dtype=np.float64)
    cheb[:, 0] = 1.0
    if moments > 1:
        cheb[:, 1] = x
    for n in range(2, moments):
        cheb[:, n] = 2.0 * x * cheb[:, n - 1] - cheb[:, n - 2]

    kernel = np.sinh(6.0 * (1.0 - np.arange(moments) / moments)) / np.sinh(6.0)
    weights = np.ones(moments, dtype=np.float64)
    weights[0] = 0.5
    gamma = np.empty((len(energies), moments, moments), dtype=np.complex128)
    for n in range(moments):
        cn = (x - 1j * n * root) * np.exp(1j * n * theta)
        for m in range(moments):
            cm = (x + 1j * m * root) * np.exp(-1j * m * theta)
            gamma[:, n, m] = ((cn * cheb[:, m] + cm * cheb[:, n]) /
                              ((1.0 - x * x) ** 2) * kernel[n] * kernel[m] *
                              weights[n] * weights[m])

    packed = np.empty((moments * moments, mu.shape[1]), dtype=np.complex128,
                      order="F")
    for l in range(mu.shape[1]):
        for m in range(moments):
            for n in range(moments):
                packed[n + moments * m, l] = mu[l, l, n, m]
    gamma_matrix = gamma.reshape((len(energies), moments * moments), order="F")
    return factor * gamma_matrix @ packed


def jackson(n: int) -> np.ndarray:
    i = np.arange(n, dtype=float)
    theta = np.pi * i / (n + 1.0)
    return ((n - i + 1) * np.cos(theta) + np.sin(theta) * (1 / np.tan(np.pi / (n + 1)))) / (n + 1)


def cheb_transfer(n: int, energies: np.ndarray, a: float, b: float) -> np.ndarray:
    x = (energies - b) / a
    theta = np.arccos(x)
    pref = 1 / np.sqrt(a * a - (energies - b) ** 2)
    i = np.arange(n, dtype=float)[:, None]
    c = np.where(i == 0, 1.0, 2.0)
    return jackson(n)[:, None] * c * (-1j) * np.exp(-1j * i * theta[None, :]) * pref[None, :]


def block_reference(a_b: np.ndarray, b2_b: np.ndarray, a_inf: np.ndarray,
                    b_inf: np.ndarray, energies: np.ndarray, eta_re: np.ndarray,
                    eta_im: np.ndarray, sym: int) -> np.ndarray:
    steps, natoms = a_b.shape[2], a_b.shape[3]
    out = np.zeros((nb, nb, len(energies), natoms), dtype=complex, order="F")
    for n in range(natoms):
        ad = 0.5 * (a_inf[0, n] + a_inf[9, n]) if nb >= 10 else a_inf[0, n]
        bd = 0.5 * (b_inf[0, n] + b_inf[9, n]) if nb >= 10 else b_inf[0, n]
        for e, energy in enumerate(energies):
            q = np.zeros((nb, nb), dtype=complex)
            for i in range(nb):
                ai = ad if sym else a_inf[i, n]
                bi = bd if sym else b_inf[i, n]
                f = 1.025 if (not sym and i == 0) else 1.0
                det = (energy - ai - 2 * bi * f) * (energy - ai + 2 * bi * f)
                root = np.sqrt(complex(det, 0))
                q[i, i] = 0.5 * ((energy + eta_re[e] - ai - root.real)
                                 + 1j * (eta_im[e] - root.imag))
            for lc in range(steps - 1, 0, -1):
                qin = (energy + eta_re[e] + 1j * eta_im[e]) * np.eye(nb) - a_b[:, :, lc - 1, n] - q
                g = np.linalg.inv(qin)
                q = b2_b[:, :, lc - 1, n].conj().T @ g @ b2_b[:, :, lc - 1, n]
            out[:, :, e, n] = q
    return out


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.max(np.abs(actual - expected)) / max(float(np.max(np.abs(expected))), 1e-300))


def check(name: str, actual: np.ndarray, expected: np.ndarray, tol: float = 2e-8) -> None:
    err = relerr(actual, expected)
    print(f"{name}: rel.err = {err:.2e} (tol {tol:.1e})")
    if not np.isfinite(err) or err > tol:
        raise AssertionError(f"{name} exceeds tolerance: {err:.3e}")


def ptr(a: np.ndarray) -> ct.c_void_p:
    return a.ctypes.data_as(ct.c_void_p)


def iptr(a: np.ndarray):
    return a.ctypes.data_as(ct.POINTER(ct.c_int))


def dptr(a: np.ndarray):
    return a.ctypes.data_as(ct.POINTER(ct.c_double))


def configure(lib: ct.CDLL) -> None:
    v, i, d = ct.c_void_p, ct.c_int, ct.c_double
    pi, pd = ct.POINTER(i), ct.POINTER(d)
    lib.rsrec_create.argtypes = [i, i, i, i, i, i]
    lib.rsrec_create.restype = v
    lib.rsrec_destroy.argtypes = [v]
    lib.rsrec_last_error.restype = ct.c_char_p
    lib.rsrec_set_hamiltonian.argtypes = [v, v, v, v, pi, pi, v, v, v]
    lib.rsrec_set_hamiltonian.restype = i
    lib.rsrec_set_velocity.argtypes = [v, v, v, v, v]
    lib.rsrec_set_velocity.restype = i
    lib.rsrec_set_precision.argtypes = [v, i]
    lib.rsrec_set_precision.restype = i
    lib.rsrec_op_apply.argtypes = [v, i, v, v, i, d, d]
    lib.rsrec_op_apply.restype = i
    lib.rsrec_chebyshev_moments.argtypes = [v, v, i, d, d, v]
    lib.rsrec_chebyshev_moments.restype = i
    lib.rsrec_block_lanczos.argtypes = [v, v, i, v, v, i]
    lib.rsrec_block_lanczos.restype = i
    lib.rsrec_scalar_lanczos.argtypes = [v, i, i, pd, pd]
    lib.rsrec_scalar_lanczos.restype = i
    lib.rsrec_stochastic_moments.argtypes = [v, v, i, d, d, v]
    lib.rsrec_stochastic_moments.restype = i
    lib.rsrec_stochastic_moments_resident.argtypes = [v, v, i, d, d, i, i, v]
    lib.rsrec_stochastic_moments_resident.restype = i
    lib.rsrec_clear_resident_moments.argtypes = [v]
    lib.rsrec_clear_resident_moments.restype = i
    lib.rsrec_resident_count.argtypes = [v, pi]
    lib.rsrec_resident_count.restype = i
    ll = ct.POINTER(ct.c_longlong)
    lib.rsrec_reconstruct_conductivity.argtypes = [
        v, pd, i, i, d, d, d, i, i, v, pd, pd, pd, pd, pd,
        ll, ll, ll, pi,
    ]
    lib.rsrec_reconstruct_conductivity.restype = i
    lib.rsrec_stochastic_profile.argtypes = [v, pd, pd, pd, pd, pd, ll, ll, ll]
    lib.rsrec_stochastic_profile.restype = i
    lib.rsrec_orbital_moments.argtypes = [v, v, v, i, d, d, v]
    lib.rsrec_orbital_moments.restype = i
    lib.rsrec_chebyshev_dos.argtypes = [v, v, i, i, pd, i, d, d, v]
    lib.rsrec_chebyshev_dos.restype = i
    lib.rsrec_chebyshev_gf_eta.argtypes = [v, v, i, i, v, i, v]
    lib.rsrec_chebyshev_gf_eta.restype = i
    lib.rsrec_block_dos.argtypes = [v, v, v, pd, pd, pd, i, d, d, i, i, i, v]
    lib.rsrec_block_dos.restype = i
    lib.rsrec_block_gf_eta.argtypes = [v, v, v, pd, pd, d, pd, pd, i, i, i, i, v]
    lib.rsrec_block_gf_eta.restype = i


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", default="./librsrec_gpu.so", type=Path)
    args = parser.parse_args()
    lib = ct.CDLL(str(args.library.resolve()))
    configure(lib)

    def error() -> str:
        msg = lib.rsrec_last_error()
        return msg.decode() if msg else "unknown rsrec error"

    def call(name: str, rc: int) -> None:
        if rc:
            raise RuntimeError(f"{name} failed: {error()}")

    nn_f = np.asfortranarray(nn)
    ctx = lib.rsrec_create(kk, nb, nnmax, ntype, nmax, 0)
    if not ctx:
        raise RuntimeError(error())
    try:
        call("set_hamiltonian", lib.rsrec_set_hamiltonian(
            ctx, ptr(ee), ptr(hall), ptr(lsham), iptr(nn_f), iptr(iz), None, None, None))
        call("set_velocity", lib.rsrec_set_velocity(ctx, ptr(va), ptr(vb), None, None))
        call("set_precision", lib.rsrec_set_precision(ctx, 1))

        x = np.asfortranarray(rng.normal(size=(nb, nb, kk)) + 1j * rng.normal(size=(nb, nb, kk)))
        y = np.zeros_like(x, order="F")
        for which in (0, 1, 2):
            call("op_apply", lib.rsrec_op_apply(ctx, which, ptr(x), ptr(y), nb, 2.3, -0.4))
            check(f"op_apply which={which}", y, np_apply(x, which, 2.3, -0.4))

        a_s, b_s = 9.0, 0.1
        psi = np.zeros((nb, nb, kk), dtype=complex, order="F")
        for l in range(nb):
            psi[l, l, 7] = 1
        mu = np.zeros((nb, nb, 2 * lld + 2), dtype=complex, order="F")
        call("chebyshev_moments(site)", lib.rsrec_chebyshev_moments(ctx, ptr(psi), lld, a_s, b_s, ptr(mu)))
        check("chebyshev_moments site", mu, np_cheb(psi, lld, a_s, b_s))
        psi[:, :, :] = 0
        for l in range(nb):
            psi[l, l, 2] = 1 / np.sqrt(2)
            psi[l, l, 30] = 1j / np.sqrt(2)
        call("chebyshev_moments(pair)", lib.rsrec_chebyshev_moments(ctx, ptr(psi), lld, a_s, b_s, ptr(mu)))
        check("chebyshev_moments pair", mu, np_cheb(psi, lld, a_s, b_s))

        psi[:, :, :] = 0
        for l in range(nb):
            psi[l, l, 9] = 1
        aa = np.zeros((nb, nb, lld), dtype=complex, order="F")
        bb = np.zeros_like(aa)
        call("block_lanczos", lib.rsrec_block_lanczos(ctx, ptr(psi), lld, ptr(aa), ptr(bb), 1))
        aa_ref, bb_ref = np_block_lanczos(psi, lld)
        check("block_lanczos a", aa, aa_ref)
        check("block_lanczos b2", bb, bb_ref)

        ao = np.zeros((lld, nb), dtype=np.float64, order="F")
        bo = np.zeros_like(ao)
        call("scalar_lanczos", lib.rsrec_scalar_lanczos(ctx, 10, lld, dptr(ao), dptr(bo)))
        ao_ref, bo_ref = np_scalar_lanczos(10, lld)
        check("scalar_lanczos a", ao, ao_ref)
        check("scalar_lanczos b2", bo, bo_ref)

        psiref = np.zeros((nb, nb, kk), dtype=complex, order="F")
        for k in range(kk):
            psiref[:, :, k] = np.eye(nb) * np.exp(2j * np.pi * rng.random()) / np.sqrt(kk)
        steps_s = 6
        mu_nm = np.zeros((nb, nb, steps_s, steps_s), dtype=complex, order="F")
        call("stochastic_moments", lib.rsrec_stochastic_moments(ctx, ptr(psiref), steps_s, a_s, b_s, ptr(mu_nm)))
        check("stochastic_moments", mu_nm, np_stochastic(psiref, steps_s, a_s, b_s), 5e-8)
        # A second request reuses the context-lifetime left/moment/recurrence
        # workspace.  The result must remain independent of the reuse.
        mu_repeat = np.zeros_like(mu_nm)
        call("stochastic_moments(reuse)", lib.rsrec_stochastic_moments(
            ctx, ptr(psiref), steps_s, a_s, b_s, ptr(mu_repeat)))
        check("stochastic_moments reuse", mu_repeat, mu_nm, 5e-8)

        # Exercise the explicit production FP32 route against the same
        # deterministic reference vector and physics.
        call("set_precision(fp32)", lib.rsrec_set_precision(ctx, 0))
        mu_fp32 = np.zeros_like(mu_nm)
        call("stochastic_moments(fp32)", lib.rsrec_stochastic_moments(
            ctx, ptr(psiref), steps_s, a_s, b_s, ptr(mu_fp32)))
        check("stochastic_moments fp32", mu_fp32,
              np_stochastic(psiref, steps_s, a_s, b_s), 5e-5)
        call("set_precision(fp64)", lib.rsrec_set_precision(ctx, 1))

        # G1.3: retain only packed diagonal moments, reconstruct with a
        # non-divisible energy tile, and prove that the optimized request did
        # not perform the full-moment D2H transfer.
        energies_g13 = np.linspace(-2.0, 2.0, 5).astype(np.float64)
        energies_g13 = np.ascontiguousarray(energies_g13)
        c64 = np.zeros((len(energies_g13), nb, 1), dtype=np.complex128, order="F")
        resident_mu = np.zeros_like(mu_nm)
        call("stochastic_moments_resident(fp64 diagnostic)", lib.rsrec_stochastic_moments_resident(
            ctx, ptr(psiref), steps_s, a_s, b_s, 1, 1, ptr(resident_mu)))
        check("resident moments fp64", resident_mu, mu_nm, 5e-8)
        call("stochastic_moments_resident(fp64)", lib.rsrec_stochastic_moments_resident(
            ctx, ptr(psiref), steps_s, a_s, b_s, 1, 0, ptr(resident_mu)))
        resident_count = np.zeros(1, dtype=np.int32)
        call("resident_count", lib.rsrec_resident_count(ctx, iptr(resident_count)))
        if int(resident_count[0]) != 1:
            raise AssertionError(f"expected one resident trace, got {resident_count[0]}")
        h2d_s = np.zeros(1, dtype=np.float64)
        cheb_s = np.zeros(1, dtype=np.float64)
        d2h_s = np.zeros(1, dtype=np.float64)
        conversion_s = np.zeros(1, dtype=np.float64)
        pack_s = np.zeros(1, dtype=np.float64)
        h2d_b = np.zeros(1, dtype=np.int64)
        d2h_b = np.zeros(1, dtype=np.int64)
        pack_b = np.zeros(1, dtype=np.int64)
        ll = ct.POINTER(ct.c_longlong)
        call("stochastic_profile(fp64 resident)", lib.rsrec_stochastic_profile(
            ctx, dptr(h2d_s), dptr(cheb_s), dptr(d2h_s), dptr(conversion_s),
            dptr(pack_s), h2d_b.ctypes.data_as(ll), d2h_b.ctypes.data_as(ll),
            pack_b.ctypes.data_as(ll)))
        if int(d2h_b[0]) != 0:
            raise AssertionError(f"resident route performed full-moment D2H: {d2h_b[0]} bytes")
        gamma_s = np.zeros(1, dtype=np.float64)
        gamma_basis_s = np.zeros(1, dtype=np.float64)
        gamma_fill_s = np.zeros(1, dtype=np.float64)
        gemm_s = np.zeros(1, dtype=np.float64)
        result_d2h_s = np.zeros(1, dtype=np.float64)
        gamma_h2d_b = np.zeros(1, dtype=np.int64)
        gamma_block_b = np.zeros(1, dtype=np.int64)
        result_d2h_b = np.zeros(1, dtype=np.int64)
        actual_be = np.zeros(1, dtype=np.int32)
        call("reconstruct_conductivity(fp64)", lib.rsrec_reconstruct_conductivity(
            ctx, dptr(energies_g13), len(energies_g13), steps_s, a_s, b_s, 1.7,
            1, 3, ptr(c64), dptr(gamma_s), dptr(gamma_basis_s), dptr(gamma_fill_s),
            dptr(gemm_s), dptr(result_d2h_s), gamma_h2d_b.ctypes.data_as(ll),
            gamma_block_b.ctypes.data_as(ll), result_d2h_b.ctypes.data_as(ll),
            iptr(actual_be)))
        g13_ref64 = transport_reconstruction_reference(energies_g13, steps_s, a_s, b_s,
                                                       1.7, mu_nm)
        check("G1.3 reconstruction fp64", c64[:, :, 0], g13_ref64, 2e-8)
        if int(actual_be[0]) != 3:
            raise AssertionError(f"requested BE=3, got {actual_be[0]}")
        print(f"G1.3 fp64: BE={actual_be[0]}, Gamma={gamma_s[0]:.4f}s, "
              f"GEMM={gemm_s[0]:.4f}s, result D2H={result_d2h_b[0]} bytes")
        call("clear_resident_moments(fp64)", lib.rsrec_clear_resident_moments(ctx))

        call("set_precision(fp32 reconstruction)", lib.rsrec_set_precision(ctx, 0))
        mu_fp32 = np.zeros_like(mu_nm)
        call("stochastic_moments(fp32 reference)", lib.rsrec_stochastic_moments(
            ctx, ptr(psiref), steps_s, a_s, b_s, ptr(mu_fp32)))
        c32 = np.zeros((len(energies_g13), nb, 1), dtype=np.complex64, order="F")
        call("stochastic_moments_resident(fp32)", lib.rsrec_stochastic_moments_resident(
            ctx, ptr(psiref), steps_s, a_s, b_s, 1, 0, ptr(resident_mu)))
        call("reconstruct_conductivity(fp32)", lib.rsrec_reconstruct_conductivity(
            ctx, dptr(energies_g13), len(energies_g13), steps_s, a_s, b_s, 1.7,
            1, 3, ptr(c32), dptr(gamma_s), dptr(gamma_basis_s), dptr(gamma_fill_s),
            dptr(gemm_s), dptr(result_d2h_s), gamma_h2d_b.ctypes.data_as(ll),
            gamma_block_b.ctypes.data_as(ll), result_d2h_b.ctypes.data_as(ll),
            iptr(actual_be)))
        check("G1.3 reconstruction fp32", c32[:, :, 0],
              transport_reconstruction_reference(energies_g13, steps_s, a_s, b_s,
                                                  1.7, mu_fp32), 5e-4)
        call("clear_resident_moments(fp32)", lib.rsrec_clear_resident_moments(ctx))
        call("set_precision(fp64 final)", lib.rsrec_set_precision(ctx, 1))

        left = np.asfortranarray(rng.normal(size=psiref.shape) + 1j * rng.normal(size=psiref.shape))
        mu_o = np.zeros((nb, nb, lld), dtype=complex, order="F")
        call("orbital_moments", lib.rsrec_orbital_moments(ctx, ptr(left), ptr(psiref), lld, a_s, b_s, ptr(mu_o)))
        check("orbital_moments", mu_o, np_orbital(left, psiref, lld, a_s, b_s))

        # Reconstruction contracts use the exact output layouts documented in rsrec.h.
        natoms, nmom, nv = 2, 2 * lld + 2, 5
        mu_r = np.asfortranarray(rng.normal(size=(nb, nb, nmom, natoms)) + 1j * rng.normal(size=(nb, nb, nmom, natoms)))
        energies = np.linspace(-2.0, 2.0, nv)
        g_dos = np.zeros((nb, nb, nv, natoms), dtype=complex, order="F")
        energies_c = np.ascontiguousarray(energies, dtype=np.float64)
        call("chebyshev_dos", lib.rsrec_chebyshev_dos(ctx, ptr(mu_r), nmom, natoms, dptr(energies_c), nv, a_s, b_s, ptr(g_dos)))
        check("chebyshev_dos", g_dos, np.einsum("ijon,oe->ijen", mu_r, cheb_transfer(nmom, energies, a_s, b_s)), 1e-7)

        eta_re = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        eta_im = np.array([0.03, 0.07, 0.12], dtype=np.float64)
        ef = 0.25
        zenergies = ef + eta_re + 1j * eta_im
        ztransfer = np.empty((nmom, len(zenergies)), dtype=complex, order="F")
        for e, z in enumerate(zenergies):
            ztransfer[:, e] = cheb_transfer(nmom, np.array([z], dtype=complex), a_s, b_s)[:, 0]
        g_eta = np.zeros((nb, nb, len(zenergies), natoms), dtype=complex, order="F")
        call("chebyshev_gf_eta", lib.rsrec_chebyshev_gf_eta(ctx, ptr(mu_r), nmom, natoms, ptr(ztransfer), len(zenergies), ptr(g_eta)))
        check("chebyshev_gf_eta", g_eta, np.einsum("ijon,oe->ijen", mu_r, ztransfer), 1e-7)

        steps_b, nat_b = 3, 2
        ab = np.zeros((nb, nb, steps_b, nat_b), dtype=complex, order="F")
        b2 = np.zeros_like(ab)
        for n in range(nat_b):
            for l in range(steps_b):
                ab[:, :, l, n] = np.diag(0.1 * (l + 1) + 0.02 * np.arange(nb))
                b2[:, :, l, n] = np.diag(0.8 + 0.03 * np.arange(nb))
        a_inf = np.asfortranarray(0.1 * np.arange(nb)[:, None] + 0.01 * np.arange(nat_b)[None, :])
        b_inf = np.asfortranarray(0.9 + 0.02 * np.arange(nb)[:, None] + 0.01 * np.arange(nat_b)[None, :])
        energies_b = np.linspace(-1.0, 1.0, 4)
        g_block = np.zeros((nb, nb, len(energies_b), nat_b), dtype=complex, order="F")
        energies_b_c = np.ascontiguousarray(energies_b, dtype=np.float64)
        call("block_dos", lib.rsrec_block_dos(ctx, ptr(ab), ptr(b2), dptr(a_inf), dptr(b_inf), dptr(energies_b_c), len(energies_b), 0.0, 0.05, nat_b, steps_b, 0, ptr(g_block)))
        check("block_dos", g_block, block_reference(ab, b2, a_inf, b_inf, energies_b, np.zeros(4), np.full(4, 0.05), 0), 1e-7)
        g_block_eta = np.zeros((nb, nb, len(eta_im), nat_b), dtype=complex, order="F")
        call("block_gf_eta", lib.rsrec_block_gf_eta(ctx, ptr(ab), ptr(b2), dptr(a_inf), dptr(b_inf), ef, dptr(eta_re), dptr(eta_im), len(eta_im), nat_b, steps_b, 0, ptr(g_block_eta)))
        check("block_gf_eta", g_block_eta, block_reference(ab, b2, a_inf, b_inf, np.full(3, ef), eta_re, eta_im, 0), 1e-7)
    finally:
        lib.rsrec_destroy(ctx)
    print("ALL LOW-LEVEL CUDA ROUTES VALIDATED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
