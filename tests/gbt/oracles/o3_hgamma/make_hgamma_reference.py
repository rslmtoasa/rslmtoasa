#!/usr/bin/env python3
"""O3 (draft, for developer review): exact-ASA H^gamma reference for bcc Fe on the fixed fe_bcc potential. numpy only.

One site, local frames as in the GBT restart spec (R1). Inputs are plain files exported by the code (see README.md):
  <export>/gbt_q0/gbt_bonds.out  S^alpha bonds (gbt_single_q run, q = 0, theta = 0)
  <export>/gbt_q0/Fe.nml, input.nml, lattice.nml
  <export>/nc_first, nc_second   periodic_nc spin-resolved eigenvalue exports (check P0 only)
Construction (Turek et al. 1997, Andersen & Jepsen 1984):
  S^gamma = S^alpha [1 - (gamma - alpha) S^alpha]^-1,   H^gamma = C + dele S^gamma dele,
  spiral: S_rot(k) = U^dag diag(S(k - q/2), S(k + q/2)) U in the 18x18 local frame, H_rot likewise with the local-frame
  diagonal gamma - alpha, C, dele.  Omega = Mermin grand potential at fixed mu (mu solved once at q = 0, theta = 0).
usage: make_hgamma_reference.py <export dir> <output .nml> [--workers N] [--skip-values] [--skip-checks]
"""
import argparse
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

import numpy as np

KB_RY_PER_K = 6.3336814e-6              # reciprocal_mod kB_Ry_per_K
ANG2AU = 1.8897259886                   # math.f90:92
ALPHA = np.array([0.3485, 0.0530, 0.0107])   # strux default screening, lattice.f90:48 (l = 0, 1, 2)
TEMPERATURE_K = 300.0
KT = TEMPERATURE_K * KB_RY_PER_K
N_ELECTRONS = 8.0
Q_LIST = [(0.125, -0.125, 0.125), (0.25, -0.25, 0.25)]       # direct
THETA_DEG = (0.25, 0.5, 1.0)
ROT_THETA_DEG = (5.0, 30.0)                  # rotated eigenvalues (b) and the driver-level theta (c): 5 deg
C1_C2_TOL = 1e-10                            # "ok" thresholds set by this generator (not by the developer): round-off scale
MOMENT_MESHES = (4, 8, 16)
M_TOL = 1e-10                                # O3b self-checks, round-off scale (set by this generator)
C3_REL_TOL = 1e-4                            # tail-fit and theta^4 truncation scale

# ------------------------------------------------------------------ input reading


def read_numbers(text):
    return [float(v.replace("D", "E").replace("d", "e")) for v in re.split(r"[,\s]+", text.strip()) if v]


def read_potential(path):
    par = Path(path).read_text().split("&par")[1].split("\n/")[0]
    out = {}
    for line in par.splitlines():
        m = re.match(r"\s*(\w+)\s*(?:\(([^)]*)\))?\s*=\s*(.*)", line)
        if not m:
            continue
        name, idx, val = m.groups()
        if idx is None:
            out[name] = read_numbers(val)[0] if len(read_numbers(val)) == 1 else np.array(read_numbers(val))
        elif idx.replace(" ", "").startswith(":,"):
            out.setdefault(name, [None, None])[int(idx.split(",")[1]) - 1] = np.array(read_numbers(val))
    return out


def read_bonds(path):
    d, s = [], []
    for line in Path(path).read_text().splitlines():
        f = line.split()
        if not f or f[0].startswith("#") or f[0] == "link":
            continue
        if f[0] == "bond":
            d.append([float(x) for x in f[f.index("d=") + 1:f.index("d=") + 4]])
            s.append(np.zeros((9, 9), complex))
        elif f[0] == "S":
            s[-1][int(f[1]) - 1, int(f[2]) - 1] = complex(float(f[3]), float(f[4]))
    return np.array(d), np.array(s)


def read_scalar(path, key):
    return float(re.search(rf"^\s*{key}\s*=\s*([-+0-9.eEdD]+)", Path(path).read_text(), re.M).group(1).replace("d", "e"))


def read_lattice(path):
    rows = [read_numbers(re.search(rf"a\(:,\s*{i}\)\s*=\s*(.*)", Path(path).read_text()).group(1)) for i in (1, 2, 3)]
    return np.array(rows).T                                           # columns are the lattice vectors, units of alat


def orbital_l(d_cart, s_bonds):
    """l of each of the 9 orbital indices from S blocks only: S_ij(-d) = (-1)^(l_i + l_j) S_ij(d). The index with a unique
    onsite S value is s (even); parity relative to s then separates p (odd) from d (even)."""
    ons = np.real(np.diag(s_bonds[int(np.argmin(np.linalg.norm(d_cart, axis=1)))]))
    uniq = [i for i in range(9) if np.sum(np.abs(ons - ons[i]) < 1e-8) == 1]
    assert len(uniq) == 1, f"s not unique in the onsite block: {ons}"
    parity = {uniq[0]: 1}
    changed = True
    while changed and len(parity) < 9:
        changed = False
        for m in range(len(d_cart)):
            mm = [x for x in range(len(d_cart)) if np.allclose(d_cart[x], -d_cart[m], atol=1e-8)][0]
            for i in range(9):
                for j in range(9):
                    if i in parity and j not in parity and abs(s_bonds[m][i, j]) > 1e-6:
                        parity[j] = parity[i]*int(np.sign(np.real(s_bonds[mm][i, j]/s_bonds[m][i, j])))
                        changed = True
    assert len(parity) == 9
    p = np.array([parity[i] for i in range(9)])
    l = np.where(p < 0, 1, 2)
    l[uniq[0]] = 0
    assert [int(np.sum(l == x)) for x in (0, 1, 2)] == [1, 3, 5], "orbital count per l is not 1, 3, 5"
    return l, ons


class Model:
    def __init__(self, export):
        g = Path(export) / "gbt_q0"
        self.d_cart, self.s_bonds = read_bonds(g / "gbt_bonds.out")
        self.alat = read_scalar(g / "input.nml", "alat")
        wav = read_scalar(g / "input.nml", "wav")
        self.a_cols = read_lattice(g / "lattice.nml")
        self.n_direct = np.linalg.solve(self.a_cols, (self.d_cart/self.alat).T).T
        assert np.allclose(self.n_direct, np.rint(self.n_direct), atol=1e-6), "bond vectors are not lattice vectors"
        self.lmap, self.onsite_s = orbital_l(self.d_cart, self.s_bonds)
        pot = read_potential(g / "Fe.nml")
        self.pot = pot
        wow = wav*ANG2AU/pot["ws_r"]
        l, i = self.lmap, self.lmap + 1
        sel = lambda name: np.array([pot[name][s][l] for s in (0, 1)])         # (spin, orbital)
        self.c, self.enu, self.srdel, self.qpar = sel("c"), sel("enu"), sel("srdel"), sel("qpar")
        self.vmad = pot["vmad"]
        self.dele = self.srdel*wow**(0.5 - i)
        self.gamma = self.qpar*wow**(1 - 2*i)
        self.cc = self.c + self.vmad                      # C of H^gamma
        self.enu_v = self.enu + self.vmad                 # E_nu
        self.gam_a = self.gamma - ALPHA[l]                # gamma - alpha
        x = 1 - self.gam_a*(self.c - self.enu)/self.dele**2
        self.center = (self.c - self.enu)*x + self.enu + self.vmad
        self.shifted = (self.c - self.enu)*x
        self.width = self.dele*x
        self.obar = self.gam_a/((self.c - self.enu)*self.gam_a - self.dele**2)
        self.stored_dev = max(np.abs(self.center - np.array([pot["center_band"][s][l] for s in (0, 1)])).max(),
                              np.abs(self.width - np.array([pot["width_band"][s][l] for s in (0, 1)])).max(),
                              np.abs(self.obar - np.array([pot["obar"][s][l] for s in (0, 1)])).max())

    def s_of_k(self, k):
        """S^alpha(k), k direct, shape (n, 3) -> (n, 9, 9); phase exp(+2 pi i k.d_direct) as in the code."""
        return np.einsum("nm,mij->nij", np.exp(2j*np.pi*(k @ self.n_direct.T)), self.s_bonds)

    # ---- collinear Hamiltonians, spin sigma = 0 (majority) or 1
    def h_first(self, s, sg):
        w = self.width[sg]
        return np.diag(self.center[sg])[None] + w[None, :, None]*s*w[None, None, :]

    def h_second(self, s, sg):
        w = self.width[sg]
        h = np.diag(self.shifted[sg])[None] + w[None, :, None]*s*w[None, None, :]
        return h - h*self.obar[sg][None, None, :] @ h + np.diag(self.enu_v[sg])[None]

    def h_gamma(self, s, sg):
        gam = self.gam_a[sg]
        sgam = s @ np.linalg.inv(np.eye(9)[None] - gam[None, :, None]*s)
        w = self.dele[sg]
        return np.diag(self.cc[sg])[None] + w[None, :, None]*sgam*w[None, None, :]

    # ---- rotated (18 x 18, local frame, index (s - 1)*9 + i)
    def s_rot(self, k, q, theta):
        u = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]])
        sk = [self.s_of_k(k - 0.5*np.asarray(q)), self.s_of_k(k + 0.5*np.asarray(q))]
        out = np.zeros((len(k), 18, 18), complex)
        for s in range(2):
            for t in range(2):
                out[:, 9*s:9*s + 9, 9*t:9*t + 9] = sum(np.conj(u[sg, s])*u[sg, t]*sk[sg] for sg in range(2))
        return out

    def loc(self, arr):
        return np.concatenate([arr[0], arr[1]])

    def h_rot(self, k, q, theta):
        sr = self.s_rot(k, q, theta)
        gam = self.loc(self.gam_a)
        sgam = sr @ np.linalg.inv(np.eye(18)[None] - gam[None, :, None]*sr)
        w = self.loc(self.dele)
        return np.diag(self.loc(self.cc))[None] + w[None, :, None]*sgam*w[None, None, :], sr

    def p_alpha(self, z, sg):
        """P^alpha_sigma(z) per orbital: (P^alpha)^-1 = (P^gamma)^-1 + gamma - alpha, P^gamma = (z - C)/dele^2."""
        d = z - self.cc[sg]
        return d/(self.dele[sg]**2 + self.gam_a[sg]*d)


MODEL = None


def init_worker(export):
    global MODEL
    MODEL = Model(export)


# ------------------------------------------------------------------ mesh work (parallel)


def mesh_k(n, lo, hi):
    idx = np.arange(lo, hi)
    return np.stack([idx // (n*n), (idx // n) % n, idx % n], axis=1)/n


def task_eigs(args):
    n, lo, hi, q, theta = args
    h, _ = MODEL.h_rot(mesh_k(n, lo, hi), q, theta)
    return np.linalg.eigvalsh(h)


def task_sums(args):
    n, lo, hi, q, theta, mu = args
    e = task_eigs((n, lo, hi, q, theta))
    x = (e - mu)/KT
    return float(np.sum(-KT*np.logaddexp(0.0, -x))), float(np.sum(0.5*(1.0 - np.tanh(0.5*x))))


def chunks(n, size=2048):
    return [(lo, min(lo + size, n**3)) for lo in range(0, n**3, size)]


def omega_probe(pool, n, q, theta, mu):
    res = pool.map(task_sums, [(n, lo, hi, q, theta, mu) for lo, hi in chunks(n)])
    return sum(r[0] for r in res)/n**3, sum(r[1] for r in res)/n**3


def solve_mu(pool, n):
    e = np.concatenate(pool.map(task_eigs, [(n, lo, hi, (0, 0, 0), 0.0) for lo, hi in chunks(n)]))
    lo_, hi_ = e.min() - 1.0, e.max() + 1.0
    for _ in range(80):
        mu = 0.5*(lo_ + hi_)
        if np.sum(0.5*(1.0 - np.tanh(0.5*(e - mu)/KT)))/n**3 < N_ELECTRONS:
            lo_ = mu
        else:
            hi_ = mu
    return mu


def neville0(x, f):
    f = list(f)
    for m in range(1, len(x)):
        f = [(x[i + m]*f[i] - x[i]*f[i + 1])/(x[i + m] - x[i]) for i in range(len(f) - 1)]
    return f[0]


def richardson_j(pool, n, mu):
    th = np.radians(THETA_DEG)
    x = np.sin(th)**2
    om0 = [omega_probe(pool, n, (0, 0, 0), t, mu)[0] for t in th]
    out = []
    for q in Q_LIST:
        f = [(omega_probe(pool, n, q, t, mu)[0] - o)/xx for t, o, xx in zip(th, om0, x)]
        out.append(neville0(x, f))
    return out, om0


# ------------------------------------------------------------------ LKAG Matsubara (alpha representation)


_S_CACHE = {}


def task_matsubara(args):
    n, iw, mu, qs = args
    if n not in _S_CACHE:
        _S_CACHE[n] = MODEL.s_of_k(mesh_k(n, 0, n**3))
    s = _S_CACHE[n]
    z = mu + 1j*(2*iw + 1)*np.pi*KT
    pu, pd = MODEL.p_alpha(z, 0), MODEL.p_alpha(z, 1)
    return mats_t(n, s, pu, pd, qs)


def mats_t(n, s, pu, pd, qs):
    au = np.linalg.inv(np.diag(pu)[None] - s).reshape(n, n, n, 9, 9)
    ad = np.linalg.inv(np.diag(pd)[None] - s).reshape(n, n, n, 9, 9)
    delta = pu - pd
    out = []
    for q in qs:
        sh = [int(round(0.5*qi*n)) for qi in q]
        assert all(abs(0.5*qi*n - sh_i) < 1e-9 for qi, sh_i in zip(q, sh)), "q/2 is not on the mesh"
        a = np.roll(au, sh, axis=(0, 1, 2))                      # g_up(k - q/2)
        b = np.roll(ad, [-v for v in sh], axis=(0, 1, 2))        # g_dn(k + q/2)
        out.append(np.einsum("i,xyzij,j,xyzji->", delta, a, delta, b)/n**3)
    return np.array(out)


def lkag_matsubara(pool, n, mu, cutoffs):
    qs = [(0, 0, 0)] + Q_LIST
    nmax = int(max(cutoffs)/(2*np.pi*KT)) + 1
    t_n = np.array(pool.map(task_matsubara, [(n, iw, mu, qs) for iw in range(nmax)]))      # (n, nq)
    k = mesh_k(n, 0, n**3)
    s = MODEL.s_of_k(k)
    big = np.array([1.0/MODEL.gam_a[0], 1.0/MODEL.gam_a[1]])
    t_inf = mats_t(n, s, big[0], big[1], qs)
    # T(z_bar) = conj(T(z)) after the k sum (q -> -q symmetric, bcc): checked on n = 0 below
    omega_n = (2*np.arange(nmax) + 1)*np.pi*KT
    z_n = mu + 1j*omega_n
    dd_all = (t_n[:, 1:] - t_n[:, :1]) - (t_inf[None, 1:] - t_inf[None, :1])          # complex, per q
    res, res_direct = {}, {}
    for wc in cutoffs:
        m = int(wc/(2*np.pi*KT))
        direct = np.sum(np.real(dd_all[:m]), axis=0)
        lo = m//3                                      # fit window [wc/3, wc]: dd(z) = sum_p a_p z^-p, a_p real
        zw = z_n[lo:m]
        tail, a1 = [], []
        for iq in range(dd_all.shape[1]):
            cols = np.stack([zw**(-p) for p in range(1, 6)], axis=1)
            a = np.linalg.lstsq(np.vstack([cols.real, cols.imag]), np.concatenate([dd_all[lo:m, iq].real, dd_all[lo:m, iq].imag]), rcond=None)[0]
            a1.append(a[0])
            acc = 0.0
            for start in range(m, 4_000_000, 500_000):
                zz = mu + 1j*(2*np.arange(start, start + 500_000) + 1)*np.pi*KT
                acc += np.sum(np.real(sum(a[p - 1]*zz**(-p) for p in range(1, 6))))
            tail.append(acc)
        # kT sum_n e^{i w_n eta} a1/z_n = a1 f(-mu), whereas the symmetric real-part sum gives a1 (f(-mu) - 1/2):
        # the convergence factor of Omega = -kT sum_n e^{i w_n eta} ln det(z_n - H) adds a1/2 (terms 1/z^p, p >= 2, converge).
        res_direct[wc] = 0.25*KT*2.0*direct
        res[wc] = 0.25*(KT*2.0*(direct + np.array(tail)) + 0.5*np.array(a1))
    z = mu - 1j*omega_n[0]
    pu, pd = MODEL.p_alpha(z, 0), MODEL.p_alpha(z, 1)
    sym = np.abs(mats_t(n, s, pu, pd, qs) - np.conj(t_n[0])).max()
    return res, res_direct, sym


# ------------------------------------------------------------------ checks


def sorted_eigs(h):
    return np.sort(np.linalg.eigvalsh(h), axis=-1)


def check_p0(model, export, ks_idx):
    meta = (Path(export)/"nc_first"/"fermi_surface.meta").read_text()
    nk = int(re.search(r"n_kpoints = (\d+)", meta).group(1))
    nb = int(re.search(r"n_bands = (\d+)", meta).group(1))
    kk = np.fromfile(Path(export)/"nc_first"/"fermi_surface.bin", dtype="<f8")[:3*nk].reshape(3, nk, order="F").T
    res = {}
    for order, fn in (("first", model.h_first), ("second", model.h_second)):
        sp = np.fromfile(Path(export)/f"nc_{order}"/"fermi_surface.spin.bin", dtype="<f8", count=2*(nb//2)*nk).reshape(2, nb//2, nk, order="F")
        dev = 0.0
        for ik in ks_idx:
            for sg in (0, 1):
                dev = max(dev, np.abs(sorted_eigs(fn(model.s_of_k(kk[ik:ik + 1]), sg))[0] - np.sort(sp[sg, :, ik])).max())
        res[order] = dev
    return res, kk


def check_c1(model, ks, thetas_q):
    worst_col, worst_rot = 0.0, 0.0
    for sg in (0, 1):
        s = model.s_of_k(ks)
        ev = np.linalg.eigvalsh(model.h_gamma(s, sg))
        for ik in range(len(ks)):
            for e in ev[ik]:
                m = np.diag(model.p_alpha(e, sg)) - s[ik]
                sv = np.linalg.svd(m, compute_uv=False)
                worst_col = max(worst_col, sv[-1]/sv[0])
    for q, th in thetas_q:
        h, sr = model.h_rot(ks, q, th)
        ev = np.linalg.eigvalsh(h)
        for ik in range(len(ks)):
            for e in ev[ik]:
                p = np.concatenate([model.p_alpha(e, 0), model.p_alpha(e, 1)])
                sv = np.linalg.svd(np.diag(p) - sr[ik], compute_uv=False)
                worst_rot = max(worst_rot, sv[-1]/sv[0])
    return worst_col, worst_rot


def check_c2(model, ks, q):
    h, _ = model.h_rot(ks, q, 0.0)
    e_rot = np.linalg.eigvalsh(h)
    e_ref = np.sort(np.concatenate([np.linalg.eigvalsh(model.h_gamma(model.s_of_k(ks - 0.5*np.asarray(q)), 0)),
                                    np.linalg.eigvalsh(model.h_gamma(model.s_of_k(ks + 0.5*np.asarray(q)), 1))], axis=1), axis=1)
    return np.abs(e_rot - e_ref).max()


# ------------------------------------------------------------------ O3b: l-resolved moments of H^gamma (collinear)


def hgamma_eig(model, n, rot=None):
    """Eigenpairs of H^gamma_sigma on the Gamma-centred n^3 mesh, both spins; rot = 9x9 unitary V gives H -> V^dag H V."""
    s = model.s_of_k(mesh_k(n, 0, n**3))
    out = []
    for sg in (0, 1):
        h = model.h_gamma(s, sg)
        if rot is not None:
            h = rot.conj().T[None] @ h @ rot[None]
        out.append(np.linalg.eigh(h))
    return out


def solve_mu_eigs(eig, n):
    e = np.concatenate([x[0] for x in eig], axis=1)
    lo_, hi_ = e.min() - 1.0, e.max() + 1.0
    for _ in range(80):
        mu = 0.5*(lo_ + hi_)
        if np.sum(0.5*(1.0 - np.tanh(0.5*(e - mu)/KT)))/n**3 < N_ELECTRONS:
            lo_ = mu
        else:
            hi_ = mu
    return mu


def moments(model, eig, n, mu, shift=0.0):
    """n_q(l, sigma) = sum_k w_k sum_j f(e_j - mu) e_j^q sum_{m in l} |u_{m j}|^2, q = 0, 1, 2, w_k = 1/n^3; shape (3 q, 3 l, 2 sigma).
    shift adds a constant to all eigenvalues and to mu."""
    sel = np.zeros((3, 9))
    sel[model.lmap, np.arange(9)] = 1.0
    out = np.zeros((3, 3, 2))
    for sg in (0, 1):
        e, u = eig[sg]
        e = e + shift
        f = 0.5*(1.0 - np.tanh(0.5*(e - (mu + shift))/KT))
        wl = np.einsum("lm,kmj->klj", sel, np.abs(u)**2)
        for q in range(3):
            out[q, :, sg] = np.einsum("kj,klj->l", f*e**q, wl)/n**3
    return out


def central(nq):
    centre = nq[1]/nq[0]
    return centre, nq[2] - 2*centre*nq[1] + centre**2*nq[0]


def random_l_unitary(model, rng):
    v = np.zeros((9, 9), complex)
    for l in range(3):
        idx = np.where(model.lmap == l)[0]
        g = rng.normal(size=(len(idx), len(idx))) + 1j*rng.normal(size=(len(idx), len(idx)))
        v[np.ix_(idx, idx)] = np.linalg.qr(g)[0]
    return v


# ------------------------------------------------------------------ output


def fmt(v):
    return ", ".join(f"{x:.16e}" for x in np.ravel(v))


def fmtf(v):
    """Fortran (column-major) flattening, first index fastest, so that a Fortran array declared with the stated shape reads correctly."""
    return fmt(np.ravel(v, order="F"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("export")
    ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--skip-values", action="store_true")
    ap.add_argument("--skip-checks", action="store_true")
    ap.add_argument("--meshes", type=int, nargs="*", default=[32, 48, 64, 96])
    ap.add_argument("--check-meshes", type=int, nargs="*", default=[16, 32])
    a = ap.parse_args()
    t_start = time.time()
    model = Model(a.export)
    global MODEL
    MODEL = model
    print(f"orbital l (from S parity): {model.lmap.tolist()}; onsite S diag: {np.round(model.onsite_s, 6).tolist()}")
    print(f"stored center_band/width_band/obar vs recomputed (strux alpha): max |dev| = {model.stored_dev:.2e}")
    rng = np.random.default_rng(1)
    kh = rng.random((7, 3))
    herm = np.abs(model.s_of_k(kh) - np.conj(np.transpose(model.s_of_k(kh), (0, 2, 1)))).max()
    print(f"S^alpha(k) Hermiticity: max |S - S^dag| = {herm:.2e}")
    lines = [f"! O3 draft (developer review). Generated by make_hgamma_reference.py; T = {TEMPERATURE_K} K, kT = {KT:.16e} Ry",
             "&o3_inputs",
             f"  orbital_l = {', '.join(str(int(x)) for x in model.lmap)}",
             f"  alpha = {fmt(ALPHA)}",
             f"  stored_param_dev = {model.stored_dev:.6e}",
             f"  s_hermiticity = {herm:.6e}", "/"]
    new = []                                                  # groups added by the 2026-10 amendment (a), (b), (c); appended at the end
    new += ["&o3_kernel",
            "  ! (a) plain-array kernel inputs, one site, as the generator uses them. Arrays are written in Fortran order (first index fastest).",
            "  ! bond_s_re/bond_s_im(i,j,m): S^alpha block of bond m, raw_orb(i,j) = sbar(j,i,m) as written in gbt_bonds.out; S(k) = sum_m exp(+2 pi i k_direct . bond_n_direct(:,m)) S(:,:,m).",
            "  ! bond_n_direct(3,m): bond vector in direct (integer) coordinates of the bccFe_prim lattice (lattice.nml); bond 1 is onsite.",
            "  ! cc(i,sigma) = c + vmad, dele(i,sigma) = sqrt(Delta) = srdel*wow^(0.5-(l+1)), gamma_minus_alpha(i,sigma) = qpar*wow^(1-2(l+1)) - alpha(l); l = orbital_l(i), sigma = 1 majority.",
            "  ! wow = wav*ang2au/ws_r with the 1-based index l+1 of symbolic_atom.f90:297-298; reproduces the stored center_band/width_band (stored_param_dev, check P0).",
            f"  n_bond = {len(model.n_direct)}",
            f"  bond_n_direct = {fmtf(model.n_direct.T)}",
            f"  bond_s_re = {fmtf(np.transpose(model.s_bonds.real, (1, 2, 0)))}",
            f"  bond_s_im = {fmtf(np.transpose(model.s_bonds.imag, (1, 2, 0)))}",
            f"  cc = {fmtf(model.cc.T)}", f"  dele = {fmtf(model.dele.T)}", f"  gamma_minus_alpha = {fmtf(model.gam_a.T)}", "/"]
    status = {}
    ks_idx, ks = [0, 7, 19, 36, 53], None
    prior = Path(a.out).read_text() if Path(a.out).exists() else ""          # previous output, for the 5-k comparison (M3)
    ctx = mp.get_context("fork")
    pool = ctx.Pool(a.workers, initializer=init_worker, initargs=(a.export,))
    if not a.skip_checks:
        p0, kk = check_p0(model, a.export, ks_idx)
        print(f"P0 max |eig(H) - code periodic_nc| at 5 k: first {p0['first']:.2e}, second {p0['second']:.2e} Ry (target 1e-10)")
        if max(p0.values()) > 1e-10:
            print("P0 NOT MET: stop"); sys.exit(2)
        ks = kk[ks_idx]
        c1 = check_c1(model, ks, [((0.25, 0.0, 0.5), np.radians(30.0)), ((0.125, -0.125, 0.125), np.radians(5.0))])
        print(f"C1 max smallest/largest singular value of P^alpha(eps) - S^alpha: collinear {c1[0]:.2e}, rotated {c1[1]:.2e}")
        c2 = max(check_c2(model, ks, q) for q in Q_LIST)
        print(f"C2 theta = 0: max |eig(H_rot) - eig(H^gamma_up(k-q/2)) U eig(H^gamma_dn(k+q/2))| = {c2:.2e} Ry")
        lines += ["&o3_checks", f"  p0_first_ry = {p0['first']:.6e}", f"  p0_second_ry = {p0['second']:.6e}",
                  f"  c1_collinear = {c1[0]:.6e}", f"  c1_rotated = {c1[1]:.6e}", f"  c2_ry = {c2:.6e}", "/"]
        mu_ref = -0.069305
        lines.append("&o3_kpoints")
        lines.append(f"  mu_context_ry = {mu_ref}")
        for ik, k in zip(ks_idx, ks):
            for sg in (0, 1):
                s = model.s_of_k(k[None])
                e1, e2, eg = (sorted_eigs(f(s, sg))[0] for f in (model.h_first, model.h_second, model.h_gamma))
                near = np.abs(eg - mu_ref) < 0.1
                lines += [f"  k{ik}_direct = {fmt(k)}", f"  k{ik}_spin{sg + 1}_h_first = {fmt(e1)}", f"  k{ik}_spin{sg + 1}_h_second = {fmt(e2)}",
                          f"  k{ik}_spin{sg + 1}_h_gamma = {fmt(eg)}", f"  k{ik}_spin{sg + 1}_near_mu = {' '.join('T' if x else 'F' for x in near)}",
                          f"  k{ik}_spin{sg + 1}_dev_first_minus_gamma_near_mu = {fmt((e1 - eg)[near])}",
                          f"  k{ik}_spin{sg + 1}_dev_second_minus_gamma_near_mu = {fmt((e2 - eg)[near])}"]
                print(f"k{ik} spin{sg + 1}: within 0.1 Ry of mu: first-order dev {np.round((e1 - eg)[near], 5).tolist()}, second-order dev {np.round((e2 - eg)[near], 5).tolist()}")
        lines.append("/")
        new += ["&o3_rotated", "  ! (b) sorted eigenvalues [Ry] of H_rot (18 values) at the k-points of &o3_kpoints; q1 = (0.125,-0.125,0.125), q2 = (0.25,-0.25,0.25) direct, phi = 0",
                f"  k_index = {', '.join(str(i) for i in ks_idx)}", f"  q1_direct = {fmt(Q_LIST[0])}", f"  q2_direct = {fmt(Q_LIST[1])}"]
        for ik, k in zip(ks_idx, ks):
            for iq, q in enumerate(Q_LIST):
                for th in ROT_THETA_DEG:
                    e = np.sort(np.linalg.eigvalsh(model.h_rot(k[None], q, np.radians(th))[0])[0])
                    new.append(f"  k{ik}_q{iq + 1}_theta{th:g}_h_rot = {fmt(e)}")
        new.append("/")
        status["P0"] = max(p0.values()) <= 1e-10
        status["C1"] = max(c1) < C1_C2_TOL
        status["C2"] = c2 < C1_C2_TOL
        cutoffs = [10.0, 20.0, 40.0]
        for n in a.check_meshes:
            t0 = time.time()
            mu = solve_mu(pool, n)
            jr, _ = richardson_j(pool, n, mu)
            res, res_direct, sym = lkag_matsubara(pool, n, mu, cutoffs)
            print(f"C3 mesh {n}^3 mu = {mu:.12f} Ry  ({time.time() - t0:.0f} s); T(z_bar) = conj T(z): {sym:.1e}")
            lines += [f"&o3_c3_mesh{n}", f"  mu_ry = {mu:.16e}", f"  j0_minus_jq_hrot_ry = {fmt(jr)}"]
            for wc in cutoffs:
                print(f"   cutoff {wc:5.1f} Ry: symmetric direct sum {np.round(res_direct[wc]*1e3, 6).tolist()}  tail + convergence-factor term {np.round(res[wc]*1e3, 6).tolist()} mRy")
                lines += [f"  j0_minus_jq_lkag_direct_cutoff_{wc:g}_ry = {fmt(res_direct[wc])}", f"  j0_minus_jq_lkag_tail_cutoff_{wc:g}_ry = {fmt(res[wc])}"]
            ext = res[max(cutoffs)]
            lines += [f"  j0_minus_jq_lkag_ry = {fmt(ext)}"]
            rel = (ext - np.array(jr))/np.array(jr)
            print(f"   LKAG (cutoff {max(cutoffs):g} Ry with fitted tail) {np.round(ext*1e3, 6).tolist()} mRy; H_rot Richardson {np.round(np.array(jr)*1e3, 6).tolist()} mRy; relative difference {rel.tolist()}")
            lines += [f"  relative_difference = {fmt(rel)}", "/"]
            status["C3"] = status.get("C3", True) and bool(np.max(np.abs(rel)) < C3_REL_TOL)
            if n == 32:
                om00, n00 = omega_probe(pool, n, (0, 0, 0), 0.0, mu)
                new += [f"&o3_driver_mesh{n}", "  ! (c) direct sums on the Gamma-centred mesh, 300 K, mu solved for (q,theta) = (0,0); Omega in Ry per cell, N electrons per cell; theta = 5 deg, no Richardson",
                        f"  mu_ry = {mu:.16e}", f"  omega_00_ry = {om00:.16e}", f"  n_00 = {n00:.16e}"]
                for iq, q in enumerate(Q_LIST):
                    o5, n5 = omega_probe(pool, n, q, np.radians(5.0), mu)
                    o0, _ = omega_probe(pool, n, q, 0.0, mu)
                    new += [f"  omega_q{iq + 1}_theta5_ry = {o5:.16e}", f"  n_q{iq + 1}_theta5 = {n5:.16e}", f"  omega_q{iq + 1}_theta0_ry = {o0:.16e}"]
                new.append("/")
    if not a.skip_values:
        for n in a.meshes:
            t0 = time.time()
            mu = solve_mu(pool, n)
            jr, om0 = richardson_j(pool, n, mu)
            print(f"values mesh {n}^3: mu = {mu:.12f} Ry, J(0)-J(q) = {np.round(np.array(jr)*1e3, 6).tolist()} mRy ({time.time() - t0:.0f} s)")
            lines += [f"&o3_values_mesh{n}", f"  mu_ry = {mu:.16e}", f"  j0_minus_jq_ry = {fmt(jr)}",
                      f"  omega_q0_theta_ry = {fmt(om0)}", "/"]
    # ---- O3b: l-resolved moments of H^gamma (collinear, per spin), Gamma-centred full meshes, 300 K, mu for N_ELECTRONS on each mesh
    rng_m = np.random.default_rng(2)
    e_shift = 0.37
    dev_n0, dev_rot, dev_shift, dev_eig = 0.0, 0.0, 0.0, None
    for n in MOMENT_MESHES:
        eig = hgamma_eig(model, n)
        mu = solve_mu_eigs(eig, n)
        nq = moments(model, eig, n, mu)
        centre, cq2 = central(nq)
        dev_n0 = max(dev_n0, abs(nq[0].sum() - N_ELECTRONS))
        v = random_l_unitary(model, rng_m)
        dev_rot = max(dev_rot, np.abs(moments(model, hgamma_eig(model, n, v), n, mu) - nq).max())
        nqs = moments(model, eig, n, mu, e_shift)
        cs, cq2s = central(nqs)
        dev_shift = max(dev_shift, np.abs(nqs[0] - nq[0]).max(), np.abs(cs - (centre + e_shift)).max(), np.abs(cq2s - cq2).max())
        print(f"moments mesh {n}^3: mu = {mu:.12f} Ry, sum n0 = {nq[0].sum():.12f}")
        for l in range(3):
            for sg in (0, 1):
                print(f"   l={l} spin{sg + 1}: n0 {nq[0, l, sg]:.9f}  n1 {nq[1, l, sg]:.9f}  n2 {nq[2, l, sg]:.9f}  centre {centre[l, sg]:.9f}  central q2 {cq2[l, sg]:.9f}")
        if n == 4 and ks is not None and prior:
            dev_eig = 0.0
            for ik, k in zip(ks_idx, ks):
                m = int(np.argmin(np.abs(((mesh_k(4, 0, 64) - k + 0.5) % 1.0 - 0.5)).sum(axis=1)))
                for sg in (0, 1):
                    old = re.search(rf"k{ik}_spin{sg + 1}_h_gamma = (.*)", prior)
                    if old:
                        dev_eig = max(dev_eig, np.abs(eig[sg][0][m] - np.array(read_numbers(old.group(1)))).max())
                    else:
                        dev_eig = None
                        break
                if dev_eig is None:
                    break
        new += [f"&o3b_moments_mesh{n}",
                "  ! (O3b) l-resolved moments of H^gamma per spin, Gamma-centred full mesh n^3, 300 K, mu solved for 8 electrons (both spins), eigenvalues include vmad. Arrays (l = 0, 1, 2; sigma = 1 majority, 2 minority), Fortran order.",
                "  ! n0 [electrons], n1 = sum f e [Ry electrons], n2 = sum f e^2 [Ry^2 electrons], centre = n1/n0 [Ry], central_q2 = n2 - 2 centre n1 + centre^2 n0 [Ry^2 electrons]; mu_ry [Ry]",
                f"  mu_ry = {mu:.16e}", f"  n0 = {fmtf(nq[0])}", f"  n1 = {fmtf(nq[1])}", f"  n2 = {fmtf(nq[2])}",
                f"  centre_ry = {fmtf(centre)}", f"  central_q2_ry2 = {fmtf(cq2)}", "/"]
    print(f"M1 |sum n0 - 8| = {dev_n0:.2e}; M2 per-l moments under random unitary within each l shell: max |change| = {dev_rot:.2e}; "
          f"M3 eigenvalues vs previous 5-k values: {'n/a' if dev_eig is None else f'{dev_eig:.2e} Ry'}; M4 energy zero shifted by {e_shift} Ry: max |change| = {dev_shift:.2e}")
    status["M1"] = dev_n0 < M_TOL
    status["M2"] = dev_rot < M_TOL
    status["M3"] = dev_eig is not None and dev_eig < M_TOL
    status["M4"] = dev_shift < M_TOL
    Path(a.out).write_text("\n".join(lines + new) + "\n")
    for name, ok in status.items():
        print(f"self-check {name}: {'ok' if ok else 'DEVIATION (see values above)'}")
    print(f"total wall time {time.time() - t_start:.0f} s; written {a.out}")


if __name__ == "__main__":
    main()
