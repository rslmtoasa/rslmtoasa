#!/usr/bin/env python3
"""Stage A toy study: linear prediction (LP) of Chebyshev moments (numpy only).

Run: python3 tests/validation/val_cheb_lp_toy.py
All references (closed forms, quadratures, k-mesh sums, Lanczos) are computed
in this script. No production code is involved. Output is tables only.
Units: scaled energy x = (E-b)/a, band half-width W = 0.85 (code scaling).
"""
import numpy as np

LLD = 150
N = 2 * LLD + 2            # known moments, as in the triad Fe deck
F = 4
NEXT = F * N               # extended length
W = 0.85
RULES = [("drop", 0.0), ("unit", 0.0), ("shrink", 1e-3), ("shrink", 1e-2)]
DELTAS = (1e-12, 1e-8, 1e-6)
PFS = (8, 4, 2)            # p = n_known // pf
NFS = (4, 2)               # n_lo = n_known // nf


def rname(r):
    return r[0] if r[0] != "shrink" else "sh%.0e" % r[1]


# ----------------------------------------------------------------- LP core
def lp_extend(x, n_ext, p, n_lo, delta, rule):
    """x: (C, n_known) real channels sharing one real coefficient vector.
    Returns ext (C, n_ext), info dict. Training window n in [n_lo, n_known-1]."""
    x = np.atleast_2d(np.asarray(x, dtype=float))
    c, n = x.shape
    assert n_lo >= p and n_ext > n
    a_mat = np.stack([x[:, n_lo - i:n - i] for i in range(1, p + 1)], axis=-1).reshape(-1, p)
    y = x[:, n_lo:n].reshape(-1)
    rr = a_mat.T @ a_mat
    r = a_mat.T @ y
    s, u = np.linalg.eigh(rr)
    keep = s > delta * s.max()
    a = -(u[:, keep] @ ((u[:, keep].T @ r) / s[keep]))
    res = np.linalg.norm(a_mat @ a + y) / np.linalg.norm(y)
    m = np.zeros((p, p))
    m[0, :] = -a
    m[1:, :-1] = np.eye(p - 1)
    lam, v = np.linalg.eig(m)
    mod = np.abs(lam)
    bad = mod >= 1.0
    lam2 = lam.copy()
    kind, eps = rule
    if kind == "drop":
        lam2[bad] = 0.0
    elif kind == "unit":
        lam2[bad] = lam[bad] / mod[bad]
    else:
        lam2[bad] = lam[bad] / mod[bad] * (1.0 - eps)
    s0 = x[:, n - 1::-1][:, :p]                       # x(n-1) ... x(n-p)
    wc = np.linalg.solve(v, s0.T.astype(complex))     # (p, C)
    k = np.arange(1, n_ext - n + 1)
    out = ((lam2[None, :] ** k[:, None]) * v[0, :][None, :]) @ wc
    ext = np.concatenate([x, out.real.T], axis=1)
    info = dict(rank=int(keep.sum()), nunst=int(bad.sum()), lam_max=float(mod.max()),
                fit_res=float(res), imag=float(np.abs(out.imag).max()))
    return ext, info


def lp_try(x, n_ext, p, n_lo, delta, rule):
    try:
        return lp_extend(x, n_ext, p, n_lo, delta, rule)
    except np.linalg.LinAlgError:
        return None, None


def configs(n_known):
    for pf in PFS:
        for nf in NFS:
            p, nlo = n_known // pf, n_known // nf
            if nlo < p:
                continue
            for d in DELTAS:
                for rule in RULES:
                    yield (pf, nf, d, rule), p, nlo


# --------------------------------------------------- reconstruction / metrics
def jackson(n):
    k = np.arange(n)
    th = np.pi * k / (n + 1)
    return ((n - k + 1) * np.cos(th) + np.sin(th) / np.tan(np.pi / (n + 1))) / (n + 1)


XG = np.linspace(-0.98, 0.98, 3920)       # even count: x = 0 is not on the grid
COS = np.cos(np.outer(np.arange(NEXT), np.arccos(XG)))
INT = np.abs(XG) < 0.75
EDGE = (np.abs(XG) >= 0.75) & (np.abs(XG) <= 0.95)


def dos(mu):
    n = len(mu)
    c = jackson(n) * mu
    c[1:] *= 2.0
    return (c @ COS[:n]) / (np.pi * np.sqrt(1.0 - XG ** 2))


def errs(rho, ref):
    d = rho - ref
    return (np.sqrt(np.mean(d[INT] ** 2)), np.abs(d[INT]).max(),
            np.sqrt(np.mean(d[EDGE] ** 2)), np.abs(d[EDGE]).max())


def fmt_e(e):
    return " ".join("%9.2e" % v for v in e)


# ------------------------------------------------------------------ fixtures
def semicircle_exact(x):
    return 2.0 / (np.pi * W ** 2) * np.sqrt(np.maximum(W ** 2 - x ** 2, 0.0))


def semicircle_moments(nmom, m=6000):
    """Gauss-Chebyshev (2nd kind) quadrature, exact to degree 2m-1 in y=x/W."""
    j = np.arange(1, m + 1)
    y = np.cos(j * np.pi / (m + 1))
    wt = 2.0 / (m + 1) * np.sin(j * np.pi / (m + 1)) ** 2
    x = W * y
    mu = np.empty(nmom)
    t0, t1 = np.ones_like(x), x.copy()
    mu[0], mu[1] = wt.sum(), wt @ t1
    for n in range(2, nmom):
        t0, t1 = t1, 2 * x * t1 - t0
        mu[n] = wt @ t1
    return mu


def ellipk(m):
    """Complete elliptic integral K(m) (parameter m) by AGM."""
    a, b = np.ones_like(m), np.sqrt(1.0 - m)
    for _ in range(40):
        a, b = 0.5 * (a + b), np.sqrt(a * b)
    return np.pi / (2.0 * a)


def square_exact(x):
    c = W / 2.0
    out = np.zeros_like(x)
    ins = np.abs(x) < W
    out[ins] = ellipk(1.0 - (x[ins] / W) ** 2) / (np.pi ** 2 * c)
    return out


def square_moments(nmom, m=1300):
    """Infinite 2D square lattice, eps = c (cos kx + cos ky), c = W/2.
    Uniform m x m mesh: trapezoid exact for trig polynomials of degree < m."""
    assert m > nmom and m % 2 == 0
    c = W / 2.0
    j = np.arange(m // 2 + 1)
    wt1 = np.where((j == 0) | (j == m // 2), 1.0, 2.0) / m
    ck = np.cos(2 * np.pi * j / m)
    eps = c * (ck[:, None] + ck[None, :])
    wt = wt1[:, None] * wt1[None, :]
    mu = np.empty(nmom)
    t0, t1 = np.ones_like(eps), eps.copy()
    mu[0], mu[1] = wt.sum(), (wt * t1).sum()
    for n in range(2, nmom):
        t0, t1 = t1, 2 * eps * t1 - t0
        mu[n] = (wt * t1).sum()
    return mu


# ------------------------------------------------------------------- A1
def a1_signal(n):
    rho = np.array([0.975, 0.982, 0.988, 0.991, 0.994, 0.997])
    th = np.array([0.35, 0.80, 1.25, 1.70, 2.20, 2.75])
    w = np.array([1.0, 0.7, 1.3, 0.5, 0.9, 0.6])
    k = np.arange(n)[:, None]
    return (w * rho ** k * np.cos(k * th)).sum(axis=1)


def run_a1(store):
    print("\n=== A1 damped poles (K=6, rho<1), N=%d -> %d ===" % (N, NEXT))
    print("err = max_{n>=N}|x_LP - x_exact| / max_n|x_exact| (n<N); noise = relative to rms of x(n<N)")
    ex = a1_signal(NEXT)
    scale = np.abs(ex[:N]).max()
    rng = np.random.default_rng(12345)
    unit_noise = rng.standard_normal(N)
    for noise in (0.0, 1e-12, 1e-7):
        xk = ex[:N] + noise * np.sqrt(np.mean(ex[:N] ** 2)) * unit_noise
        print("\nnoise=%g" % noise)
        print("pf nf   delta  | rank nunst   lam_max   fit_res |" + "".join("%11s" % rname(r) for r in RULES))
        seen = {}
        for key, p, nlo in configs(N):
            pf, nf, d, rule = key
            ext, info = lp_try(xk, NEXT, p, nlo, d, rule)
            e = np.inf if ext is None else np.abs(ext[0, N:] - ex[N:]).max() / scale
            store.setdefault(key, {})[noise] = e
            seen.setdefault((pf, nf, d), []).append((e, info))
        for (pf, nf, d), lst in seen.items():
            info = lst[0][1]
            head = "%2d %2d %7.0e | " % (pf, nf, d)
            head += ("%4d %5d %9.5f %9.2e |" % (info["rank"], info["nunst"], info["lam_max"], info["fit_res"])
                     if info else "   fail |")
            print(head + "".join("%11.2e" % e for e, _ in lst))


# ------------------------------------------------------------ A2 / A3 / A4
def dos_table(title, mu_known, mu_ref_full, rho_ref, n_known, store, tag, noise=0.0, extra=None):
    """Baselines + LP grid for one fixture. mu_known: moments (n_known,);
    mu_ref_full: reference moments (NEXT,) for Jackson@NEXT and moment errors."""
    print("\n--- %s ---" % title)
    print("method                         p  nlo  delta  rule  nunst  lam_max |  int_rms  int_max  edge_rms edge_max" +
          ("  mom_err" if mu_ref_full is not None else ""))
    base = errs(dos(mu_known), rho_ref)
    print("%-28s %s" % ("Jackson@%d (given moments)" % n_known, " " * 30 + fmt_e(base)))
    if mu_ref_full is not None:
        print("%-28s %s" % ("Jackson@%d (reference)" % NEXT, " " * 30 + fmt_e(errs(dos(mu_ref_full), rho_ref))))
    for key, p, nlo in configs(n_known):
        pf, nf, d, rule = key
        ext, info = lp_try(mu_known, NEXT, p, nlo, d, rule)
        if ext is None:
            print("LP %3d %4d %7.0e %-6s  fail" % (p, nlo, d, rname(rule)))
            store.setdefault(key, {})[tag] = np.inf
            continue
        e = errs(dos(ext[0]), rho_ref)
        store.setdefault(key, {})[tag] = max(e[0] / base[0], e[2] / base[2])
        mom = ""
        if mu_ref_full is not None:
            mom = " %9.2e" % np.abs(ext[0, n_known:] - mu_ref_full[n_known:NEXT]).max()
        print("LP  %-24s %3d %4d %7.0e %-6s %5d %8.5f | %s%s" %
              ("", p, nlo, d, rname(rule), info["nunst"], info["lam_max"], fmt_e(e), mom))


def run_a2(store, mu_sc):
    print("\n=== A2 semicircle, half-width %.2f ===" % W)
    print("sanity: mu_2 = %.12f (closed form W^2/2-1 = %.12f)" % (mu_sc[2], W ** 2 / 2 - 1))
    rho_ref = semicircle_exact(XG)
    rng = np.random.default_rng(2)
    unit_noise = rng.standard_normal(N)
    for noise in (0.0, 1e-7):
        xk = mu_sc[:N] + noise * np.sqrt(np.mean(mu_sc[:N] ** 2)) * unit_noise
        dos_table("A2 semicircle noise=%g (ratios in store use noise 0 only)" % noise, xk, mu_sc, rho_ref, N,
                  store, "A2" if noise == 0.0 else "A2n")


def run_a3(store, mu_sq):
    print("\n=== A3 2D square lattice nn (infinite), band |x|<=%.2f ===" % W)
    print("sanity: mu_2 = %.12f (closed form 2c^2-1 = %.12f); odd moments max|mu_odd| = %.1e" %
          (mu_sq[2], 2 * (W / 2) ** 2 - 1, np.abs(mu_sq[1::2]).max()))
    rho_ref = square_exact(XG)
    print("sanity: trapezoid int rho_exact dx over grid = %.4f (log singularity at 0, grid cuts at 0.98)" %
          (rho_ref.sum() * (XG[1] - XG[0])))
    rng = np.random.default_rng(3)
    unit_noise = rng.standard_normal(N)
    for noise in (0.0, 1e-7):
        xk = mu_sq[:N] + noise * np.sqrt(np.mean(mu_sq[:N] ** 2)) * unit_noise
        dos_table("A3 square noise=%g" % noise, xk, mu_sq, rho_ref, N, store, "A3" if noise == 0.0 else "A3n")


def cluster_moments(radius, napps):
    """Diamond cluster |dx|+|dy|<=radius, H = (W/4) sum_nn; doubling trick."""
    g = 2 * radius + 3
    mask = np.zeros((g, g))
    ctr = radius + 1
    ix = np.arange(g) - ctr
    mask[(np.abs(ix)[:, None] + np.abs(ix)[None, :]) <= radius] = 1.0
    h = W / 4.0

    def ap(v):
        return h * mask * (np.roll(v, 1, 0) + np.roll(v, -1, 0) + np.roll(v, 1, 1) + np.roll(v, -1, 1))
    phi = [np.zeros((g, g)) for _ in range(napps + 1)]
    phi[0][ctr, ctr] = 1.0
    phi[1] = ap(phi[0])
    for n in range(1, napps):
        phi[n + 1] = 2 * ap(phi[n]) - phi[n - 1]
    mu = np.zeros(2 * napps)
    mu0 = np.sum(phi[0] * phi[0])
    mu1 = np.sum(phi[1] * phi[0])
    for n in range(napps + 1):
        if 2 * n < len(mu):
            mu[2 * n] = 2 * np.sum(phi[n] * phi[n]) - mu0
        if 2 * n + 1 < len(mu) and n + 1 <= napps:
            mu[2 * n + 1] = 2 * np.sum(phi[n + 1] * phi[n]) - mu1
    return mu, int(mask.sum())


def run_a4(store, mu_sq, napps, radius):
    nb, na = 2 * radius + 2, 2 * napps
    print("\n=== A4 cluster boundary: diamond R=%d hops, %d H-applications ===" % (radius, napps))
    mu_cl, nsites = cluster_moments(radius, napps)
    print("sites=%d; cluster moments n=0..%d. Bulk-exact expected for n<=2R+1=%d" % (nsites, na - 1, nb - 1))
    d = np.abs(mu_cl - mu_sq[:na])
    print("max|mu_cl-mu_bulk| for n<%d: %.2e ; for %d<=n<%d: %.2e (boundary effect)" %
          (nb, d[:nb].max(), nb, na, d[nb:].max()))
    rho_ref = square_exact(XG)
    dos_table("A4 before boundary: train on exact-bulk n<%d, extend to %d" % (nb, NEXT),
              mu_cl[:nb], mu_sq, rho_ref, nb, store, "A4b")
    dos_table("A4 after boundary: train on cluster moments n<%d (include boundary-affected), extend to %d" % (na, NEXT),
              mu_cl, mu_sq, rho_ref, na, store, "A4a")
    return mu_cl


# -------------------------------------------------------------------- A5
def lanczos(napps, radius):
    g = 2 * radius + 3
    mask = np.zeros((g, g))
    ctr = radius + 1
    ix = np.arange(g) - ctr
    mask[(np.abs(ix)[:, None] + np.abs(ix)[None, :]) <= radius] = 1.0
    h = W / 4.0

    def ap(v):
        v = v.reshape(g, g)
        return (h * mask * (np.roll(v, 1, 0) + np.roll(v, -1, 0) + np.roll(v, 1, 1) + np.roll(v, -1, 1))).ravel()
    vs = [np.zeros(g * g)]
    vs[0][ctr * g + ctr] = 1.0
    a, b = [], [0.0]
    for n in range(napps):
        w = ap(vs[n])
        a.append(vs[n] @ w)
        w = w - a[n] * vs[n] - (b[n] * vs[n - 1] if n > 0 else 0.0)
        for _ in range(2):                       # full re-orthogonalisation
            for v in vs:
                w = w - (v @ w) * v
        b.append(np.linalg.norm(w))
        vs.append(w / b[-1])
    return np.array(a), np.array(b[1:])          # a_0..a_{L-1}, b_1..b_L


def lanczos_bp_dos(a, b, ntail=10):
    """Continued fraction with Beer-Pettifor constant terminator, z = x + i*eta."""
    ll = len(a)
    ainf, binf = a[-ntail:].mean(), b[-ntail:].mean()
    z = XG + 1e-9j
    zp = z - ainf
    sq = np.sqrt(zp * zp - 4 * binf ** 2)
    t = (zp - sq) / (2 * binf ** 2)
    t = np.where(t.imag > 0, (zp + sq) / (2 * binf ** 2), t)
    g = t
    for n in range(ll - 1, -1, -1):
        g = 1.0 / (z - a[n] - b[n] ** 2 * g)
    return -g.imag / np.pi, ainf, binf


def jacobi_moments(a, b, nmom):
    ll = len(a)
    j = np.diag(a) + np.diag(b[:ll - 1], 1) + np.diag(b[:ll - 1], -1)
    v0 = np.zeros(ll)
    v0[0] = 1.0
    v1 = j @ v0
    mu = np.empty(nmom)
    mu[0], mu[1] = 1.0, v1[0]
    for n in range(2, nmom):
        v0, v1 = v1, 2 * (j @ v1) - v0
        mu[n] = v1[0]
    return mu


def best_per_rule(store):
    out = {}
    for rule in RULES:
        cand = [(score(v), k) for k, v in store.items() if k[3] == rule]
        out[rname(rule)] = min(cand)[1]
    return out


def run_a5(store, mu_sq, mu_cl, napps, radius):
    print("\n=== A5 matched cost: L=%d H-applications on the R=%d cluster; scored vs infinite lattice ===" % (napps, radius))
    rho_ref = square_exact(XG)
    a, b = lanczos(napps, radius)
    rho_l, ainf, binf = lanczos_bp_dos(a, b)
    print("Lanczos: a_inf=%.5f b_inf=%.5f (mean of last 10); tail spread a: %.1e b: %.1e; band edge a_inf+2b_inf=%.4f" %
          (ainf, binf, np.ptp(a[-10:]), np.ptp(b[-10:]), ainf + 2 * binf))
    print("method                                    |  int_rms  int_max  edge_rms edge_max")
    print("%-42s| %s" % ("Lanczos + BP terminator", fmt_e(errs(rho_l, rho_ref))))
    na = 2 * napps
    print("%-42s| %s" % ("Chebyshev Jackson@%d" % na, fmt_e(errs(dos(mu_cl), rho_ref))))
    print("%-42s| %s" % ("Jackson@%d bulk-exact (limit)" % NEXT, fmt_e(errs(dos(mu_sq), rho_ref))))
    mu_j = jacobi_moments(a, b, na + 4)
    dm = np.abs(mu_j[:na + 4] - np.concatenate([mu_cl, [np.nan] * 4]))
    ok = np.where(dm[:na] > 1e-9)[0]
    print("Lanczos->Chebyshev moments vs doubling moments: max|diff| n<%d = %.2e; first n with diff>1e-9: %s" %
          (na, np.nanmax(dm[:na]), ok[0] if len(ok) else "none"))
    for name, k in best_per_rule(store).items():
        pf, nf, d, rule = k
        p, nlo = na // pf, na // nf
        for lab, mu in (("Cheb+LP", mu_cl), ("Lanczos->Cheb+LP", mu_j[:na])):
            ext, info = lp_try(mu, NEXT, p, nlo, d, rule)
            e = errs(dos(ext[0]), rho_ref) if ext is not None else (np.nan,) * 4
            print("%-42s| %s" % ("%s pf=%d nf=%d d=%.0e %s" % (lab, pf, nf, d, name), fmt_e(e)))


# ---------------------------------------------------------------- ranking
def score(v):
    """Worst ratio (LP error / Jackson@N_known error, RMS, interior and edge) over the
    clean-moment fixtures A2, A3, A4b. A4a (boundary-contaminated training) is reported separately."""
    t = [v.get(k, np.inf) for k in ("A2", "A3", "A4b")]
    return max(t)


def summary(store_dos, store_a1):
    print("\n=== Ranking: score = max over A2,A3,A4b (clean moments) of max(int_rms, edge_rms)_LP / (same)_Jackson@N_known ===")
    print("(score < 1: LP beats Jackson@N on every fixture. A1 columns: max rel extrapolation error at noise 0 / 1e-12 / 1e-7)")
    print("rank  pf nf   delta  rule     score |     A2     A3    A4b    A4a |   A1(0)  A1(1e-12)  A1(1e-7)")
    rows = sorted(store_dos.items(), key=lambda kv: score(kv[1]))
    for i, (k, v) in enumerate(rows[:12]):
        pf, nf, d, rule = k
        a1 = store_a1.get(k, {})
        print("%4d  %2d %2d %7.0e  %-6s %8.3f | %6.3f %6.3f %6.3f %6.3f | %8.1e %9.1e %9.1e" %
              (i + 1, pf, nf, d, rname(rule), score(v), v.get("A2", np.nan), v.get("A3", np.nan),
               v.get("A4b", np.nan), v.get("A4a", np.nan), a1.get(0.0, np.nan), a1.get(1e-12, np.nan), a1.get(1e-7, np.nan)))
    print("\nbest per rule:")
    for name, k in best_per_rule(store_dos).items():
        pf, nf, d, rule = k
        v = store_dos[k]
        a1 = store_a1.get(k, {})
        print("  %-6s pf=%d nf=%d delta=%.0e score=%.3f  A1(0/1e-12/1e-7)=%.1e/%.1e/%.1e" %
              (name, pf, nf, d, score(v), a1.get(0.0, np.nan), a1.get(1e-12, np.nan), a1.get(1e-7, np.nan)))
    k4 = min(store_dos, key=lambda k: store_dos[k].get("A4a", np.inf))
    print("A4a (training window crosses the cluster boundary): best ratio %.3f at pf=%d nf=%d delta=%.0e %s; "
          "median %.3f; max %.3g" % (store_dos[k4]["A4a"], k4[0], k4[1], k4[2], rname(k4[3]),
                                      np.median([v["A4a"] for v in store_dos.values()]),
                                      max(v["A4a"] for v in store_dos.values())))
    nfail = sum(1 for v in store_dos.values() if score(v) == np.inf)
    print("configs with a failed linear solve: %d of %d" % (nfail, len(store_dos)))


def main():
    np.set_printoptions(linewidth=160)
    print("N=%d (lld=%d), F=%d, N_ext=%d, W=%.2f; p=N//pf, n_lo=N//nf; rules: %s" %
          (N, LLD, F, NEXT, W, ", ".join(rname(r) for r in RULES)))
    store_a1, store_dos = {}, {}
    run_a1(store_a1)
    mu_sc = semicircle_moments(NEXT)
    mu_sq = square_moments(NEXT)
    run_a2(store_dos, mu_sc)
    run_a3(store_dos, mu_sq)
    napps = N // 2                    # L: 2L = N moments by doubling
    radius = 75
    mu_cl = run_a4(store_dos, mu_sq, napps, radius)
    summary(store_dos, store_a1)
    run_a5(store_dos, mu_sq, mu_cl, napps, radius)


if __name__ == "__main__":
    main()
