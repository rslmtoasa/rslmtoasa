#!/usr/bin/env python3
"""Model oracles for the k-space GBT kernels (draft for developer review).

O1  One s-like orbital on a simple cubic lattice (lattice constant 1), nearest-neighbour
    hopping t, W_up = W_dn = 1, local-frame onsite C = -Delta/2 (majority) and +Delta/2.
    Writes o1_one_band.nml: closed-form cone eigenvalues for one site per cell, and the
    eigenvalues of a two-site cell (sc doubled along x), parallel and antiparallel.

Conventions (spec, section "Physics"):
  q Cartesian in 2 pi/alat, k and q below are DIRECT coordinates of the cell they are given for;
  the moment at position r advances its azimuth by +q.r; D(a) = diag(e^{-ia/2}, e^{+ia/2});
  site frame U_a = D(phi_a) R_y(theta_a) D(-phi_a); Bloch phases use the full bond vector
  d = R + tau_b - tau_a (atomic gauge).
  R1:  H_{a s,b s'}(k) = C_{a s} delta_ab delta_ss' + W_a s [U_a^dag diag(S_ab(k-q/2), S_ab(k+q/2)) U_b]_{s s'} W_b s'
  One site: eps_pm = eps_bar +- sqrt(delta^2 + Delta^2/4 - delta*Delta*cos(theta)),
            eps_bar = [s(k-q/2) + s(k+q/2)]/2, delta = [s(k-q/2) - s(k+q/2)]/2, s(k) = 2 t sum_i cos(2 pi k_i).
  Antiparallel second site: the frame alone carries the orientation (W and C as for a majority
  site): theta_2 = pi - theta_ss, phi_2 = phi_ss + pi; at theta_ss = 0 this is reference theta = pi.

Self-checks (independent of the GBT algebra; each prints ok or FAIL with the deviation):
  1. closed form vs numeric diagonalisation of the R1 operator (one site, all cases);
  2. for commensurate q (1/2, 1/3, 1/4 along x, in each cell's direct coordinates) the explicit
     lab-frame supercell with spin-independent hopping and onsite -(Delta/2) sigma.m_i,
     m_i = (sin th cos(ph + q.r_i), sin th sin(ph + q.r_i), cos th), at the Bloch vector K that the
     spiral's translation phase fixes (K_x = k_x + M/(2L), M = q_x L, in 2 pi/alat units), has as
     its spectrum the union of the GBT eigenvalues at the L/n_sub folded k-points
     k_x + j/L, j = 0..L/n_sub-1 (n_sub = sites per cell). One site: the closed form is used;
     two sites, parallel and antiparallel: the GBT 4 x 4 eigenvalues.

O2  Force theorem on the O1 one-site model (spec R7, R8, R11): Gamma-centred 24^3 mesh,
    k_B T = 0.002 Ry, mu solved once for (q, theta) = (0, 0) of each set. Writes o2_force_theorem.nml:
    Omega(q,theta) = -k_B T sum_k w ln(1 + e^{-(eps-mu)/k_B T}) at that fixed mu, N(q,theta), the band
    energy at fixed mu without the -mu N term (the mutation set), and the model LKAG exchange as a band
    (Lehmann) sum, with H = -sum_{i/=j} J_ij e_i.e_j:
      J(q) = (Delta^2/4) (1/N_k) sum_k [f(eps_up(k)) - f(eps_dn(k+q))] / [eps_dn(k+q) - eps_up(k)],
      eps_up = s(k) - Delta/2 + shift, eps_dn = s(k) + Delta/2 + shift (limit -f' for a degenerate pair).
    Model form of the J(q) expression accompanying Eq. (7) of Pajda et al., PRB 64, 174402 (2001), with
    P_up - P_dn -> Delta. The paper places spin-up at k+q; the two placements are equal here because
    eps(k) = eps(-k). The finite-T Fermi-function form is our extension; the paper uses a T = 0 contour.
    Sets: 1) N = 0.6 per site; 2) N = 0.9 (a different mu); 3) the set-1 model shifted rigidly by +0.5 Ry
    (both C values), N = 0.6, so mu moves by the same amount.
    Self-checks: Omega(q,0) = Omega(0,0) (R2); Richardson extrapolation (Neville in sin^2 theta, the three
    smallest angles 0.25, 0.5, 1 deg) of [Omega(q,theta) - Omega(0,theta)]/sin^2(theta) equals J(0) - J(q)
    (R8); set 3 equals set 1 in Omega, N and J and its band energy exceeds set 1 by 0.5 N (R7, R11);
    the mutation violates the relation in every set and its violation is printed per set.

  3. the doubled-cell spectrum at K equals the union of the one-site closed forms at the two k that fold
     onto K, with the same Cartesian q (k_cart and k_cart + (1/2,0,0) in 2 pi/alat), all 72 cases.
  Every case also carries Cartesian k and q (2 pi/alat) next to the direct coordinates of its cell.

Usage: python3 make_gbt_references.py [outdir]   (numpy only)
"""
import sys
from itertools import product
from pathlib import Path

import numpy as np

T = -0.05                      # nearest-neighbour hopping (Ry)
DELTA = 0.15                   # exchange splitting (Ry)
KS = [(0.0, 0.0, 0.0), (0.1, 0.2, 0.3), (0.25, -0.1, 0.4)]               # direct
QS = [(0.5, 0.0, 0.0), (1.0 / 3.0, 0.0, 0.0), (0.1, 0.2, 0.05)]          # direct
THETAS = [0.0, 30.0, 90.0, 150.0]                                        # deg
PHIS = [0.0, 40.0]                                                       # deg
SC_ONE = dict(A=np.eye(3), tau=np.zeros((1, 3)))
SC_TWO = dict(A=np.diag([2.0, 1.0, 1.0]), tau=np.array([[0.0, 0, 0], [1.0, 0, 0]]))
TOL_CHECK = 1.0e-12

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def dmat(a):
    return np.diag([np.exp(-0.5j * a), np.exp(0.5j * a)])


def ry(t):
    c, s = np.cos(0.5 * t), np.sin(0.5 * t)
    return np.array([[c, -s], [s, c]], complex)


def frame(theta, phi):
    return dmat(phi) @ ry(theta) @ dmat(-phi)


def cart(k_direct, cell):
    return np.asarray(k_direct, float) @ np.linalg.inv(cell["A"]).T


def bonds(cell):
    """Nearest-neighbour bonds (a, b, R integer, d Cartesian), |d| = 1."""
    out = []
    n = len(cell["tau"])
    for a, b in product(range(n), repeat=2):
        for r in product(range(-3, 4), repeat=3):
            d = np.array(r, float) @ cell["A"] + cell["tau"][b] - cell["tau"][a]
            if abs(np.linalg.norm(d) - 1.0) < 1e-9:
                out.append((a, b, r, d))
    return out


def gbt_h1(cell, k_dir, q_dir, thetas, phis):
    """R1 operator for one k; basis (site, spin) site-major, W = 1, local C = -+Delta/2."""
    n = len(cell["tau"])
    kc, qc = cart(k_dir, cell), cart(q_dir, cell)
    u = [frame(thetas[a], phis[a]) for a in range(n)]
    h = np.zeros((2 * n, 2 * n), complex)
    for a in range(n):
        h[2 * a:2 * a + 2, 2 * a:2 * a + 2] += np.diag([-0.5 * DELTA, 0.5 * DELTA])
    for a, b, _, d in bonds(cell):
        s_minus = T * np.exp(2j * np.pi * np.dot(kc - 0.5 * qc, d))   # first (spin-up) entry: k - q/2
        s_plus = T * np.exp(2j * np.pi * np.dot(kc + 0.5 * qc, d))    # second (spin-down) entry: k + q/2
        h[2 * a:2 * a + 2, 2 * b:2 * b + 2] += u[a].conj().T @ np.diag([s_minus, s_plus]) @ u[b]
    return h


def band(kc):
    return 2.0 * T * np.sum(np.cos(2.0 * np.pi * np.asarray(kc)))


def closed_form(k_dir, q_dir, theta):
    """One site per cell; k, q direct = Cartesian (2 pi/alat)."""
    k, q = np.asarray(k_dir, float), np.asarray(q_dir, float)
    sm, sp = band(k - 0.5 * q), band(k + 0.5 * q)
    bar, dl = 0.5 * (sm + sp), 0.5 * (sm - sp)
    root = np.hypot(dl * np.sin(theta), dl * np.cos(theta) - 0.5 * DELTA)   # = sqrt(dl^2 + Delta^2/4 - dl Delta cos)
    return np.array([bar - root, bar + root])


def two_site_angles(theta, phi, antiparallel):
    if antiparallel:
        return [theta, np.pi - theta], [phi, phi + np.pi]
    return [theta, theta], [phi, phi]


def supercell_spectrum(k_cart, q_cart, length, types, thetas, phis):
    """Explicit lab-frame supercell, `length` sites along x (1 x 1 in y, z); spin-independent
    hopping T, onsite -(Delta/2) sigma.m_i. Returns eigenvalues at the spiral-fixed Bloch vector."""
    kx = k_cart[0] + round(q_cart[0] * length) / (2.0 * length)
    ky_kz = 2.0 * T * (np.cos(2 * np.pi * k_cart[1]) + np.cos(2 * np.pi * k_cart[2]))
    h = np.zeros((2 * length, 2 * length), complex)
    for i in range(length):
        a = types[i]
        ph = phis[a] + 2.0 * np.pi * q_cart[0] * i
        m = (np.sin(thetas[a]) * np.cos(ph), np.sin(thetas[a]) * np.sin(ph), np.cos(thetas[a]))
        h[2 * i:2 * i + 2, 2 * i:2 * i + 2] += ky_kz * np.eye(2) - 0.5 * DELTA * (m[0] * SX + m[1] * SY + m[2] * SZ)
        j = (i + 1) % length
        phase = np.exp(2j * np.pi * kx * length) if i == length - 1 else 1.0   # hop i -> i+1 across the wrap
        h[2 * i:2 * i + 2, 2 * j:2 * j + 2] += T * phase * np.eye(2)
        h[2 * j:2 * j + 2, 2 * i:2 * i + 2] += T * np.conj(phase) * np.eye(2)
    return np.linalg.eigvalsh(h)


def folded_one_site(k_dir, n, theta):
    ev = []
    for j in range(n):
        ev.extend(closed_form(np.array(k_dir) + np.array([j / n, 0, 0]), np.array([1.0 / n, 0, 0]), theta))
    return np.sort(ev)


def folded_two_site(k_dir, n, theta, phi, antiparallel):
    ev = []
    th, ph = two_site_angles(theta, phi, antiparallel)
    for j in range(n):          # k_cart shift j/L with L = 2n, direct shift 2 j/L = j/n
        ev.extend(np.linalg.eigvalsh(gbt_h1(SC_TWO, np.array(k_dir) + np.array([j / n, 0, 0]),
                                            [1.0 / n, 0, 0], th, ph)))
    return np.sort(ev)


def fmt(x):
    return f"{x:.17e}"


def arr(name, values, per_line=4):
    flat = [fmt(v) for v in np.asarray(values, float).ravel(order="F")]
    lines = [", ".join(flat[i:i + per_line]) for i in range(0, len(flat), per_line)]
    return f"  {name} =\n    " + ",\n    ".join(lines) + "\n"


def compute_sets():
    nk, nq, nt, npf = len(KS), len(QS), len(THETAS), len(PHIS)
    one = np.zeros((2, nk, nq, nt, npf))
    par = np.zeros((4, nk, nq, nt, npf))
    anti = np.zeros((4, nk, nq, nt, npf))
    for (ik, k), (iq, q), (it, t), (ip, p) in product(enumerate(KS), enumerate(QS), enumerate(THETAS), enumerate(PHIS)):
        th, ph = np.radians(t), np.radians(p)
        one[:, ik, iq, it, ip] = closed_form(k, q, th)
        par[:, ik, iq, it, ip] = np.linalg.eigvalsh(gbt_h1(SC_TWO, k, q, *two_site_angles(th, ph, False)))
        anti[:, ik, iq, it, ip] = np.linalg.eigvalsh(gbt_h1(SC_TWO, k, q, *two_site_angles(th, ph, True)))
    return one, par, anti


def self_checks(one):
    ok = True

    def report(label, dev):
        nonlocal ok
        good = dev < TOL_CHECK
        ok &= good
        print(f"  {label:<62s} max dev {dev:.3e} Ry  {'ok' if good else 'FAIL'}")

    dev = 0.0
    for (ik, k), (iq, q), (it, t), (ip, p) in product(enumerate(KS), enumerate(QS), enumerate(THETAS), enumerate(PHIS)):
        num = np.linalg.eigvalsh(gbt_h1(SC_ONE, k, q, [np.radians(t)], [np.radians(p)]))
        dev = max(dev, np.abs(num - one[:, ik, iq, it, ip]).max())
    report("closed form vs R1 operator, one site, 72 cases", dev)

    dev = 0.0
    for k, q, t, p in product(KS, QS, THETAS, PHIS):
        th, ph = np.radians(t), np.radians(p)
        k2, q2 = np.array(SC_TWO["A"]) @ np.array(k), np.array(SC_TWO["A"]) @ np.array(q)    # same Cartesian k, q
        two = np.linalg.eigvalsh(gbt_h1(SC_TWO, k2, q2, *two_site_angles(th, ph, False)))
        fold = np.sort(np.concatenate([closed_form(k, q, th), closed_form(np.array(k) + [0.5, 0, 0], q, th)]))
        dev = max(dev, np.abs(two - fold).max())
    report("doubled cell vs two folded one-site spectra, same Cartesian q, 72 cases", dev)

    for n in (2, 3, 4):
        dev1 = dev_p = dev_a = 0.0
        for k, t, p in product(KS, THETAS, PHIS):
            th, ph = np.radians(t), np.radians(p)
            kc = cart(k, SC_ONE)
            ev = supercell_spectrum(kc, [1.0 / n, 0, 0], n, [0] * n, [th], [ph])
            dev1 = max(dev1, np.abs(ev - folded_one_site(k, n, th)).max())
            for anti, which in ((False, "p"), (True, "a")):
                ths, phs = two_site_angles(th, ph, anti)
                qc = [1.0 / (2 * n), 0, 0]                     # q = (1/n,0,0) direct in the doubled cell
                ev = supercell_spectrum(cart(k, SC_TWO), qc, 2 * n, [i % 2 for i in range(2 * n)], ths, phs)
                d2 = np.abs(ev - folded_two_site(k, n, th, ph, anti)).max()
                if anti:
                    dev_a = max(dev_a, d2)
                else:
                    dev_p = max(dev_p, d2)
        report(f"supercell vs closed form, one site, q = (1/{n},0,0), 24 (k,theta,phi)", dev1)
        report(f"supercell vs GBT, two sites parallel, q = (1/{n},0,0)", dev_p)
        report(f"supercell vs GBT, two sites antiparallel, q = (1/{n},0,0)", dev_a)
    return ok


def write_o1(outdir, one, par, anti):
    nk, nq, nt, npf = len(KS), len(QS), len(THETAS), len(PHIS)
    s = "! O1: one-orbital spiral model, k-space GBT operator (spec R1). Generated by make_gbt_references.py; do not edit.\n"
    s += "&o1_model\n"
    s += f"  t = {fmt(T)}\n  delta = {fmt(DELTA)}\n  w_up = 1.0\n  w_dn = 1.0\n"
    s += f"  c_majority = {fmt(-0.5 * DELTA)}\n  c_minority = {fmt(0.5 * DELTA)}\n"
    s += "  s_rule = 'S_ab(k) = t sum over nearest-neighbour bonds d (|d| = 1, d = R + tau_b - tau_a) of exp(2 pi i k.d); k Cartesian in 2 pi/alat'\n"
    s += "  frame_rule = 'U_a = D(phi_a) R_y(theta_a) D(-phi_a), D(x) = diag(exp(-ix/2), exp(+ix/2)); H_aa carries C = diag(c_majority, c_minority) in the local frame'\n"
    s += "  eigenvalue_order = 'ascending, Ry'\n/\n"
    s += "&o1_cases\n"
    s += f"  nk = {nk}\n  nq = {nq}\n  ntheta = {nt}\n  nphi = {npf}\n"
    s += "  k_rule = 'k and q are direct coordinates of the cell the eigenvalues belong to (one-site: simple cubic, a = 1; two-site: doubled cell below)'\n"
    s += arr("k_direct", np.array(KS).T) + arr("q_direct", np.array(QS).T)
    s += arr("theta_deg", THETAS) + arr("phi_deg", PHIS) + "/\n"
    s += "&o1_one_site\n  eig_shape = '2, nk, nq, ntheta, nphi (Fortran order)'\n"
    s += "  nsite = 1\n  a_cell = 1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0\n  tau = 0.0, 0.0, 0.0\n"
    s += "  cart_rule = 'k_cart, q_cart in 2 pi/alat, one column per k or q'\n"
    s += arr("k_cart", np.array([cart(k, SC_ONE) for k in KS]).T) + arr("q_cart", np.array([cart(q, SC_ONE) for q in QS]).T)
    s += arr("eig", one) + "/\n"
    bl = bonds(SC_TWO)
    s += "&o1_two_site\n  eig_shape = '4, nk, nq, ntheta, nphi (Fortran order)'\n  nsite = 2\n"
    s += "  a_cell = 2.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0\n  tau = 0.0, 0.0, 0.0,  1.0, 0.0, 0.0\n"
    s += "  site_rule_parallel = 'theta_a = theta_ss, phi_a = phi_ss on both sites'\n"
    s += "  site_rule_antiparallel = 'site 1: theta_ss, phi_ss; site 2: theta = pi - theta_ss, phi = phi_ss + pi; W and C as for a majority site'\n"
    s += f"  nbond = {len(bl)}\n  bond_s = {fmt(T)}\n"
    s += "  bond_site_a = " + ", ".join(str(b[0] + 1) for b in bl) + "\n"
    s += "  bond_site_b = " + ", ".join(str(b[1] + 1) for b in bl) + "\n"
    s += "  bond_d = " + ",  ".join(", ".join(f"{x:.1f}" for x in b[3]) for b in bl) + "\n"
    s += "  cart_rule = 'k_cart, q_cart in 2 pi/alat, one column per k or q; the direct coordinates in &o1_cases are those of this doubled cell'\n"
    s += arr("k_cart", np.array([cart(k, SC_TWO) for k in KS]).T) + arr("q_cart", np.array([cart(q, SC_TWO) for q in QS]).T)
    s += arr("eig_parallel", par) + arr("eig_antiparallel", anti) + "/\n"
    Path(outdir, "o1_one_band.nml").write_text(s)


NK_MESH = 24
KT = 0.002
O2_SETS = [(0.6, 0.0), (0.9, 0.0), (0.6, 0.5)]      # (electrons per site, rigid shift of both C values in Ry)
O2_QS = [(1 / 12, 0, 0), (1 / 6, 0, 0), (1 / 4, 0, 0), (1 / 2, 0, 0), (1 / 6, 1 / 6, 0)]
O2_THETAS = [0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 20.0]  # deg
N_RICH = 3                                           # Richardson uses the three smallest angles
FT_CHECK_REL = 1.0e-6                                # residual / max_q |J(0) - J(q)|
MESH_K = np.meshgrid(*[np.arange(NK_MESH) / NK_MESH] * 3, indexing="ij")


def band_mesh(k):
    return 2.0 * T * sum(np.cos(2.0 * np.pi * c) for c in k)


def eig_mesh(q, theta, shift):
    sm = band_mesh([c - 0.5 * x for c, x in zip(MESH_K, q)])
    sp = band_mesh([c + 0.5 * x for c, x in zip(MESH_K, q)])
    bar, dl = 0.5 * (sm + sp), 0.5 * (sm - sp)
    root = np.hypot(dl * np.sin(theta), dl * np.cos(theta) - 0.5 * DELTA)
    return bar + shift - root, bar + shift + root


def fermi(e, mu):
    return 0.5 * (1.0 - np.tanh((e - mu) / (2.0 * KT)))


def ft_quantities(q, theta, mu, shift):
    """Per site: Omega (Mermin, fixed mu), N, band energy at fixed mu without -mu N."""
    om = nn = eb = 0.0
    for e in eig_mesh(q, theta, shift):
        om += -KT * np.logaddexp(0.0, -(e - mu) / KT).sum()
        nn += fermi(e, mu).sum()
        eb += (fermi(e, mu) * e).sum()
    n_k = NK_MESH ** 3
    return om / n_k, nn / n_k, eb / n_k


def solve_mu(n_target, shift):
    lo, hi = shift - 1.0, shift + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ft_quantities((0, 0, 0), 0.0, mid, shift)[1] < n_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lkag_j(q, mu, shift):
    ea = band_mesh(MESH_K) - 0.5 * DELTA + shift
    eb = band_mesh([c + x for c, x in zip(MESH_K, q)]) + 0.5 * DELTA + shift
    fa, fb, d = fermi(ea, mu), fermi(eb, mu), eb - ea
    safe = np.abs(d) > 1.0e-9
    term = np.where(safe, (fa - fb) / np.where(safe, d, 1.0), fa * (1.0 - fa) / KT)
    return 0.25 * DELTA ** 2 * term.sum() / NK_MESH ** 3


def neville(x, y):
    p = list(y)
    for m in range(1, len(x)):
        p = [(x[i + m] * p[i] - x[i] * p[i + 1]) / (x[i + m] - x[i]) for i in range(len(p) - 1)]
    return p[0]


def compute_o2(n_target, shift):
    mu = solve_mu(n_target, shift)
    th = np.radians(O2_THETAS)
    nq, nt = len(O2_QS), len(th)
    omega, nel, eband = (np.zeros((nq, nt)) for _ in range(3))
    for iq, q in enumerate(O2_QS):
        for it, t in enumerate(th):
            omega[iq, it], nel[iq, it], eband[iq, it] = ft_quantities(q, t, mu, shift)
    ref = np.array([ft_quantities((0, 0, 0), t, mu, shift) for t in th])
    omega_q0 = np.array([ft_quantities(q, 0.0, mu, shift)[0] for q in O2_QS])
    j0 = lkag_j((0, 0, 0), mu, shift)
    jdiff = np.array([j0 - lkag_j(q, mu, shift) for q in O2_QS])
    return dict(mu=mu, shift=shift, omega=omega, nel=nel, eband=eband, omega_ref=ref[:, 0], nel_ref=ref[:, 1],
                eband_ref=ref[:, 2], omega_q0=omega_q0, omega_00=ft_quantities((0, 0, 0), 0.0, mu, shift)[0],
                jdiff=jdiff)


def richardson(res, key, ref_key, idx):
    s2 = np.sin(np.radians(O2_THETAS)) ** 2
    ratio = (res[key] - res[ref_key][None, :]) / s2[None, :]
    return np.array([neville(s2[idx], ratio[iq, idx]) for iq in range(len(O2_QS))])


def o2_checks(sets):
    ok = True
    small, mid, wide = list(range(N_RICH)), [2, 3, 4], list(range(len(O2_THETAS)))
    for idx_set, ((n_target, shift), res) in enumerate(zip(O2_SETS, sets)):
        tag = f"set {idx_set + 1} (N = {n_target}, shift {shift})"
        dev = max(abs(res["omega_q0"] - res["omega_00"]))
        good = dev < 1.0e-14
        ok &= good
        print(f"  {tag}: Omega(q,0) - Omega(0,0), 5 q: max dev {dev:.3e} Ry  {'ok' if good else 'FAIL'}")
        print(f"  {tag}: Omega(0,theta) - Omega(0,0): {np.abs(res['omega_ref'] - res['omega_00']).max():.3e} Ry  (by construction, R6, not counted)")
        scale = np.abs(res["jdiff"]).max()
        r = richardson(res, "omega", "omega_ref", small)
        resid = np.abs(r - res["jdiff"])
        good = resid.max() / scale < FT_CHECK_REL
        ok &= good
        print(f"  {tag}: Richardson ({O2_THETAS[0]}, {O2_THETAS[1]}, {O2_THETAS[2]} deg) vs J(0)-J(q): residual/max|J0-Jq| {resid.max() / scale:.3e};"
              f" relative to J0-Jq per q: {', '.join(f'{x:.1e}' for x in resid / np.abs(res['jdiff']))}  {'ok' if good else 'FAIL'}")
        for idx, label in ((mid, "2, 4, 10 deg (information)"), (wide, "all seven angles (information)")):
            rr = np.abs(richardson(res, "omega", "omega_ref", idx) - res["jdiff"]) / scale
            print(f"  {tag}: Richardson ({label}): residual/max|J0-Jq| {rr.max():.3e}")
        mut = np.abs(richardson(res, "eband", "eband_ref", small) - res["jdiff"]).max() / scale
        good = mut > 10.0 * FT_CHECK_REL
        ok &= good
        print(f"  {tag}: mutation (band energy at fixed mu, no -mu N): violation/max|J0-Jq| {mut:.3e}  "
              f"{'detected, ok' if good else 'NOT DETECTED, FAIL'}")
    one, three = sets[0], sets[2]
    dev_o = max(np.abs(three["omega"] - one["omega"]).max(), np.abs(three["omega_ref"] - one["omega_ref"]).max())
    dev_n = np.abs(three["nel"] - one["nel"]).max()
    dev_j = np.abs(three["jdiff"] - one["jdiff"]).max()
    dev_e = max(np.abs(three["eband"] - one["eband"] - 0.5 * one["nel"]).max(),
                np.abs(three["eband_ref"] - one["eband_ref"] - 0.5 * one["nel_ref"]).max())
    for label, dev, tol in (("Omega", dev_o, 1.0e-14), ("N", dev_n, 1.0e-14), ("J(0)-J(q)", dev_j, 1.0e-12),
                            ("band energy minus 0.5 N", dev_e, 1.0e-14)):
        good = dev < tol
        ok &= good
        print(f"  set 3 vs set 1, {label}: max dev {dev:.3e}  {'ok' if good else 'FAIL'}")
    return ok


def write_o2(outdir, sets):
    s = "! O2: force theorem on the O1 one-site model (spec R7, R8, R11). Generated by make_gbt_references.py; do not edit.\n"
    s += "&o2_settings\n"
    s += f"  nk_mesh = {NK_MESH}, {NK_MESH}, {NK_MESH}\n  k_offset = 0.0, 0.0, 0.0\n  kt_ry = {fmt(KT)}\n"
    s += "  mesh_rule = 'Gamma-centred full mesh, weights 1/nk^3, no symmetry reduction; k direct = Cartesian (2 pi/alat); all q/2 are mesh points'\n"
    s += f"  nq = {len(O2_QS)}\n  ntheta = {len(O2_THETAS)}\n  n_richardson = {N_RICH}\n"
    s += arr("q_direct", np.array(O2_QS).T) + arr("q_cart", np.array(O2_QS).T) + arr("theta_deg", O2_THETAS)
    s += "  richardson_rule = 'Neville extrapolation in sin^2(theta) of [omega(q,theta) - omega_ref(theta)]/sin^2(theta) over the first n_richardson angles'\n"
    s += "  array_rule = 'omega, n, eband_fixed_mu: (nq, ntheta) Fortran order; *_q0: Omega(q, 0) per q; *_ref: q = 0 per theta; Ry per site'\n"
    s += "  lkag_rule = 'j0_minus_jq = J(0) - J(q), H = -sum_{i/=j} J_ij e_i.e_j, band (Lehmann) sum of the LKAG product, same mesh, kT and mu'\n"
    s += "  mutation_rule = 'eband_fixed_mu = sum_k w f(eps - mu) eps without the -mu N term (the fixed-E_F band energy of finding F2)'\n"
    s += "  set_rule = 'set 1: N = 0.6; set 2: N = 0.9; set 3: set 1 with both C values shifted by c_shift_ry (Omega, N and J(0)-J(q) unchanged, eband_fixed_mu larger by c_shift_ry N)'\n/\n"
    for idx, ((n_target, shift), res) in enumerate(zip(O2_SETS, sets)):
        s += f"&o2_set{idx + 1}\n  n_electrons = {n_target}\n  c_shift_ry = {fmt(shift)}\n  mu_ry = {fmt(res['mu'])}\n"
        s += arr("omega", res["omega"]) + arr("n", res["nel"]) + arr("eband_fixed_mu", res["eband"])
        s += arr("omega_ref", res["omega_ref"]) + arr("n_ref", res["nel_ref"]) + arr("eband_fixed_mu_ref", res["eband_ref"])
        s += arr("omega_q0", res["omega_q0"]) + f"  omega_00 = {fmt(res['omega_00'])}\n"
        s += arr("j0_minus_jq", res["jdiff"]) + "/\n"
    Path(outdir, "o2_force_theorem.nml").write_text(s)


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    one, par, anti = compute_sets()
    print("O1 self-checks:")
    if not self_checks(one):
        sys.exit("O1 self-checks FAILED; nothing written")
    sets = [compute_o2(n, shift) for n, shift in O2_SETS]
    print("O2 self-checks:")
    for i, ((n, shift), res) in enumerate(zip(O2_SETS, sets)):
        print(f"  set {i + 1}: N = {n}, shift {shift}: mu = {res['mu']:.12f} Ry")
    if not o2_checks(sets):
        sys.exit("self-checks FAILED; nothing written")
    write_o1(outdir, one, par, anti)
    write_o2(outdir, sets)
    print(f"wrote {outdir / 'o1_one_band.nml'} and {outdir / 'o2_force_theorem.nml'}")


if __name__ == "__main__":
    main()
