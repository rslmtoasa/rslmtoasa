#!/usr/bin/env python3
"""Model references O1-O4 for the Stage 1 BdG kernels (bdg_assemble, bdg_accumulate, bdg_gap).

Developer-owned oracle: drafted on request, reviewed by the developer,
committed under tests/bdg/oracles/, read-only to AI sessions afterwards.

Where the expected values come from
-----------------------------------
Every reference value is computed WITHOUT the spec's doubled 4-component
Nambu construction:
  * gap and Tc values: scalar (or 2x2 multi-channel) BCS gap equations,
    solved by bisection / Newton on the same discrete k mesh;
  * electron counts, spin densities and anomalous densities: the reduced
    2x2 Nambu basis Phi_k = (c_{k,up}, c^dag_{-k,dn}), diagonalized per k;
  * O3 spectra: closed form.
The spec's doubled construction appears only in self-check 'spec', which
confirms that the spec formulas (section 1) reproduce these references.
If that self-check fails, the spec is wrong, not the references.

Conventions (spec section 1, restated)
--------------------------------------
  energies in Ry; k fractional; mesh k_n = (n - 1/2)/nk, n = 1..nk, w = 1/nk,
  closed under k -> -k (index n -> nk + 1 - n);
  pairing term  sum_iL Delta_iL c^dag_{iL up} c^dag_{iL dn} + h.c.;
  chi_iL = <c_{iL up} c_{iL dn}>;   Delta_iL = sum_L' K_LL' chi_iL',  K > 0 attractive;
  with a real positive seed, Delta stays real and positive in every model here.
  Code basis per site: index = L + (sigma - 1)*norb, sigma = 1 up, 2 dn.

Models (one site per cell, 1-D chain, no inter-orbital hopping)
---------------------------------------------------------------
  O1   single band, eps(k) = -2 t cos(2 pi k), K = g.  Delta at kT = 0, 0.5 Tc,
       0.9 Tc; Tc; at 1.1 Tc Delta = 0 and the linearized contraction rate.
  O1b  three degenerate copies of the O1 band, K = (g/3) * ones(3,3).
       Exact: every Delta_L equals O1's kT = 0 value.  K applied as g/3 * I
       gives Delta = 0 (g/3 is below the mesh's critical coupling).
  O2   O1 with a Zeeman field: eps_up = eps - h, eps_dn = eps + h.
       Case a: kT = 0, h = Delta0/2.  Exact: Delta = Delta0, N = N0, m = 0.
       Case b: kT = 0.5 Tc, h = Delta0/2, from the field gap equation.
  O3   k-odd, time-reversal-symmetric spin splitting:
       eps_up(k) = eps(k) + a sin(2 pi k), eps_dn(k) = eps(k) - a sin(2 pi k).
       Correct hole block (-H(-k)*) pairs eps_up(k) with eps_dn(-k) = eps_up(k);
       spectrum +-E'(k), +-E''(k) with xi' = eps + a s - mu, xi'' = eps - a s - mu.
       A hole block built from -H(k)* gives a s +- sqrt(xi^2 + Delta^2) instead
       (written as informational 'wrong' eigenvalues; spectra are taken at the mesh
       points nearest the two Fermi points, where the difference is O(alpha)).
       Self-consistently, the wrong block breaks the pairing: Delta -> 0.
  O4   two decoupled bands a, b (norb = 2), K = [[g_aa, g_ab], [g_ba, 0]].
       Case a: symmetric K.  Case b: asymmetric K (kernel check only; catches
       a transposed K).  Band b has no own coupling: Delta_b is induced.

Output: one Fortran namelist file per model, values in Ry.  Only numpy.

Usage: python3 make_references.py [outdir]
"""
import sys
from pathlib import Path

import numpy as np

NK = 2000
T_HOP = 0.05
MU = -0.03
G1 = 0.07            # O1 coupling: Delta0 ~ 5 mRy on this mesh
ALPHA_O3 = 0.02      # O3 k-odd splitting
DELTA_O3_FIXED = 0.004   # O3 fixed-Delta spectrum points
K_O3_GENERIC = 501      # 1-based mesh index of a k point away from the Fermi points
T_B, C_B = 0.03, 0.01    # O4 band b: eps_b = -2 t_b cos + c_b
G_AA, G_AB_SYM = 0.07, 0.02
G_AB_ASYM, G_BA_ASYM = 0.03, 0.01

CHECKS = []


def check(name, ok, detail):
    CHECKS.append((name, bool(ok), detail))


# ----------------------------------------------------------------------------
# mesh and bands
# ----------------------------------------------------------------------------
def mesh(nk=NK):
    k = (np.arange(1, nk + 1) - 0.5) / nk
    return k, np.full(nk, 1.0 / nk)


def band(k, t=T_HOP, c=0.0):
    return -2.0 * t * np.cos(2.0 * np.pi * k) + c


def minus_k(arr):
    """Value at -k for an array on the symmetric mesh."""
    return arr[::-1]


def theta(e, h, kt):
    """1 - f(E+h) - f(E-h) = (tanh((E+h)/2kT) + tanh((E-h)/2kT))/2; kT = 0 -> step."""
    if kt == 0.0:
        return 0.5 * (np.sign(e + h) + np.sign(e - h))
    return 0.5 * (np.tanh(0.5 * (e + h) / kt) + np.tanh(0.5 * (e - h) / kt))


def fermi(e, kt):
    if kt == 0.0:
        return np.where(e < 0.0, 1.0, np.where(e > 0.0, 0.0, 0.5))
    return 0.5 * (1.0 - np.tanh(0.5 * e / kt))


# ----------------------------------------------------------------------------
# scalar / multi-channel gap equations (primary references)
# ----------------------------------------------------------------------------
def chi_scalar(xi, w, delta, kt, h=0.0):
    """chi = sum_k w Delta/(2E) theta(E, h, kT) for one channel (spin-singlet pair)."""
    e = np.sqrt(xi * xi + delta * delta)
    return float(np.sum(w * delta / (2.0 * e) * theta(e, h, kt)))


def solve_gap(xis, w, kmat, kt, h=0.0, seed=None, tol=1e-15, maxit=200):
    """Solve Delta = K chi(Delta) for decoupled channels xis, by Newton with a numeric Jacobian."""
    n = len(xis)
    d = np.full(n, 0.005) if seed is None else np.array(seed, float)
    for _ in range(maxit):
        chi = np.array([chi_scalar(xis[i], w, d[i], kt, h) for i in range(n)])
        r = d - kmat @ chi
        if np.max(np.abs(r)) < tol:
            break
        jac = np.eye(n)
        for j in range(n):
            dd = d.copy()
            step = 1e-7 * max(abs(d[j]), 1e-6)
            dd[j] += step
            chi2 = np.array([chi_scalar(xis[i], w, dd[i], kt, h) for i in range(n)])
            jac[:, j] = ((dd - kmat @ chi2) - r) / step
        d = d - np.linalg.solve(jac, r)
    chi = np.array([chi_scalar(xis[i], w, d[i], kt, h) for i in range(n)])
    return d, chi, float(np.max(np.abs(d - kmat @ chi)))


def tc_single(xi, w, g):
    """kTc from 1 = g sum_k w tanh(xi/2kT)/(2 xi), by bisection."""
    def lin(kt):
        return g * float(np.sum(w * np.tanh(0.5 * xi / kt) / (2.0 * xi)))
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if lin(mid) > 1.0 else (lo, mid)
    return 0.5 * (lo + hi), lin


# ----------------------------------------------------------------------------
# reduced 2x2 Nambu basis (independent observables)
# ----------------------------------------------------------------------------
def reduced_observables(eps_up, eps_dn, w, delta, mu, kt):
    """Phi_k = (c_{k up}, c^dag_{-k dn}); H = [[eps_up(k)-mu, D], [D, -(eps_dn(-k)-mu)]].

    Returns n_up, n_dn, chi = <c_up c_dn> (on-site), all per site.
    """
    a = eps_up - mu
    b = -(minus_k(eps_dn) - mu)
    n_up = n_dn = chi = 0.0
    for i in range(len(w)):
        hm = np.array([[a[i], delta], [delta, b[i]]])
        ev, vec = np.linalg.eigh(hm)
        f = fermi(ev, kt)
        u, v = vec[0, :], vec[1, :]
        n_up += w[i] * np.sum(np.abs(u) ** 2 * f)
        n_dn += w[i] * (1.0 - np.sum(np.abs(v) ** 2 * f))     # state -k dn
        chi += w[i] * (-np.sum(u * np.conj(v) * f)).real       # <c_up c_dn> = -<c_dn c_up>
    return n_up, n_dn, chi


# ----------------------------------------------------------------------------
# spec construction (self-check only)
# ----------------------------------------------------------------------------
def spec_observables(hk_list, hmk_list, w, deltas, mu, kt, norb):
    """Spec section 1: doubled Nambu, D = Delta (x) i sigma_y; chi and N by the spec formulas."""
    nb = 2 * norb
    chi = np.zeros(norb, complex)
    nel = 0.0
    eigs = []
    for hk, hmk, wk in zip(hk_list, hmk_list, w):
        d = np.zeros((nb, nb), complex)
        for l in range(norb):
            d[l, l + norb] = deltas[l]
            d[l + norb, l] = -deltas[l]
        hb = np.block([[hk - mu * np.eye(nb), d], [d.conj().T, -hmk.conj() + mu * np.eye(nb)]])
        ev, vec = np.linalg.eigh(hb)
        eigs.append(ev)
        u, v = vec[:nb, :], vec[nb:, :]
        f = fermi(ev, kt)
        for l in range(norb):
            chi[l] += wk * np.sum(u[l, :] * np.conj(v[l + norb, :]) * (1.0 - f))
        nel += wk * np.sum(f * np.sum(np.abs(u) ** 2, axis=0))
    return chi, nel, eigs


def diag_h(eps_by_index):
    """Spin-diagonal H(k) list from per-(orbital, spin) band arrays, code index order."""
    nb = len(eps_by_index)
    nk = len(eps_by_index[0])
    return [np.diag([eps_by_index[a][i] for a in range(nb)]).astype(complex) for i in range(nk)]


# ----------------------------------------------------------------------------
# namelist writer
# ----------------------------------------------------------------------------
def nml_value(v):
    if isinstance(v, str):
        return "'" + v + "'"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.17e}"
    arr = np.asarray(v, dtype=float).ravel(order="F")
    lines, row = [], []
    for x in arr:
        row.append(f"{x:.17e}")
        if len(row) == 4:
            lines.append("    " + ", ".join(row) + ",")
            row = []
    if row:
        lines.append("    " + ", ".join(row) + ",")
    lines[-1] = lines[-1].rstrip(",")
    return "\n" + "\n".join(lines)


def write_nml(path, groups):
    with open(path, "w") as fh:
        for gname, items in groups:
            fh.write(f"&{gname}\n")
            for key, val in items:
                fh.write(f"  {key} = {nml_value(val)}\n")
            fh.write("/\n\n")


# ----------------------------------------------------------------------------
# models
# ----------------------------------------------------------------------------
def model_o1(k, w):
    eps = band(k)
    xi = eps - MU
    ktc, lin = tc_single(xi, w, G1)
    out = {"ktc": ktc}
    for tag, kt in (("t0", 0.0), ("t05", 0.5 * ktc), ("t09", 0.9 * ktc)):
        d, chi, res = solve_gap([xi], w, np.array([[G1]]), kt)
        n_up, n_dn, chi_red = reduced_observables(eps, eps, w, d[0], MU, kt)
        n0u, n0d, _ = reduced_observables(eps, eps, w, 0.0, MU, kt)
        check(f"o1_{tag}_residual", res < 1e-13, f"{res:.1e}")
        check(f"o1_{tag}_chi_reduced", abs(chi_red - chi[0]) < 1e-12 * max(1, abs(chi[0])),
              f"{abs(chi_red - chi[0]):.1e}")
        check(f"o1_{tag}_spin_symmetric", abs(n_up - n_dn) < 1e-12, f"{abs(n_up - n_dn):.1e}")
        out[tag] = dict(kt=kt, delta=d[0], chi=chi[0], nel=n_up + n_dn, nel_normal=n0u + n0d)
    kt11 = 1.1 * ktc
    rate = lin(kt11)
    d, _, _ = solve_gap([xi], w, np.array([[G1]]), kt11, seed=[0.005])
    check("o1_t11_normal", abs(d[0]) < 1e-12 and rate < 1.0, f"Delta={d[0]:.1e}, rate={rate:.4f}")
    out["t11"] = dict(kt=kt11, delta=0.0, lin_rate=rate)
    check("o1_gap_ratio_finite_band", 3.0 < 2 * out["t0"]["delta"] / ktc < 4.0,
          f"2Delta0/kTc = {2 * out['t0']['delta'] / ktc:.4f} (weak-coupling 3.528)")
    return eps, out


def model_o1b(k, w, delta0):
    xi = band(k) - MU
    kmat = np.full((3, 3), G1 / 3.0)
    d, chi, res = solve_gap([xi, xi, xi], w, kmat, 0.0, seed=[0.003, 0.004, 0.005])
    check("o1b_equals_o1", np.max(np.abs(d - delta0)) < 1e-12, f"{np.max(np.abs(d - delta0)):.1e}")
    dbad, _, _ = solve_gap([xi], w, np.array([[G1 / 3.0]]), 0.0, seed=[0.005])
    check("o1b_diagonal_bug_gives_zero", abs(dbad[0]) < 1e-12, f"{dbad[0]:.1e}")
    return d, chi


def model_o2(k, w, o1):
    eps = band(k)
    xi = eps - MU
    d0 = o1["t0"]["delta"]
    h = 0.5 * d0
    res = {}
    # case a: kT = 0
    d, chi, r = solve_gap([xi], w, np.array([[G1]]), 0.0, h=h, seed=[0.5 * d0])
    nu, nd, chi_red = reduced_observables(eps - h, eps + h, w, d[0], MU, 0.0)
    check("o2a_delta_equals_delta0", abs(d[0] - d0) < 1e-13, f"{abs(d[0] - d0):.1e}")
    check("o2a_m_zero", abs(nu - nd) < 1e-12, f"{nu - nd:.1e}")
    check("o2a_n_equals_n0", abs(nu + nd - o1["t0"]["nel"]) < 1e-12, f"{abs(nu + nd - o1['t0']['nel']):.1e}")
    check("o2a_chi_reduced", abs(chi_red - chi[0]) < 1e-12, f"{abs(chi_red - chi[0]):.1e}")
    res["a"] = dict(kt=0.0, h=h, delta=d[0], chi=chi[0], n_up=nu, n_dn=nd)
    # case b: kT = 0.5 Tc; count nonzero roots on a Delta grid
    kt = 0.5 * o1["ktc"]
    grid = np.linspace(1e-6, 2.0 * d0, 4001)
    f = np.array([x - G1 * chi_scalar(xi, w, x, kt, h) for x in grid])
    roots = int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
    d, chi, r = solve_gap([xi], w, np.array([[G1]]), kt, h=h, seed=[d0])
    nu, nd, chi_red = reduced_observables(eps - h, eps + h, w, d[0], MU, kt)
    check("o2b_single_nonzero_root", roots == 1, f"{roots} sign change(s) on (0, 2 Delta0]")
    check("o2b_residual", r < 1e-13, f"{r:.1e}")
    check("o2b_chi_reduced", abs(chi_red - chi[0]) < 1e-12, f"{abs(chi_red - chi[0]):.1e}")
    check("o2b_m_positive", nu - nd > 0.0, f"m = {nu - nd:.3e}")
    res["b"] = dict(kt=kt, h=h, delta=d[0], chi=chi[0], n_up=nu, n_dn=nd)
    return res


def model_o3(k, w):
    eps = band(k)
    s = np.sin(2.0 * np.pi * k)
    eps_up, eps_dn = eps + ALPHA_O3 * s, eps - ALPHA_O3 * s
    check("o3_time_reversal", np.max(np.abs(minus_k(eps_dn) - eps_up)) < 1e-15, "eps_dn(-k) = eps_up(k)")
    xi1 = eps_up - MU
    d, chi, r = solve_gap([xi1], w, np.array([[G1]]), 0.0)
    nu, nd, chi_red = reduced_observables(eps_up, eps_dn, w, d[0], MU, 0.0)
    check("o3_residual", r < 1e-13, f"{r:.1e}")
    check("o3_chi_reduced", abs(chi_red - chi[0]) < 1e-12, f"{abs(chi_red - chi[0]):.1e}")
    # spectrum points: nearest mesh points to xi' = 0 and to xi0 = 0 (first half), and one generic point
    half = len(k) // 2
    k_o3 = (int(np.argmin(np.abs(xi1[:half]))) + 1, int(np.argmin(np.abs((eps - MU)[:half]))) + 1, K_O3_GENERIC)
    spectra, wrong = [], []
    dfix = DELTA_O3_FIXED
    for n in k_o3:
        i = n - 1
        xp = eps[i] + ALPHA_O3 * s[i] - MU
        xm = eps[i] - ALPHA_O3 * s[i] - MU
        ep, em = np.hypot(xp, dfix), np.hypot(xm, dfix)
        spectra.append(np.sort([-ep, ep, -em, em]))
        x0 = eps[i] - MU
        e0 = np.hypot(x0, dfix)
        wrong.append(np.sort([ALPHA_O3 * s[i] - e0, ALPHA_O3 * s[i] + e0,
                              -ALPHA_O3 * s[i] - e0, -ALPHA_O3 * s[i] + e0]))
    diff = max(np.max(np.abs(a - b)) for a, b in zip(spectra, wrong))
    check("o3_spectra_differ_from_wrong_block", diff > 1e-3, f"max {diff:.2e} Ry at the Fermi points")
    # self-consistent Delta with the wrong block: k-odd pair breaking, h_k = alpha sin(2 pi k) on xi0
    dw, _, _ = solve_gap([eps - MU], w, np.array([[G1]]), 0.0, h=ALPHA_O3 * s, seed=[d[0]])
    check("o3_wrong_block_delta_differs", abs(dw[0] - d[0]) > 1e-3, f"Delta_wrong = {dw[0]:.3e}, Delta = {d[0]:.3e}")
    return eps_up, eps_dn, dict(delta=d[0], chi=chi[0], n_up=nu, n_dn=nd, k_index=k_o3, delta_wrong=dw[0],
                                spectra=np.array(spectra).T, wrong=np.array(wrong).T)


def model_o4(k, w):
    xa = band(k) - MU
    xb = band(k, T_B, C_B) - MU
    out = {}
    for tag, gab, gba in (("a", G_AB_SYM, G_AB_SYM), ("b", G_AB_ASYM, G_BA_ASYM)):
        kmat = np.array([[G_AA, gab], [gba, 0.0]])
        d, chi, r = solve_gap([xa, xb], w, kmat, 0.0, seed=[0.005, 0.001])
        check(f"o4{tag}_residual", r < 1e-13, f"{r:.1e}")
        check(f"o4{tag}_induced", d[1] > 1e-4, f"Delta_b = {d[1]:.3e}")
        if tag == "b":
            dt, _, _ = solve_gap([xa, xb], w, kmat.T, 0.0, seed=[0.005, 0.001])
            check("o4b_transpose_differs", np.max(np.abs(dt - d)) > 1e-4, f"{np.max(np.abs(dt - d)):.2e}")
        out[tag] = dict(kmat=kmat, delta=d, chi=chi)
    return band(k), band(k, T_B, C_B), out


# ----------------------------------------------------------------------------
# spec self-check: doubled Nambu (spec section 1) vs the references
# ----------------------------------------------------------------------------
def spec_selfcheck(k, w, o1, o2, o3, eps_up3, eps_dn3, o4, eps_a, eps_b):
    eps = band(k)
    # O1, kT = 0.5 Tc
    r = o1["t05"]
    hk = diag_h([eps, eps])
    chi, nel, _ = spec_observables(hk, hk[::-1], w, [r["delta"]], MU, r["kt"], 1)
    check("spec_o1_chi", abs(chi[0] - r["chi"]) < 1e-12, f"{abs(chi[0] - r['chi']):.1e}")
    check("spec_o1_nel", abs(nel - r["nel"]) < 1e-11, f"{abs(nel - r['nel']):.1e}")
    # O2 case b (Zeeman, finite T)
    r = o2["b"]
    hk = diag_h([eps - r["h"], eps + r["h"]])
    chi, nel, _ = spec_observables(hk, hk[::-1], w, [r["delta"]], MU, r["kt"], 1)
    check("spec_o2b_chi", abs(chi[0] - r["chi"]) < 1e-12, f"{abs(chi[0] - r['chi']):.1e}")
    check("spec_o2b_nel", abs(nel - (r["n_up"] + r["n_dn"])) < 1e-11, f"{abs(nel - r['n_up'] - r['n_dn']):.1e}")
    # O3 spectra at fixed Delta, with the correct -H(-k)* hole block
    hk = diag_h([eps_up3, eps_dn3])
    _, _, eigs = spec_observables(hk, hk[::-1], w, [DELTA_O3_FIXED], MU, 0.0, 1)
    err = max(np.max(np.abs(eigs[n - 1] - o3["spectra"][:, j])) for j, n in enumerate(o3["k_index"]))
    check("spec_o3_spectra", err < 1e-13, f"{err:.1e}")
    _, _, eigs_bad = spec_observables(hk, hk, w, [DELTA_O3_FIXED], MU, 0.0, 1)
    err = max(np.max(np.abs(eigs_bad[n - 1] - o3["wrong"][:, j])) for j, n in enumerate(o3["k_index"]))
    check("spec_o3_wrong_block_reproduces_wrong", err < 1e-13, f"{err:.1e}")
    # O4 case b, two orbitals: code index order (a up, b up, a dn, b dn)
    r = o4["b"]
    hk = diag_h([eps_a, eps_b, eps_a, eps_b])
    chi, _, _ = spec_observables(hk, hk[::-1], w, r["delta"], MU, 0.0, 2)
    err = np.max(np.abs(chi - r["chi"]))
    check("spec_o4b_chi", err < 1e-12, f"{err:.1e}")


# ----------------------------------------------------------------------------
def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)
    k, w = mesh()
    eps, o1 = model_o1(k, w)
    d1b, chi1b = model_o1b(k, w, o1["t0"]["delta"])
    o2 = model_o2(k, w, o1)
    eps_up3, eps_dn3, o3 = model_o3(k, w)
    eps_a, eps_b, o4 = model_o4(k, w)
    spec_selfcheck(k, w, o1, o2, o3, eps_up3, eps_dn3, o4, eps_a, eps_b)

    common = [("nk", NK), ("mu", MU), ("mesh_rule", "k_n = (n - 1/2)/nk, n = 1..nk, w = 1/nk; -k_n = k_(nk+1-n)"),
              ("basis_rule", "one site; code index = L + (sigma-1)*norb, sigma 1 = up; H(k) diagonal")]
    write_nml(outdir / "o1_single_band.nml", [
        ("o1_inputs", common + [("norb", 1), ("t", T_HOP), ("g", G1),
                                ("band_rule", "eps(k) = -2 t cos(2 pi k), same for up and dn"),
                                ("k", k), ("eps", eps)]),
        ("o1_reference", [("ktc", o1["ktc"])] + [
            (f"{key}_{tag}", o1[tag][key]) for tag in ("t0", "t05", "t09")
            for key in ("kt", "delta", "chi", "nel", "nel_normal")] +
            [("kt_t11", o1["t11"]["kt"]), ("delta_t11", 0.0), ("lin_rate_t11", o1["t11"]["lin_rate"])]),
    ])
    write_nml(outdir / "o1b_uniform_k.nml", [
        ("o1b_inputs", common + [("norb", 3), ("kmat", np.full((3, 3), G1 / 3.0)), ("kt", 0.0),
                                 ("band_rule", "three degenerate copies of the O1 band (use o1 eps)")]),
        ("o1b_reference", [("delta", d1b), ("chi", chi1b)]),
    ])
    write_nml(outdir / "o2_zeeman.nml", [
        ("o2_inputs", common + [("norb", 1), ("g", G1),
                                ("band_rule", "eps_up = eps - h, eps_dn = eps + h, eps from o1")]),
        ("o2_reference", [(f"{key}_{c}", o2[c][key]) for c in ("a", "b")
                          for key in ("kt", "h", "delta", "chi", "n_up", "n_dn")]),
    ])
    write_nml(outdir / "o3_kodd.nml", [
        ("o3_inputs", common + [("norb", 1), ("g", G1), ("kt", 0.0), ("alpha", ALPHA_O3),
                                ("band_rule", "eps_up = eps + alpha sin(2 pi k), eps_dn = eps - alpha sin(2 pi k)"),
                                ("eps_up", eps_up3), ("eps_dn", eps_dn3),
                                ("delta_fixed", DELTA_O3_FIXED), ("spectrum_k_index", np.array(o3["k_index"], float))]),
        ("o3_reference", [("delta", o3["delta"]), ("chi", o3["chi"]), ("n_up", o3["n_up"]), ("n_dn", o3["n_dn"]),
                          ("spectra", o3["spectra"]), ("spectra_dims", "(4, n_index): sorted BdG eigenvalues at delta_fixed"),
                          ("wrong_spectra", o3["wrong"]),
                          ("wrong_delta", o3["delta_wrong"]),
                          ("wrong_note", "informational: hole block from -H(k)* instead of -H(-k)*")]),
    ])
    write_nml(outdir / "o4_two_band.nml", [
        ("o4_inputs", common + [("norb", 2), ("kt", 0.0), ("t_b", T_B), ("c_b", C_B),
                                ("band_rule", "orbital 1: eps (o1), orbital 2: -2 t_b cos(2 pi k) + c_b; no hybridization"),
                                ("eps_a", eps_a), ("eps_b", eps_b),
                                ("kmat_a", o4["a"]["kmat"]), ("kmat_b", o4["b"]["kmat"]),
                                ("kmat_note", "kmat_b is asymmetric on purpose: kernel check only")]),
        ("o4_reference", [("delta_a", o4["a"]["delta"]), ("chi_a", o4["a"]["chi"]),
                          ("delta_b", o4["b"]["delta"]), ("chi_b", o4["b"]["chi"])]),
    ])

    width = max(len(c[0]) for c in CHECKS)
    for name, ok, detail in CHECKS:
        print(f"{name:<{width}}  {'ok' if ok else 'FAIL'}  {detail}")
    print(f"Delta0 = {o1['t0']['delta']:.10e} Ry, kTc = {o1['ktc']:.10e} Ry")
    sys.exit(0 if all(c[1] for c in CHECKS) else 1)


if __name__ == "__main__":
    main()
