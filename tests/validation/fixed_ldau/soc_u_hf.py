#!/usr/bin/env python3
"""nsp=2 SOC+U Hellmann-Feynman tests (bcc Fe, flag on, local or spectral). Run with the repository .venv interpreter.

  soc_u_hf.py run   --work DIR --kind ks|rs NAME...     NAME = w1_{bp,bm,s1p,s1m} | U<hundredths> | Us<hundredths> | J<hundredths>
  soc_u_hf.py analyze --work DIR --kind ks|rs

Decks are fresh copies of tests/scf/cases (bulk/bccFe or k_space_scf/bccFe) with the listed edits. W1 runs are U=J=0
with an identity test operator on the d block (both spins, or spin 1 for s1*), lambda = +-0.005 Ry. F is candidate_b_F (KS,
free energy) or candidate_b (RS) plus lambda*Tr(W P) with Tr(W P) from the d occupations ql (the logged e_external is zero at U=J=0).
Target dE_Hub/dX at fixed P: hubbard_spinor_functional at X +- eps through soc_u_fixed_p_derivative (build it with
gfortran -Ibuild/modules soc_u_fixed_p_derivative.f90 build/source/librslmto.a -framework Accelerate -o WORK/fixedp).
Allowance A_X = 2 sum_s |g_s dq_s/dX| + |D(h1) - D(h2)|, g_s from R^T g = B (R[s,a] = dq_s/dlambda_a, B_a = dF/dlambda_a - q_a).
"""
import argparse, json, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
import ldau_regression as R

RY = 13.605703976
LAM = 0.005


def spec(name):
    """name -> (nsp-2 keys: u, j, mode, lambda, spin)"""
    if name.startswith("w1_"):
        tag = name[3:]
        return 0.0, 0.0, "local", (LAM if tag.endswith("p") else -LAM), (1 if tag.startswith("s1") else 0)
    m = re.fullmatch(r"(Us|U|J)(\d+)", name)
    v = int(m.group(2)) / 100.0
    return {"U": (v, 0.0, "local"), "Us": (v, 0.0, "spectral"), "J": (2.0, v, "local")}[m.group(1)] + (0.0, 0)


def deck(kind):
    def edit_for(u, j, mode, lam, spin):
        def edit(t):
            if kind == "ks":
                t = R.use_kspace(t).replace("&control", "&control\nlmax = 2", 1)
                t = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 400", t)
                t = re.sub(r"(?m)^fermi\s*=.*", "fermi = -0.082950", t)
            else:
                t = R.set_nstep(t, 150)
                t = re.sub(r"(?m)^fermi\s*=.*", "fermi = -0.085033", t)
                t = re.sub(r"(?m)^beta\s*=.*", "beta = 0.005", t)
                t = t.replace("&self", "&self\nconv_thr = 1.0e-30", 1)
            t = re.sub(r"(?m)^nsp\s*=\s*\d", "nsp = 2", t)
            keys = (f"hubbard_u_general(1,3) = {u}\nhubbard_j_general(1,3) = {j}\nhubbard_u_potential_form = 'liechtenstein'\n"
                    f"hubbard_radial_moments = '{mode}'\nhubbard_soc_experimental = .true.\n")
            if lam:
                keys += f"hubbard_test_operator = 'w1'\nhubbard_test_lambda = {lam}\nhubbard_test_spin = {spin}\n"
            return R.in_hamiltonian(t, keys) if "&hamiltonian" in t else t + "\n&hamiltonian\n" + keys + "/\n"
        return edit
    return edit_for


def parse(d):
    t = R.log_text(d)
    last = re.findall(r"HUBBARD_MOMENT_ENERGY_TOTAL (.*)", t)[-1]
    out = {k: float(re.search(rf"\b{k}=\s*([-0-9.EeD+]+)", last).group(1).replace("D", "E")) for k in ("E_Hub",)}
    for k in ("candidate_b", "candidate_b_F", "candidate_c", "candidate_c_F"):
        m = re.search(rf"\b{k}=\s*([-0-9.EeD+]+)", last)
        out[k] = float(m.group(1)) if m else None
    out["diff"] = float(re.findall(r"Diff=\s*([0-9.Ee+-]+)", t)[-1])
    ql = {}
    for line in (d / "Fe_out.nml").read_text().splitlines():
        m = re.match(r"\s*ql\((\d+),\s*:\s*,\s*(\d)\)\s*=\s*(.*)", line)
        if m:
            ql[(int(m.group(1)), int(m.group(2)))] = [float(x) for x in m.group(3).replace(",", " ").split()]
    out["d_up"], out["d_dn"] = ql[(1, 1)][2], ql[(1, 2)][2]
    return out


def run_all(args):
    work = Path(args.work)

    def one(name):
        d = R.prepare(R.CASES / ("bulk/bccFe" if args.kind == "rs" else "k_space_scf/bccFe"), work / f"{args.kind}_{name}",
                      deck(args.kind)(*spec(name)))
        rc = R.run(Path(args.binary), d, "Fe.nml" if args.kind == "rs" else None)
        return name, ({"rc": rc} if rc else {"rc": 0, **parse(d)})

    with ThreadPoolExecutor(4) as ex:
        res = dict(ex.map(one, args.names))
    f = work / f"{args.kind}_results.json"
    old = json.loads(f.read_text()) if f.exists() else {}
    old.update(res)
    f.write_text(json.dumps(old, indent=1))
    for k, v in res.items():
        print(k, v)


def sidecar_p(d):
    rows = {}
    for line in (d / "Fe_dm.nml").read_text().splitlines():
        m = re.match(r"\s*dm_(re|im)\(3,(\d),(\d),(\d),:\)\s*=\s*(.*)", line)
        if m:
            rows[(m.group(1), *map(int, m.group(2, 3, 4)))] = np.array([float(x) for x in m.group(5).replace(",", " ").split()])
    p = np.zeros((10, 10), complex)
    for sr in (1, 2):
        for sc in (1, 2):
            for m1 in range(1, 6):
                p[(sr - 1) * 5 + m1 - 1, (sc - 1) * 5:sc * 5] = rows[("re", sr, sc, m1)] + 1j * rows[("im", sr, sc, m1)]
    return p


def analyze(args):
    work = Path(args.work)
    r = json.loads((work / f"{args.kind}_results.json").read_text())
    kind = args.kind
    fkey, ckey = ("candidate_b_F", "candidate_c_F") if kind == "ks" else ("candidate_b", "candidate_c")

    def fixed_p(name, u_ev, j_ev):
        f = work / f"{kind}_{name}" / "P.txt"
        f.write_text("\n".join(f"{z.real:.17e} {z.imag:.17e}" for z in sidecar_p(work / f"{kind}_{name}").ravel()) + "\n")
        o = subprocess.run([str(work / "fixedp"), str(f), str(u_ev / RY), str(j_ev / RY), "1e-4"], capture_output=True, text=True, check=True).stdout.split()
        return float(o[1]), float(o[3]) / RY, float(o[5]) / RY  # E_Hub (Ry), dE/dU, dE/dJ (Ry/eV)

    q = lambda n: np.array([r[n]["d_up"], r[n]["d_dn"]])
    F = lambda n, lam, spins: r[n][fkey] + lam * sum(q(n)[s] for s in spins)
    q0 = 0.5 * (q("w1_bp") + q("w1_bm"))
    dF = {t: (F(f"w1_{t}p", LAM, s) - F(f"w1_{t}m", -LAM, s)) / (2 * LAM) for t, s in (("b", (0, 1)), ("s1", (0,)))}
    dq = {t: (q(f"w1_{t}p") - q(f"w1_{t}m")) / (2 * LAM) for t in ("b", "s1")}
    print(f"W1 both spins: dF/dlambda {dF['b']:.6f}  Tr(WP) {q0.sum():.6f}  relative {(dF['b'] - q0.sum()) / q0.sum():+.4%}")
    resp = np.array([dq["s1"], dq["b"] - dq["s1"]]).T  # resp[s,a] = dq_s/dlambda_a
    g = np.linalg.solve(resp.T, np.array([dF["s1"] - q0[0], dF["b"] - dF["s1"] - q0[1]]))
    print(f"g_s (Ry/electron): {g[0]:.6f} {g[1]:.6f}   dq_s/dlambda_a: {resp.tolist()}")

    def allowance(lo, hi, h):
        return 2 * np.sum(np.abs(g * (q(hi) - q(lo)) / (2 * h)))

    for X, lo_hi in (("U", [(0.05, "U195", "U205"), (0.10, "U190", "U210")]), ("J", [(0.10, "J040", "J060")])):
        D, T, A = [], [], []
        for h, lo, hi in lo_hi:
            if lo not in r or hi not in r:
                continue
            xs = (2 - h, 2 + h) if X == "U" else (0.5 - h, 0.5 + h)
            args_ = ((xs[0], 0.0), (xs[1], 0.0)) if X == "U" else ((2.0, xs[0]), (2.0, xs[1]))
            D.append((r[hi][fkey] - r[lo][fkey]) / (2 * h))
            T.append(0.5 * (fixed_p(lo, *args_[0])[1 if X == "U" else 2] + fixed_p(hi, *args_[1])[1 if X == "U" else 2]))
            A.append(allowance(lo, hi, h))
        if not D:
            continue
        a_x = np.mean(A) + (abs(D[0] - D[1]) if len(D) == 2 else 0.0)
        print(f"{X} b: steps {len(D)}  D {np.mean(D):.8f}  target {np.mean(T):.8f}  miss {np.mean(D) - np.mean(T):+.3e}  A_X {a_x:.3e}"
              f"{'' if len(D) == 2 else ' (single step: no |D(h1)-D(h2)| term)'}  b within A_X: {abs(np.mean(D) - np.mean(T)) <= a_x}")
        if X == "U" and "Us190" in r and "Us210" in r:
            Dc = (r["Us210"][ckey] - r["Us190"][ckey]) / 0.2
            Tc = 0.5 * (fixed_p("Us190", 1.9, 0.0)[1] + fixed_p("Us210", 2.1, 0.0)[1])
            print(f"U c (spectral, h=0.10): D {Dc:.8f}  target {Tc:.8f}  miss {Dc - Tc:+.3e}  miss/A_X {abs(Dc - Tc) / a_x:.1f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("run", "analyze"))
    ap.add_argument("--work", required=True)
    ap.add_argument("--kind", required=True, choices=("ks", "rs"))
    ap.add_argument("--binary", default=str(ROOT / "build/bin/rslmto.x"))
    ap.add_argument("names", nargs="*")
    a = ap.parse_args()
    run_all(a) if a.cmd == "run" else analyze(a)
