"""analyze.py : per (mesh, eta, xi) local maxima of tr L, bare-continuum onset, zero crossings of Re(1+U chi0) below onset."""
import glob, os, re, sys
import numpy as np
SP = os.getcwd()
RY = 13605.693123  # meV per Ry
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../spin_response/oracles/fe_lswt.dat")
rows = np.array([[float(x) for x in l.split()[:9]] for l in open(REF) if l.strip() and not l.startswith("#") and l.split()[0][0] in "-0123456789"])
gh = rows[:21]                      # Gamma-H rows: xi = 2*q_direct(1)
def ref(xi):
    return np.interp(xi, gh[:, 4], gh[:, 7])

def load(d):
    state = open(f"{d}/spin_response_state.dat").read().split("\n")
    uj = float(state[[i for i, l in enumerate(state) if l.startswith("# site u_juelich")][0] + 1].split()[1])
    ef = float([l for l in state if l.startswith("fermi_used")][0].split()[1])
    disp = np.loadtxt(f"{d}/spin_response_dispersion.dat", usecols=range(7), ndmin=2)
    xi = np.sqrt((disp[:, 3:6] ** 2).sum(1))
    return uj, ef, xi

def analyse(f, uj, ons_frac=0.01):
    a = np.loadtxt(f)
    w, c0r, c0i, trl = a[:, 0], a[:, 1], a[:, 2], a[:, 5]
    mx = []
    for i in range(1, len(w) - 1):
        if trl[i] > trl[i - 1] and trl[i] >= trl[i + 1]:
            den = trl[i - 1] - 2 * trl[i] + trl[i + 1]
            h = 0.5 * (w[i + 1] - w[i - 1])
            wp = w[i] + h * (trl[i - 1] - trl[i + 1]) / (2 * den)
            mx.append((wp, trl[i]))
    def onset(fr):
        s = -c0i
        idx = np.nonzero(s >= fr * s.max())[0]
        return w[idx[0]] if len(idx) else np.nan
    on1, on5 = onset(0.01), onset(0.05)
    g = 1 + uj * c0r
    cr = []
    for i in range(len(w) - 1):
        if g[i] * g[i + 1] < 0:
            cr.append(w[i] - g[i] * (w[i + 1] - w[i]) / (g[i + 1] - g[i]))
    return mx, on1, on5, cr

def run(sel=None):
    out = []
    for d in sorted(glob.glob(f"{SP}/runs/m*_e*"), key=lambda p: (int(re.search(r"m(\d+)_", p).group(1)), float(p.split("_e")[1]))):
        if not os.path.exists(f"{d}/spin_response_dispersion.dat"): continue
        N = int(re.search(r"m(\d+)_", d).group(1)); eta = float(d.split("_e")[1])
        uj, ef, xis = load(d)
        for iq, xi in enumerate(xis):
            mx, on1, on5, cr = analyse(f"{d}/spin_response_q{iq+1:03d}.dat", uj)
            out.append(dict(N=N, eta=eta, xi=xi, uj=uj, ef=ef, mx=mx, on1=on1, on5=on5, cr=cr, ref=ref(xi)))
    return out

if __name__ == "__main__":
    res = run()
    m = lambda x: f"{x*RY:8.2f}"
    for r in res:
        print(f"N={r['N']:2d} eta={r['eta']:.0e} xi={r['xi']:.4f} ref={r['ref']:.1f}meV U={r['uj']:.6f} EF={r['ef']:.8f} onset1%={m(r['on1'])} onset5%={m(r['on5'])} meV")
        print("   maxima (meV, trL): " + "; ".join(f"{w*RY:.1f}({h:.3g})" for w, h in r['mx']))
        print("   crossings (meV): " + ", ".join(f"{c*RY:.1f}" for c in r['cr']) + f"   [below onset1%: " + ", ".join(f"{c*RY:.1f}" for c in r['cr'] if c < r['on1']) + "]")
