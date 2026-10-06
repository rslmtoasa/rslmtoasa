"""analyze3.py : round-3 scan tables and the xi plot. Peak = highest interior local maximum of tr L; FWHM from half of each maximum's own height."""
import glob, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A
RY, SP = A.RY, A.SP
ROWS = A.gh                                         # Gamma-H rows of fe_lswt.dat: col 4 = xi, col 7 = omega (meV)
SENS = [(0.15, .33), (1/6, .32), (0.25, .20), (1/3, .07), (0.5, .01)]   # header: cutoff sensitivity

def reference(xi):
    k = np.nonzero(abs(ROWS[:, 4] - xi) < 1e-6)[0]
    if len(k): return ROWS[k[0], 7], "read"
    i = np.argsort(abs(ROWS[:, 4] - xi))[:4]
    return float(np.polyval(np.polyfit(ROWS[i, 4], ROWS[i, 7], 3), xi)), "interp"

def sens(xi): return np.interp(xi, *zip(*SENS))

def fwhm(w, y, i):
    half = y[i] / 2; l = r = i
    while l > 0 and y[l] >= half: l -= 1
    while r < len(w) - 1 and y[r] >= half: r += 1
    if y[l] >= half or y[r] >= half: return None, False
    wl = w[l] + (half - y[l]) * (w[l + 1] - w[l]) / (y[l + 1] - y[l])
    wr = w[r - 1] + (half - y[r - 1]) * (w[r] - w[r - 1]) / (y[r] - y[r - 1])
    merged = any(y[j] > y[j - 1] and y[j] >= y[j + 1] for j in range(l + 1, r) if j != i)
    return (wr - wl) * RY, merged

def scan(prefix):
    out = []
    for d in sorted(glob.glob(f"{SP}/runs/{prefix}*")):
        uj, ef, xis = A.load(d)
        st = open(f"{d}/spin_response_state.dat").read().split("\n")
        u_m = float(st[[i for i, l in enumerate(st) if l.startswith("# site u_juelich")][0] + 1].split()[2])
        meth = [l for l in st if l.startswith("method")][0].split()[1]
        U = u_m if meth == "mills" else uj
        for iq, xi in enumerate(xis):
            a = np.loadtxt(f"{d}/spin_response_q{iq+1:03d}.dat"); w, trl = a[:, 0], a[:, 5]
            mx = [(i, w[i] * RY, trl[i]) for i in range(1, len(w) - 1) if trl[i] > trl[i - 1] and trl[i] >= trl[i + 1]]
            # parabolic refinement as in analyze.analyse
            ref_mx = []
            for i, _, h in mx:
                den = trl[i - 1] - 2 * trl[i] + trl[i + 1]; hh = 0.5 * (w[i + 1] - w[i - 1])
                ref_mx.append((i, (w[i] + hh * (trl[i - 1] - trl[i + 1]) / (2 * den)) * RY, h))
            ref_mx.sort(key=lambda m: -m[2])
            fw = {m[0]: fwhm(w, trl, m[0]) for m in ref_mx}
            out.append(dict(xi=xi, U=U, uj=uj, ef=ef, mx=ref_mx, fw=fw, delta=1 + U * a[0, 1], om0=w[0] * RY))
    return sorted(out, key=lambda r: r["xi"])

def show(tag, res):
    print(f"\n== {tag}: U = {res[0]['U']:.8f} Ry (E_F {res[0]['ef']:.8f}); peak (meV), height, FWHM; maxima >= 10% of largest")
    for r in res:
        ref, kind = reference(r["xi"]); top = r["mx"][0]
        s = "; ".join(f"{m[1]:.1f}({m[2]:.0f}," + ("n/a" if r['fw'][m[0]][0] is None else f"{r['fw'][m[0]][0]:.0f}" + ("m" if r['fw'][m[0]][1] else "")) + ")"
                      for m in r["mx"] if m[2] >= 0.1 * top[2])
        print(f"xi={r['xi']:.4f} ref={ref:.1f}({kind}) peak={top[1]:.1f} ({100*(top[1]-ref)/ref:+.1f}%) delta_q(w=0,eta=1e-3)={r['delta']:.4f} | {s}")

if __name__ == "__main__":
    s1, s2 = scan("m60_e0.001_s1"), scan("m60_e0.001_s2")
    show("scan 1 juelich", s1); show("scan 2 coefficient + Goldstone U", s2)
    print("\nsummary xi | ref | s1 peak | s2 peak | s2-s1 (meV) | delta_q s1 | delta_q s2")
    for a, b in zip(s1, s2):
        ref, k = reference(a["xi"])
        print(f"{a['xi']:.4f} | {ref:.1f}{'*' if k=='interp' else ''} | {a['mx'][0][1]:.1f} | {b['mx'][0][1]:.1f} | {b['mx'][0][1]-a['mx'][0][1]:+.1f} | {a['delta']:.4f} | {b['delta']:.4f}")
    # plot
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    xx = np.linspace(0.15, 0.55, 81); inside = (ROWS[:, 4] >= 0.15) & (ROWS[:, 4] <= 0.55)
    cur = np.interp(xx, ROWS[:, 4], ROWS[:, 7]); s = sens(np.clip(xx, 0.15, 0.5))
    ax.fill_between(xx, cur * (1 - s), cur * (1 + s), color="C3", alpha=0.18, label="LSWT cutoff-sensitivity band (only)")
    ax.plot(ROWS[inside, 4], ROWS[inside, 7], "-", color="C3", lw=1.5, label="LSWT (fe_lswt.dat)")
    for res, col, lab, mk in [(s1, "k", "Juelich (item 1)", "o"), (s2, "C0", "coefficient + own Goldstone U (item 2)", "s")]:
        ax.plot([r["xi"] for r in res], [r["mx"][0][1] for r in res], mk + "-", color=col, ms=6, lw=1.3, label=f"dominant peak, {lab}")
        for r in res:
            top = r["mx"][0][2]
            for m in r["mx"][1:]:
                if m[2] >= 0.1 * top:
                    ax.scatter(r["xi"], m[1], s=6 + 80 * m[2] / top, facecolors="none" if col != "k" else col, edgecolors=col, alpha=0.5, lw=0.8)
    ax.set_xlabel("xi (2pi/a, direction (0, 1, 0))"); ax.set_ylabel("omega (meV)"); ax.set_xlim(0.18, 0.52); ax.set_ylim(0, 720)
    ax.set_title("60^3, eta = 1e-3 Ry. Dots: other tr L maxima >= 10% of the dominant, area ~ height")
    ax.legend(fontsize=8, loc="upper left")
    fig.text(0.01, 0.005, "The band shows only the reference's own cutoff sensitivity (header values, linearly interpolated); it is not an error bar on the response.", fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{SP}/xi_scan.png", dpi=130)
