"""plots4.py : round-4 outputs. eta=2e-3 scan tables, jump comparison, xi-scan plot with a third series, intensity maps of tr L(xi, omega)."""
import glob, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A, analyze3 as B
RY, SP = A.RY, A.SP

def columns(prefix):
    """(xi array, omega meV, trL[nxi, nw]) of a juelich scan, sorted by xi."""
    xis, rows = [], []
    for d in sorted(glob.glob(f"{SP}/runs/{prefix}*")):
        _, _, xi = A.load(d)
        for iq, x in enumerate(xi):
            a = np.loadtxt(f"{d}/spin_response_q{iq+1:03d}.dat"); xis.append(x); rows.append(a[:, 5]); w = a[:, 0] * RY
    o = np.argsort(xis)
    return np.array(xis)[o], w, np.array(rows)[o]

if __name__ == "__main__":
    s1 = B.scan("m60_e0.001_s1"); s3 = B.scan("m60_e0.002_r4"); s2 = B.scan("m60_e0.001_s2")
    B.show("scan 3 juelich eta=2e-3", s3)
    print("\nxi | ref | peak eta1e-3 | peak eta2e-3 | step to next xi (1e-3) | step to next xi (2e-3)")
    for i, (a, c) in enumerate(zip(s1, s3)):
        ref, k = B.reference(a["xi"]); pa, pc = a["mx"][0][1], c["mx"][0][1]
        na = f"{s1[i+1]['mx'][0][1]-pa:+.1f}" if i + 1 < len(s1) else "-"; nc = f"{s3[i+1]['mx'][0][1]-pc:+.1f}" if i + 1 < len(s3) else "-"
        print(f"{a['xi']:.4f} | {ref:.1f}{'*' if k=='interp' else ''} | {pa:.1f} | {pc:.1f} | {na} | {nc}")
    # xi-scan plot with the eta = 2e-3 series
    ROWS = B.ROWS
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    xx = np.linspace(0.15, 0.55, 81); inside = (ROWS[:, 4] >= 0.15) & (ROWS[:, 4] <= 0.55)
    cur = np.interp(xx, ROWS[:, 4], ROWS[:, 7]); s = B.sens(np.clip(xx, 0.15, 0.5))
    ax.fill_between(xx, cur * (1 - s), cur * (1 + s), color="C3", alpha=0.18, label="LSWT cutoff-sensitivity band (only)")
    ax.plot(ROWS[inside, 4], ROWS[inside, 7], "-", color="C3", lw=1.5, label="LSWT (fe_lswt.dat)")
    for res, col, lab, mk, fill in [(s1, "k", "Juelich, eta 1e-3", "o", True), (s3, "C2", "Juelich, eta 2e-3", "D", True),
                                    (s2, "C0", "coefficient + own Goldstone U, eta 1e-3", "s", False)]:
        ax.plot([r["xi"] for r in res], [r["mx"][0][1] for r in res], mk + "-", color=col, ms=6, lw=1.3, label=f"dominant peak, {lab}")
        for r in res:
            top = r["mx"][0][2]
            for m in r["mx"][1:]:
                if m[2] >= 0.1 * top:
                    ax.scatter(r["xi"], m[1], s=6 + 80 * m[2] / top, facecolors=col if fill else "none", edgecolors=col, alpha=0.4, lw=0.8)
    ax.set_xlabel("xi (2pi/a, direction (0, 1, 0))"); ax.set_ylabel("omega (meV)"); ax.set_xlim(0.18, 0.52); ax.set_ylim(0, 720)
    ax.set_title("60^3. Dots: other tr L maxima >= 10% of the dominant, area ~ height")
    ax.legend(fontsize=7, loc="upper left")
    fig.text(0.01, 0.005, "The band shows only the reference's own cutoff sensitivity (header values, linearly interpolated); it is not an error bar on the response.", fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{SP}/xi_scan_r4.png", dpi=130); plt.close(fig)
    # intensity maps
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, prefix, eta in [(axs[0], "m60_e0.001_s1", "1e-3"), (axs[1], "m60_e0.002_r4", "2e-3")]:
        xi, w, t = columns(prefix)
        t = t / t.max(axis=1, keepdims=True)
        dx = np.diff(xi).mean(); edges = np.append(xi - dx / 2, xi[-1] + dx / 2)
        wd = np.diff(w).mean(); wedges = np.append(w - wd / 2, w[-1] + wd / 2)
        m = ax.pcolormesh(edges, wedges, t.T, cmap="viridis", vmin=0, vmax=1, shading="flat")
        ax.fill_between(xx, cur * (1 - s), cur * (1 + s), color="w", alpha=0.25)
        ax.plot(ROWS[inside, 4], ROWS[inside, 7], "-", color="w", lw=1.3, label="LSWT")
        ax.plot(xx, cur * (1 - s), "w:", lw=0.8); ax.plot(xx, cur * (1 + s), "w:", lw=0.8, label="LSWT +/- cutoff sensitivity (only)")
        ax.set_xlim(edges[0], edges[-1]); ax.set_ylim(0, 700); ax.set_xlabel("xi (2pi/a)"); ax.set_title(f"60^3, eta = {eta} Ry, Juelich")
        ax.legend(fontsize=7, loc="upper left")
    axs[0].set_ylabel("omega (meV)")
    fig.colorbar(m, ax=axs, label="tr L / max over omega at each xi")
    fig.text(0.01, 0.005, "Each column (one xi) is normalized to its own maximum. The band shows only the reference's own cutoff sensitivity; it is not an error bar on the response.", fontsize=7)
    fig.savefig(f"{SP}/intensity_map.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    print("plots written")
