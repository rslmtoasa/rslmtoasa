"""val18_table.py : the VAL-18 results table from ./runs (60^3 Juelich scans at eta 1e-3 and 2e-3) and tolerances.nml.
Run in the work directory that holds runs/. Region I: xi <= 0.25, II: 0.25 < xi < 0.45 (relative peak deviation vs region1_rel,
region2_rel), III: xi >= 0.45 (LSWT energy inside the contiguous interval around the dominant peak where tr L >= half its height)."""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze3 as B
import plots4 as P
RY = B.RY
TOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../spin_response/oracles/tolerances.nml")

def key(name):
    m = re.search(rf"^\s*{name}\s*=\s*([0-9.eE+-]+)", open(TOL).read(), flags=re.M)
    return float(m.group(1)) if m else None

def half_span(w, y):
    i = int(np.argmax(y[1:-1])) + 1
    l = r = i
    while l > 0 and y[l] >= y[i] / 2: l -= 1
    while r < len(w) - 1 and y[r] >= y[i] / 2: r += 1
    lo = w[l] + (y[i] / 2 - y[l]) * (w[l + 1] - w[l]) / (y[l + 1] - y[l]) if y[l] < y[i] / 2 else np.nan
    hi = w[r - 1] + (y[i] / 2 - y[r - 1]) * (w[r] - w[r - 1]) / (y[r] - y[r - 1]) if y[r] < y[i] / 2 else np.nan
    return lo, hi

if __name__ == "__main__":
    t1, t2 = key("region1_rel"), key("region2_rel")
    runs = {"1e-3": ("m60_e0.001_s1", B.scan("m60_e0.001_s1")), "2e-3": ("m60_e0.002_r4", B.scan("m60_e0.002_r4"))}
    cols = {k: P.columns(v[0]) for k, v in runs.items()}
    print("| xi | LSWT (meV) | peak eta 1e-3 (meV) | peak eta 2e-3 (meV) | region | criterion | deviation eta 1e-3 / 2e-3 | tolerance | within tolerance eta 1e-3 / 2e-3 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for i, r1 in enumerate(runs["1e-3"][1]):
        xi = r1["xi"]; ref, kind = B.reference(xi); r2 = runs["2e-3"][1][i]
        p1, p2 = r1["mx"][0][1], r2["mx"][0][1]
        if xi <= 0.25 + 1e-9: reg, crit, tol = "I", "region1_rel", t1
        elif xi < 0.45: reg, crit, tol = "II", "region2_rel", t2
        else: reg, crit, tol = "III", "half-maximum span", None
        if tol is not None:
            dev = [(p - ref) / ref for p in (p1, p2)]
            dtxt = " / ".join(f"{d:+.3f}" for d in dev); ttxt = f"{tol}"
            ok = " / ".join("yes" if abs(d) <= tol else "no" for d in dev)
        else:
            sp = [half_span(cols[k][1], cols[k][2][np.argmin(abs(cols[k][0] - xi))]) for k in ("1e-3", "2e-3")]
            dtxt = " / ".join(f"span {lo:.0f}-{hi:.0f} meV" for lo, hi in sp)
            ttxt = "n/a"; ok = " / ".join("yes" if lo <= ref <= hi else "no" for lo, hi in sp)
        print(f"| {xi:.4f} | {ref:.1f}{'*' if kind=='interp' else ''} | {p1:.1f} | {p2:.1f} | {reg} | {crit} | {dtxt} | {ttxt} | {ok} |")
