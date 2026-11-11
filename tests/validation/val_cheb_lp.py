#!/usr/bin/env python3
"""Stage C production check of linear prediction (LP) of Chebyshev moments.

Frozen potential, post-processing only. Per deck three runs with the same
reconstruction and the same total number of moments for (ii) and (iii):
  (i)   lld = L,              Jackson@N                  (baseline)
  (ii)  lld = L, cheb_lp = T, factor F                   (LP)
  (iii) lld = F*(L+1) - 1,    Jackson                    (computed reference)
(ii) vs (iii) is the holdout test: same reconstruction, they differ only in
the moments n >= N. Decks: example/bulk/bccFe (on-site DOS; the potential is
converged once with the deck's block recursion) and
tests/regression/triad_bccFe_exchange (J_ij, shells 1-3).
Usage: val_cheb_lp.py --binary build/bin/rslmto.x --work <scratch dir> [--L 150] [--F 4]
Prints tables only.
"""
import argparse
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def set_key(text, nml, key, value):
    """Set `key = value` inside namelist &nml (replace the line or append it)."""
    m = re.search(r"(&%s\b.*?\n)(.*?)(\n/)" % nml, text, re.S | re.I)
    body = m.group(2)
    line = "%s = %s" % (key, value)
    pat = re.compile(r"^\s*%s\s*=.*$" % re.escape(key), re.M | re.I)
    body = pat.sub(line, body, count=1) if pat.search(body) else body + "\n" + line
    return text[:m.start(2)] + body + text[m.end(2):]


def run(binary, d):
    """Run the code in d; returns wall time (s), peak RSS (MB), exit code."""
    t0 = time.time()
    with open(d / "run.out", "w") as o, open(d / "run.err", "w") as e:
        p = subprocess.Popen([str(binary)], cwd=d, stdout=o, stderr=e)
        _, status, ru = os.wait4(p.pid, 0)
    return time.time() - t0, ru.ru_maxrss / 1024.0, os.waitstatus_to_exitcode(status)


def log_of(d):
    return ANSI.sub("", (d / "run.out").read_text(errors="replace"))


def variants(L, F):
    return [("i", dict(lld=L, lp=False)), ("ii", dict(lld=L, lp=True)),
            ("iii", dict(lld=F * (L + 1) - 1, lp=False))]


def prepare(src, dst, v, F, extra):
    """Copy src to dst and patch input.nml for one variant."""
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    p = dst / "input.nml"
    t = p.read_text()
    t = set_key(t, "control", "recur", "'chebyshev'")
    t = set_key(t, "control", "lld", v["lld"])
    if v["lp"]:
        t = set_key(t, "control", "cheb_lp", ".true.")
        t = set_key(t, "control", "cheb_lp_factor", F)
    for nml, key, val in extra:
        t = set_key(t, nml, key, val)
    p.write_text(t)


def execute(binary, dst):
    wall, rss, code = run(binary, dst)
    lg = log_of(dst)
    info = dict(wall=wall, rss=rss, code=code)
    f = re.findall(r"Free Fermi energy:\s*(\S+)", lg)
    info["ef"] = float(f[-1]) if f else float("nan")
    mem = re.findall(r"mu_n\+mu_ng = (\S+) MB", lg)
    info["mem"] = mem[-1] if mem else "n/a"
    info["pos"] = sum(int(x) for x in re.findall(r"local_blocks=(\d+)", lg))
    info["lp_lines"] = len(re.findall(r"Chebyshev LP ", lg))
    lp = re.findall(r"Chebyshev LP .*?: (p=.*)", lg)
    info["lp_first"] = lp[0] if lp else "n/a"
    return info


def top_peaks(e, y, n=2, sep=0.05):
    idx = [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1]]
    out = []
    for i in sorted(idx, key=lambda k: -y[k]):
        if all(abs(e[i] - e[j]) > sep for j in out):
            out.append(i)
        if len(out) == n:
            break
    return [(y[i], e[i]) for i in out]


def bulk(binary, work, L, F):
    print("\n=== bcc Fe bulk (example/bulk/bccFe): on-site DOS, frozen potential ===")
    scf = work / "bulk_scf"
    if not (scf / "Fe_out.nml").exists():
        if scf.exists():
            shutil.rmtree(scf)
        shutil.copytree(ROOT / "example/bulk/bccFe", scf)
        w, _, code = run(binary, scf)
        print("converged potential from the deck's block SCF: exit=%d, %.0f s" % (code, w))
    nsites = re.search(r"II\s*=\s*(\d+)", (scf / "clust").read_text()).group(1)
    print("cluster: %s sites; L=%d N=%d, F=%d N_ext=%d; Chebyshev window -1.2..1.2 Ry, channels_ldos=2500" %
          (nsites, L, 2 * L + 2, F, F * (2 * L + 2)))
    res = {}
    for name, v in variants(L, F):
        dst = work / ("bulk_" + name)
        prepare(scf, dst, v, F, [("self", "nstep", 1), ("energy", "energy_min", -1.2),
                                 ("energy", "energy_max", 1.2), ("energy", "channels_ldos", 2500)])
        info = execute(binary, dst)
        info["dos"] = np.loadtxt(dst / "Fe_dos.out", comments="#")
        res[name] = info
    ref = res["iii"]
    e = ref["dos"][:, 0]
    ef = ref["ef"]
    print("\nrun   lld   LP  exit  E_F(Ry)    wall(s) maxrss(MB) mu(MB) posproj  LP log (first atom)")
    for name, v in variants(L, F):
        r = res[name]
        print("%-4s %4d  %-3s %4d  %9.6f %7.1f %9.0f %7s %6d   %s" % (name, v["lld"], "T" if v["lp"] else "F", r["code"], r["ef"],
              r["wall"], r["rss"], r["mem"], r["pos"], r["lp_first"]))
    print("\nDOS differences vs (iii), states/Ry; occupied band = E<=E_F(iii)=%.6f, window = |E-E_F(iii)|<=0.1" % ef)
    print("run   spin     occ_max    occ_rms   ef_max     ef_rms")
    for name in ("i", "ii"):
        for col, spin in ((2, "up"), (3, "dw")):
            d = res[name]["dos"][:, col] - ref["dos"][:, col]
            occ, win = e <= ef, np.abs(e - ef) <= 0.1
            print("%-4s  %-4s  %9.3e  %9.3e  %9.3e  %9.3e" % (name, spin, np.abs(d[occ]).max(), np.sqrt(np.mean(d[occ] ** 2)),
                  np.abs(d[win]).max(), np.sqrt(np.mean(d[win] ** 2))))
    print("\nrun   DOS_total(E_F own)  DOS_total(E_F(iii))   two highest peaks of total DOS: height@E(Ry)")
    for name, _ in variants(L, F):
        r = res[name]
        tot = r["dos"][:, 1]
        pk = top_peaks(e, tot)
        print("%-4s  %12.4f  %14.4f        %s" % (name, np.interp(r["ef"], e, tot), np.interp(ef, e, tot),
              "  ".join("%.3f@%.3f" % p for p in pk)))


def triad(binary, work, L, F):
    print("\n=== triad bcc Fe exchange (tests/regression/triad_bccFe_exchange): J_ij shells 1-3 ===")
    src = ROOT / "tests/regression/triad_bccFe_exchange"
    pairs = [(335, "shell 1"), (336, "shell 2"), (235, "shell 3")]
    extra = [("lattice", "njij", len(pairs)), ("energy", "channels_ldos", 2500)] + \
            [("lattice", "ijpair(%d, :)" % (k + 1), "1, %d" % p[0]) for k, p in enumerate(pairs)]
    res = {}
    for name, v in variants(L, F):
        dst = work / ("triad_" + name)
        prepare(src, dst, v, F, extra)
        info = execute(binary, dst)
        info["jij"] = np.loadtxt(dst / "jij.out")[:, 5]
        res[name] = info
    print("run   lld   LP  exit   wall(s) maxrss(MB) mu(MB) posproj  " + "  ".join("J(%s)" % p[1] for p in pairs))
    for name, v in variants(L, F):
        r = res[name]
        print("%-4s %4d  %-3s %4d  %7.1f %9.0f %7s %6d   %s" % (name, v["lld"], "T" if v["lp"] else "F", r["code"], r["wall"], r["rss"],
              r["mem"], r["pos"], "  ".join("%9.5f" % x for x in r["jij"])))
    print("(J column 6 of jij.out, in the units the code writes)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--L", type=int, default=150)
    ap.add_argument("--F", type=int, default=4)
    ap.add_argument("--deck", choices=["bulk", "triad", "all"], default="all")
    a = ap.parse_args()
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    binary = Path(a.binary).resolve()
    if a.deck in ("bulk", "all"):
        bulk(binary, work, a.L, a.F)
    if a.deck in ("triad", "all"):
        triad(binary, work, a.L, a.F)


if __name__ == "__main__":
    main()
