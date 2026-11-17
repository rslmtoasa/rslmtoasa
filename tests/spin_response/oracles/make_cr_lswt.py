#!/usr/bin/env python3
"""Record the Cr adiabatic magnon oracle: LKAG J(r) + two-sublattice linear spin waves.

Developer-owned, next to make_fe_lswt.py. Three subcommands:

  qpath     write the q path of the simple cubic two-atom cell (UppASD qfile
            format: count, then q_direct and label per line), Gamma-X-M-Gamma-R-X|M-R,
            in steps of 1/n (direct = Cartesian in units of 2 pi / a for this cell).
            Use a response mesh that is a multiple of n.
  lswt      compute omega(q) from a J(r) shell file and write cr_lswt.dat. If an
            UppASD AMS file is given, compare with it and record the deviation.
  selftest  check this script: (1) a J1-J2 CsCl antiferromagnet against the
            closed form; (2) optionally the committed fe_lswt.dat from its J(r).

Convention (the one the Fe oracle verified): H = -sum_{i != j} J_ij e_i . e_j,
J in mRy, ferromagnet omega = (4/M) [J(0) - J(q)]. For collinear moments
e_i = s_i z (s = +1 or -1), linearizing the precession gives
    omega phi_a = sum_b D_ab phi_b,
    D_ab = (4 / M_a) [ delta_ab sum_c J_ac(0) s_c - s_a J_ab(q) ],
    J_ab(q) = sum_{j in b} J_ij exp(i q . (r_j - r_i)).
For the antiferromagnet the two eigenvalues are +w and -w, the SAME signed
convention as the response's adiabatic dispersion (the minority sublattice
precesses at negative omega in the global +- channel).

J(r) file: one line per symmetry shell, as lswt_inputs/jfile:
    itype jtype  rx ry rz  J[mRy]  |r|      (r Cartesian, units of a)
The J of the file are those of a ferromagnetic reference state (the UppASD jfile
convention): J_ab refers to parallel moments. The lswt subcommand converts them to
the global spin frame, J_ab -> s_a s_b J_ab, so that J_12 and J_21 change sign for the
antiferromagnet; the selftest models are already given in the global frame.
Each shell is expanded by the 48 operations of the cubic group, valid because
every site of the CsCl structure has full cubic site symmetry. The number of
neighbours after expansion is printed; check it against your exchange run.

Usage:
  python3 make_cr_lswt.py qpath --n 32 --out qfile.kpath
  python3 make_cr_lswt.py lswt --jfile jfile --qfile qfile.kpath --moment 0.923 \
        --alat 2.95 --out cr_lswt.dat [--ams ams.rs_expor.out] [--note "provenance line"] ...
  python3 make_cr_lswt.py selftest [--fe-jfile lswt_inputs/jfile --fe-lswt fe_lswt.dat]
"""
import argparse
import itertools
import sys

import numpy as np

MEV_PER_MRY = 13.605693
MEV_PER_RY = 13605.693
SIGNS = (+1, -1)                       # sublattice 1 up, sublattice 2 down (global frame)
CORNERS = {"GAMMA": (0, 0, 0), "X": (0, 0.5, 0), "M": (0.5, 0.5, 0), "R": (0.5, 0.5, 0.5)}
PATH = [("GAMMA", "X"), ("X", "M"), ("M", "GAMMA"), ("GAMMA", "R"), ("R", "X"), None, ("M", "R")]


# ----------------------------------------------------------------------------
# q path
# ----------------------------------------------------------------------------
def build_qpath(n):
    """Points in steps of 1/n along the path; a break (None) starts a new segment."""
    pts, labels = [], []
    start_new = True
    for seg in PATH:
        if seg is None:
            start_new = True
            continue
        a, b = (np.array(CORNERS[s], float) for s in seg)
        steps = int(round(np.max(np.abs(b - a)) * n))
        for i in range(0 if start_new else 1, steps + 1):
            q = a + (b - a) * i / steps
            if not np.allclose(q * n, np.round(q * n), atol=1e-12):
                raise SystemExit(f"q {q} not on the 1/{n} grid")
            lab = seg[0] if i == 0 else (seg[1] if i == steps else "")
            pts.append(q)
            labels.append(lab)
        start_new = False
    return np.array(pts), labels


def write_qfile(path, q, labels):
    with open(path, "w") as fh:
        fh.write(f"{len(q):10d}\n")
        for qq, lab in zip(q, labels):
            fh.write(" ".join(f"{x:11.8f}" for x in qq) + f"    {lab}\n")


def read_qfile(path):
    lines = open(path).read().split("\n")
    n = int(lines[0].split()[0])
    q, lab = [], []
    for ln in lines[1:1 + n]:
        p = ln.split()
        q.append([float(x) for x in p[:3]])
        lab.append(p[3] if len(p) > 3 else "")
    return np.array(q), lab


# ----------------------------------------------------------------------------
# J(r) and the spin-wave matrix
# ----------------------------------------------------------------------------
def cubic_ops():
    ops = []
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            m = np.zeros((3, 3))
            for i in range(3):
                m[i, perm[i]] = sg[i]
            ops.append(m)
    return ops


def expand_shells(rows, expand=True):
    """rows: (itype, jtype, r(3), J). Returns {(a, b): [(r, J), ...]} with the full stars."""
    ops = cubic_ops()
    pairs = {}
    for it, jt, r, j in rows:
        stars = {tuple(np.round(op @ r, 10)) for op in ops} if expand else {tuple(np.round(r, 10))}
        lst = pairs.setdefault((it, jt), [])
        seen = {tuple(np.round(v, 10)) for v, _ in lst}
        for v in stars:
            if v in seen:
                raise SystemExit(f"shell {r} of pair {(it, jt)} overlaps another line: "
                                 "the file lists more than one line per shell (use --no-expand?)")
            lst.append((np.array(v), j))
    return pairs


def read_jfile(path):
    rows = []
    for ln in open(path):
        p = ln.split()
        if len(p) < 6 or p[0].startswith("#"):
            continue
        rows.append((int(p[0]), int(p[1]), np.array([float(x) for x in p[2:5]]), float(p[5])))
    return rows


def jq(pairs, a, b, q_cart):
    """J_ab(q), q Cartesian in units of 2 pi / a, r in units of a."""
    lst = pairs.get((a, b), [])
    if not lst:
        return 0.0
    r = np.array([v for v, _ in lst])
    j = np.array([x for _, x in lst])
    return np.sum(j * np.exp(2j * np.pi * (r @ q_cart)))


def omega_signed(pairs, q_cart, moments, signs):
    """Eigenvalues of D(q) in mRy, sorted ascending; also max |Im| / max |Re|."""
    n = len(moments)
    j0 = np.array([[jq(pairs, a + 1, b + 1, np.zeros(3)) for b in range(n)] for a in range(n)])
    jqm = np.array([[jq(pairs, a + 1, b + 1, q_cart) for b in range(n)] for a in range(n)])
    d = np.zeros((n, n), dtype=complex)
    for a in range(n):
        for b in range(n):
            d[a, b] = -signs[a] * jqm[a, b]
        d[a, a] += np.sum(j0[a, :] * np.array(signs))
        d[a, :] *= 4.0 / moments[a]
    ev = np.linalg.eigvals(d)
    imrel = np.max(np.abs(ev.imag)) / max(np.max(np.abs(ev.real)), 1e-300)
    return np.sort(ev.real), imrel


# ----------------------------------------------------------------------------
# subcommands
# ----------------------------------------------------------------------------
def cmd_qpath(args):
    q, lab = build_qpath(args.n)
    write_qfile(args.out, q, lab)
    print(f"{len(q)} q points in steps of 1/{args.n} written to {args.out}; "
          f"use a response mesh that is a multiple of {args.n}")


def cmd_lswt(args):
    rows = [(a, b, r, j * SIGNS[a - 1] * SIGNS[b - 1]) for a, b, r, j in read_jfile(args.jfile)]
    pairs = expand_shells(rows, expand=not args.no_expand)
    for key in sorted(pairs):
        print(f"  pair {key}: {len(pairs[key])} neighbours, sum J = {sum(x for _, x in pairs[key]):+.6f} mRy")
    q, lab = read_qfile(args.qfile)
    moments = [args.moment, args.moment] if args.moment2 is None else [args.moment, args.moment2]
    w = np.zeros((len(q), 2))
    worst_im = 0.0
    for i, qq in enumerate(q):
        ev, imrel = omega_signed(pairs, qq, moments, SIGNS)
        w[i] = ev
        if np.max(np.abs(ev)) > 1e-9:
            worst_im = max(worst_im, imrel)
    print(f"  max |Im omega| / |Re omega| over the path: {worst_im:.2e} (stable in LSWT if ~0)")
    notes = list(args.note or [])
    if args.ams:
        ams = np.loadtxt(args.ams)
        if ams.shape[0] != len(q):
            raise SystemExit("q count and AMS rows differ")
        ams_mev = np.sort(np.abs(ams[:, 1:1 + 2]), axis=1) if ams.shape[1] > 3 else np.abs(ams[:, 1:2])
        mine = np.sort(np.abs(w), axis=1) * MEV_PER_MRY
        if ams_mev.shape[1] == 1:
            dev = np.max(np.abs(mine[:, 1] - ams_mev[:, 0]))
        else:
            dev = np.max(np.abs(mine - ams_mev))
        print(f"  this LSWT vs the UppASD AMS: max |delta omega| = {dev:.3e} meV")
        notes.append(f"recomputed from J(r) with make_cr_lswt.py, it reproduces the UppASD AMS to {dev:.1e} meV")
    qabs = np.linalg.norm(q, axis=1) * 2 * np.pi / args.alat
    with open(args.out, "w") as fh:
        fh.write("# Cr bcc (CsCl two-atom cell) antiferromagnet, adiabatic magnons from LKAG J(r) + two-sublattice LSWT.\n")
        fh.write("# Developer-owned oracle, made by make_cr_lswt.py.\n")
        fh.write(f"# alat = {args.alat} A; |M| = {moments[0]} / {moments[1]} mu_B (sublattice 1 up, 2 down).\n")
        fh.write("# convention: H = -sum_{i!=j} J_ij e_i.e_j, J in mRy; D_ab = (4/M_a)[delta_ab sum_c J_ac(0) s_c - s_a J_ab(q)];\n")
        fh.write("#   eigenvalues +w and -w, signed as in the response's adiabatic output.\n")
        fh.write("# the J(r) file is in the ferromagnetic-reference (UppASD jfile) convention; converted to the global frame J_ab -> s_a s_b J_ab.\n")
        for k in sorted(pairs):
            fh.write(f"# J(r) pair {k}: {len(pairs[k])} neighbours (shells expanded by the cubic group).\n")
        for line in notes:
            fh.write(f"# {line}\n")
        fh.write("# columns: q_direct(1:3) = q_cart [2pi/a] | |q| [1/A] | omega_minus omega_plus [meV] | same [Ry] | label\n")
        for i in range(len(q)):
            fh.write(" ".join(f"{x:12.8f}" for x in q[i]) + f"  {qabs[i]:12.8f}  " +
                     " ".join(f"{x * MEV_PER_MRY:14.9f}" for x in w[i]) + "  " +
                     " ".join(f"{x * MEV_PER_MRY / MEV_PER_RY:16.9e}" for x in w[i]) + f"  {lab[i]}\n")
    print(f"written to {args.out}")


def cmd_selftest(args):
    ok = True
    # (1) J1-J2 CsCl antiferromagnet against the closed form
    j1, j2, m = -0.35, 0.05, 0.9
    rows = [(1, 2, np.array([0.5, 0.5, 0.5]), j1), (2, 1, np.array([0.5, 0.5, 0.5]), j1),
            (1, 1, np.array([1.0, 0.0, 0.0]), j2), (2, 2, np.array([1.0, 0.0, 0.0]), j2)]
    pairs = expand_shells(rows)
    counts = {k: len(v) for k, v in pairs.items()}
    print(f"  J1-J2 model neighbour counts: {counts} (expected 8 for A-B, 6 for A-A)")
    ok &= counts == {(1, 2): 8, (2, 1): 8, (1, 1): 6, (2, 2): 6}
    q, _ = build_qpath(32)
    worst = 0.0
    for qq in q:
        ev, _ = omega_signed(pairs, qq, [m, m], SIGNS)
        jab = 8 * j1 * np.prod(np.cos(np.pi * qq))
        jaa = 2 * j2 * np.sum(np.cos(2 * np.pi * qq))
        jab0, jaa0 = 8 * j1, 6 * j2
        wexact = 4 / m * np.sqrt(max((jaa0 - jaa - jab0) ** 2 - jab ** 2, 0.0))
        worst = max(worst, abs(ev[1] - wexact), abs(ev[0] + wexact))
    print(f"  J1-J2 antiferromagnet vs closed form, max |delta| = {worst:.2e} mRy")
    ok &= worst < 1e-12
    qs = [np.array([0, 0, x]) for x in (1 / 64, 2 / 64, 4 / 64)]
    ws = [omega_signed(pairs, qq, [m, m], SIGNS)[0][1] for qq in qs]
    p = np.log2(np.array(ws[1:]) / np.array(ws[:-1]))
    print(f"  small-q exponents along Gamma-X: {p[0]:.4f}, {p[1]:.4f} (linear expected)")
    ok &= bool(np.all(np.abs(p - 1) < 0.05))
    # (2) optionally the committed Fe oracle, from its own J(r) (single sublattice)
    if args.fe_jfile and args.fe_lswt:
        pairs_fe = expand_shells(read_jfile(args.fe_jfile))
        nfe = len(pairs_fe[(1, 1)])
        data = np.loadtxt(args.fe_lswt, usecols=(3, 4, 5, 7))
        mfe = 2.087188
        dev = max(abs(omega_signed(pairs_fe, row[:3], [mfe], (+1,))[0][0] * MEV_PER_MRY - row[3])
                  for row in data)
        print(f"  Fe: {nfe} neighbours (header says 1066); fe_lswt.dat reproduced to {dev:.2e} meV")
        ok &= nfe == 1066 and dev < 1e-5
    print("selftest", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("qpath"); a.add_argument("--n", type=int, default=32); a.add_argument("--out", default="qfile.kpath")
    b = sub.add_parser("lswt")
    b.add_argument("--jfile", required=True); b.add_argument("--qfile", required=True)
    b.add_argument("--moment", type=float, required=True, help="|M| of sublattice 1 (mu_B)")
    b.add_argument("--moment2", type=float, help="|M| of sublattice 2 if different")
    b.add_argument("--alat", type=float, required=True); b.add_argument("--out", default="cr_lswt.dat")
    b.add_argument("--ams", help="UppASD AMS output to compare with (index, omega[meV], ...)")
    b.add_argument("--no-expand", action="store_true", help="the J(r) file already lists every neighbour")
    b.add_argument("--note", action="append", help="provenance line for the header (repeatable)")
    c = sub.add_parser("selftest"); c.add_argument("--fe-jfile"); c.add_argument("--fe-lswt")
    args = ap.parse_args()
    return {"qpath": cmd_qpath, "lswt": cmd_lswt, "selftest": cmd_selftest}[args.cmd](args) or 0


if __name__ == "__main__":
    sys.exit(main())
