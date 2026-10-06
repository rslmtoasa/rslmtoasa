"""S3: bcc Fe as a three-atom cell against the one-atom cell, through gauge-invariant quantities.

Usage: s3_supercell.py <rslmto.x> <deck directory> <scratch directory> <tolerances.nml>
The deck directory (tests/spin_response/fe_bcc) is only read. Two scratch decks are run with post_processing = 'spin_response':
  prim: the deck as it is, mesh 4 x 4 x 12;
  sc3 : the same crystal as three atoms per cell, A1 = a1, A2 = a2, A3 = 3 a3 (a_i the bcc primitive vectors), the same Fe.nml
        on every site, mesh 4 x 4 x 4. Its reciprocal basis is b1, b2, b3/3, so its k mesh and the folded set
        {k, k + b3/3, k + 2 b3/3} are exactly the 4 x 4 x 12 mesh of the one-atom cell.
A supercell q with direct coordinates (d1, d2, d3) is the one-atom q = (d1, d2, (d3 + n)/3), n = 0, 1, 2. For identical sites and
a site-independent U, chi0 of the cell is chi0_ij = (1/3) sum_n exp(i G_n.(tau_i - tau_j)) chi0_prim(q + G_n), so
  tr chi0, tr chi = sum_n of the one-atom values, the eigenvalues of I + chi0 U are 1 + U chi0_prim(q + G_n),
and U_Juelich is the one-atom value on every site. Compared, as max|a - b|/max|b| (U per site: relative), against s3_rel:
tr chi0, tr chi, the smallest |eig(I + chi0 U)| (the 7th column of the q files) and U per site.

Known blind spots, none of which this test can see:
  - faults that act identically on both cells (a factor in chi0, the sign of eta, potential-parameter or occupation errors):
    both cells run through the same binary;
  - a transposed vertex (chi0 -> chi0^T leaves traces and eigenvalues unchanged);
  - a wrong site index between identical sites (the sites have identical potentials); this closes with non-identical sites in Stage 3.
With two atoms the Bloch phase exp(+-2 pi i G.tau_i) is +-1, so the sign of the phase would be invisible; three atoms make it complex.
"""
import pathlib
import re
import shutil
import subprocess
import sys

M = 4
Q_SC = [(-0.25, 0.25, 0.75), (0.25, 0.0, 0.25)]   # supercell direct coordinates
TAU_CART = [(0, 0, 0), (0.5, 0.5, -0.5), (1.0, 1.0, -1.0)]   # cartesian, units of alat


def deck(binary, source, scratch, name, nk, extra, qs):
    d = pathlib.Path(scratch, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(source, d)
    text = pathlib.Path(source, "input.nml").read_text()
    old = "pre_processing = 'bravais'"
    if old not in text:
        sys.exit("s3_supercell.py: " + old + " not found in " + source + "/input.nml")
    text = text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'")
    for key, n in zip(("nk1", "nk2", "nk3"), nk):
        text, c = re.subn(rf"^{key} = 12$", f"{key} = {n}", text, flags=re.M)
        if c != 1:
            sys.exit("s3_supercell.py: " + key + " = 12 not found in " + source + "/input.nml")
    if extra:
        for old, new in (("crystal_sym = 'bcc'", "crystal_sym = 'file'"), ("ntype = 1", "ntype = 3"), ("ct(:) = 4.0d0", "ct(:) = 3*4.0d0"),
                         ("label(1) = 'Fe'", "label(1) = 'Fe'\nlabel(2) = 'Fe'\nlabel(3) = 'Fe'")):
            if old not in text:
                sys.exit("s3_supercell.py: " + old + " not found in " + source + "/input.nml")
            text = text.replace(old, new)
        lat = "&lattice\n nbulk_bulk = 3\n ntot = 3\n nbas = 3\n nrec = 3\n a(:, 1) = -0.5, 0.5, 0.5\n a(:, 2) = 0.5, -0.5, 0.5\n a(:, 3) = 1.5, 1.5, -1.5\n"
        lat += "".join(f" crd(:, {i + 1}) = {t[0]}, {t[1]}, {t[2]}\n" for i, t in enumerate(TAU_CART))
        lat += "".join(f" {k}({i}) = {i}\n" for k in ("izp", "no", "iu", "ib", "irec") for i in (1, 2, 3)) + "/\n"
        pathlib.Path(d, "lattice.nml").write_text(lat)
    q = "".join(f"q_list(:, {i + 1}) = {t[0]!r}, {t[1]!r}, {t[2]!r}\n" for i, t in enumerate(qs))
    block = f"\n&spin_response\nmethod = 'juelich'\nn_q = {len(qs)}\n{q}omega_min = 0.0\nomega_max = 0.03\nn_omega = 41\neta = 3.0e-3\n/\n"
    nml = pathlib.Path(d, "input_driver.nml")
    nml.write_text(text + block)
    subprocess.run([binary, nml.name], cwd=d, check=True, stdout=open(pathlib.Path(d, "driver.log"), "w"), stderr=subprocess.STDOUT)
    return d


def qfile(d, i):
    rows = [[float(x) for x in l.split()] for l in pathlib.Path(d, f"spin_response_q{i:03d}.dat").read_text().splitlines() if not l.startswith("#")]
    return ([complex(r[1], r[2]) for r in rows], [complex(r[3], r[4]) for r in rows], [r[6] for r in rows])


def uj(d):
    lines = pathlib.Path(d, "spin_response_state.dat").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# site u_juelich"))
    return [float(l.split()[1]) for l in lines[i + 1:] if l.strip()]


def rel(a, b):
    return max(abs(x - y) for x, y in zip(a, b)) / max(abs(y) for y in b)


def main():
    binary, source, scratch, tolfile = sys.argv[1:5]
    m = re.search(r"^\s*s3_rel\s*=\s*([0-9.eEdD+-]+)", pathlib.Path(tolfile).read_text(), flags=re.M)
    tol = float(m.group(1).replace("d", "e").replace("D", "e")) if m else -1.0
    if not tol > 0:
        sys.exit(f"{tolfile}: key s3_rel missing or <= 0")
    qs_prim = [(d1, d2, (d3 + n) / 3) for d1, d2, d3 in Q_SC for n in (0, 1, 2)]
    prim = deck(binary, source, scratch, "prim", (M, M, 3 * M), False, qs_prim)
    sc3 = deck(binary, source, scratch, "sc3", (M, M, M), True, Q_SC)
    u_prim, u_sc = uj(prim)[0], uj(sc3)
    results = [("U_Juelich per site", max(abs(u - u_prim) for u in u_sc) / u_prim)]
    vacuous = []
    for i in range(len(Q_SC)):
        c0, c, me = qfile(sc3, i + 1)
        parts = [qfile(prim, 3 * i + n + 1) for n in range(3)]
        sum0 = [sum(p[0][w] for p in parts) for w in range(len(c0))]
        sum1 = [sum(p[1][w] for p in parts) for w in range(len(c0))]
        eig = [min(abs(1 + u_prim * p[0][w]) for p in parts) for w in range(len(c0))]
        results += [(f"q{i + 1} tr chi0", rel(c0, sum0)), (f"q{i + 1} tr chi", rel(c, sum1)), (f"q{i + 1} min |eig(I + chi0 U)|", rel(me, eig))]
        diag = [abs(1 + u_prim * z / 3) for z in c0]      # what dropping the off-diagonal elements would give
        vacuous.append((f"q{i + 1} non-vacuity: min |eig| differs from the diagonal-only value by more than 1e-3", rel(me, diag) > 1e-3))
    print("S3 check | measured | tolerance")
    failed = False
    for label, value in results:
        print(f"{label} | {value:.3e} | {tol:.1e} | {'ok' if value <= tol else 'FAIL'}")
        failed = failed or not value <= tol
    for label, ok in vacuous:
        print(f"{label} | {'ok' if ok else 'FAIL'}")
        failed = failed or not ok
    sys.exit(1 if failed else 0)


main()
