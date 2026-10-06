#!/usr/bin/env python3
"""Convert the developer's UppASD AMS output into tests/spin_response/oracles/fe_lswt.dat.

Inputs: qfile.kpath (direct q, first line = count) and ams.rs_expor.out
(index, omega, path coordinate).  Assumes omega in meV (see header).
"""
import sys
import numpy as np

ALAT = 2.86120                       # Angstrom, from tests/spin_response/fe_bcc/input.nml
MEV_PER_RY = 13605.693
cell = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]])   # rows a1..a3, units of a
recip = np.linalg.inv(cell).T        # rows b1..b3, units of 2*pi/a

def read_q(path):
    lines = open(path).read().split("\n")
    n = int(lines[0].split()[0])
    q, lab = [], []
    for ln in lines[1:1 + n]:
        p = ln.split()
        q.append([float(x) for x in p[:3]]); lab.append(p[3] if len(p) > 3 else "")
    return np.array(q), lab

def main(qpath, amspath, out):
    qd, lab = read_q(qpath)
    ams = np.loadtxt(amspath)
    assert ams.shape[0] == qd.shape[0], "q count and AMS rows differ"
    w_mev = ams[:, 1]
    qc = qd @ recip                                  # Cartesian, units 2*pi/a
    qabs = np.linalg.norm(qc, axis=1) * 2 * np.pi / ALAT
    with open(out, "w") as fh:
        fh.write("# Fe bcc adiabatic magnon spectrum, LKAG J(r) + LSWT (UppASD AMS). Developer-owned oracle.\n")
        fh.write("# ground state: tests/spin_response/fe_bcc (Fe.nml), level h; J(r): real-space block recursion,\n")
        fh.write("#   lld = 31, 2000 energy points (empirical baseline); M = 2.087188 mu_B (total spin moment)\n")
        fh.write("# J(r): 46 shells, cutoff 5a (1066 neighbours), plain undamped sum; J in mRy; omega in meV.\n")
        fh.write("# convention verified: omega = (4/M)[J(0) - J(q)], H = -sum_{i!=j} J_ij e_i.e_j; recomputed from\n")
        fh.write("#   J(r) it reproduces this AMS at all 91 points to 3e-7 meV.\n")
        fh.write("# cutoff sensitivity along Gamma-H (|omega(4a) - omega(5a)| / omega(5a)): xi 0.05 38%, 0.15 33%,\n")
        fh.write("#   1/6 32%, 0.25 20%, 1/3 7%, 0.5 1%, H 6%.  Small q is NOT converged in the shell cutoff.\n")
        fh.write("# columns: q_direct(1:3) | q_cart(1:3) [2pi/a] | |q| [1/A] | omega [meV] | omega [Ry] | label\n")
        for i in range(qd.shape[0]):
            fh.write(" ".join(f"{x:12.8f}" for x in qd[i]) + "  " +
                     " ".join(f"{x:12.8f}" for x in qc[i]) +
                     f"  {qabs[i]:12.8f}  {w_mev[i]:14.9f}  {w_mev[i]/MEV_PER_RY:16.9e}  {lab[i]}\n")
    return qd, qc, qabs, w_mev

if __name__ == "__main__":
    main(*sys.argv[1:4])
