"""Shared reader for <prefix>_state.dat and <prefix>_q<NNN>.dat of one run directory."""
import numpy as np

RY_MEV = 13605.693122994


def read_state(d, prefix="spin_response"):
    """Return dict with m_juelich, m_mills, u_juelich, u_mills (single site) and the mesh size."""
    out, rows, mode = {}, [], None
    for line in open(f"{d}/{prefix}_state.dat"):
        t = line.split()
        if line.startswith("# site m_mills"): mode = "m"; continue
        if line.startswith("# site u_juelich"): mode = "u"; continue
        if t and not line.startswith("#") and mode == "m": out["m_mills"], out["m_juelich"] = float(t[1]), float(t[2])
        elif t and not line.startswith("#") and mode == "u": out["u_juelich"], out["u_mills"] = float(t[1]), float(t[2])
        elif t and t[0] == "mesh_nk": out["nk"] = int(t[1])
        elif t and t[0] in ("fermi_used", "kT"): out[t[0]] = float(t[1])
    return out


def read_q(d, iq, prefix="spin_response"):
    """Columns of <prefix>_q<iq>.dat: omega, re/im tr chi0, re/im tr chi, tr L, min|eig| (Ry units)."""
    return np.loadtxt(f"{d}/{prefix}_q{iq:03d}.dat")


def read_dispersion(d, prefix="spin_response"):
    """q direct (3), q cartesian (3, 2pi/a), |q| (1/Angstrom) per q, as an (nq, 7) array."""
    return np.array([[float(x) for x in l.replace("n/a", " n/a").split()[:7]] for l in open(f"{d}/{prefix}_dispersion.dat") if not l.startswith("#")])
