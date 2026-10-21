"""Shared helpers for the AppGbtKspace* tests: scratch decks from tests/spin_response, the production binary and
the fermi_surface eigenvalue export. Decks are copied; nothing in the repository is modified."""
import re
import shutil
import subprocess
from pathlib import Path

import f90nml
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SR = ROOT / "tests" / "spin_response"
FIX = ROOT / "tests" / "gbt" / "fixtures"


def tolerance(key):
    return f90nml.read(ROOT / "tests/gbt/oracles/tolerances.nml")["gbt_tolerances"][key]


def lattice_vectors(kind):
    return np.array(f90nml.read(FIX / kind / "lattice.nml")["lattice"]["a"])      # rows are the lattice vectors


def make_deck(dst, kind, rep, nk, q_dir=(0.0, 0.0, 0.0), order="first", spin_resolved=False, theta_deg=0.0, ccor=False):
    """kind: bccFe_prim, bccFe_sc2 (Fe, same potential on every site) or cr_afm (cr_bcc, Cr2 = Cr1 with mom = -z).
    q_dir is in direct coordinates of the cell; a gbt_single_q deck gets it as Cartesian q_ss with the cone angle theta_deg."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "cr_bcc", dst)
    nml = f90nml.read(str(SR / "cr_bcc" / "input.nml"))
    if kind == "cr_afm":
        text = (SR / "cr_bcc" / "Cr1.nml").read_text().replace("symbol = 'Cr1'", "symbol = 'Cr2'")
        (dst / "Cr2.nml").write_text(re.sub(r"^ mom = .*$", " mom = 0.0, 0.0, -1.0", text, flags=re.M))
    else:
        for f in ("Cr1.nml", "Cr2.nml"):
            (dst / f).unlink()
        shutil.copy(SR / "fe_bcc" / "Fe.nml", dst / "Fe.nml")
        shutil.copy(FIX / kind / "lattice.nml", dst / "lattice.nml")
        nsite = 1 if kind == "bccFe_prim" else 2
        nml["lattice"].update(alat=2.8612, wav=1.4088, ntype=nsite, ct=[4.0] * nsite)
        nml["atoms"] = {"database": "./", "label": ["Fe"] * nsite}
    return finish_deck(dst, nml, rep, nk, order, spin_resolved, q_dir=q_dir, kind=kind, theta_deg=theta_deg, ccor=ccor)


def make_nio_deck(dst, nid, rep, nk, q_dir=(0.0, 0.0, 0.0), order="second", spin_resolved=False, theta_deg=0.0, nstep=0, vmt_mode=None, nid_file=None):
    """AFM-II NiO, four sites (Niu, Nid, O1, O2), hoh + ccor_2c (strict) as the deck of the fixture nio_af. nid is 'mirrored' (the
    deck's Nid: mom = +z, spin 1 holds the minority values) or 'b1' (Nid_B1: mom = -z, spin 1 holds the majority values).
    nstep = 0 probes the stored potentials unchanged. vmt_mode sets ccor_vmt_mode, nid_file replaces Nid.nml."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    for f in ("lattice.nml", "Niu.nml", "Nid.nml", "O1.nml", "O2.nml"):
        shutil.copy(FIX / "nio_af" / f, dst / f)
    if nid == "b1":
        shutil.copy(nid_file or FIX / "nio_af" / "Nid_B1.nml", dst / "Nid.nml")
    nml = f90nml.read(str(FIX / "nio_af" / "input.nml"))
    nml["self"]["nstep"] = nstep
    return finish_deck(dst, nml, rep, nk, order, spin_resolved, q_dir=q_dir, kind="nio_af", theta_deg=theta_deg, ccor=True, vmt_mode=vmt_mode)


def make_nio_scf_deck(dst, nid, gbt, nk, q_dir=(0.0, 0.0, 0.0), order="second", nstep=250):
    """Self-consistent run of the nio_af deck (hoh + ccor_2c strict): nsp = 1 collinear (nid 'mirrored' only) or, with gbt, an
    unconstrained nsp = 3 gbt_single_q at theta = 0 with q_dir (direct coordinates). txc = 1 (native Barth-Hedin) replaces the
    fixture's libXC selector, which the CI configuration does not have."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    for f in ("lattice.nml", "Niu.nml", "Nid.nml", "O1.nml", "O2.nml"):
        shutil.copy(FIX / "nio_af" / f, dst / f)
    if nid == "b1":
        shutil.copy(FIX / "nio_af" / "Nid_B1.nml", dst / "Nid.nml")
    nml = f90nml.read(str(FIX / "nio_af" / "input.nml"))
    nml["self"]["nstep"] = nstep
    nml["control"].update(txc=1, nsp=3 if gbt else 1)
    nml["reciprocal"].update(nk1=nk[0], nk2=nk[1], nk3=nk[2], kspace_ham_order=order)
    nml["hamiltonian"] = {"hoh": order == "second", "ccor_2c": True, "ccor_strict": True}
    if gbt:
        q_cart = np.linalg.inv(lattice_vectors("nio_af")) @ np.array(q_dir)
        nml["hamiltonian"].update(magnetic_representation="gbt_single_q", q_ss=list(q_cart), theta_ss=0.0)
    nml.write(str(dst / "input.nml"), force=True)
    return dst


def finish_deck(dst, nml, rep, nk, order, spin_resolved, q_dir=(0.0, 0.0, 0.0), kind="cr_afm", theta_deg=0.0, ccor=False, vmt_mode=None):
    nml["calculation"].update(pre_processing="none", post_processing="fermi_surface")
    nml["control"].update(nsp=1 if spin_resolved else 3, recur="block")     # the spin-resolved export needs nsp = 1
    nml["hamiltonian"] = {"hoh": order == "second", "magnetic_representation": rep}
    if rep == "gbt_single_q":
        q_cart = np.linalg.inv(lattice_vectors(kind)) @ np.array(q_dir) if kind != "cr_afm" else np.array(q_dir)
        nml["hamiltonian"].update(q_ss=list(q_cart), theta_ss=theta_deg)
    nml["lattice"]["strux_want_sdot"] = False
    if ccor:       # the two-centre CCOR needs Sdot and a screened structure; strict turns a silent Hcc = 0 into a stop
        nml["lattice"].update(strux_want_sdot=True, screening="fitted")
        nml["hamiltonian"].update(ccor_2c=True, ccor_strict=True)
        if vmt_mode:
            nml["hamiltonian"]["ccor_vmt_mode"] = vmt_mode
    nml["reciprocal"].update(nk1=nk[0], nk2=nk[1], nk3=nk[2], use_symmetry_reduction=False, use_time_reversal=False,
                             auto_find_fermi=True, kspace_ham_order=order, write_spin_resolved_eigenpairs=spin_resolved)
    nml.write(str(dst / "input.nml"), force=True)
    return dst


def run(binary, deck):
    subprocess.run([str(binary), "input.nml"], cwd=deck, check=True, stdout=open(deck / "run.log", "w"),
                   stderr=subprocess.STDOUT)


def read_export(deck, spin=False):
    """(k_points (nk,3), eigenvalues (nbands,nk)) or, for spin=True, spin sectors (2,nbasis,nk)."""
    meta = (deck / "fermi_surface.meta").read_text()
    nk, nb = (int(re.search(rf"n_{n} = (\d+)", meta).group(1)) for n in ("kpoints", "bands"))
    raw = np.fromfile(deck / "fermi_surface.bin", dtype="<f8")
    k = raw[:3 * nk].reshape(3, nk, order="F").T
    if not spin:
        return k, raw[4 * nk:4 * nk + nb * nk].reshape(nb, nk, order="F")
    nbasis = nb // 2
    sp = np.fromfile(deck / "fermi_surface.spin.bin", dtype="<f8", count=2 * nbasis * nk)
    return k, sp.reshape(2, nbasis, nk, order="F")


def mesh_index(k, nk):
    """Map direct k (nk points, folded or not) to mesh integer coordinates."""
    m = np.rint(k * np.array(nk)).astype(int) % np.array(nk)
    return {tuple(row): i for i, row in enumerate(m)}


def make_super_deck(dst, fixture, positions_z, q_z, theta_deg, order, nk, gbt=False, ccor=False):
    """periodic_nc supercell of Fe: one type per site, every type the Fe potential with its moment rotated to
    m(R + tau) = (sin th cos(2 pi q.r), sin th sin(2 pi q.r), cos th), q = (0, 0, q_z) in 2 pi/alat. With gbt=True
    the deck is a q = 0, theta = 0 gbt_single_q one (every mom = +z), used only for its bond dump."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "cr_bcc", dst)
    for f in ("Cr1.nml", "Cr2.nml"):
        (dst / f).unlink()
    shutil.copy(FIX / fixture / "lattice.nml", dst / "lattice.nml")
    th = np.radians(0.0 if gbt else theta_deg)
    for i, z in enumerate(positions_z, start=1):
        phi = 0.0 if gbt else 2.0 * np.pi * q_z * z
        mom = ", ".join(repr(float(x)) for x in (np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), np.cos(th)))
        text = (SR / "fe_bcc" / "Fe.nml").read_text().replace("symbol = 'Fe'", f"symbol = 'Fe{i}'")
        (dst / f"Fe{i}.nml").write_text(re.sub(r"^ mom = .*$", f" mom = {mom}", text, flags=re.M))
    nml = f90nml.read(str(SR / "cr_bcc" / "input.nml"))
    n = len(positions_z)
    n3 = -(-20 // round(positions_z[-1] + 0.5))     # cells along z so that the cluster still reaches 20 alat
    nml["lattice"].update(alat=2.8612, wav=1.4088, ntype=n, ct=[4.0] * n, n3=n3, ndim=20 * 20 * n3 * n + 1000)
    nml["atoms"] = {"database": "./", "label": [f"Fe{i}" for i in range(1, n + 1)]}
    return finish_deck(dst, nml, "gbt_single_q" if gbt else "periodic_nc", nk, order, False, ccor=ccor)


def read_bonds(deck):
    """Raw S blocks from the gbt_bonds.out dump of a gbt_single_q run: [(atom type, d Cartesian, S(i, j))]."""
    bonds = []
    for line in (deck / "gbt_bonds.out").read_text().splitlines():
        f = line.split()
        if f and f[0] == "bond":
            bonds.append((int(f[1]), np.array([float(x) for x in f[-3:]]), np.zeros((9, 9), complex)))
        elif f and f[0] == "S":
            bonds[-1][2][int(f[1]) - 1, int(f[2]) - 1] = complex(float(f[3]), float(f[4]))
    return bonds


def run_orthogonal(binary, dst, nk, temperature, theta_deg, qs):
    """frozen_magnon (mft) with kspace_ham_order = 'orthogonal' on the fixed tests/spin_response/fe_bcc potential, zero SCF
    steps, Gamma-centred nk^3 mesh. Returns (mu, M, DeltaOmega/sin^2(theta) [mRy], omega [meV]) with DeltaOmega = Omega(q, theta) - Omega(0, 0)."""
    dst = Path(dst)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(SR / "fe_bcc", dst)
    nml = f90nml.read(str(SR / "fe_bcc" / "input.nml"))
    nml["calculation"].update(post_processing="frozen_magnon")
    nml["self"].update(nstep=0)
    nml["control"].update(nsp=3)
    nml["reciprocal"].update(nk1=nk, nk2=nk, nk3=nk, temperature=temperature, auto_find_fermi=True, kspace_ham_order="orthogonal")
    nml["hamiltonian"] = {"hoh": False, "theta_ss": theta_deg}
    nml["frozen_magnon"] = {"mode": "mft", "q_file": "q.dat", "q_coordinates": "direct"}
    nml.write(str(dst / "input.nml"), force=True)
    (dst / "q.dat").write_text(f"{len(qs) + 1}\n0 0 0\n" + "\n".join(" ".join(map(str, q)) for q in qs) + "\n")
    run(binary, dst)
    header = (dst / "frozen_magnon_kspace.dat").read_text().splitlines()[0]
    mu, moment = (float(re.search(rf"{k}\s*=\s*([-+0-9.eE]+)", header).group(1)) for k in (r"mu \[Ry\]", r"M \[muB\]"))
    row = np.loadtxt(dst / "frozen_magnon_kspace.dat", ndmin=2)     # columns: see the file header; 7 Omega(q,theta), 8 Omega(q,0), 11 omega [meV]
    dom = 1.0e3 * (row[1:, 7] - row[0, 8]) / np.sin(np.radians(theta_deg)) ** 2
    return mu, moment, dom, row[1:, 11]
