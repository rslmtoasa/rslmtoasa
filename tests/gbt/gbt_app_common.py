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


def make_deck(dst, kind, rep, nk, q_dir=(0.0, 0.0, 0.0), order="first", spin_resolved=False):
    """kind: bccFe_prim, bccFe_sc2 (Fe, same potential on every site) or cr_afm (cr_bcc, Cr2 = Cr1 with mom = -z).
    q_dir is in direct coordinates of the cell; a gbt_single_q deck gets it as Cartesian q_ss with theta_ss = 0."""
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
    nml["calculation"].update(pre_processing="none", post_processing="fermi_surface")
    nml["control"].update(nsp=1 if spin_resolved else 3, recur="block")     # the spin-resolved export needs nsp = 1
    nml["hamiltonian"] = {"hoh": order == "second", "magnetic_representation": rep}
    if rep == "gbt_single_q":
        q_cart = np.linalg.inv(lattice_vectors(kind)) @ np.array(q_dir) if kind != "cr_afm" else np.array(q_dir)
        nml["hamiltonian"].update(q_ss=list(q_cart), theta_ss=0.0)
    nml["lattice"]["strux_want_sdot"] = False
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
