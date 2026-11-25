"""Exercise external TB import with lattice.nml and no RS export."""
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

import f90nml
from test_effective_workflow import compare_tables, metadata, numbers, run


def make_deck(fixture, work, periodic):
    work.mkdir()
    for name in ("Fe.nml", "Co.nml", "lattice.nml"):
        shutil.copy2(fixture / name, work)
    deck = f90nml.read(fixture / "input.nml")
    deck["calculation"] = dict(post_processing="exchange_import", exchange_delta="onsite")
    deck["lattice"].update(rc=4, pbc=periodic, b1=periodic, b2=periodic, b3=periodic,
                          n1=1, n2=1, n3=1, njij=1, ijpair=[[1], [2]])
    deck["control"].update(recur="chebyshev", lld=24, cond_ll=16, cheb_backend="legacy",
                           transport_backend="legacy")
    deck["energy"].update(fermi=0.01, channels_ldos=400, energy_min=-1.0, energy_max=1.0)
    deck["hamiltonian"] = dict(effective_import="external")
    deck.write(work / "input.nml")
    return deck


def write_tb(work, nsites, image, reverse=False):
    # Magnetic sites have different onsite splittings. Add an off-diagonal
    # orbital term to the up-spin onsite block and connect two s orbitals.
    records = []
    for site, up, down in ((0, -0.8, 0.65), (1, -0.45, 0.95)):
        first = 9*site + 1
        spin_down = first + 9*nsites
        records.extend([(0, first, first, up), (0, spin_down, spin_down, down),
                        (0, first, first + 1, 0.15), (0, first + 1, first, 0.15)])
    for offset in (0, 9*nsites):
        records.extend([(image, 1 + offset, 10 + offset, -0.25),
                        (-image, 10 + offset, 1 + offset, -0.25)])
    n = 18*nsites
    if reverse:
        records = [(r, n + 1 - i, n + 1 - j, val) for r, i, j, val in records]
        (work / "rs_tb_orbital_map.dat").write_text("".join(f"{i} {n + 1 - i}\n" for i in range(1, n + 1)))
    (work / "paoham.dat").write_text("".join(f"{r} 0 0 {i} {j} {val} 0\n" for r, i, j, val in records))


def action_reference(probe, work, image):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    result = subprocess.run([probe, "import"], cwd=work, env=env,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout
    data = (work / "import.action").read_bytes()
    nb, nsites = struct.unpack_from("=ii", data)
    raw = struct.unpack_from(f"={2*5*3*nb*nsites}d", data, offset=8)
    values = list(map(complex, raw[::2], raw[1::2]))
    count = 3*nb*nsites
    psi = values[:count]
    lines = (work / "rs_tb.meta").read_text().splitlines()
    start = lines.index("positions_cart") + 1
    a = list(map(float, lines[lines.index("lattice_vectors_cart") + 1].split()))
    left, right = [list(map(float, line.split()[1:4])) for line in lines[start:start+2]]
    hopping = -0.25/float(metadata(work)["ry_to_eV"])

    for axis in range(3):
        displ = right[axis] - left[axis] + image*a[axis]
        expected = [0j]*count
        for rhs in range(3):
            for spin in (0, 9):
                i = spin + nb*rhs
                j = i + 3*nb
                expected[i] = 1j * displ * hopping * psi[j]
                expected[j] = -1j * displ * hopping * psi[i]
        actual = values[(axis+2)*count:(axis+3)*count]
        assert max(abs(x-y) for x,y in zip(expected, actual)) < 1e-13, (axis, image)


def main(binary, probe, fixture, root):
    with tempfile.TemporaryDirectory(prefix="external-tb-", dir=root) as tmp:
        for periodic in (False, True):
            work = Path(tmp) / str(periodic)
            deck = make_deck(fixture, work, periodic)
            # A diagonal TB record is sufficient for the first call to
            # generate the cluster site map. No native export is involved.
            (work / "paoham.dat").write_text("0 0 0 1 1 -0.4 0\n")
            deck["calculation"] = dict(post_processing="dos_import")
            deck.write(work / "input.nml", force=True)
            run(binary, work)
            nsites = int(metadata(work)["n_atoms"])
            assert nsites >= 2
            (work / "rs_tb.meta").unlink()

            image = -int(periodic)
            write_tb(work, nsites, image)
            action_reference(probe, work, image)
            deck["calculation"] = dict(post_processing="exchange_import", exchange_delta="onsite")
            deck.write(work / "input.nml", force=True)
            run(binary, work)
            jij = numbers(work / "jij.out")
            assert len(jij) == 1 and math.isfinite(jij[0][-2]) and abs(jij[0][-2]) > 1e-8

            mapped = Path(tmp) / f"mapped-{periodic}"
            make_deck(fixture, mapped, periodic)
            write_tb(mapped, nsites, image, reverse=True)
            action_reference(probe, mapped, image)
            run(binary, mapped)
            assert "orbital_index_map" in (mapped / "rs_tb.meta").read_text()
            for name in ("jij.out", "dij.out", "aij.out", "jijparts.out"):
                compare_tables(work / name, mapped / name)
            (mapped / "rs_tb_orbital_map.dat").unlink()
            deck["hamiltonian"]["effective_import"] = "legacy"
            deck.write(mapped / "input.nml", force=True)
            run(binary, mapped)
            for name in ("jij.out", "dij.out", "aij.out", "jijparts.out"):
                compare_tables(work / name, mapped / name)

            if periodic:
                assert nsites == 2
                # The one-cell hopping table also defines the compact model.
                sidecar = mapped / "rs_tb.meta"
                cluster_meta = sidecar.read_text()
                sidecar.write_text(cluster_meta.replace("scope cluster", "scope periodic", 1))
                run(binary, mapped)
                for name in ("jij.out", "dij.out", "aij.out", "jijparts.out"):
                    compare_tables(work / name, mapped / name)
                sidecar.write_text(cluster_meta)

            deck["calculation"] = dict(post_processing="conductivity_import")
            for path, mode in ((work, "external"), (mapped, "legacy")):
                deck["hamiltonian"]["effective_import"] = mode
                deck.write(path / "input.nml", force=True)
                run(binary, path)
            compare_tables(work / "response_total.out", mapped / "response_total.out")
            assert any(abs(row[1]) > 1e-8 for row in numbers(work / "response_total.out"))

            deck["hamiltonian"]["effective_import"] = "external"
            deck.write(mapped / "input.nml", force=True)
            (mapped / "rs_tb_orbital_map.dat").write_text("1 1\n2 1\n")
            run(binary, mapped, "Invalid orbital map")


if __name__ == "__main__":
    main(*(Path(arg).resolve() for arg in sys.argv[1:]))
