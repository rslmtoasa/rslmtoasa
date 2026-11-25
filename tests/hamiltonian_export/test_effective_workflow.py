"""Compare native and imported real-space post-processing through the CLI."""
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import f90nml


def metadata(directory):
    lines = (directory / "rs_tb.meta").read_text().splitlines()
    assert lines[0] == "# rs-tb-meta-v3"
    return dict(line.split(maxsplit=1) for line in lines if line and line[0].isalpha() and " " in line)


def run(binary, directory, error=None):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    result = subprocess.run([binary], cwd=directory, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory / "test.log").write_text(result.stdout)
    if error is None:
        assert result.returncode == 0, result.stdout
    else:
        assert result.returncode != 0 and error in result.stdout, result.stdout
    return result.stdout


def numbers(path):
    return [list(map(float, line.split())) for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]


def compare_tables(native, imported):
    ref = numbers(native)
    actual = numbers(imported)
    assert len(ref) == len(actual) and ref
    for first, second in zip(ref, actual):
        assert len(first) == len(second)
        assert all(math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-9) for x, y in zip(first, second)), (
            native.name, first, second
        )


def main(binary, fixture, root, mode):
    with tempfile.TemporaryDirectory(prefix=f"{mode}-", dir=root) as work:
        work = Path(work)
        base = f90nml.read(fixture / "input.nml")
        base["calculation"].update(pre_processing="none", post_processing="conductivity")
        base["lattice"].update(pbc=True, b1=True, b2=True, b3=False, n1=4, n2=4, n3=4)
        base["hamiltonian"].update(hoh=True, export="python", export_scope="auto")
        base["control"].update(nsp=1, recur="chebyshev", lld=16, cond_ll=16,
                               linear_in="charge", linear_out="charge",
                               cheb_backend="legacy", transport_backend="legacy")
        base["energy"].update(fermi=-0.3, energy_min=-2.0, energy_max=2.5, channels_ldos=400,
                              fix_fermi=True)

        def setup(name, post_processing, recur=None):
            directory = work / name
            directory.mkdir()
            shutil.copy2(fixture / "Cu.nml", directory)
            base["calculation"]["post_processing"] = post_processing
            if recur:
                base["control"]["recur"] = recur
            base.write(directory / "input.nml")
            return directory

        if mode == "compact-sheet":
            base["lattice"].update(n1=6, n2=6, n3=1)
            base["hamiltonian"].update(export="rs2pao", export_scope="cluster")
            cluster = setup("cluster", "hamiltonian_export")
            run(binary, cluster)
            base["hamiltonian"]["export_scope"] = "periodic"
            compact = setup("compact", "hamiltonian_export")
            # A rerun can publish over outputs left by an earlier run.
            for name in ("rs_tb.meta", "rs_tb.tb", "rs2paoham.dat"):
                (compact / name).write_text("stale export\n")
            run(binary, compact)
            assert not list(compact.glob("rs*.tmp"))
            fields = metadata(compact)
            assert fields["scope"] == "periodic" and fields["periodic"].split() == ["T", "T", "F"]
            assert int(fields["n_atoms"]) == 1 and int(metadata(cluster)["n_atoms"]) == 36
            translations = [tuple(map(int, line.split()[:3]))
                            for line in (compact / "rs2paoham.dat").read_text().splitlines()]
            assert any(x or y for x, y, _ in translations)
            assert all(z == 0 for _, _, z in translations)

            base["hamiltonian"].update(export="none", effective_import="legacy")
            outputs = []
            for source in (cluster, compact):
                directory = setup(f"replay-{source.name}", "conductivity_import")
                for name in ("rs_tb.meta", "rs2paoham.dat"):
                    shutil.copy2(source / name, directory)
                run(binary, directory)
                outputs.append(directory)
            for name in ("response_total.out", "response_orbitals_real.out",
                         "response_orbitals_imag.out", "Cu_dos.out"):
                compare_tables(outputs[0] / name, outputs[1] / name)

            # Retile the same frozen matrix. Independent native builds at
            # other replica sizes can produce different screened hoppings.
            base["lattice"]["n1"] = 8
            larger = setup("larger-replica", "dos_import")
            for name in ("rs_tb.meta", "rs2paoham.dat"):
                shutil.copy2(compact / name, larger)
            run(binary, larger)
            assert int((larger / "fort.805").read_text().splitlines()[0]) == 48
            assert any(row[1] > 0 for row in numbers(larger / "Cu_dos.out"))

            directory = setup("boundary-mismatch", "conductivity_import")
            shutil.copy2(compact / "rs2paoham.dat", directory)
            sidecar = (compact / "rs_tb.meta").read_text()
            record = next(line for line in sidecar.splitlines() if line.startswith("periodic "))
            (directory / "rs_tb.meta").write_text(sidecar.replace(record, "periodic T T T", 1))
            run(binary, directory, "Effective Hamiltonian boundary flags differ from input.")
            return

        if mode == "block-dos":
            for recur in ("chebyshev", "block", "lanczos"):
                directory = setup(recur, "hamiltonian_export", recur)
                run(binary, directory)
                assert metadata(directory)["scope"] == "cluster"
            reference = (work / "chebyshev/rs_tb.tb").read_bytes()
            for recur in ("block", "lanczos"):
                assert (work / recur / "rs_tb.tb").read_bytes() == reference
            assert metadata(work / "block")["fermi_provenance"] == "frozen_input"
            assert (work / "block/rs_tb_paoham.dat").is_file()
            base["hamiltonian"].update(export="rs2pao", export_scope="cluster")
            base["lattice"].update(b1=False, b2=False, b3=False)
            directory = setup("open", "hamiltonian_export")
            run(binary, directory)
            fields = metadata(directory)
            assert fields["scope"] == "cluster" and fields["periodic"].split() == ["F", "F", "F"]
            assert int(fields["n_atoms"]) > 1
            with (directory / "rs2paoham.dat").open() as payload:
                first = next(payload).split()
                assert len(first) == 7 and first[:3] == ["0", "0", "0"]
                assert all(line.split(maxsplit=3)[:3] == ["0", "0", "0"] for line in payload)
            base["hamiltonian"].update(export="python", export_scope="auto")
            base["lattice"].update(b1=True, b2=True)
            directory = setup("dos", "dos_frozen", "block")
            log = run(binary, directory)
            assert "Block recursion on progress" in log
            table = numbers(directory / "Cu_dos.out")
            assert table and all(math.isfinite(value) for row in table for value in row)
            assert any(row[1] > 0 for row in table)
            base["hamiltonian"].update(export="none")
            replay = setup("block-import", "dos_import")
            for name in ("rs_tb.meta", "rs_tb_paoham.dat"):
                shutil.copy2(work / "block" / name, replay)
            run(binary, replay)
            compare_tables(directory / "Cu_dos.out", replay / "Cu_dos.out")
            return

        base["hamiltonian"]["export"] = "rs2pao"
        native = setup("native", "conductivity")
        run(binary, native)
        # Export settings cannot trigger a writer in native conductivity.
        for name in ("rs_tb.meta", "rs_tb.tb", "rs2paoham.dat", "rs_tb.kspace"):
            assert not (native / name).exists()

        exported = setup("export", "hamiltonian_export")
        run(binary, exported)
        fields = metadata(exported)
        assert fields["scope"] == "cluster" and fields["periodic"].split() == ["T", "T", "F"]
        assert fields["fermi_provenance"] == "frozen_input"
        assert float(fields["fermi_Ry"]) == base["energy"]["fermi"]
        with (exported / "rs2paoham.dat").open() as payload:
            has_wrapped_bond = False
            count = 0
            for count, line in enumerate(payload, 1):
                record = line.split()
                assert len(record) == 7 and int(record[2]) == 0
                has_wrapped_bond |= int(record[0]) != 0 or int(record[1]) != 0
        assert has_wrapped_bond and count == int(fields["legacy_records"])

        base["hamiltonian"]["export"] = "none"
        imported = setup("import", "conductivity_import")
        for name in ("rs_tb.meta", "rs_tb.tb"):
            shutil.copy2(exported / name, imported)
        run(binary, imported)
        output_names = ("response_total.out", "response_orbitals_real.out", "response_orbitals_imag.out", "Cu_dos.out")
        for name in output_names:
            compare_tables(native / name, imported / name)

        legacy = setup("seven-column", "conductivity_import")
        for name in ("rs_tb.meta", "rs2paoham.dat"):
            shutil.copy2(exported / name, legacy)
        run(binary, legacy)
        for name in output_names:
            compare_tables(native / name, legacy / name)

        meta = imported / "rs_tb.meta"
        original = meta.read_text()
        # Older paired files can specify a separate DOS mesh reference.
        fermi = float(fields["fermi_Ry"])
        delta = 0.075
        fermi_line = next(line for line in original.splitlines() if line.startswith("fermi_Ry "))
        meta.write_text(original.replace(fermi_line + "\n",
                                         f"fermi_Ry {fermi + delta:.17e}\nmesh_fermi_Ry {fermi:.17e}\n", 1))
        run(binary, imported)
        old_dos = numbers(native / "Cu_dos.out")
        new_dos = numbers(imported / "Cu_dos.out")
        for old, new in zip(old_dos, new_dos):
            assert abs((old[0] - new[0]) - delta) < 1.1e-5, (old, new)
        meta.write_text(original.replace("spin_then_l_then_m_increasing", "incompatible"))
        run(binary, imported, "Invalid or incomplete rs_tb.meta")
        meta.write_text(original)
        payload = imported / "rs_tb.tb"
        lines = payload.read_text().splitlines(keepends=True)
        index = next(i for i, line in enumerate(lines) if not line.startswith("#"))
        fields = lines[index].split()
        fields[-2] = "1.25E-4"
        lines[index] = " ".join(fields) + "\n"
        payload.write_text("".join(lines))
        run(binary, imported, "payload does not match rs_tb.meta")


def magnetic_round_trip(binary, fixture, root):
    with tempfile.TemporaryDirectory(prefix="magnetic-", dir=root) as work:
        work = Path(work)
        base = f90nml.read(fixture / "input.nml")
        base["calculation"].update(pre_processing="none", processing="none")
        base["lattice"].update(rc=6)
        base["control"].update(nsp=4, recur="chebyshev", lld=18, cond_ll=14,
                               cond_calctype="per_type", cheb_backend="legacy", transport_backend="legacy")
        base["energy"].update(energy_min=-3.5, energy_max=3.5, channels_ldos=400, fix_fermi=True)
        base["hamiltonian"] = dict(hoh=True, export_scope="cluster")

        def setup(name, post_processing, source=None):
            directory = work / name
            directory.mkdir()
            shutil.copy2(fixture / "Fe.nml", directory)
            base["calculation"]["post_processing"] = post_processing
            base["hamiltonian"]["export"] = "rs2pao" if post_processing == "hamiltonian_export" else "none"
            base.write(directory / "input.nml")
            if source:
                for name in ("rs_tb.meta", "rs2paoham.dat"):
                    shutil.copy2(source / name, directory)
            return directory

        exported = setup("export", "hamiltonian_export")
        run(binary, exported)
        assert "onsite_soc_Ry " in (exported / "rs_tb.meta").read_text()
        fields = metadata(exported)
        assert int(fields["two_hop_paths"]) > int(fields["original_bond_records"])
        sites = [line.split() for line in (exported / "clust").read_text().splitlines()[1:]]
        # Select an offsite nearest neighbor and include an onsite pair.
        distances = [sum(float(x)**2 for x in site[:3]) for site in sites]
        neighbor = min((i for i, d in enumerate(distances) if d > 0), key=distances.__getitem__) + 1
        base["lattice"].update(njij=2, ijpair=[[1, 1], [neighbor, 1]])
        for recur in ("chebyshev", "block"):
            base["control"]["recur"] = recur
            native = setup(f"native-{recur}", "exchange")
            replay = setup(f"replay-{recur}", "exchange_import", exported)
            run(binary, native)
            run(binary, replay)
            for name in ("jij.out", "dij.out", "aij.out", "jijparts.out"):
                compare_tables(native / name, replay / name)
            pairs = numbers(replay / "jij.out")
            assert len(pairs) == 2 and abs(pairs[1][-2]) > 1e-6

        base["lattice"]["njij"] = 0
        del base["lattice"]["ijpair"]
        base["control"].update(recur="chebyshev", linear_in="spin_torque", linear_out="spin_soc_torque")
        base["hamiltonian"].update(pol_alpha="x", pol_beta="x")
        native = setup("native-torque", "conductivity")
        replay = setup("replay-torque", "conductivity_import", exported)
        run(binary, native)
        run(binary, replay)
        for name in ("response_total.out", "response_orbitals_real.out", "response_orbitals_imag.out"):
            compare_tables(native / name, replay / name)


def surface_round_trip(binary, fixture, root):
    with tempfile.TemporaryDirectory(prefix="surface-", dir=root) as work:
        work = Path(work)
        base = f90nml.read(fixture / "input.nml")
        base["calculation"].update(pre_processing="none", processing="none")
        base["lattice"].update(rc=8, nlay=2, ntype=3, ct=[2.7]*3, r2=9.0, strux_backend="legacy")
        base["atoms"]["label"] = ["Cu", "ES", "Cu-S"]
        base["control"].update(recur="chebyshev", lld=16, cond_ll=12, cond_calctype="per_type",
                               cheb_backend="legacy", transport_backend="legacy",
                               linear_in="charge", linear_out="orbital_accumulation")
        base["energy"].update(energy_min=-3.5, energy_max=3.5, channels_ldos=350, fix_fermi=True)
        base["hamiltonian"] = dict(hoh=True, export_scope="cluster", pol_alpha="y", v_beta=[1.0, 0.0, 0.0])

        def setup(route):
            directory = work / route
            directory.mkdir()
            for label in base["atoms"]["label"]:
                shutil.copy2(fixture / f"{label}.nml", directory)
            base["calculation"]["post_processing"] = route
            base["hamiltonian"]["export"] = "rs2pao" if route == "hamiltonian_export" else "none"
            base.write(directory / "input.nml")
            return directory

        exported = setup("hamiltonian_export")
        run(binary, exported)
        fields = metadata(exported)
        assert int(fields["effective_blocks"]) > int(fields["n_atoms"])
        native = setup("conductivity")
        replay = setup("conductivity_import")
        for name in ("rs_tb.meta", "rs2paoham.dat"):
            shutil.copy2(exported / name, replay)
        run(binary, native)
        run(binary, replay)
        for name in ("response_total.out", "response_orbitals_real.out", "response_orbitals_imag.out"):
            compare_tables(native / name, replay / name)
        assert any(abs(value) > 1e-6 for row in numbers(replay / "response_total.out") for value in row[1:])


if __name__ == "__main__":
    binary, fixture, root = (Path(arg).resolve() for arg in sys.argv[1:4])
    if sys.argv[4] == "magnetic":
        magnetic_round_trip(binary, fixture, root)
    elif sys.argv[4] == "surface":
        surface_round_trip(binary, fixture, root)
    else:
        main(binary, fixture, root, sys.argv[4])
