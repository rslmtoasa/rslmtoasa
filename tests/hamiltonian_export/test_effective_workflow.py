"""Run mixed-boundary export/import and frozen block DOS through the CLI."""
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
                              fix_fermi=False)

        def setup(name, post_processing, recur=None):
            directory = work / name
            directory.mkdir()
            shutil.copy2(fixture / "Cu.nml", directory)
            base["calculation"]["post_processing"] = post_processing
            if recur:
                base["control"]["recur"] = recur
            base.write(directory / "input.nml")
            return directory

        if mode == "block-dos":
            for recur in ("chebyshev", "block", "lanczos"):
                directory = setup(recur, "hamiltonian_export", recur)
                run(binary, directory)
                assert metadata(directory)["scope"] == "cluster"
            reference = (work / "chebyshev/rs_tb.tb").read_bytes()
            for recur in ("block", "lanczos"):
                assert (work / recur / "rs_tb.tb").read_bytes() == reference
            assert metadata(work / "block")["fermi_provenance"] == "frozen_input"
            directory = setup("dos", "dos_frozen", "block")
            log = run(binary, directory)
            assert "Block recursion on progress" in log
            table = numbers(directory / "Cu_dos.out")
            assert table and all(math.isfinite(value) for row in table for value in row)
            assert any(row[1] > 0 for row in table)
            return

        native = setup("native", "conductivity")
        log = run(binary, native)
        fields = metadata(native)
        assert fields["scope"] == "cluster" and fields["periodic"].split() == ["T", "T", "F"]
        assert fields["fermi_provenance"] == "conductivity_resolved"
        reported = next(float(line.split("Free Fermi energy:")[1])
                        for line in log.splitlines() if "Free Fermi energy:" in line)
        assert abs(float(fields["fermi_Ry"]) - reported) <= 5.1e-7
        assert float(fields["mesh_fermi_Ry"]) == base["energy"]["fermi"]

        base["hamiltonian"]["export"] = "none"
        imported = setup("import", "conductivity_effective")
        for name in ("rs_tb.meta", "rs_tb.tb"):
            shutil.copy2(native / name, imported)
        run(binary, imported)
        for name in ("response_total.out", "response_orbitals_real.out", "response_orbitals_imag.out", "Cu_dos.out"):
            compare_tables(native / name, imported / name)

        meta = imported / "rs_tb.meta"
        original = meta.read_text()
        # The small native conductivity fixture can retain the input Fermi
        # level. Force distinct values in metadata to check independent mesh
        # alignment and response energies on the effective route.
        fermi = float(fields["fermi_Ry"])
        delta = 0.075
        meta.write_text(original.replace(fields["fermi_Ry"] + "\n", f"{fermi + delta:.17e}\n", 1))
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


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve(), sys.argv[4])
