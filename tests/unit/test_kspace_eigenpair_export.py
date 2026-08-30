"""Contracts and reader round-trip checks for final k-space eigenpair output."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.read_kspace_eigenpairs import read_kspace_eigenpairs


def test_fermi_surface_postprocessing_and_export_controls_are_registered() -> None:
    namelist = (ROOT / "source/include_codes/namelists/reciprocal.f90").read_text(encoding="utf-8")
    reciprocal = (ROOT / "source/reciprocal.f90").read_text(encoding="utf-8")
    calculation = (ROOT / "source/calculation.f90").read_text(encoding="utf-8")
    calculation_reciprocal = (ROOT / "source/calculation_reciprocal.f90").read_text(encoding="utf-8")
    export = (ROOT / "source/reciprocal_export.f90").read_text(encoding="utf-8")
    self_source = (ROOT / "source/self.f90").read_text(encoding="utf-8")
    cmake = (ROOT / "source/CMakeLists.txt").read_text(encoding="utf-8")

    for name in (
        "fs_nk1",
        "fs_nk2",
        "fs_nk3",
        "fs_use_symmetry_reduction",
        "write_eigenpair_projections",
        "write_spin_resolved_eigenpairs",
        "eigenpair_output_file",
    ):
        assert name in namelist
        assert name in reciprocal
    assert "case ('fermi_surface')" in calculation
    assert "post_processing_fermi_surface" in calculation_reciprocal
    assert "write_kspace_eigenpairs" in calculation_reciprocal
    assert "site_" in export
    assert "_species =" in export
    assert "global_spin_blocks" in export
    assert "H/S has nonzero up/down coupling" in export
    assert "require nsp=1" in export
    assert "write_kspace_eigenpairs" not in self_source
    assert "write_eigenpairs" not in namelist
    assert "write_eigenpairs" not in reciprocal
    assert "reciprocal_export.f90" in cmake


def test_reader_preserves_fortran_shapes_and_complex_vectors(tmp_path: Path) -> None:
    base = tmp_path / "pairs"
    nk, nbands, nbasis, nsites = 3, 2, 4, 1
    k_points = np.arange(3 * nk, dtype=np.float64).reshape((3, nk), order="F")
    k_weights = np.full(nk, 1.0 / nk)
    eigenvalues = np.arange(nbands * nk, dtype=np.float64).reshape((nbands, nk), order="F")
    eigenvectors_real = np.arange(nbasis * nbands * nk, dtype=np.float64).reshape((nbasis, nbands, nk), order="F")
    eigenvectors_imag = -eigenvectors_real
    projections = np.arange(nsites * 4 * 2 * nbands * nk, dtype=np.float64).reshape(
        (nsites, 4, 2, nbands, nk), order="F"
    )

    (base.with_suffix(".meta")).write_text(
        "\n".join(
            [
                "n_kpoints = 3",
                "n_bands = 2",
                "n_basis = 4",
                "n_sites = 1",
                "projection_weights = enabled",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with base.with_suffix(".bin").open("wb") as stream:
        for array in (k_points, k_weights, eigenvalues, eigenvectors_real, eigenvectors_imag, projections):
            np.asarray(array, dtype=np.float64, order="F").ravel(order="F").tofile(stream)

    actual = read_kspace_eigenpairs(base)
    np.testing.assert_array_equal(actual["k_points"], k_points)
    np.testing.assert_array_equal(actual["k_weights"], k_weights)
    np.testing.assert_array_equal(actual["eigenvalues"], eigenvalues)
    np.testing.assert_array_equal(actual["eigenvectors"], eigenvectors_real + 1j * eigenvectors_imag)
    np.testing.assert_array_equal(actual["projection_weights"], projections)


def test_reader_loads_independently_solved_spin_payload(tmp_path: Path) -> None:
    base = tmp_path / "pairs"
    nk, nbands, nbasis, nsites = 2, 1, 2, 1
    k_points = np.zeros((3, nk), dtype=np.float64)
    k_weights = np.full(nk, 0.5, dtype=np.float64)
    eigenvalues = np.zeros((nbands, nk), dtype=np.float64)
    eigenvectors_real = np.ones((nbasis, nbands, nk), dtype=np.float64)
    eigenvectors_imag = np.zeros_like(eigenvectors_real)
    spin_eigenvalues = np.array(
        [[[1.0, 2.0]], [[3.0, 4.0]]],
        dtype=np.float64,
    )
    spin_vectors_real = np.arange(2 * nbasis * nbands * nk, dtype=np.float64).reshape(
        (2, nbasis, nbands, nk), order="F"
    )
    spin_vectors_imag = -spin_vectors_real

    (base.with_suffix(".meta")).write_text(
        "\n".join(
            [
                "n_kpoints = 2",
                "n_bands = 1",
                "n_basis = 2",
                "n_sites = 1",
                "spin_resolved_eigenpairs = enabled",
                "spin_resolved_binary_file = pairs.spin.bin",
                "spin_resolved_n_spin = 2",
                "spin_resolved_n_bands = 1",
                "spin_resolved_n_basis = 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with base.with_suffix(".bin").open("wb") as stream:
        for array in (k_points, k_weights, eigenvalues, eigenvectors_real, eigenvectors_imag):
            np.asarray(array, dtype=np.float64, order="F").ravel(order="F").tofile(stream)
    with (tmp_path / "pairs.spin.bin").open("wb") as stream:
        for array in (spin_eigenvalues, spin_vectors_real, spin_vectors_imag):
            np.asarray(array, dtype=np.float64, order="F").ravel(order="F").tofile(stream)

    actual = read_kspace_eigenpairs(base)
    np.testing.assert_array_equal(actual["spin_eigenvalues"], spin_eigenvalues)
    np.testing.assert_array_equal(actual["spin_eigenvectors"], spin_vectors_real + 1j * spin_vectors_imag)


if __name__ == "__main__":
    test_fermi_surface_postprocessing_and_export_controls_are_registered()
    print("k-space eigenpair export source contract PASS")
