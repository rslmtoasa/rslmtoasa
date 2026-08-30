"""Pure-Python checks for the regular-grid preparation used by the viewer."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.plot_fermi_surface import (
    _cartesian_coordinates,
    _basis_projection_weights,
    _crossing_bands,
    _nesting_bands,
    _print_nesting_report,
    _nesting_sheet_data,
    _nesting_spectrum,
    _cross_spin_nesting_spectrum,
    _save_nesting_map,
    _spectral_refine_grid,
    _periodic_local_maximum_mask,
    _parser,
    _projection_values,
    _radial_nesting_average,
    _refine_axes,
    _render_window_size,
    _spin_sector_projection_weights,
    _site_species,
    _top_nesting_vectors,
    _unique_species,
    _values_on_grid,
)


def test_values_are_scattered_by_coordinates_not_file_order() -> None:
    dimensions = (2, 3, 2)
    coordinates = np.array(
        [
            [0, 0, 0],
            [1, 2, 1],
            [0, 1, 0],
            [1, 0, 1],
            [0, 2, 1],
            [1, 1, 0],
            [0, 0, 1],
            [1, 2, 0],
            [0, 1, 1],
            [1, 0, 0],
            [0, 2, 0],
            [1, 1, 1],
        ],
        dtype=float,
    ).T
    values = coordinates[0] + 10.0 * coordinates[1] + 100.0 * coordinates[2]
    actual = _values_on_grid(values, coordinates, dimensions)
    expected = np.empty(dimensions)
    for ix in range(dimensions[0]):
        for iy in range(dimensions[1]):
            for iz in range(dimensions[2]):
                expected[ix, iy, iz] = ix + 10.0 * iy + 100.0 * iz
    np.testing.assert_array_equal(actual, expected)


def test_fractional_mesh_is_affinely_transformed_to_cartesian_space() -> None:
    fractional = np.zeros((3, 2, 1, 1))
    fractional[:, 1, 0, 0] = [1.0, 0.0, 0.0]
    fractional[:, 0, 0, 0] = [0.0, 1.0, 0.0]
    reciprocal_vectors = np.array(
        [[2.0, 10.0, 20.0], [3.0, 11.0, 21.0], [4.0, 12.0, 22.0]],
    )
    cartesian = _cartesian_coordinates(fractional, reciprocal_vectors)
    np.testing.assert_array_equal(cartesian[:, 1, 0, 0], reciprocal_vectors[:, 0])
    np.testing.assert_array_equal(cartesian[:, 0, 0, 0], reciprocal_vectors[:, 1])


def test_cartesian_coordinates_are_the_viewer_default() -> None:
    parser = _parser()
    args = parser.parse_args(["fermi_surface"])
    assert args.coordinate_system == "cartesian"
    assert args.close_bz_seams
    assert not parser.parse_args(["fermi_surface", "--no-close-bz-seams"]).close_bz_seams
    assert parser.parse_args(["fermi_surface", "--periodic"]).close_bz_seams


def test_spin_sheet_modes_are_explicit_and_combined_is_the_default() -> None:
    parser = _parser()
    assert parser.parse_args(["fermi_surface"]).spin_sheets == "combined"
    assert parser.parse_args(["fermi_surface", "--spin-sheets", "up"]).spin_sheets == "up"
    assert parser.parse_args(["fermi_surface", "--spin-sheets", "both"]).spin_sheets == "both"


def test_nesting_pair_and_sheet_options_are_explicit() -> None:
    parser = _parser()
    args = parser.parse_args(
        [
            "fermi_surface",
            "--nesting-pairs",
            "interband",
            "--nesting-sheets",
            "up",
            "--nesting-output",
            "nesting.npz",
            "--nesting-only",
        ]
    )
    assert args.nesting_pairs == "interband"
    assert args.nesting_sheets == "up"
    assert args.nesting_output == Path("nesting.npz")
    assert args.nesting_only
    assert parser.parse_args(["fermi_surface", "--nesting-pairs", "cross-spin"]).nesting_pairs == "cross-spin"

    radial = parser.parse_args(
        ["fermi_surface", "--nesting-radial-output", "nesting.png", "--nesting-radial-bins", "12"]
    )
    assert radial.nesting_radial_output == Path("nesting.png")
    assert radial.nesting_radial_bins == 12
    assert parser.parse_args(["fermi_surface"]).nesting_interpolation == 1


def test_nesting_sheet_selection_uses_spin_payload_when_requested() -> None:
    main = np.arange(2.0 * 3.0).reshape((2, 3))
    spin = np.stack((main + 10.0, main + 20.0))
    data = {"eigenvalues": main, "spin_eigenvalues": spin, "k_points": np.zeros((3, 3))}

    selected, labels = _nesting_sheet_data(data, "down")
    np.testing.assert_array_equal(selected, spin[1])
    assert labels == ("down:1", "down:2")


def test_high_quality_render_scale_defaults_and_overrides() -> None:
    parser = _parser()
    high_quality = parser.parse_args(["fermi_surface", "--high-quality"])
    assert _render_window_size(high_quality) == (2400, 1800)

    custom = parser.parse_args(["fermi_surface", "--high-quality", "--render-scale", "1.5", "--window-size", "100", "101"])
    assert _render_window_size(custom) == (150, 152)


def test_smooth_interpolation_refines_each_cell_and_preserves_endpoints() -> None:
    axes = (np.array([0.0, 1.0, 2.0]), np.array([-1.0, 1.0]), np.array([0.0]))
    refined = _refine_axes(axes, 2)
    assert tuple(axis.size for axis in refined) == (5, 3, 1)
    assert refined[0][0] == axes[0][0]
    assert refined[0][-1] == axes[0][-1]


def test_projection_modes_sum_the_expected_channels() -> None:
    projections = np.zeros((2, 4, 2, 1, 3))
    projections[:, :, :, 0, :] = np.arange(projections[:, :, :, 0, :].size).reshape((2, 4, 2, 3))
    assert np.array_equal(_projection_values(projections, 0, "site", 1, 0, 0), projections[1].sum(axis=(0, 1))[0])
    assert np.array_equal(
        _projection_values(projections, 0, "orbital", 0, 2, 0), projections[:, 2].sum(axis=(0, 1))[0]
    )
    spin_up = projections[:, :, 0, 0].sum(axis=(0, 1))
    spin_down = projections[:, :, 1, 0].sum(axis=(0, 1))
    expected_polarization = (spin_up - spin_down) / (spin_up + spin_down)
    np.testing.assert_allclose(
        _projection_values(projections, 0, "spin", 0, 0, 1),
        expected_polarization,
    )
    assert np.array_equal(_projection_values(projections, 0, "channel", 1, 2, 1), projections[1, 2, 1, 0])


def test_spin_projection_has_fixed_signed_up_down_limits() -> None:
    projections = np.zeros((1, 1, 2, 1, 3))
    projections[0, 0, 0, 0] = [1.0, 0.0, 0.25]
    projections[0, 0, 1, 0] = [0.0, 1.0, 0.25]

    np.testing.assert_array_equal(
        _projection_values(projections, 0, "spin", 0, 0, 0),
        [1.0, -1.0, 0.0],
    )


def test_basis_projection_weights_follow_recorded_site_spin_orbital_order() -> None:
    eigenvectors = np.zeros((36, 1, 1), dtype=complex)
    eigenvectors[0, 0, 0] = 0.5  # site 1, spin-up, s
    eigenvectors[1, 0, 0] = 0.5  # site 1, spin-up, p
    eigenvectors[18 + 9 + 4, 0, 0] = np.sqrt(0.5)  # site 2, spin-down, d
    metadata = {
        "n_sites": "2",
        "basis_layout": "site_major_spin_blocked",
        "basis_site_block": "18",
        "basis_orbitals_per_spin": "9",
    }

    projections = _basis_projection_weights(eigenvectors, metadata)
    assert projections is not None
    np.testing.assert_allclose(projections[0, 0, 0, 0, 0], 0.25)
    np.testing.assert_allclose(projections[0, 1, 0, 0, 0], 0.25)
    np.testing.assert_allclose(projections[1, 2, 1, 0, 0], 0.5)
    np.testing.assert_allclose(projections.sum(axis=(0, 1, 2))[0, 0], 1.0)


def test_spin_sector_projection_weights_keep_the_selected_global_spin() -> None:
    eigenvectors = np.zeros((18, 1, 1), dtype=complex)
    eigenvectors[0, 0, 0] = 0.5
    eigenvectors[4, 0, 0] = np.sqrt(0.75)
    metadata = {
        "n_sites": "2",
        "basis_layout": "site_major_spin_blocked",
        "basis_orbitals_per_spin": "9",
    }

    projections = _spin_sector_projection_weights(eigenvectors, metadata, 1)
    assert projections is not None
    np.testing.assert_allclose(projections[0, 0, 1, 0, 0], 0.25)
    np.testing.assert_allclose(projections[0, 2, 1, 0, 0], 0.75)
    np.testing.assert_allclose(projections[:, :, 0], 0.0)


def test_all_band_selection_and_species_projection() -> None:
    eigenvalues = np.array(
        [[-1.0, 1.0, 0.0], [2.0, 3.0, 4.0], [-3.0, -2.0, -1.0]],
    )
    np.testing.assert_array_equal(_crossing_bands(eigenvalues, 0.0), np.array([0]))

    projections = np.zeros((3, 4, 2, 1, 2))
    projections[0, :, :, 0, :] = 1.0
    projections[1, :, :, 0, :] = 2.0
    projections[2, :, :, 0, :] = 4.0
    labels = _site_species({"site_1_species": "Fe", "site_2_species": "Ni", "site_3_species": "Fe"}, 3)
    assert labels == ("Fe", "Ni", "Fe")
    assert _unique_species(labels) == ("Fe", "Ni")
    np.testing.assert_array_equal(
        _projection_values(projections, 0, "species", 0, 0, 0, "Fe", labels),
        np.full(2, 40.0),
    )


def test_nesting_band_selection_keeps_five_sigma_shells_separate_from_fs_selection() -> None:
    eigenvalues = np.array([[0.02, 0.02], [0.2, 0.3], [-0.2, -0.3]])
    np.testing.assert_array_equal(_crossing_bands(eigenvalues, 0.0), np.array([], dtype=int))
    np.testing.assert_array_equal(_nesting_bands(eigenvalues, 0.0, 0.01), np.array([0]))


def test_nesting_estimate_uses_periodic_fractional_mesh() -> None:
    dimensions = (4, 4, 4)
    axes = [np.arange(size, dtype=float) / size - 0.5 + 0.5 / size for size in dimensions]
    mesh = np.stack(np.meshgrid(*axes, indexing="ij"), axis=0)
    k_points = mesh.reshape((3, -1), order="F")
    eigenvalues = mesh[0].reshape((1, -1), order="F")

    nesting, q_axes, bands = _nesting_spectrum(eigenvalues, k_points, dimensions, 0.0, 0.05)
    assert nesting.shape == dimensions
    assert bands.tolist() == [0]
    assert all(axis.shape == (4,) for axis in q_axes)
    top = _top_nesting_vectors(nesting, q_axes, 3)
    assert top
    assert all(np.isfinite(score) for _, score in top)
    assert all(not np.allclose(vector, 0.0) for vector, _ in top)


def test_nesting_pair_channels_are_normalized_at_q_zero() -> None:
    dimensions = (4, 4, 4)
    axes = [np.arange(size, dtype=float) / size - 0.5 + 0.5 / size for size in dimensions]
    mesh = np.stack(np.meshgrid(*axes, indexing="ij"), axis=0)
    k_points = mesh.reshape((3, -1), order="F")
    eigenvalues = np.vstack(
        (
            mesh[0].reshape((1, -1), order="F"),
            mesh[1].reshape((1, -1), order="F"),
        )
    )

    all_pairs, _, _ = _nesting_spectrum(eigenvalues, k_points, dimensions, 0.0, 0.05, "all")
    intraband, _, _ = _nesting_spectrum(eigenvalues, k_points, dimensions, 0.0, 0.05, "intraband")
    interband, _, _ = _nesting_spectrum(eigenvalues, k_points, dimensions, 0.0, 0.05, "interband")
    center = tuple(size // 2 for size in dimensions)
    assert all_pairs[center] == 1.0
    assert intraband[center] == 1.0
    assert interband[center] == 1.0
    assert not np.allclose(all_pairs, intraband)
    assert np.min(interband) >= -1.0e-12


def test_spectral_interpolation_preserves_periodic_grid_samples() -> None:
    axis = np.arange(4, dtype=float) / 4.0
    mesh = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=0)
    values = np.sin(2.0 * np.pi * mesh[0]) + 0.5 * np.cos(2.0 * np.pi * mesh[1])
    refined = _spectral_refine_grid(values, 2)
    np.testing.assert_allclose(refined[::2, ::2, ::2], values, atol=1.0e-12)


def test_cross_spin_nesting_keeps_the_ordered_q_direction() -> None:
    dimensions = (4, 1, 1)
    axes = (np.arange(4, dtype=float) / 4.0, np.array([0.0]), np.array([0.0]))
    mesh = np.stack(np.meshgrid(*axes, indexing="ij"), axis=0)
    k_points = mesh.reshape((3, -1), order="F")
    spin_eigenvalues = np.array(
        [
            [[0.0], [1.0], [1.0], [1.0]],
            [[1.0], [0.0], [1.0], [1.0]],
        ]
    ).reshape((2, 1, 4))

    nesting, _, _, _ = _cross_spin_nesting_spectrum(
        spin_eigenvalues,
        k_points,
        dimensions,
        0.0,
        0.01,
    )
    assert nesting[3, 0, 0] == 1.0
    assert nesting[1, 0, 0] != nesting[3, 0, 0]


def test_cross_spin_report_keeps_up_and_down_band_labels(capsys) -> None:
    dimensions = (2, 2, 2)
    axis = np.arange(2, dtype=float) / 2.0
    mesh = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=0)
    k_points = mesh.reshape((3, -1), order="F")
    eigenvalues = mesh[0].reshape((1, -1), order="F")
    data = {
        "metadata": {
            "k_mesh": "2 2 2",
            "fermi_level": "0.0",
            "reciprocal_b1": "1 0 0",
            "reciprocal_b2": "0 1 0",
            "reciprocal_b3": "0 0 1",
        },
        "k_points": k_points,
        "eigenvalues": eigenvalues,
        "spin_eigenvalues": np.stack((eigenvalues, eigenvalues)),
    }
    args = _parser().parse_args(
        [
            "fermi_surface",
            "--nesting-pairs",
            "cross-spin",
            "--nesting-sheets",
            "both",
            "--nesting-width",
            "0.5",
            "--nesting-top",
            "0",
        ]
    )

    _print_nesting_report(args, data)
    assert "bands contributing near EF=up:1,down:1" in capsys.readouterr().out


def test_periodic_local_maximum_mask_wraps_across_the_bz_boundary() -> None:
    values = np.zeros((4, 1, 1))
    values[0, 0, 0] = 2.0
    values[1, 0, 0] = 1.0
    values[3, 0, 0] = 1.0
    mask = _periodic_local_maximum_mask(values)
    assert mask[0, 0, 0]
    assert not mask[1, 0, 0]
    assert not mask[3, 0, 0]


def test_zero_nesting_map_has_no_reportable_peaks() -> None:
    nesting = np.zeros((2, 2, 2))
    q_axes = (np.array([-0.5, 0.0]),) * 3
    assert _top_nesting_vectors(nesting, q_axes, 3) == []


def test_radial_nesting_average_uses_cartesian_q_lengths() -> None:
    q_axes = (np.array([0.0, 1.0]), np.array([0.0]), np.array([0.0]))
    nesting = np.array([[[2.0]], [[4.0]]])
    radii, values, counts = _radial_nesting_average(nesting, q_axes, np.diag([2.0, 3.0, 4.0]), 2)

    np.testing.assert_allclose(radii, [0.0, 2.0])
    np.testing.assert_allclose(values, [2.0, 4.0])
    np.testing.assert_array_equal(counts, [1, 1])


def test_nesting_map_export_contains_aligned_fractional_and_cartesian_grids(tmp_path: Path) -> None:
    q_axes = (np.array([-0.5, 0.0]), np.array([-0.25, 0.25]), np.array([0.0]))
    nesting = np.arange(4.0).reshape((2, 2, 1))
    reciprocal_vectors = np.diag([2.0, 3.0, 4.0])
    output = tmp_path / "nesting.npz"

    _save_nesting_map(
        output,
        nesting,
        q_axes,
        np.array([0, 3]),
        ("up:1", "up:2", "down:1", "down:2"),
        0.0,
        0.05,
        "all",
        "both",
        reciprocal_vectors,
    )

    with np.load(output) as saved:
        assert saved["nesting"].shape == (2, 2, 1)
        assert saved["q_fractional"].shape == (3, 2, 2, 1)
        assert saved["q_cartesian"].shape == (3, 2, 2, 1)
        assert saved["crossing_band_labels"].tolist() == ["up:1", "down:2"]
        assert saved["contributing_band_labels"].tolist() == ["up:1", "down:2"]
        np.testing.assert_array_equal(saved["q_cartesian"][:, 0, 1, 0], [-1.0, 0.75, 0.0])


def test_nonfinite_bands_are_not_selected_for_the_fs() -> None:
    eigenvalues = np.array([[np.nan, 1.0], [np.inf, -np.inf], [-1.0, 1.0]])
    np.testing.assert_array_equal(_crossing_bands(eigenvalues, 0.0), np.array([2]))


if __name__ == "__main__":
    test_values_are_scattered_by_coordinates_not_file_order()
    test_projection_modes_sum_the_expected_channels()
    print("Fermi-surface plot grid helpers PASS")
