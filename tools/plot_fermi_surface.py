#!/usr/bin/env python3
"""Interactive PyVista viewer for an RSLMTO Fermi-surface export.

The input is the base name written by the ``fermi_surface`` post-processing
route, for example ``fermi_surface`` for ``fermi_surface.bin`` and
``fermi_surface.meta``.  The viewer contours one selected band, or every band
crossing the chosen Fermi level, on the regular k-mesh.  Optional site,
orbital, spin, and species projections can be used for colouring.
Collinear exports with a ``.spin.bin`` companion also support independently
solved ``up``, ``down``, or ``both`` global spin sheets.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from typing import Any

import numpy as np

try:
    from .read_kspace_eigenpairs import read_kspace_eigenpairs
except ImportError:  # direct execution: python tools/plot_fermi_surface.py ...
    from read_kspace_eigenpairs import read_kspace_eigenpairs


ORBITALS = ("s", "p", "d", "f")
SPINS = ("up", "down")


def _mesh_dimensions(metadata: dict[str, str], nk: int) -> tuple[int, int, int]:
    """Return the regular mesh dimensions recorded by the Fortran writer."""

    try:
        dimensions = tuple(int(value) for value in metadata["k_mesh"].split())
    except (KeyError, ValueError) as exc:
        raise ValueError("metadata is missing a valid 'k_mesh = nx ny nz' entry") from exc
    if len(dimensions) != 3 or any(value < 1 for value in dimensions):
        raise ValueError(f"invalid k_mesh dimensions: {dimensions}")
    if int(np.prod(dimensions)) != nk:
        raise ValueError(
            f"the viewer needs a complete regular mesh: k_mesh={dimensions} but payload has {nk} points"
        )
    return dimensions


def _regular_grid_index_map(
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
) -> tuple[tuple[np.ndarray, ...], np.ndarray]:
    """Return regular-grid axes and a cached Fortran-order scatter map."""

    axes: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    for coordinate, size in zip(k_points, dimensions):
        axis = np.unique(np.round(coordinate, decimals=12))
        if axis.size != size:
            raise ValueError(
                f"coordinate axis has {axis.size} unique values, expected {size}; "
                "the export is not a regular full mesh"
            )
        index = np.searchsorted(axis, np.round(coordinate, decimals=12))
        if not np.allclose(coordinate, axis[index], rtol=0.0, atol=1.0e-10):
            raise ValueError("k-point coordinates do not map cleanly onto a regular mesh")
        axes.append(axis)
        indices.append(index)

    ix, iy, iz = indices
    nx, ny, nz = dimensions
    flat_indices = ix + nx * (iy + ny * iz)
    if np.unique(flat_indices).size != nx * ny * nz:
        raise ValueError("k-point payload contains duplicate or missing regular-grid points")
    return tuple(axes), flat_indices.astype(np.int64, copy=False)


def _regular_axes(k_points: np.ndarray, dimensions: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    """Find and validate the three coordinate axes, independent of point order."""

    axes, _ = _regular_grid_index_map(k_points, dimensions)
    return axes


def _values_on_grid(
    values: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    grid_indices: np.ndarray | None = None,
) -> np.ndarray:
    """Scatter values stored as (k-point,) onto an (nx, ny, nz) NumPy grid."""

    if grid_indices is None:
        _, grid_indices = _regular_grid_index_map(k_points, dimensions)
    result = np.empty(int(np.prod(dimensions)), dtype=np.float64)
    result[grid_indices] = np.asarray(values, dtype=np.float64)
    return result.reshape(dimensions, order="F")


def _refine_axes(axes: tuple[np.ndarray, ...], factor: int) -> tuple[np.ndarray, ...]:
    """Refine regular axes while preserving their endpoints."""

    if factor < 1:
        raise ValueError("smooth interpolation factor must be at least one")
    if factor == 1:
        return axes
    refined: list[np.ndarray] = []
    for axis in axes:
        if axis.size == 1:
            refined.append(axis.copy())
            continue
        size = (axis.size - 1) * factor + 1
        refined.append(np.linspace(axis[0], axis[-1], size))
    return tuple(refined)


def _smooth_grid(values: np.ndarray, factor: int) -> np.ndarray:
    """Cubicly interpolate a scalar grid onto a uniformly refined grid."""

    if factor == 1:
        return values
    try:
        from scipy.ndimage import map_coordinates
    except ImportError as exc:
        raise RuntimeError(
            "--smooth-interpolation needs SciPy; install it with 'python -m pip install scipy'"
        ) from exc

    old_shape = values.shape
    new_shape = tuple((size - 1) * factor + 1 for size in old_shape)
    coordinates = np.meshgrid(
        *(np.linspace(0.0, size - 1.0, new_size) for size, new_size in zip(old_shape, new_shape)),
        indexing="ij",
    )
    order = min(3, *(size - 1 for size in old_shape))
    return map_coordinates(values, coordinates, order=order, mode="nearest", prefilter=order > 1)


def _spectral_refine_grid(values: np.ndarray, factor: int) -> np.ndarray:
    """Periodically interpolate a scalar grid by zero-padding Fourier modes."""

    if factor < 1:
        raise ValueError("nesting interpolation factor must be at least one")
    if factor == 1:
        return np.asarray(values, dtype=np.float64)
    old_shape = tuple(int(size) for size in values.shape)
    if any(size < 2 for size in old_shape):
        raise ValueError("nesting spectral interpolation needs at least two k-points per direction")
    new_shape = tuple(size * factor for size in old_shape)
    shifted = np.fft.fftshift(np.fft.fftn(np.asarray(values, dtype=np.float64)))
    padded = np.zeros(new_shape, dtype=np.complex128)
    slices = tuple(
        slice((new_size - old_size) // 2, (new_size - old_size) // 2 + old_size)
        for old_size, new_size in zip(old_shape, new_shape)
    )
    padded[slices] = shifted
    scale = float(np.prod(new_shape)) / float(np.prod(old_shape))
    return (np.fft.ifftn(np.fft.ifftshift(padded)).real * scale).astype(np.float64, copy=False)


def _spectral_refine_nesting_mesh(
    eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    factor: int,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int]]:
    """Spectrally refine all band energies and construct the refined k-grid."""

    if factor < 1:
        raise ValueError("nesting interpolation factor must be at least one")
    if factor == 1:
        return eigenvalues, k_points, dimensions
    axes, grid_indices = _regular_grid_index_map(k_points, dimensions)
    refined_axes: list[np.ndarray] = []
    for axis in axes:
        if axis.size < 2:
            raise ValueError("nesting spectral interpolation needs at least two k-points per direction")
        step = float(axis[1] - axis[0])
        if not np.allclose(np.diff(axis), step, rtol=0.0, atol=1.0e-10):
            raise ValueError("nesting spectral interpolation requires uniformly spaced k-point axes")
        refined_axes.append(axis[0] + np.arange(axis.size * factor, dtype=np.float64) * step / factor)
    refined_dimensions = tuple(axis.size for axis in refined_axes)
    refined_grid = np.empty((eigenvalues.shape[0], *refined_dimensions), dtype=np.float64)
    for band in range(eigenvalues.shape[0]):
        grid_values = _values_on_grid(eigenvalues[band], k_points, dimensions, grid_indices)
        refined_grid[band] = _spectral_refine_grid(grid_values, factor)
    refined_mesh = np.stack(np.meshgrid(*refined_axes, indexing="ij"), axis=0)
    refined_k_points = refined_mesh.reshape((3, -1), order="F")
    return refined_grid.reshape((eigenvalues.shape[0], -1), order="F"), refined_k_points, refined_dimensions


def _periodic_extend(values: np.ndarray, axes: tuple[np.ndarray, ...]) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    """Add a translated endpoint in each direction to close periodic seams."""

    dimensions = values.shape
    extended_indices = [np.arange(size + 1) % size for size in dimensions]
    extended_values = values[np.ix_(*extended_indices)]
    extended_axes: list[np.ndarray] = []
    for axis in axes:
        if axis.size == 1:
            step = 1.0
        else:
            step = float(axis[1] - axis[0])
            if not np.allclose(np.diff(axis), step, rtol=0.0, atol=1.0e-10):
                raise ValueError("periodic plotting requires uniformly spaced k-point axes")
        extended_axes.append(np.concatenate((axis, [axis[0] + axis.size * step])))
    return extended_values, tuple(extended_axes)


def _metadata_reciprocal_vectors(metadata: dict[str, str]) -> np.ndarray | None:
    """Read reciprocal vectors as columns, returning None for older sidecars."""

    vectors = []
    for index in range(1, 4):
        key = f"reciprocal_b{index}"
        if key not in metadata:
            return None
        values = np.fromstring(metadata[key], sep=" ")
        if values.size != 3:
            return None
        vectors.append(values)
    return np.column_stack(vectors)


def _basis_projection_weights(
    eigenvectors: np.ndarray,
    metadata: dict[str, str],
) -> np.ndarray | None:
    """Project eigenvectors using the recorded site-major spin-blocked basis.

    The normal-state basis is ordered per site as
    ``[s,p,d,f orbitals for spin up, s,p,d,f orbitals for spin down]``.
    This gives direct basis-spin weights even when optional local-axis
    projection weights were not written by the Fortran exporter.
    """

    try:
        nsites = int(metadata["n_sites"])
        site_block = int(metadata.get("basis_site_block", "0"))
        norb_per_spin = int(metadata.get("basis_orbitals_per_spin", "0"))
    except (KeyError, ValueError):
        return None
    if metadata.get("basis_layout", "site_major_spin_blocked") != "site_major_spin_blocked":
        return None
    if nsites < 1:
        return None
    nbasis = eigenvectors.shape[0]
    if site_block < 1:
        if nbasis % nsites != 0:
            return None
        site_block = nbasis // nsites
    if norb_per_spin < 1:
        if site_block % 2 != 0:
            return None
        norb_per_spin = site_block // 2
    if site_block != 2 * norb_per_spin or site_block * nsites != nbasis:
        return None

    amplitudes = np.abs(eigenvectors) ** 2
    projections = np.zeros(
        (nsites, 4, 2, eigenvectors.shape[1], eigenvectors.shape[2]),
        dtype=np.float64,
    )
    orbital_first = (0, 1, 4, 9)
    orbital_last = (1, 4, 9, 16)
    for site in range(nsites):
        site_offset = site * site_block
        for orbital, (first, last) in enumerate(zip(orbital_first, orbital_last)):
            if first >= norb_per_spin:
                continue
            last = min(last, norb_per_spin)
            for spin in range(2):
                start = site_offset + spin * norb_per_spin + first
                stop = site_offset + spin * norb_per_spin + last
                projections[site, orbital, spin] = amplitudes[start:stop].sum(axis=0)
    return projections


def _spin_sector_projection_weights(
    eigenvectors: np.ndarray,
    metadata: dict[str, str],
    spin: int,
) -> np.ndarray | None:
    """Build site/orbital/spin weights for one reduced global spin sector."""

    if spin not in (0, 1):
        raise ValueError("spin sector must be 0 (up) or 1 (down)")
    try:
        nsites = int(metadata["n_sites"])
        norb_per_spin = int(metadata.get("basis_orbitals_per_spin", "0"))
    except (KeyError, ValueError):
        return None
    if metadata.get("basis_layout", "site_major_spin_blocked") != "site_major_spin_blocked":
        return None
    if nsites < 1:
        return None
    if norb_per_spin < 1:
        if eigenvectors.shape[0] % nsites != 0:
            return None
        norb_per_spin = eigenvectors.shape[0] // nsites
    if nsites * norb_per_spin != eigenvectors.shape[0]:
        return None

    amplitudes = np.abs(eigenvectors) ** 2
    projections = np.zeros(
        (nsites, 4, 2, eigenvectors.shape[1], eigenvectors.shape[2]),
        dtype=np.float64,
    )
    orbital_first = (0, 1, 4, 9)
    orbital_last = (1, 4, 9, 16)
    for site in range(nsites):
        site_offset = site * norb_per_spin
        for orbital, (first, last) in enumerate(zip(orbital_first, orbital_last)):
            if first >= norb_per_spin:
                continue
            last = min(last, norb_per_spin)
            start = site_offset + first
            stop = site_offset + last
            projections[site, orbital, spin] = amplitudes[start:stop].sum(axis=0)
    return projections


def _cartesian_coordinates(fractional_mesh: np.ndarray, reciprocal_vectors: np.ndarray) -> np.ndarray:
    """Transform fractional reciprocal coordinates using basis vectors as columns."""

    vectors = np.asarray(reciprocal_vectors, dtype=np.float64)
    if vectors.shape != (3, 3):
        raise ValueError(f"reciprocal vectors must have shape (3, 3), got {vectors.shape}")
    return np.einsum("ij,j...->i...", vectors, fractional_mesh)


def _site_species(metadata: dict[str, str], nsites: int) -> tuple[str, ...]:
    """Return species labels for projection sites, with an old-file fallback."""

    labels = []
    for index in range(1, nsites + 1):
        label = metadata.get(f"site_{index}_species", "")
        labels.append(label.strip() or f"site {index}")
    return tuple(labels)


def _unique_species(site_species: tuple[str, ...]) -> tuple[str, ...]:
    """Keep species in first-occurrence order for stable UI selection."""

    return tuple(dict.fromkeys(site_species))


def _crossing_bands(eigenvalues: np.ndarray, fermi: float) -> np.ndarray:
    """Return zero-based bands whose sampled energy range brackets ``fermi``."""

    finite = np.all(np.isfinite(eigenvalues), axis=1)
    lower = np.min(eigenvalues, axis=1)
    upper = np.max(eigenvalues, axis=1)
    return np.flatnonzero(finite & (lower <= fermi) & (upper >= fermi))


def _nesting_bands(eigenvalues: np.ndarray, fermi: float, width: float) -> np.ndarray:
    """Return bands within five Gaussian widths of EF on the sampled mesh."""

    if not np.isfinite(width) or width <= 0.0:
        raise ValueError("nesting width must be positive")
    finite = np.all(np.isfinite(eigenvalues), axis=1)
    lower = np.min(eigenvalues, axis=1)
    upper = np.max(eigenvalues, axis=1)
    margin = 5.0 * width
    return np.flatnonzero(finite & (lower <= fermi + margin) & (upper >= fermi - margin))


def _snap_integer_slider(value: float, widget: Any, minimum: int, maximum: int) -> int:
    """Snap a PyVista slider to an integer and return its bounded value."""

    integer = int(np.clip(np.rint(value), minimum, maximum))
    representation = widget.GetRepresentation()
    if representation.GetValue() != float(integer):
        representation.SetValue(float(integer))
    return integer


def _nesting_weights(
    eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    fermi: float,
    width: float,
) -> tuple[np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
    """Return Gaussian-delta weights for every sampled band crossing EF."""

    if width <= 0.0:
        raise ValueError("nesting width must be positive")
    axes, grid_indices = _regular_grid_index_map(k_points, dimensions)
    bands = _nesting_bands(eigenvalues, fermi, width)
    normalization = width * np.sqrt(2.0 * np.pi)
    weights = np.empty((bands.size, *dimensions), dtype=np.float64)
    for index, band in enumerate(bands):
        energy = _values_on_grid(eigenvalues[band], k_points, dimensions, grid_indices)
        weights[index] = np.exp(-0.5 * ((energy - fermi) / width) ** 2) / normalization
    return bands, weights, axes


def _typical_nesting_energy_step(
    eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    fermi: float,
    width: float,
) -> float | None:
    """Estimate the energy change between neighbouring raw mesh points."""

    bands = _nesting_bands(eigenvalues, fermi, width)
    if bands.size == 0:
        return None
    _, grid_indices = _regular_grid_index_map(k_points, dimensions)
    steps: list[np.ndarray] = []
    for band in bands:
        energy = _values_on_grid(eigenvalues[band], k_points, dimensions, grid_indices)
        for axis in range(3):
            difference = np.abs(np.roll(energy, -1, axis=axis) - energy).ravel()
            finite = difference[np.isfinite(difference) & (difference > 0.0)]
            if finite.size:
                steps.append(finite)
    if not steps:
        return None
    return float(np.median(np.concatenate(steps)))


def _warn_if_nesting_width_is_narrow(
    eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    fermi: float,
    width: float,
) -> None:
    """Warn when Gaussian shells are substantially narrower than mesh steps."""

    typical_step = _typical_nesting_energy_step(eigenvalues, k_points, dimensions, fermi, width)
    if typical_step is not None and width < 0.5 * typical_step:
        warnings.warn(
            f"nesting Gaussian width ({width:.6g} Ry) is smaller than half the "
            f"typical raw-mesh energy step ({typical_step:.6g} Ry); increase "
            "--nesting-width or use --nesting-interpolation to reduce mesh dependence",
            RuntimeWarning,
            stacklevel=3,
        )


def _sum_rfft_weights(weights: np.ndarray) -> np.ndarray:
    """Accumulate real-grid weights in the compact real-FFT representation."""

    dimensions = weights.shape[1:]
    transformed = np.zeros((dimensions[0], dimensions[1], dimensions[2] // 2 + 1), dtype=np.complex128)
    for weight in weights:
        transformed += np.fft.rfftn(weight)
    return transformed


def _nesting_from_weights(
    weights: np.ndarray,
    dimensions: tuple[int, int, int],
    pair_mode: str,
) -> np.ndarray:
    """Build an all-, same-index-, or different-index pair correlation."""

    total = _sum_rfft_weights(weights)
    if pair_mode == "all":
        spectrum = np.conj(total) * total
    else:
        same_index = np.zeros_like(total.real)
        for weight in weights:
            transformed = np.fft.rfftn(weight)
            same_index += np.abs(transformed) ** 2
        if pair_mode == "intraband":
            spectrum = same_index
        else:
            spectrum = np.conj(total) * total - same_index
    return np.fft.irfftn(spectrum, s=dimensions, axes=(0, 1, 2)).real


def _normalize_nesting(nesting: np.ndarray) -> np.ndarray:
    """Normalize a nesting map at q=0, with a maximum fallback."""

    nesting = np.asarray(nesting, dtype=np.float64)
    tolerance = 1.0e-12 * max(1.0, float(np.max(np.abs(nesting))))
    nesting[(nesting < 0.0) & (nesting > -tolerance)] = 0.0
    q0 = float(nesting[(0,) * len(nesting.shape)]) if nesting.size else 0.0
    if q0 > tolerance:
        nesting /= q0
    elif nesting.size:
        maximum = float(np.max(nesting))
        if maximum > 0.0:
            nesting /= maximum
    return np.fft.fftshift(nesting)


def _nesting_q_axes(axes: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    """Return the FFT q-grid axes corresponding to sampled k-grid axes."""

    q_axes: list[np.ndarray] = []
    for axis in axes:
        if axis.size == 1:
            step = 1.0
        else:
            step = float(axis[1] - axis[0])
        q_axes.append(np.fft.fftshift(np.fft.fftfreq(axis.size) * axis.size * step))
    return tuple(q_axes)


def _q_magnitude_grid(
    q_axes: tuple[np.ndarray, ...],
    reciprocal_vectors: np.ndarray,
) -> np.ndarray:
    """Return Cartesian q magnitudes on the nesting grid."""

    q_fractional = np.stack(np.meshgrid(*q_axes, indexing="ij"), axis=0)
    q_cartesian = _cartesian_coordinates(q_fractional, reciprocal_vectors)
    return np.sqrt(np.sum(q_cartesian**2, axis=0))


def _angular_nesting_baseline(
    nesting: np.ndarray,
    q_axes: tuple[np.ndarray, ...],
    reciprocal_vectors: np.ndarray,
    bins: int = 64,
) -> np.ndarray:
    """Return the Cartesian-shell average used to remove the small-q envelope."""

    if bins < 1:
        raise ValueError("nesting radial bin count must be positive")
    magnitudes = _q_magnitude_grid(q_axes, reciprocal_vectors)
    if magnitudes.shape != nesting.shape:
        raise ValueError("nesting map and q-grid have different shapes")
    flat_magnitudes = magnitudes.ravel()
    flat_nesting = np.asarray(nesting, dtype=np.float64).ravel()
    maximum = float(np.max(flat_magnitudes)) if flat_magnitudes.size else 0.0
    if maximum <= 0.0:
        return np.full(nesting.shape, float(np.mean(flat_nesting)), dtype=np.float64)
    edges = np.linspace(0.0, maximum, bins + 1)
    indices = np.floor(flat_magnitudes / maximum * bins).astype(np.int64)
    indices = np.clip(indices, 0, bins - 1)
    totals = np.bincount(indices, weights=flat_nesting, minlength=bins)
    counts = np.bincount(indices, minlength=bins)
    shell_average = np.divide(
        totals,
        counts,
        out=np.ones_like(totals),
        where=counts > 0,
    )
    return shell_average[indices].reshape(nesting.shape)


def _periodic_local_maximum_mask(values: np.ndarray) -> np.ndarray:
    """Find local maxima on a periodic grid without a SciPy dependency."""

    values = np.asarray(values, dtype=np.float64)
    mask = np.isfinite(values) & (values > 0.0)
    axes = tuple(range(values.ndim))
    for offset_index in np.ndindex(*(3 for _ in axes)):
        offset = tuple(value - 1 for value in offset_index)
        if all(value == 0 for value in offset):
            continue
        mask &= values >= np.roll(values, offset, axis=axes)
    return mask


def _nesting_spectrum(
    eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    fermi: float,
    width: float,
    pair_mode: str = "all",
    interpolation_factor: int = 1,
) -> tuple[np.ndarray, tuple[np.ndarray, ...], np.ndarray]:
    """Estimate the FS-JDOS with a periodic FFT cross-correlation.

    For ``pair_mode='all'`` this is

    ``N(q) = sum_(n,m,k) delta_width(E_n(k)-EF) * delta_width(E_m(k+q)-EF)``.

    ``intraband`` retains only ``n=m`` and ``interband`` retains ``n!=m``.
    These two channels use the same sorted eigenvalue index at each k-point;
    they are not connectivity-tracked bands through avoided crossings.
    The result is normalized by its q=0 value when possible, so a nonzero
    peak reports its relative strength rather than an arbitrary FFT scale.
    The q-grid is expressed in fractional reciprocal-lattice coordinates;
    Cartesian q-vectors can be obtained by multiplying by the reciprocal basis.
    """

    if pair_mode not in ("all", "intraband", "interband"):
        raise ValueError("nesting pair mode must be 'all', 'intraband', or 'interband'")
    if interpolation_factor < 1:
        raise ValueError("nesting interpolation factor must be at least one")
    _warn_if_nesting_width_is_narrow(eigenvalues, k_points, dimensions, fermi, width)
    eigenvalues, k_points, dimensions = _spectral_refine_nesting_mesh(
        eigenvalues,
        k_points,
        dimensions,
        interpolation_factor,
    )
    bands, weights, axes = _nesting_weights(eigenvalues, k_points, dimensions, fermi, width)
    nesting = _nesting_from_weights(weights, dimensions, pair_mode)
    return _normalize_nesting(nesting), _nesting_q_axes(axes), bands


def _cross_spin_nesting_spectrum(
    spin_eigenvalues: np.ndarray,
    k_points: np.ndarray,
    dimensions: tuple[int, int, int],
    fermi: float,
    width: float,
    interpolation_factor: int = 1,
) -> tuple[np.ndarray, tuple[np.ndarray, ...], np.ndarray, np.ndarray]:
    """Estimate the ordered up-to-down cross-spin FS-JDOS."""

    spin_eigenvalues = np.asarray(spin_eigenvalues, dtype=np.float64)
    if spin_eigenvalues.ndim != 3 or spin_eigenvalues.shape[0] != 2:
        raise ValueError("cross-spin nesting requires an up/down eigenvalue payload")
    if spin_eigenvalues.shape[2] != k_points.shape[1]:
        raise ValueError("spin-resolved eigenvalues and k-point payload have different sizes")
    _warn_if_nesting_width_is_narrow(
        np.concatenate((spin_eigenvalues[0], spin_eigenvalues[1]), axis=0),
        k_points,
        dimensions,
        fermi,
        width,
    )
    up_values, refined_k_points, refined_dimensions = _spectral_refine_nesting_mesh(
        spin_eigenvalues[0], k_points, dimensions, interpolation_factor
    )
    down_values, _, _ = _spectral_refine_nesting_mesh(
        spin_eigenvalues[1], k_points, dimensions, interpolation_factor
    )
    up_bands, up_weights, axes = _nesting_weights(
        up_values, refined_k_points, refined_dimensions, fermi, width
    )
    down_bands, down_weights, _ = _nesting_weights(
        down_values, refined_k_points, refined_dimensions, fermi, width
    )
    up_transform = _sum_rfft_weights(up_weights)
    down_transform = _sum_rfft_weights(down_weights)
    nesting = np.fft.irfftn(
        np.conj(up_transform) * down_transform,
        s=refined_dimensions,
        axes=(0, 1, 2),
    ).real
    return (
        _normalize_nesting(nesting),
        _nesting_q_axes(axes),
        up_bands,
        down_bands,
    )


def _top_nesting_vectors(
    nesting: np.ndarray,
    q_axes: tuple[np.ndarray, ...],
    count: int,
    ranking: np.ndarray | None = None,
    deduplicate_opposites: bool = True,
) -> list[tuple[np.ndarray, float]]:
    """Return strongest local-maxima q-vectors ranked by an optional map."""

    if count < 1:
        return []
    ranking = np.asarray(nesting if ranking is None else ranking, dtype=np.float64)
    if ranking.shape != nesting.shape:
        raise ValueError("nesting ranking map and nesting map have different shapes")
    local_maxima = _periodic_local_maximum_mask(ranking)
    candidates = np.flatnonzero(local_maxima.ravel())
    candidates = candidates[np.argsort(ranking.ravel()[candidates])[::-1]]
    results: list[tuple[np.ndarray, float]] = []
    for flat_index in candidates:
        index = np.unravel_index(flat_index, nesting.shape)
        rank = float(ranking[index])
        if not np.isfinite(rank) or rank <= 0.0:
            break
        vector = np.array([q_axes[axis][index[axis]] for axis in range(3)], dtype=np.float64)
        if np.allclose(vector, 0.0, rtol=0.0, atol=1.0e-12):
            continue
        if deduplicate_opposites and any(
            np.allclose(vector, -previous, rtol=0.0, atol=1.0e-12) for previous, _ in results
        ):
            continue
        results.append((vector, float(nesting[index])))
        if len(results) == count:
            break
    return results


def _nesting_sheet_data(
    data: dict[str, Any],
    sheet_mode: str,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Select the eigenvalue sheets used by the nesting calculation."""

    if sheet_mode == "combined":
        eigenvalues = np.asarray(data["eigenvalues"], dtype=np.float64)
        labels = tuple(str(index + 1) for index in range(eigenvalues.shape[0]))
        return eigenvalues, labels
    if sheet_mode not in ("both", "up", "down"):
        raise ValueError("nesting sheet mode must be 'combined', 'both', 'up', or 'down'")

    spin_eigenvalues = data.get("spin_eigenvalues")
    if spin_eigenvalues is None:
        raise ValueError(
            "--nesting-sheets up/down/both requires a spin-resolved eigenpair export"
        )
    eigenvalues = np.asarray(spin_eigenvalues, dtype=np.float64)
    if eigenvalues.ndim != 3 or eigenvalues.shape[0] != 2:
        raise ValueError("spin-resolved nesting requires an up/down eigenvalue payload")
    k_points = np.asarray(data["k_points"])
    if eigenvalues.shape[2] != k_points.shape[1]:
        raise ValueError("spin-resolved eigenvalues and k-point payload have different sizes")

    sectors = {"up": (0,), "down": (1,), "both": (0, 1)}[sheet_mode]
    selected = np.concatenate([eigenvalues[sector] for sector in sectors], axis=0)
    labels = tuple(
        f"{SPINS[sector]}:{band + 1}"
        for sector in sectors
        for band in range(eigenvalues.shape[1])
    )
    return selected, labels


def _save_nesting_map(
    path: Path,
    nesting: np.ndarray,
    q_axes: tuple[np.ndarray, ...],
    bands: np.ndarray,
    band_labels: tuple[str, ...],
    fermi: float,
    width: float,
    pair_mode: str,
    sheet_mode: str,
    reciprocal_vectors: np.ndarray | None,
    angular_baseline: np.ndarray | None = None,
    angular_ratio: np.ndarray | None = None,
    interpolation_factor: int = 1,
) -> None:
    """Save the normalized nesting map and its fractional q-grid."""

    q_fractional = np.stack(np.meshgrid(*q_axes, indexing="ij"), axis=0)
    payload: dict[str, np.ndarray] = {
        "nesting": np.asarray(nesting, dtype=np.float64),
        "q_fractional": q_fractional,
        "q1_fractional": np.asarray(q_axes[0], dtype=np.float64),
        "q2_fractional": np.asarray(q_axes[1], dtype=np.float64),
        "q3_fractional": np.asarray(q_axes[2], dtype=np.float64),
        "crossing_bands": np.asarray(bands, dtype=np.int64),
        "contributing_bands": np.asarray(bands, dtype=np.int64),
        "crossing_band_labels": np.asarray([band_labels[band] for band in bands]),
        "contributing_band_labels": np.asarray([band_labels[band] for band in bands]),
        "fermi_level": np.asarray(fermi, dtype=np.float64),
        "gaussian_width": np.asarray(width, dtype=np.float64),
        "band_selection_margin": np.asarray(5.0 * width, dtype=np.float64),
        "pair_mode": np.asarray(pair_mode),
        "sheet_mode": np.asarray(sheet_mode),
        "interpolation_factor": np.asarray(interpolation_factor, dtype=np.int64),
    }
    if reciprocal_vectors is not None:
        payload["reciprocal_vectors"] = np.asarray(reciprocal_vectors, dtype=np.float64)
        payload["q_cartesian"] = _cartesian_coordinates(q_fractional, reciprocal_vectors)
    if angular_baseline is not None:
        payload["angular_baseline"] = np.asarray(angular_baseline, dtype=np.float64)
    if angular_ratio is not None:
        payload["angular_ratio"] = np.asarray(angular_ratio, dtype=np.float64)
    np.savez_compressed(path, **payload)


def _radial_nesting_average(
    nesting: np.ndarray,
    q_axes: tuple[np.ndarray, ...],
    reciprocal_vectors: np.ndarray,
    bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return shell-averaged nesting versus Cartesian ``|q|``."""

    if bins < 1:
        raise ValueError("--nesting-radial-bins must be positive")
    magnitudes = _q_magnitude_grid(q_axes, reciprocal_vectors).ravel()
    values = np.asarray(nesting, dtype=np.float64).ravel()
    if magnitudes.size != values.size:
        raise ValueError("nesting map and q-grid have different sizes")
    if magnitudes.size == 0:
        return np.empty(0), np.empty(0), np.empty(0, dtype=np.int64)

    maximum = float(np.max(magnitudes))
    if maximum <= 0.0:
        return (
            np.array([0.0]),
            np.array([float(np.mean(values))]),
            np.array([values.size], dtype=np.int64),
        )

    edges = np.linspace(0.0, maximum, bins + 1)
    indices = np.floor(magnitudes / maximum * bins).astype(np.int64)
    indices = np.clip(indices, 0, bins - 1)
    sums = np.bincount(indices, weights=values, minlength=bins)
    counts = np.bincount(indices, minlength=bins)
    radii = np.bincount(indices, weights=magnitudes, minlength=bins)
    occupied = counts > 0
    return radii[occupied] / counts[occupied], sums[occupied] / counts[occupied], counts[occupied]


def _save_radial_nesting_plot(
    path: Path,
    radii: np.ndarray,
    values: np.ndarray,
    pair_mode: str,
    sheet_mode: str,
    fermi: float,
    width: float,
) -> None:
    """Save a high-resolution shell-averaged nesting plot."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "--nesting-radial-output needs Matplotlib; install it with "
            "'python -m pip install matplotlib'"
        ) from exc

    figure, axis = plt.subplots(figsize=(8.0, 5.0), dpi=160)
    axis.plot(radii, values, color="#b2182b", linewidth=2.0)
    axis.scatter(radii, values, color="#b2182b", s=12, zorder=3)
    axis.set_xlabel(r"$|\mathbf{q}|$ (Cartesian reciprocal units)")
    axis.set_ylabel(r"$S(|\mathbf{q}|)=\langle N(\mathbf{q})\rangle_{|\mathbf{q}|}$")
    axis.set_title(f"Radially averaged FS nesting ({pair_mode}, {sheet_mode})")
    axis.set_xlim(left=0.0)
    axis.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
    figure.text(0.99, 0.01, f"EF={fermi:.6g} Ry, Gaussian width={width:.6g} Ry", ha="right", fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def _print_nesting_report(args: argparse.Namespace, data: dict[str, Any]) -> None:
    """Print strongest approximate nesting vectors for the current EF."""

    metadata = data["metadata"]
    k_points = np.asarray(data["k_points"], dtype=np.float64)
    dimensions = _mesh_dimensions(metadata, k_points.shape[1])
    try:
        fermi = float(metadata["fermi_level"]) if args.fermi is None else float(args.fermi)
    except KeyError as exc:
        raise ValueError("metadata is missing 'fermi_level'") from exc

    if args.nesting_pairs == "cross-spin":
        if args.nesting_sheets != "both":
            raise ValueError("--nesting-pairs cross-spin requires --nesting-sheets both")
        spin_eigenvalues = data.get("spin_eigenvalues")
        if spin_eigenvalues is None:
            raise ValueError("--nesting-pairs cross-spin requires a spin-resolved eigenpair export")
        spin_eigenvalues = np.asarray(spin_eigenvalues, dtype=np.float64)
        nesting, q_axes, up_bands, down_bands = _cross_spin_nesting_spectrum(
            spin_eigenvalues,
            k_points,
            dimensions,
            fermi,
            args.nesting_width,
            args.nesting_interpolation,
        )
        bands = np.concatenate((up_bands, down_bands + spin_eigenvalues.shape[1]))
        band_labels = tuple(
            f"{SPINS[spin]}:{band + 1}"
            for spin in range(2)
            for band in range(spin_eigenvalues.shape[1])
        )
    else:
        eigenvalues, band_labels = _nesting_sheet_data(data, args.nesting_sheets)
        nesting, q_axes, bands = _nesting_spectrum(
            eigenvalues,
            k_points,
            dimensions,
            fermi,
            args.nesting_width,
            args.nesting_pairs,
            args.nesting_interpolation,
        )
    if args.reciprocal_vectors is not None:
        reciprocal_vectors = np.asarray(args.reciprocal_vectors, dtype=np.float64).reshape((3, 3)).T
    else:
        reciprocal_vectors = _metadata_reciprocal_vectors(metadata)
    if reciprocal_vectors is None:
        warnings.warn(
            "nesting peak ranking is falling back to fractional q shells because no "
            "Cartesian reciprocal vectors are available",
            RuntimeWarning,
            stacklevel=2,
        )
        ranking = nesting
    else:
        baseline = _angular_nesting_baseline(nesting, q_axes, reciprocal_vectors)
        ranking = nesting / np.maximum(baseline, np.finfo(np.float64).tiny)
    top = _top_nesting_vectors(
        nesting,
        q_axes,
        args.nesting_top,
        ranking=ranking,
        deduplicate_opposites=args.nesting_pairs != "cross-spin",
    )

    print(
        f"Nesting FS-JDOS: pairs={args.nesting_pairs}; sheets={args.nesting_sheets}; "
        f"interpolation={args.nesting_interpolation}x; "
        f"Gaussian width={args.nesting_width:.6g} Ry; "
        f"bands contributing near EF={','.join(band_labels[band] for band in bands) or 'none'}"
    )
    if args.nesting_output is not None:
        _save_nesting_map(
            args.nesting_output,
            nesting,
            q_axes,
            bands,
            band_labels,
            fermi,
            args.nesting_width,
            args.nesting_pairs,
            args.nesting_sheets,
            reciprocal_vectors,
            angular_baseline=baseline if reciprocal_vectors is not None else None,
            angular_ratio=ranking if reciprocal_vectors is not None else None,
            interpolation_factor=args.nesting_interpolation,
        )
        print(f"  saved normalized nesting map: {args.nesting_output}")
    if args.nesting_radial_output is not None:
        if reciprocal_vectors is None:
            raise ValueError(
                "radial nesting needs Cartesian reciprocal vectors in the sidecar or "
                "nine values supplied with --reciprocal-vectors"
            )
        radii, radial_values, counts = _radial_nesting_average(
            nesting,
            q_axes,
            reciprocal_vectors,
            args.nesting_radial_bins,
        )
        _save_radial_nesting_plot(
            args.nesting_radial_output,
            radii,
            radial_values,
            args.nesting_pairs,
            args.nesting_sheets,
            fermi,
            args.nesting_width,
        )
        print(
            f"  saved radial nesting plot: {args.nesting_radial_output} "
            f"({radii.size} occupied bins, {int(np.sum(counts))} q-points)"
        )
    if not top:
        print("  no nonzero nesting vectors found")
        return
    for index, (q_fractional, score) in enumerate(top, start=1):
        direct = " ".join(f"{value:.6g}" for value in q_fractional)
        q_index = tuple(int(np.argmin(np.abs(axis - value))) for axis, value in zip(q_axes, q_fractional))
        angular_ratio = float(ranking[q_index])
        line = (
            f"  {index:2d}: q_direct=({direct})  relative_overlap={score:.6g}"
            f"  angular_ratio={angular_ratio:.6g}"
        )
        if reciprocal_vectors is not None:
            q_cartesian = reciprocal_vectors @ q_fractional
            cartesian = " ".join(f"{value:.6g}" for value in q_cartesian)
            line += f"  q_cartesian=({cartesian})"
        print(line)


def _render_window_size(args: argparse.Namespace) -> tuple[int, int]:
    """Return the output size after applying the optional export scale."""

    scale = args.render_scale
    if scale is None:
        scale = 2.0 if args.high_quality else 1.0
    if scale <= 0.0:
        raise ValueError("--render-scale must be positive")
    return tuple(max(1, int(round(size * scale))) for size in args.window_size)


def _configure_high_quality_renderer(plotter: Any, args: argparse.Namespace) -> None:
    """Configure VTK's quality features for a polished still-image export."""

    if args.background is not None:
        plotter.set_background(args.background)
    elif args.high_quality:
        plotter.set_background("white")

    if not args.high_quality:
        return
    if args.ao_radius <= 0.0:
        raise ValueError("--ao-radius must be positive")
    if args.ao_bias < 0.0:
        raise ValueError("--ao-bias must be non-negative")
    if args.ao_kernel_size < 1:
        raise ValueError("--ao-kernel-size must be positive")

    # FXAA is applied after rasterization, while the larger render window
    # supplies additional geometric detail for the final still image.
    try:
        plotter.enable_anti_aliasing("fxaa")
    except (AttributeError, RuntimeError) as exc:
        warnings.warn(f"FXAA is unavailable; continuing without it ({exc})", RuntimeWarning, stacklevel=2)

    try:
        plotter.render_window.SetAlphaBitPlanes(1)
        depth_peeling_supported = plotter.enable_depth_peeling(
            number_of_peels=32,
            occlusion_ratio=0.0,
        )
        if not depth_peeling_supported:
            warnings.warn("depth peeling is unavailable; translucent sheets may sort less accurately", RuntimeWarning, stacklevel=2)
    except (AttributeError, RuntimeError) as exc:
        warnings.warn(f"depth peeling is unavailable ({exc})", RuntimeWarning, stacklevel=2)

    try:
        plotter.enable_ssao(
            radius=args.ao_radius,
            bias=args.ao_bias,
            kernel_size=args.ao_kernel_size,
            blur=True,
        )
    except (AttributeError, RuntimeError) as exc:
        warnings.warn(f"SSAO is unavailable; continuing without it ({exc})", RuntimeWarning, stacklevel=2)


def _projection_values(
    projections: np.ndarray,
    band: int,
    mode: str,
    site: int,
    orbital: int,
    spin: int,
    species: str | None = None,
    site_species: tuple[str, ...] | None = None,
) -> np.ndarray:
    """Select or sum the stored (site, orbital, spin, band, k) projections."""

    selected = projections[:, :, :, band, :]
    if mode == "site":
        return selected[site, :, :, :].sum(axis=(0, 1))
    if mode == "orbital":
        return selected[:, orbital, :, :].sum(axis=(0, 1))
    if mode == "spin":
        up = selected[:, :, 0, :].sum(axis=(0, 1))
        down = selected[:, :, 1, :].sum(axis=(0, 1))
        total = up + down
        return np.divide(up - down, total, out=np.zeros_like(total), where=total > 1.0e-14)
    if mode == "channel":
        return selected[site, orbital, spin, :]
    if mode == "species":
        if species is None or site_species is None:
            raise ValueError("species colouring needs a selected species and site labels")
        if len(site_species) != selected.shape[0]:
            raise ValueError("site species metadata does not match projection dimensions")
        mask = np.asarray(site_species) == species
        if not np.any(mask):
            raise ValueError(f"unknown species: {species}")
        return selected[mask, :, :, :].sum(axis=(0, 1, 2))
    raise ValueError(f"unknown projection mode: {mode}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", help="export base name, with or without .bin/.meta")
    parser.add_argument(
        "--all-bands",
        action="store_true",
        help="plot every band whose sampled range crosses EF (default: one selected band)",
    )
    parser.add_argument(
        "--spin-sheets",
        choices=("combined", "both", "up", "down"),
        default="combined",
        help="plot the full combined FS or independently solved collinear spin sheets (default: combined)",
    )
    parser.add_argument("--band", type=int, default=None, help="initial band number, 1-based (default: closest to EF)")
    parser.add_argument("--fermi", type=float, default=None, help="override the Fermi level in the metadata, in Ry")
    parser.add_argument(
        "--fermi-window",
        type=float,
        default=0.25,
        help="half-width of the interactive Fermi-level slider in Ry (default: 0.25)",
    )
    parser.add_argument(
        "--color-by",
        choices=("none", "band", "site", "orbital", "spin", "species", "channel"),
        default="none",
        help="initial colouring mode (default: none; press p to toggle projection colouring)",
    )
    parser.add_argument("--site", type=int, default=1, help="initial site number, 1-based (default: 1)")
    parser.add_argument(
        "--orbital",
        choices=ORBITALS,
        default="d",
        help="initial orbital channel for orbital/channel colouring (default: d)",
    )
    parser.add_argument("--spin", choices=SPINS, default="up", help="initial spin channel (default: up)")
    parser.add_argument(
        "--spin-frame",
        choices=("basis", "local"),
        default="basis",
        help="spin projection frame: eigenvector basis up/down (default) or exported local moment axis",
    )
    parser.add_argument(
        "--species",
        default=None,
        help="initial species for species colouring (default: first species in the export)",
    )
    parser.add_argument(
        "--coordinate-system",
        choices=("fractional", "cartesian"),
        default="cartesian",
        help="plot Cartesian reciprocal coordinates (default) or fractional/direct coordinates",
    )
    parser.add_argument(
        "--reciprocal-vectors",
        type=float,
        nargs=9,
        metavar=("B1X", "B1Y", "B1Z", "B2X", "B2Y", "B2Z", "B3X", "B3Y", "B3Z"),
        help="optional Cartesian reciprocal vectors, used with --coordinate-system cartesian",
    )
    seam_group = parser.add_mutually_exclusive_group()
    seam_group.add_argument(
        "--close-bz-seams",
        dest="close_bz_seams",
        action="store_true",
        help="close equivalent Brillouin-zone faces before contouring (default)",
    )
    seam_group.add_argument(
        "--no-close-bz-seams",
        dest="close_bz_seams",
        action="store_false",
        help="leave the sampled BZ faces open for diagnostic plots",
    )
    parser.set_defaults(close_bz_seams=True)
    # Backwards-compatible spelling; the operation is a visualization seam
    # closure, not a switch that enables physical periodic boundary conditions.
    parser.add_argument("--periodic", dest="close_bz_seams", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--nesting",
        action="store_true",
        help="print an FFT-based Fermi-surface joint-density-of-states nesting analysis",
    )
    parser.add_argument(
        "--nesting-only",
        action="store_true",
        help="print/save the nesting analysis and exit without opening the PyVista viewer",
    )
    parser.add_argument(
        "--nesting-pairs",
        choices=("all", "intraband", "interband", "cross-spin"),
        default="all",
        help="pair channel: all, same sorted index, different sorted index, or up-to-down cross-spin",
    )
    parser.add_argument(
        "--nesting-sheets",
        choices=("combined", "both", "up", "down"),
        default="combined",
        help="eigenvalue sheets for --nesting: combined or spin sectors (default: combined)",
    )
    parser.add_argument(
        "--nesting-output",
        type=Path,
        default=None,
        metavar="FILE",
        help="save the normalized nesting map and q-grids as a compressed .npz file",
    )
    parser.add_argument(
        "--nesting-radial-output",
        type=Path,
        default=None,
        metavar="IMAGE",
        help="save a shell-averaged S(|q|) nesting plot (requires Cartesian reciprocal vectors)",
    )
    parser.add_argument(
        "--nesting-radial-bins",
        type=int,
        default=64,
        metavar="N",
        help="number of radial bins for --nesting-radial-output (default: 64)",
    )
    parser.add_argument(
        "--nesting-interpolation",
        type=int,
        default=1,
        metavar="FACTOR",
        help="periodically interpolate energies by this factor before nesting (default: 1)",
    )
    parser.add_argument(
        "--nesting-width",
        type=float,
        default=0.01,
        metavar="RY",
        help="Gaussian energy width for --nesting (default: 0.01 Ry)",
    )
    parser.add_argument(
        "--nesting-top",
        type=int,
        default=10,
        metavar="N",
        help="number of non-equivalent strongest nesting vectors to print (default: 10)",
    )
    parser.add_argument(
        "--high-quality",
        action="store_true",
        help="use a polished still-image renderer: 2x resolution, FXAA, SSAO, and depth peeling",
    )
    parser.add_argument(
        "--render-scale",
        type=float,
        default=None,
        metavar="FACTOR",
        help="multiply --window-size for output (default: 1, or 2 with --high-quality)",
    )
    parser.add_argument(
        "--background",
        default=None,
        help="override the background colour; high-quality mode otherwise uses solid white",
    )
    parser.add_argument(
        "--ao-radius",
        type=float,
        default=0.6,
        metavar="RADIUS",
        help="SSAO neighbour radius in high-quality mode (default: 0.6)",
    )
    parser.add_argument(
        "--ao-bias",
        type=float,
        default=0.012,
        metavar="BIAS",
        help="SSAO depth bias in high-quality mode (default: 0.012)",
    )
    parser.add_argument(
        "--ao-kernel-size",
        type=int,
        default=256,
        metavar="N",
        help="SSAO sample count in high-quality mode (default: 256)",
    )
    parser.add_argument(
        "--smooth-interpolation",
        type=int,
        default=1,
        metavar="FACTOR",
        help="refine each mesh cell by this factor with cubic interpolation (default: disabled)",
    )
    parser.add_argument("--screenshot", type=Path, help="save a screenshot instead of opening an interactive window")
    parser.add_argument("--window-size", type=int, nargs=2, default=(1200, 900), metavar=("WIDTH", "HEIGHT"))
    return parser


def _build_plotter(args: argparse.Namespace, data: dict[str, Any]):
    try:
        import pyvista as pv
    except ImportError as exc:
        raise RuntimeError("PyVista is required; install it with 'python -m pip install pyvista'") from exc

    metadata = data["metadata"]
    k_points = np.asarray(data["k_points"], dtype=np.float64)
    eigenvalues = np.asarray(data["eigenvalues"], dtype=np.float64)
    eigenvectors = np.asarray(data["eigenvectors"])
    spin_eigenvalues = data.get("spin_eigenvalues")
    spin_eigenvectors = data.get("spin_eigenvectors")
    if args.spin_sheets == "combined":
        selected_sheets = (0,)
        sheet_labels = ("combined",)
        eigenvalues_by_sheet = (eigenvalues,)
        eigenvectors_by_sheet = (eigenvectors,)
    else:
        if spin_eigenvalues is None or spin_eigenvectors is None:
            raise ValueError(
                "spin-resolved sheets were requested, but the export has no spin payload; "
                "rerun with write_spin_resolved_eigenpairs=.true."
            )
        spin_eigenvalues = np.asarray(spin_eigenvalues, dtype=np.float64)
        spin_eigenvectors = np.asarray(spin_eigenvectors)
        if spin_eigenvalues.ndim != 3 or spin_eigenvalues.shape[0] != 2:
            raise ValueError("spin-resolved eigenvalues must have shape (2,n_bands,n_kpoints)")
        if spin_eigenvectors.ndim != 4 or spin_eigenvectors.shape[0] != 2:
            raise ValueError("spin-resolved eigenvectors must have shape (2,n_basis,n_bands,n_kpoints)")
        selected_sheets = {"up": (0,), "down": (1,), "both": (0, 1)}[args.spin_sheets]
        sheet_labels = tuple(SPINS[index] for index in selected_sheets)
        eigenvalues_by_sheet = tuple(spin_eigenvalues[index] for index in selected_sheets)
        eigenvectors_by_sheet = tuple(spin_eigenvectors[index] for index in selected_sheets)
    dimensions = _mesh_dimensions(metadata, k_points.shape[1])
    axes, grid_indices = _regular_grid_index_map(k_points, dimensions)
    for sheet_values in eigenvalues_by_sheet:
        if sheet_values.shape[1] != k_points.shape[1]:
            raise ValueError("eigenvalue and k-point counts differ")
    n_bands = eigenvalues_by_sheet[0].shape[0]
    if any(sheet_values.shape[0] != n_bands for sheet_values in eigenvalues_by_sheet):
        raise ValueError("spin-resolved sheets have different band counts")

    try:
        fermi = float(metadata["fermi_level"]) if args.fermi is None else float(args.fermi)
    except KeyError as exc:
        raise ValueError("metadata is missing 'fermi_level'") from exc
    if args.fermi_window <= 0.0:
        raise ValueError("--fermi-window must be positive")

    if args.band is None:
        distances = np.concatenate(
            [np.where(np.isfinite(values), np.abs(values - fermi), np.inf) for values in eigenvalues_by_sheet],
            axis=0,
        )
        band = int(np.argmin(np.min(distances, axis=1)) % n_bands)
    else:
        if args.band < 1 or args.band > n_bands:
            raise ValueError(f"--band must be between 1 and {n_bands}")
        band = args.band - 1

    exported_projections = data.get("projection_weights")
    if exported_projections is not None:
        exported_projections = np.asarray(exported_projections, dtype=np.float64)
    basis_projections = _basis_projection_weights(eigenvectors, metadata)
    projections = basis_projections if basis_projections is not None else exported_projections
    if args.spin_frame == "local" and exported_projections is None:
        raise ValueError("--spin-frame local needs exported projection_weights")
    if args.spin_sheets != "combined" and args.spin_frame == "local":
        raise ValueError("--spin-frame local is unavailable for independently solved global spin sheets")

    if args.spin_sheets == "combined":
        projections_by_sheet = (projections,)
        spin_projections_by_sheet = (exported_projections if args.spin_frame == "local" else projections,)
    else:
        projections_by_sheet = tuple(
            _spin_sector_projection_weights(vector, metadata, spin_index)
            for vector, spin_index in zip(eigenvectors_by_sheet, selected_sheets)
        )
        if any(projection is None for projection in projections_by_sheet):
            projections_by_sheet = tuple(None for _ in projections_by_sheet)
        projections = projections_by_sheet[0]
        spin_projections_by_sheet = projections_by_sheet
    if projections is not None:
        projections = np.asarray(projections, dtype=np.float64)
        nsites = projections.shape[0]
        if args.site < 1 or args.site > nsites:
            raise ValueError(f"--site must be between 1 and {nsites}")
    elif args.color_by in ("site", "orbital", "spin", "species", "channel"):
        raise ValueError("projection colouring was requested, but no compatible eigenvector projections are available")
    site = args.site - 1
    orbital = ORBITALS.index(args.orbital)
    spin = SPINS.index(args.spin)
    site_species = _site_species(metadata, nsites) if projections is not None else ()
    species_options = _unique_species(site_species)
    if args.color_by == "species":
        if not species_options:
            raise ValueError("species colouring needs site metadata and a compatible eigenvector basis")
        species = args.species or species_options[0]
        if species not in species_options:
            raise ValueError(f"--species must be one of: {', '.join(species_options)}")
    else:
        species = args.species if args.species in species_options else (species_options[0] if species_options else "")

    reciprocal_vectors = None
    if args.coordinate_system == "cartesian":
        if args.reciprocal_vectors is not None:
            reciprocal_vectors = np.asarray(args.reciprocal_vectors, dtype=np.float64).reshape((3, 3)).T
        else:
            reciprocal_vectors = _metadata_reciprocal_vectors(metadata)
        if reciprocal_vectors is None:
            raise ValueError(
                "Cartesian plotting needs reciprocal vectors in the sidecar or "
                "nine values supplied with --reciprocal-vectors"
            )

    plot_axes = axes
    if args.close_bz_seams:
        # The scalar arrays are extended together with the coordinate axes in
        # update_surface; this flag only changes the initial grid dimensions.
        plot_axes = tuple(
            np.concatenate((axis, [axis[0] + axis.size * (axis[1] - axis[0] if axis.size > 1 else 1.0)]))
            for axis in axes
        )
    if args.smooth_interpolation < 1:
        raise ValueError("--smooth-interpolation must be at least one")
    plot_axes = _refine_axes(plot_axes, args.smooth_interpolation)

    fractional_mesh = np.stack(np.meshgrid(*plot_axes, indexing="ij"), axis=0)
    if reciprocal_vectors is None:
        coordinate_mesh = fractional_mesh
        coordinate_label = "fractional reciprocal coordinates"
    else:
        coordinate_mesh = _cartesian_coordinates(fractional_mesh, reciprocal_vectors)
        coordinate_label = "Cartesian reciprocal coordinates"
    grid = pv.StructuredGrid(coordinate_mesh[0], coordinate_mesh[1], coordinate_mesh[2])
    window_size = _render_window_size(args)
    clean_export = args.high_quality and args.screenshot is not None

    state: dict[str, Any] = {
        "actors": [],
        "scalar_bar_title": None,
        "status": None,
        "band": band,
        "fermi": fermi,
        "opacity": 0.85,
        "site": site,
        "orbital": orbital,
        "spin": spin,
        "species": species,
        "color_mode": args.color_by if args.color_by != "none" else "site",
        "color_enabled": args.color_by != "none",
        "use_projection": args.color_by in ("site", "orbital", "spin", "species", "channel"),
        "bands": [],
        "axes_visible": True,
        "labels_visible": True,
        "grid_visible": True,
    }

    axes_widget = None
    cube_axes_actor = None
    coordinate_annotation = None
    help_annotation = None
    axis_titles: dict[str, str] = {}

    def apply_axes_visibility() -> None:
        if axes_widget is not None:
            axes_widget.SetEnabled(int(state["axes_visible"]))
        if cube_axes_actor is not None:
            visible = int(state["axes_visible"])
            for axis in "XYZ":
                getattr(cube_axes_actor, f"Set{axis}AxisVisibility")(visible)

    def apply_grid_visibility() -> None:
        if cube_axes_actor is None:
            return
        visible = int(state["grid_visible"])
        for axis in "XYZ":
            getattr(cube_axes_actor, f"SetDraw{axis}Gridlines")(visible)

    def apply_label_visibility() -> None:
        visible = bool(state["labels_visible"])
        for actor in (state["status"], coordinate_annotation, help_annotation):
            if actor is not None:
                actor.SetVisibility(int(visible))
        if cube_axes_actor is not None:
            for axis in "XYZ":
                getattr(cube_axes_actor, f"Set{axis}Title")(axis_titles[axis] if visible else "")
                getattr(cube_axes_actor, f"Set{axis}AxisLabelVisibility")(int(visible))
                getattr(cube_axes_actor, f"Set{axis}AxisTickVisibility")(int(visible))
        for actor in getattr(plotter, "_scalar_bar_actors", {}).values():
            actor.SetVisibility(int(visible))

    def active_color_mode() -> str:
        if not state["color_enabled"]:
            return "none"
        if state["color_mode"] == "band":
            return "band"
        if state["use_projection"] and projections is not None:
            return state["color_mode"]
        return "none"

    def update_status(plotter) -> None:
        if clean_export:
            return
        if state["status"] is not None:
            plotter.remove_actor(state["status"], render=False)
        if args.all_bands:
            if state["bands"]:
                band_text = "bands " + ",".join(str(value + 1) for value in state["bands"])
            else:
                band_text = "no band crosses EF"
        else:
            band_text = f"band {state['band'] + 1}/{n_bands}"
        if args.spin_sheets != "combined":
            band_text = f"{','.join(sheet_labels)} sheets; {band_text}"

        mode = active_color_mode()
        if mode == "none":
            color_text = "off"
        elif mode == "band":
            color_text = "band index"
        else:
            color_text = mode
            if mode in ("site", "channel"):
                color_text += f", site {state['site'] + 1}"
            if mode in ("orbital", "channel"):
                color_text += f", {ORBITALS[state['orbital']]}"
            if mode == "channel":
                color_text += f", {SPINS[state['spin']]}"
            elif mode == "spin":
                color_text += ", up/down polarization"
            elif mode == "species":
                color_text += f", {state['species']}"
        text = f"{band_text}   EF={state['fermi']:.8f} Ry   {color_text}"
        state["status"] = plotter.add_text(text, position="upper_left", font_size=11, name="fs_status")
        apply_label_visibility()

    def update_surface(plotter) -> None:
        if args.all_bands:
            bands = sorted(
                {
                    int(band)
                    for sheet_values in eigenvalues_by_sheet
                    for band in _crossing_bands(sheet_values, state["fermi"])
                }
            )
        else:
            bands = [state["band"]]
        state["bands"] = bands

        for actor in state["actors"]:
            plotter.remove_actor(actor, render=False)
        state["actors"] = []
        if state["scalar_bar_title"] is not None:
            plotter.remove_scalar_bar(state["scalar_bar_title"])
            state["scalar_bar_title"] = None

        color_mode = active_color_mode()
        scalar_bar_title = {
            "band": "band index",
            "spin": "spin polarization (up-down)",
            "species": f"{state['species']} weight",
        }.get(color_mode, "projection")
        scalar_bar_args = {
            "title": scalar_bar_title,
            "vertical": True,
            "position_x": 0.88,
            "position_y": 0.25,
            "height": 0.4,
            "width": 0.08,
        }
        if args.high_quality:
            scalar_bar_args["color"] = "black"
        band_ranks = {band: rank for rank, band in enumerate(bands, start=1)}
        band_lookup_table = None
        if color_mode == "band":
            number_of_bands = max(1, len(bands))
            # The annotations carry the original (possibly non-contiguous)
            # band numbers; suppress the auxiliary rank ticks on the bar.
            scalar_bar_args["n_labels"] = 0
            band_lookup_table = pv.LookupTable(
                cmap="tab20",
                n_values=number_of_bands,
                scalar_range=(0.5, number_of_bands + 0.5),
                annotations={float(rank): str(current_band + 1) for current_band, rank in band_ranks.items()},
            )

        actor_index = 0
        for sheet_index, (sheet_values, sheet_projections) in enumerate(
            zip(eigenvalues_by_sheet, projections_by_sheet)
        ):
            for current_band in bands:
                band_values = _values_on_grid(sheet_values[current_band], k_points, dimensions, grid_indices)
                scalar_values = band_values
                if args.close_bz_seams:
                    scalar_values, _ = _periodic_extend(band_values, axes)
                scalar_values = _smooth_grid(scalar_values, args.smooth_interpolation)

                projection_values = None
                if color_mode not in ("none", "band"):
                    projection_source = (
                        spin_projections_by_sheet[sheet_index]
                        if color_mode == "spin"
                        else sheet_projections
                    )
                    if projection_source is None:
                        raise ValueError("the selected spin sheet has no compatible projection basis")
                    raw_values = _projection_values(
                        projection_source,
                        current_band,
                        color_mode,
                        state["site"],
                        state["orbital"],
                        state["spin"],
                        state["species"],
                        site_species,
                    )
                    projection_values = _values_on_grid(raw_values, k_points, dimensions, grid_indices)
                    if args.close_bz_seams:
                        projection_values, _ = _periodic_extend(projection_values, axes)
                    projection_values = _smooth_grid(projection_values, args.smooth_interpolation)
                    if color_mode == "spin":
                        projection_values = np.clip(projection_values, -1.0, 1.0)
                    else:
                        projection_values = np.clip(projection_values, 0.0, 1.0)

                # Grid coordinates are fixed, but scalars change with every
                # widget. Reuse the same structured-grid topology for each band.
                grid.point_data["energy"] = scalar_values.ravel(order="F")
                if color_mode == "band":
                    grid.point_data["band_index"] = np.full(
                        scalar_values.shape, float(band_ranks[current_band]), dtype=np.float64
                    ).ravel(order="F")
                elif projection_values is not None:
                    grid.point_data["projection"] = projection_values.ravel(order="F")
                surface = grid.contour(isosurfaces=[state["fermi"]], scalars="energy")
                if not surface.n_points:
                    continue

                mesh_kwargs: dict[str, Any] = {
                    "opacity": state["opacity"],
                    "smooth_shading": True,
                    "name": f"fermi_surface_{actor_index}",
                }
                actor_index += 1
                if args.high_quality:
                    mesh_kwargs.update(
                        ambient=0.24,
                        diffuse=0.76,
                        specular=0.12,
                        specular_power=20.0,
                    )
                if color_mode == "band":
                    mesh_kwargs.update(
                        scalars="band_index",
                        cmap=band_lookup_table,
                    )
                elif projection_values is not None:
                    mesh_kwargs.update(
                        scalars="projection",
                        cmap="Reds" if color_mode == "spin" else "viridis",
                        clim=(-1.0, 1.0) if color_mode == "spin" else (0.0, 1.0),
                    )
                else:
                    mesh_kwargs["color"] = "cornflowerblue"

                if color_mode != "none" and not state["actors"]:
                    mesh_kwargs["scalar_bar_args"] = scalar_bar_args
                    state["scalar_bar_title"] = scalar_bar_title
                elif color_mode != "none":
                    mesh_kwargs["show_scalar_bar"] = False
                state["actors"].append(plotter.add_mesh(surface, **mesh_kwargs))
        update_status(plotter)

    plotter = pv.Plotter(
        window_size=window_size,
        off_screen=args.screenshot is not None,
        title="RSLMTO Fermi surface",
    )
    _configure_high_quality_renderer(plotter, args)
    update_surface(plotter)
    coordinate_extent = np.ptp(coordinate_mesh.reshape(3, -1), axis=1)
    # Newer PyVista releases expose an explicit box-aspect setter. Keep the
    # call optional because older releases preserve world-coordinate aspect
    # directly through the VTK renderer and do not provide this method.
    set_box_aspect = getattr(plotter, "set_box_aspect", None)
    if callable(set_box_aspect):
        set_box_aspect(tuple(np.maximum(coordinate_extent, np.finfo(float).eps)))
    plotter.add_axes(
        color="black" if args.high_quality else None,
        x_color="black" if args.high_quality else None,
        y_color="black" if args.high_quality else None,
        z_color="black" if args.high_quality else None,
    )
    axes_widget = getattr(plotter.renderer, "axes_widget", None)
    cube_axes_actor = plotter.show_bounds(
        grid="back",
        location="outer",
        ticks="outside",
        xtitle="kx",
        ytitle="ky",
        ztitle="kz",
        color="black" if args.high_quality else None,
    )
    axis_titles = {
        "X": cube_axes_actor.GetXTitle(),
        "Y": cube_axes_actor.GetYTitle(),
        "Z": cube_axes_actor.GetZTitle(),
    }
    coordinate_annotation = plotter.add_text(
        coordinate_label,
        position="upper_right",
        font_size=10,
        color="black" if args.high_quality else None,
        name="fs_coordinates",
    )
    apply_axes_visibility()
    apply_grid_visibility()
    apply_label_visibility()
    if not clean_export:
        help_text = "EF/opacity sliders below" if args.all_bands else "Band/EF/opacity sliders below"
        if args.all_bands:
            help_text += "; all bands crossing EF are shown"
        if projections is not None:
            help_text += "; P toggles projection; 1=site 2=orbital 3=spin 4=channel 5=species"
        help_annotation = plotter.add_text(
            help_text,
            position="lower_left",
            font_size=9,
            color="black" if args.high_quality else None,
            name="fs_help",
        )

        def toggle_axes(visible: bool) -> None:
            state["axes_visible"] = bool(visible)
            apply_axes_visibility()

        def toggle_labels(visible: bool) -> None:
            state["labels_visible"] = bool(visible)
            apply_label_visibility()

        def toggle_grid(visible: bool) -> None:
            state["grid_visible"] = bool(visible)
            apply_grid_visibility()

        # Keep the controls below the top annotations.  VTK uses a
        # bottom-left pixel origin for button widgets.
        button_y = max(10, window_size[1] - 85)
        for index, (title, callback) in enumerate(
            (("Axes", toggle_axes), ("Labels", toggle_labels), ("Grid", toggle_grid))
        ):
            button_x = 10 + index * 100
            plotter.add_checkbox_button_widget(
                callback,
                value=True,
                position=(button_x, button_y),
                size=24,
                border_size=2,
                color_on="black",
                color_off="lightgray",
                background_color="white",
            )
            plotter.add_text(
                title,
                position=(button_x + 30, button_y + 5),
                font_size=9,
                color="black" if args.high_quality else None,
                name=f"toggle_label_{title.lower()}",
            )

    def band_callback(value: float, widget: Any) -> None:
        band = _snap_integer_slider(value, widget, 1, n_bands) - 1
        if band != state["band"]:
            state["band"] = band
            update_surface(plotter)

    def fermi_callback(value: float) -> None:
        state["fermi"] = float(value)
        update_surface(plotter)

    def opacity_callback(value: float) -> None:
        state["opacity"] = float(value)
        for actor in state["actors"]:
            actor.GetProperty().SetOpacity(state["opacity"])

    if not args.all_bands and not clean_export:
        plotter.add_slider_widget(
            band_callback,
            (1.0, float(n_bands)),
            value=float(band + 1),
            title="Band",
            pointa=(0.03, 0.08),
            pointb=(0.27, 0.08),
            style="modern",
            interaction_event="always",
            pass_widget=True,
        )
    fermi_start = 0.03 if args.all_bands else 0.32
    fermi_end = 0.27 if args.all_bands else 0.56
    opacity_start = 0.32 if args.all_bands else 0.61
    opacity_end = 0.56 if args.all_bands else 0.84
    if not clean_export:
        plotter.add_slider_widget(
            fermi_callback,
            (fermi - args.fermi_window, fermi + args.fermi_window),
            value=fermi,
            title="EF (Ry)",
            pointa=(fermi_start, 0.08),
            pointb=(fermi_end, 0.08),
            style="modern",
            interaction_event="end",
        )
        plotter.add_slider_widget(
            opacity_callback,
            (0.1, 1.0),
            value=state["opacity"],
            title="Opacity",
            pointa=(opacity_start, 0.08),
            pointb=(opacity_end, 0.08),
            style="modern",
            interaction_event="end",
        )

    if projections is not None and not clean_export:
        nsites = projections.shape[0]

        def site_callback(value: float, widget: Any) -> None:
            site = _snap_integer_slider(value, widget, 1, nsites) - 1
            if site != state["site"]:
                state["site"] = site
                update_surface(plotter)

        def orbital_callback(value: float, widget: Any) -> None:
            orbital = _snap_integer_slider(value, widget, 1, 4) - 1
            if orbital != state["orbital"]:
                state["orbital"] = orbital
                update_surface(plotter)

        def spin_callback(value: float, widget: Any) -> None:
            spin = _snap_integer_slider(value, widget, 1, 2) - 1
            if spin != state["spin"]:
                state["spin"] = spin
                update_surface(plotter)

        plotter.add_slider_widget(
            site_callback,
            (1.0, float(nsites)),
            value=float(site + 1),
            title="Site",
            pointa=(0.03, 0.15),
            pointb=(0.27, 0.15),
            style="modern",
            interaction_event="always",
            pass_widget=True,
        )
        plotter.add_slider_widget(
            orbital_callback,
            (1.0, 4.0),
            value=float(orbital + 1),
            title="Orbital",
            pointa=(0.32, 0.15),
            pointb=(0.56, 0.15),
            style="modern",
            interaction_event="always",
            pass_widget=True,
        )
        plotter.add_slider_widget(
            spin_callback,
            (1.0, 2.0),
            value=float(spin + 1),
            title="Spin",
            pointa=(0.61, 0.15),
            pointb=(0.84, 0.15),
            style="modern",
            interaction_event="always",
            pass_widget=True,
        )

        if len(species_options) > 1:

            def species_callback(value: float, widget: Any) -> None:
                species_index = _snap_integer_slider(value, widget, 1, len(species_options)) - 1
                species = species_options[species_index]
                if species != state["species"]:
                    state["species"] = species
                    update_surface(plotter)

            plotter.add_slider_widget(
                species_callback,
                (1.0, float(len(species_options))),
                value=float(species_options.index(state["species"]) + 1),
                title="Species",
                pointa=(0.03, 0.22),
                pointb=(0.27, 0.22),
                style="modern",
                interaction_event="always",
                pass_widget=True,
            )

        def toggle_projection() -> None:
            if state["color_mode"] == "band":
                state["color_enabled"] = not state["color_enabled"]
            else:
                state["color_enabled"] = True
                state["use_projection"] = not state["use_projection"]
            update_surface(plotter)

        plotter.add_key_event("p", toggle_projection)

        def set_projection_mode(mode: str) -> None:
            state["color_mode"] = mode
            state["color_enabled"] = True
            state["use_projection"] = True
            update_surface(plotter)

        for key, mode in (
            ("1", "site"),
            ("2", "orbital"),
            ("3", "spin"),
            ("4", "channel"),
            ("5", "species"),
        ):
            plotter.add_key_event(key, lambda mode=mode: set_projection_mode(mode))

    def set_band_coloring() -> None:
        state["color_mode"] = "band"
        state["color_enabled"] = True
        state["use_projection"] = False
        update_surface(plotter)

    plotter.add_key_event("b", set_band_coloring)

    return plotter, args.screenshot


def main() -> None:
    args = _parser().parse_args()
    try:
        data = read_kspace_eigenpairs(args.base)
        if args.nesting or args.nesting_only or args.nesting_radial_output is not None:
            _print_nesting_report(args, data)
        if args.nesting_only:
            return
        plotter, screenshot = _build_plotter(args, data)
    except (OSError, ValueError, RuntimeError) as exc:
        _parser().error(str(exc))

    if screenshot is not None:
        plotter.show(screenshot=str(screenshot), auto_close=True)
        print(f"saved {screenshot}")
    else:
        plotter.show()


if __name__ == "__main__":
    main()
