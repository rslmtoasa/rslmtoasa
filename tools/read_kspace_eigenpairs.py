#!/usr/bin/env python3
"""Read an RSLMTO dense k-space eigenpair export.

The returned arrays use the same names and shapes listed in the ``.meta``
sidecar.  All k-points are fractional reciprocal-lattice coordinates and all
energies are in Ry.  Eigenvectors are returned as a complex array with shape
``(n_basis, n_bands, n_kpoints)``.
When present, independently solved collinear sectors are returned as
``spin_eigenvalues`` with shape ``(2, n_spin_bands, n_kpoints)`` and
``spin_eigenvectors`` with shape ``(2, n_spin_basis, n_spin_bands, n_kpoints)``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def _read_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        metadata[key.strip()] = value.strip()
    return metadata


def read_kspace_eigenpairs(base: str | Path) -> dict[str, object]:
    """Read ``base + '.bin'`` and ``base + '.meta'`` into NumPy arrays."""

    base_path = Path(base)
    binary_path = base_path.with_suffix(".bin") if base_path.suffix else Path(str(base_path) + ".bin")
    metadata_path = base_path.with_suffix(".meta") if base_path.suffix else Path(str(base_path) + ".meta")
    metadata = _read_metadata(metadata_path)

    nk = int(metadata["n_kpoints"])
    nbands = int(metadata["n_bands"])
    nbasis = int(metadata["n_basis"])
    nsites = int(metadata["n_sites"])
    values = np.fromfile(binary_path, dtype=np.dtype("=f8"))
    cursor = 0

    def take(shape: tuple[int, ...]) -> np.ndarray:
        nonlocal cursor
        count = int(np.prod(shape, dtype=np.int64))
        raw = values[cursor : cursor + count]
        if raw.size != count:
            raise ValueError(f"truncated eigenpair payload in {binary_path}")
        cursor += count
        return raw.reshape(shape, order="F").copy()

    k_points = take((3, nk))
    k_weights = take((nk,))
    eigenvalues = take((nbands, nk))
    eigenvectors_real = take((nbasis, nbands, nk))
    eigenvectors_imag = take((nbasis, nbands, nk))
    result: dict[str, object] = {
        "metadata": metadata,
        "k_points": k_points,
        "k_weights": k_weights,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors_real + 1j * eigenvectors_imag,
    }

    if metadata.get("projection_weights", "disabled") == "enabled":
        result["projection_weights"] = take((nsites, 4, 2, nbands, nk))
    if cursor != values.size:
        raise ValueError(f"unexpected trailing data in {binary_path}")

    # Independently solved collinear spin sectors live in a separate stream so
    # the original payload remains byte-for-byte backwards compatible.
    if metadata.get("spin_resolved_eigenpairs", "disabled") == "enabled":
        try:
            spin_binary_name = metadata["spin_resolved_binary_file"]
            nspin = int(metadata["spin_resolved_n_spin"])
            spin_nbands = int(metadata["spin_resolved_n_bands"])
            spin_nbasis = int(metadata["spin_resolved_n_basis"])
        except (KeyError, ValueError) as exc:
            raise ValueError("incomplete spin-resolved eigenpair metadata") from exc
        if nspin < 1 or spin_nbands < 1 or spin_nbasis < 1:
            raise ValueError("invalid spin-resolved eigenpair dimensions")

        # The Fortran writer records the same path spelling used for the main
        # payload. Deriving it from the resolved ``.bin`` path also handles a
        # base name containing directories (``out/fermi_surface``).
        derived_spin_binary_path = binary_path.with_suffix(".spin.bin")
        recorded_spin_binary_path = Path(spin_binary_name)
        if recorded_spin_binary_path.is_absolute():
            spin_binary_path = recorded_spin_binary_path
        elif derived_spin_binary_path.exists():
            spin_binary_path = derived_spin_binary_path
        else:
            spin_binary_path = metadata_path.parent / recorded_spin_binary_path
            if not spin_binary_path.exists():
                spin_binary_path = recorded_spin_binary_path
        spin_values = np.fromfile(spin_binary_path, dtype=np.dtype("=f8"))
        spin_cursor = 0

        def take_spin(shape: tuple[int, ...]) -> np.ndarray:
            nonlocal spin_cursor
            count = int(np.prod(shape, dtype=np.int64))
            raw = spin_values[spin_cursor : spin_cursor + count]
            if raw.size != count:
                raise ValueError(f"truncated spin-resolved payload in {spin_binary_path}")
            spin_cursor += count
            return raw.reshape(shape, order="F").copy()

        spin_eigenvalues = take_spin((nspin, spin_nbands, nk))
        spin_eigenvectors_real = take_spin((nspin, spin_nbasis, spin_nbands, nk))
        spin_eigenvectors_imag = take_spin((nspin, spin_nbasis, spin_nbands, nk))
        if spin_cursor != spin_values.size:
            raise ValueError(f"unexpected trailing data in {spin_binary_path}")
        result["spin_eigenvalues"] = spin_eigenvalues
        result["spin_eigenvectors"] = spin_eigenvectors_real + 1j * spin_eigenvectors_imag
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", help="export base name, with or without .bin/.meta")
    args = parser.parse_args()
    data = read_kspace_eigenpairs(args.base)
    metadata = data["metadata"]
    print(
        f"Read {metadata['n_kpoints']} k-points, {metadata['n_bands']} bands, "
        f"{metadata['n_basis']} basis coefficients; EF={metadata['fermi_level']} Ry"
    )
    print("arrays:", ", ".join(key for key in data if key != "metadata"))


if __name__ == "__main__":
    main()
