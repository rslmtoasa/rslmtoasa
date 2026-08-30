#!/usr/bin/env python3
"""Reusable matched primitive-operator construction for GBT WP04.

The utility deliberately keeps the two representations separate:

* ``build_primitive_gbt`` applies the endpoint link to the raw directed
  orbital block before the local collinear LMTO factors;
* ``build_explicit_supercell`` uses the same directed records to rotate the
  local endpoint factors into the lab frame and tiles the finite operator.

The records retain three-component cell and physical displacements.  The
production convention is q in Cartesian units of ``2*pi/alat`` and the
ordinary Fourier factor is ``exp(+2*pi*i*k.R)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


COMPLEX = np.complex128
IDENTITY_SPIN = np.eye(2, dtype=COMPLEX)


def su2_frame(theta: float, phi: float) -> np.ndarray:
    """Return the repository's reference spin frame U(theta, phi)."""
    c = np.cos(0.5 * theta)
    s = np.sin(0.5 * theta)
    ephi = np.exp(1.0j * phi)
    return np.array([[c, -np.conj(ephi) * s], [ephi * s, c]], dtype=COMPLEX)


def z_rotation(alpha: float) -> np.ndarray:
    """Return exp(-i alpha sigma_z / 2), independently assembled."""
    return np.diag(np.exp(-0.5j * alpha * np.array([1.0, -1.0]))).astype(COMPLEX)


def spin_lift(orbital: np.ndarray, spin: np.ndarray = IDENTITY_SPIN) -> np.ndarray:
    """Lift an orbital block in the production spin-block basis."""
    return np.kron(spin, orbital).astype(COMPLEX)


def endpoint_factor(w0: np.ndarray, w1: np.ndarray) -> np.ndarray:
    """Build W = diag(w0+w1, w0-w1) in spin-block ordering."""
    w0 = np.asarray(w0, dtype=COMPLEX)
    w1 = np.asarray(w1, dtype=COMPLEX)
    zeros = np.zeros((len(w0), len(w0)), dtype=COMPLEX)
    return np.block([[np.diag(w0 + w1), zeros], [zeros, np.diag(w0 - w1)]])


@dataclass(frozen=True)
class DirectedBond:
    """One immutable row of the canonical finite primitive directed table."""

    source: int
    target: int
    cell_displacement: tuple[int, int, int]
    physical_displacement: tuple[float, float, float]
    structure: np.ndarray
    sdot: np.ndarray | None = None

    def __post_init__(self) -> None:
        structure = np.array(self.structure, dtype=COMPLEX, copy=True)
        structure.setflags(write=False)
        object.__setattr__(self, "structure", structure)
        if self.sdot is not None:
            sdot = np.array(self.sdot, dtype=COMPLEX, copy=True)
            sdot.setflags(write=False)
            object.__setattr__(self, "sdot", sdot)

    @property
    def is_onsite(self) -> bool:
        return self.source == self.target and self.cell_displacement == (0, 0, 0)


@dataclass
class MatchedPrimitiveOperator:
    """Finite primitive operator shared by both sides of the WP04 oracle."""

    alat: float
    bonds: tuple[DirectedBond, ...]
    w0: tuple[np.ndarray, ...]
    w1: tuple[np.ndarray, ...]
    c0: tuple[np.ndarray, ...]
    c1: tuple[np.ndarray, ...]
    frames: tuple[np.ndarray, ...]

    def __post_init__(self) -> None:
        if self.alat <= 0.0:
            raise ValueError("alat must be positive")
        if not self.bonds:
            raise ValueError("the matched operator needs at least one bond")
        nsite = len(self.w0)
        if not (len(self.w1) == len(self.c0) == len(self.c1) == len(self.frames) == nsite):
            raise ValueError("site-local arrays have inconsistent lengths")
        norb = len(self.w0[0])
        for site in range(nsite):
            for values in (self.w0[site], self.w1[site], self.c0[site], self.c1[site]):
                if len(values) != norb:
                    raise ValueError("all fixture sites must have the same orbital dimension")
            if self.frames[site].shape != (2, 2):
                raise ValueError("reference frames must be 2x2")
        for index, bond in enumerate(self.bonds):
            if not (0 <= bond.source < nsite and 0 <= bond.target < nsite):
                raise ValueError(f"bond {index} has an invalid endpoint")
            if bond.structure.shape != (norb, norb):
                raise ValueError(f"bond {index} has an invalid orbital block shape")
            if bond.sdot is not None and bond.sdot.shape != (norb, norb):
                raise ValueError(f"bond {index} has an invalid Sdot block shape")
            expected = self.alat * np.asarray(bond.cell_displacement, dtype=float)
            if not np.allclose(bond.physical_displacement, expected, atol=1.0e-14, rtol=0.0):
                raise ValueError(
                    "this exact one-cell fixture requires physical displacement "
                    "to equal alat times the cell displacement"
                )

    @property
    def nsite(self) -> int:
        return len(self.w0)

    @property
    def norb(self) -> int:
        return len(self.w0[0])

    @property
    def block_size(self) -> int:
        return 2 * self.norb

    def basis_slice(self, site: int) -> slice:
        """Return the deterministic spin-block basis slice for one site."""
        if site < 0 or site >= self.nsite:
            raise IndexError("site is outside the primitive basis")
        start = site * self.block_size
        return slice(start, start + self.block_size)

    def basis_index(self, site: int, spin: int, orbital: int) -> int:
        """Return a zero-based deterministic (site, spin, orbital) index."""
        if spin not in (0, 1) or orbital < 0 or orbital >= self.norb:
            raise IndexError("spin or orbital is outside the primitive basis")
        return site * self.block_size + spin * self.norb + orbital

    def canonical_bond_table(self) -> tuple[DirectedBond, ...]:
        """Return a deterministic view of the single shared bond table."""
        return tuple(
            sorted(
                self.bonds,
                key=lambda bond: (bond.source, bond.target, bond.cell_displacement),
            )
        )

    def _link(self, bond: DirectedBond, q: np.ndarray) -> np.ndarray:
        alpha = 2.0 * np.pi * float(np.dot(q, np.asarray(bond.physical_displacement))) / self.alat
        return self.frames[bond.source].conj().T @ z_rotation(alpha) @ self.frames[bond.target]

    def _local_correction(self, site: int) -> np.ndarray:
        return endpoint_factor(self.c0[site], self.c1[site])

    def build_primitive_gbt(self, k: np.ndarray, q: np.ndarray) -> np.ndarray:
        """Build H_GBT(k,q) from the canonical records."""
        k = np.asarray(k, dtype=float)
        q = np.asarray(q, dtype=float)
        result = np.zeros((self.nsite * self.block_size, self.nsite * self.block_size), dtype=COMPLEX)
        for bond in self.canonical_bond_table():
            left = endpoint_factor(self.w0[bond.source], self.w1[bond.source])
            right = endpoint_factor(self.w0[bond.target], self.w1[bond.target])
            linked = spin_lift(bond.structure, self._link(bond, q))
            block = left @ linked @ right
            if bond.is_onsite:
                block = block + self._local_correction(bond.source)
            phase = np.exp(2.0j * np.pi * float(np.dot(k, bond.cell_displacement)))
            i0 = bond.source * self.block_size
            j0 = bond.target * self.block_size
            result[i0 : i0 + self.block_size, j0 : j0 + self.block_size] += phase * block
        return result

    def _site_frame(self, cell: int, site: int, q: np.ndarray) -> np.ndarray:
        cell_vector = np.array([cell, 0.0, 0.0], dtype=float)
        return z_rotation(2.0 * np.pi * float(np.dot(q, cell_vector))) @ self.frames[site]

    def build_explicit_supercell(self, nperiod: int, k_super: np.ndarray, q: np.ndarray) -> np.ndarray:
        """Build the explicit lab-frame N-cell operator from the same records.

        ``k_super`` is expressed in the primitive reciprocal coordinates, so
        the boundary phase for a translated directed record is exp(+2*pi*i*
        k_super.delta).  This makes the fold map explicit and avoids hiding a
        second reciprocal-coordinate convention in the utility.
        """
        if nperiod < 2:
            raise ValueError("a commensurate supercell needs at least two cells")
        k_super = np.asarray(k_super, dtype=float)
        q = np.asarray(q, dtype=float)
        dimension = nperiod * self.nsite * self.block_size
        result = np.zeros((dimension, dimension), dtype=COMPLEX)

        for cell in range(nperiod):
            for bond in self.canonical_bond_table():
                delta = np.asarray(bond.cell_displacement, dtype=int)
                if delta[1] != 0 or delta[2] != 0:
                    raise ValueError("the one-dimensional supercell fixture has y/z cell extent zero")
                target_cell_unwrapped = cell + int(delta[0])
                target_cell = target_cell_unwrapped % nperiod
                source_spin = self._site_frame(cell, bond.source, q)
                target_spin = self._site_frame(target_cell_unwrapped, bond.target, q)
                source_p = np.kron(source_spin, np.eye(self.norb, dtype=COMPLEX))
                target_p = np.kron(target_spin, np.eye(self.norb, dtype=COMPLEX))
                source_w = endpoint_factor(self.w0[bond.source], self.w1[bond.source])
                target_w = endpoint_factor(self.w0[bond.target], self.w1[bond.target])
                source_lab = source_p @ source_w @ source_p.conj().T
                target_lab = target_p @ target_w @ target_p.conj().T
                block = source_lab @ spin_lift(bond.structure) @ target_lab
                phase = np.exp(2.0j * np.pi * float(np.dot(k_super, delta)))
                i0 = (cell * self.nsite + bond.source) * self.block_size
                j0 = (target_cell * self.nsite + bond.target) * self.block_size
                result[i0 : i0 + self.block_size, j0 : j0 + self.block_size] += phase * block

            site = cell * self.nsite
            for local_site in range(self.nsite):
                p = np.kron(self._site_frame(cell, local_site, q), np.eye(self.norb, dtype=COMPLEX))
                i0 = (site + local_site) * self.block_size
                result[i0 : i0 + self.block_size, i0 : i0 + self.block_size] += p @ self._local_correction(local_site) @ p.conj().T
        return result

    def folded_kpoints(self, nperiod: int, k_super: np.ndarray) -> tuple[np.ndarray, ...]:
        """Return the spinor-periodic fold k_l=k_super+(l+1/2)/N along x."""
        k_super = np.asarray(k_super, dtype=float)
        return tuple(
            k_super + np.array([(ell + 0.5) / nperiod, 0.0, 0.0], dtype=float)
            for ell in range(nperiod)
        )

    def folding_unitary(self, nperiod: int, k_super: np.ndarray, q: np.ndarray) -> np.ndarray:
        """Map folded primitive GBT blocks into the explicit supercell basis."""
        dimension = nperiod * self.nsite * self.block_size
        result = np.zeros((dimension, dimension), dtype=COMPLEX)
        for ell, _ in enumerate(self.folded_kpoints(nperiod, k_super)):
            for cell in range(nperiod):
                row = cell * self.nsite * self.block_size
                col = ell * self.nsite * self.block_size
                p = np.zeros((self.nsite * self.block_size, self.nsite * self.block_size), dtype=COMPLEX)
                for site in range(self.nsite):
                    i0 = site * self.block_size
                    site_p = np.kron(self._site_frame(cell, site, q), np.eye(self.norb, dtype=COMPLEX))
                    p[i0 : i0 + self.block_size, i0 : i0 + self.block_size] = site_p
                p *= np.exp(2.0j * np.pi * (ell + 0.5) * cell / nperiod) / np.sqrt(nperiod)
                result[row : row + p.shape[0], col : col + p.shape[1]] = p
        return result


def block_diagonal(blocks: Iterable[np.ndarray]) -> np.ndarray:
    """Small dependency-free block diagonal helper for test diagnostics."""
    blocks = tuple(blocks)
    dimension = sum(block.shape[0] for block in blocks)
    result = np.zeros((dimension, dimension), dtype=COMPLEX)
    offset = 0
    for block in blocks:
        size = block.shape[0]
        result[offset : offset + size, offset : offset + size] = block
        offset += size
    return result
