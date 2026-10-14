#!/usr/bin/env python3
"""Compare production KS H(k) against a Fourier reconstruction of RS hoppings.

This checker consumes the opt-in fixed-U trace after the trace writer emits
`rs_hopping` and `rs_bond_map` rows. It independently applies the production
Fourier convention to the captured real-space neighbor blocks and compares the
result with captured `h_total` matrices at every k point. The test is only a
periodic operator-assembly comparison; it does not validate SCF convergence,
occupation integration, or the radial charge metric.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


FIELDS = (
    "kind", "run_id", "lambda", "iteration", "stage", "site", "type", "l", "spin",
    "k", "band", "row", "col", "order", "weight", "energy", "scalar", "re", "im",
)


def integer(row: dict[str, str], key: str) -> int:
    value = row[key]
    if value == "":
        raise ValueError(f"missing {key} in {row['kind']} record")
    return int(value)


def complex_value(row: dict[str, str]) -> complex:
    return complex(float(row["re"] or 0.0), float(row["im"] or 0.0))


def read_records(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"unexpected CSV header in {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"no records in {path}")
    return rows


def filtered(records: list[dict[str, str]], run_id: str, iteration: int, lam: float) -> list[dict[str, str]]:
    selected = [r for r in records if r["run_id"] == run_id
                and integer(r, "iteration") == iteration
                and math.isclose(float(r["lambda"]), lam, rel_tol=0.0, abs_tol=1e-14)]
    if not selected:
        raise ValueError("no records match the requested run, iteration, and lambda")
    return selected


def _extract_basis_map(rows: list[dict[str, str]], stage: str) -> dict[tuple[int, int], tuple[int, int, int, int]]:
    """Map (site, local orbital-spin index) to type/l/spin/m metadata."""
    result: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    for r in rows:
        if r["kind"] != "basis_map" or r["stage"] != stage:
            continue
        site, row = integer(r, "site"), integer(r, "row")
        key = (site, row)
        value = (integer(r, "type"), integer(r, "l"), integer(r, "spin"), integer(r, "col"))
        if key in result:
            raise ValueError(f"duplicate {stage} basis map entry {key}")
        result[key] = value
    return result


def _extract_realspace(rows: list[dict[str, str]], nsite: int, nblock: int) -> tuple[
        dict[tuple[int, int], np.ndarray], dict[tuple[int, int], tuple[int, int]],
        dict[tuple[int, int], np.ndarray]]:
    """Read full ee blocks and per-site bond maps from one backend capture."""
    hopping: dict[tuple[int, int], np.ndarray] = {}
    seen_h: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for r in rows:
        if r["kind"] != "rs_hopping" or r["stage"] != "RS_OPERATOR":
            continue
        key = (integer(r, "type"), integer(r, "band"))
        ij = (integer(r, "row"), integer(r, "col"))
        if not (0 <= ij[0] < nblock and 0 <= ij[1] < nblock):
            raise ValueError(f"bad real-space block index {ij}; local block dimension is {nblock}")
        if ij in seen_h[key]:
            raise ValueError(f"duplicate real-space hopping element for key={key}, index={ij}")
        seen_h[key].add(ij)
        hopping.setdefault(key, np.zeros((nblock, nblock), dtype=np.complex128))[ij] = complex_value(r)
    if not hopping:
        raise ValueError("missing RS_OPERATOR rs_hopping records")
    for key, seen in seen_h.items():
        if len(seen) != nblock * nblock:
            raise ValueError(f"incomplete real-space block {key}: {len(seen)}/{nblock*nblock}")

    # A bond row stores source site, source type, neighbor slot, target site,
    # fractional R component (order=0..2), and that component in energy.
    bond_targets: dict[tuple[int, int], tuple[int, int]] = {}
    bond_parts: dict[tuple[int, int], dict[int, float]] = defaultdict(dict)
    for r in rows:
        if r["kind"] != "rs_bond_map" or r["stage"] != "RS_OPERATOR":
            continue
        key = (integer(r, "site"), integer(r, "band"))
        source_type, target = integer(r, "type"), integer(r, "col")
        coord = integer(r, "order")
        if coord not in (0, 1, 2) or not (0 <= key[0] < nsite) or not (0 <= target < nsite):
            raise ValueError(f"invalid real-space bond map record: {r}")
        if key in bond_targets and bond_targets[key] != (source_type, target):
            raise ValueError(f"inconsistent bond map metadata for source/neighbor {key}")
        bond_targets[key] = (source_type, target)
        if coord in bond_parts[key]:
            raise ValueError(f"duplicate R component for source/neighbor {key}")
        bond_parts[key][coord] = float(r["energy"])
    if not bond_targets:
        raise ValueError("missing RS_OPERATOR rs_bond_map records")
    for key in bond_targets:
        if set(bond_parts[key]) != {0, 1, 2}:
            raise ValueError(f"incomplete fractional bond vector for {key}")
    bonds = {key: meta for key, meta in bond_targets.items()}
    vectors = {key: np.array([bond_parts[key][i] for i in range(3)]) for key in bond_targets}
    return hopping, bonds, vectors


def _compare_realspace(rs: tuple[dict[tuple[int, int], np.ndarray],
                                   dict[tuple[int, int], tuple[int, int]],
                                   dict[tuple[int, int], np.ndarray]],
                      ks: tuple[dict[tuple[int, int], np.ndarray],
                                   dict[tuple[int, int], tuple[int, int]],
                                   dict[tuple[int, int], np.ndarray]],
                      tolerance: float) -> float:
    rs_hop, rs_bonds, rs_vectors = rs
    ks_hop, ks_bonds, ks_vectors = ks
    if rs_hop.keys() != ks_hop.keys():
        raise ValueError("RS and KS capture different type/neighbor hopping blocks")
    if rs_bonds.keys() != ks_bonds.keys():
        raise ValueError("RS and KS capture different source-site/neighbor maps")
    block_err = max(float(np.max(np.abs(rs_hop[key] - ks_hop[key]))) for key in rs_hop)
    vector_err = 0.0
    for key in rs_bonds:
        if rs_bonds[key] != ks_bonds[key]:
            raise ValueError(f"RS/KS bond mapping differs for {key}: {rs_bonds[key]} != {ks_bonds[key]}")
        vector_err = max(vector_err, float(np.max(np.abs(rs_vectors[key] - ks_vectors[key]))))
    mismatch = max(block_err, vector_err)
    if mismatch > tolerance:
        raise AssertionError(f"RS and KS real-space payloads differ by {mismatch:.3e} (tol {tolerance:.3e})")
    return mismatch


def operator_check(ks_rows: list[dict[str, str]], rs_rows: list[dict[str, str]],
                   tolerance: float, cross_backend: bool) -> tuple[float, float, int, int]:
    # Basis map is site-major; its full range gives the matrix dimension.
    basis = [r for r in ks_rows if r["kind"] == "basis_map" and r["stage"] == "EIGENPAIRS"]
    if not basis:
        raise ValueError("missing EIGENPAIRS basis_map records")
    basis_rows = {integer(r, "row") for r in basis}
    if basis_rows != set(range(len(basis_rows))):
        raise ValueError("basis_map row indices are not contiguous from zero")
    nsite = 1 + max(integer(r, "site") for r in basis)
    if len(basis_rows) % nsite:
        raise ValueError("basis_map dimension is not an integer number of site blocks")
    nblock = len(basis_rows) // nsite
    if nblock != 18:
        raise ValueError(f"expected scalar spd two-spin block size 18, found {nblock}")
    for r in basis:
        row_index = integer(r, "row")
        if integer(r, "site") != row_index // nblock:
            raise ValueError("basis_map is not packed site-major")

    ks_realspace = _extract_realspace(ks_rows, nsite, nblock)
    ks_basis = {(integer(r, "site"), integer(r, "row") % nblock):
                (integer(r, "type"), integer(r, "l"), integer(r, "spin"), integer(r, "col"))
                for r in basis}
    if cross_backend:
        rs_basis = _extract_basis_map(rs_rows, "RS_GREEN")
        if set(rs_basis) != set(ks_basis):
            raise ValueError("RS and KS local basis maps cover different site/orbital indices")
        if rs_basis != ks_basis:
            raise ValueError("RS and KS local basis-map conventions differ")

    # The same real-space arrays must be observed in both solver runs. For a
    # single capture this is only an assembler check; pass --rs-trace to make
    # it an actual RS-vs-KS capture comparison.
    rs_realspace = _extract_realspace(rs_rows, nsite, nblock)
    realspace_error = _compare_realspace(rs_realspace, ks_realspace, tolerance) if cross_backend else 0.0
    hopping, bonds, vectors = rs_realspace
    types_by_site = {site: integer(next(r for r in basis if integer(r, "site") == site), "type") - 1
                     for site in range(nsite)}
    for site, type_id in types_by_site.items():
        expected = {neighbor for block_type, neighbor in hopping if block_type == type_id}
        found = {neighbor for source, neighbor in bonds if source == site}
        if found != expected:
            raise ValueError(f"site {site} bond slots {sorted(found)} do not cover type {type_id} blocks {sorted(expected)}")
    mapped_blocks = {(source_type, neighbor)
                     for (source, neighbor), (source_type, _target) in bonds.items()}
    if mapped_blocks != set(hopping):
        raise ValueError("hopping block set and source-site bond map do not match")

    # K coordinates use scalar values keyed by k and row=Cartesian component.
    kcoords: dict[int, dict[int, float]] = defaultdict(dict)
    for r in ks_rows:
        if r["kind"] == "k_coordinate" and r["stage"] == "EIGENPAIRS":
            k = integer(r, "k")
            component = integer(r, "row")
            if component in kcoords[k]:
                raise ValueError(f"duplicate k coordinate component k={k}, component={component}")
            kcoords[k][component] = float(r["scalar"])
    if not kcoords or any(set(v) != {0, 1, 2} for v in kcoords.values()):
        raise ValueError("missing or incomplete EIGENPAIRS k_coordinate records")

    hk: dict[int, np.ndarray] = {}
    seen_k: dict[int, set[tuple[int, int]]] = defaultdict(set)
    for r in ks_rows:
        if r["kind"] != "h_total" or r["stage"] != "EIGENPAIRS":
            continue
        k = integer(r, "k")
        ij = (integer(r, "row"), integer(r, "col"))
        if not (0 <= ij[0] < nblock * nsite and 0 <= ij[1] < nblock * nsite):
            raise ValueError(f"bad H(k) matrix index {ij}")
        if ij in seen_k[k]:
            raise ValueError(f"duplicate H(k) element k={k}, index={ij}")
        seen_k[k].add(ij)
        hk.setdefault(k, np.zeros((nblock * nsite, nblock * nsite), dtype=np.complex128))[ij] = complex_value(r)
    if not hk:
        raise ValueError("missing EIGENPAIRS h_total matrices")
    if set(hk) != set(kcoords):
        raise ValueError("k-coordinate indices and captured H(k) indices differ")

    worst = 0.0
    for k, actual in hk.items():
        if k not in kcoords:
            raise ValueError(f"missing k coordinates for H(k) index {k}")
        if len(seen_k[k]) != actual.size:
            raise ValueError(f"incomplete H(k) k={k}: {len(seen_k[k])}/{actual.size}")
        kval = np.array([kcoords[k][i] for i in range(3)])
        reconstructed = np.zeros_like(actual)
        for (source, neighbor), (type_id, target) in bonds.items():
            block = hopping.get((type_id, neighbor))
            if block is None:
                raise ValueError(f"missing hopping block type={type_id}, neighbor={neighbor}")
            rvec = vectors[(source, neighbor)]
            phase = np.exp(2j * np.pi * np.dot(kval, rvec))
            i0, j0 = source * nblock, target * nblock
            reconstructed[i0:i0+nblock, j0:j0+nblock] += phase * block
        err = float(np.max(np.abs(reconstructed - actual)))
        worst = max(worst, err)
        if err > tolerance:
            raise AssertionError(f"RS Fourier reconstruction differs from production KS H(k={k}) by {err:.3e} (tol {tolerance:.3e})")
    return worst, realspace_error, len(hk), nsite


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records_csv", type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--lambda", dest="lam", type=float, required=True)
    parser.add_argument("--tolerance", type=float, default=2e-11)
    parser.add_argument("--rs-trace", type=Path,
                         help="separate RS-run records.csv; required for a cross-backend claim")
    parser.add_argument("--rs-run-id")
    parser.add_argument("--rs-iteration", type=int)
    parser.add_argument("--rs-lambda", type=float)
    args = parser.parse_args()
    ks_rows = filtered(read_records(args.records_csv), args.run_id, args.iteration, args.lam)
    cross_backend = args.rs_trace is not None
    if cross_backend:
        if args.rs_run_id is None or args.rs_iteration is None or args.rs_lambda is None:
            parser.error("--rs-trace requires --rs-run-id, --rs-iteration, and --rs-lambda")
        rs_rows = filtered(read_records(args.rs_trace), args.rs_run_id, args.rs_iteration, args.rs_lambda)
    else:
        if any(v is not None for v in (args.rs_run_id, args.rs_iteration, args.rs_lambda)):
            parser.error("RS identifiers require --rs-trace")
        rs_rows = ks_rows
    error, realspace_error, nk, nsite = operator_check(ks_rows, rs_rows, args.tolerance, cross_backend)
    label = "cross-backend" if cross_backend else "same-capture assembly"
    print(f"periodic RS-to-KS operator {label} PASS: sites={nsite}, kpoints={nk}, "
          f"max_RS_realspace_difference={realspace_error:.6e}, max_RS_fourier_vs_KS_Hk={error:.6e}, "
          f"tolerance={args.tolerance:.1e}")


if __name__ == "__main__":
    main()
