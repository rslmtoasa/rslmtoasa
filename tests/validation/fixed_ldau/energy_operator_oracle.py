#!/usr/bin/env python3
"""Independent KS finite-H oracle for RSLMTO_FIXED_U_TRACE captures.

The capture is supplied by the normal production calculation.  This script
does not call Fortran kernels: it reconstructs the operator, projectors,
moments, and finite-temperature Hellmann-Feynman slopes with NumPy.

Required directory contents are manifest.json and records.csv.  The manifest
schema and the CSV record vocabulary are described by --help and in
REMEDIATION_GATES.md.  At present generalized overlap matrices are accepted
for eigensolves; local projector/moment comparison uses the captured basis
coefficients and must be interpreted with the separately recorded projector
metric.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment diagnostic
    raise SystemExit("energy_operator_oracle.py requires NumPy") from exc


SCHEMA = "rslmto-fixed-u-trace-v1"
CSV_FIELDS = (
    "kind", "run_id", "lambda", "iteration", "stage", "site", "type", "l", "spin",
    "k", "band", "row", "col", "order", "weight", "energy", "scalar", "re", "im",
)


def _int(row: dict[str, str], key: str, default: int | None = None) -> int | None:
    value = row.get(key, "")
    return default if value == "" else int(value)


def _float(row: dict[str, str], key: str, default: float | None = None) -> float | None:
    value = row.get(key, "")
    return default if value == "" else float(value)


def _complex(row: dict[str, str]) -> complex:
    return complex(float(row.get("re", "0") or 0.0), float(row.get("im", "0") or 0.0))


def read_capture(directory: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    manifest_path = directory / "manifest.json"
    records_path = directory / "records.csv"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"unsupported trace schema: {manifest.get('schema')!r}")
    raw = records_path.read_bytes()
    expected_hash = manifest.get("csv_sha256")
    actual_hash = hashlib.sha256(raw).hexdigest()
    if expected_hash and expected_hash != actual_hash:
        raise ValueError(f"records.csv SHA256 mismatch: manifest={expected_hash}, actual={actual_hash}")
    with records_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != CSV_FIELDS:
            raise ValueError(f"records.csv header must be exactly {','.join(CSV_FIELDS)}")
        records = list(reader)
    if not records:
        raise ValueError("records.csv has no records")
    run_id = str(manifest.get("run_id", ""))
    if any(r["run_id"] != run_id for r in records):
        raise ValueError("manifest run_id does not match every CSV row")
    return manifest, records


def read_manifest(directory: Path) -> dict[str, Any]:
    """Read and hash-check a capture manifest without loading its CSV into RAM."""
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"unsupported trace schema: {manifest.get('schema')!r}")
    expected_hash = manifest.get("csv_sha256")
    if expected_hash:
        digest = hashlib.sha256()
        with (directory / "records.csv").open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != expected_hash:
            raise ValueError(f"records.csv SHA256 mismatch: manifest={expected_hash}, actual={digest.hexdigest()}")
    return manifest


def _complex_matrix(records: list[dict[str, str]], kind: str, k: int, n: int,
                    stage: str | None = None) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.complex128)
    seen: set[tuple[int, int]] = set()
    for row in records:
        if row["kind"] != kind or _int(row, "k", -1) not in (-1, k):
            continue
        if stage is not None and row["stage"] != stage:
            continue
        i, j = _int(row, "row"), _int(row, "col")
        if i is None or j is None or not (0 <= i < n and 0 <= j < n):
            raise ValueError(f"bad {kind} matrix index at k={k}: row={i}, col={j}")
        if (i, j) in seen:
            raise ValueError(f"duplicate {kind} element k={k}, row={i}, col={j}")
        matrix[i, j] = _complex(row)
        seen.add((i, j))
    if not seen:
        raise ValueError(f"missing {kind} matrix for k={k}")
    # Hermitian matrices are expected as full records; sparse Hubbard V is
    # permitted and is completed by Hermitian conjugation only when absent.
    if kind == "v_hub":
        for i in range(n):
            for j in range(n):
                if (i, j) not in seen and (j, i) in seen:
                    matrix[i, j] = matrix[j, i].conjugate()
    elif len(seen) != n * n:
        raise ValueError(f"{kind} matrix has {len(seen)} of {n*n} elements at k={k}")
    residual = float(np.max(np.abs(matrix - matrix.conj().T)))
    if residual > 1e-9:
        raise ValueError(f"{kind} is not Hermitian at k={k}: residual {residual:.3e}")
    return 0.5 * (matrix + matrix.conj().T)


def _embedded_hubbard_potential(records: list[dict[str, str]], basis: list[dict[str, int]], n: int) -> np.ndarray:
    """Embed type-local spin-orbital V_U blocks into the Bloch basis."""
    local_blocks: dict[int, np.ndarray] = {}
    seen: dict[int, set[tuple[int, int]]] = defaultdict(set)
    nper_spin = max(((int(v["l"]) + 1) ** 2 for v in basis), default=0)
    local_n = 2 * nper_spin
    for row in records:
        if row["kind"] != "v_hub" or row["stage"] != "Pin_H":
            continue
        type_id = int(_int(row, "type")); i = int(_int(row, "row")); j = int(_int(row, "col"))
        block = local_blocks.setdefault(type_id, np.zeros((local_n, local_n), dtype=np.complex128))
        key = (i, j)
        if key in seen[type_id]:
            raise ValueError(f"duplicate type-local V_U entry type={type_id}, index={key}")
        seen[type_id].add(key)
        if not (0 <= i < local_n and 0 <= j < local_n):
            raise ValueError(f"V_U index {key} outside local spin-orbital dimension {local_n}")
        block[i, j] = _complex(row)
    if not local_blocks:
        raise ValueError("trace lacks Pin_H v_hub matrix records")

    result = np.zeros((n, n), dtype=np.complex128)
    for i, bi in enumerate(basis):
        li = int(bi["spin"]) * nper_spin + int(bi["l"]) ** 2 + int(bi["m"]) + int(bi["l"])
        for j, bj in enumerate(basis):
            if (bi["site"], bi["type"]) != (bj["site"], bj["type"]):
                continue
            lj = int(bj["spin"]) * nper_spin + int(bj["l"]) ** 2 + int(bj["m"]) + int(bj["l"])
            block = local_blocks.get(int(bi["type"]))
            if block is not None:
                result[i, j] = block[li, lj]
    residual = float(np.max(np.abs(result - result.conj().T)))
    if residual > 1e-9:
        raise ValueError(f"embedded V_U is not Hermitian: residual {residual:.3e}")
    return 0.5 * (result + result.conj().T)


def _generalized_hermitian_eigh(h: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solve H C = S C e with C†SC=I using a Cholesky transform."""
    lmat = np.linalg.cholesky(s)
    left = np.linalg.solve(lmat, h)
    transformed = np.linalg.solve(lmat.conj(), left.T).T
    transformed = 0.5 * (transformed + transformed.conj().T)
    evals, y = np.linalg.eigh(transformed)
    coeff = np.linalg.solve(lmat.conj().T, y)
    return evals, coeff


def _fermi_occupations(evals: np.ndarray, kweights: np.ndarray, nelectrons: float, temp: float) -> tuple[np.ndarray, float]:
    """Solve the fixed-N Fermi occupations; states are along axis 1."""
    if temp <= 0.0:
        flat = sorted((float(evals[k, b]), float(kweights[k]), k, b)
                      for k in range(len(kweights)) for b in range(evals.shape[1]))
        occ = np.zeros_like(evals)
        remaining = nelectrons
        mu = flat[0][0]
        for energy, weight, k, band in flat:
            fill = min(1.0, max(0.0, remaining / weight))
            occ[k, band] = fill
            remaining -= fill * weight
            mu = energy
            if remaining <= 1e-12:
                break
        if remaining > 1e-8:
            raise ValueError(f"electron count {nelectrons} exceeds available states")
        return occ, mu

    if nelectrons <= 0 or nelectrons >= float(len(kweights) * evals.shape[1]):
        raise ValueError("finite-T electron count must lie strictly between zero and capacity")
    lo = float(np.min(evals) - 80.0 * temp - 1.0)
    hi = float(np.max(evals) + 80.0 * temp + 1.0)

    def occupations(mu: float) -> np.ndarray:
        x = np.clip((evals - mu) / temp, -700.0, 700.0)
        return 1.0 / (1.0 + np.exp(x))

    for _ in range(240):
        mid = 0.5 * (lo + hi)
        count = float(np.sum(occupations(mid) * kweights[:, None]))
        if count < nelectrons:
            lo = mid
        else:
            hi = mid
        if hi - lo < 2e-14:
            break
    mu = 0.5 * (lo + hi)
    return occupations(mu), mu


def _electronic_free_energy(evals: np.ndarray, kweights: np.ndarray, occ: np.ndarray, temp: float) -> tuple[float, float]:
    band_energy = float(np.sum(kweights[:, None] * occ * evals))
    if temp <= 0:
        return band_energy, band_energy
    clipped = np.clip(occ, 1e-300, 1.0 - 1e-16)
    entropy_term = temp * float(np.sum(kweights[:, None] * (
        clipped * np.log(clipped) + (1.0 - clipped) * np.log(1.0 - clipped))))
    return band_energy + entropy_term, band_energy


def _basis_map(records: list[dict[str, str]], n: int) -> list[dict[str, int]]:
    mapping: dict[int, dict[str, int]] = {}
    for row in records:
        if row["kind"] != "basis_map":
            continue
        index = _int(row, "row")
        m = _int(row, "col")
        values = {key: _int(row, key) for key in ("site", "type", "l", "spin")}
        if index is None or m is None or any(v is None for v in values.values()):
            raise ValueError("basis_map records require row, col=m, site, type, l, and spin")
        if index in mapping:
            raise ValueError(f"duplicate basis_map row {index}")
        mapping[index] = {"m": m, **values}  # type: ignore[arg-type]
    if set(mapping) != set(range(n)):
        raise ValueError(f"basis_map covers indices {sorted(mapping)}, expected 0..{n-1}")
    return [mapping[i] for i in range(n)]


def _group_k(records: list[dict[str, str]], kind: str) -> dict[int, list[dict[str, str]]]:
    result: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        if row["kind"] == kind:
            k = _int(row, "k")
            if k is None:
                raise ValueError(f"{kind} record lacks k index")
            result[k].append(row)
    return result


def analyze_rs_capture(directory: Path, tolerance: float) -> dict[str, Any]:
    """Independently integrate captured retarded Green matrices for RS LDM.

    This validates the actual RS occupation contraction against its raw G(E)
    payload and the recorded production Simpson weights. It is distinct from
    the finite-H exact-eigensystem/Hellmann-Feynman oracle below.
    """
    manifest = read_manifest(directory)
    if str(manifest.get("backend", "")).upper() != "RS":
        raise ValueError("RS Green oracle requires manifest backend='RS'")
    iteration = int(manifest["iteration"])
    records_path = directory / "records.csv"
    with records_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != CSV_FIELDS:
            raise ValueError(f"records.csv header must be exactly {','.join(CSV_FIELDS)}")
        energies: dict[tuple[int, int], tuple[float, float, float]] = {}
        by_site_energy: dict[tuple[int, int], np.ndarray] = {}
        seen: set[tuple[int, int, int, int]] = set()
        basis_rows: list[dict[str, str]] = []
        ldm_rows: list[dict[str, str]] = []
        moment_rows: list[dict[str, str]] = []
        vhub_rows: list[dict[str, str]] = []
        sum_ev_rows: list[dict[str, str]] = []
        run_id = str(manifest.get("run_id", ""))
        for row in reader:
            if row["run_id"] != run_id:
                raise ValueError("manifest run_id does not match every CSV row")
            if int(_int(row, "iteration")) != iteration:
                continue
            kind = row["kind"]
            if kind == "fermi_step":
                site, ie = int(_int(row, "site")), int(_int(row, "k"))
                energies[(site, ie)] = (float(_float(row, "energy")), float(_float(row, "weight")),
                                         float(_float(row, "scalar")))
            elif kind == "basis_map":
                basis_rows.append(row)
            elif kind == "ldm":
                ldm_rows.append(row)
            elif kind == "v_hub" and row["stage"] == "Pin_H":
                vhub_rows.append(row)
            elif kind == "sum_ev" and row["stage"] == "POST_SCF_RADIAL":
                sum_ev_rows.append(row)
            elif kind == "radial_moment" and row["stage"] == "P_OUT":
                moment_rows.append(row)
            elif kind == "green_function":
                site, ie = int(_int(row, "site")), int(_int(row, "k"))
                i, j = int(_int(row, "row")), int(_int(row, "col"))
                n = int(manifest.get("basis_size", max(i, j) + 1))
                matrix = by_site_energy.setdefault((site, ie), np.zeros((n, n), dtype=np.complex128))
                key = (site, ie, i, j)
                if key in seen:
                    raise ValueError(f"duplicate RS Green element {key}")
                seen.add(key)
                matrix[i, j] = _complex(row)

    if not energies or not by_site_energy:
        raise ValueError("RS trace needs fermi_step and green_function records")
    basis_map_captured = bool(basis_rows)
    if not basis_rows:
        # Older opt-in captures predate explicit basis_map rows. Reconstruct
        # only the documented scalar-collinear l,m,spin packing, and label it
        # as inferred in the result. This is not accepted for the full H/k
        # oracle, which requires a captured basis map.
        lmax = manifest.get("lmax")
        norb = manifest.get("norb_per_spin")
        if lmax is None or norb is None:
            raise ValueError("RS trace lacks basis_map; supply manifest lmax and norb_per_spin for the legacy source convention")
        sites = sorted({site for site, _ in by_site_energy})
        types = {}
        with records_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["kind"] == "green_function":
                    types[int(_int(row, "site"))] = int(_int(row, "type"))
        for site in sites:
            for spin in range(2):
                for orbital in range(int(norb)):
                    ell = math.isqrt(orbital)
                    if ell > int(lmax):
                        continue
                    m = orbital - ell * ell - ell
                    basis_rows.append({"site": str(site), "type": str(types[site]), "l": str(ell),
                                       "spin": str(spin), "row": str(spin * int(norb) + orbital),
                                       "col": str(m)})

    p_expected: dict[tuple[int, int, int, int], np.ndarray] = {}
    moment_expected: dict[tuple[int, int, int, int, int], float] = defaultdict(float)
    full_projectors: dict[int, np.ndarray] = {}
    global_spectral_m0 = 0.0
    global_spectral_m1 = 0.0
    for (site, ie), g in by_site_energy.items():
        if (site, ie) not in energies:
            raise ValueError(f"missing RS quadrature/Fermi row for site={site}, energy={ie}")
        energy, weight, fermi_step = energies[(site, ie)]
        spectral = 1j * (g - g.conj().T) / (2.0 * math.pi)
        if np.max(np.abs(spectral - spectral.conj().T)) > 1e-9:
            raise ValueError(f"RS spectral matrix is not Hermitian at site={site}, energy={ie}")
        # Simpson weights already encode the partial final interval to EF;
        # do not multiply by fermi_step again.
        full_projectors.setdefault(site, np.zeros_like(g))
        full_projectors[site] += weight * spectral
        global_spectral_m0 += weight * float(np.trace(spectral).real)
        global_spectral_m1 += weight * energy * float(np.trace(spectral).real)
        for map_row in basis_rows:
            if int(_int(map_row, "site")) != site:
                continue
            index = int(_int(map_row, "row")); ell = int(_int(map_row, "l"))
            spin = int(_int(map_row, "spin")); type_id = int(_int(map_row, "type"))
            m = int(_int(map_row, "col"))
            key = (site, type_id, ell, spin)
            dim = 2 * ell + 1
            if key not in p_expected:
                p_expected[key] = np.zeros((dim, dim), dtype=np.complex128)
            for map_col in basis_rows:
                if (int(_int(map_col, "site")), int(_int(map_col, "type")), int(_int(map_col, "l")),
                    int(_int(map_col, "spin"))) != key:
                    continue
                j = int(_int(map_col, "row")); m2 = int(_int(map_col, "col"))
                p_expected[key][m + ell, m2 + ell] += weight * spectral[index, j]
            # The recorded order-0 quantity is the spectral occupation Q0.
            # Higher recorded radial_moment fields are Q0*gravity_center
            # and Q0*(gravity_center**2+central_variance), not raw spectral
            # M1/M2; keep them out of this independent Green comparison.
            moment_expected[(*key, 0)] += weight * spectral[index, index].real

    if not ldm_rows:
        raise ValueError("RS trace lacks P_OUT LDM records")
    max_ldm = 0.0
    for row in ldm_rows:
        if row["stage"] != "P_OUT":
            continue
        # Legacy LDM rows used one-based site/spin labels while Green/basis
        # matrices used zero-based labels. Normalize before comparison.
        site = int(_int(row, "site")) - (0 if basis_map_captured else 1)
        # LDM shell labels remain one-based even after local site labels were
        # corrected to the zero-based Green/basis-map convention.
        spin = int(_int(row, "spin")) - 1
        key = (site, int(_int(row, "type")), int(_int(row, "l")), spin)
        i, j = int(_int(row, "row")), int(_int(row, "col"))
        if key not in p_expected:
            raise ValueError(f"no Green-derived projector for LDM key={key}")
        max_ldm = max(max_ldm, float(abs(p_expected[key][i, j] - _complex(row))))
    if max_ldm > tolerance:
        raise ValueError(f"RS Green-vs-P_OUT projector residual {max_ldm:.3e} exceeds {tolerance:.3e}")
    if not moment_rows:
        raise ValueError("RS trace lacks radial_moment rows")
    max_moment = 0.0
    for row in moment_rows:
        if int(_int(row, "order")) != 0:
            continue
        key = (int(_int(row, "site")), int(_int(row, "type")), int(_int(row, "l")),
               int(_int(row, "spin")) - 1, 0)
        if key in moment_expected:
            max_moment = max(max_moment, abs(float(_float(row, "scalar")) - moment_expected[key]))
    tr_vp_by_state: dict[str, dict[str, float]] = {}
    if vhub_rows:
        local_n = int(manifest.get("basis_size", 0))
        v_by_type: dict[int, np.ndarray] = {}
        for row in vhub_rows:
            type_id = int(_int(row, "type"))
            v = v_by_type.setdefault(type_id, np.zeros((local_n, local_n), dtype=np.complex128))
            v[int(_int(row, "row")), int(_int(row, "col"))] = _complex(row)
        # Reconstruct each local spin-orbital density from shell LDMs using
        # the same inferred/captured basis map as the Green projector.
        labels = {}
        for row in basis_rows:
            labels[(int(_int(row, "site")), int(_int(row, "type")), int(_int(row, "l")),
                    int(_int(row, "spin")), int(_int(row, "col")))] = int(_int(row, "row"))
        stages = ("Pin_H", "P_OUT", "P_MIX")
        p_records: dict[str, dict[int, np.ndarray]] = {stage: {} for stage in stages}
        for row in ldm_rows:
            stage = row["stage"]
            if stage not in p_records:
                continue
            # Pre-basis-map traces label these two fields one-based; newer
            # traces use the captured zero-based local site/spin convention.
            site = int(_int(row, "site")) - (0 if basis_map_captured else 1)
            spin = int(_int(row, "spin")) - 1
            type_id, ell = int(_int(row, "type")), int(_int(row, "l"))
            matrix = p_records[stage].setdefault(site, np.zeros((local_n, local_n), dtype=np.complex128))
            for local_row, local_col in ((int(_int(row, "row")), int(_int(row, "col"))),):
                i = labels.get((site, type_id, ell, spin, local_row - ell))
                j = labels.get((site, type_id, ell, spin, local_col - ell))
                if i is not None and j is not None:
                    matrix[i, j] = _complex(row)
        for stage in stages:
            for site, pmat in p_records[stage].items():
                type_ids = {int(_int(r, "type")) for r in ldm_rows if r["stage"] == stage and
                            int(_int(r, "site")) - (0 if basis_map_captured else 1) == site}
                if len(type_ids) == 1 and next(iter(type_ids)) in v_by_type:
                    tr = np.trace(v_by_type[next(iter(type_ids))] @ pmat)
                    tr_vp_by_state.setdefault(stage, {})[str(site)] = float(np.real(tr))
    return {"schema": SCHEMA, "run_id": manifest["run_id"], "backend": "RS",
            "source_sha": manifest["source_sha"], "executable_sha256": manifest["executable_sha256"],
            "capture_sha256": manifest.get("csv_sha256"), "iteration": iteration,
            "max_green_projector_residual": max_ldm, "max_projected_moment_residual": max_moment,
            "tr_vhub_p_by_stage_ry": tr_vp_by_state,
            "green_spectral_m0": global_spectral_m0,
            "green_spectral_m1_ry": global_spectral_m1,
            "radial_sumev_ry": sum(float(_float(r, "scalar")) for r in sum_ev_rows),
            "green_m1_minus_radial_sumev_ry": global_spectral_m1 -
                sum(float(_float(r, "scalar")) for r in sum_ev_rows),
            "status": "PASS_RS_GREEN_PROJECTOR" if max_ldm <= tolerance and max_moment <= tolerance else "FAIL_RS_GREEN_PROJECTOR",
            "basis_map_source": "captured" if basis_map_captured else "inferred_legacy_scalar_collinear",
            "scope_note": "RS Green/projector and global spectral M1 are independently integrated; comparison to radial SUMEV is diagnostic only and does not establish the physical ETOT double-counting correction."}


def analyze_capture(directory: Path, fd_steps: list[float], tolerance: float) -> dict[str, Any]:
    manifest, records = read_capture(directory)
    n = int(manifest["basis_size"])
    nk = int(manifest["kpoint_count"])
    lam0 = float(manifest["lambda"])
    temp = float(manifest["temperature_ry"])
    ne = float(manifest["electron_count"])
    basis = _basis_map(records, n)
    # Index matrix rows once. Scanning the 1.2M-row trace for every k point
    # makes the independent oracle needlessly O(nk * nrecords).
    hs = _group_k([r for r in records if r["kind"] == "h_total" and r["stage"] == "EIGENPAIRS"], "h_total")
    ss = _group_k(records, "s_overlap")
    eig_records = _group_k(records, "eigenvalue")
    vec_records = _group_k(records, "eigenvector")
    occ_records = _group_k(records, "occupation")
    kw_rows = {int(_int(r, "k")): float(_float(r, "weight")) for r in records if r["kind"] == "k_weight"}
    if set(kw_rows) != set(range(nk)):
        raise ValueError(f"k weights cover {sorted(kw_rows)}, expected 0..{nk-1}")
    kweights = np.array([kw_rows[k] for k in range(nk)], dtype=float)
    if not np.isclose(np.sum(kweights), 1.0, atol=1e-10, rtol=0.0):
        raise ValueError(f"k weights sum to {np.sum(kweights):.16g}, expected 1")

    hfull, vhub, overlap = [], [], []
    evals_capture = np.zeros((nk, n), dtype=float)
    occ_capture = np.zeros((nk, n), dtype=float)
    coeff_capture = np.zeros((nk, n, n), dtype=np.complex128)
    v_embedded = _embedded_hubbard_potential(records, basis, n)
    for k in range(nk):
        hfull.append(_complex_matrix(hs.get(k, []), "h_total", k, n, stage="EIGENPAIRS"))
        vhub.append(v_embedded)
        if manifest.get("metric", "identity") == "overlap":
            overlap.append(_complex_matrix(records, "s_overlap", k, n))
        else:
            overlap.append(np.eye(n, dtype=np.complex128))
        if k not in eig_records or k not in vec_records or k not in occ_records:
            raise ValueError(f"missing eigenvalue/eigenvector/occupation records for k={k}")
        for row in eig_records[k]:
            band = _int(row, "band")
            if band is None or not 0 <= band < n:
                raise ValueError(f"bad band index in eigenvalue row: {band}")
            evals_capture[k, band] = float(_float(row, "energy"))
        for row in occ_records[k]:
            band = _int(row, "band")
            if band is None or not 0 <= band < n:
                raise ValueError(f"bad band index in occupation row: {band}")
            occ_capture[k, band] = float(_float(row, "scalar"))
        for row in vec_records[k]:
            i, b = _int(row, "row"), _int(row, "band")
            if i is None or b is None or not (0 <= i < n and 0 <= b < n):
                raise ValueError(f"bad eigenvector index at k={k}: ({i},{b})")
            coeff_capture[k, i, b] = _complex(row)

    if np.any(occ_capture < -1e-10) or np.any(occ_capture > 1.0 + 1e-10):
        raise ValueError("captured scalar occupations outside [0,1]")
    count_capture = float(np.sum(kweights[:, None] * occ_capture))
    if abs(count_capture - ne) > 1e-7:
        raise ValueError(f"captured occupations give N={count_capture:.12g}, manifest says {ne:.12g}")

    max_eig_resid = 0.0
    max_orth_resid = 0.0
    max_occ_resid = 0.0
    p_local: dict[tuple[int, int, int, int], np.ndarray] = {}
    moment_local: dict[tuple[int, int, int, int, int], float] = defaultdict(float)
    tr_vp = 0.0
    for k in range(nk):
        h, v, s = hfull[k], vhub[k], overlap[k]
        if np.linalg.norm(h - h.conj().T, ord=np.inf) > 1e-9:
            raise ValueError(f"H is not Hermitian at k={k}")
        c = coeff_capture[k]
        residual = h @ c - (s @ c) * evals_capture[k][None, :]
        max_eig_resid = max(max_eig_resid, float(np.max(np.abs(residual))))
        gram = c.conj().T @ s @ c
        max_orth_resid = max(max_orth_resid, float(np.max(np.abs(gram - np.eye(n)))))
        # Check recorded occupations against the manifest Fermi state.  This
        # explicitly includes the recorded chemical potential and temperature.
        if temp > 0:
            mu = float(manifest["fermi_ry"])
            target = 1.0 / (1.0 + np.exp(np.clip((evals_capture[k] - mu) / temp, -700, 700)))
            max_occ_resid = max(max_occ_resid, float(np.max(np.abs(target - occ_capture[k]))))
        p = (c * occ_capture[k][None, :]) @ c.conj().T
        tr_vp += kweights[k] * float(np.real(np.trace(p @ v)))
        for i, bmap_i in enumerate(basis):
            for j, bmap_j in enumerate(basis):
                if (bmap_i["site"], bmap_i["type"], bmap_i["l"], bmap_i["spin"]) == \
                   (bmap_j["site"], bmap_j["type"], bmap_j["l"], bmap_j["spin"]):
                    key = (bmap_i["site"], bmap_i["type"], bmap_i["l"], bmap_i["spin"])
                    if key not in p_local:
                        dim = 2 * key[2] + 1
                        p_local[key] = np.zeros((dim, dim), dtype=np.complex128)
                    a, b = bmap_i["m"] + bmap_i["l"], bmap_j["m"] + bmap_j["l"]
                    p_local[key][a, b] += kweights[k] * p[i, j]
        # Projected radial moments from exactly the same captured states.
        for band in range(n):
            f = occ_capture[k, band] * kweights[k]
            e = evals_capture[k, band]
            for idx, bmap_i in enumerate(basis):
                key = (bmap_i["site"], bmap_i["type"], bmap_i["l"], bmap_i["spin"])
                q = abs(c[idx, band]) ** 2
                for order in range(3):
                    moment_local[(*key, order)] += f * q * e**order

    if max_eig_resid > tolerance:
        raise ValueError(f"captured eigensystem residual {max_eig_resid:.3e} exceeds {tolerance:.3e}")
    if max_orth_resid > tolerance:
        raise ValueError(f"captured C†SC residual {max_orth_resid:.3e} exceeds {tolerance:.3e}")
    if max_occ_resid > max(tolerance, 1e-8):
        raise ValueError(f"captured fermi occupation residual {max_occ_resid:.3e} exceeds tolerance")

    # Compare production LDM and projected radial moments.  Missing records
    # are reported explicitly rather than silently treated as a pass.
    ldm_stages = {r["stage"] for r in records if r["kind"] == "ldm"}
    if "P_OUT" not in ldm_stages:
        raise ValueError("trace lacks Pout LDM records")
    ldm_max = {stage: 0.0 for stage in ldm_stages}
    for row in records:
        if row["kind"] != "ldm":
            continue
        key = (int(_int(row, "site")), int(_int(row, "type")), int(_int(row, "l")), int(_int(row, "spin")) - 1)
        i, j = _int(row, "row"), _int(row, "col")
        pred = p_local.get(key)
        if pred is None or i is None or j is None:
            continue
        residual = abs(pred[i, j] - _complex(row))
        ldm_max[row["stage"]] = max(ldm_max[row["stage"]], float(residual))
    if ldm_max.get("P_OUT", 0.0) > tolerance:
        raise ValueError(f"independent KS projector vs captured P_OUT residual {ldm_max['P_OUT']:.3e} exceeds {tolerance:.3e}")

    moment_rows = [r for r in records if r["kind"] == "radial_moment"]
    max_moment_resid = 0.0
    for row in moment_rows:
        # The capture stores Q0, Q0*gravity_center, and
        # Q0*(gravity_center**2 + central_variance).  Only order zero is the
        # raw projected spectral moment available from eigenpairs; orders 1/2
        # are radial reconstruction parameters, not sum(f*epsilon**p).
        if row["stage"] != "P_OUT" or int(_int(row, "order")) != 0:
            continue
        key = (int(_int(row, "site")), int(_int(row, "type")), int(_int(row, "l")),
               int(_int(row, "spin")) - 1, 0)
        expected = moment_local.get(key)
        if expected is None:
            raise ValueError(f"no eigenstate-derived projected moment for {key}")
        max_moment_resid = max(max_moment_resid, abs(expected - float(_float(row, "scalar"))))
    if not moment_rows:
        raise ValueError("trace lacks raw radial_moment records; cannot validate the spectral projector mapping")
    if max_moment_resid > tolerance:
        raise ValueError(f"independent projected moment residual {max_moment_resid:.3e} exceeds {tolerance:.3e}")

    # Independent finite-H oracle: H(lambda)=H(lambda0)+(lambda-lambda0)V/lambda0.
    # The fixed-N Helmholtz free energy includes the entropy term; its
    # derivative is Tr(PV_unit) at finite T. Band energy alone has an extra
    # occupation-redistribution term, which is reported rather than mistaken
    # for the Hellmann-Feynman expectation.
    v_unit = [v / lam0 for v in vhub]
    fermi = []
    band = []
    slopes = []
    for step in fd_steps:
        energies_plus, energies_minus = [], []
        occ_plus, occ_minus = [], []
        for sign, target in ((1.0, energies_plus), (-1.0, energies_minus)):
            vals = []
            for k in range(nk):
                hk = hfull[k] + sign * step * v_unit[k]
                if manifest.get("metric", "identity") == "overlap":
                    ev, _ = _generalized_hermitian_eigh(hk, overlap[k])
                else:
                    ev = np.linalg.eigvalsh(0.5 * (hk + hk.conj().T))
                vals.append(ev)
            target.extend(vals)
        ep, em = np.array(energies_plus), np.array(energies_minus)
        fp, mu_p = _fermi_occupations(ep, kweights, ne, temp)
        fm, mu_m = _fermi_occupations(em, kweights, ne, temp)
        free_p, band_p = _electronic_free_energy(ep, kweights, fp, temp)
        free_m, band_m = _electronic_free_energy(em, kweights, fm, temp)
        slope_free = (free_p - free_m) / (2.0 * step)
        slope_band = (band_p - band_m) / (2.0 * step)
        slopes.append({"step": step, "free_slope": slope_free, "band_slope": slope_band,
                       "expected_tr_pv": tr_vp / lam0, "free_abs_error": abs(slope_free - tr_vp / lam0),
                       "mu_plus": mu_p, "mu_minus": mu_m})
        fermi.append((mu_p, mu_m))
        band.append((band_p, band_m))

    max_fd_error = max(item["free_abs_error"] for item in slopes)
    # Central-difference truncation decreases as h^2; require a resolved
    # convergence trend, and an endpoint error consistent with tolerance.
    if len(slopes) >= 3 and not (slopes[-1]["free_abs_error"] <= slopes[0]["free_abs_error"] * 1.05 + tolerance):
        raise ValueError("finite-T free-energy Hellmann-Feynman errors did not decrease across the chosen steps")
    report = {
        "schema": SCHEMA,
        "run_id": manifest["run_id"],
        "backend": manifest["backend"],
        "lambda": lam0,
        "source_sha": manifest["source_sha"],
        "executable_sha256": manifest["executable_sha256"],
        "capture_sha256": manifest.get("csv_sha256"),
        "count_from_occupations": count_capture,
        "trace_vp": tr_vp,
        "max_eigensystem_residual": max_eig_resid,
        "max_metric_orthogonality_residual": max_orth_resid,
        "max_fermi_occupation_residual": max_occ_resid,
        "max_pout_projector_residual": ldm_max["P_OUT"],
        "max_projected_moment_residual": max_moment_resid,
        "finite_difference_slopes": slopes,
        "max_free_energy_slope_error": max_fd_error,
        "status": "PASS_FINITE_H_ORACLE" if max_fd_error <= tolerance else "FAIL_FINITE_H_ORACLE",
        "scope_note": "Finite-H/operator and projected spectral-moment oracle only; not a radial ETOT identity or converged-material result.",
    }
    if max_fd_error > tolerance:
        raise ValueError(f"finite-temperature free-energy slope error {max_fd_error:.3e} exceeds {tolerance:.3e}")
    return report


def self_test() -> dict[str, Any]:
    """Independent small Hermitian test of fixed-N finite-T Hellmann-Feynman."""
    # Closed 2x2 fixture. H0 and V do not commute, ensuring off-diagonal
    # response is present. The reference derivative is computed as Tr(PV)
    # from the eigensystem, not by the finite-difference routine under test.
    h0 = np.array([[0.1, 0.07j], [-0.07j, 0.6]], dtype=np.complex128)
    v = np.array([[0.3, 0.11 + 0.04j], [0.11 - 0.04j, -0.2]], dtype=np.complex128)
    lam0, temp, ne = 0.8, 0.03, 0.85
    e0, c0 = np.linalg.eigh(h0 + lam0 * v)
    occ0, _ = _fermi_occupations(e0[None, :], np.array([1.0]), ne, temp)
    p0 = (c0 * occ0[0][None, :]) @ c0.conj().T
    expected = float(np.real(np.trace(p0 @ v)))
    rows = []
    for step in (1e-2, 5e-3, 2.5e-3):
        free = []
        for sign in (-1, 1):
            e, _ = np.linalg.eigh(h0 + (lam0 + sign * step) * v)
            f, _ = _fermi_occupations(e[None, :], np.array([1.0]), ne, temp)
            free.append(_electronic_free_energy(e[None, :], np.array([1.0]), f, temp)[0])
        slope = (free[1] - free[0]) / (2.0 * step)
        rows.append({"step": step, "slope": slope, "expected_tr_pv": expected, "abs_error": abs(slope - expected)})
    errors = [r["abs_error"] for r in rows]
    if not (errors[-1] < errors[0] and errors[-1] < 2e-5):
        raise ValueError(f"self-test finite-difference convergence failed: {errors}")
    return {"status": "PASS_SELF_TEST_ONLY", "temperature_ry": temp, "electron_count": ne,
            "expected_tr_pv": expected, "finite_difference_slopes": rows,
            "scope_note": "Parser/math self-test only; not production evidence."}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", nargs="?", type=Path, help="directory containing manifest.json and records.csv")
    parser.add_argument("--self-test", action="store_true", help="run the independent finite-T oracle self-test")
    parser.add_argument("--fd-steps", nargs="+", type=float, default=[1e-3, 1e-4, 1e-5],
                        help="dimensionless central-difference steps in lambda")
    parser.add_argument("--tolerance", type=float, default=1e-7,
                        help="absolute residual tolerance for matrix and energy identities (Ry for slopes)")
    parser.add_argument("--json", type=Path, help="also write the JSON report to this path")
    args = parser.parse_args()
    if args.self_test:
        result = self_test()
    elif args.capture is not None:
        manifest = read_manifest(args.capture)
        if str(manifest.get("backend", "")).upper() == "RS":
            result = analyze_rs_capture(args.capture, args.tolerance)
        else:
            result = analyze_capture(args.capture, args.fd_steps, args.tolerance)
    else:
        parser.error("provide a capture directory or --self-test")
    encoded = json.dumps(result, indent=2, sort_keys=True)
    print(encoded)
    if args.json:
        args.json.write_text(encoded + "\n")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(f"energy_operator_oracle: FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
