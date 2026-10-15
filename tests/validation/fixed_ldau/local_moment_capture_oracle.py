#!/usr/bin/env python3
"""Independent, fail-closed local-moment arithmetic for production traces.

This checker validates the captured fixed-U pairing and independently computes
the candidate local d-shell moments.  It does not call the Fortran correction,
calculate_pl, or NEWRHO routines; captured replay results must be checked
separately against these predictions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "rslmto-fixed-u-trace-v1"
FIELDS = (
    "kind", "run_id", "lambda", "iteration", "stage", "site", "type", "l", "spin",
    "k", "band", "row", "col", "order", "weight", "energy", "scalar", "re", "im",
)
TRACE_IDENTITY_TOL = 1.0e-10  # Ry, same-call production M1 and q0*c
HERMITIAN_TOL = 1.0e-6       # production LDM diagnostics' declared scale
PSD_TOL = 1.0e-3             # production occupation-eigenvalue diagnostic scale
Q2_NEGATIVE_TOL = 1.0e-10    # Ry^2; only roundoff-sized negativity is accepted
# Same-input copied NEWRHO/ATOMSC replay should reproduce same-call outputs to
# roundoff; these are separate from the later operator-response/floor budget.
RADIAL_REPLAY_RHO_TOL = 1.0e-10
RADIAL_REPLAY_ENERGY_TOL_RY = 1.0e-8


def _i(row: dict[str, str], key: str) -> int:
    value = row[key]
    if value == "":
        raise ValueError(f"missing {key} in {row['kind']} row")
    return int(value)


def _f(row: dict[str, str], key: str) -> float:
    value = row[key]
    if value == "":
        raise ValueError(f"missing {key} in {row['kind']} row")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {key} in {row['kind']} row")
    return result


def _z(row: dict[str, str]) -> complex:
    return complex(_f(row, "re"), _f(row, "im"))


def load_capture(directory: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"unsupported trace schema {manifest.get('schema')!r}")
    raw = (directory / "records.csv").read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if manifest.get("csv_sha256") != digest:
        raise ValueError(f"records.csv hash mismatch: manifest={manifest.get('csv_sha256')}, actual={digest}")
    with (directory / "records.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("unexpected CSV field order/schema")
        rows = list(reader)
    run_id = manifest.get("run_id")
    if not rows or any(row["run_id"] != run_id for row in rows):
        raise ValueError("empty capture or inconsistent run_id")
    if int(manifest.get("nsp", -1)) != 1 or int(manifest.get("lmax", -1)) != 2:
        raise ValueError("local-moment checker is restricted to scalar-relativistic nsp=1, lmax=2")
    return manifest, rows


def matrix(rows: list[dict[str, str]], *, kind: str, stage: str, iteration: int,
           site: int, ell: int, spin: int, dim: int) -> np.ndarray:
    selected = [r for r in rows if r["kind"] == kind and r["stage"] == stage
                and _i(r, "iteration") == iteration and _i(r, "site") == site
                and _i(r, "l") == ell and _i(r, "spin") == spin]
    out = np.zeros((dim, dim), dtype=np.complex128)
    seen: set[tuple[int, int]] = set()
    for row in selected:
        ij = (_i(row, "row"), _i(row, "col"))
        if ij in seen or not (0 <= ij[0] < dim and 0 <= ij[1] < dim):
            raise ValueError(f"duplicate/out-of-range {kind} element at {ij}")
        seen.add(ij)
        out[ij] = _z(row)
    if len(seen) != dim * dim:
        raise ValueError(f"incomplete {kind}/{stage} matrix site={site} l={ell} spin={spin}: "
                         f"{len(seen)} of {dim*dim} entries")
    anti = float(np.max(np.abs(out - out.conj().T)))
    if anti > HERMITIAN_TOL:
        raise ValueError(f"{kind} is non-Hermitian: max residual {anti:.6e}")
    return 0.5 * (out + out.conj().T)


def scalar(rows: list[dict[str, str]], *, kind: str, stage: str, iteration: int,
           site: int, ell: int, spin: int, order: int, field: str) -> float:
    found = [r for r in rows if r["kind"] == kind and r["stage"] == stage
             and _i(r, "iteration") == iteration and _i(r, "site") == site
             and _i(r, "l") == ell and _i(r, "spin") == spin
             and _i(r, "order") == order]
    if len(found) != 1:
        raise ValueError(f"expected one {kind}/{stage} order={order} record, got {len(found)}")
    return _f(found[0], field)


def local_correction(p: np.ndarray, m1: np.ndarray, v: np.ndarray,
                     q0: float, c: float, q2: float) -> dict[str, float | bool]:
    if q0 <= 0.0:
        raise ValueError(f"correlated-shell q0 must be positive, got {q0}")
    t = float(np.trace(v @ p).real)
    m1_centered = m1 - c * p
    c_local = c - t / q0
    q2_local = q2 - 2.0 * float(np.trace(v @ m1_centered).real) \
        + float(np.trace(v @ v @ p).real) - t * t / q0
    p_eigs = np.linalg.eigvalsh(0.5 * (p + p.conj().T))
    v_norm = float(np.max(np.abs(v)))
    if q2_local < -Q2_NEGATIVE_TOL:
        raise ValueError(f"computed q2_local is negative ({q2_local:.9e} Ry^2)")
    if float(np.min(p_eigs)) < -PSD_TOL or float(np.max(p_eigs)) > 1.0 + PSD_TOL:
        raise ValueError(f"P eigenvalues outside [0,1] tolerance: [{p_eigs.min():.8g},{p_eigs.max():.8g}]")

    # Independent algebra controls applied to the captured P/M1/q0/c/q2.
    zero = np.zeros_like(v)
    zc = local_correction_no_controls(p, m1, zero, q0, c, q2)
    if zc["c_local_ry"] != c or zc["q2_local_ry2"] != q2:
        raise AssertionError("V=0 no-op control was not exact")
    shift = 0.137
    scalar_v = shift * np.eye(v.shape[0], dtype=np.complex128)
    sc = local_correction_no_controls(p, m1, scalar_v, q0, c, q2)
    scalar_err = abs(sc["q2_local_ry2"] - q2)
    scalar_tol = 2.0e-12 * max(1.0, abs(q2), abs(shift * shift * q0))
    if scalar_err > scalar_tol or abs(sc["c_local_ry"] - (c - shift)) > 1.0e-12:
        raise AssertionError("scalar-shift identity failed")
    return {
        "tr_vp_ry": t,
        "c_local_ry": c_local,
        "q2_local_ry2": q2_local,
        "p_min_eigenvalue": float(p_eigs.min()),
        "p_max_eigenvalue": float(p_eigs.max()),
        "v_max_abs_ry": v_norm,
        "scalar_shift_control_q2_error_ry2": scalar_err,
        "zero_v_noop_exact": True,
    }


def local_correction_no_controls(p: np.ndarray, m1: np.ndarray, v: np.ndarray,
                                 q0: float, c: float, q2: float) -> dict[str, float]:
    if not np.any(v):
        return {"c_local_ry": c, "q2_local_ry2": q2}
    t = float(np.trace(v @ p).real)
    m1c = m1 - c * p
    return {
        "c_local_ry": c - t / q0,
        "q2_local_ry2": q2 - 2.0 * float(np.trace(v @ m1c).real)
        + float(np.trace(v @ v @ p).real) - t * t / q0,
    }


def independent_spectral_moments(manifest: dict[str, Any], rows: list[dict[str, str]]) -> tuple[
        dict[tuple[int, int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray | None]], float]:
    """Rebuild shell P/M1 from raw RS G(E) or exact KS eigenpairs.

    The optional third result is the KS coefficient-space local H0 moment,
    computed from captured H_total - embedded V_in and its actual eigenvectors.
    """
    backend = str(manifest.get("backend", "")).upper()
    result: dict[tuple[int, int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray | None]] = {}
    if backend == "RS":
        maps: dict[tuple[int, int, int], dict[int, int]] = defaultdict(dict)
        for row in rows:
            if row["kind"] != "basis_map" or row["stage"] != "RS_GREEN":
                continue
            key = (_i(row, "site"), _i(row, "l"), _i(row, "spin"))
            m, idx = _i(row, "col"), _i(row, "row")
            if m in maps[key]:
                raise ValueError(f"duplicate RS basis-map m={m} for {key}")
            maps[key][m] = idx
        if not maps:
            raise ValueError("RS capture lacks RS_GREEN basis_map")
        g_rows: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            if row["kind"] == "green_function" and row["stage"] == "RS_GREEN":
                g_rows[(_i(row, "site"), _i(row, "k"))].append(row)
        if not g_rows:
            raise ValueError("RS capture lacks raw RS_GREEN matrices")
        for (site, energy_index), block in g_rows.items():
            g = np.zeros((18, 18), dtype=np.complex128)
            seen: set[tuple[int, int]] = set()
            for row in block:
                ij = (_i(row, "row"), _i(row, "col"))
                if ij in seen:
                    raise ValueError(f"duplicate RS Green entry at site={site}, E-index={energy_index}, {ij}")
                seen.add(ij)
                g[ij] = _z(row)
            if len(seen) != 18 * 18:
                raise ValueError(f"incomplete RS Green matrix at site={site}, E-index={energy_index}")
            energy_set = {_f(r, "energy") for r in block}
            weight_set = {_f(r, "weight") for r in block}
            if len(energy_set) != 1 or len(weight_set) != 1:
                raise ValueError("inconsistent energy/weight within RS Green matrix")
            energy, weight = next(iter(energy_set)), next(iter(weight_set))
            spectral = 1j * (g - g.conj().T) / (2.0 * math.pi)
            if float(np.max(np.abs(spectral - spectral.conj().T))) > 1e-9:
                raise ValueError("RS Green spectral matrix is non-Hermitian")
            for (map_site, ell, spin), idx_by_m in maps.items():
                if map_site != site:
                    continue
                expected_m = list(range(-ell, ell + 1))
                if sorted(idx_by_m) != expected_m:
                    raise ValueError(f"incomplete RS basis map for site={site}, l={ell}, spin={spin}")
                idx = [idx_by_m[m] for m in expected_m]
                key = (site, ell, spin, 0)
                dim = 2 * ell + 1
                if key not in result:
                    result[key] = (np.zeros((dim, dim), complex), np.zeros((dim, dim), complex), None)
                p_acc, m1_acc, _ = result[key]
                p_acc += weight * spectral[np.ix_(idx, idx)]
                m1_acc += weight * energy * spectral[np.ix_(idx, idx)]
        return result, 0.0
    if backend != "KS":
        raise ValueError(f"unsupported backend {backend!r}")

    n = int(manifest.get("basis_size", 0))
    nk = int(manifest.get("kpoint_count", 0))
    if n != 18 or nk < 1:
        raise ValueError(f"unexpected scalar spd KS basis/mesh size n={n}, nk={nk}")
    basis: dict[int, dict[str, int]] = {}
    for row in rows:
        if row["kind"] == "basis_map" and row["stage"] == "EIGENPAIRS":
            i = _i(row, "row")
            if i in basis:
                raise ValueError(f"duplicate KS basis_map index {i}")
            basis[i] = {field: _i(row, field) for field in ("site", "type", "l", "spin", "col")}
    if set(basis) != set(range(n)):
        raise ValueError("KS basis_map is incomplete")
    weights: dict[int, float] = {}
    eig: dict[int, np.ndarray] = {}
    occ: dict[int, np.ndarray] = {}
    vectors: dict[int, np.ndarray] = {}
    hamiltonians: dict[int, np.ndarray] = {}
    for row in rows:
        if row["kind"] == "k_weight":
            k = _i(row, "k")
            if k in weights:
                raise ValueError(f"duplicate k weight {k}")
            weights[k] = _f(row, "weight")
        elif row["kind"] in ("eigenvalue", "occupation") and row["stage"] == "EIGENPAIRS":
            k, band = _i(row, "k"), _i(row, "band")
            target = eig if row["kind"] == "eigenvalue" else occ
            target.setdefault(k, np.full(n, np.nan))[band] = _f(row, "energy" if row["kind"] == "eigenvalue" else "scalar")
        elif row["kind"] == "eigenvector" and row["stage"] == "EIGENPAIRS":
            k, i, b = _i(row, "k"), _i(row, "row"), _i(row, "band")
            vectors.setdefault(k, np.full((n, n), np.nan + 1j*np.nan))[i, b] = _z(row)
        elif row["kind"] == "h_total" and row["stage"] == "EIGENPAIRS":
            k, i, j = _i(row, "k"), _i(row, "row"), _i(row, "col")
            hamiltonians.setdefault(k, np.full((n, n), np.nan + 1j*np.nan))[i, j] = _z(row)
    if set(weights) != set(range(nk)):
        raise ValueError("KS k_weights do not cover the declared mesh")
    if not math.isclose(sum(weights.values()), 1.0, abs_tol=1e-10, rel_tol=0.0):
        raise ValueError("KS k_weights do not sum to one")

    # Map each captured site to the symbolic-atom type index used by v_hub.
    site_type: dict[int, int] = {}
    for row in rows:
        if row["kind"] == "m0_carrier_input" and row["stage"] == "FROZEN_REPLAY" and _i(row, "spin") == 1:
            site = _i(row, "site")
            itype = _i(row, "type") + 1
            if site in site_type and site_type[site] != itype:
                raise ValueError("inconsistent M0 site-to-type map")
            site_type[site] = itype
    if not site_type:
        raise ValueError("KS capture lacks M0 site-to-symbolic-type mapping")
    v_by_type: dict[int, np.ndarray] = {}
    for itype in sorted({_i(r, "type") for r in rows if r["kind"] == "v_hub" and r["stage"] == "Pin_H"}):
        v = np.zeros((18, 18), complex)
        seen: set[tuple[int, int]] = set()
        for row in rows:
            if row["kind"] != "v_hub" or row["stage"] != "Pin_H" or _i(row, "type") != itype:
                continue
            ij = (_i(row, "row"), _i(row, "col"))
            if ij in seen:
                raise ValueError(f"duplicate Vin row for type {itype}")
            seen.add(ij); v[ij] = _z(row)
        if len(seen) != 18*18:
            raise ValueError(f"incomplete Vin matrix for type {itype}")
        v_by_type[itype] = v
    idx_by_key: dict[tuple[int, int, int], list[tuple[int, int]]] = defaultdict(list)
    for idx, labels in basis.items():
        key = (labels["site"], labels["l"], labels["spin"])
        idx_by_key[key].append((labels["col"], idx))
    p_m1: dict[tuple[int, int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray | None]] = {}
    direct_sums: dict[tuple[int, int, int], np.ndarray] = defaultdict(lambda: np.zeros(3, dtype=float))
    max_eig_resid = 0.0
    for k in range(nk):
        if k not in eig or k not in occ or k not in vectors or k not in hamiltonians:
            raise ValueError(f"KS capture missing eigensystem/H at k={k}")
        e, f, cmat, h = eig[k], occ[k], vectors[k], hamiltonians[k]
        if np.any(~np.isfinite(e)) or np.any(~np.isfinite(f)) or np.any(~np.isfinite(cmat)) or np.any(~np.isfinite(h)):
            raise ValueError(f"KS capture has incomplete eigensystem/H matrix at k={k}")
        site_v = np.zeros((n, n), complex)
        for i, li in basis.items():
            for j, lj in basis.items():
                if li["site"] != lj["site"] or li["spin"] != lj["spin"] or li["l"] != 2 or lj["l"] != 2:
                    continue
                itype = site_type.get(li["site"])
                if itype not in v_by_type:
                    raise ValueError(f"missing Vin for site {li['site']} type {itype}")
                ii = li["spin"]*9 + 4 + (li["col"] + 2)
                jj = lj["spin"]*9 + 4 + (lj["col"] + 2)
                site_v[i, j] = v_by_type[itype][ii, jj]
        resid = float(np.max(np.abs(h @ cmat - cmat * e[None, :])))
        max_eig_resid = max(max_eig_resid, resid)
        if resid > 1e-8:
            raise ValueError(f"KS captured eigenpair residual exceeds 1e-8 at k={k}: {resid:.3e}")
        for (site, ell, spin), mapped in idx_by_key.items():
            ordered = sorted(mapped)
            if [m for m, _ in ordered] != list(range(-ell, ell + 1)):
                raise ValueError(f"invalid KS orbital map at site={site} l={ell} spin={spin}")
            idx = [i for _, i in ordered]
            coeff = cmat[idx, :]
            key = (site, ell, spin, 0)
            dim = 2*ell + 1
            if key not in p_m1:
                p_m1[key] = (np.zeros((dim, dim), complex), np.zeros((dim, dim), complex), None)
            p_acc, m1_acc, _ = p_m1[key]
            weighted_f = weights[k] * f
            p_acc += (coeff * weighted_f[None, :]) @ coeff.conj().T
            m1_acc += (coeff * (weighted_f * e)[None, :]) @ coeff.conj().T
            if ell == 2:
                a = coeff
                h0a = (h - site_v) @ cmat
                selected_h0a = h0a[idx, :]
                norms = np.sum(np.abs(a)**2, axis=0)
                cross = np.sum(a.conj() * selected_h0a, axis=0).real
                h0_norms = np.sum(np.abs(selected_h0a)**2, axis=0)
                acc = direct_sums[(site, ell, spin)]
                acc[0] += float(np.sum(weighted_f * norms))
                acc[1] += float(np.sum(weighted_f * cross))
                acc[2] += float(np.sum(weighted_f * h0_norms))
    for site, ell, spin, _ in p_m1:
        moment = p_m1[(site, ell, spin, 0)]
        direct = None
        if ell == 2:
            p_shell, m1_shell, _ = moment
            q0 = float(np.trace(p_shell).real)
            center = float(np.trace(m1_shell).real / q0)
            itype = site_type[site]
            vfull = v_by_type[itype]
            spin0 = spin
            v_shell = vfull[4 + 9*spin0:9 + 9*spin0, 4 + 9*spin0:9 + 9*spin0]
            t = float(np.trace(v_shell @ p_shell).real)
            c_local = center - t/q0
            s0, s1, s2 = direct_sums[(site, ell, spin)]
            direct_q2 = s2 - 2.0*c_local*s1 + c_local*c_local*s0
            direct = np.array([[direct_q2]], dtype=float)
        result[(site, ell, spin, 0)] = moment[0], moment[1], direct
    return result, max_eig_resid


def _case_rows(rows: list[dict[str, str]], kind: str, stage: str, site: int, case: str) -> list[dict[str, str]]:
    return [r for r in rows if r["kind"] == kind and r["stage"] == case
            and _i(r, "site") == site]


def _single_by_order(rows: list[dict[str, str]], kind: str, case: str, site: int,
                     ell: int, spin: int, order: int) -> dict[str, str]:
    found = [r for r in rows if r["kind"] == kind and r["stage"] == case
             and _i(r, "site") == site and _i(r, "l") == ell and _i(r, "spin") == spin
             and _i(r, "order") == order]
    if len(found) != 1:
        raise ValueError(f"expected one {kind}/{case} site={site},l={ell},spin={spin},order={order}; got {len(found)}")
    return found[0]


def _production_calculate_pl(gravity: float, c: float, qpar: float, srdel: float, old_pl: float) -> float:
    """Independent transcription of the scalar phase-shift map, with explicit inputs."""
    ell = 2
    rq = 1.0 / qpar
    delta2 = srdel * srdel
    dc = c - gravity
    dnu = ell + (2.0 * ell + 1.0) / (
        rq * dc / (2.0 * (2.0 * ell + 1.0) * (dc - delta2 * rq)) - 1.0
    )
    return -math.atan(dnu) / math.pi + 0.5 + math.trunc(old_pl)


def _atomsc_radial_charge(rho: list[float], rofi: list[float], a: float, b: float) -> float:
    """AT0MSC's full-NR Simpson integral of rho on the logarithmic mesh."""
    if len(rho) != len(rofi) or len(rho) < 3 or len(rho) % 2 != 1:
        raise ValueError("AT0MSC radial charge requires matching odd-length rho/rofi arrays")
    n = len(rho)
    total = 0.0
    for ir0, value in enumerate(rho):
        ir = ir0 + 1 # Fortran one-based index
        weight = 1.0 / 3.0 if ir == 1 or ir == n else (4.0 / 3.0 if ir % 2 == 0 else 2.0 / 3.0)
        total += weight * a * (rofi[ir0] + b) * value
    return total


def validate_production_radial_replay(rows: list[dict[str, str]], sites: list[int]) -> dict[str, float]:
    """Fail closed unless the copied production NEWRHO replay reproduces its source state.

    `radial_replay_*_max_abs` records compare replay with live same-call outputs.
    `LIVE_STATE` component deltas are exact same-call checks. RAW_PL/LOCAL_PL
    component deltas use the P_OUT QL snapshot, which can precede mixed QL
    consumed by live ATOMSC, so those are reported separately only.
    """
    limits = {
        "radial_replay_rho_max_abs": RADIAL_REPLAY_RHO_TOL,
        "radial_replay_sumec_max_abs_ry": RADIAL_REPLAY_ENERGY_TOL_RY,
        "radial_replay_sumev_max_abs_ry": RADIAL_REPLAY_ENERGY_TOL_RY,
        "radial_replay_ec_max_abs_ry": RADIAL_REPLAY_ENERGY_TOL_RY,
        "radial_replay_ev_max_abs_ry": RADIAL_REPLAY_ENERGY_TOL_RY,
    }
    maxima: dict[str, float] = {}
    for label, limit in limits.items():
        matches = [r for r in rows if r["kind"] == label and r["stage"] == "RADIAL_REPLAY"]
        if len(matches) != len(sites):
            raise ValueError(f"frozen production replay requires one {label} row per site; got {len(matches)} for {len(sites)} sites")
        by_site = {_i(r, "site"): _f(r, "scalar") for r in matches}
        if set(by_site) != set(sites):
            raise ValueError(f"{label} site map does not match production site map")
        value = max(abs(x) for x in by_site.values())
        if not math.isfinite(value) or value > limit:
            raise ValueError(f"frozen production replay {label}={value:.6g} exceeds {limit:.1e}")
        maxima[label] = value

    match_rows = [r for r in rows if r["kind"] == "m0_replay_match" and r["stage"] == "FROZEN_INPUT"]
    if len(match_rows) != len(sites) or {_i(r, "site") for r in match_rows} != set(sites):
        raise ValueError("frozen production replay requires one m0_replay_match row per site")
    pl_error = max(abs(_f(r, "energy")) for r in match_rows)
    ql_stage_gap = max(abs(_f(r, "scalar")) for r in match_rows)
    if pl_error > RADIAL_REPLAY_ENERGY_TOL_RY:
        raise ValueError(f"raw PL replay does not reproduce same-call production PL: {pl_error:.6g}")
    maxima["m0_replay_pl_max_abs"] = pl_error
    maxima["m0_ql_unmixed_vs_live_mixed_stage_gap"] = ql_stage_gap

    all_component_rows = [r for r in rows if r["kind"] == "m0_replay_component_delta"]
    live_component_rows = [r for r in all_component_rows if r["stage"] == "LIVE_STATE"]
    if all_component_rows and not live_component_rows:
        raise ValueError("component replay rows are present but exact LIVE_STATE component baseline is missing")
    if live_component_rows:
        for site in sites:
            site_rows = [r for r in live_component_rows if _i(r, "site") == site]
            by_order = {_i(r, "order"): r for r in site_rows}
            if set(by_order) != set(range(10)):
                raise ValueError(f"site {site} LIVE_STATE component replay must include ATOMSC orders 0..9")
            deltas = []
            for order, row in by_order.items():
                delta, scratch, live = _f(row, "scalar"), _f(row, "re"), _f(row, "im")
                if abs(delta - (scratch - live)) > 2.0e-12:
                    raise ValueError(f"site {site} LIVE_STATE component order {order} delta != scratch-live")
                deltas.append(abs(delta))
            value = max(deltas)
            if value > RADIAL_REPLAY_COMPONENT_TOL_RY:
                raise ValueError(f"site {site} LIVE_STATE ATOMSC component residual {value:.6g} Ry exceeds "
                                 f"{RADIAL_REPLAY_COMPONENT_TOL_RY:.1e} Ry")
            maxima[f"site_{site}_live_component_max_abs_ry"] = value
        for case in ("RAW_PL", "LOCAL_PL"):
            case_rows = [r for r in all_component_rows if r["stage"] == case]
            if case_rows:
                for site in sites:
                    values = [_f(r, "scalar") for r in case_rows if _i(r, "site") == site]
                    if len(values) != 10:
                        raise ValueError(f"site {site} {case} component diagnostic must include orders 0..9")
                    maxima[f"site_{site}_{case.lower()}_component_delta_max_abs_ry"] = max(map(abs, values))
    return maxima


def validate_frozen_replay(rows: list[dict[str, str]], sites: list[int],
                           source_moments: dict[tuple[int, int], dict[str, float]],
                           tr_vp_by_site_spin: dict[tuple[int, int], float],
                           m02_floor: dict[str, Any] | None = None,
                           backend: str = "") -> dict[str, Any]:
    """Check RAW_PL/LOCAL_PL are paired real NEWRHO replays of the same saved state.

    This validates the replay plumbing and reports, but does not enforce, a
    physical total-energy coefficient.  The scalar SUMEV slope is observed
    evidence only; the captured atom-centered ql and Bloch projector live in
    different representations until an independently justified map is proven.
    """
    exact_fields = ("weight", "energy", "scalar", "re", "im")
    production_radial_replay = validate_production_radial_replay(rows, sites)
    identity_kinds = ("m0_replay_rofi", "m0_replay_v_input", "m0_replay_rho_input",
                      "m0_replay_ec_seed", "m0_replay_ev_seed", "m0_replay_control",
                      "m0_replay_radial_parameter")
    result: list[dict[str, Any]] = []
    max_frozen_mismatch = 0.0
    max_pl_formula_error = 0.0
    max_raw_pl_live_error = 0.0
    max_ql_match_error = 0.0
    for site in sites:
        cases: dict[str, dict[str, list[dict[str, str]]]] = {}
        for case in ("RAW_PL", "LOCAL_PL"):
            selected = [r for r in rows if r["stage"] == case and _i(r, "site") == site
                        and (r["kind"].startswith("m0_replay_") or r["kind"] == "m0_scalar_v_identity")]
            if not selected:
                raise ValueError(f"site {site} lacks {case} NEWRHO replay records")
            grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
            for r in selected:
                grouped[r["kind"]].append(r)
            cases[case] = grouped
            etot_rows = grouped.get("m0_replay_etot_available", [])
            if len(etot_rows) != 1 or _f(etot_rows[0], "scalar") != 0.0:
                raise ValueError(f"site {site} {case} must explicitly report NEWRHO ETOT unavailable")
        raw, local = cases["RAW_PL"], cases["LOCAL_PL"]
        for kind in identity_kinds:
            left = raw.get(kind, [])
            right = local.get(kind, [])
            if not left or len(left) != len(right):
                raise ValueError(f"site {site} replay pair has missing/different {kind} rows")
            key = lambda r: (_i(r, "l"), _i(r, "spin"), _i(r, "row"), _i(r, "col"), _i(r, "order"))
            lmap, rmap = {key(r): r for r in left}, {key(r): r for r in right}
            if kind == "m0_replay_control":
                # Order 3 describes the explicit scalar-V intervention and is
                # expected to differ; orders 0..2 define the frozen solver.
                lmap = {k: v for k, v in lmap.items() if k[-1] <= 2}
                rmap = {k: v for k, v in rmap.items() if k[-1] <= 2}
            if lmap.keys() != rmap.keys():
                raise ValueError(f"site {site} replay {kind} index sets differ")
            for ix in lmap:
                if any(lmap[ix][f] != rmap[ix][f] for f in exact_fields):
                    raise ValueError(f"site {site} frozen {kind} differs between RAW_PL/LOCAL_PL at {ix}")

        # Replay controls: nr/ncore/nval and radial solver settings must match.
        ctl = { _i(r, "order"): r for r in raw["m0_replay_control"] }
        if set(ctl) != {0, 1, 2}:
            raise ValueError(f"site {site} replay lacks complete solver controls")
        nr = round(_f(ctl[2], "energy"))
        ncore = round(_f(ctl[2], "scalar"))
        nval = round(_f(ctl[2], "re"))
        if nr < 3 or ncore < 0 or nval < 0:
            raise ValueError(f"site {site} has invalid replay NR/core/valence counts")

        # Independently reproduce calculate_pl for every shell and spin.
        pl_errors: dict[str, float] = {}
        for case in ("RAW_PL", "LOCAL_PL"):
            for spin in (1, 2):
                for ell in range(3):
                    pvals = [_f(_single_by_order(rows, "m0_replay_pl_input", case, site, ell, spin, o), "re")
                             for o in range(5)]
                    predicted = _production_calculate_pl(*pvals)
                    actual = _f(_single_by_order(rows, "m0_replay_pl_return", case, site, ell, spin, 0), "re")
                    err = abs(predicted - actual)
                    pl_errors[f"{case}:l{ell}:spin{spin}"] = err
                    max_pl_formula_error = max(max_pl_formula_error, err)
                    if err > 2.0e-12:
                        raise ValueError(f"site {site} {case} calculate_pl formula mismatch l={ell} spin={spin}: {err:.3e}")
        # RAW_PL must use and reproduce the actual production PL call inputs/output.
        for spin in (1, 2):
            for ell in range(3):
                for order in range(5):
                    raw_input = _f(_single_by_order(rows, "m0_replay_pl_input", "RAW_PL", site, ell, spin, order), "re")
                    actual_input = scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=_i(ctl[0], "iteration"),
                                          site=site, ell=ell, spin=spin, order=order, field="re")
                    max_raw_pl_live_error = max(max_raw_pl_live_error, abs(raw_input - actual_input))
                raw_return = _f(_single_by_order(rows, "m0_replay_pl_return", "RAW_PL", site, ell, spin, 0), "re")
                live_return = scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=_i(ctl[0], "iteration"),
                                     site=site, ell=ell, spin=spin, order=5, field="re")
                max_raw_pl_live_error = max(max_raw_pl_live_error, abs(raw_return - live_return))
        if max_raw_pl_live_error > 2.0e-12:
            raise ValueError(f"site {site} RAW_PL does not reproduce the live calculate_pl snapshot")

        # Only d-shell q2 and its phase-shift gravity center may change.
        ql_differences = []
        for spin in (1, 2):
            for ell in range(3):
                for order in range(3):
                    qraw = _f(_single_by_order(rows, "m0_replay_ql_input", "RAW_PL", site, ell, spin, order), "re")
                    qloc = _f(_single_by_order(rows, "m0_replay_ql_input", "LOCAL_PL", site, ell, spin, order), "re")
                    delta = qloc - qraw
                    ql_differences.append((ell, spin, order, delta))
                    if ell != 2 or order != 2:
                        max_ql_match_error = max(max_ql_match_error, abs(delta))
                if ell == 2:
                    source = source_moments[(site, spin)]
                    qraw0 = _f(_single_by_order(rows, "m0_replay_ql_input", "RAW_PL", site, ell, spin, 0), "re")
                    qraw2 = _f(_single_by_order(rows, "m0_replay_ql_input", "RAW_PL", site, ell, spin, 2), "re")
                    qloc2 = _f(_single_by_order(rows, "m0_replay_ql_input", "LOCAL_PL", site, ell, spin, 2), "re")
                    if abs(qraw0 - source["q0"]) > 1.0e-10 or abs(qraw2 - source["q2"]) > 1.0e-10:
                        raise ValueError(f"site {site} RAW_PL QL does not match same-iteration source moments")
                    if abs(qloc2 - source["q2_local"]) > 1.0e-10:
                        raise ValueError(f"site {site} LOCAL_PL Q2 does not match independent corrected moment")
        if max_ql_match_error > 1.0e-12:
            raise ValueError(f"site {site} replay changed non-q2 or non-d-shell QL input")

        # The corrected center uses the same raw M1/P/V snapshot and vmad frame.
        for spin in (1, 2):
            base = source_moments[(site, spin)]
            local_input = _f(_single_by_order(rows, "m0_replay_pl_input", "LOCAL_PL", site, 2, spin, 0), "re")
            # radial_replay_vmad metadata is emitted once per radial snapshot.
            vmad_rows = [r for r in rows if r["kind"] == "radial_replay_vmad"
                         and r["stage"] == "RADIAL_REPLAY" and _i(r, "site") == site]
            if not vmad_rows:
                raise ValueError(f"site {site} lacks captured radial_replay_vmad metadata")
            vmad = _f(vmad_rows[0], "scalar")
            expected_gc = base["c_local"] - vmad
            if abs(local_input - expected_gc) > 1.0e-10:
                raise ValueError(f"site {site} LOCAL_PL gravity-center mismatch for spin {spin}")

        # Report real NEWRHO outputs. These are measured outcomes, not claims
        # that the spectral first-moment carrier is the complete DFT energy.
        energies: dict[str, dict[int, tuple[float, float]]] = {}
        rhos: dict[str, dict[tuple[int, int], float]] = {}
        for case in ("RAW_PL", "LOCAL_PL"):
            energies[case] = {}
            for r in cases[case].get("m0_replay_energy", []):
                spin = _i(r, "spin")
                energies[case][spin] = (_f(r, "scalar"), _f(r, "energy")) # SUMEV, SUMEC
            if set(energies[case]) != {1, 2}:
                raise ValueError(f"site {site} {case} lacks two-spin NEWRHO energies")
            rhos[case] = {}
            for r in cases[case].get("m0_replay_rho", []):
                key = (_i(r, "row"), _i(r, "spin"))
                if key in rhos[case]:
                    raise ValueError(f"duplicate {case} NEWRHO density row {key}")
                rhos[case][key] = _f(r, "scalar")
            if len(rhos[case]) != 2 * nr:
                raise ValueError(f"site {site} {case} NEWRHO density has {len(rhos[case])} entries, expected {2*nr}")
        sum_ev_raw = sum(x[0] for x in energies["RAW_PL"].values())
        sum_ev_local = sum(x[0] for x in energies["LOCAL_PL"].values())
        sum_ec_raw = sum(x[1] for x in energies["RAW_PL"].values())
        sum_ec_local = sum(x[1] for x in energies["LOCAL_PL"].values())
        rho_max = max(abs(rhos["RAW_PL"][k] - rhos["LOCAL_PL"][k]) for k in rhos["RAW_PL"])
        rofi_rows = sorted(raw["m0_replay_rofi"], key=lambda r: _i(r, "row"))
        rofi = [_f(r, "scalar") for r in rofi_rows]
        a = _f(ctl[0], "energy")
        b = _f(ctl[1], "energy")
        if len(rofi) != nr:
            raise ValueError(f"site {site} replay mesh has {len(rofi)} points, expected NR={nr}")
        rho_integrals = {}
        for case in ("RAW_PL", "LOCAL_PL"):
            per_spin = {}
            for spin in (1, 2):
                values = [rhos[case][(ir, spin)] for ir in range(nr)]
                per_spin[spin] = _atomsc_radial_charge(values, rofi, a, b)
            rho_integrals[case] = per_spin
        raw_spin_moment = rho_integrals["RAW_PL"][1] - rho_integrals["RAW_PL"][2]
        local_spin_moment = rho_integrals["LOCAL_PL"][1] - rho_integrals["LOCAL_PL"][2]
        trvp = sum(tr_vp_by_site_spin[(site, spin)] for spin in (1, 2))
        result.append({
            "site": site, "nr": nr, "ncore": ncore, "nval": nval,
            "max_frozen_input_difference": max_frozen_mismatch,
            "max_calculate_pl_formula_error": max(pl_errors.values()),
            "max_raw_pl_vs_live_snapshot_error": max_raw_pl_live_error,
            "max_unchanged_ql_input_error": max_ql_match_error,
            "q2_input_differences_l_spin_order": ql_differences,
            "sum_ev_raw_ry": sum_ev_raw, "sum_ev_local_ry": sum_ev_local,
            "sum_ec_raw_ry": sum_ec_raw, "sum_ec_local_ry": sum_ec_local,
            "sum_ev_raw_minus_local_ry": sum_ev_raw - sum_ev_local,
            "sum_tr_vp_out_ry": trvp,
            "observed_sum_ev_ratio": ((sum_ev_raw - sum_ev_local) / trvp if abs(trvp) > 1e-14 else None),
            "max_rho_output_difference": rho_max,
            "atomsc_radial_charge_by_spin_raw": rho_integrals["RAW_PL"],
            "atomsc_radial_charge_by_spin_local": rho_integrals["LOCAL_PL"],
            "atomsc_radial_charge_delta_local_minus_raw": {
                spin: rho_integrals["LOCAL_PL"][spin] - rho_integrals["RAW_PL"][spin] for spin in (1, 2)
            },
            "atomsc_sphere_spin_moment_raw": raw_spin_moment,
            "atomsc_sphere_spin_moment_local": local_spin_moment,
            "atomsc_sphere_spin_moment_delta_local_minus_raw": local_spin_moment - raw_spin_moment,
            "status": "replay plumbing validated; energy coefficient not inferred",
        })
        # Additional, deliberately different-operator controls are replayed
        # from the same frozen radial inputs. Their finite differences are
        # reported to establish discrimination, never used to tune the main
        # RAW/LOCAL result.
        for case in ("SCALAR_V", "VOUT_CONTROL", "P_TRANSPOSE_CONTROL"):
            crows = [r for r in rows if r["stage"] == case and _i(r, "site") == site
                     and (r["kind"].startswith("m0_replay_") or r["kind"] == "m0_scalar_v_identity")]
            if not crows:
                raise ValueError(f"site {site} lacks required {case} replay/negative-control records")
            grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
            for r in crows:
                grouped[r["kind"]].append(r)
            for kind in ("m0_replay_rofi", "m0_replay_v_input", "m0_replay_rho_input",
                         "m0_replay_ec_seed", "m0_replay_ev_seed", "m0_replay_radial_parameter"):
                ref = {key(r): r for r in raw.get(kind, [])}
                got = {key(r): r for r in grouped.get(kind, [])}
                if not ref or ref.keys() != got.keys():
                    raise ValueError(f"site {site} {case} lacks the same frozen {kind} payload as RAW_PL")
                if any(ref[k][f] != got[k][f] for k in ref for f in exact_fields):
                    raise ValueError(f"site {site} {case} changed frozen {kind} input")
            rctl = {_i(r, "order"): r for r in grouped.get("m0_replay_control", [])}
            rawctl = {_i(r, "order"): r for r in raw["m0_replay_control"]}
            row_moment = next((r for r in grouped.get("m0_replay_control", []) if _i(r, "order") == 3), None)
            if not {0, 1, 2}.issubset(rctl):
                raise ValueError(f"site {site} {case} lacks A/B/NR/core/valence controls")
            for order in (0, 1, 2):
                if any(rctl[order][f] != rawctl[order][f] for f in exact_fields):
                    raise ValueError(f"site {site} {case} changed frozen control order {order}")
            for spin in (1, 2):
                for ell in range(3):
                    pins = [_f(_single_by_order(rows, "m0_replay_pl_input", case, site, ell, spin, o), "re")
                            for o in range(5)]
                    pcalc = _production_calculate_pl(*pins)
                    pret = _f(_single_by_order(rows, "m0_replay_pl_return", case, site, ell, spin, 0), "re")
                    if abs(pcalc - pret) > 2.0e-12:
                        raise ValueError(f"site {site} {case} calculate_pl formula mismatch l={ell} spin={spin}")
                if case in ("VOUT_CONTROL", "P_TRANSPOSE_CONTROL"):
                    carrier = [r for r in rows if r["kind"] == "m0_carrier_control" and r["stage"] == case
                               and _i(r, "site") == site and _i(r, "spin") == spin]
                    if len(carrier) != 1:
                        raise ValueError(f"site {site} {case} lacks independent carrier-control row")
                    vmad_rows = [r for r in rows if r["kind"] == "radial_replay_vmad"
                                 and r["stage"] == "RADIAL_REPLAY" and _i(r, "site") == site]
                    if not vmad_rows:
                        raise ValueError(f"site {site} lacks vmad for {case} phase-shift check")
                    got_center = _f(_single_by_order(rows, "m0_replay_pl_input", case, site, 2, spin, 0), "re")
                    expected_center = _f(carrier[0], "scalar") - _f(vmad_rows[0], "scalar")
                    got_q0 = _f(_single_by_order(rows, "m0_replay_ql_input", case, site, 2, spin, 0), "re")
                    got_q2 = _f(_single_by_order(rows, "m0_replay_ql_input", case, site, 2, spin, 2), "re")
                    if (abs(got_center - expected_center) > 1.0e-10
                            or abs(got_q0 - _f(carrier[0], "re")) > 1.0e-10
                            or abs(got_q2 - _f(carrier[0], "im")) > 1.0e-10):
                        raise ValueError(f"site {site} {case} replay does not match its independent carrier row")
                if case == "SCALAR_V":
                    scalar_v = _f(row_moment, "energy") if row_moment is not None else 0.01
                    graw = _f(_single_by_order(rows, "m0_replay_pl_input", "RAW_PL", site, 2, spin, 0), "re")
                    gtest = _f(_single_by_order(rows, "m0_replay_pl_input", case, site, 2, spin, 0), "re")
                    if abs((gtest - graw) + scalar_v) > 1.0e-12:
                        raise ValueError(f"site {site} SCALAR_V must subtract positive v from the raw center")
                    for order in (0, 1, 2):
                        qraw = _f(_single_by_order(rows, "m0_replay_ql_input", "RAW_PL", site, 2, spin, order), "re")
                        qtest = _f(_single_by_order(rows, "m0_replay_ql_input", case, site, 2, spin, order), "re")
                        if abs(qtest - qraw) > 1.0e-12:
                            raise ValueError(f"site {site} SCALAR_V changed QL input order {order}")
            e_rows = {_i(r, "spin"): r for r in grouped.get("m0_replay_energy", [])}
            d_rows = {_i(r, "spin"): r for r in grouped.get("m0_replay_delta", [])}
            r_energy = {_i(r, "spin"): r for r in raw.get("m0_replay_energy", [])}
            rho_case_rows = grouped.get("m0_replay_rho", [])
            if set(e_rows) != {1, 2} or set(d_rows) != {1, 2} or len(rho_case_rows) != 2 * nr:
                raise ValueError(f"site {site} {case} missing two-spin energies/deltas or radial rho")
            for spin in (1, 2):
                de = _f(e_rows[spin], "energy") - _f(r_energy[spin], "energy")
                dv = _f(e_rows[spin], "scalar") - _f(r_energy[spin], "scalar")
                if abs(de - _f(d_rows[spin], "energy")) > 2.0e-12 or abs(dv - _f(d_rows[spin], "scalar")) > 2.0e-12:
                    raise ValueError(f"site {site} {case} explicit replay delta disagrees with energy rows")
            rho_by_key = {(_i(r, "row"), _i(r, "spin")): _f(r, "scalar") for r in rho_case_rows}
            if len(rho_by_key) != 2 * nr:
                raise ValueError(f"site {site} {case} has duplicate radial rho output")
            rho_charge = {spin: _atomsc_radial_charge([rho_by_key[(ir, spin)] for ir in range(nr)], rofi, a, b)
                           for spin in (1, 2)}
            rho_ref = {spin: _atomsc_radial_charge([rhos["RAW_PL"][(ir, spin)] for ir in range(nr)], rofi, a, b)
                       for spin in (1, 2)}
            sphere_moment = rho_charge[1] - rho_charge[2]
            ref_moment = rho_ref[1] - rho_ref[2]
            d_ev = sum(_f(d_rows[spin], "scalar") for spin in (1, 2))
            control_info: dict[str, Any] = {}
            etot_rows = grouped.get("m0_replay_etot_available", [])
            if len(etot_rows) != 1 or _f(etot_rows[0], "scalar") != 0.0:
                raise ValueError(f"site {site} {case} must explicitly report NEWRHO ETOT unavailable")
            if case == "SCALAR_V":
                if row_moment is None:
                    raise ValueError("SCALAR_V replay lacks signed intervention metadata")
                scalar_v = _f(row_moment, "energy")
                control_info["scalar_v_ry"] = scalar_v
                control_info["expected_delta_sumev_ry"] = scalar_v * sum(source_moments[(site, s)]["q0"] for s in (1, 2))
                control_info["measured_minus_expected_delta_sumev_ry"] = (
                    d_ev - control_info["expected_delta_sumev_ry"])
                identity_rows = [r for r in grouped.get("m0_scalar_v_identity", []) if _i(r, "spin") in (1, 2)]
                if len(identity_rows) != 2:
                    raise ValueError(f"site {site} SCALAR_V lacks per-spin production identity records")
                expected_from_source = sum(_f(r, "energy") for r in identity_rows)
                measured_from_source = sum(_f(r, "scalar") for r in identity_rows)
                if (abs(expected_from_source - control_info["expected_delta_sumev_ry"]) > 1.0e-12
                        or abs(measured_from_source - d_ev) > 2.0e-12):
                    raise ValueError(f"site {site} SCALAR_V identity rows disagree with independent integrated output")
                # Positive v is removed from the local center: c_local=c-v.
                # Hence SUMEV_raw-SUMEV_local=+v*q0 for V=vI.
                # The 3e-4 Ry M0.1 quadrature allowance is not an M0.2
                # stationary-state inversion floor. Never use it as PASS.
                discrepancy = abs(control_info["measured_minus_expected_delta_sumev_ry"])
                control_info["engineering_3e-4_ry_flag_only"] = (
                    "within" if discrepancy <= 3.0e-4 else "exceeded")
                if m02_floor is None:
                    control_info["acceptance_status"] = "UNASSESSED_M02_FLOOR_MISSING"
                else:
                    if str(m02_floor.get("backend", "")).upper() != backend.upper():
                        raise ValueError(f"M0.2 floor backend {m02_floor.get('backend')} does not match {backend}")
                    floor = float(m02_floor["max_sum_ev_floor_ry"])
                    uncertainty = float(m02_floor["uncertainty_ry"])
                    if floor < 0.0 or uncertainty < 0.0:
                        raise ValueError("M0.2 floor and uncertainty must be nonnegative")
                    allowance = floor + uncertainty
                    control_info["m02_floor_ry"] = floor
                    control_info["m02_floor_uncertainty_ry"] = uncertainty
                    control_info["m02_propagated_allowance_ry"] = allowance
                    control_info["acceptance_status"] = "PASS" if discrepancy <= allowance else "FAIL"
            result[-1].setdefault("operator_controls", {})[case] = {
                **control_info,
                "delta_sumev_ry": d_ev,
                "delta_sumec_ry": sum(_f(d_rows[s], "energy") for s in (1, 2)),
                "delta_atomsc_charge_by_spin": {s: rho_charge[s] - rho_ref[s] for s in (1, 2)},
                "atomsc_radial_charge_by_spin": rho_charge,
                "atomsc_sphere_spin_moment": sphere_moment,
                "delta_atomsc_sphere_spin_moment": sphere_moment - ref_moment,
                "sum_tr_vp_raw_ry": sum(tr_vp_by_site_spin[(site, s)] for s in (1, 2)),
                "measured_delta_over_raw_tr_vp": (
                    d_ev / sum(tr_vp_by_site_spin[(site, s)] for s in (1, 2))
                    if abs(sum(tr_vp_by_site_spin[(site, s)] for s in (1, 2))) > 1.0e-14 else None),
                "status": "measured operator-control response; not an energy closure claim",
            }
    return {"sites": result, "production_radial_replay": production_radial_replay,
            "max_calculate_pl_formula_error": max_pl_formula_error,
            "max_raw_pl_vs_live_snapshot_error": max_raw_pl_live_error,
            "max_unchanged_ql_input_error": max_ql_match_error}


def run_no_u_calibration(manifest: dict[str, Any], rows: list[dict[str, str]], iteration: int,
                         m02_floor: dict[str, Any] | None = None) -> dict[str, Any]:
    availability = [r for r in rows if r["kind"] == "ldm_availability" and r["stage"] == "P_OUT"]
    if len(availability) != 1 or _f(availability[0], "scalar") != 0.0:
        raise ValueError("M0.2 no-U requires explicit P_OUT LDM availability=0")
    if any(r["kind"] == "ldm" and r["stage"] == "P_OUT" for r in rows):
        raise ValueError("no-U capture contains P_OUT LDM rows; refusing stale projector state")
    if any(abs(_z(r)) > 1.0e-14 for r in rows if r["kind"] == "v_hub" and r["stage"] == "Pin_H"):
        raise ValueError("M0.2 no-U capture contains nonzero V_in")
    independent, max_eig_resid = independent_spectral_moments(manifest, rows)
    production_sites = sorted({key[0] for key in independent})
    radial_replay = validate_production_radial_replay(rows, production_sites)
    m1_rows = [r for r in rows if r["kind"] == "ldm_m1" and r["stage"] == "M1_OUT"]
    if not m1_rows:
        raise ValueError("no-U capture lacks nonmutating M1_OUT")
    backend = str(manifest.get("backend", "")).upper()
    tolerance = 5.0e-12 if backend == "KS" else 3.0e-4
    keys = sorted({(_i(r, "site"), _i(r, "l"), _i(r, "spin")) for r in m1_rows})
    m1_error = trace_error = q0_error = band_m1 = 0.0
    site_summary: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for site, ell, spin in keys:
        dim = 2 * ell + 1
        m1 = matrix(rows, kind="ldm_m1", stage="M1_OUT", iteration=iteration,
                    site=site, ell=ell, spin=spin, dim=dim)
        key = (site, ell, spin - 1, 0)
        if key not in independent:
            raise ValueError(f"independent no-U spectrum missing shell {key}")
        p_ind, m1_ind, _ = independent[key]
        err = float(np.max(np.abs(m1 - m1_ind)))
        m1_error = max(m1_error, err)
        if err > tolerance:
            raise ValueError(f"no-U production M1 disagrees with independent {backend} spectrum at {key}: {err:.3e}")
        q0 = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration,
                    site=site, ell=ell, spin=spin, order=0, field="scalar")
        q0c = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration,
                     site=site, ell=ell, spin=spin, order=1, field="scalar")
        trerr = abs(float(np.trace(m1).real) - q0c)
        trace_error = max(trace_error, trerr)
        if trerr > TRACE_IDENTITY_TOL:
            raise ValueError(f"no-U M1 trace/q0c identity fails at {key}: {trerr:.3e} Ry")
        p_err = abs(float(np.trace(p_ind).real) - q0)
        q0_error = max(q0_error, p_err)
        if p_err > 1.0e-4:
            raise ValueError(f"no-U independently reconstructed projector trace differs from q0 at {key}: {p_err:.3e}")
        band_m1 += q0c
        site_summary[site]["band_m1_sum_ry"] += q0c
    atomsc = [r for r in rows if r["kind"] == "energy_component_sumev" and r["stage"] == "ATOMSC_ETOT"]
    replay = [r for r in rows if r["kind"] == "radial_replay_sumev" and r["stage"] == "RADIAL_REPLAY"]
    if not atomsc or not replay:
        raise ValueError("M0.2 capture lacks production ATOMSC and frozen single-NEWRHO SUMEV")
    sum_atomsc = sum(_f(r, "scalar") for r in atomsc)
    sum_replay = sum(_f(r, "scalar") for r in replay)
    scalar_rows = [r for r in rows if r["kind"] == "m0_scalar_v_identity" and r["stage"] == "SCALAR_V"]
    if len(scalar_rows) < 2:
        raise ValueError("M0.2 no-U lacks scalar-shift identity rows")
    scalar_expected = sum(_f(r, "energy") for r in scalar_rows)
    scalar_measured = sum(_f(r, "scalar") for r in scalar_rows)
    discrepancy = abs(scalar_expected - scalar_measured)
    floor_status = "UNASSESSED_M02_FLOOR_MISSING"
    allowance = None
    if m02_floor is not None:
        if str(m02_floor.get("backend", "")).upper() != backend:
            raise ValueError("M0.2 floor backend does not match capture")
        floor = float(m02_floor["max_sum_ev_floor_ry"])
        uncertainty = float(m02_floor["uncertainty_ry"])
        if floor < 0.0 or uncertainty < 0.0:
            raise ValueError("M0.2 floor and uncertainty must be nonnegative")
        allowance = floor + uncertainty
        floor_status = "PASS" if discrepancy <= allowance else "FAIL"
    return {
        "stage": "M0.2 no-U calibration", "evidence_class": "diagnostic; no local-moment correction claim",
        "run_id": manifest["run_id"], "backend": backend, "iteration": iteration,
        "manifest_csv_sha256": manifest["csv_sha256"], "max_eigenpair_residual": max_eig_resid,
        "max_production_independent_m1_error_ry": m1_error, "predeclared_m1_tolerance_ry": tolerance,
        "max_m1_trace_q0c_error_ry": trace_error, "max_independent_projector_q0_error": q0_error,
        "band_m1_sum_q0c_ry": band_m1, "production_atomsc_sumev_ry": sum_atomsc,
        "single_frozen_newrho_sumev_ry": sum_replay,
        "u0_asa_inversion_relaxation_floor_ry": band_m1 - sum_atomsc,
        "newrho_replay_difference_from_atomsc_ry": sum_replay - sum_atomsc,
        "production_radial_replay": radial_replay,
        "scalar_v_expected_raw_minus_local_ry": scalar_expected,
        "scalar_v_measured_raw_minus_local_ry": scalar_measured,
        "scalar_v_discrepancy_ry": discrepancy, "scalar_v_floor_status": floor_status,
        "m02_floor_allowance_ry": allowance, "sites": {s: dict(v) for s, v in site_summary.items()},
        "limitations": ["No production P/LDM is serialized; M1 uses the nonmutating diagnostic path.",
                        "Scalar control is unassessed until backend-matched M0.2 floor evidence is supplied.",
                        "No-U calibration alone does not establish M0.3 carrier or total-energy closure."],
    }


def run(directory: Path, selected_iteration: int | None = None,
        m02_floor: dict[str, Any] | None = None, mode: str = "m03") -> dict[str, Any]:
    manifest, rows = load_capture(directory)
    iteration = selected_iteration if selected_iteration is not None else int(manifest["iteration"])
    if iteration != int(manifest["iteration"]):
        raise ValueError("requested iteration must equal hash-manifest iteration")
    # Filter early: trace files may contain a full SCF history, while the
    # manifest pins one accepted iteration.
    rows = [r for r in rows if _i(r, "iteration") == iteration]
    if mode == "m02-no-u":
        return run_no_u_calibration(manifest, rows, iteration, m02_floor)
    if mode != "m03":
        raise ValueError(f"unsupported mode {mode!r}")
    stages = {r["stage"] for r in rows if _i(r, "iteration") == iteration}
    for needed in ("Pin_H", "P_OUT", "M1_OUT"):
        if needed not in stages:
            raise ValueError(f"capture missing required stage {needed} at iteration {iteration}")
    availability = [r for r in rows if r["kind"] == "ldm_availability" and r["stage"] == "P_OUT"
                    and _i(r, "iteration") == iteration]
    if len(availability) != 1 or _f(availability[0], "scalar") != 1.0:
        raise ValueError("P_OUT LDM is absent/stale; M0.3 requires the actual output projector")
    independent, max_eigenpair_residual = independent_spectral_moments(manifest, rows)

    # One global type-local 18x18 spin-orbital Hubbard potential is emitted
    # at Pin_H. Restrict its up/down d blocks according to PROJECTOR_CONTRACT.
    vrows = [r for r in rows if r["kind"] == "v_hub" and r["stage"] == "Pin_H"
             and _i(r, "iteration") == iteration]
    if not vrows:
        raise ValueError("missing Pin_H V_in matrix")
    type_ids = sorted({_i(r, "type") for r in vrows})
    v_by_type: dict[int, np.ndarray] = {}
    v_anti = 0.0
    for type_id in type_ids:
        vfull = np.zeros((18, 18), dtype=np.complex128)
        seen: set[tuple[int, int]] = set()
        for row in vrows:
            if _i(row, "type") != type_id:
                continue
            ij = (_i(row, "row"), _i(row, "col"))
            if ij in seen or not all(0 <= x < 18 for x in ij):
                raise ValueError(f"duplicate/out-of-range Pin_H V_in element type={type_id} {ij}")
            seen.add(ij)
            vfull[ij] = _z(row)
        if len(seen) != 18 * 18:
            raise ValueError(f"incomplete Pin_H V_in type={type_id}: {len(seen)} of 324 entries")
        anti = float(np.max(np.abs(vfull - vfull.conj().T)))
        v_anti = max(v_anti, anti)
        if anti > HERMITIAN_TOL:
            raise ValueError(f"Pin_H V_in type={type_id} is non-Hermitian: {anti:.3e}")
        v_by_type[type_id] = 0.5 * (vfull + vfull.conj().T)
    correlated_mask = np.zeros((18, 18), dtype=bool)
    for spin in (0, 1):
        correlated_mask[4 + 9*spin:9 + 9*spin, 4 + 9*spin:9 + 9*spin] = True
    off_shell_max = max(float(np.max(np.abs(v[~correlated_mask]))) for v in v_by_type.values())
    if off_shell_max > 1.0e-10:
        raise ValueError(f"fixed-U scope expects d-only V_in; max s/p/off-spin element {off_shell_max:.3e} Ry")

    sites = sorted({_i(r, "site") for r in rows if r["kind"] == "ldm" and r["stage"] == "P_OUT"
                    and _i(r, "iteration") == iteration and _i(r, "l") == 2 and _i(r, "spin") in (1, 2)})
    if not sites:
        raise ValueError("no P_OUT correlated d-shell sites")
    output: list[dict[str, Any]] = []
    trace_max = 0.0
    pl_max = 0.0
    source_formula_max = 0.0
    production_p_max = 0.0
    production_m1_max = 0.0
    direct_coefficient_max = 0.0
    source_moments: dict[tuple[int, int], dict[str, float]] = {}
    tr_vp_by_site_spin: dict[tuple[int, int], float] = {}
    backend_tol = 5.0e-12 if str(manifest.get("backend", "")).upper() == "KS" else 3.0e-4
    for site in sites:
        site_rows = [r for r in rows if _i(r, "iteration") == iteration and _i(r, "site") == site]
        site_type_rows = [r for r in site_rows if r["kind"] == "ldm" and r["stage"] == "P_OUT" and _i(r, "l") == 2]
        if not site_type_rows or len({_i(r, "type") for r in site_type_rows}) != 1:
            raise ValueError(f"site {site} has no unique correlated-element label")
        type_map_rows = [r for r in site_rows if r["kind"] == "m0_carrier_input" and _i(r, "spin") == 1]
        if len(type_map_rows) != 1:
            raise ValueError(f"site {site} lacks unique FROZEN_REPLAY type-index mapping")
        # The historical `ldm.type` field is atomic Z, while v_hub.type is
        # the 1-based symbolic-atom-array index. The explicit M0 map is
        # zero-based symbolic index, so convert it rather than conflating keys.
        type_index0 = _i(type_map_rows[0], "type")
        type_id = type_index0 + 1
        if type_id not in v_by_type:
            raise ValueError(f"site {site} maps to missing Pin_H V_in type index {type_id}")
        vfull = v_by_type[type_id]
        spins: list[dict[str, Any]] = []
        for spin in (1, 2):
            shell_identities = []
            for ell in range(3):
                shell_dim = 2 * ell + 1
                p_shell = matrix(rows, kind="ldm", stage="P_OUT", iteration=iteration, site=site,
                                 ell=ell, spin=spin, dim=shell_dim)
                m1_shell = matrix(rows, kind="ldm_m1", stage="M1_OUT", iteration=iteration, site=site,
                                  ell=ell, spin=spin, dim=shell_dim)
                independent_key = (site, ell, spin - 1, 0)
                if independent_key not in independent:
                    raise ValueError(f"independent raw spectrum lacks site/l/spin key {independent_key}")
                p_independent, m1_independent, _ = independent[independent_key]
                p_err = float(np.max(np.abs(p_shell - p_independent)))
                m1_err = float(np.max(np.abs(m1_shell - m1_independent)))
                production_p_max = max(production_p_max, p_err)
                production_m1_max = max(production_m1_max, m1_err)
                if p_err > backend_tol or m1_err > backend_tol:
                    raise ValueError(f"raw {manifest['backend']} spectral reconstruction mismatch at "
                                     f"site={site},l={ell},spin={spin}: P={p_err:.3e}, M1={m1_err:.3e}, "
                                     f"limit={backend_tol:.3e}")
                q0_shell = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration,
                                  site=site, ell=ell, spin=spin, order=0, field="scalar")
                q0c_shell = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration,
                                   site=site, ell=ell, spin=spin, order=1, field="scalar")
                shell_trace_error = abs(float(np.trace(m1_shell).real) - q0c_shell)
                trace_max = max(trace_max, shell_trace_error)
                if shell_trace_error > TRACE_IDENTITY_TOL:
                    raise ValueError(f"M1/q0c identity fails at site={site} l={ell} spin={spin}: "
                                     f"{shell_trace_error:.3e} Ry")
                if abs(float(np.trace(p_shell).real) - q0_shell) > 1.0e-4:
                    raise ValueError(f"P trace and radial q0 disagree at site={site} l={ell} spin={spin}")
                shell_identities.append({"l": ell, "q0": q0_shell,
                                         "trace_p": float(np.trace(p_shell).real),
                                         "trace_m1_minus_q0c_ry": shell_trace_error})
            p = matrix(rows, kind="ldm", stage="P_OUT", iteration=iteration, site=site, ell=2, spin=spin, dim=5)
            m1 = matrix(rows, kind="ldm_m1", stage="M1_OUT", iteration=iteration, site=site, ell=2, spin=spin, dim=5)
            q0 = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration, site=site,
                        ell=2, spin=spin, order=0, field="scalar")
            q0c = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration, site=site,
                         ell=2, spin=spin, order=1, field="scalar")
            center_first = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration, site=site,
                                  ell=2, spin=spin, order=1, field="energy")
            center2_q2 = scalar(rows, kind="radial_moment", stage="P_OUT", iteration=iteration, site=site,
                                ell=2, spin=spin, order=2, field="energy")
            if q0 <= 0.0:
                raise ValueError(f"site {site} spin {spin} has nonpositive d-shell q0")
            center = center_first
            q2 = center2_q2 - center * center
            identity = abs(float(np.trace(m1).real) - q0c)
            trace_max = max(trace_max, identity)
            if identity > TRACE_IDENTITY_TOL:
                raise ValueError(f"M1/q0c identity fails at site={site} spin={spin}: {identity:.3e} Ry")
            q0p = float(np.trace(p).real)
            if abs(q0p - q0) > 1.0e-4:
                raise ValueError(f"P trace and radial q0 disagree at site={site} spin={spin}: {q0p} vs {q0}")
            spin_zero = spin - 1
            vi = vfull[4 + 9*spin_zero:9 + 9*spin_zero, 4 + 9*spin_zero:9 + 9*spin_zero]
            corr = local_correction(p, m1, vi, q0, center, q2)
            source_moments[(site, spin)] = {
                "q0": q0, "c": center, "q2": q2,
                "t": float(corr["tr_vp_ry"]), "c_local": float(corr["c_local_ry"]),
                "q2_local": float(corr["q2_local_ry2"]),
            }
            tr_vp_by_site_spin[(site, spin)] = float(corr["tr_vp_ry"])
            if str(manifest.get("backend", "")).upper() == "KS":
                direct_q2 = float(independent[(site, 2, spin - 1, 0)][2][0, 0])
                direct_err = abs(float(corr["q2_local_ry2"]) - direct_q2)
                direct_coefficient_max = max(direct_coefficient_max, direct_err)
                if direct_err > 1.0e-10:
                    raise ValueError(f"independent KS H-Vin coefficient q2 mismatch at site={site}, spin={spin}: "
                                     f"{direct_err:.3e} Ry^2")
            source_in = [r for r in site_rows if r["kind"] == "m0_carrier_input" and _i(r, "spin") == spin]
            source_local = [r for r in site_rows if r["kind"] == "m0_carrier_local" and _i(r, "spin") == spin]
            source_identity = [r for r in site_rows if r["kind"] == "m0_carrier_identity" and _i(r, "spin") == spin]
            if len(source_in) != 1 or len(source_local) != 1 or len(source_identity) != 1:
                raise ValueError(f"site {site} spin {spin} lacks complete source M0 carrier snapshot")
            sin, slocal, sident = source_in[0], source_local[0], source_identity[0]
            source_vals = {
                "q0": _f(sin, "scalar"), "c": _f(sin, "energy"), "q2": _f(sin, "re"),
                "t": _f(slocal, "weight"), "c_local": _f(slocal, "energy"),
                "q2_local": _f(slocal, "re"),
            }
            independent_vals = {
                "q0": q0, "c": center, "q2": q2,
                "t": float(corr["tr_vp_ry"]), "c_local": float(corr["c_local_ry"]),
                "q2_local": float(corr["q2_local_ry2"]),
            }
            source_formula_error = max(abs(source_vals[k] - independent_vals[k]) for k in source_vals)
            source_formula_max = max(source_formula_max, source_formula_error)
            if source_formula_error > 1.0e-10:
                raise ValueError(f"independent/source local-moment mismatch site={site} spin={spin}: "
                                 f"{source_formula_error:.3e}")
            if _f(sident, "scalar") > 1.0e-10:
                raise ValueError(f"source complex-contraction imaginary residual at site={site} spin={spin} "
                                 f"is {_f(sident, 'scalar'):.3e}")
            # Legacy P_OUT gravity_center and calculate_pl input must be the
            # same call-time value. Their equality is independent of later vmad.
            radial_gc = scalar(rows, kind="radial_parameter", stage="P_OUT", iteration=iteration,
                               site=site, ell=2, spin=spin, order=-1, field="energy")
            pl_gc = scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                           site=site, ell=2, spin=spin, order=0, field="scalar")
            pl_err = abs(radial_gc - pl_gc)
            pl_max = max(pl_max, pl_err)
            if pl_err > 1.0e-10:
                raise ValueError(f"CALCULATE_PL center is not paired with P_OUT radial center: {pl_err:.3e} Ry")
            spins.append({
                "spin_zero_based": spin_zero, "q0": q0, "trace_p": q0p,
                "all_spd_shell_trace_identities": shell_identities,
                "center_ry": center, "q2_ry2": q2,
                "trace_m1_minus_q0c_ry": identity,
                "pl_input_gravity_center_ry": pl_gc,
                "pl_map_input": {
                    "C": scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                                site=site, ell=2, spin=spin, order=1, field="re"),
                    "QPAR": scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                                   site=site, ell=2, spin=spin, order=2, field="re"),
                    "SRDEL": scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                                    site=site, ell=2, spin=spin, order=3, field="re"),
                    "oldPL": scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                                    site=site, ell=2, spin=spin, order=4, field="re"),
                    "newPL": scalar(rows, kind="pl_map_input", stage="CALCULATE_PL", iteration=iteration,
                                    site=site, ell=2, spin=spin, order=5, field="re"),
                },
                **corr,
                "wrong_transposed_p_tr_vp_ry": float(np.trace(vi @ p.T).real),
                "wrong_transposed_p_center_ry": center - float(np.trace(vi @ p.T).real) / q0,
                "wrong_transposed_p_q2_ry2": local_correction_no_controls(p.T, m1, vi, q0, center, q2)["q2_local_ry2"],
                "source_formula_max_abs_error": source_formula_error,
            })
        output.append({"site": site, "type": type_id, "spins": spins})
    if pl_max > 1.0e-10:
        raise ValueError("PL center pairing tolerance exceeded")
    replay = validate_frozen_replay(rows, sites, source_moments, tr_vp_by_site_spin,
                                    m02_floor, str(manifest.get("backend", "")))
    return {
        "stage": "M0.3 independent trace arithmetic only",
        "evidence_class": "diagnostic; not a replay or Gate M0 pass",
        "run_id": manifest["run_id"], "backend": manifest["backend"], "iteration": iteration,
        "manifest_csv_sha256": manifest["csv_sha256"],
        "executable_sha256": manifest.get("executable_sha256"),
        "input_sha256": manifest.get("input_sha256"),
        "operator_convention": "P_ab=psi_a psi_b*, M1_ab=sum f epsilon psi_a psi_b*, Tr(VP)=Re Tr(V@P)",
        "energy_units": "Ry; c before vmad; q2 in Ry^2",
        "v_hermiticity_max_abs": v_anti,
        "v_off_d_max_abs_ry": off_shell_max,
        "max_same_call_trace_identity_ry": trace_max,
        "max_pl_input_pairing_error_ry": pl_max,
        "max_independent_vs_source_formula_error": source_formula_max,
        "raw_spectrum_reconstruction": {
            "max_P_residual": production_p_max,
            "max_M1_residual_ry": production_m1_max,
            "predeclared_backend_tolerance": backend_tol,
            "ks_max_eigenpair_residual": max_eigenpair_residual,
            "ks_direct_HminusVin_q2_local_error_ry2": direct_coefficient_max if str(manifest.get("backend", "")).upper() == "KS" else None,
        },
        "sites": output,
        "frozen_replay": replay,
        "controls": {
            "zero_v_exact_no_op": True,
            "scalar_shift": "tested per site/spin with v=0.137 Ry",
            "wrong_transpose": "reported as a negative control; do not claim discrimination unless its induced moment differs beyond the declared M0.3 carrier floor",
        },
        "limitations": [
            "This script computes candidate corrected moments but does not invoke production calculate_pl/newrho.",
            "The caller must separately validate paired frozen RAW_PL/LOCAL_PL production replay outputs.",
            "A wrong-stage V_out or transpose result can be nondiscriminating at a converged fixed point and must be measured, not assumed to fail.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture_dir", type=Path, help="directory containing hash-matched manifest.json and records.csv")
    parser.add_argument("--iteration", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--m02-floor-json", type=Path,
                        help="backend-matched M0.2 evidence with max_sum_ev_floor_ry and uncertainty_ry")
    parser.add_argument("--mode", choices=("m02-no-u", "m03"), default="m03")
    args = parser.parse_args()
    floor = json.loads(args.m02_floor_json.read_text(encoding="utf-8")) if args.m02_floor_json else None
    result = run(args.capture_dir, args.iteration, floor, args.mode)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
