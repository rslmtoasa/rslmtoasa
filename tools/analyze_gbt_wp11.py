#!/usr/bin/env python3
"""Analyse the WP11 LKAG ``J(q)`` versus frozen-spiral response.

The program deliberately treats the exchange file as a *real-space pair
table*, not as a ready-made ``J(0)``.  ``jij.out`` contains one representative
bond per requested pair record, so shell degeneracies must be supplied when a
record stands for a symmetry shell (for example 8 and 6 for the first two
bcc shells).  This prevents an under-supported LKAG table from being silently
presented as a converged Fourier sum.

The public helpers are standard-library-only and are useful from tests or
small notebooks:

``read_jij_file``
    Read the production ``jij.out`` layout.
``construct_jq``
    Construct scalar one-sublattice ``J(q)`` in Ry.
``construct_exchange_matrix``
    Construct a directed multi-sublattice ``J_ab(q)`` matrix.
``analyze_wp11``
    Compare one or more LKAG range sets with WP09/WP10 sweep rows and return
    a JSON-serializable evidence record.

The scalar centrosymmetric convention is

    J(q) = sum_R m_R J_R cos(2*pi*q.R),

where ``m_R`` is the total number of equivalent vectors represented by the
record.  ``q`` is the code's Cartesian ``2*pi/alat`` coordinate and the bond
vectors printed by the current bcc fixture are in units of ``alat``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
WP09_TOOL = ROOT / "tools" / "analyze_gbt_wp09.py"
_WP09 = None

RY_PER_UNIT = {
    "ry": 1.0,
    "mry": 1.0e-3,
    "ev": 1.0 / 13.605693122994,
    "mev": 1.0e-3 / 13.605693122994,
}
Q_TOLERANCE = 1.0e-8
DEFAULT_CURVATURE_RELATIVE_TOLERANCE = 0.10
DEFAULT_RANGE_RELATIVE_TOLERANCE = 0.05
DEFAULT_GOLDSTONE_ABSOLUTE_TOLERANCE = 1.0e-10


def _number(value: str | float | int) -> float:
    return float(str(value).strip().replace("D", "E").replace("d", "e"))


def _q_tuple(q: Sequence[float]) -> tuple[float, float, float]:
    if len(q) != 3:
        raise ValueError(f"q must have three components, got {q!r}")
    return tuple(float(value) for value in q)  # type: ignore[return-value]


def _q_key(q: Sequence[float]) -> tuple[float, float, float]:
    return tuple(round(float(value), 8) for value in q)  # type: ignore[return-value]


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in vector))


def _unit(vector: Sequence[float]) -> tuple[float, float, float]:
    length = _norm(vector)
    if length <= 0.0:
        raise ValueError("direction must be nonzero")
    return tuple(float(value) / length for value in vector)  # type: ignore[return-value]


def _phase(q: Sequence[float], bond: Sequence[float]) -> float:
    return 2.0 * math.pi * sum(float(a) * float(b) for a, b in zip(q, bond))


def _load_wp09():
    global _WP09
    if _WP09 is not None:
        return _WP09
    spec = importlib.util.spec_from_file_location("gbt_wp09_analysis_for_wp11", WP09_TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {WP09_TOOL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _WP09 = module
    return module


def read_jij_file(
    path: str | Path,
    *,
    multiplicities: Sequence[float] | None = None,
    energy_unit: str = "mry",
    bond_scale: float = 1.0,
) -> list[dict[str, Any]]:
    """Read production ``jij.out`` records into a normalized pair table.

    The current writer emits seven fields::

        type_i type_j Rx Ry Rz Jij distance

    An optional eighth field is accepted as an inline shell multiplicity.  An
    explicit ``multiplicities`` argument takes precedence and is recommended
    for the unmodified production output.  Multiplicity means the *total*
    number of equivalent directed vectors represented by the row; it is not a
    factor to apply again for the ``cos`` form of a centrosymmetric shell.
    """

    unit = energy_unit.strip().lower()
    if unit not in RY_PER_UNIT:
        raise ValueError(f"unsupported Jij energy unit {energy_unit!r}")
    if bond_scale <= 0.0:
        raise ValueError("bond_scale must be positive")

    source = Path(path)
    records: list[dict[str, Any]] = []
    for line_number, raw in enumerate(source.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 7:
            raise ValueError(f"{source}:{line_number}: expected at least 7 fields, got {len(fields)}")
        try:
            species_i = int(round(_number(fields[0])))
            species_j = int(round(_number(fields[1])))
            bond = tuple(bond_scale * _number(fields[index]) for index in (2, 3, 4))
            jij = _number(fields[5])
            distance = _number(fields[6])
            has_inline_multiplicity = len(fields) >= 8
            inline_multiplicity = _number(fields[7]) if has_inline_multiplicity else 1.0
        except ValueError as exc:
            raise ValueError(f"{source}:{line_number}: malformed Jij record") from exc
        if species_i < 1 or species_j < 1:
            raise ValueError(f"{source}:{line_number}: species indices must be positive")
        if inline_multiplicity <= 0.0:
            raise ValueError(f"{source}:{line_number}: multiplicity must be positive")
        records.append({
            "index": len(records) + 1,
            "species_i": species_i,
            "species_j": species_j,
            "R": bond,
            "J_input": jij,
            "J_Ry": jij * RY_PER_UNIT[unit],
            "J_mRy": jij * RY_PER_UNIT[unit] * 1.0e3,
            "distance": distance,
            "multiplicity": inline_multiplicity,
            "multiplicity_source": "inline" if has_inline_multiplicity else "default",
            "energy_unit": unit,
            "source": str(source),
        })

    if not records:
        raise ValueError(f"{source}: no Jij records")
    if multiplicities is not None:
        if len(multiplicities) != len(records):
            raise ValueError(
                f"{source}: expected {len(records)} multiplicities, got {len(multiplicities)}"
            )
        for record, multiplicity in zip(records, multiplicities):
            value = float(multiplicity)
            if value <= 0.0:
                raise ValueError(f"{source}: multiplicities must be positive")
            record["multiplicity"] = value
            record["multiplicity_source"] = "explicit"
    else:
        for record in records:
            # ``read_jij_file`` records already carry the correct source tag.
            # Keep this branch explicit for callers that construct compatible
            # records themselves in the future.
            record.setdefault("multiplicity_source", "default")
    return records


def parse_multiplicities(text: str) -> list[float]:
    """Parse ``8,6`` or whitespace-separated shell multiplicities."""
    values = [item for item in re.split(r"[\s,]+", text.strip()) if item]
    if not values:
        raise ValueError("multiplicity list is empty")
    return [float(item) for item in values]


def read_multiplicity_file(path: str | Path) -> list[float]:
    """Read one multiplicity per data line, or ``record_index multiplicity``."""
    values: list[float] = []
    source = Path(path)
    for line_number, raw in enumerate(source.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        try:
            value = _number(fields[-1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{source}:{line_number}: malformed multiplicity") from exc
        if value <= 0.0:
            raise ValueError(f"{source}:{line_number}: multiplicity must be positive")
        values.append(value)
    if not values:
        raise ValueError(f"{source}: no multiplicities")
    return values


def _record_energy(record: dict[str, Any]) -> float:
    if "J_Ry" not in record:
        raise ValueError("Jij record is missing normalized J_Ry")
    return float(record["J_Ry"])


def construct_jq(
    records: Sequence[dict[str, Any]],
    q: Sequence[float],
    *,
    centrosymmetric: bool = True,
) -> complex:
    """Construct scalar one-sublattice ``J(q)`` in Ry.

    With ``centrosymmetric=True`` each row is a shell representative whose
    multiplicity already counts both members of every ``+R/-R`` pair.  Set it
    to ``False`` for a fully directed table, in which case the complex Fourier
    phase is retained.
    """
    q_value = _q_tuple(q)
    result = 0.0j
    for record in records:
        phase = _phase(q_value, record["R"])
        weight = float(record.get("multiplicity", 1.0)) * _record_energy(record)
        result += weight * (math.cos(phase) if centrosymmetric else complex(math.cos(phase), math.sin(phase)))
    return result


def construct_exchange_matrix(
    records: Sequence[dict[str, Any]],
    q: Sequence[float],
    *,
    n_sublattices: int,
    centrosymmetric: bool = False,
) -> list[list[complex]]:
    """Construct directed multi-sublattice ``J_ab(q)`` in Ry.

    Multi-sublattice records are directed: a row ``a b R`` contributes to
    ``J_ab(q)``.  The default therefore retains ``exp(i q.R)`` and expects the
    reverse pair to be present when Hermiticity is required.  The optional
    centrosymmetric form is available only for a table whose rows are already
    shell-summed for each ordered sublattice pair.
    """
    if n_sublattices < 1:
        raise ValueError("n_sublattices must be positive")
    q_value = _q_tuple(q)
    matrix = [[0.0j for _ in range(n_sublattices)] for _ in range(n_sublattices)]
    for record in records:
        a = int(record["species_i"]) - 1
        b = int(record["species_j"]) - 1
        if not (0 <= a < n_sublattices and 0 <= b < n_sublattices):
            raise ValueError(f"species index outside n_sublattices: {record}")
        phase = _phase(q_value, record["R"])
        factor: complex = math.cos(phase) if centrosymmetric else complex(math.cos(phase), math.sin(phase))
        matrix[a][b] += float(record.get("multiplicity", 1.0)) * _record_energy(record) * factor
    return matrix


def _symmetric_eigenvalues(matrix: Sequence[Sequence[float]]) -> list[float]:
    """Small Jacobi eigensolver for the dependency-free matrix path."""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("eigenvalue matrix must be non-empty and square")
    a = [[float(value) for value in row] for row in matrix]
    for _ in range(max(20, 12 * n * n)):
        p, q = 0, 1 if n > 1 else 0
        largest = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > largest:
                    largest = abs(a[i][j])
                    p, q = i, j
        if largest <= 1.0e-14:
            break
        angle = 0.5 * math.atan2(2.0 * a[p][q], a[p][p] - a[q][q])
        c, s = math.cos(angle), math.sin(angle)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        a[p][q] = a[q][p] = 0.0
        for i in range(n):
            if i in (p, q):
                continue
            aip, aiq = a[i][p], a[i][q]
            a[i][p] = a[p][i] = c * aip - s * aiq
            a[i][q] = a[q][i] = s * aip + c * aiq
    return sorted(a[i][i] for i in range(n))


def construct_acoustic_frequencies(
    records: Sequence[dict[str, Any]],
    q: Sequence[float],
    moments: Sequence[float],
    *,
    centrosymmetric: bool = False,
) -> list[float]:
    """Return classical acoustic/optic frequencies in Ry for a multi-sublattice table.

    For ``H = -1/2 sum_(a,b,R) J_ab(R) e_a . e_b`` the quadratic matrix is
    ``K = diag(J(0) 1) - Re J(q)`` and the corresponding classical frequency
    matrix is ``2 M^(-1/2) K M^(-1/2)``.  This reduces to the scalar WP11
    relation ``omega = 2 [J(0)-J(q)] / M``.
    """
    n = len(moments)
    if n < 1 or any(float(moment) <= 0.0 for moment in moments):
        raise ValueError("all sublattice moments must be positive")
    j0 = construct_exchange_matrix(records, (0.0, 0.0, 0.0), n_sublattices=n, centrosymmetric=centrosymmetric)
    jq = construct_exchange_matrix(records, q, n_sublattices=n, centrosymmetric=centrosymmetric)
    row_sum = [sum(j0[a][b].real for b in range(n)) for a in range(n)]
    k_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            raw = (row_sum[a] if a == b else 0.0) - jq[a][b].real
            reverse = (row_sum[b] if b == a else 0.0) - jq[b][a].real
            k_matrix[a][b] = 0.5 * (raw + reverse)
    scaled = [[2.0 * k_matrix[a][b] / math.sqrt(float(moments[a]) * float(moments[b]))
               for b in range(n)] for a in range(n)]
    return _symmetric_eigenvalues(scaled)


def _fit_zero_intercept(q_values: Sequence[float], y_values: Sequence[float], quartic: bool = True) -> dict[str, Any]:
    if len(q_values) != len(y_values) or len(q_values) < (3 if quartic else 1):
        return {"status": "insufficient_points", "n_points": len(q_values)}
    x_values = [float(q) ** 2 for q in q_values]
    if not quartic:
        denominator = sum(x * x for x in x_values)
        if denominator <= 0.0:
            return {"status": "singular"}
        coefficient = sum(x * y for x, y in zip(x_values, y_values)) / denominator
        residuals = [y - coefficient * x for x, y in zip(x_values, y_values)]
        return {"status": "fit", "A": coefficient, "B": 0.0,
                "rmse": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
                "n_points": len(q_values)}

    z_values = [x * x for x in x_values]
    sxx = sum(x * x for x in x_values)
    sxz = sum(x * z for x, z in zip(x_values, z_values))
    szz = sum(z * z for z in z_values)
    syx = sum(y * x for x, y in zip(x_values, y_values))
    syz = sum(y * z for z, y in zip(z_values, y_values))
    determinant = sxx * szz - sxz * sxz
    if abs(determinant) <= 1.0e-30:
        return {"status": "singular", "n_points": len(q_values)}
    coefficient_a = (syx * szz - syz * sxz) / determinant
    coefficient_b = (sxx * syz - sxz * syx) / determinant
    residuals = [y - coefficient_a * x - coefficient_b * z for x, z, y in zip(x_values, z_values, y_values)]
    return {"status": "fit", "A": coefficient_a, "B": coefficient_b,
            "rmse": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
            "n_points": len(q_values)}


def _load_gbt_rows(paths: Iterable[str | Path]) -> list[dict[str, Any]]:
    analyzer = _load_wp09()
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.extend(analyzer.read_sweep_file(path))
    return rows


def _row_omega(row: dict[str, Any]) -> float | None:
    delta = row.get("delta_gauge")
    moment = row.get("mtot")
    sin2 = row.get("sin2theta")
    if delta is not None and moment is not None and sin2 is not None and abs(float(moment) * float(sin2)) > 0.0:
        return 4.0 * float(delta) / (float(moment) * float(sin2))
    value = row.get("omega")
    return None if value is None else float(value)


def _is_zero_q(q: Sequence[float]) -> bool:
    return _norm(q) <= Q_TOLERANCE


def _unique_nonzero_q(rows: Sequence[dict[str, Any]]) -> list[tuple[float, float, float]]:
    result: list[tuple[float, float, float]] = []
    seen: set[tuple[float, float, float]] = set()
    for row in rows:
        q = _q_tuple(row["q"])
        if _is_zero_q(q):
            continue
        key = _q_key(q)
        if key not in seen:
            seen.add(key)
            result.append(q)
    return result


def _q_radius(q: Sequence[float], direction: Sequence[float]) -> float:
    unit = _unit(direction)
    radial = sum(float(a) * float(b) for a, b in zip(q, unit))
    transverse = _norm([float(q[index]) - radial * unit[index] for index in range(3)])
    if transverse > 1.0e-6:
        raise ValueError(f"q={q!r} is not on requested direction {direction!r}")
    return abs(radial)


def _range_convergence(
    exchange_sets: Sequence[dict[str, Any]],
    q_points: Sequence[Sequence[float]],
    *,
    centrosymmetric: bool,
    relative_tolerance: float,
) -> dict[str, Any]:
    if len(exchange_sets) < 2:
        return {
            "status": "not_assessed",
            "reason": "at least two nested LKAG interaction ranges are required",
            "n_ranges": len(exchange_sets),
        }
    if not q_points:
        return {
            "status": "not_assessed",
            "reason": "at least one nonzero GBT q point is required",
            "n_ranges": len(exchange_sets),
        }
    comparisons: list[dict[str, Any]] = []
    for previous, current in zip(exchange_sets, exchange_sets[1:]):
        points: list[dict[str, Any]] = []
        max_relative = 0.0
        for q in q_points:
            old_delta = (construct_jq(previous["records"], (0.0, 0.0, 0.0), centrosymmetric=centrosymmetric).real -
                         construct_jq(previous["records"], q, centrosymmetric=centrosymmetric).real)
            new_delta = (construct_jq(current["records"], (0.0, 0.0, 0.0), centrosymmetric=centrosymmetric).real -
                         construct_jq(current["records"], q, centrosymmetric=centrosymmetric).real)
            scale = max(abs(old_delta), abs(new_delta), 1.0e-12)
            relative = abs(new_delta - old_delta) / scale
            max_relative = max(max_relative, relative)
            points.append({"q": list(_q_tuple(q)), "old_delta_J_Ry": old_delta,
                           "new_delta_J_Ry": new_delta, "relative_difference": relative})
        comparisons.append({
            "from": previous["name"], "to": current["name"],
            "max_relative_difference": max_relative,
            "tolerance": relative_tolerance,
            "status": "converged" if max_relative <= relative_tolerance else "not_converged",
            "points": points,
        })
    return {"status": "converged" if all(item["status"] == "converged" for item in comparisons) else "not_converged",
            "n_ranges": len(exchange_sets), "comparisons": comparisons,
            "tolerance": relative_tolerance}


def _wp10_summary(analysis: dict[str, Any] | None, tolerance: float) -> dict[str, Any]:
    if analysis is None:
        return {"status": "not_assessed", "reason": "no WP10 analysis JSON supplied"}
    sensitivity = analysis.get("sensitivity", {}).get("sin2", [])
    mft = [item for item in sensitivity if item.get("mode") == "mft"]
    if not mft:
        return {"status": "not_assessed", "reason": "WP10 JSON has no bare-MFT sin2 sensitivity summary"}
    spread = float(mft[0].get("A_relative_spread", math.inf))
    odd_statuses = [item.get("reversal", {}).get("odd_component", {}).get("status")
                    for item in analysis.get("groups", []) if item.get("mode") == "mft"]
    odd_pass = bool(odd_statuses) and all(status == "within_tolerance" for status in odd_statuses)
    status = "converged" if spread <= tolerance else "not_converged"
    return {"status": status, "A_relative_spread": spread, "tolerance": tolerance,
            "q_reversal_status": "pass" if odd_pass else "not_passed",
            "source": "WP10 analysis JSON"}


def analyze_wp11(
    exchange_sets: Sequence[dict[str, Any]],
    gbt_rows: Sequence[dict[str, Any]] | None = None,
    *,
    direction: Sequence[float] = (0.0, 0.0, 1.0),
    moment: float | None = None,
    n_sublattices: int = 1,
    centrosymmetric: bool = True,
    wp10_analysis: dict[str, Any] | None = None,
    wp10_relative_tolerance: float = DEFAULT_CURVATURE_RELATIVE_TOLERANCE,
    range_relative_tolerance: float = DEFAULT_RANGE_RELATIVE_TOLERANCE,
    curvature_relative_tolerance: float = DEFAULT_CURVATURE_RELATIVE_TOLERANCE,
    goldstone_absolute_tolerance: float = DEFAULT_GOLDSTONE_ABSOLUTE_TOLERANCE,
) -> dict[str, Any]:
    """Build the WP11 evidence record from normalized exchange and GBT rows."""
    if not exchange_sets:
        raise ValueError("at least one exchange set is required")
    direction = _unit(direction)
    rows = list(gbt_rows or [])
    q_points = _unique_nonzero_q(rows)

    sets_payload: list[dict[str, Any]] = []
    for exchange_set in exchange_sets:
        records = exchange_set["records"]
        j0 = construct_jq(records, (0.0, 0.0, 0.0), centrosymmetric=centrosymmetric).real
        entry: dict[str, Any] = {
            "name": exchange_set["name"], "n_records": len(records),
            "multiplicities": [float(record.get("multiplicity", 1.0)) for record in records],
            "J0_Ry": j0, "J0_mRy": 1.0e3 * j0,
            "records": [{"index": int(record["index"]), "species_i": int(record["species_i"]),
                         "species_j": int(record["species_j"]), "R": list(record["R"]),
                         "J_mRy": float(record["J_mRy"]), "distance": float(record["distance"]),
                         "multiplicity": float(record.get("multiplicity", 1.0))} for record in records],
        }
        if n_sublattices == 1:
            entry["Jq"] = []
            for q in q_points:
                jq = construct_jq(records, q, centrosymmetric=centrosymmetric).real
                entry["Jq"].append({"q": list(q), "Jq_Ry": jq, "J0_minus_Jq_Ry": j0 - jq,
                                    "q_radius": _q_radius(q, direction)})
        sets_payload.append(entry)

    comparison: list[dict[str, Any]] = []
    gbt_groups: dict[tuple[str, tuple[int, int, int] | None, float | None], list[dict[str, Any]]] = {}
    for row in rows:
        key = (str(row.get("mode", "unknown")), row.get("mesh"), row.get("theta_deg"))
        gbt_groups.setdefault(key, []).append(row)
    for key, group in sorted(gbt_groups.items(), key=lambda item: str(item[0])):
        if n_sublattices != 1:
            continue
        for row in group:
            q = _q_tuple(row["q"])
            gbt_omega = _row_omega(row)
            if gbt_omega is None:
                continue
            radius = 0.0 if _is_zero_q(q) else _q_radius(q, direction)
            for exchange_set, payload in zip(exchange_sets, sets_payload):
                jq = construct_jq(exchange_set["records"], q, centrosymmetric=centrosymmetric).real
                j0 = payload["J0_Ry"]
                current_moment = float(moment if moment is not None else row.get("mtot") or 0.0)
                if current_moment <= 0.0:
                    raise ValueError("a positive total moment is required for scalar LKAG comparison")
                expected = 2.0 * (j0 - jq) / current_moment
                residual = gbt_omega - expected
                sin2 = row.get("sin2theta")
                if sin2 is not None:
                    gbt_delta = float(row.get("delta_gauge")) if row.get("delta_gauge") is not None else (
                        gbt_omega * current_moment * float(sin2) / 4.0
                    )
                    lkag_delta = 0.5 * (j0 - jq) * float(sin2)
                    delta_residual = gbt_delta - lkag_delta
                    relative_delta_residual = abs(delta_residual) / max(
                        abs(gbt_delta), abs(lkag_delta), 1.0e-12
                    )
                else:
                    gbt_delta = None
                    lkag_delta = None
                    delta_residual = None
                    relative_delta_residual = None
                comparison.append({
                    "mode": key[0], "mesh": list(key[1]) if key[1] is not None else None,
                    "theta_deg": key[2], "q": list(q), "q_radius": radius,
                    "gbt_omega_Ry": gbt_omega, "lkag_omega_Ry": expected,
                    "residual_Ry": residual,
                    "relative_residual": abs(residual) / max(abs(gbt_omega), abs(expected), 1.0e-8),
                    "gbt_delta_Ry": gbt_delta,
                    "lkag_delta_Ry": lkag_delta,
                    "delta_residual_Ry": delta_residual,
                    "relative_delta_residual": relative_delta_residual,
                    "exchange_set": exchange_set["name"],
                })

    curvature: list[dict[str, Any]] = []
    if n_sublattices == 1:
        for key, group in sorted(gbt_groups.items(), key=lambda item: str(item[0])):
            paired: dict[float, list[float]] = {}
            for row in group:
                omega = _row_omega(row)
                if omega is None or _is_zero_q(row["q"]):
                    continue
                radius = _q_radius(row["q"], direction)
                paired.setdefault(round(radius, 8), []).append(omega)
            q_fit = sorted(paired)
            y_fit = [sum(values) / len(values) for radius, values in sorted(paired.items())]
            if len(q_fit) >= 2:
                gbt_fit = _fit_zero_intercept(q_fit, y_fit, quartic=len(q_fit) >= 3)
            else:
                gbt_fit = {"status": "insufficient_points", "n_points": len(q_fit)}
            curvature.append({"mode": key[0], "mesh": list(key[1]) if key[1] is not None else None,
                              "theta_deg": key[2], "q_points": q_fit, "gbt_omega_fit": gbt_fit})
            for exchange_set in exchange_sets:
                j_values = []
                j_q = []
                for radius in q_fit:
                    q = tuple(direction[index] * radius for index in range(3))
                    j0 = construct_jq(exchange_set["records"], (0.0, 0.0, 0.0), centrosymmetric=centrosymmetric).real
                    jq = construct_jq(exchange_set["records"], q, centrosymmetric=centrosymmetric).real
                    current_moment = float(moment if moment is not None else (group[0].get("mtot") or 0.0))
                    j_q.append(2.0 * (j0 - jq) / current_moment)
                    j_values.append(radius)
                lkag_fit = _fit_zero_intercept(j_values, j_q, quartic=len(j_values) >= 3)
                item = curvature[-1]
                item.setdefault("lkag_omega_fits", []).append({"exchange_set": exchange_set["name"], "fit": lkag_fit})
                if gbt_fit.get("status") == "fit" and lkag_fit.get("status") == "fit":
                    gbt_a = float(gbt_fit["A"])
                    lkag_a = float(lkag_fit["A"])
                    relative = abs(gbt_a - lkag_a) / max(abs(gbt_a), abs(lkag_a), 1.0e-12)
                    item.setdefault("curvature_comparisons", []).append({
                        "exchange_set": exchange_set["name"], "gbt_A_Ry": gbt_a,
                        "lkag_A_Ry": lkag_a, "relative_difference": relative,
                        "tolerance": curvature_relative_tolerance,
                        "status": "within_tolerance" if relative <= curvature_relative_tolerance else "not_within_tolerance",
                    })

    gamma_rows = [item for item in comparison if _is_zero_q(item["q"])]
    max_gamma = max((abs(item["gbt_omega_Ry"]) for item in gamma_rows), default=None)
    curve_max_relative = max((item["relative_residual"] for item in comparison if not _is_zero_q(item["q"])), default=None)
    curve_status = "not_assessed" if curve_max_relative is None else ("within_tolerance" if curve_max_relative <= curvature_relative_tolerance else "not_within_tolerance")
    curvature_comparisons = [item for item in curvature for item in item.get("curvature_comparisons", [])]
    curvature_status = "not_assessed" if not curvature_comparisons else ("within_tolerance" if all(item["status"] == "within_tolerance" for item in curvature_comparisons) else "not_within_tolerance")
    range_status = _range_convergence(exchange_sets, q_points, centrosymmetric=centrosymmetric,
                                      relative_tolerance=range_relative_tolerance)
    wp10_status = _wp10_summary(wp10_analysis, wp10_relative_tolerance)
    goldstone_status = "not_assessed" if max_gamma is None else ("within_tolerance" if max_gamma <= goldstone_absolute_tolerance else "not_within_tolerance")
    acceptance = {
        "wp10_dependency": wp10_status,
        "exchange_range": range_status,
        "q0_goldstone": {"status": goldstone_status, "max_abs_omega_Ry": max_gamma,
                         "tolerance_Ry": goldstone_absolute_tolerance},
        "curve": {"status": curve_status, "max_relative_residual": curve_max_relative,
                   "tolerance": curvature_relative_tolerance},
        "curvature": {"status": curvature_status, "comparisons": curvature_comparisons,
                       "tolerance": curvature_relative_tolerance},
    }
    required = [wp10_status["status"], range_status["status"], goldstone_status, curve_status, curvature_status]
    overall = "pass" if all(status in {"converged", "within_tolerance"} for status in required) else "not_qualified"
    return {
        "schema": "gbt_wp11_lkag_jq_v1",
        "convention": {
            "exchange_hamiltonian": "H = -1/2 sum_(i,R) J_R e_i dot e_(i+R)",
            "positive_J": "ferromagnetic",
            "pair_counting": "directed R table; 1/2 removes +/- double counting",
            "frozen_magnon_formula": "omega_GBT = 4 DeltaE / (M sin(theta)^2)",
            "derived_relation": "omega_LKAG = 2 [J(0)-J(q)] / M",
            "energy_unit": "Ry internally; jij.out is written in mRy",
            "q_unit": "Cartesian 2*pi/alat; R from jij.out in alat units",
            "centrosymmetric": bool(centrosymmetric),
        },
        "direction_internal": list(direction),
        "n_sublattices": n_sublattices,
        "exchange_sets": sets_payload,
        "comparison_count": len(comparison),
        "comparisons": comparison,
        "curvature": curvature,
        "acceptance": acceptance,
        "overall_status": overall,
    }


def format_report(result: dict[str, Any]) -> str:
    acceptance = result["acceptance"]
    lines = [
        "GBT WP11 LKAG J(q) analysis",
        f"status: {result['overall_status']}",
        f"exchange sets: {', '.join(item['name'] for item in result['exchange_sets'])}",
        f"comparison rows: {result['comparison_count']}",
        "",
        "acceptance:",
    ]
    for name, item in acceptance.items():
        lines.append(f"  {name}: {item.get('status')}")
    for item in result["exchange_sets"]:
        lines.append(f"  {item['name']}: J(0) = {item['J0_mRy']:.9g} mRy, records = {item['n_records']}")
    if acceptance["curve"].get("max_relative_residual") is not None:
        lines.append(f"  max curve relative residual: {acceptance['curve']['max_relative_residual']:.6g}")
    return "\n".join(lines)


def _parse_q_point(value: str) -> tuple[float, float, float]:
    fields = [item for item in re.split(r"[\s,]+", value.strip()) if item]
    if len(fields) != 3:
        raise ValueError(f"q-point must have three values, got {value!r}")
    return _q_tuple([_number(item) for item in fields])


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jij", action="append", required=True, help="production jij.out; repeat for nested ranges")
    parser.add_argument("--sweep", action="append", help="WP09/WP10 sweep CSV or frozen-magnon output")
    parser.add_argument("--wp10-analysis", type=Path, help="WP10 analysis JSON used as a convergence dependency")
    parser.add_argument("--multiplicities", action="append", help="one comma-separated multiplicity list per --jij")
    parser.add_argument("--multiplicity-file", action="append", type=Path, help="one multiplicity file per --jij")
    parser.add_argument("--jij-name", action="append", help="name one exchange set; defaults to filename stem")
    parser.add_argument("--energy-unit", default="mry", choices=sorted(RY_PER_UNIT))
    parser.add_argument("--bond-scale", type=float, default=1.0)
    parser.add_argument("--direction", nargs=3, type=float, default=(0.0, 0.0, 1.0), metavar=("QX", "QY", "QZ"))
    parser.add_argument("--q-point", action="append", help="q point if no sweep is supplied")
    parser.add_argument("--moment", type=float, help="fixed one-sublattice moment in mu_B")
    parser.add_argument("--n-sublattices", type=int, default=1)
    parser.add_argument("--directed", action="store_true", help="retain exp(i q.R) instead of centrosymmetric cos")
    parser.add_argument("--json", action="store_true", help="print complete JSON")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.multiplicities and len(args.multiplicities) not in {1, len(args.jij)}:
            raise ValueError("--multiplicities must be supplied once or once per --jij")
        if args.multiplicity_file and len(args.multiplicity_file) not in {1, len(args.jij)}:
            raise ValueError("--multiplicity-file must be supplied once or once per --jij")
        if args.jij_name and len(args.jij_name) not in {1, len(args.jij)}:
            raise ValueError("--jij-name must be supplied once or once per --jij")
        sets: list[dict[str, Any]] = []
        for index, path in enumerate(args.jij):
            multiplicities = None
            if args.multiplicities:
                value = args.multiplicities[0 if len(args.multiplicities) == 1 else index]
                multiplicities = parse_multiplicities(value)
            if args.multiplicity_file:
                values = read_multiplicity_file(args.multiplicity_file[0 if len(args.multiplicity_file) == 1 else index])
                if multiplicities is not None:
                    raise ValueError("use either --multiplicities or --multiplicity-file, not both")
                multiplicities = values
            records = read_jij_file(path, multiplicities=multiplicities,
                                    energy_unit=args.energy_unit, bond_scale=args.bond_scale)
            name = args.jij_name[0 if len(args.jij_name) == 1 else index] if args.jij_name else Path(path).stem
            sets.append({"name": name, "records": records})
        rows = _load_gbt_rows(args.sweep or []) if args.sweep else []
        if not rows and args.q_point:
            rows = [{"q": _parse_q_point(value), "mode": "synthetic", "mtot": args.moment,
                     "omega": None, "delta_gauge": None, "sin2theta": None} for value in args.q_point]
        wp10 = json.loads(args.wp10_analysis.read_text()) if args.wp10_analysis else None
        result = analyze_wp11(sets, rows, direction=args.direction, moment=args.moment,
                              n_sublattices=args.n_sublattices, centrosymmetric=not args.directed,
                              wp10_analysis=wp10)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(payload + "\n")
    print(payload if args.json else format_report(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
