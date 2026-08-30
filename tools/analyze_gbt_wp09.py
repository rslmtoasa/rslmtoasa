#!/usr/bin/env python3
"""Analyse raw GBT cone-angle and k-grid sweep files.

The utility deliberately does not fit a dispersion or apply an energy shift.
It reads the raw per-q records emitted by
``frozen_magnon_harmonic_diagnostics.dat`` and reports the quantities needed
by WP09:

* raw and same-q gauge-subtracted energies;
* ``DeltaE / sin(theta)**2`` plateau candidates;
* convergence of each energy definition as the k mesh is refined;
* pure-gauge drift, Fermi/electron-count consistency, and q/2 mesh mapping.

The preferred input schema is the whitespace-delimited file whose header
contains ``# columns = ...``.  A CSV/TSV file with the same column names is
also accepted, which makes it convenient to collate separate MFT and
constrained-MFT runs.  The old ``frozen_magnon_diagnostics.dat`` and
``frozen_magnon.dat`` layouts are accepted as a migration path, but they
cannot provide fields that those historical files never wrote.

Example::

    python tools/analyze_gbt_wp09.py run_nk8 run_nk12 run_nk16 \
        --json-out wp09_analysis.json

This module is intentionally standard-library-only so it can run in the lean
unit-test environment as well as on a production calculation directory.

Mesh convergence is reported as ``converged`` when the full mesh spread is
within either the default absolute tolerance of 1e-6 Ry or the default
relative tolerance of 1 percent.  The tolerances are command-line options so
the evidence report can state exactly which gate was applied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence


MISSING = {"", "-", "--", "na", "n/a", "none", "nan", "null"}
Q_TOLERANCE = 1.0e-8
ANGLE_TOLERANCE = 1.0e-8
KGRID_ABSOLUTE_TOLERANCE = 1.0e-6
KGRID_RELATIVE_TOLERANCE = 1.0e-2
DEFAULT_STABILITY_RELATIVE_TOLERANCE = 0.05


def _number(value: str | float | int | None) -> float | None:
    if value is None:
        return None
    if isinstance(value, (float, int)):
        return float(value)
    text = value.strip()
    if text.lower() in MISSING:
        return None
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(f"expected a number, got {value!r}") from exc


def _boolean(value: str | float | int | bool | None) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"t", ".true.", "true", "yes", "y", "1", "1.0"}:
        return True
    if text in {"f", ".false.", "false", "no", "n", "0", "0.0"}:
        return False
    return None


def _normal_name(value: str) -> str:
    """Normalize schema names without changing their semantic spelling."""
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def _q_key(q: Sequence[float]) -> tuple[float, float, float]:
    return tuple(round(float(value), 8) for value in q)  # type: ignore[return-value]


def _mesh_key(mesh: Sequence[int] | None) -> tuple[int, int, int] | None:
    if mesh is None:
        return None
    return tuple(int(value) for value in mesh)  # type: ignore[return-value]


def _parse_metadata(lines: Iterable[str]) -> tuple[dict[str, str], list[str] | None, list[str]]:
    metadata: dict[str, str] = {}
    columns: list[str] | None = None
    data_lines: list[str] = []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            comment = line[1:].strip()
            if "=" in comment:
                key, value = comment.split("=", 1)
                key = _normal_name(key)
                metadata[key] = value.strip()
                if key == "columns":
                    columns = _split_columns(value)
            continue
        data_lines.append(raw.rstrip("\n"))

    if data_lines and columns is None:
        candidate = data_lines[0].strip()
        tokens = _split_columns(candidate)
        # A header is recognized by a q/energy field name, not merely by the
        # absence of a successful float conversion.  This keeps malformed
        # numeric rows from being silently discarded as headers.
        if any(_normal_name(token) in {"q1", "q_1", "e_q_theta", "mode"} for token in tokens):
            columns = tokens
            data_lines = data_lines[1:]

    return metadata, columns, data_lines


def _split_columns(line: str) -> list[str]:
    if "," in line:
        return [item.strip() for item in next(csv.reader([line]))]
    if "\t" in line:
        return [item.strip() for item in line.split("\t")]
    return line.split()


def _metadata_float(metadata: dict[str, str], *keys: str) -> float | None:
    for key in keys:
        value = metadata.get(_normal_name(key))
        if value is not None:
            return _number(value)
    return None


def _metadata_mesh(metadata: dict[str, str]) -> tuple[int, int, int] | None:
    values: list[int] = []
    for axis in ("nk1", "nk2", "nk3"):
        value = _metadata_float(metadata, axis)
        if value is None:
            values = []
            break
        values.append(int(round(value)))
    if len(values) == 3:
        return _mesh_key(values)
    for key in ("k_mesh", "mesh", "kgrid", "k_grid"):
        raw = metadata.get(_normal_name(key))
        if raw:
            numbers = re.findall(r"[-+]?\d+", raw)
            if len(numbers) == 3:
                return _mesh_key([int(item) for item in numbers])
    return None


def _field_index(columns: Sequence[str] | None) -> dict[str, int]:
    if columns is None:
        return {}
    aliases = {
        "q_1": "q1", "q_2": "q2", "q_3": "q3",
        "theta": "theta_deg", "theta_ss": "theta_deg",
        "e_q": "e_q_theta", "e_q_theta_": "e_q_theta",
        "e_0": "e_qref_theta", "e_ref_theta": "e_qref_theta",
        "e_q_zero": "e_q0", "e_q_theta_zero": "e_q0",
        "e_ref_zero": "e_qref0", "e_q_reference_zero": "e_qref0",
        "delta_e_raw": "delta_raw", "delta_e_gauge": "delta_gauge",
        "delta_e_pure": "delta_pure", "sin2": "sin2theta",
        "m_tot": "mtot", "n": "electron_count", "electrons": "electron_count",
        "electron_number": "electron_count", "dn": "electron_error",
        "target_n": "target_electrons", "fermi": "fermi_level",
        "ef": "fermi_level", "weight": "weight_sum",
        "gauge": "gauge_available", "q_half_maps": "q_half_maps_to_mesh",
        "q_half_mapping": "q_half_maps_to_mesh",
    }
    result: dict[str, int] = {}
    for index, column in enumerate(columns):
        name = _normal_name(column)
        result[aliases.get(name, name)] = index
    return result


def _raw_value(tokens: Sequence[str], fields: dict[str, int], name: str) -> float | str | None:
    index = fields.get(name)
    if index is None or index >= len(tokens):
        return None
    value = tokens[index]
    if name == "mode":
        return value
    return _number(value)


def _legacy_fields(path: Path, tokens: Sequence[str]) -> dict[str, float | str | None]:
    """Map the two pre-WP09 output layouts into the common row schema."""
    values: dict[str, float | str | None] = {}
    if path.name == "frozen_magnon_diagnostics.dat":
        # q1 q2 q3 E(q,theta) E(qref,theta) DeltaE_raw E(q,0)
        # DeltaE_gauge sin2theta Mtot omega
        if len(tokens) < 11:
            raise ValueError(f"{path}: expected at least 11 legacy diagnostic columns")
        values.update({
            "q1": _number(tokens[0]), "q2": _number(tokens[1]), "q3": _number(tokens[2]),
            "e_q_theta": _number(tokens[3]), "e_qref_theta": _number(tokens[4]),
            "delta_raw": _number(tokens[5]), "e_q0": _number(tokens[6]),
            "delta_gauge": _number(tokens[7]), "sin2theta": _number(tokens[8]),
            "mtot": _number(tokens[9]), "omega": _number(tokens[10]),
            "mode": "mft", "gauge_available": True,
        })
        return values
    if path.name == "frozen_magnon.dat":
        # q1 q2 q3 etot eband mtot_1..mtot_nrec omega.  The first row is the
        # reference row; its reference energy is filled after all rows load.
        if len(tokens) < 6:
            raise ValueError(f"{path}: expected a frozen_magnon.dat data row")
        values.update({
            "q1": _number(tokens[0]), "q2": _number(tokens[1]), "q3": _number(tokens[2]),
            "e_q_theta": _number(tokens[4]), "mtot": _number(tokens[5]),
            "omega": _number(tokens[-1]), "mode": "mft",
        })
        return values
    raise ValueError(f"{path}: no column header and no supported legacy filename")


def read_sweep_file(
    path: str | Path,
    *,
    mode: str | None = None,
    theta_deg: float | None = None,
    mesh: Sequence[int] | None = None,
    target_electrons: float | None = None,
) -> list[dict[str, Any]]:
    """Read one WP09 sweep file into normalized per-q dictionaries.

    ``mode``, ``theta_deg``, ``mesh``, and ``target_electrons`` are explicit
    overrides for old files that predate the WP09 metadata header.
    """
    source = Path(path)
    metadata, columns, data_lines = _parse_metadata(source.read_text().splitlines())
    fields = _field_index(columns)
    file_mode = mode or metadata.get("mode")
    file_theta = theta_deg if theta_deg is not None else _metadata_float(metadata, "theta_deg", "theta")
    file_mesh = _mesh_key(mesh) if mesh is not None else _metadata_mesh(metadata)
    file_alat = _metadata_float(metadata, "alat_angstrom", "alat")
    target = target_electrons if target_electrons is not None else _metadata_float(
        metadata, "target_electrons", "total_electrons", "electron_target"
    )

    rows: list[dict[str, Any]] = []
    for line_number, line in enumerate(data_lines, 1):
        tokens = _split_columns(line.strip())
        try:
            values = _legacy_fields(source, tokens) if not fields else {
                name: (_raw_value(tokens, fields, name)) for name in fields
            }
        except ValueError as exc:
            raise ValueError(f"{source}:{line_number}: {exc}") from exc
        if not values:
            continue
        try:
            q = tuple(_number(values.get(axis)) for axis in ("q1", "q2", "q3"))
        except ValueError as exc:
            raise ValueError(f"{source}:{line_number}: malformed q row") from exc
        if any(value is None for value in q):
            raise ValueError(f"{source}:{line_number}: q1, q2, q3 are required")
        row_mesh = file_mesh
        if row_mesh is None:
            row_mesh_values = [_number(values.get(axis)) for axis in ("nk1", "nk2", "nk3")]
            if all(value is not None for value in row_mesh_values):
                row_mesh = _mesh_key([int(round(float(value))) for value in row_mesh_values])
        row: dict[str, Any] = {
            "source": str(source),
            "line": line_number,
            "mode": str(file_mode or values.get("mode") or "unknown").strip().lower(),
            "q": tuple(float(value) for value in q),
            "theta_deg": _number(values.get("theta_deg")) if values.get("theta_deg") is not None else file_theta,
            "mesh": row_mesh,
            "alat_angstrom": file_alat,
        }
        for name in (
            "e_q_theta", "e_q0", "e_qref_theta", "e_qref0", "delta_raw", "delta_gauge",
            "delta_pure", "sin2theta", "mtot", "omega", "fermi_level", "electron_count",
            "electron_error", "target_electrons", "weight_sum", "nk_total",
        ):
            value = values.get(name)
            row[name] = _number(value) if value is not None else None
        row["target_electrons"] = row["target_electrons"] if row["target_electrons"] is not None else target
        gauge = _boolean(values.get("gauge_available"))
        row["gauge_available"] = gauge
        mapping = _boolean(values.get("q_half_maps_to_mesh"))
        row["q_half_maps_to_mesh"] = mapping
        row["q_half_mapping_source"] = "file_metadata" if mapping is not None else None
        if row["theta_deg"] is None and row["sin2theta"] is not None:
            sin2 = float(row["sin2theta"])
            if 0.0 <= sin2 <= 1.0:
                row["theta_deg"] = math.degrees(math.asin(math.sqrt(sin2)))
        rows.append(row)

    if not rows:
        raise ValueError(f"{source}: no data rows")

    # The reference row is a property of the sweep file, not of a parser
    # convention.  Fill omitted reference fields from row 1 and derive the
    # raw/gauge/pure differences from the retained raw energies.
    reference = rows[0]
    for row in rows:
        if row["e_qref_theta"] is None:
            row["e_qref_theta"] = reference["e_q_theta"]
        if row["e_qref0"] is None:
            row["e_qref0"] = reference["e_q0"]
        if row["delta_raw"] is None and row["e_q_theta"] is not None and row["e_qref_theta"] is not None:
            row["delta_raw"] = row["e_q_theta"] - row["e_qref_theta"]
        if row["gauge_available"] is not False:
            if row["delta_gauge"] is None and row["e_q_theta"] is not None and row["e_q0"] is not None:
                row["delta_gauge"] = row["e_q_theta"] - row["e_q0"]
            if row["delta_pure"] is None and row["e_q0"] is not None and row["e_qref0"] is not None:
                row["delta_pure"] = row["e_q0"] - row["e_qref0"]
        if row["sin2theta"] is None and row["theta_deg"] is not None:
            row["sin2theta"] = math.sin(math.radians(float(row["theta_deg"]))) ** 2
        if row["mesh"] is None:
            row["mesh"] = _mesh_key(mesh)
        if row["q_half_maps_to_mesh"] is None:
            row["q_half_maps_to_mesh"] = infer_q_half_mesh_mapping(row["q"], row["mesh"])
            if row["q_half_maps_to_mesh"] is not None:
                row["q_half_mapping_source"] = "inferred_uniform_cartesian"
    return rows


def infer_q_half_mesh_mapping(q: Sequence[float], mesh: Sequence[int] | None, tolerance: float = 1.0e-8) -> bool | None:
    """Infer q/2 mesh closure for a uniform Cartesian 2*pi/alat mesh.

    This is exact for the cubic/Cartesian fixtures used by WP09.  For a
    non-cubic lattice callers should provide ``q_half_maps_to_mesh`` from a
    lattice-aware driver rather than relying on this fallback.
    """
    if mesh is None:
        return None
    return all(abs(float(value) * int(nk) / 2.0 - round(float(value) * int(nk) / 2.0)) <= tolerance
               for value, nk in zip(q, mesh))


def _group(rows: Iterable[dict[str, Any]], keys: Sequence[str]) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        values: list[Any] = []
        for key in keys:
            value = row[key]
            if key == "q":
                value = _q_key(value)
            elif key == "mesh":
                value = _mesh_key(value)
            elif key == "theta_deg" and value is not None:
                value = round(float(value), 8)
            values.append(value)
        grouped[tuple(values)].append(row)
    return grouped


def _spread(values: Sequence[float], energy_floor: float = 1.0e-30) -> tuple[float, float, float]:
    if not values:
        return math.nan, math.nan, math.nan
    mean = sum(values) / len(values)
    absolute = max(values) - min(values)
    relative = absolute / max(abs(mean), energy_floor)
    return absolute, relative, mean


def _definition_value(row: dict[str, Any], definition: str) -> float | None:
    if definition == "raw":
        delta = row.get("delta_raw")
    elif definition == "gauge":
        delta = row.get("delta_gauge")
    elif definition == "pure":
        delta = row.get("delta_pure")
    else:
        raise ValueError(f"unknown energy definition {definition!r}")
    if delta is None:
        return None
    return float(delta)


def _plateau_for_group(rows: list[dict[str, Any]], definition: str, relative_tolerance: float,
                       min_points: int, energy_floor: float) -> dict[str, Any] | None:
    points: list[tuple[float, float]] = []
    for row in rows:
        angle = row.get("theta_deg")
        sin2 = row.get("sin2theta")
        delta = _definition_value(row, definition)
        if angle is None or sin2 is None or delta is None or float(sin2) <= 0.0:
            continue
        points.append((float(angle), delta / float(sin2)))
    by_angle: dict[float, float] = {}
    for angle, value in points:
        by_angle.setdefault(round(angle, 8), value)
    ordered = sorted(by_angle.items())
    if len(ordered) < min_points:
        return None

    candidates: list[tuple[int, float, float, int, int]] = []
    for start in range(len(ordered)):
        for end in range(start + min_points, len(ordered) + 1):
            values = [value for _, value in ordered[start:end]]
            absolute, relative, _ = _spread(values, energy_floor)
            if relative <= relative_tolerance:
                # Longest window first, then smallest spread, then the
                # lowest maximum angle.  The tie-break makes the admitted
                # range deterministic when all points happen to be flat.
                candidates.append((end - start, relative, absolute, start, end))
    if not candidates:
        return {
            "angles_deg": [angle for angle, _ in ordered],
            "admitted_angles_deg": [],
            "excluded_angles_deg": [angle for angle, _ in ordered],
            "spread_abs": None,
            "spread_rel": None,
            "status": "no_plateau",
        }
    _, _, _, start, end = min(candidates, key=lambda item: (-item[0], item[1], item[2], ordered[item[4] - 1][0]))
    admitted = ordered[start:end]
    values = [value for _, value in admitted]
    absolute, relative, mean = _spread(values, energy_floor)
    admitted_angles = [angle for angle, _ in admitted]
    return {
        "angles_deg": [angle for angle, _ in ordered],
        "admitted_angles_deg": admitted_angles,
        "excluded_angles_deg": [angle for angle, _ in ordered if angle not in admitted_angles],
        "spread_abs": absolute,
        "spread_rel": relative,
        "mean_k": mean,
        "status": "plateau",
    }


def analyze_rows(rows: Sequence[dict[str, Any]], relative_tolerance: float = 0.05,
                 min_plateau_points: int = 2, energy_floor: float = 1.0e-30,
                 kgrid_absolute_tolerance: float = KGRID_ABSOLUTE_TOLERANCE,
                 kgrid_relative_tolerance: float = KGRID_RELATIVE_TOLERANCE) -> dict[str, Any]:
    """Return JSON-serializable WP09 analysis results for normalized rows."""
    rows = [dict(row) for row in rows]
    for row in rows:
        for definition in ("raw", "gauge", "pure"):
            delta = _definition_value(row, definition)
            sin2 = row.get("sin2theta")
            row[f"k_{definition}"] = None if delta is None or not sin2 else delta / float(sin2)

    plateaus: list[dict[str, Any]] = []
    for (mode, q, mesh), group in _group(rows, ("mode", "q", "mesh")).items():
        definitions = ["gauge", "raw"] if any(_definition_value(item, "gauge") is not None for item in group) else ["raw"]
        for definition in definitions:
            result = _plateau_for_group(group, definition, relative_tolerance, min_plateau_points, energy_floor)
            if result is not None:
                plateaus.append({"mode": mode, "q": list(q), "mesh": list(mesh) if mesh else None,
                                 "definition": definition, **result})

    kgrid: list[dict[str, Any]] = []
    for mode in sorted({row["mode"] for row in rows}):
        for q in sorted({_q_key(row["q"]) for row in rows if row["mode"] == mode}):
            for theta in sorted({row.get("theta_deg") for row in rows if row["mode"] == mode and _q_key(row["q"]) == q and row.get("theta_deg") is not None}):
                selected = [row for row in rows if row["mode"] == mode and _q_key(row["q"]) == q and row.get("theta_deg") == theta]
                for definition in ("raw", "gauge", "pure"):
                    by_mesh: dict[tuple[int, int, int] | None, float] = {}
                    for row in selected:
                        value = _definition_value(row, definition)
                        if value is None:
                            continue
                        mesh = _mesh_key(row.get("mesh"))
                        by_mesh.setdefault(mesh, value)
                    if not by_mesh:
                        continue
                    ordered = sorted(by_mesh.items(), key=lambda item: (math.prod(item[0]) if item[0] else -1, item[0] or ()))
                    values = [value for _, value in ordered]
                    absolute, relative, mean = _spread(values, energy_floor)
                    successive = [abs(values[index] - values[index - 1]) for index in range(1, len(values))]
                    if len(values) == 1:
                        status = "single_mesh"
                    elif absolute <= kgrid_absolute_tolerance or relative <= kgrid_relative_tolerance:
                        status = "converged"
                    else:
                        status = "not_converged"
                    kgrid.append({
                        "mode": mode, "q": list(q), "theta_deg": theta, "definition": definition,
                        "meshes": [{"mesh": list(mesh) if mesh else None, "delta": value} for mesh, value in ordered],
                        "spread_abs": absolute, "spread_rel": relative, "mean_delta": mean,
                        "successive_delta_abs": successive,
                        "status": status,
                    })

    pure_drift: list[dict[str, Any]] = []
    for (mode, theta, mesh), group in _group(rows, ("mode", "theta_deg", "mesh")).items():
        gamma = next((row for row in group if max(abs(value) for value in row["q"]) <= Q_TOLERANCE), None)
        for row in group:
            if row.get("e_q0") is None:
                continue
            reference_zero = gamma.get("e_q0") if gamma else row.get("e_qref0")
            if reference_zero is None:
                continue
            pure_drift.append({
                "mode": mode, "theta_deg": theta, "mesh": list(mesh) if mesh else None,
                "q": list(_q_key(row["q"])), "drift": float(row["e_q0"]) - float(reference_zero),
                "reference": "gamma" if gamma else "q_ref",
            })

    electron_consistency: list[dict[str, Any]] = []
    for (mode, q, theta), group in _group(rows, ("mode", "q", "theta_deg")).items():
        errors: list[float] = []
        fermi: list[float] = []
        for row in group:
            error = row.get("electron_error")
            if error is None and row.get("electron_count") is not None and row.get("target_electrons") is not None:
                error = float(row["electron_count"]) - float(row["target_electrons"])
            if error is not None:
                errors.append(abs(float(error)))
            if row.get("fermi_level") is not None:
                fermi.append(float(row["fermi_level"]))
        if errors or fermi:
            electron_consistency.append({
                "mode": mode, "q": list(q), "theta_deg": theta,
                "max_abs_electron_error": max(errors) if errors else None,
                "fermi_min": min(fermi) if fermi else None,
                "fermi_max": max(fermi) if fermi else None,
                "samples": len(group),
            })

    mapping: list[dict[str, Any]] = []
    for row in rows:
        mapping.append({
            "mode": row["mode"], "q": list(_q_key(row["q"])), "theta_deg": row.get("theta_deg"),
            "mesh": list(row["mesh"]) if row.get("mesh") else None,
            "q_half_maps_to_mesh": row.get("q_half_maps_to_mesh"),
            "method": row.get("q_half_mapping_source") or "unavailable",
        })

    return {
        "schema": "gbt_wp09_analysis_v1",
        "input_rows": len(rows),
        "modes": sorted({row["mode"] for row in rows}),
        "q_points": [list(q) for q in sorted({_q_key(row["q"]) for row in rows})],
        "meshes": [list(mesh) for mesh in sorted({_mesh_key(row["mesh"]) for row in rows if row.get("mesh")},
                                                   key=lambda item: (math.prod(item), item))],
        "plateau_relative_tolerance": relative_tolerance,
        "kgrid_absolute_tolerance_ry": kgrid_absolute_tolerance,
        "kgrid_relative_tolerance": kgrid_relative_tolerance,
        "plateaus": plateaus,
        "k_grid_convergence": kgrid,
        "pure_gauge_drift": pure_drift,
        "electron_consistency": electron_consistency,
        "q_half_mapping": mapping,
        "rows": rows,
    }


def analyze(paths: Sequence[str | Path], **kwargs: Any) -> dict[str, Any]:
    """Read and analyse multiple files, retaining every raw row."""
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.extend(read_sweep_file(path))
    return analyze_rows(rows, **kwargs)


def _wp10_unit(vector: Sequence[float]) -> tuple[float, float, float]:
    length = math.sqrt(sum(float(value) ** 2 for value in vector))
    if length <= Q_TOLERANCE:
        raise ValueError("q direction must be nonzero")
    return tuple(float(value) / length for value in vector)  # type: ignore[return-value]


def _wp10_direction(rows: Sequence[dict[str, Any]]) -> tuple[float, float, float]:
    first = next((row["q"] for row in rows if math.sqrt(sum(float(value) ** 2 for value in row["q"])) > Q_TOLERANCE), None)
    if first is None:
        raise ValueError("the sweep contains no nonzero q point")
    axis = max(range(3), key=lambda index: abs(float(first[index])))
    direction = [0.0, 0.0, 0.0]
    direction[axis] = 1.0 if float(first[axis]) >= 0.0 else -1.0
    return tuple(direction)  # type: ignore[return-value]


def _wp10_energy(row: dict[str, Any], definition: str) -> float | None:
    value = row.get("delta_gauge" if definition == "gauge" else "delta_raw")
    return None if value is None else float(value)


def _wp10_fit(points: Sequence[tuple[float, float]], quartic: bool) -> dict[str, Any]:
    """Fit y=A*q^2 (+ B*q^4), fixing the q=0 intercept at zero."""
    npar = 2 if quartic else 1
    design = [(q * q, q**4) if quartic else (q * q,) for q, _ in points]
    if len(points) < npar:
        return {"status": "insufficient_points", "n_points": len(points), "coefficients": {"A": None, "B": None}}
    if quartic:
        s40 = sum(row[0] * row[0] for row in design)
        s60 = sum(row[0] * row[1] for row in design)
        s80 = sum(row[1] * row[1] for row in design)
        sy2 = sum(row[0] * y for row, (_, y) in zip(design, points))
        sy4 = sum(row[1] * y for row, (_, y) in zip(design, points))
        determinant = s40 * s80 - s60 * s60
        if abs(determinant) <= 1.0e-300:
            return {"status": "singular", "n_points": len(points), "coefficients": {"A": None, "B": None}}
        coefficients = [(sy2 * s80 - sy4 * s60) / determinant, (s40 * sy4 - s60 * sy2) / determinant]
    else:
        denominator = sum(row[0] * row[0] for row in design)
        coefficients = [sum(row[0] * y for row, (_, y) in zip(design, points)) / denominator]
    residuals = [y - sum(c * x for c, x in zip(coefficients, row)) for row, (_, y) in zip(design, points)]
    rss = sum(value * value for value in residuals)
    dof = len(points) - npar
    sigma = math.sqrt(rss / dof) if dof > 0 else None
    if sigma is None:
        uncertainties = [None] * npar
    elif quartic:
        uncertainties = [sigma * math.sqrt(s80 / determinant), sigma * math.sqrt(s40 / determinant)]
    else:
        uncertainties = [sigma / math.sqrt(denominator)]
    return {
        "status": "fit", "n_points": len(points), "degrees_of_freedom": dof,
        "coefficients": {"A": coefficients[0], "B": coefficients[1] if quartic else None},
        "uncertainties": {"A": uncertainties[0], "B": uncertainties[1] if quartic else None},
        "rss": rss, "rmse": math.sqrt(rss / len(points)), "residuals": residuals,
        "q_points_internal": [q for q, _ in points],
    }


def _wp10_physical(fit: dict[str, Any], q_scale: float | None) -> dict[str, Any]:
    result = dict(fit)
    if q_scale is None or fit.get("coefficients", {}).get("A") is None:
        result["physical_units"] = None
        return result
    coefficients = fit["coefficients"]
    uncertainties = fit["uncertainties"]
    result["physical_units"] = {
        "q": "1/angstrom",
        "A_Ry_angstrom2": coefficients["A"] / q_scale**2,
        "B_Ry_angstrom4": None if coefficients["B"] is None else coefficients["B"] / q_scale**4,
        "A_uncertainty_Ry_angstrom2": None if uncertainties["A"] is None else uncertainties["A"] / q_scale**2,
        "B_uncertainty_Ry_angstrom4": None if uncertainties["B"] is None else uncertainties["B"] / q_scale**4,
    }
    return result


def _wp10_windows(points: Sequence[tuple[float, float]], q_scale: float | None) -> dict[str, Any]:
    magnitudes = sorted({round(abs(q), 10) for q, _ in points})
    windows: dict[str, list[dict[str, Any]]] = {"quadratic": [], "quadratic_quartic": []}
    for max_q in magnitudes:
        selected = [(q, value) for q, value in points if abs(q) <= max_q + Q_TOLERANCE]
        if len(selected) >= 2:
            fit = _wp10_physical(_wp10_fit(selected, quartic=False), q_scale)
            fit["max_q_internal"] = max_q
            windows["quadratic"].append(fit)
        if len(selected) >= 3:
            fit = _wp10_physical(_wp10_fit(selected, quartic=True), q_scale)
            fit["max_q_internal"] = max_q
            windows["quadratic_quartic"].append(fit)
    result: dict[str, Any] = {}
    for name, values in windows.items():
        result[name + "_windows"] = values
        result[name] = values[-1] if values else {"status": "insufficient_points", "n_points": 0}
        coefficients = [item["coefficients"]["A"] for item in values if item["coefficients"].get("A") is not None]
        if coefficients:
            reference = coefficients[-1]
            spread = (max(coefficients) - min(coefficients)) / max(abs(reference), 1.0e-300)
            result[name + "_A_stability"] = {
                "status": "stable" if spread <= DEFAULT_STABILITY_RELATIVE_TOLERANCE else "unstable",
                "relative_spread": spread, "tolerance": DEFAULT_STABILITY_RELATIVE_TOLERANCE,
                "reference_A": reference, "reference_max_q_internal": values[-1]["max_q_internal"],
            }
        else:
            result[name + "_A_stability"] = {"status": "unavailable", "relative_spread": None}
    return result


def _wp10_group(rows: Sequence[dict[str, Any]], direction: Sequence[float], definition: str,
                q_scale: float | None, odd_absolute_tolerance: float,
                odd_relative_tolerance: float) -> dict[str, Any]:
    samples: dict[tuple[float, float, float], list[float]] = defaultdict(list)
    for row in rows:
        value = _wp10_energy(row, definition)
        if value is not None:
            samples[tuple(round(float(value), 10) for value in row["q"])].append(value)
    projected: dict[float, dict[int, float]] = {}
    transverse = 0.0
    for q, values in samples.items():
        signed = sum(float(a) * float(b) for a, b in zip(q, direction))
        transverse = max(transverse, math.sqrt(sum((float(q[i]) - signed * direction[i])**2 for i in range(3))))
        magnitude = round(abs(signed), 10)
        sign = 0 if magnitude <= Q_TOLERANCE else (1 if signed > 0.0 else -1)
        projected.setdefault(magnitude, {})[sign] = sum(values) / len(values)
    zero = projected.get(0.0, {}).get(0)
    pairs: list[dict[str, Any]] = []
    points: list[tuple[float, float]] = []
    for magnitude in sorted(value for value in projected if value > Q_TOLERANCE):
        positive = projected[magnitude].get(1)
        negative = projected[magnitude].get(-1)
        if positive is None or negative is None:
            pairs.append({"q_abs_internal": magnitude, "status": "missing_positive" if positive is None else "missing_negative",
                          "q_abs_inv_angstrom": None if q_scale is None else magnitude * q_scale,
                          "energy_plus": positive, "energy_minus": negative})
            continue
        even = 0.5 * (positive + negative)
        odd = 0.5 * (positive - negative)
        delta_even = None if zero is None else even - zero
        item = {"q_abs_internal": magnitude,
                "q_abs_inv_angstrom": None if q_scale is None else magnitude * q_scale,
                "status": "paired" if delta_even is not None else "missing_gamma",
                "energy_plus": positive, "energy_minus": negative, "energy_even": even, "energy_odd": odd,
                "delta_even": delta_even, "delta_odd": odd}
        pairs.append(item)
        if delta_even is not None:
            points.append((magnitude, delta_even))
    complete = [item for item in pairs if item["status"] == "paired"]
    odd_values = [float(item["delta_odd"]) for item in complete]
    even_scale = max((abs(float(item["delta_even"])) for item in complete), default=0.0)
    odd_max = max((abs(value) for value in odd_values), default=None)
    odd_limit = odd_absolute_tolerance + odd_relative_tolerance * even_scale
    sin2 = rows[0].get("sin2theta")
    normalized = [] if sin2 is None or float(sin2) <= 0.0 else [(q, value / float(sin2)) for q, value in points]
    return {
        "energy_definition": definition, "zero_energy": zero, "pairs": pairs,
        "complete_pair_count": len(complete), "transverse_residual_internal": transverse,
        "odd_component": {"max_abs_Ry": odd_max,
                          "rms_Ry": math.sqrt(sum(value * value for value in odd_values) / len(odd_values)) if odd_values else None,
                          "tolerance_Ry": odd_limit,
                          "absolute_tolerance_Ry": odd_absolute_tolerance,
                          "relative_tolerance": odd_relative_tolerance,
                          "status": "within_tolerance" if odd_values and odd_max <= odd_limit else ("red_flag" if odd_values else "unavailable")},
        "fits": {"raw": _wp10_windows(points, q_scale), "sin2": _wp10_windows(normalized, q_scale)},
    }


def analyze_wp10_rows(rows: Sequence[dict[str, Any]], *, alat_angstrom: float | None = None,
                      direction: Sequence[float] | None = None, energy_definition: str = "auto",
                      odd_absolute_tolerance: float = 1.0e-10,
                      odd_relative_tolerance: float = 1.0e-6) -> dict[str, Any]:
    """Return WP10 q-reversal, fit-window, and mesh/cone sensitivity data."""
    if not rows:
        raise ValueError("no sweep rows supplied")
    if alat_angstrom is None:
        alat_values = {round(float(row["alat_angstrom"]), 10) for row in rows if row.get("alat_angstrom") is not None}
        if len(alat_values) == 1:
            alat_angstrom = alat_values.pop()
    if alat_angstrom is not None and alat_angstrom <= 0.0:
        raise ValueError("alat must be positive")
    direction = _wp10_unit(direction) if direction is not None else _wp10_direction(rows)
    q_scale = None if alat_angstrom is None else 2.0 * math.pi / float(alat_angstrom)
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        theta = row.get("theta_deg")
        grouped[(row.get("mode", "unknown"), _mesh_key(row.get("mesh")),
                 None if theta is None else round(float(theta), 8))].append(dict(row))
    groups: list[dict[str, Any]] = []
    for (mode, mesh, theta), group in sorted(grouped.items(), key=lambda item: (str(item[0][0]), item[0][1] or (), item[0][2] if item[0][2] is not None else -1.0)):
        definition = energy_definition
        if definition == "auto":
            definition = "gauge" if any(_wp10_energy(row, "gauge") is not None for row in group) else "raw"
        groups.append({"mode": mode, "mesh": list(mesh) if mesh else None, "theta_deg": theta,
                       "direction_internal": list(direction),
                       "reversal": _wp10_group(group, direction, definition, q_scale,
                                                odd_absolute_tolerance, odd_relative_tolerance)})
    sensitivity: dict[str, list[dict[str, Any]]] = {"raw": [], "sin2": []}
    for quantity in sensitivity:
        by_mode: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for group in groups:
            fit = group["reversal"]["fits"][quantity]["quadratic_quartic"]
            coefficient = fit.get("coefficients", {}).get("A")
            if coefficient is not None:
                by_mode[group["mode"]].append({"mesh": group["mesh"], "theta_deg": group["theta_deg"],
                                                "A": coefficient, "A_uncertainty": fit.get("uncertainties", {}).get("A")})
        for mode, values in sorted(by_mode.items()):
            coefficients = [item["A"] for item in values]
            sensitivity[quantity].append({"mode": mode, "samples": values,
                                          "A_min": min(coefficients), "A_max": max(coefficients),
                                          "A_spread": max(coefficients) - min(coefficients),
                                          "A_relative_spread": (max(coefficients) - min(coefficients)) / max(abs(sum(coefficients) / len(coefficients)), 1.0e-300)})
    return {"schema": "gbt_wp10_small_q_v1", "input_rows": len(rows),
            "raw_rows": [dict(row) for row in rows], "q_units_internal": "Cartesian 2*pi/alat",
            "direction_internal": list(direction), "alat_angstrom": alat_angstrom,
            "q_scale_inv_angstrom": q_scale, "physical_q_units": "1/angstrom" if q_scale is not None else None,
            "groups": groups, "sensitivity": sensitivity,
            "odd_component_status": [{"mode": group["mode"], "mesh": group["mesh"], "theta_deg": group["theta_deg"],
                                       **group["reversal"]["odd_component"]} for group in groups]}


def analyze_wp10(paths: Sequence[str | Path], **kwargs: Any) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.extend(read_sweep_file(path))
    return analyze_wp10_rows(rows, **kwargs)


def format_wp10_report(result: dict[str, Any]) -> str:
    lines = ["GBT WP10 small-q quadratic and q-reversal analysis",
             f"raw rows: {result['input_rows']}; q units: {result['q_units_internal']}",
             f"direction: {tuple(result['direction_internal'])}; q scale: {result['q_scale_inv_angstrom'] or 'n/a'} 1/angstrom", "",
             "Q-REVERSAL, QUADRATIC+QUARTIC FITS, AND WINDOW STABILITY"]
    for group in result["groups"]:
        reversal = group["reversal"]; odd = reversal["odd_component"]
        lines.append(f"  mode={group['mode']} mesh={group['mesh']} theta={group['theta_deg']} definition={reversal['energy_definition']} "
                     f"pairs={reversal['complete_pair_count']} odd={odd['status']} max|odd|={odd['max_abs_Ry']}")
        for quantity in ("raw", "sin2"):
            fit = reversal["fits"][quantity]; q2 = fit["quadratic"]; q4 = fit["quadratic_quartic"]
            q4_physical = q4.get("physical_units") or {}
            lines.append(f"    {quantity}: A(q2)={q2.get('coefficients', {}).get('A')} +/- {q2.get('uncertainties', {}).get('A')}; "
                         f"A(q2+q4)={q4.get('coefficients', {}).get('A')} +/- {q4.get('uncertainties', {}).get('A')}, "
                         f"B={q4.get('coefficients', {}).get('B')}; stability={fit['quadratic_quartic_A_stability']['status']}; "
                         f"physical A={q4_physical.get('A_Ry_angstrom2', 'n/a')} Ry Angstrom^2, "
                         f"B={q4_physical.get('B_Ry_angstrom4', 'n/a')} Ry Angstrom^4")
    lines.append("\nMESH/CONE SENSITIVITY (q2+q4 A)")
    for quantity, entries in result["sensitivity"].items():
        for item in entries:
            lines.append(f"  {quantity} mode={item['mode']}: samples={len(item['samples'])}, A range=[{item['A_min']}, {item['A_max']}], "
                         f"relative spread={item['A_relative_spread']}")
    return "\n".join(lines)


def _fmt(value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.6e}"
    return str(value)


def format_report(result: dict[str, Any]) -> str:
    """Format the analysis as a compact human-readable report."""
    lines = [
        "GBT WP09 harmonic cone-angle and k-grid analysis",
        f"raw rows: {result['input_rows']}; modes: {', '.join(result['modes']) or 'none'}",
        "",
        "HARMONIC PLATEAUS (K = DeltaE/sin(theta)^2)",
    ]
    if not result["plateaus"]:
        lines.append("  no angle groups contain enough usable points")
    for item in result["plateaus"]:
        admitted = item["admitted_angles_deg"]
        excluded = item["excluded_angles_deg"]
        lines.append(
            f"  mode={item['mode']} q={tuple(item['q'])} mesh={item['mesh']} {item['definition']}: "
            f"{item['status']}, admitted={admitted or 'none'}, excluded={excluded or 'none'}, "
            f"spread={_fmt(item['spread_rel'])} relative ({_fmt(item['spread_abs'])} absolute)"
        )
    lines.extend(["", "K-GRID CONVERGENCE (raw delta before any fitting)"])
    for item in result["k_grid_convergence"]:
        lines.append(
            f"  mode={item['mode']} q={tuple(item['q'])} theta={item['theta_deg']} "
            f"{item['definition']}: meshes={len(item['meshes'])}, spread={_fmt(item['spread_abs'])} Ry "
            f"({_fmt(item['spread_rel'])} relative), status={item['status']}"
        )
    lines.extend(["", "PURE-GAUGE DRIFT E(q,0)-E(Gamma,0)"])
    for item in result["pure_gauge_drift"]:
        lines.append(
            f"  mode={item['mode']} q={tuple(item['q'])} theta={item['theta_deg']} mesh={item['mesh']}: "
            f"{_fmt(item['drift'])} Ry ({item['reference']})"
        )
    lines.extend(["", "FERMI/ELECTRON-COUNT CONSISTENCY"])
    for item in result["electron_consistency"]:
        lines.append(
            f"  mode={item['mode']} q={tuple(item['q'])} theta={item['theta_deg']}: "
            f"max|dN|={_fmt(item['max_abs_electron_error'])}, "
            f"EF=[{_fmt(item['fermi_min'])}, {_fmt(item['fermi_max'])}] Ry"
        )
    controls = result["q_half_mapping"]
    known = sum(item["q_half_maps_to_mesh"] is not None for item in controls)
    lines.extend(["", f"q/2 MESH-MAPPING CONTROL: {known}/{len(controls)} rows classified"])
    if known:
        lines.append("  classification is file metadata or the uniform-Cartesian fallback; use a lattice-aware value for non-cubic cells")
    else:
        lines.append("  no mapping metadata available")
    return "\n".join(lines)


def _expand_inputs(inputs: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for value in inputs:
        path = Path(value)
        if path.is_dir():
            candidates = sorted(path.rglob("frozen_magnon_harmonic_diagnostics.dat"))
            candidates += sorted(path.rglob("*.csv"))
            if not candidates:
                raise ValueError(f"{path}: directory contains no WP09 sweep files")
            paths.extend(candidates)
        else:
            paths.append(path)
    return paths


def main_wp10(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Analyse symmetric small-q GBT frozen-magnon sweeps")
    parser.add_argument("inputs", nargs="+", help="WP09 sweep files or directories")
    parser.add_argument("--alat", type=float, help="lattice parameter in Angstrom")
    parser.add_argument("--direction", type=float, nargs=3, metavar=("QX", "QY", "QZ"))
    parser.add_argument("--energy-definition", choices=("auto", "raw", "gauge"), default="auto")
    parser.add_argument("--odd-absolute-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--odd-relative-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        result = analyze_wp10(_expand_inputs(args.inputs), alat_angstrom=args.alat, direction=args.direction,
                              energy_definition=args.energy_definition, odd_absolute_tolerance=args.odd_absolute_tolerance,
                              odd_relative_tolerance=args.odd_relative_tolerance)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(payload + "\n")
    print(payload if args.json else format_wp10_report(result))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if "--wp10" in argv:
        return main_wp10([arg for arg in argv if arg != "--wp10"])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="WP09 sweep files or directories")
    parser.add_argument("--json", action="store_true", help="print the complete JSON result")
    parser.add_argument("--json-out", type=Path, help="write the complete JSON result to this file")
    parser.add_argument("--plateau-relative-tolerance", type=float, default=0.05)
    parser.add_argument("--min-plateau-points", type=int, default=2)
    parser.add_argument("--kgrid-absolute-tolerance", type=float, default=KGRID_ABSOLUTE_TOLERANCE,
                        help="absolute mesh-spread tolerance in Ry")
    parser.add_argument("--kgrid-relative-tolerance", type=float, default=KGRID_RELATIVE_TOLERANCE,
                        help="relative mesh-spread tolerance")
    args = parser.parse_args(argv)
    try:
        paths = _expand_inputs(args.inputs)
        rows: list[dict[str, Any]] = []
        for path in paths:
            rows.extend(read_sweep_file(path))
        result = analyze_rows(rows, relative_tolerance=args.plateau_relative_tolerance,
                              min_plateau_points=args.min_plateau_points,
                              kgrid_absolute_tolerance=args.kgrid_absolute_tolerance,
                              kgrid_relative_tolerance=args.kgrid_relative_tolerance)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(payload + "\n")
    if args.json:
        print(payload)
    else:
        print(format_report(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())
