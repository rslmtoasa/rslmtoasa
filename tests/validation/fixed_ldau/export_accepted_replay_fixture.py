#!/usr/bin/env python3
"""Export an exact, hash-bound production radial replay fixture.

The fixture is derived only from the immutable accepted A1 no-U iteration-23
capture. CSV numeric lexemes are retained as JSON strings so a loader can
round-trip the source doubles without decimal reformatting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


EXPECTED_CSV = "55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4"
STAGES = ("RAW_PL", *(f"A1_SHIFT_{i:02d}" for i in range(1, 9)))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exactly(rows: list[dict[str, str]], kind: str, stage: str, *, order: int | None = None,
            row: int | None = None, col: int | None = None, ell: int | None = None,
            spin: int | None = None, band: int | None = None) -> dict[str, str]:
    found = [r for r in rows if r["kind"] == kind and r["stage"] == stage
             and (order is None or int(r["order"]) == order)
             and (row is None or int(r["row"]) == row)
             and (col is None or int(r["col"]) == col)
             and (ell is None or int(r["l"]) == ell)
             and (spin is None or int(r["spin"]) == spin)
             and (band is None or int(r["band"]) == band)]
    if len(found) != 1:
        raise ValueError(f"expected one {kind}/{stage} row at order,row,col,l,spin="
                         f"{order},{row},{col},{ell},{spin}; found {len(found)}")
    return found[0]


def array(rows: list[dict[str, str]], kind: str, stage: str, *, field: str,
          key_fields: tuple[str, ...], expected_keys: list[tuple[int, ...]]) -> list[str]:
    found: dict[tuple[int, ...], str] = {}
    for r in rows:
        if r["kind"] != kind or r["stage"] != stage:
            continue
        key = tuple(int(r[k]) for k in key_fields)
        if key in found:
            raise ValueError(f"duplicate {kind}/{stage} array key {key}")
        found[key] = r[field]
    if set(found) != set(expected_keys):
        missing = sorted(set(expected_keys) - set(found))[:3]
        extra = sorted(set(found) - set(expected_keys))[:3]
        raise ValueError(f"incomplete {kind}/{stage}: missing={missing}, extra={extra}")
    return [found[k] for k in expected_keys]


def parse_fe(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    element = text.split("&element", 1)[1].split("/", 1)[0]
    def match(pattern: str) -> str:
        m = re.search(pattern, element, flags=re.I)
        if not m:
            raise ValueError(f"missing Fe element field matching {pattern!r}")
        return m.group(1)
    return {"symbol": match(r"\bsymbol\s*=\s*'([^']+)'"),
            "f_core": int(match(r"\bf_core\s*=\s*([+-]?\d+)")),
            "atomic_number": match(r"\batomic_number\s*=\s*([+\-0-9.Ee]+)"),
            "core_electrons": match(r"\bcore\s*=\s*([+\-0-9.Ee]+)"),
            "lmax": int(re.search(r"\blmax\s*=\s*(\d+)", text, flags=re.I).group(1))}


def export(capture: Path) -> dict[str, Any]:
    csv_path = capture / "records.csv"
    manifest_path = capture / "manifest.json"
    fe_path = capture / "Fe.nml"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    csv_hash = sha(csv_path)
    if csv_hash != EXPECTED_CSV or manifest.get("post_run_sha256", {}).get("records.csv") != EXPECTED_CSV:
        raise ValueError("accepted capture CSV hash is not the immutable approved A1 CSV")
    if manifest.get("selected_trace_iteration") != "rolling final (TRACE_ITERATION=0)":
        raise ValueError("accepted capture does not identify its rolling-final iteration")
    if manifest.get("backend") != "KS" or manifest.get("expected_accepted_iteration") != 23:
        raise ValueError("fixture source is not the accepted scalar KS iteration-23 capture")
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    iters = {int(r["iteration"]) for r in rows}
    if iters != {23}:
        raise ValueError(f"expected only accepted iteration 23, got {sorted(iters)}")
    nr, lmax, nsp = 495, 2, 2
    controls = {i: exactly(rows, "m0_replay_control", "RAW_PL", order=i)
                for i in range(6)}
    params = {"a": controls[0]["energy"], "z": controls[0]["scalar"],
              "b": controls[1]["energy"], "tol": controls[1]["scalar"],
              "nr": int(float(controls[2]["energy"])),
              "ncore": int(float(controls[2]["scalar"])),
              "nval": int(float(controls[2]["re"])),
              "control_rows_exact": [{k: r[k] for k in ("order", "energy", "scalar", "re", "im")}
                                     for r in controls.values()]}
    if (params["nr"], params["ncore"], params["nval"]) != (nr, 10, 6):
        raise ValueError("captured production radial extents differ from the verified schema")

    radial = {
        "rofi": array(rows, "m0_replay_rofi", "RAW_PL", field="scalar",
                      key_fields=("row",), expected_keys=[(i,) for i in range(nr)]),
        "v_by_row_spin": array(rows, "m0_replay_v_input", "RAW_PL", field="scalar",
                      key_fields=("row", "spin"), expected_keys=[(i, s) for i in range(nr) for s in (1, 2)]),
        "rho_input_by_row_spin": array(rows, "m0_replay_rho_input", "RAW_PL", field="scalar",
                      key_fields=("row", "spin"), expected_keys=[(i, s) for i in range(nr) for s in (1, 2)]),
        "ec_seed": array(rows, "m0_replay_ec_seed", "RAW_PL", field="scalar",
                      key_fields=("row",), expected_keys=[(i,) for i in range(params["ncore"])]),
        "ev_seed": array(rows, "m0_replay_ev_seed", "RAW_PL", field="scalar",
                      key_fields=("row",), expected_keys=[(i,) for i in range(params["nval"])]),
    }
    pl_inputs: dict[str, Any] = {}
    for s in (1, 2):
        pl_inputs[f"spin_{s}"] = {}
        for ell in range(lmax + 1):
            pl_inputs[f"spin_{s}"][f"l_{ell}"] = {
                str(order): exactly(rows, "m0_replay_pl_input", "RAW_PL", ell=ell,
                                    spin=s, order=order)["scalar"] for order in range(5)}
    ql = {f"spin_{s}": {f"l_{ell}": [
                exactly(rows, "m0_replay_ql_input", "RAW_PL", ell=ell, spin=s, order=m)["scalar"]
                for m in range(3)] for ell in range(lmax + 1)} for s in (1, 2)}
    raw_pl = {f"spin_{s}": [exactly(rows, "m0_replay_pl_return", "RAW_PL", ell=ell,
                                  spin=s, order=0)["scalar"] for ell in range(lmax + 1)]
              for s in (1, 2)}
    historical_a1_pl: dict[str, Any] = {}
    for stage in STAGES[1:]:
        historical_a1_pl[stage] = {f"spin_{s}": [
            exactly(rows, "m0_replay_pl_return", stage, ell=ell, spin=s, order=0)["scalar"]
            for ell in range(lmax + 1)] for s in (1, 2)}
    raw_energy = {f"spin_{s}": {
        "sumec": exactly(rows, "m0_replay_energy", "RAW_PL", spin=s, order=0)["energy"],
        "sumev": exactly(rows, "m0_replay_energy", "RAW_PL", spin=s, order=0)["scalar"]}
        for s in (1, 2)}
    # Keep the independently emitted A1 case summary beside the inputs; the
    # standalone loader must reproduce its raw case before any shift is read.
    for s in (1, 2):
        ref = exactly(rows, "m0_a11_case", "RAW_PL", spin=s, ell=2)
        if float(raw_energy[f"spin_{s}"]["sumev"]) != float(ref["energy"]):
            raise ValueError(f"RAW_PL SUMEV disagrees between replay and A1 baseline row for spin {s}")
    component_controls = {
        "vhrmax": controls[4]["energy"], "b_fsm": controls[5]["energy"],
    }
    raw_ec_ev_seeds = {
        "ec_seed": radial["ec_seed"], "ev_seed": radial["ev_seed"],
        "ec_output_by_spin": {f"spin_{s}": exactly(rows, "m0_replay_energy", "RAW_PL",
                                      spin=s, order=0)["energy"] for s in (1, 2)},
        "ev_output_by_spin": {f"spin_{s}": exactly(rows, "m0_replay_energy", "RAW_PL",
                                      spin=s, order=0)["scalar"] for s in (1, 2)},
    }
    raw_rho = array(rows, "m0_replay_rho", "RAW_PL", field="scalar",
                    key_fields=("row", "spin"), expected_keys=[(i, s) for i in range(nr) for s in (1, 2)])
    radial_verification = {}
    for index, name in enumerate(("rho_max_abs", "sumec_max_abs_ry", "sumev_max_abs_ry",
                                  "ec_max_abs_ry", "ev_max_abs_ry")):
        radial_verification[name] = exactly(rows, f"radial_replay_{name}", "RADIAL_REPLAY",
                                            order=index)["scalar"]
    old_a1: dict[str, Any] = {}
    for stage in STAGES:
        old_a1[stage] = []
        for spin in (1, 2):
            r = exactly(rows, "m0_a11_case", stage, spin=spin, ell=2)
            old_a1[stage].append({k: r[k] for k in ("weight", "energy", "scalar", "re", "im")})
    raw_valence_ev = []
    for state in range(params["nval"]):
        r = exactly(rows, "m0_a11_ev", "RAW_PL", band=state)
        raw_valence_ev.append({k: r[k] for k in ("l", "spin", "band", "energy", "scalar", "re")})

    return {
        "schema": "rslmto-accepted-iter23-private-radial-replay-fixture-v1",
        "provenance": {
            "capture_dir": str(capture), "capture_csv_sha256": csv_hash,
            "capture_manifest_sha256": sha(manifest_path), "capture_self_sha256": manifest["self_f90_sha256"],
            "release_binary_sha256": manifest["binary_sha256"],
            "input_pre_sha256": manifest["input_sha256"]["input.nml.pre"],
            "fe_pre_sha256": manifest["input_sha256"]["Fe.nml.pre"],
            "fe_nml_sha256": sha(fe_path), "input_nml_postrun_sha256": sha(capture / "input.nml"),
            "accepted_iteration": 23, "accepted_fe_out_sha256": manifest["post_run_sha256"]["Fe_out.nml"],
            "numeric_encoding": "all source CSV numeric values preserved as strings",
        },
        "static_atom": parse_fe(fe_path),
        "dimensions": {"lmax": lmax, "nsp": nsp, "nr": nr,
                       "ql_axes": ["spin 1..2", "l 0..2", "moment q0/q1/q2"],
                       "radial_axes": ["row 0..494", "spin 1..2"],
                       "pl_axes": ["spin 1..2", "l 0..2", "operand order center,c,qpar,srdel,old_pl"]},
        "parameters": params,
        "asa_component_parameters": component_controls,
        "hyperfine": {"value": False, "basis": "source default source/control.f90:588; not set in input.nml"},
        "radial": radial, "pl_input_by_spin_l": pl_inputs, "ql_by_spin_l": ql,
        "raw_pl_return_by_spin": raw_pl, "historical_a1_pl_return_by_stage_spin": historical_a1_pl,
        "raw_energy_by_spin": raw_energy,
        "raw_rho_output_by_row_spin": raw_rho, "raw_neWRHO_reference": raw_ec_ev_seeds,
        "raw_live_neWRHO_verification_error": radial_verification,
        "raw_valence_ev_reference_by_state": raw_valence_ev,
        "historical_a1_case_rows_by_stage_spin": old_a1,
        "stage_names": list(STAGES),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("capture", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    result = export(args.capture)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(json.dumps({"output": str(args.output),
                      "fixture_sha256": hashlib.sha256(payload.encode()).hexdigest(),
                      "capture_csv_sha256": result["provenance"]["capture_csv_sha256"],
                      "nr": result["dimensions"]["nr"], "stages": result["stage_names"]}, sort_keys=True))


if __name__ == "__main__":
    main()
