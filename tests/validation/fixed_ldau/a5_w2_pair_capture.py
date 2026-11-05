#!/usr/bin/env python3
"""Verify and combine paired full/moments-only W2 endpoint captures.

Only direct P_OUT density/q0/Tr(WP) rows are copied from moments-only into
the full trace. Both immutable raw capture directories remain untouched.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_capture(path: Path) -> tuple[dict, list[dict]]:
    manifest = json.loads((path / "manifest.json").read_text())
    records = path / "records.csv"
    if sha(records) != manifest["records_sha256"]:
        raise ValueError(f"records hash mismatch: {records}")
    with records.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if {int(row["iteration"]) for row in rows} != {1}:
        raise ValueError(f"capture is not a single selected iteration: {path}")
    return manifest, rows


def thermo(path: Path) -> tuple[float, float]:
    log = (path / "stdout.log").read_text(errors="replace")
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    matches = re.findall(rf"entropy=\s*({number}).*?TS=\s*({number})", log)
    if not matches:
        raise ValueError(f"finite entropy/TS absent from {path}/stdout.log")
    return tuple(float(value.replace("D", "E").replace("d", "e")) for value in matches[-1])


def keyed(rows: list[dict], kind: str, stage: str | None = None) -> dict[tuple, dict]:
    selected = [r for r in rows if r["kind"] == kind and (stage is None or r["stage"] == stage)]
    return {(r["site"], r["l"], r["spin"], r["order"]): r for r in selected}


def close(a: float, b: float, label: str, atol: float = 2e-12) -> None:
    if not (math.isfinite(a) and math.isfinite(b)) or abs(a - b) > atol + 2e-13 * max(abs(a), abs(b)):
        raise ValueError(f"paired capture mismatch {label}: {a:.17g} vs {b:.17g}")


def verify(full_dir: Path, moments_dir: Path, merged: Path) -> dict:
    fm, full = read_capture(full_dir)
    mm, moments = read_capture(moments_dir)
    for key in ("binary_sha256", "lambda_ry", "converged_endpoint_restart_sha256",
                "capture_restart_pre_sha256", "input_pre_sha256"):
        if fm[key] != mm[key]:
            raise ValueError(f"paired captures differ in {key}: {fm[key]} vs {mm[key]}")
    if fm["moments_only"] or not mm["moments_only"]:
        raise ValueError("capture modes are not full/moments-only")
    fk, mk = {r["kind"] for r in full}, {r["kind"] for r in moments}
    radial = {"radial_replay_mesh", "radial_replay_v", "radial_replay_rho_in",
              "radial_replay_ec_seed", "radial_replay_eigen_seed", "radial_replay_wave"}
    direct = {"test_density_matrix", "test_q0_trace", "test_trace_wp"}
    if not radial <= fk or radial & mk or not direct <= mk or direct & fk:
        raise ValueError("paired mode output contract failed")

    f_energy = keyed(full, "etot_total"); m_energy = keyed(moments, "etot_total")
    if set(f_energy) != set(m_energy):
        raise ValueError("ETOT row keys differ")
    for key in f_energy:
        close(float(f_energy[key]["scalar"]), float(m_energy[key]["scalar"]), f"ETOT {key}")
    f_entropy, f_ts = thermo(full_dir); m_entropy, m_ts = thermo(moments_dir)
    close(f_entropy, m_entropy, "entropy"); close(f_ts, m_ts, "TS")

    f_ql = keyed(full, "ql", "P_OUT")
    m_ql = keyed(moments, "ql", "P_OUT")
    f_moment = keyed(full, "radial_moment", "P_OUT")
    m_moment = keyed(moments, "radial_moment", "P_OUT")
    f_param = keyed(full, "radial_parameter", "P_OUT")
    m_param = keyed(moments, "radial_parameter", "P_OUT")
    for label, left, right in (("ql", f_ql, m_ql), ("radial moment", f_moment, m_moment),
                               ("radial parameter", f_param, m_param)):
        if set(left) != set(right):
            raise ValueError(f"P_OUT {label} keys differ")
        for key in left:
            close(float(left[key]["scalar"]), float(right[key]["scalar"]), f"{label} {key}")

    ql_d = {(r["spin"], int(r["l"])): float(r["scalar"])
            for r in full if r["kind"] == "ql" and r["stage"] == "P_OUT" and r["l"] == "2"}
    p_q0 = {r["spin"]: float(r["scalar"]) for r in moments if r["kind"] == "test_q0_trace"}
    if set(p_q0) != {"0", "1"}:
        raise ValueError(f"direct physical-spin q0 coverage incorrect: {set(p_q0)}")
    for spin, q0 in p_q0.items():
        close(q0, ql_d[(str(int(spin) + 1), 2)], f"P_OUT direct q0 vs ql spin={spin}")
    p_rows = [r for r in moments if r["kind"] in direct and r["stage"] == "P_OUT"]
    matrix = [r for r in p_rows if r["kind"] == "test_density_matrix"]
    for spin in ("0", "1"):
        if len([r for r in matrix if r["spin"] == spin]) != 25:
            raise ValueError(f"direct P block incomplete, spin={spin}")

    # Merge only the direct P_OUT rows; all scientific/radial operands come
    # from the full trace. Keep moments-only raw file and manifest unchanged.
    merged.parent.mkdir(parents=True, exist_ok=True)
    with merged.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(full[0].keys()), lineterminator="\n")
        writer.writeheader(); writer.writerows(full); writer.writerows(p_rows)
    receipt = {
        "status": "A5.4_W2_PAIRED_CAPTURE_VERIFIED",
        "lambda_ry": fm["lambda_ry"],
        "full_mode_manifest_sha256": sha(full_dir / "manifest.json"),
        "moments_mode_manifest_sha256": sha(moments_dir / "manifest.json"),
        "full_records_sha256": sha(full_dir / "records.csv"),
        "moments_records_sha256": sha(moments_dir / "records.csv"),
        "capture_runner_sha256": sha(Path(__file__).resolve().parent / "a5_w2_endpoint_capture.py"),
        "pair_verifier_sha256": sha(Path(__file__).resolve()),
        "merged_records_sha256": sha(merged),
        "merged_records_path": str(merged),
        "copied_direct_rows": len(p_rows),
        "direct_matrix_rows": len(matrix),
        "verified_same_endpoint_input_binary": True,
        "capture_mode_environment_delta": {
            "full_RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "0",
            "moments_RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1",
            "all_other_capture_environment_values": "identical",
        },
        "same_etot_entropy_ts_and_pout_ql_moments_parameters": True,
        "full_has_radial_seed_and_wave_operands": True,
        "moments_only_has_no_radial_replay_rows": True,
        "capture_completeness_check_only": True,
        "entropy": f_entropy,
        "TS_ry": f_ts,
        "etot_ry": float(next(iter(f_energy.values()))["scalar"]),
        "q0_by_physical_spin": p_q0,
    }
    (merged.parent / "pair_manifest.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", type=Path, required=True)
    ap.add_argument("--moments", type=Path, required=True)
    ap.add_argument("--merged", type=Path, required=True)
    args = ap.parse_args()
    print(json.dumps(verify(args.full, args.moments, args.merged), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
