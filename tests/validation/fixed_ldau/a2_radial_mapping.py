#!/usr/bin/env python3
"""Non-gating A2.4 D(E)->EVAL and PHDFSR derivative diagnostic.

The capture CSV and JSON manifest are inputs only.  Each radial solve is made
by a fresh process calling self_mod.legacy_radial_fixture (RSEQSR/PHDFSR).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path

H_RY = (0.006, 0.003, 0.0015, 0.00075)
C_AU = 274.074


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def value(row: dict[str, str], key: str = "scalar") -> float:
    return float(row[key])


def exactly(rows: list[dict[str, str]], pred, label: str) -> dict[str, str]:
    found = [r for r in rows if pred(r)]
    if len(found) != 1:
        raise ValueError(f"expected one {label} row, found {len(found)}")
    return found[0]


def solve(driver: Path, rofi: list[float], pot: list[float], z: float, a: float,
          b: float, boundary: float) -> dict:
    with tempfile.TemporaryDirectory(prefix="a2-radial-") as td:
        inp = Path(td) / "radial.in"
        with inp.open("w") as f:
            f.write(f"{len(rofi)} 2 {z:.17g} {a:.17g} {b:.17g} {boundary:.17g}\n")
            for r, v in zip(rofi, pot):
                f.write(f"{r:.17g} {v:.17g}\n")
        out = subprocess.run([str(driver), str(inp)], check=True, text=True,
                             capture_output=True).stdout.splitlines()
    energy = float(out[0].split()[1])
    arr = [[float(x) for x in line.split()] for line in out[1:]]
    if len(arr) != len(rofi) or any(len(row) != 4 for row in arr):
        raise ValueError("radial driver returned an incomplete G/GP mesh")
    return {"eval": energy, "g": [r[:2] for r in arr], "gp": [r[2:] for r in arr]}


def d_of_e(e: float, c: float, qpar: float, srdel: float) -> float:
    # calculate_pl's unrounded dnu; the integer PL branch is deliberately out.
    ell = 2
    rq = 1.0 / qpar
    delta2 = srdel * srdel
    dc = c - e
    return ell + (2 * ell + 1) / (rq * dc / (2 * (2 * ell + 1) * (dc - delta2 * rq)) - 1)


def gintsr_inner(x, y, rofi, pot, a, b, z, e):
    # Exact legacy GINTSR Simpson ordering and scalar-relativistic metric.
    nr = len(rofi)
    def term(i):
        r = rofi[i]
        tmc = C_AU - (pot[i] - 2 * z / r - e) / C_AU
        gfac = 1 + 6 / (tmc * r) ** 2
        return (r + b) * (x[i][0] * y[i][0] * gfac + x[i][1] * y[i][1])
    s = 4 * sum(term(i) for i in range(1, nr - 1, 2))
    s += 2 * sum(term(i) for i in range(2, nr - 2, 2))
    s += term(nr - 1)
    return a * s / 3


def nodes(g):
    vals = [r[0] for r in g if abs(r[0]) > 1.0e-14]
    return sum(1 for x, y in zip(vals, vals[1:]) if x*y < 0)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--capture", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--iteration", type=int, required=True)
    p.add_argument("--site", type=int, required=True)
    p.add_argument("--spin", type=int, required=True, help="CSV zero-based spin index")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = list(csv.DictReader(args.capture.open(newline="")))
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("backend", "").upper() != "KS":
        raise ValueError("A2.4 helper expects the accepted KS capture")
    bound_csv = manifest.get("csv_sha256")
    if bound_csv and bound_csv != sha(args.capture):
        raise ValueError("capture CSV hash does not match the supplied manifest")
    it, site, spin = args.iteration, args.site, args.spin
    common = lambda r: int(r["iteration"]) == it and int(r["site"]) == site
    mesh = sorted([r for r in rows if common(r) and r["kind"] == "radial_replay_mesh"
                   and int(r["spin"]) == spin], key=lambda r: int(r["row"]))
    vrows = sorted([r for r in rows if common(r) and r["kind"] == "radial_replay_v"
                    and int(r["spin"]) == spin], key=lambda r: int(r["row"]))
    if len(mesh) < 35 or len(mesh) % 2 != 1 or len(vrows) != len(mesh):
        raise ValueError("missing, mismatched, or unsupported captured radial mesh/potential")
    rofi = [value(r, "energy") for r in mesh]
    pot = [value(r) for r in vrows]
    metadata = {}
    for key in ("radial_A", "radial_B", "radial_Z", "radial_NR"):
        metadata[key] = value(exactly(rows, lambda r: common(r) and r["kind"] == key, key))
    if int(metadata["radial_NR"]) != len(rofi):
        raise ValueError("captured radial_NR does not match mesh length")
    enu = value(exactly(rows, lambda r: common(r) and r["kind"] == "radial_replay_enu"
                          and int(r["l"]) == 2 and int(r["spin"]) == spin, "d-shell Enu"))
    plrows = [r for r in rows if common(r) and r["kind"] == "pl_map_input" and
              r["stage"] == "CALCULATE_PL" and int(r["l"]) == 2 and int(r["spin"]) == spin + 1]
    pl = {int(r["order"]): value(r) for r in plrows}
    if not all(i in pl for i in range(6)):
        raise ValueError("incomplete six-field calculate_pl input tuple")
    center, c, qpar, srdel, pl_saved = (pl[0], pl[1], pl[2], pl[3], pl[4])
    d_base = d_of_e(enu, c, qpar, srdel)
    # Check the saved phase branch independently; only the unrounded D enters RSEQSR.
    d_saved = d_of_e(center, c, qpar, srdel)
    pl_calc = 0.5 - math.atan(d_saved) / math.pi + math.trunc(pl_saved)
    result = {"diagnostic": "A2.4_non_gating", "iteration": it, "site": site, "spin": spin,
              "source_sha256": sha(Path(__file__)), "driver_sha256": sha(args.driver),
              "capture_sha256": sha(args.capture), "manifest_sha256": sha(args.manifest),
              "driver_path": str(args.driver.resolve()), "capture_path": str(args.capture.resolve()),
              "manifest_path": str(args.manifest.resolve()), "mesh_points": len(rofi),
              "mesh_sha256": hashlib.sha256("\n".join(f"{x:.17g}" for x in rofi).encode()).hexdigest(),
              "captured": metadata, "enu_ry": enu, "potpar": {"center": center, "c": c,
              "qpar": qpar, "srdel": srdel, "pl_saved": pl_saved},
              "D_enu": d_base, "D_plus_one_boundary": d_base + 1,
              "D_at_saved_center": d_saved,
              "PL_branch_reconstruction": pl_calc,
              "PL_branch_difference": pl_calc - pl_saved, "phdfsr_dele_ry": 0.003,
              "h_sweep_ry": list(H_RY), "points": []}
    center_solve = solve(args.driver, rofi, pot, metadata["radial_Z"], metadata["radial_A"],
                         metadata["radial_B"], d_base + 1)
    result["center_solve"] = {"input_E": enu, "EVAL": center_solve["eval"],
                              "f_E": center_solve["eval"] - enu,
                              "large_component_nodes": nodes(center_solve["g"])}
    for h in H_RY:
        plus = solve(args.driver, rofi, pot, metadata["radial_Z"], metadata["radial_A"],
                     metadata["radial_B"], d_of_e(enu + h, c, qpar, srdel) + 1)
        minus = solve(args.driver, rofi, pot, metadata["radial_Z"], metadata["radial_A"],
                      metadata["radial_B"], d_of_e(enu - h, c, qpar, srdel) + 1)
        fprime = (plus["eval"] - minus["eval"]) / (2 * h) - 1
        # Differentiate the normalized radial solution against actual EVAL spacing.
        overlap_plus = gintsr_inner(center_solve["g"], plus["g"], rofi, pot,
                                    metadata["radial_A"], metadata["radial_B"], metadata["radial_Z"], enu)
        overlap_minus = gintsr_inner(center_solve["g"], minus["g"], rofi, pot,
                                     metadata["radial_A"], metadata["radial_B"], metadata["radial_Z"], enu)
        phase_plus = 1.0 if overlap_plus >= 0 else -1.0
        phase_minus = 1.0 if overlap_minus >= 0 else -1.0
        gp_fd = [[(phase_plus*plus["g"][i][j] - phase_minus*minus["g"][i][j]) /
                   (plus["eval"] - minus["eval"]) for j in range(2)] for i in range(len(rofi))]
        diff = [[center_solve["gp"][i][j] - gp_fd[i][j]
                 for j in range(2)] for i in range(len(rofi))]
        gp_norm = gintsr_inner(center_solve["gp"], center_solve["gp"], rofi, pot,
                               metadata["radial_A"], metadata["radial_B"], metadata["radial_Z"], enu)
        err_norm = gintsr_inner(diff, diff, rofi, pot, metadata["radial_A"], metadata["radial_B"],
                                metadata["radial_Z"], enu)
        maxerr = max(abs(v) for row in diff for v in row)
        result["points"].append({"h_ry": h, "input_spacing_ry": 2*h,
            "eval_spacing_ry": plus["eval"] - minus["eval"], "fprime": fprime,
            "plus_eval": plus["eval"], "minus_eval": minus["eval"],
            "plus_nodes": nodes(plus["g"]), "minus_nodes": nodes(minus["g"]),
            "phase_alignment": {"plus": phase_plus, "minus": phase_minus},
            "root_deviation_from_center_max_ry": max(abs(plus["eval"] - center_solve["eval"]),
                                                     abs(minus["eval"] - center_solve["eval"])),
            "gp_vs_fd_max_abs": maxerr,
            "gp_vs_fd_weighted_l2": math.sqrt(max(err_norm, 0.0)),
            "gp_weighted_l2": math.sqrt(max(gp_norm, 0.0))})
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "output_sha256": sha(args.output),
                      "fprime": [x["fprime"] for x in result["points"]]}, sort_keys=True))


if __name__ == "__main__":
    main()
