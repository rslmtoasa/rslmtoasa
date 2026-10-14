#!/usr/bin/env python3
"""Independently check captured production NEWRHO moments and radial density.

This validates the captured source-connected radial calculation; it does not
establish the missing Hubbard total-energy functional identity.
"""

import argparse
import csv
import json
import math
from collections import defaultdict


def read_trace(path):
    with open(path, newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trace_csv")
    parser.add_argument("--site", type=int, default=0)
    parser.add_argument("--c-light", type=float, default=274.074)
    parser.add_argument("--tolerance", type=float, default=2e-10)
    args = parser.parse_args()

    rows = read_trace(args.trace_csv)
    replay = [r for r in rows if r["stage"] == "RADIAL_REPLAY" and int(r["site"]) == args.site]
    if not replay:
        raise SystemExit("no RADIAL_REPLAY rows for requested site")

    meta = {r["kind"]: number(r, "scalar") for r in replay if r["kind"] in {"radial_A", "radial_B", "radial_tolerance", "radial_NR", "radial_Z"}}
    required_meta = {"radial_A", "radial_B", "radial_tolerance", "radial_NR", "radial_Z"}
    if required_meta - meta.keys():
        raise SystemExit(f"missing replay metadata: {sorted(required_meta - meta.keys())}")

    by_kind = defaultdict(list)
    for row in replay:
        by_kind[row["kind"]].append(row)

    # Radial mesh and per-spin inputs/outputs are indexed from zero in the CSV.
    radii = {int(r["row"]): number(r, "energy") for r in by_kind["radial_replay_mesh"]}
    potentials = {(int(r["row"]), int(r["spin"])): number(r, "scalar") for r in by_kind["radial_replay_v"]}
    rho_in = {(int(r["row"]), int(r["spin"])): number(r, "scalar") for r in by_kind["radial_replay_rho_in"]}
    rho_out = {(int(r["row"]), int(r["spin"])): number(r, "scalar") for r in by_kind["radial_replay_rho_out"]}
    rho_core = {(int(r["row"]), int(r["spin"])): number(r, "scalar") for r in by_kind["radial_replay_rho_core"]}
    ql = {(int(r["l"]), int(r["spin"]), int(r["order"])): number(r, "scalar") for r in by_kind["radial_replay_ql"]}
    eigenvalues = {(int(r["l"]), int(r["spin"])): number(r, "scalar") for r in by_kind["radial_replay_eigenvalue"]}
    nre = {(int(r["l"]), int(r["spin"])): int(round(number(r, "scalar"))) for r in by_kind["radial_replay_nre"]}
    replay_residuals = {
        kind.removeprefix("radial_replay_"): number(entries[0], "scalar")
        for kind, entries in by_kind.items()
        if kind.startswith("radial_replay_") and kind.endswith(("_max_abs", "_max_abs_ry"))
    }
    wave = {}
    for row in by_kind["radial_replay_wave"]:
        key = (int(row["l"]), int(row["spin"]), int(row["row"]), int(row["col"]))
        wave[key] = number(row, "scalar")

    if not nre:
        raise SystemExit("missing radial_replay_nre rows; NEWRHO integrates only through NRE")
    if not eigenvalues or not ql or not wave or not rho_core:
        raise SystemExit("incomplete NEWRHO replay payload (need QL, eigenvalues, waves, and core density)")
    for (ell, spin), cutoff in nre.items():
        if (ell, spin) not in eigenvalues or any((ell, spin, order) not in ql for order in range(3)):
            raise SystemExit(f"missing channel inputs for l={ell}, spin={spin}")
        for ir0 in range(1, cutoff):
            for component in range(6):
                if (ell, spin, ir0, component) not in wave:
                    raise SystemExit(f"missing wave row l={ell}, spin={spin}, radial_index={ir0}, component={component}")

    spins = sorted({s for _, s in rho_out})
    max_rho_error = 0.0
    max_sumev_error = 0.0
    per_spin = {}
    for spin in spins:
        predicted = defaultdict(float)
        predicted_sumev = 0.0
        for ell in sorted(l for l, s, _ in ql if s == spin and _ == 0):
            q0 = ql[(ell, spin, 0)]
            q1 = ql[(ell, spin, 1)]
            q2 = ql[(ell, spin, 2)]
            if q0 < 1e-5:
                continue
            energy = eigenvalues[(ell, spin)]
            predicted_sumev += energy * q0 + q1
            last_index = nre[(ell, spin)] - 1
            for ir0 in range(1, last_index + 1):
                radius = radii[ir0]
                v = potentials[(ir0, spin)]
                tmc = args.c_light - (v - 2.0 * meta["radial_Z"] / radius - energy) / args.c_light
                gfac = 1.0 + ell * (ell + 1.0) / (tmc * radius) ** 2
                g = wave[(ell, spin, ir0, 0)]
                gp = wave[(ell, spin, ir0, 1)]
                gpp = wave[(ell, spin, ir0, 2)]
                gs = wave[(ell, spin, ir0, 3)]
                gps = wave[(ell, spin, ir0, 4)]
                gpps = wave[(ell, spin, ir0, 5)]
                predicted[ir0] += (
                    q0 * (gfac * g * g + gs * gs)
                    + 2.0 * q1 * (gfac * g * gp + gs * gps)
                    + q2 * (gfac * (gp * gp + g * gpp) + gps * gps + gs * gpps)
                )

        for ir0, valence_density in predicted.items():
            actual = rho_out[(ir0, spin)] - rho_core[(ir0, spin)]
            max_rho_error = max(max_rho_error, abs(actual - valence_density))
        actual_sumev = number(next(r for r in by_kind["radial_replay_sumev"] if int(r["spin"]) == spin), "scalar")
        max_sumev_error = max(max_sumev_error, abs(actual_sumev - predicted_sumev))
        per_spin[str(spin)] = {"sum_ev_actual_ry": actual_sumev, "sum_ev_reconstructed_ry": predicted_sumev,
                               "nre_by_l": {str(l): nre[(l, spin)] for l, s in nre if s == spin}}

    result = {"site": args.site, "source_connected_capture": True,
              "max_valence_density_abs_error": max_rho_error,
              "max_sumev_abs_error_ry": max_sumev_error,
              "same_process_ne_wrho_replay_residuals": replay_residuals,
              "tolerance": args.tolerance, "per_spin": per_spin,
              "claim_limit": "Checks NEWRHO source equations against captured production arrays; does not prove Hubbard ETOT variational accounting."}
    print(json.dumps(result, indent=2, sort_keys=True))
    if max_rho_error > args.tolerance or max_sumev_error > args.tolerance:
        raise SystemExit("captured radial density/SUMEV reconstruction failed")
    if not replay_residuals:
        raise SystemExit("missing in-process production NEWRHO replay residuals")
    if max(replay_residuals.values()) > args.tolerance:
        raise SystemExit("same-state production NEWRHO replay did not reproduce its outputs")


if __name__ == "__main__":
    main()
