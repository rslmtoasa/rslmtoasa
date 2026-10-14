#!/usr/bin/env python3
"""Check the scalar-relativistic radial eigen-equation on production captures.

This is an independent numerical contraction of the radial fields emitted by
the production NEWRHO/RSEQSR/PHDFSR path.  It is not a total-energy oracle.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def val(row: dict[str, str], key: str) -> float:
    return float(row[key])


def simpson_log(values: list[float], radius: list[float], a: float, b: float, nre: int) -> float:
    """Production GINTSR/ATOMSC Simpson weights, through Fortran NRE."""
    if nre < 3 or nre > len(values) or nre % 2 != 1:
        raise ValueError(f"invalid odd NRE={nre} for {len(values)} radial samples")
    total = 0.0
    # Fortran IR=2,4,...,NRE-1 has Simpson coefficient 4; IR=3,5,... has 2.
    for i in range(1, nre - 1):
        coefficient = 4.0 if (i + 1) % 2 == 0 else 2.0
        total += coefficient * (radius[i] + b) * values[i]
    total += (radius[nre - 1] + b) * values[nre - 1]
    return a * total / 3.0


def simpson_log_stride2(values: list[float], radius: list[float], a: float, b: float, nre: int) -> float:
    """Coarsened quadrature on every other captured log-mesh point.

    The interval count after decimation is odd for the odd NRE values used by
    RSEQSR, so Simpson 1/3 is applied on all but the final three intervals;
    Simpson 3/8 closes at the same captured sphere boundary.
    """
    if nre < 7 or nre > len(values) or nre % 2 != 1:
        raise ValueError(f"invalid odd NRE={nre} for stride-2 quadrature")
    indices = list(range(0, nre, 2))
    if indices[-1] != nre - 1:
        raise ValueError("stride-2 quadrature must include sphere endpoint")
    vals = [(radius[i] + b) * values[i] for i in indices]
    intervals = len(indices) - 1
    prefix_intervals = intervals - 3
    if prefix_intervals % 2:
        raise ValueError("stride-2 Simpson prefix must have an even interval count")
    total = vals[0] + vals[prefix_intervals]
    total += 4.0 * sum(vals[i] for i in range(1, prefix_intervals, 2))
    total += 2.0 * sum(vals[i] for i in range(2, prefix_intervals, 2))
    total *= (2.0 * a) / 3.0
    # Simpson 3/8 over the final three coarsened intervals.
    total += (3.0 * 2.0 * a / 8.0) * (
        vals[prefix_intervals] + 3.0 * vals[prefix_intervals + 1] +
        3.0 * vals[prefix_intervals + 2] + vals[prefix_intervals + 3]
    )
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace_csv", type=Path)
    parser.add_argument("--site", type=int, default=0)
    parser.add_argument("--c-light", type=float, default=274.074)
    parser.add_argument("--identity-tolerance-ry", type=float, default=3e-6)
    args = parser.parse_args()

    rows = [r for r in read_rows(args.trace_csv) if r["stage"] == "RADIAL_REPLAY" and int(r["site"]) == args.site]
    if not rows:
        raise SystemExit(f"no RADIAL_REPLAY rows for site {args.site}")
    by_kind: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_kind[row["kind"]].append(row)
    meta = {r["kind"]: val(r, "scalar") for kind in (
        "radial_A", "radial_B", "radial_Z", "radial_NR"
    ) for r in by_kind[kind]}
    required = {"radial_A", "radial_B", "radial_Z", "radial_NR"}
    if required - meta.keys():
        raise SystemExit(f"missing metadata {sorted(required - meta.keys())}")
    a, b, z = (meta["radial_A"], meta["radial_B"], meta["radial_Z"])
    nr = int(round(meta["radial_NR"]))
    radius = [math.nan] * nr
    for row in by_kind["radial_replay_mesh"]:
        radius[int(row["row"])] = val(row, "energy")
    if any(not math.isfinite(x) for x in radius):
        raise SystemExit("radial mesh capture is incomplete")

    potential = {(int(r["row"]), int(r["spin"])): val(r, "scalar")
                 for r in by_kind["radial_replay_v"]}
    ql = {(int(r["l"]), int(r["spin"]), int(r["order"])): val(r, "scalar")
          for r in by_kind["radial_replay_ql"]}
    evals = {(int(r["l"]), int(r["spin"])): val(r, "scalar")
             for r in by_kind["radial_replay_eigenvalue"]}
    nre = {(int(r["l"]), int(r["spin"])): int(round(val(r, "scalar")))
           for r in by_kind["radial_replay_nre"]}
    wave = {(int(r["l"]), int(r["spin"]), int(r["row"]), int(r["col"])): val(r, "scalar")
            for r in by_kind["radial_replay_wave"]}
    rho_core = {(int(r["row"]), int(r["spin"])): val(r, "scalar")
                for r in by_kind["radial_replay_rho_core"]}
    rho_out = {(int(r["row"]), int(r["spin"])): val(r, "scalar")
               for r in by_kind["radial_replay_rho_out"]}
    sumev_rows = {int(r["spin"]): val(r, "scalar") for r in by_kind["radial_replay_sumev"]}
    active = [(l, s) for (l, s), q0n in ((k[:2], k[2]) for k in ql) if q0n == 0 and ql[(l, s, 0)] >= 1e-5]
    if not active or not nre or not rho_out or not rho_core:
        raise SystemExit("capture lacks active channels, NRE, or radial density")

    per_channel: list[dict[str, float | int]] = []
    lhs_weighted = rhs_weighted = sumev_q0 = sum_q1 = 0.0
    q0_potential_integral = 0.0
    q0_kinetic = q0_surface = q0_norm_correction = q0_potential_metric_correction = 0.0
    moment_density = {(i, spin): 0.0 for i in range(nr) for spin in {s for _, s in nre}}
    for l, spin in sorted(set(active)):
        cutoff = nre[(l, spin)]
        if (l, spin) not in evals or any((l, spin, order) not in ql for order in range(3)):
            raise SystemExit(f"missing radial inputs for l={l}, spin={spin}")
        u = [wave[(l, spin, i, 0)] for i in range(nr)]
        f = [wave[(l, spin, i, 3)] for i in range(nr)]
        energy = evals[(l, spin)]
        angular = float(l * (l + 1))
        wpot: list[float] = []
        tmc: list[float] = []
        for i, r in enumerate(radius):
            w = potential[(i, spin)] - 2.0 * z / r if r > 0.0 else 0.0
            wpot.append(w)
            tmc.append(args.c_light + (energy - w) / args.c_light)

        # From RSEQSR/RSQSR1/RSQSR2 (the source's first-order SR system):
        # u' = t f + u/r; f' = -f/r + [l(l+1)/(t*r^2) - (E-W)/c]u.
        # Add u*f' and f*u' and integrate over the finite ASA sphere.
        norm_dirac = simpson_log([u[i] ** 2 - f[i] ** 2 for i in range(nr)], radius, a, b, cutoff)
        source_norm = simpson_log([
            (1.0 + angular / (tmc[i] * radius[i]) ** 2) * u[i] ** 2 + f[i] ** 2
            if radius[i] > 0.0 else u[i] ** 2 + f[i] ** 2 for i in range(nr)
        ], radius, a, b, cutoff)
        potential_dirac = simpson_log([wpot[i] * (u[i] ** 2 - f[i] ** 2) for i in range(nr)],
                                      radius, a, b, cutoff)
        potential_source_metric = simpson_log([wpot[i] * (
            (1.0 + angular / (tmc[i] * radius[i]) ** 2) * u[i] ** 2 + f[i] ** 2
            if radius[i] > 0.0 else u[i] ** 2 + f[i] ** 2
        ) for i in range(nr)], radius, a, b, cutoff)
        kinetic_sr = simpson_log([
            args.c_light ** 2 * f[i] ** 2 +
            (args.c_light * angular * u[i] ** 2 / (tmc[i] * radius[i] ** 2) if radius[i] > 0.0 else 0.0)
            for i in range(nr)
        ], radius, a, b, cutoff)
        surface = -args.c_light * u[cutoff - 1] * f[cutoff - 1]
        lhs = energy * norm_dirac
        rhs = potential_dirac + kinetic_sr + surface
        residual = lhs - rhs

        q0 = ql[(l, spin, 0)]
        q1 = ql[(l, spin, 1)]
        q2 = ql[(l, spin, 2)]
        sumev_q0 += energy * q0
        sum_q1 += q1
        lhs_weighted += q0 * lhs
        rhs_weighted += q0 * rhs
        # This is the Q0-only potential matrix element; production rho also
        # contains Q2 terms, reported separately below.
        q0_potential_integral += q0 * potential_source_metric
        q0_kinetic += q0 * kinetic_sr
        q0_surface += q0 * surface
        q0_norm_correction += q0 * energy * (source_norm - norm_dirac)
        q0_potential_metric_correction += q0 * (potential_dirac - potential_source_metric)
        q1_cross = [wave[(l, spin, i, 0)] * wave[(l, spin, i, 1)] * (
            1.0 + angular / (tmc[i] * radius[i]) ** 2
        ) + wave[(l, spin, i, 3)] * wave[(l, spin, i, 4)] if radius[i] > 0.0 else
            wave[(l, spin, i, 0)] * wave[(l, spin, i, 1)] + wave[(l, spin, i, 3)] * wave[(l, spin, i, 4)]
            for i in range(nr)]
        q2_cross = [
            (1.0 + angular / (tmc[i] * radius[i]) ** 2) * (
                wave[(l, spin, i, 1)] ** 2 + wave[(l, spin, i, 0)] * wave[(l, spin, i, 2)]
            ) + wave[(l, spin, i, 4)] ** 2 + wave[(l, spin, i, 3)] * wave[(l, spin, i, 5)]
            if radius[i] > 0.0 else
            wave[(l, spin, i, 1)] ** 2 + wave[(l, spin, i, 0)] * wave[(l, spin, i, 2)] +
            wave[(l, spin, i, 4)] ** 2 + wave[(l, spin, i, 3)] * wave[(l, spin, i, 5)]
            for i in range(nr)
        ]
        q0_cross = [(1.0 + angular / (tmc[i] * radius[i]) ** 2) * u[i] ** 2 + f[i] ** 2
                    if radius[i] > 0.0 else u[i] ** 2 + f[i] ** 2 for i in range(nr)]
        for i in range(cutoff):
            moment_density[(i, spin)] += q0 * q0_cross[i] + 2.0 * q1 * q1_cross[i] + q2 * q2_cross[i]
        per_channel.append({
            "l": l, "spin_zero_based": spin, "q0": q0, "q1_ry": q1,
            "energy_ry": energy, "source_charge_norm": source_norm,
            "u2_minus_f2_norm": norm_dirac, "identity_lhs_ry": lhs,
            "identity_rhs_ry": rhs, "identity_residual_ry": residual,
            "kinetic_term_ry": kinetic_sr, "potential_term_with_u2_minus_f2_ry": potential_dirac,
            "potential_term_with_production_metric_ry": potential_source_metric,
            "finite_sphere_surface_term_ry": surface,
        })

    # Integrate the actual captured production valence density (rho_out-rho_core)
    # against W=V-2Z/r, using the ATOMSC Simpson convention on all NR points.
    actual_w_rho = 0.0
    reconstructed_w_rho = 0.0
    for spin in sorted(sumev_rows):
        radial = []
        reconstructed = []
        for i, r in enumerate(radius):
            w = potential[(i, spin)] - 2.0 * z / r if r > 0.0 else 0.0
            radial.append(w * (rho_out[(i, spin)] - rho_core[(i, spin)]))
            reconstructed.append(w * moment_density[(i, spin)])
        actual_w_rho += simpson_log(radial, radius, a, b, nr)
        reconstructed_w_rho += simpson_log(reconstructed, radius, a, b, nr)

    q0_radial_difference = sumev_q0 - q0_potential_integral
    q0_identity_difference = (q0_kinetic + q0_surface + q0_norm_correction +
                               q0_potential_metric_correction)

    max_residual = max(abs(float(x["identity_residual_ry"])) for x in per_channel)
    coarse_residuals = []
    for l, spin in sorted(set(active)):
        cutoff = nre[(l, spin)]
        energy = evals[(l, spin)]
        angular = float(l * (l + 1))
        u = [wave[(l, spin, i, 0)] for i in range(nr)]
        f = [wave[(l, spin, i, 3)] for i in range(nr)]
        wpot = [potential[(i, spin)] - 2.0 * z / radius[i] if radius[i] > 0.0 else 0.0
                for i in range(nr)]
        tmc = [args.c_light + (energy - wpot[i]) / args.c_light for i in range(nr)]
        lhs_coarse = energy * simpson_log_stride2([u[i] ** 2 - f[i] ** 2 for i in range(nr)], radius, a, b, cutoff)
        potential_coarse = simpson_log_stride2([wpot[i] * (u[i] ** 2 - f[i] ** 2) for i in range(nr)],
                                                radius, a, b, cutoff)
        kinetic_coarse = simpson_log_stride2([
            args.c_light ** 2 * f[i] ** 2 +
            (args.c_light * angular * u[i] ** 2 / (tmc[i] * radius[i] ** 2) if radius[i] > 0.0 else 0.0)
            for i in range(nr)
        ], radius, a, b, cutoff)
        boundary = -args.c_light * u[cutoff - 1] * f[cutoff - 1]
        coarse_residuals.append(lhs_coarse - potential_coarse - kinetic_coarse - boundary)
    source_sumev = sum(sumev_rows.values())
    result = {
        "capture": str(args.trace_csv), "site": args.site, "source_equations": "RSEQSR/RSQSR1/RSQSR2",
        "radial_mesh_points": nr, "mesh_a": a, "mesh_b": b, "nuclear_z": z,
        "speed_of_light_ry_units": args.c_light,
        "identity": "E*int(u^2-f^2)=int[W(u^2-f^2)+c^2 f^2+c*l(l+1)u^2/(t*r^2)]-c*u(R)f(R)",
        "max_channel_identity_residual_ry": max_residual,
        "max_channel_stride2_quadrature_identity_residual_ry": max(abs(x) for x in coarse_residuals),
        "max_channel_full_vs_stride2_identity_residual_change_ry": max(
            abs(coarse_residuals[i] - float(per_channel[i]["identity_residual_ry"]))
            for i in range(len(per_channel))
        ),
        "sum_e_q0_ry": sumev_q0, "sum_q1_ry": sum_q1,
        "captured_sumev_ry": source_sumev,
        "sumev_minus_q0_e_minus_q1_ry": source_sumev - sumev_q0 - sum_q1,
        "q0_weighted_identity_lhs_ry": lhs_weighted,
        "q0_weighted_identity_rhs_ry": rhs_weighted,
        "q0_weighted_source_metric_potential_ry": q0_potential_integral,
        "q0_weighted_scalar_relativistic_kinetic_term_ry": q0_kinetic,
        "q0_weighted_finite_sphere_surface_term_ry": q0_surface,
        "q0_weighted_norm_metric_correction_ry": q0_norm_correction,
        "q0_weighted_potential_metric_correction_ry": q0_potential_metric_correction,
        "q0_sum_ev_minus_source_metric_potential_ry": q0_radial_difference,
        "q0_eigen_equation_decomposition_residual_ry": q0_radial_difference - q0_identity_difference,
        "actual_valence_density_integral_Wrho_ry": actual_w_rho,
        "moment_reconstructed_valence_density_integral_Wrho_ry": reconstructed_w_rho,
        "captured_vs_reconstructed_Wrho_residual_ry": actual_w_rho - reconstructed_w_rho,
        "sumev_minus_actual_valence_Wrho_ry": source_sumev - actual_w_rho,
        "q1_q2_density_potential_increment_ry": actual_w_rho - q0_potential_integral,
        "per_channel": per_channel,
        "claim_limit": (
            "Checks the source scalar-relativistic radial eigen-equation identity with its finite-sphere boundary term. "
            "The eigen-equation norm (u^2-f^2) differs from the production charge metric "
            "(GFAC*u^2+f^2); production density also contains Q2 energy-derivative terms. The stride-2 "
            "quadrature is a coarsening check on captured wavefunctions, not a radial-solver mesh convergence test. "
            "This does not identify SUMEV minus the actual density-potential integral with a complete kinetic functional, "
            "and does not settle the Hubbard double-counting correction."
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if max_residual > args.identity_tolerance_ry:
        raise SystemExit(f"radial eigen-equation identity residual {max_residual:.3e} exceeds {args.identity_tolerance_ry:.3e} Ry")
    if abs(result["sumev_minus_q0_e_minus_q1_ry"]) > 2e-10:
        raise SystemExit("captured SUMEV does not reconstruct from EVAL*Q0+Q1")


if __name__ == "__main__":
    main()
