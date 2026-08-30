#!/usr/bin/env python3
"""Run and plot the XCR-01 fixed-density XC reconciliation sweep.

The Fortran executable is the production-path driver.  This script only
orchestrates it, summarizes the CSV evidence, and makes diagnostic plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


PAIR_LABELS = {
    "BH_VBH": "legacy TXC=1 / libXC TXC=101",
    "PW92_PW": "legacy TXC=5 / libXC TXC=105",
    "VWN_VWN": "legacy TXC=4 / libXC TXC=106",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, skipinitialspace=True))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key].strip())


def make_plots(output_dir: Path, comparison: list[dict[str, str]], response: list[dict[str, str]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    regular = [row for row in comparison if row["evaluation_status"].strip() == "REGULAR"]
    pairs = list(PAIR_LABELS)
    rs_values = sorted({number(row, "rs") for row in regular})

    fig, axes = plt.subplots(len(pairs), 1, figsize=(7.5, 10.5), sharex=True)
    if len(pairs) == 1:
        axes = [axes]
    for axis, pair in zip(axes, pairs):
        rows = [row for row in regular if row["functional_pair"].strip() == pair]
        for rs in rs_values:
            curve = sorted((row for row in rows if number(row, "rs") == rs), key=lambda r: number(r, "zeta"))
            axis.plot(
                [number(row, "zeta") for row in curve],
                [number(row, "bxc_legacy") for row in curve],
                marker="o",
                ms=2.5,
                label=fr"$r_s={rs:g}$",
            )
            axis.plot(
                [number(row, "zeta") for row in curve],
                [number(row, "bxc_libxc") for row in curve],
                linestyle="--",
                linewidth=0.8,
                color=axis.lines[-1].get_color(),
            )
        axis.set_ylabel(r"$B_{xc}$ (Ry)")
        axis.set_title(PAIR_LABELS[pair])
        axis.grid(alpha=0.25)
    axes[-1].set_xlabel(r"polarization $\zeta$")
    axes[0].legend(ncol=3, fontsize="small")
    fig.tight_layout()
    fig.savefig(output_dir / "bxc_vs_zeta.png", dpi=180)
    plt.close(fig)

    response = [row for row in response if row["functional_pair"].strip() in pairs]
    available_deltas = {number(row, "delta") for row in response if number(row, "delta") > 0}
    delta = 1.0e-4 if 1.0e-4 in available_deltas else min(available_deltas)
    fig, axes = plt.subplots(len(pairs), 1, figsize=(7.5, 10.5), sharex=True)
    if len(pairs) == 1:
        axes = [axes]
    for axis, pair in zip(axes, pairs):
        rows = [
            row
            for row in response
            if row["functional_pair"].strip() == pair and abs(number(row, "delta") - delta) < 1.0e-14
        ]
        rows.sort(key=lambda row: number(row, "rs"))
        axis.plot([number(row, "rs") for row in rows], [number(row, "i_xc_legacy") for row in rows], "o-", label="legacy")
        axis.plot([number(row, "rs") for row in rows], [number(row, "i_xc_libxc") for row in rows], "s--", label="libXC")
        axis.set_ylabel(r"$I_{xc}$ (Ry bohr$^3$)")
        axis.set_title(f"{PAIR_LABELS[pair]}, symmetric delta={delta:g}")
        axis.grid(alpha=0.25)
    axes[-1].set_xlabel(r"$r_s$ (bohr)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(output_dir / "stoner_response.png", dpi=180)
    plt.close(fig)


def summarize(
    comparison: list[dict[str, str]],
    response: list[dict[str, str]],
    exchange: list[dict[str, str]],
    energy_gradient: list[dict[str, str]],
) -> dict:
    regular = [row for row in comparison if row["evaluation_status"].strip() == "REGULAR"]
    fully_polarized = [row for row in comparison if row["evaluation_status"].strip() != "REGULAR"]
    summary: dict[str, object] = {
        "comparison_rows": len(comparison),
        "regular_comparison_rows": len(regular),
        "fully_polarized_rows": len(comparison) - len(regular),
        "exchange_validation": {
            "max_absolute_error": max(number(row, "max_abs_error") for row in exchange),
            "max_relative_error": max(number(row, "max_relative_error") for row in exchange),
        },
        "fully_polarized_limits": {},
        "pairs": {},
    }
    pair_summary = summary["pairs"]
    assert isinstance(pair_summary, dict)
    limit_summary = summary["fully_polarized_limits"]
    assert isinstance(limit_summary, dict)
    for pair in PAIR_LABELS:
        rows = [row for row in regular if row["functional_pair"].strip() == pair]
        limit_rows = [row for row in fully_polarized if row["functional_pair"].strip() == pair]
        limit_summary[pair] = {
            "rows": len(limit_rows),
            "max_exc_absolute_difference": max(number(row, "exc_abs_difference") for row in limit_rows),
            "max_vup_absolute_difference": max(number(row, "vup_abs_difference") for row in limit_rows),
            "max_vdn_absolute_difference": max(number(row, "vdn_abs_difference") for row in limit_rows),
            "max_bxc_absolute_difference": max(number(row, "bxc_abs_difference") for row in limit_rows),
            "classifications": sorted({row["classification"].strip() for row in limit_rows}),
        }
        pair_summary[pair] = {
            "rows": len(rows),
            "max_exc_absolute_difference": max(number(row, "exc_abs_difference") for row in rows),
            "max_vup_absolute_difference": max(number(row, "vup_abs_difference") for row in rows),
            "max_vdn_absolute_difference": max(number(row, "vdn_abs_difference") for row in rows),
            "max_bxc_absolute_difference": max(number(row, "bxc_abs_difference") for row in rows),
            "max_relative_difference": max(number(row, "relative_difference") for row in rows),
            "classifications": sorted({row["classification"].strip() for row in rows}),
        }
        pair_response = [row for row in response if row["functional_pair"].strip() == pair]
        pair_summary[pair]["stoner_at_delta_1e-4"] = [
            {
                "rs": number(row, "rs"),
                "i_xc_legacy": number(row, "i_xc_legacy"),
                "i_xc_libxc": number(row, "i_xc_libxc"),
                "absolute_difference": number(row, "absolute_difference"),
                "relative_difference": number(row, "relative_difference"),
            }
            for row in pair_response
            if abs(number(row, "delta") - 1.0e-4) < 1.0e-14
        ]
        pair_gradient = [row for row in energy_gradient if row["functional_pair"].strip() == pair]
        pair_summary[pair]["max_legacy_energy_gradient_residual"] = max(
            number(row, "max_residual") for row in pair_gradient
        )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path, help="libXC-enabled UnitXcLdaReconciliation executable")
    parser.add_argument("--output-dir", required=True, type=Path, help="directory for CSV, JSON, plots, and run log")
    parser.add_argument("--no-plots", action="store_true", help="skip matplotlib diagnostic plots")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / "run.log"
    with log_path.open("w") as log:
        completed = subprocess.run(
            [str(args.binary), str(args.output_dir)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise SystemExit(f"{args.binary} failed; see {log_path}")

    comparison = read_csv(args.output_dir / "xc_lda_reconciliation.csv")
    response = read_csv(args.output_dir / "xc_lda_stoner_response.csv")
    exchange = read_csv(args.output_dir / "xc_lda_exchange_validation.csv")
    energy_gradient = read_csv(args.output_dir / "xc_lda_energy_gradient.csv")
    summary = summarize(comparison, response, exchange, energy_gradient)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not args.no_plots:
        make_plots(args.output_dir, comparison, response)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
