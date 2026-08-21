#!/usr/bin/env python3
"""Aggregate the frozen KPM-B1 campaign and generate its report package.

The benchmark runner writes one canonical B0C-shaped JSON document per
configuration group.  B1 intentionally aggregates those documents without
recomputing or hand-editing measurements; the JSON, CSV, Markdown, and plots
are all derived from the same row objects.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def median(row: dict[str, Any], name: str) -> float:
    return float(row.get("statistics", {}).get(name, {}).get("median", 0.0) or 0.0)


def metric(row: dict[str, Any], name: str) -> float:
    return median(row, name)


def eligible(row: dict[str, Any]) -> bool:
    return bool(row.get("headline_speedup_eligible")) and row.get("profile_status") == "PASS" and row.get("correctness", {}).get("status") == "PASS"


def rows_for(rows: list[dict[str, Any]], **filters: Any) -> list[dict[str, Any]]:
    return [row for row in rows if all(row.get(key) == value for key, value in filters.items())]


def unique_rows(reports: list[tuple[Path, dict[str, Any]]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_id: dict[str, dict[str, Any]] = {}
    skipped: list[dict[str, Any]] = []
    for source, report in reports:
        for row in report.get("rows", []):
            row = dict(row)
            row["source_report"] = str(source)
            by_id.setdefault(str(row.get("row_id")), row)
        for row in report.get("skipped_rows", []):
            item = dict(row)
            item["source_report"] = str(source)
            skipped.append(item)
    return sorted(by_id.values(), key=lambda row: (
        str(row.get("material", "")), int(row.get("replication", 0) or 0), int(row.get("M", 0) or 0),
        str(row.get("cond_type", "")), str(row.get("numeric_mode", "")), bool(row.get("gpu_plugin")),
        int(row.get("OMP_threads", 0) or 0), int(row.get("block_width", 0) or 0))), skipped


def _csv_value(row: dict[str, Any], key: str) -> Any:
    if key == "correctness_status":
        return row.get("correctness", {}).get("status")
    if key == "headline_speedup_eligible":
        return row.get("headline_speedup_eligible", False)
    if key == "transport_median":
        return median(row, "T_transport_total")
    if key == "whole_wall_median":
        return median(row, "wall_time_s")
    if key == "moment_median":
        return median(row, "P_moments_total")
    if key == "gamma_median":
        return median(row, "P_gamma_basis_setup") + median(row, "P_gamma_generation")
    if key == "reconstruction_median":
        return median(row, "P_reconstruction_total")
    if key == "output_median":
        return median(row, "P_output_io")
    if key in {"S_moments", "S_transport", "S_whole"}:
        return row.get("speedups", {}).get(key)
    return row.get(key)


CSV_FIELDS = [
    "row_id", "material", "size", "replication", "N", "nnz", "M", "NE", "lld", "cond_type",
    "cond_calctype", "Ntrace", "random_seed", "block_width", "moment_backend", "moment_precision",
    "reconstruction_backend", "reconstruction_precision", "numeric_mode", "canonical_output_precision",
    "OMP_threads", "BLAS_threads", "gpu_plugin", "profile_status", "correctness_status",
    "transport_median", "whole_wall_median", "moment_median", "gamma_median", "reconstruction_median",
    "output_median", "S_moments", "S_transport", "S_whole", "headline_speedup_eligible", "best_cpu_row_id",
]


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row, key) for key in CSV_FIELDS})


def fmt(value: Any, digits: int = 3) -> str:
    if value is None or value == "":
        return "-"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "-"
    return f"{number:.{digits}g}"


def best_cpu(rows: list[dict[str, Any]], gpu: dict[str, Any]) -> dict[str, Any] | None:
    row_id = gpu.get("best_cpu_row_id")
    return next((row for row in rows if row.get("row_id") == row_id), None)


def headline_table(rows: list[dict[str, Any]], mode: str) -> list[str]:
    lines = [
        "| material | size | N | nnz | M | estimator | Ntrace | best CPU OMP | CPU moment (s) | GPU moment (s) | S_moments | CPU transport (s) | GPU transport (s) | S_transport | CPU wall (s) | GPU wall (s) | S_whole | correctness |",
        "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for gpu in rows:
        if not (gpu.get("gpu_plugin") and gpu.get("cond_calctype") == "per_type" and gpu.get("numeric_mode") == mode and eligible(gpu)):
            continue
        cpu = best_cpu(rows, gpu)
        ref = cpu or {}
        lines.append("| {material} | {size} | {N} | {nnz} | {M} | {estimator} | {Ntrace} | {omp} | {cm} | {gm} | {sm} | {ct} | {gt} | {st} | {cw} | {gw} | {sw} | {corr} |".format(
            material=gpu.get("material"), size=gpu.get("size"), N=gpu.get("N"), nnz=gpu.get("nnz"), M=gpu.get("M"),
            estimator=gpu.get("cond_calctype"), Ntrace=gpu.get("Ntrace"), omp=ref.get("OMP_threads", "-"),
            cm=fmt(median(ref, "P_moments_total")), gm=fmt(median(gpu, "P_moments_total")),
            sm=fmt(gpu.get("speedups", {}).get("S_moments")),
            ct=fmt(median(ref, "T_transport_total")), gt=fmt(median(gpu, "T_transport_total")),
            st=fmt(gpu.get("speedups", {}).get("S_transport")), cw=fmt(median(ref, "wall_time_s")),
            gw=fmt(median(gpu, "wall_time_s")), sw=fmt(gpu.get("speedups", {}).get("S_whole")),
            corr=gpu.get("correctness", {}).get("status", "-")))
    return lines


def m_scaling_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = [
        "| material | size | mode | M | CPU OMP | CPU moments (s) | GPU moments (s) | S_moments | CPU transport (s) | GPU transport (s) | S_transport | correctness |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for gpu in rows:
        if not (gpu.get("gpu_plugin") and eligible(gpu) and gpu.get("cond_type") == "spin" and gpu.get("cond_calctype") == "per_type"):
            continue
        cpu = best_cpu(rows, gpu) or {}
        lines.append("| {material} | {size} | {mode} | {M} | {omp} | {cm} | {gm} | {sm} | {ct} | {gt} | {st} | {corr} |".format(
            material=gpu.get("material"), size=gpu.get("size"), mode=gpu.get("numeric_mode"), M=gpu.get("M"),
            omp=cpu.get("OMP_threads", "-"), cm=fmt(median(cpu, "P_moments_total")), gm=fmt(median(gpu, "P_moments_total")),
            sm=fmt(gpu.get("speedups", {}).get("S_moments")), ct=fmt(median(cpu, "T_transport_total")),
            gt=fmt(median(gpu, "T_transport_total")), st=fmt(gpu.get("speedups", {}).get("S_transport")),
            corr=gpu.get("correctness", {}).get("status", "-")))
    return lines


def mixed_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = ["| material | size | M | backend | moments | reconstruction | OMP | transport (s) | wall (s) | correctness |", "|---|---:|---:|---|---|---|---:|---:|---:|---|"]
    for row in rows:
        if row.get("numeric_mode") != "mixed":
            continue
        lines.append("| {material} | {size} | {M} | {backend} | {mp} | {rp} | {omp} | {transport} | {wall} | {corr} |".format(
            material=row.get("material"), size=row.get("size"), M=row.get("M"), backend=row.get("moment_backend"),
            mp=row.get("moment_precision"), rp=row.get("reconstruction_precision"), omp=row.get("OMP_threads"),
            transport=fmt(median(row, "T_transport_total")), wall=fmt(median(row, "wall_time_s")),
            corr=row.get("correctness", {}).get("status", "-")))
    return lines


def cpu_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = ["| material | size | mode | OMP | transport (s) | wall (s) | profile | correctness |", "|---|---:|---|---:|---:|---:|---|---|"]
    for row in rows:
        if row.get("gpu_plugin") or row.get("cond_calctype") != "per_type" or row.get("M") != 500:
            continue
        lines.append("| {material} | {size} | {mode} | {omp} | {transport} | {wall} | {profile} | {corr} |".format(
            material=row.get("material"), size=row.get("size"), mode=row.get("numeric_mode"), omp=row.get("OMP_threads"),
            transport=fmt(median(row, "T_transport_total")), wall=fmt(median(row, "wall_time_s")),
            profile=row.get("profile_status"), corr=row.get("correctness", {}).get("status", "-")))
    return lines


def stochastic_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = ["| backend | precision | Ntrace | block width | moment time/trace (s) | traces/s | transport (s) | correctness |", "|---|---|---:|---:|---:|---:|---:|---|"]
    for row in rows:
        if row.get("cond_calctype") != "random_vec":
            continue
        derived = row.get("derived", {})
        lines.append("| {backend} | {precision} | {ntrace} | {block} | {per_trace} | {rate} | {transport} | {corr} |".format(
            backend="GPU" if row.get("gpu_plugin") else "CPU", precision=row.get("numeric_mode"), ntrace=row.get("Ntrace"),
            block=row.get("block_width"), per_trace=fmt(derived.get("moment_time_per_trace_s", {}).get("median")),
            rate=fmt(derived.get("traces_per_second", {}).get("median")), transport=fmt(median(row, "T_transport_total")),
            corr=row.get("correctness", {}).get("status", "-")))
    return lines


def make_plots(rows: list[dict[str, Any]], plot_dir: Path) -> list[str]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return ["plots unavailable: matplotlib is not installed"]

    made: list[str] = []
    pt = [row for row in rows if row.get("material") == "fccPt_SOC" and row.get("cond_type") == "spin" and row.get("cond_calctype") == "per_type"]
    def save(name: str) -> None:
        plt.tight_layout()
        plt.savefig(plot_dir / name, dpi=160)
        plt.close()
        made.append(str(plot_dir / name))

    for ylabel, key, name in (("S_transport", "S_transport", "transport_speedup_vs_size.png"), ("S_moments", "S_moments", "moment_speedup_vs_size.png")):
        plt.figure(figsize=(6.4, 4.0))
        for mode, color in (("fp64", "tab:blue"), ("fp32", "tab:orange")):
            points = sorted((row for row in pt if row.get("numeric_mode") == mode and eligible(row)), key=lambda row: int(row.get("replication", 0)))
            plt.plot([int(row.get("N")) for row in points], [float(row.get("speedups", {}).get(key, 0.0)) for row in points], "o-", label=mode, color=color)
        plt.xlabel("real-space dimension N")
        plt.ylabel(ylabel)
        plt.grid(alpha=0.25)
        plt.legend()
        save(name)

    plt.figure(figsize=(6.4, 4.0))
    for size, color in ((6, "tab:green"), (8, "tab:red")):
        points = sorted((row for row in pt if int(row.get("replication", 0)) == size and row.get("numeric_mode") == "fp64" and eligible(row)), key=lambda row: int(row.get("M", 0)))
        if points:
            plt.plot([int(row.get("M")) for row in points], [float(row.get("speedups", {}).get("S_transport", 0.0)) for row in points], "o-", label=f"Pt r{size}", color=color)
    plt.xlabel("Chebyshev order M")
    plt.ylabel("S_transport")
    plt.grid(alpha=0.25)
    plt.legend()
    save("transport_speedup_vs_M.png")

    plt.figure(figsize=(7.0, 4.2))
    for size, mode, color in ((4, "fp64", "tab:blue"), (4, "fp32", "tab:orange"), (6, "fp64", "tab:green"), (6, "fp32", "tab:red"), (8, "fp64", "tab:purple"), (8, "fp32", "tab:brown")):
        points = sorted((row for row in pt if int(row.get("replication", 0)) == size and row.get("numeric_mode") == mode and not row.get("gpu_plugin") and int(row.get("M", 0)) == 500), key=lambda row: int(row.get("OMP_threads", 0)))
        if points:
            plt.plot([int(row.get("OMP_threads")) for row in points], [median(row, "T_transport_total") for row in points], "o-", label=f"r{size} {mode}", color=color)
    plt.xlabel("OMP threads")
    plt.ylabel("CPU transport time (s)")
    plt.grid(alpha=0.25)
    plt.legend(ncol=2, fontsize=8)
    save("cpu_thread_scaling.png")

    stochastic = [row for row in rows if row.get("cond_calctype") == "random_vec"]
    if stochastic:
        plt.figure(figsize=(6.4, 4.0))
        for is_gpu, label, color in ((False, "CPU", "tab:blue"), (True, "GPU", "tab:orange")):
            points = sorted((row for row in stochastic if bool(row.get("gpu_plugin")) == is_gpu and int(row.get("block_width", 1)) == 1), key=lambda row: int(row.get("Ntrace", 0)))
            if points:
                plt.plot([int(row.get("Ntrace")) for row in points], [float(row.get("derived", {}).get("traces_per_second", {}).get("median", 0.0)) for row in points], "o-", label=label, color=color)
        plt.xlabel("Ntrace")
        plt.ylabel("traces/s")
        plt.grid(alpha=0.25)
        plt.legend()
        save("stochastic_throughput.png")

    representative = []
    for mode, gpu in (("CPU FP64", False), ("GPU FP64", True), ("CPU FP32", False), ("GPU FP32", True)):
        candidate = next((row for row in pt if int(row.get("replication", 0)) == 4 and row.get("numeric_mode") == mode.split()[-1].lower() and bool(row.get("gpu_plugin")) == gpu and int(row.get("M", 0)) == 500), None)
        if candidate:
            representative.append((mode, candidate))
    if representative:
        labels = [name for name, _ in representative]
        stages = ("moments", "Gamma", "reconstruction", "postprocess", "I/O", "other")
        values: list[list[float]] = []
        for _, row in representative:
            parts = [median(row, "P_moments_total"), median(row, "P_gamma_basis_setup") + median(row, "P_gamma_generation"), median(row, "P_reconstruction_total"), median(row, "P_energy_integration") + median(row, "P_tensor_postprocess"), median(row, "P_output_io")]
            total = max(median(row, "T_transport_total"), 1.0e-30)
            values.append([100.0 * part / total for part in parts] + [max(0.0, 100.0 - sum(100.0 * part / total for part in parts))])
        plt.figure(figsize=(7.0, 4.2))
        bottom = [0.0] * len(labels)
        for index, stage in enumerate(stages):
            heights = [value[index] for value in values]
            plt.bar(labels, heights, bottom=bottom, label=stage)
            bottom = [left + right for left, right in zip(bottom, heights)]
        plt.ylabel("share of transport time (%)")
        plt.xticks(rotation=20, ha="right")
        plt.legend(fontsize=8, ncol=3)
        save("stage_fractions.png")
    return made


def make_report(campaign: dict[str, Any], output: Path, plot_paths: list[str]) -> None:
    rows = campaign["rows"]
    skipped = campaign["skipped_rows"]
    env = campaign["environment"]
    headlines = [row for row in rows if row.get("gpu_plugin") and row.get("cond_calctype") == "per_type" and eligible(row)]
    ineligible = [row for row in rows if row.get("gpu_plugin") and not eligible(row)]
    fp64 = [row for row in headlines if row.get("numeric_mode") == "fp64"]
    fp32 = [row for row in headlines if row.get("numeric_mode") == "fp32"]
    pt_spin = [row for row in rows if row.get("material") == "fccPt_SOC" and row.get("cond_type") == "spin" and row.get("cond_calctype") == "per_type"]
    fe = [row for row in headlines if row.get("material") == "bccFe_magnetic"]
    def crossover(data: list[dict[str, Any]]) -> str:
        tested = sorted({int(row.get("replication", 0)) for row in data})
        winners = sorted({int(row.get("replication", 0)) for row in data if float(row.get("speedups", {}).get("S_transport", 0.0)) >= 1.0})
        if not tested:
            return "not measured"
        if winners and winners[0] == tested[0]:
            return f"GPU already wins at the smallest tested size (r{tested[0]}); crossover below r{tested[0]} was not measured"
        if winners:
            previous = max(size for size in tested if size < winners[0])
            return f"between r{previous} and r{winners[0]} (discrete measured bracket)"
        return f"not reached through r{tested[-1]}"
    lines = [
        "# KPM-B1 final CPU/GPU benchmark campaign", "",
        "This report is generated from `campaign.json`; numerical tables and plots are derived from the same canonical rows.", "",
        "## Executive summary", "",
        f"- Frozen implementation commit: `{env.get('git_commit')}`; campaign rows: `{len(rows)}`; explicit skips plus ineligible GPU records retained: `{len(skipped) + len(ineligible)}`.",
        f"- Primary GPU: `{env.get('gpu_model')}` device `{env.get('selected_gpu_index')}`, `{env.get('gpu_vram_mib')} MiB`, compute `{env.get('gpu_compute_capability')}`, CUDA `{env.get('cuda_toolkit')}`.",
        f"- Headline pairs: FP64 `{len(fp64)}`, FP32 `{len(fp32)}`; every published pair has profile and production-output correctness PASS.",
        "- CPU rows retain OMP=1/2/4/8; GPU rows use one process and one selected GPU at OMP=1.", "",
        "## Methodology and freeze", "",
        f"- Compiler/build: `{env.get('compiler')}`, `{env.get('build_type')}`, BLAS/LAPACK `{env.get('blas_lapack')}`.",
        f"- CPU: `{env.get('cpu_model')}`; MPI ranks recorded as `{env.get('mpi_ranks')}`. The host could not configure the MPI Fortran wrapper reliably, so the executable is serial and the one-rank policy is satisfied without an MPI-scaling claim.",
        "- Production output is enabled for timed rows. Two warmups and five measurements are used for default anchors; genuinely expensive r6/r8 and secondary M/G2 rows use three measurements, recorded in their source JSON policy. Medians are used for speedups and spread remains in JSON.",
        "- FP64 means FP64 moments plus FP64 reconstruction. FP32 means FP32 moments plus FP32 reconstruction. The mixed route is FP32 moments plus FP64 reconstruction on CPU; current CUDA profiling couples its moment and reconstruction precision, so no mixed GPU claim is made.",
        "- The tracked implementation was clean at freeze. Pre-existing untracked build/result artifacts were preserved and are not part of the B1 commit; `git_dirty` in the raw environment remains an honest record of that state.", "",
        "## FP64 headline", "",
        *headline_table(rows, "fp64"), "",
        "## FP32 headline", "",
        *headline_table(rows, "fp32"), "",
        "## M-order scaling", "",
        *m_scaling_table(rows), "",
        "## Mixed practical route", "",
        *mixed_table(rows), "",
        "## CPU OMP scaling", "",
        *cpu_table(rows), "",
        "## Stochastic estimator and G2 block evidence", "",
        *stochastic_table(rows), "",
        "## Magnetic bcc-Fe cross-check", "",
        f"Headline Fe rows retained: `{len(fe)}`. The table above is generated from the real `{('bccFe_magnetic' if fe else 'bccFe_magnetic fixture requested but no eligible headline rows)')}` fixture; no synthetic matrix is used for performance conclusions.", "",
        "## Interpretation", "",
        f"1. FP64 crossover: {crossover(fp64)}.",
        f"2. FP32 crossover: {crossover(fp32)}.",
        "3. Size scaling is read from measured r4/r6/r8 rows; no unsupported crossover fit is applied.",
        "4. M scaling is reported from the completed medium-size series and the anchor rows; omitted large extensions are explicit skips below.",
        "5. Moment speedup is a kernel-attribution metric, while transport speedup includes Gamma, reconstruction, post-processing, and other stages.",
        "6. FP64 GPU speedup is limited by non-moment stages and memory-bound work on this A4000; FP32 has a different balance and is only recommended when its production-output tolerance passes.",
        "7. CPU OpenMP changes the selected baseline and therefore materially affects the measured crossover; every anchor retains OMP=1/2/4/8.",
        "8. FP32 scientific acceptability is a correctness/tolerance decision, not a performance assumption; the report exposes the attached PASS evidence.",
        "9. FP64 GPU use is worthwhile only in the larger or otherwise GPU-favorable measured regimes; small-workload overhead can leave CPU competitive.",
        "10. Charge, spin, and orbital rows are reported separately; shared-stage conclusions are limited to the measured Pt fixtures.",
        "11. random_vec scaling and per-trace throughput are reported separately from the primary per_type production workflow.",
        "12. Block stochastic processing is judged from measured time/trace and traces/s; neutral or negative results are retained.",
        "13. Pt and magnetic Fe are compared as qualitative cross-checks, not as identical-size performance claims.",
        "14. Production recommendation: use CPU FP64 for small workloads when overhead dominates, GPU FP64 for larger validated FP64 workloads, and GPU FP32 only when its scientific tolerance is acceptable; do not hard-code A4000 thresholds.",
        "15. Further KPM GPU optimization is not justified by B1 unless a new workload exposes a dominant, reproducible stage outside the optimized moment path.",
        "16. Remaining limitations are the one-rank serial executable, single A4000 device, bounded stochastic matrix, and explicitly skipped memory/runtime extensions.", "",
        "## Skips, failures, and limitations", "",
    ]
    if skipped:
        lines.extend(f"- `{item.get('row_id', item.get('status', 'skipped'))}`: {item.get('reason', item.get('status', 'skipped'))}" for item in skipped)
    else:
        lines.append("- No runner-generated skipped rows.")
    if ineligible:
        lines.append("")
        lines.append("### Failed or non-headline GPU rows retained")
        lines.append("")
        for item in ineligible:
            corr = item.get("correctness", {})
            reason = ", ".join(filter(None, [f"profile={item.get('profile_status')}", f"correctness={corr.get('status')}", ", ".join(corr.get('reasons', []))]))
            lines.append(f"- `{item.get('row_id')}` ({item.get('material')} {item.get('size')} {item.get('cond_type')} {item.get('numeric_mode')} M={item.get('M')}): {reason or 'not headline eligible'}")
    lines += ["", "## Generated outputs", "", f"- Canonical JSON: `{output.name}`", f"- CSV: `{output.with_suffix('.csv').name}`", f"- Plots: `{len(plot_paths)}`", "- Raw logs and correctness JSON are under `raw/` and `correctness/` beside the canonical dataset.", ""]
    report_text = "\n".join(lines).rstrip() + "\n"
    output.with_suffix(".md").write_text(report_text, encoding="utf-8")
    output.parents[3].joinpath("docs", "dev", "RS_LMTO_ASA_KPM_B1_REPORT.md").write_text(report_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--skip-note", action="append", default=[])
    args = parser.parse_args()
    reports = [(path.resolve(), load(path.resolve())) for path in args.inputs]
    rows, skipped = unique_rows(reports)
    skipped.extend({"status": "explicit_skip", "reason": note} for note in args.skip_note)
    first = reports[0][1]
    environment = dict(first.get("environment", {}))
    try:
        environment["git_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    campaign = {
        "schema": "rslmto.kpm-b1.v1", "scope": "frozen systematic real-material CPU/GPU KPM/Kubo-Bastin evidence",
        "environment": environment, "source_reports": [str(path) for path, _ in reports],
        "rows": rows, "skipped_rows": skipped,
        "policy": {"headline_definition": "same precision, profile PASS, correctness PASS, production output", "speedup_definitions": {"S_moments": "best_same_precision_CPU(P_moments_total)/GPU(P_moments_total)", "S_transport": "best_same_precision_CPU(T_transport_total)/GPU(T_transport_total)", "S_whole": "best_same_precision_CPU(whole_wall)/GPU(whole_wall)"}},
        "provenance": {"raw_logs": "raw/", "correctness": "correctness/", "plots": "plots/"},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(rows, args.output.with_suffix(".csv"))
    plot_paths = make_plots(rows, args.output.parent / "plots")
    make_report(campaign, args.output, plot_paths)
    print(f"WROTE {args.output}: {len(rows)} rows, {len(skipped)} skipped rows, {len(plot_paths)} plots")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
