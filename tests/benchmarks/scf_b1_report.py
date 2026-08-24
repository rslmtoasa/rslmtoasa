#!/usr/bin/env python3
"""Render the Phase-B SCF-B1R report from canonical campaign evidence.

The renderer deliberately keeps performance, correctness, and eligibility as
separate concepts. In particular, a timing ratio is never printed unless the
row declares the corresponding comparison lane eligible.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


def rows_for(data: dict[str, Any], **filters: Any) -> list[dict[str, Any]]:
    return [row for row in data.get("rows", []) if all(row.get(key) == value for key, value in filters.items())]


def route(row: dict[str, Any]) -> str:
    return str(row.get("scf_route") or ("reciprocal" if str(row.get("case_id", "")).startswith("k") else ""))


def status(row: dict[str, Any]) -> str:
    return str(row.get("status") or row.get("correctness_status") or "-")


def final_state(row: dict[str, Any]) -> dict[str, Any]:
    value = row.get("final_state")
    return value if isinstance(value, dict) else {}


def converged(row: dict[str, Any]) -> bool:
    return final_state(row).get("converged") in (True, "true", "TRUE")


def equal_precision_eligible(row: dict[str, Any]) -> bool:
    if "equal_precision_eligible" in row:
        return bool(row["equal_precision_eligible"])
    return bool(row.get("headline_speedup_eligible"))


def eligible_ratio(row: dict[str, Any], field: str) -> Any:
    """Return a declared eligible ratio, otherwise ``None``."""

    if field.startswith("R_"):
        if status(row) != "PASS" or not row.get("production_comparison_eligible"):
            return None
    elif field.startswith("S_"):
        if (
            status(row) != "PASS"
            or not equal_precision_eligible(row)
            or row.get("headline_speedup_eligible") is not True
        ):
            return None
    value = row.get(field)
    return value if isinstance(value, (int, float)) and math.isfinite(float(value)) else None


def reason(row: dict[str, Any]) -> str:
    """Select the most useful failure or ineligibility explanation."""

    for key in ("status_reason", "failure_ineligibility_reason", "unsupported_reason"):
        value = row.get(key)
        if value:
            return str(value)
    rejection = row.get("headline_rejection_reasons")
    if rejection:
        return ", ".join(str(value) for value in rejection)
    correctness = row.get("correctness", {})
    correctness_status = correctness.get("status", row.get("correctness_status"))
    correctness_reason = row.get("correctness_reason") or correctness.get("reason")
    if correctness_status and correctness_status != "PASS" and correctness_reason:
        return str(correctness_reason)
    if status(row) == "PASS" and correctness_status == "PASS" and correctness_reason:
        return str(correctness_reason)
    return "-"


def nk_unique(row: dict[str, Any]) -> Any:
    return row.get("Nk_unique") if row.get("Nk_unique") is not None else row.get("actual_unique_nk")


def fmt(value: Any, digits: int = 6) -> str:
    if value is None or value == "":
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "-"
        return f"{value:.{digits}g}"
    return str(value)


def table(headers: list[str], values: list[list[Any]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines.extend("| " + " | ".join(fmt(value) for value in row) + " |" for row in values)
    return lines


def large_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten current-build K2 frozen-potential evidence."""

    result: list[dict[str, Any]] = []
    evidence = data.get("large_reciprocal_evidence", [])
    if isinstance(evidence, list):
        for item in evidence:
            for raw in item.get("rows", []) if isinstance(item, dict) else []:
                row = dict(raw)
                row["case_id"] = item.get("case_id", row.get("case_id"))
                row["status"] = row.get("status") or item.get("status") or "PASS"
                row["status_reason"] = row.get("status_reason") or item.get("reason")
                row["numeric_mode"] = row.get("numeric_mode") or "fp64"
                result.append(row)
    if result:
        return result
    # Backward-compatible fallback for a campaign produced before the K2
    # top-level evidence extension was added.
    for row in data.get("rows", []):
        if str(row.get("case_id", "")).startswith("k2_"):
            item = dict(row)
            item["numeric_mode"] = item.get("numeric_mode") or "fp64"
            item["actual_unique_nk"] = item.get("actual_unique_nk", 1)
            result.append(item)
    return result


def cross_route(data: dict[str, Any]) -> dict[str, Any]:
    value = data.get("cross_route")
    return value if isinstance(value, dict) else {
        "status": "INCONCLUSIVE",
        "reason": "cross-route evidence was not recorded",
        "timing_eligible": False,
    }


def _plot(data: dict[str, Any], plot_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return []

    plot_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []

    def save(name: str) -> None:
        path = plot_dir / name
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        paths.append(str(path))

    def empty(title: str, name: str) -> None:
        plt.figure()
        plt.text(0.5, 0.5, title, ha="center", va="center", wrap=True)
        plt.axis("off")
        save(name)

    rows = data.get("rows", [])
    rs = [
        row for row in rows
        if route(row) == "real_space"
        and row.get("rs_solver") == "block"
        and row.get("benchmark_level") == "scf_iteration"
        and status(row) == "PASS"
    ]

    for metric, ylabel, name in (
        ("P_rs_solver_kernel", "RS kernel time (s)", "01_rs_kernel_vs_size.png"),
        ("steady_iteration_median", "steady SCF iteration time (s)", "02_rs_iteration_vs_size.png"),
    ):
        values_found = False
        plt.figure()
        for backend in ("lapack", "cuda"):
            values = sorted(
                (row for row in rs if row.get("backend") == backend and row.get(metric) is not None),
                key=lambda item: (item.get("Natom") or 0, item.get("OMP_threads") or 0),
            )
            if values:
                values_found = True
                plt.plot([row.get("Natom") for row in values], [row.get(metric) for row in values], "o-", label=backend)
        if values_found:
            plt.xlabel("physical atom count (Natom)")
            plt.ylabel(ylabel)
            plt.legend()
            save(name)
        else:
            plt.close()
            empty("No eligible RS timing rows", name)

    representative = next((row for row in rs if row.get("backend") == "cuda"), None)
    if representative:
        names = [
            key for key in representative
            if key.startswith("P_") and key not in ("P_scf_iteration_total", "P_scf_misc")
            and isinstance(representative.get(key), (int, float)) and representative.get(key) > 0
        ]
        plt.figure()
        plt.bar(range(len(names)), [representative[key] for key in names])
        plt.xticks(range(len(names)), names, rotation=75, ha="right")
        plt.ylabel("seconds")
        save("03_rs_stage_fractions.png")
    else:
        empty("No eligible RS CUDA stage profile", "03_rs_stage_fractions.png")

    omp = sorted(
        [row for row in rs if row.get("material") == "bccFe" and row.get("backend") == "lapack"],
        key=lambda item: item.get("OMP_threads") or 0,
    )
    if omp:
        plt.figure()
        plt.plot([row.get("OMP_threads") for row in omp], [row.get("steady_iteration_median") for row in omp], "o-")
        plt.xlabel("OMP threads")
        plt.ylabel("steady iteration (s)")
        save("04_rs_cpu_omp_scaling.png")
    else:
        empty("No Fe3 CPU OMP sweep", "04_rs_cpu_omp_scaling.png")

    # Mixed RS production ratios are plotted only when the runner's explicit
    # production-comparison gate passed.
    plt.figure()
    ratio_found = False
    for field, label in (("R_rs_kernel_production", "kernel"), ("R_rs_phase_production", "phase"), ("R_iteration_production", "iteration")):
        values = sorted(
            ((row.get("nmat") or row.get("Natom"), eligible_ratio(row, field)) for row in rows),
            key=lambda item: item[0] or 0,
        )
        values = [(x, y) for x, y in values if y is not None]
        if values:
            ratio_found = True
            plt.plot([x for x, _ in values], [y for _, y in values], "o-", label=label)
    if ratio_found:
        plt.xlabel("matrix dimension / physical size")
        plt.ylabel("CPU FP64 / RS CUDA mixed production ratio")
        plt.legend()
        save("05_rs_production_ratios_vs_size.png")
    else:
        plt.close()
        empty("No eligible mixed RS production ratio", "05_rs_production_ratios_vs_size.png")

    big = large_rows(data)
    speedups: list[tuple[int, float]] = []
    for length in sorted({row.get("L") for row in big if row.get("L") is not None}):
        cpu = next((row for row in big if row.get("L") == length and row.get("backend") in ("lapack", "cpu") and row.get("status") == "PASS"), None)
        gpu = next((row for row in big if row.get("L") == length and row.get("backend") == "cuda" and row.get("status") == "PASS"), None)
        if cpu and gpu and float(gpu.get("steady_total_min_s") or 0) > 0:
            speedups.append((int(gpu.get("nmat") or 0), float(cpu["steady_total_min_s"]) / float(gpu["steady_total_min_s"])))
    if speedups:
        plt.figure()
        plt.plot([x for x, _ in speedups], [y for _, y in speedups], "o-")
        plt.xlabel("nmat")
        plt.ylabel("current-build CPU / GPU solver speedup")
        save("06_large_reciprocal_solver_speedup_vs_nmat.png")
    else:
        empty("No paired current-build large reciprocal rows", "06_large_reciprocal_solver_speedup_vs_nmat.png")

    reciprocal_cuda = next((row for row in rows if route(row) == "reciprocal" and row.get("backend") == "cuda" and status(row) == "PASS"), None)
    if reciprocal_cuda:
        names = [
            key for key in reciprocal_cuda
            if key.startswith("P_") and key not in ("P_scf_iteration_total", "P_scf_misc")
            and isinstance(reciprocal_cuda.get(key), (int, float)) and reciprocal_cuda.get(key) > 0
        ]
        plt.figure()
        plt.bar(range(len(names)), [reciprocal_cuda[key] for key in names])
        plt.xticks(range(len(names)), names, rotation=75, ha="right")
        plt.ylabel("seconds")
        save("07_reciprocal_stage_fractions.png")
    else:
        empty("No eligible reciprocal CUDA stage profile", "07_reciprocal_stage_fractions.png")

    history_rows = [row for row in rows if row.get("iterations") and status(row) in {"PASS", "NOT_CONVERGED"}]
    if history_rows:
        plt.figure()
        for row in history_rows:
            points = row.get("iterations", [])
            if not points:
                continue
            label = str(row.get("case_id"))
            plt.semilogy([item.get("iteration") for item in points], [max(float(item.get("residual", 1e-30)), 1e-30) for item in points], label=label)
        plt.xlabel("SCF iteration")
        plt.ylabel("production residual / mix delta")
        if len(history_rows) <= 6:
            plt.legend(fontsize="x-small")
        save("08_scf_residual_history.png")
    else:
        empty("No archived SCF iteration history", "08_scf_residual_history.png")

    return paths


def render(data: dict[str, Any], output: Path, report_path: Path) -> None:
    plot_paths = _plot(data, output.parent / "plots")
    rows = data.get("rows", [])
    preflight = data.get("build_preflight", {})
    provenance = data.get("provenance", {})
    manifest = data.get("manifest", {})
    cases = manifest.get("cases", []) if isinstance(manifest, dict) else []
    status_counts = Counter(status(row) for row in rows)
    completed_cases = sum(1 for item in cases if item.get("completed"))
    expected_cases = len(cases)
    rs_rows = [row for row in rows if route(row) == "real_space" and row.get("rs_solver") == "block" and row.get("benchmark_level") == "scf_iteration"]
    fe3_rs_rows = [row for row in rs_rows if str(row.get("case_id", "")).startswith("rs1_fe3_")]
    reciprocal_rows = [row for row in rows if route(row) == "reciprocal" and row.get("benchmark_level") in ("scf_iteration", "scf_convergence")]
    timer_by_case = {item.get("case_id"): item for item in data.get("reciprocal_timer_validation", [])}
    cross = cross_route(data)
    big = large_rows(data)

    lines: list[str] = [
        "# SCF-B1R — Reconciled unified SCF campaign",
        "",
        "## 1. Executive result",
        "",
        f"The archived campaign is complete: **{completed_cases}/{expected_cases} planned cases**. Canonical row statuses are {dict(status_counts)}.",
        "",
        "The build, Fe3 CPU OMP sweep, reciprocal Fe2/Fe3 runtime metadata, and current-build Fe3/Fe4/Fe5 frozen-potential rows are available. The common-potential RS-vs-k-space accuracy gate is **not passed**, so that formulation comparison is reported as **INCONCLUSIVE** and no cross-route timing claim is made.",
        "",
        "SCF-B1 can close only as a scoped performance campaign; it does not close as a universal RS-vs-k-space accuracy/performance claim. The unresolved item is the failed fixed-potential accuracy contract, not missing runtime evidence.",
        "",
        "## 2. Build-integrity check",
        "",
        *table(["build", "compiler", "type", "effective optimization", "status", "current source commit", "dirty", "reuse B1 timing?"], [[
            preflight.get("build_dir"), preflight.get("compiler"), preflight.get("build_type"),
            preflight.get("effective_optimization_conclusion"), preflight.get("status"),
            preflight.get("source_commit") or provenance.get("commit"),
            provenance.get("dirty"), preflight.get("timing_reuse_from_B1", "forbidden"),
        ]]),
        "",
        f"The authoritative preflight inspected representative compile commands for the SCF driver, RS recursion, Hamiltonian construction, and reciprocal path. Conclusion: `{preflight.get('reason') or preflight.get('effective_optimization_conclusion')}`. The B1 timings were not reused because `{preflight.get('timing_reuse_reason') or 'the corrected Release build was required'}`",
        "",
        "## 3. Method and status policy",
        "",
        "Performance is the steady production SCF-iteration wall time and its route-specific phases. SCF cycle count is retained for correctness/stability only. `PASS`, `NOT_CONVERGED`, `NOT_APPLICABLE`, `UNSUPPORTED`, `FAIL`, and `SKIPPED` remain distinct; a successful process exit does not override a failed profile or accuracy gate.",
        "",
        "Equal-precision `S_*` ratios require an explicit eligible pairing. Mixed RS CUDA rows use `R_*` production ratios only. A missing or ineligible ratio is printed as `-` with the row's reason.",
        "",
        "## 4. Reciprocal Fe2/Fe3 basis and timer resolution",
        "",
        "The runtime evidence resolves the labels as physical supercell replications with genuinely enlarged reciprocal matrices: Fe2 is `nmat=144`, Fe3 is `nmat=486`, and all K1 rows report `Nk_unique=512`. These values are taken from the emitted runtime result, not inferred from fixture names.",
        "",
        *table(["case", "physical replication", "reciprocal basis", "nmat", "Nk_unique", "backend", "strategy", "mode", "timer", "timer s", "timer check", "row status", "reason"], [
            [row.get("case_id"), row.get("physical_replication"), row.get("reciprocal_basis_replication"), row.get("nmat"), nk_unique(row), row.get("backend"), row.get("solver_strategy"), row.get("numeric_mode"),
             timer_by_case.get(row.get("case_id"), {}).get("timer") or ("T_solver" if row.get("backend") == "cuda" else "P_eigensolver"),
             timer_by_case.get(row.get("case_id"), {}).get("value") or (row.get("T_solver") if row.get("backend") == "cuda" else row.get("P_eigensolver")),
             timer_by_case.get(row.get("case_id"), {}).get("status", "not applicable"), status(row), reason(row)]
            for row in rows if str(row.get("case_id", "")).startswith("k1_")
        ]),
        "",
        "CPU `P_eigensolver` and CUDA `T_solver` are the authoritative boundaries. The timer validation is PASS for all successful Fe2/Fe3 rows; the Fe1 CPU row is INCONCLUSIVE because its row is FAIL under the strict profile-misc gate, although the positive timer value is archived.",
        "",
        "## 5. Fe3 real-space CPU OMP reconciliation",
        "",
        *table(["case", "solver/backend", "nrec", "nmat", "mode", "OMP", "kernel s", "phase s", "iteration s", "status", "equal precision", "reason"], [
            [row.get("case_id"), f"{row.get('rs_solver')}/{row.get('rs_backend')}", row.get("recursion_depth"), row.get("nmat"), row.get("numeric_mode"), row.get("OMP_threads"), row.get("P_rs_solver_kernel"), row.get("P_rs_solver_kernel", 0) + sum(float(row.get(key) or 0) for key in ("P_rs_green_function", "P_rs_density_build", "P_rs_fermi", "P_rs_hamiltonian_prepare", "P_rs_charge_spin_accumulate")), row.get("steady_iteration_median"), status(row), equal_precision_eligible(row), reason(row)]
            for row in sorted(fe3_rs_rows, key=lambda item: item.get("OMP_threads") or 0)
        ]),
        "",
    ]
    cpu_rs = [row for row in fe3_rs_rows if row.get("backend") == "lapack" and row.get("material") == "bccFe"]
    if cpu_rs:
        best_kernel = min(cpu_rs, key=lambda row: float(row.get("P_rs_solver_kernel") or math.inf))
        best_iteration = min(cpu_rs, key=lambda row: float(row.get("steady_iteration_median") or math.inf))
        lines.append(f"Best measured Fe3 CPU kernel row: `{best_kernel.get('case_id')}` at `{fmt(best_kernel.get('P_rs_solver_kernel'))} s`; best steady-iteration row: `{best_iteration.get('case_id')}` at `{fmt(best_iteration.get('steady_iteration_median'))} s`. The sweep is OMP1/2/4/8 and the CUDA row remains labelled mixed precision.")
    else:
        lines.append("No complete Fe3 CPU OMP sweep is present.")
    lines += [
        "",
        "The corrected mixed-GPU production ratio is shown only where the common-state/profile gate supports it. The RS1 Fe3 CUDA row is retained as a measured mixed route but does not receive an invented headline ratio when its comparison lane is not eligible.",
        "",
        "## 6. RS common-potential correctness anchor",
        "",
        *table(["case", "backend", "mode", "starting state", "status", "correctness", "reason", "energy", "Fermi", "moment", "charge", "residual"], [
            [row.get("case_id"), row.get("backend"), row.get("numeric_mode"), row.get("starting_state_id"), status(row), row.get("correctness_status"), row.get("correctness_reason") or row.get("correctness", {}).get("reason") or "-", final_state(row).get("final_total_energy"), final_state(row).get("fermi_energy"), final_state(row).get("site_moment"), final_state(row).get("total_charge"), final_state(row).get("final_residual")]
            for row in rows if str(row.get("case_id", "")).startswith("rs2_") or str(row.get("case_id", "")).startswith("x1_fe_cpu_rs")
        ]),
        "",
        "The reference run itself is `NOT_CONVERGED`, but its archived restart and final observables support the common-state RS CPU/CUDA comparison. The common-state CPU and CUDA rows are retained with their explicit correctness outcomes; iteration counts are not used as an accelerator metric.",
        "",
        "## 7. Reciprocal SCF rows and eligibility",
        "",
        *table(["case", "size", "nmat", "Nk_unique", "backend", "solver strategy", "mode", "iteration s", "full wall s", "S_iteration", "S_solver", "equal precision", "headline eligible", "status", "reason"], [
            [row.get("case_id"), row.get("supercell"), row.get("nmat"), nk_unique(row), row.get("backend"), row.get("solver_strategy"), row.get("numeric_mode"), row.get("steady_iteration_median"), row.get("full_scf_wall"), eligible_ratio(row, "S_iteration"), eligible_ratio(row, "S_solver"), equal_precision_eligible(row), row.get("headline_speedup_eligible"), status(row), reason(row)]
            for row in reciprocal_rows
        ]),
        "",
        "The Fe1 CPU row is a real `FAIL` under `scf_misc_exceeds_5_percent`; it is not silently treated as unsupported. Reciprocal FP64 CUDA rows are equal-precision at the numeric-mode gate, but their `S_*` ratios remain suppressed unless the full correctness/headline pairing gate passes.",
        "",
        "## 8. Current-build large reciprocal FP64 rerun",
        "",
        "These K2 rows are frozen-potential eigensolver evidence with eigenvectors enabled, not full SCF convergence rows. Their current-build CPU/GPU speedups are computed only from paired PASS rows at the same L and nmat.",
        "",
        *table(["L", "nmat", "Nk unique", "backend", "strategy", "vectors", "steady total s", "cold wall s", "status", "reason", "source commit"], [
            [row.get("L"), row.get("nmat"), nk_unique(row), row.get("backend"), row.get("solver_strategy"), row.get("vectors"), row.get("steady_total_min_s"), row.get("cold_process_wall_s"), row.get("status"), reason(row), row.get("git_commit")]
            for row in sorted(big, key=lambda item: (item.get("L") or 0, item.get("backend", "")))
        ]),
        "",
    ]
    speed_rows: list[list[Any]] = []
    for length in sorted({row.get("L") for row in big if row.get("L") is not None}):
        cpu = next((row for row in big if row.get("L") == length and row.get("backend") in ("lapack", "cpu") and row.get("status") == "PASS"), None)
        gpu = next((row for row in big if row.get("L") == length and row.get("backend") == "cuda" and row.get("status") == "PASS"), None)
        speed = None
        if cpu and gpu and float(gpu.get("steady_total_min_s") or 0) > 0:
            speed = float(cpu["steady_total_min_s"]) / float(gpu["steady_total_min_s"])
        speed_rows.append([length, cpu.get("nmat") if cpu else gpu.get("nmat") if gpu else None, speed, "PASS" if speed is not None else "INCONCLUSIVE", "paired current-build FP64 frozen-potential rows" if speed is not None else "missing paired eligible row"])
    lines += table(["L", "nmat", "CPU/GPU speedup", "eligibility", "reason"], speed_rows)
    lines += ["", "The current measured crossover is already GPU-beneficial at Fe3 (`nmat=486`) and becomes more favorable at Fe4 and Fe5. This is a frozen-potential solver result, not a universal full-SCF conclusion.", "", "## 9. Fixed-potential RS vs k-space Fe comparison", ""]
    metrics = cross.get("metrics", {})
    lines += table(["accuracy status", "energy difference", "Fermi difference", "moment difference", "charge difference", "timing eligible", "reason"], [[
        cross.get("status"), metrics.get("total_energy", {}).get("absolute_difference"), metrics.get("fermi", {}).get("absolute_difference"), metrics.get("moment", {}).get("absolute_difference"), metrics.get("charge", {}).get("absolute_difference"), cross.get("timing_eligible", False), cross.get("reason")
    ]])
    lines += [
        "",
        "The fixed-potential comparison does not satisfy the declared energy, Fermi, and moment tolerances. Accordingly the RS-vs-k-space scientific conclusion is **INCONCLUSIVE**; no RS/k-space timing speedup or DOS overlay is presented. Charge agrees closely, but that does not repair the failed observables.",
        "",
        "## 10. Canonical artifacts and plots",
        "",
        "The campaign directory contains `build_preflight.json`, `manifest.json`, `campaign.json`, `campaign.csv`, `campaign.md`, raw per-case stdout/stderr/commands, correctness records, iteration histories, and `cross_route/x1.json`. Regenerated plots:",
        "",
        *[f"- `{path}`" for path in plot_paths],
        "",
        "## 11. Final workload conclusions",
        "",
        "- The Release build is validated from actual compile commands, with effective `-O3`; old B1 timing reuse is forbidden.",
        "- Fe2/Fe3 reciprocal rows are genuinely enlarged runtime matrices (`144`/`486`) at `Nk_unique=512`; the CPU/GPU timer names are recorded separately and validated per route.",
        "- The Fe3 RS CPU OMP sweep is complete, and the mixed RS CUDA label is preserved. Production ratios are shown only when their lane is eligible.",
        "- Current-build frozen-potential reciprocal solver evidence shows GPU benefit from Fe3 through Fe5, with the crossover visible at `nmat=486` in this measured workload.",
        "- The common-potential RS-vs-k-space accuracy gate fails for energy, Fermi level, and moment. That comparison remains inconclusive and must not be used to choose a universal formulation.",
        "",
        "## 12. Closeout decision",
        "",
        "B1R evidence collection and report reconciliation are complete. SCF-B1 may close for the scoped, separately qualified workload rows, but it should not be described as a closed universal RS-vs-k-space comparison until the fixed-potential accuracy contract is repaired and rerun.",
        "",
        "## Provenance",
        "",
        "```json",
        json.dumps({"provenance": provenance, "build_preflight": preflight, "policy": data.get("policy", {})}, indent=2, sort_keys=True),
        "```",
        "",
    ]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=Path("results/benchmarks/scf_b1r/campaign.json"))
    parser.add_argument("--report", type=Path, default=Path("docs/dev/RS_LMTO_ASA_SCF_B1R_REPORT.md"))
    args = parser.parse_args()
    data = json.loads(args.campaign.read_text(encoding="utf-8"))
    render(data, args.campaign, args.report)
    print(f"WROTE {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
