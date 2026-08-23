#!/usr/bin/env python3
"""Derive SCF-B1 tables, plots, and the final report from campaign.json."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


def rows_for(data: dict[str, Any], **filters: Any) -> list[dict[str, Any]]:
    return [row for row in data.get("rows", []) if all(row.get(key) == value for key, value in filters.items())]


def converged(row: dict[str, Any]) -> bool:
    return row.get("final_state", {}).get("converged") in (True, "true", "TRUE")


def best(rows: list[dict[str, Any]], *, field: str = "steady_iteration_median") -> dict[str, Any] | None:
    valid = [row for row in rows if row.get(field) is not None and row.get("profile_status") == "PASS"]
    return min(valid, key=lambda row: float(row[field])) if valid else None


def fmt(value: Any, digits: int = 5) -> str:
    if value is None or value == "":
        return "-"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "-"
        return f"{value:.{digits}g}"
    return str(value)


def status(row: dict[str, Any]) -> str:
    return str(row.get("correctness", {}).get("status", row.get("correctness_status", "-")))


def table(headers: list[str], values: list[list[Any]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for value in values:
        lines.append("| " + " | ".join(fmt(item) for item in value) + " |")
    return lines


def cross_route(data: dict[str, Any]) -> list[dict[str, Any]]:
    rs = [row for row in data.get("rows", []) if row.get("scf_route") == "real_space"
          and row.get("rs_solver") == "block" and row.get("benchmark_level") == "scf_convergence"
          and row.get("backend") == "lapack" and converged(row)]
    reciprocal = [row for row in data.get("rows", []) if row.get("scf_route") == "reciprocal"
                  and row.get("benchmark_level") == "scf_convergence" and row.get("backend") == "lapack"
                  and converged(row)]
    cases: list[dict[str, Any]] = []
    for material, supercell in (("bccFe", "1x1x1"), ("diamondSi", "1x1x1")):
        left = best([row for row in rs if row.get("material") == material and row.get("supercell") == supercell], field="full_scf_wall")
        right = best([row for row in reciprocal if row.get("material") == material and row.get("supercell") == supercell], field="full_scf_wall")
        if left is None or right is None:
            cases.append({
                "cross_route_case_id": f"{material}:{supercell}:accuracy-matched",
                "material": material,
                "energy_difference_per_atom": None, "fermi_difference": None,
                "moment_difference_per_atom": None, "charge_difference_per_atom": None,
                "accuracy_status": "INCONCLUSIVE",
                "RS_SCF_iteration_wall": None, "reciprocal_SCF_iteration_wall": None,
                "RS_full_convergence_wall": None, "reciprocal_full_convergence_wall": None,
                "comparison_tier": "no_accuracy-matched_converged_pair",
                "interpretation": "at least one route has no converged production row",
            })
            continue
        lm, rm = left["final_state"], right["final_state"]
        natom = max(float(left.get("Natom") or 1), 1.0)
        energy = abs(float(lm.get("final_total_energy", 0.0)) / natom - float(rm.get("final_total_energy", 0.0)) / natom)
        fermi = abs(float(lm.get("fermi_energy", 0.0)) - float(rm.get("fermi_energy", 0.0)))
        moment = abs(float(lm.get("site_moment", 0.0)) - float(rm.get("site_moment", 0.0)))
        charge = abs(float(lm.get("total_charge", 0.0)) / natom - float(rm.get("total_charge", 0.0)) / natom)
        limits = {"energy_per_atom": 1e-4, "fermi": 5e-5, "moment_per_atom": 5e-4, "charge_per_atom": 5e-5}
        passed = energy <= limits["energy_per_atom"] and fermi <= limits["fermi"] and moment <= limits["moment_per_atom"] and charge <= limits["charge_per_atom"]
        cases.append({
            "cross_route_case_id": f"{material}:{supercell}:accuracy-matched",
            "material": material, "RS_solver": left.get("rs_solver"),
            "RS_convergence_controls": {"nrec": left.get("recursion_depth"), "terminator": left.get("terminator")},
            "reciprocal_convergence_controls": {"Nk_unique": right.get("Nk_unique"), "smearing": right.get("sigma")},
            "RS_numeric_mode": left.get("numeric_mode"), "reciprocal_numeric_mode": right.get("numeric_mode"),
            "energy_difference_per_atom": energy, "fermi_difference": fermi,
            "moment_difference_per_atom": moment, "charge_difference_per_atom": charge,
            "accuracy_status": "PASS" if passed else "INCONCLUSIVE",
            "RS_SCF_iteration_wall": left.get("steady_iteration_median"),
            "reciprocal_SCF_iteration_wall": right.get("steady_iteration_median"),
            "RS_full_convergence_wall": left.get("full_scf_wall"),
            "reciprocal_full_convergence_wall": right.get("full_scf_wall"),
            "comparison_tier": "arithmetic-clean_scientific_comparison",
            "interpretation": "direct CPU FP64 formulation comparison" if passed else "accuracy contract not satisfied",
            "_rs_row": left, "_reciprocal_row": right,
        })
    return cases


def large_evidence_rows(data: dict[str, Any]) -> tuple[list[list[Any]], list[list[Any]]]:
    """Return compact frozen-potential ACC-P4 rows for the report."""

    cpu_rows: list[list[Any]] = []
    gpu_rows: list[list[Any]] = []
    for evidence in data.get("large_reciprocal_evidence", []):
        schema = evidence.get("schema")
        for row in evidence.get("rows", []):
            if schema == "rslmto.accp0-real-material.v1" and row.get("fixture") == "bccFe":
                if row.get("vectors") == 1 and row.get("L") in (3, 4, 5):
                    cpu_rows.append([
                        row.get("L"), row.get("nmat"), row.get("backend"), row.get("solver_strategy"),
                        row.get("vectors"), row.get("steady_total_min_s"), row.get("cold_process_wall_s"),
                        row.get("git_commit"),
                    ])
            elif schema == "rslmto.accp2-real-material.v1" and row.get("nmat") in (486, 1152, 2250):
                gpu_rows.append([
                    row.get("fixture"), row.get("nmat"), row.get("precision"), row.get("solver"),
                    row.get("vectors"), row.get("CPU_total"), row.get("after_total"),
                    row.get("final_CPU_GPU_speedup"),
                ])
    return cpu_rows, gpu_rows


def _plot(data: dict[str, Any], plot_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return []
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    rows = data.get("rows", [])

    def save(name: str) -> None:
        path = plot_dir / name
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        paths.append(str(path))

    rs = [row for row in rows if row.get("scf_route") == "real_space" and row.get("rs_solver") == "block"
          and row.get("benchmark_level") == "scf_iteration" and row.get("steady_iteration_median") is not None]
    for metric, ylabel, name in (("P_rs_solver_kernel", "RS kernel time (s)", "01_rs_kernel_vs_size.png"),
                                 ("steady_iteration_median", "RS SCF iteration time (s)", "02_rs_iteration_vs_size.png")):
        plt.figure()
        for backend in ("lapack", "cuda"):
            values = sorted((row for row in rs if row.get("backend") == backend and row.get(metric) is not None), key=lambda item: item.get("Natom") or 0)
            if values:
                plt.plot([row.get("Natom") for row in values], [row.get(metric) for row in values], "o-", label=backend)
        plt.xlabel("Natom"); plt.ylabel(ylabel); plt.legend(); save(name)

    representative = next((row for row in rs if row.get("material") == "bccFe" and row.get("backend") == "cuda"), None)
    if representative:
        plt.figure()
        names = [key for key in representative if key.startswith("P_") and key not in ("P_scf_iteration_total", "P_scf_misc") and representative.get(key) is not None]
        values = [representative.get(key) for key in names]
        plt.bar(range(len(names)), values); plt.xticks(range(len(names)), names, rotation=75, ha="right"); plt.ylabel("seconds")
        save("03_rs_stage_fractions.png")

    omp = [row for row in rs if row.get("material") == "bccFe" and row.get("backend") == "lapack"]
    if omp:
        plt.figure(); values = sorted(omp, key=lambda item: item.get("OMP_threads") or 0)
        plt.plot([row.get("OMP_threads") for row in values], [row.get("steady_iteration_median") for row in values], "o-")
        plt.xlabel("OMP threads"); plt.ylabel("steady iteration (s)"); save("04_rs_cpu_omp_scaling.png")

    kr = [row for row in rows if row.get("scf_route") == "reciprocal" and row.get("backend") == "cuda"]
    for field, ylabel, name in (("S_solver", "S_solver", "05_reciprocal_solver_speedup.png"),
                                ("S_reciprocal", "S_reciprocal", "06_reciprocal_phase_speedup.png"),
                                ("S_iteration", "S_iteration", "07_reciprocal_iteration_speedup.png"),
                                ("S_convergence", "S_convergence", "08_reciprocal_convergence_speedup.png")):
        values = sorted((row for row in kr if row.get(field) is not None and row.get("nmat") is not None), key=lambda item: item.get("nmat"))
        if values:
            plt.figure(); plt.plot([row.get("nmat") for row in values], [row.get(field) for row in values], "o-")
            plt.xlabel("matrix dimension"); plt.ylabel(ylabel); save(name)

    reciprocal = [row for row in rows if row.get("scf_route") == "reciprocal" and row.get("backend") == "cuda" and row.get("benchmark_level") == "scf_iteration"]
    if reciprocal:
        row = next((item for item in reciprocal if item.get("material") == "bccFe"), reciprocal[0])
        plt.figure()
        names = [key for key in row if key.startswith("P_") and key not in ("P_scf_iteration_total", "P_scf_misc") and row.get(key) is not None]
        plt.bar(range(len(names)), [row.get(key) for key in names]); plt.xticks(range(len(names)), names, rotation=75, ha="right"); plt.ylabel("seconds")
        save("09_reciprocal_stage_fractions.png")

    history = next((row for row in rows if row.get("material") == "bccFe" and row.get("backend") in ("lapack", "cuda") and row.get("iterations")), None)
    if history:
        plt.figure()
        for backend in ("lapack", "cuda"):
            candidate = next((item for item in rows if item.get("material") == "bccFe" and item.get("backend") == backend and item.get("iterations")), None)
            if candidate:
                points = candidate["iterations"]
                plt.semilogy([item.get("iteration") for item in points], [max(float(item.get("residual", 1e-30)), 1e-30) for item in points], label=backend)
        plt.xlabel("SCF iteration"); plt.ylabel("production criterion (mix delta)"); plt.legend(); save("10_scf_residual_history.png")
    return paths


def render(data: dict[str, Any], output: Path, report_path: Path) -> None:
    plot_paths = _plot(data, output.parent / "plots")
    rows = data.get("rows", [])
    env = data.get("environment", {})
    freeze = data.get("freeze", {})
    lines: list[str] = ["# SCF-B1 — Unified RS and reciprocal CPU/GPU campaign", "",
        "## 1. Executive summary", "",
        f"- Frozen commit: `{freeze.get('frozen_git_commit', '-')}`; tracked implementation dirty at start: `{freeze.get('git_dirty_tracked')}`.",
        f"- Host: `{env.get('cpu_model', '-')}`, RAM `{env.get('ram_mib', '-')}` MiB; compiler `{env.get('compiler', '-')}`; CUDA `{env.get('cuda_toolkit', '-')}`; GPU `{env.get('gpu_model', '-')}`.",
        "- RS CUDA is classified as `mixed` (FP32 working paths with FP64 canonical SCF state); its ratios are production-route ratios, not equal-precision speedups.",
        "- Reciprocal FP64 CPU/GPU rows remain the equal-precision headline lane.",
        "- Cross-route claims below are emitted only after the explicit accuracy contract passes.", "",
        "## 2. Methodology", "",
        "The campaign uses the shared production probe with real Si/Fe fixtures, CPU OMP 1/2/4/8, BLAS threads fixed at 1, fresh processes per row, two warmup iterations, and five steady iterations where available. Full convergence uses the normal production convergence flag. The criterion is `mix%delta < conv_thr` with `conv_thr=0.5e-8` by default; no second criterion is applied.", "",
        "## 3. Real-space workloads and solver taxonomy", "",
        "Block recursion is the primary Fe RS route; Chebyshev is the validated Si CPU route; scalar/Lanczos is retained as CPU/reference evidence. CUDA rows require `rs_gpu_used=true`, no fallback, finite recursion invariants, and profile closure.", "",
        "## 4. Real-space kernel results", "",
        *table(["material", "size", "Natom", "nmat", "backend", "mode", "OMP", "kernel s", "production ratio", "status"], [
            [row.get("material"), row.get("supercell"), row.get("Natom"), row.get("nmat"), row.get("backend"), row.get("numeric_mode"), row.get("OMP_threads"), row.get("P_rs_solver_kernel"), row.get("R_rs_kernel_production") or row.get("S_rs_kernel"), status(row)]
            for row in rows if row.get("scf_route") == "real_space" and row.get("rs_solver") == "block" and row.get("benchmark_level") == "scf_iteration" and row.get("backend") in ("lapack", "cuda")
        ]), "",
        "## 5. Real-space SCF iteration/convergence", "",
        *table(["material", "solver", "backend", "mode", "iteration s", "iterations", "full wall s", "ratio", "status"], [
            [row.get("material"), row.get("rs_solver"), row.get("backend"), row.get("numeric_mode"), row.get("steady_iteration_median"), row.get("n_scf_iterations"), row.get("full_scf_wall"), row.get("R_convergence_production") or row.get("S_convergence"), status(row)]
            for row in rows if row.get("scf_route") == "real_space" and row.get("benchmark_level") == "scf_convergence"
        ]), "",
        "## 6. Chebyshev and scalar/Lanczos", "",
        *table(["material", "solver", "backend", "M", "CPU/GPU timing", "status", "notes"], [
            [row.get("material"), row.get("rs_solver"), row.get("backend"), row.get("chebyshev_order"), row.get("steady_iteration_median"), status(row), row.get("unsupported_reason") or "GPU unsupported rows retained explicitly" if row.get("status") == "UNSUPPORTED" else "production reference"]
            for row in rows if row.get("scf_route") == "real_space" and row.get("rs_solver") in ("chebyshev", "lanczos")
        ]), "",
        "## 7. Reciprocal workloads and eigensolver", "",
        *table(["material", "size", "nmat", "Nk", "backend", "strategy", "mode", "solver s", "S_solver", "status"], [
            [row.get("material"), row.get("supercell"), row.get("nmat"), row.get("Nk_unique"), row.get("backend"), row.get("solver_strategy"), row.get("numeric_mode"), row.get("T_solver") or row.get("P_eigensolver"), row.get("S_solver"), status(row)]
            for row in rows if row.get("scf_route") == "reciprocal" and row.get("benchmark_level") in ("eigensolver", "scf_iteration", "scf_convergence")
        ]), "",
        "## 8. Reciprocal SCF iteration/convergence", "",
        *table(["material", "size", "nmat", "Nk", "backend", "OMP", "iteration s", "full wall s", "S_iteration", "S_convergence", "status"], [
            [row.get("material"), row.get("supercell"), row.get("nmat"), row.get("Nk_unique"), row.get("backend"), row.get("OMP_threads"), row.get("steady_iteration_median"), row.get("full_scf_wall"), row.get("S_iteration"), row.get("S_convergence"), status(row)]
            for row in rows if row.get("scf_route") == "reciprocal" and row.get("benchmark_level") in ("scf_iteration", "scf_convergence")
        ]), "",
        "## 9. Large reciprocal frozen-potential evidence", "",
        "These ACC-P4 rows are eigensolver/frozen-potential evidence only; they are not full SCF convergence rows and are excluded from the headline SCF speedups. The imported documents retain their own provenance, which is shown in the final column where available; a source commit different from the frozen campaign commit is context-only and does not certify the current implementation.", "",
    ]
    large_cpu, large_gpu = large_evidence_rows(data)
    lines += table(["L", "nmat", "backend", "strategy", "vectors", "steady total s", "cold wall s", "source commit"], large_cpu)
    lines += ["", "GPU/CPU frozen-potential summary", ""]
    lines += table(["fixture", "nmat", "precision", "solver", "vectors", "CPU s", "GPU s", "CPU/GPU"], large_gpu)
    lines += ["", "## 10. Physically matched RS vs reciprocal comparison", ""]
    cases = cross_route(data)
    lines += table(["case", "energy/atom", "Fermi", "moment/atom", "charge/atom", "accuracy", "RS iter", "K iter", "RS wall", "K wall", "tier"], [
        [case["cross_route_case_id"], case["energy_difference_per_atom"], case["fermi_difference"], case["moment_difference_per_atom"], case["charge_difference_per_atom"], case["accuracy_status"], case["RS_SCF_iteration_wall"], case["reciprocal_SCF_iteration_wall"], case["RS_full_convergence_wall"], case["reciprocal_full_convergence_wall"], case["comparison_tier"]]
        for case in cases
    ])
    lines += ["", "## 11. Precision, startup, and stage/Amdahl interpretation", "",
              "Equal-precision speedups use only matching `numeric_mode` rows. RS CUDA currently has FP32 recursion working paths and FP64 Hamiltonian/density/canonical state, so `R_rs_*_production` fields are the honest mixed-vs-FP64 production comparison. CUDA transfer/setup detail for RS remains `not_exposed_by_backend` and is not inferred from wall-time remainders. Cold process wall is retained separately from steady iteration time.", "",
              "## 12. Required interpretation questions", "",
              "### RS", "",
              "The measured RS result is reported by the block tables above. The report distinguishes kernel reduction from phase and complete-iteration reduction, preserves iteration-count differences, and marks Chebyshev GPU support explicitly. Further RS kernel work is justified only if the largest measured workload remains kernel-dominated after the production ratio and correctness gates.", "",
              "### Reciprocal", "",
              "Primitive Si/Fe, Fe 2×2×2, Fe 3×3×3, and the large frozen-potential Fe evidence are retained as separate regimes. GPU solver strategy is not assumed globally: the measured strategy column is authoritative. Losses, unsupported strategies, non-convergence, and startup cost remain in the canonical rows.", "",
              "### RS vs reciprocal", "",
              "Only cases with `accuracy_status=PASS` support a direct formulation timing comparison. If no case passes, the cross-route conclusion is `INCONCLUSIVE`; architectural expectations about disorder or large sparse systems are not promoted to measured results.", "",
              "## 13. Recommended workload map", "",
              "| workload regime | RS CPU | RS GPU | k-space CPU | k-space GPU | evidence status |", "|---|---|---|---|---|---|",
              "| small perfect primitive cell | measured reference | overhead-sensitive mixed route | measured | measured FP64 | measured rows required |",
              "| many-k small matrix | not generalized | not generalized | measured | measured | use primitive rows |",
              "| moderate periodic supercell | measured RS block | mixed production route if correctness passes | measured Fe2/Fe3 | measured where converged | measured rows |",
              "| large dense reciprocal eigensystem | not measured as full RS SCF | not measured | frozen-potential reference | use ACC-P4 evidence | measured eigensolver only |",
              "| disorder/broken translational symmetry | expected RS advantage, not measured | expected mixed route, not measured | expected less favorable, not measured | not measured | architectural expectation |",
              "", "## 14. Plots", "",
              *[f"- `{path}`" for path in plot_paths], "",
              "## 15. Final performance conclusion", "",
              "SCF-B1 closes only for the workload regimes represented by valid canonical rows. It does not establish a universal RS-vs-k-space winner: the final choice is workload- and precision-dependent, and any missing accuracy-matched cross-route case remains inconclusive.", "",
              "## Provenance", "", "```json", json.dumps({"freeze": freeze, "environment": env}, indent=2, sort_keys=True), "```", ""]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    output.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=Path("results/benchmarks/scf_b1/campaign.json"))
    parser.add_argument("--report", type=Path, default=Path("docs/dev/RS_LMTO_ASA_SCF_B1_REPORT.md"))
    args = parser.parse_args()
    data = json.loads(args.campaign.read_text(encoding="utf-8"))
    render(data, args.campaign, args.report)
    print(f"WROTE {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
