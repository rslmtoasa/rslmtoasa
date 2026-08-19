#!/usr/bin/env python3
"""Run the ACC-P0 persistent real-material benchmark campaign.

The script intentionally separates three concerns:

* production fixtures and repeated bcc-Fe supercells;
* persistent-process measurements recorded by ``benchmark_harness``;
* a CPU-only production oracle for primitive/supercell band folding.

It does not declare a GPU win.  It records the measurements needed to decide
whether the accelerator is useful for a stated workload and mesh density.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_harness import capture_environment, make_document, run_command  # noqa: E402


REQUIRED_COLUMNS = [
    "fixture",
    "source",
    "workload",
    "L",
    "natom",
    "nmat",
    "nominal_mesh",
    "actual_unique_nk",
    "tile",
    "eigenvectors",
    "backend",
    "solver_strategy",
    "support_status",
    "unsupported_reason",
    "cold_process_wall_s",
    "cuda_context_backend_init_s",
    "first_solve_s",
    "steady_solve_median_s",
    "metric_repetitions",
    "Hk_CPU_s",
    "H2D_s",
    "solver_s",
    "D2H_s",
    "total_steady_s",
    "memory_estimate_mib",
    "memory_free_before_mib",
    "memory_total_mib",
]


def parse_mesh(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace("x", " ").replace(",", " ").split()
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("mesh must be NxNxN")
    values = tuple(int(part) for part in parts)
    if any(value < 1 for value in values):
        raise argparse.ArgumentTypeError("mesh dimensions must be positive")
    return values  # type: ignore[return-value]


def parse_int_list(text: str) -> list[int]:
    values = [int(part) for part in text.split(",") if part]
    if not values or any(value < 1 for value in values):
        raise argparse.ArgumentTypeError("expected a comma-separated list of positive integers")
    return values


def fixture_input(source: Path, destination: Path) -> Path:
    shutil.copytree(source, destination, dirs_exist_ok=True)
    return destination / "input.nml"


def make_bcc_supercell(source: Path, destination: Path, length: int) -> Path:
    """Create an explicit L^3 bcc primitive supercell from production Fe data."""

    destination.mkdir(parents=True, exist_ok=True)
    potential = source / "Fe.nml"
    if not potential.is_file():
        raise FileNotFoundError(potential)
    natom = length**3
    labels = [f"Fe{i:04d}" for i in range(1, natom + 1)]
    for label in labels:
        shutil.copy2(potential, destination / f"{label}.nml")

    input_lines = [
        "&calculation",
        "  pre_processing = 'bravais'",
        "  verbose = F",
        "/",
        "&lattice",
        "  ndim = 60000",
        "  rc = 80",
        "  alat = 2.86120",
        "  crystal_sym = 'file'",
        "  wav = 1.40880",
        "  ntype = %d" % natom,
        "  ct(:) = 3.0d0",
        "  r2 = 9.0d0",
        # The canonical bccFe input currently uses the production legacy
        # structure-constant backend.  Keep the primitive and supercell
        # paths identical so the folding oracle tests physics, not backend
        # implementation differences.
        "  strux_backend = 'legacy'",
        "/",
        "&atoms",
        "  database = './'",
    ]
    input_lines.extend(f"  label({index}) = '{label}'" for index, label in enumerate(labels, start=1))
    input_lines.extend(
        [
            "/",
            "&self",
            "  nstep = 1",
            "/",
            "&energy",
            "  fermi = -0.042265",
            "  energy_min = -2.0",
            "  energy_max = 0.8",
            "/",
            "&control",
            "  calctype = 'B'",
            "  nsp = 2",
            "  recur = 'block'",
            "  lmax = 2",
            "/",
            "&mix",
            "  beta = 0.01",
            "  mixtype = 'broyden'",
            "/",
            "&hamiltonian",
            "  hoh = .false.",
            "/",
        ]
    )
    (destination / "input.nml").write_text("\n".join(input_lines) + "\n", encoding="utf-8")

    # bcc primitive vectors from the production bcc lattice path, repeated
    # along each primitive direction.  Fractional coordinates are with
    # respect to the resulting supercell vectors.
    vectors = [
        (-0.5 * length, 0.5 * length, 0.5 * length),
        (0.5 * length, -0.5 * length, 0.5 * length),
        (0.5 * length, 0.5 * length, -0.5 * length),
    ]
    lattice_lines = [
        "&lattice",
        f"  nbulk_bulk = {natom}",
        f"  ntot = {natom}",
        f"  nbas = {natom}",
        f"  nrec = {natom}",
    ]
    for column, vector in enumerate(vectors, start=1):
        lattice_lines.append("  a(:, %d) = %0.16f, %0.16f, %0.16f" % (column, *vector))
    base_vectors = [
        (-0.5, 0.5, 0.5),
        (0.5, -0.5, 0.5),
        (0.5, 0.5, -0.5),
    ]
    index = 0
    for i in range(length):
        for j in range(length):
            for k in range(length):
                index += 1
                position = tuple(
                    i * base_vectors[0][component]
                    + j * base_vectors[1][component]
                    + k * base_vectors[2][component]
                    for component in range(3)
                )
                lattice_lines.append(
                    "  crd(:, %d) = %0.16f, %0.16f, %0.16f" % (index, *position)
                )
    for key in ("izp", "no", "iu", "ib", "irec"):
        lattice_lines.extend(f"  {key}({index}) = {index}" for index in range(1, natom + 1))
    lattice_lines.append("/")
    (destination / "lattice.nml").write_text("\n".join(lattice_lines) + "\n", encoding="utf-8")
    return destination / "input.nml"


def run_oracle(
    binary: Path,
    fixture_dir: Path,
    output: Path,
    mesh: tuple[int, int, int],
    length: int,
    backend: str = "lapack",
) -> list[float]:
    dump = output.with_suffix(".eig")
    command = [
        str(binary),
        "--mode",
        "oracle",
        "--backend",
        backend,
        "--fixture",
        "bccFe",
        "--input",
        "input.nml",
        "--mesh",
        "%dx%dx%d" % mesh,
        "--dump-eigenvalues",
        str(dump),
        "--L",
        str(length),
    ]
    completed = subprocess.run(command, cwd=fixture_dir, capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    output.write_text(text, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"{backend} eigenvalue oracle failed ({completed.returncode}):\n{text}")
    return [float(line) for line in dump.read_text(encoding="utf-8").split()]


def degeneracy_signature(values: list[float], tolerance: float = 1.0e-8) -> list[int]:
    if not values:
        return []
    ordered = sorted(values)
    groups: list[int] = []
    current = 1
    for previous, value in zip(ordered, ordered[1:]):
        if abs(value - previous) <= tolerance:
            current += 1
        else:
            groups.append(current)
            current = 1
    groups.append(current)
    return groups


def compare_oracle(
    primitive: list[float], supercell: list[float], *, require_degeneracy: bool = True
) -> dict[str, Any]:
    if len(primitive) != len(supercell):
        return {
            "passed": False,
            "primitive_eigenvalues": len(primitive),
            "supercell_eigenvalues": len(supercell),
            "max_abs_error": None,
            "degeneracy_match": False,
        }
    left = sorted(primitive)
    right = sorted(supercell)
    error = max((abs(a - b) for a, b in zip(left, right)), default=0.0)
    left_degeneracy = degeneracy_signature(left)
    right_degeneracy = degeneracy_signature(right)
    return {
        "passed": error <= 1.0e-7 and (not require_degeneracy or left_degeneracy == right_degeneracy),
        "primitive_eigenvalues": len(primitive),
        "supercell_eigenvalues": len(supercell),
        "max_abs_error": error,
        "degeneracy_match": left_degeneracy == right_degeneracy,
        "degeneracy_required": require_degeneracy,
        "primitive_degeneracy": left_degeneracy,
        "supercell_degeneracy": right_degeneracy,
    }


def run_measurement(
    *,
    binary: Path,
    build_dir: Path,
    fixture_dir: Path,
    output: Path,
    fixture: str,
    source: str,
    workload: str,
    length: int,
    mesh: tuple[int, int, int],
    backend: str,
    solver_strategy: str,
    tile: int,
    vectors: bool,
    warmups: int,
    repetitions: int,
) -> dict[str, Any]:
    command = [
        str(binary),
        "--mode",
        "persistent",
        "--backend",
        backend,
        "--solver-strategy",
        solver_strategy,
        "--fixture",
        fixture,
        "--source",
        source,
        "--workload",
        workload,
        "--input",
        "input.nml",
        "--mesh",
        "%dx%dx%d" % mesh,
        "--tile-size",
        str(tile),
        "--eigenvectors",
        "1" if vectors else "0",
        "--warmups",
        str(warmups),
        "--repetitions",
        str(repetitions),
        "--L",
        str(length),
    ]
    wall_times, profiles, last_output = run_command(
        command,
        cwd=fixture_dir,
        warmups=0,
        repetitions=1,
        persistent=True,
    )
    if len(profiles) != 1 or not profiles[0]:
        raise RuntimeError(f"ACC-P0 command emitted no parseable profile:\n{last_output}")
    for profile in profiles[0].values():
        profile.setdefault("metrics", {})["cold_process_wall_s"] = wall_times[0]
    metadata = capture_environment(ROOT, build_dir, omp_threads=None, mpi_ranks=1)
    metadata.update({"fixture": fixture, "source": source, "workload": workload, "L": length})
    name = f"accp1_{workload}_{fixture}_L{length}_{mesh[0]}x{mesh[1]}x{mesh[2]}_{backend}_{solver_strategy}_tile{tile}_v{int(vectors)}"
    document = make_document(
        name=name,
        benchmark_class="component",
        labels=["performance", "reciprocal", "accp1", workload, backend, solver_strategy],
        command=command,
        metadata=metadata,
        wall_times=wall_times,
        profile_samples=profiles,
        last_output=last_output,
        warmups=warmups,
        persistent=True,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    record = document["profile_records"][0]
    row = dict(record.get("metadata", {}))
    row.update(record.get("samples", [{}])[0])
    row.update({"name": name, "backend": backend, "solver_strategy": solver_strategy,
                "vectors": int(vectors), "nominal_mesh": "x".join(map(str, mesh))})
    row.setdefault("support_status", "supported")
    row.setdefault("unsupported_reason", "")
    return row


def preflight(binary: Path, fixture_dir: Path, fixture: str, length: int, mesh: tuple[int, int, int], backend: str, solver_strategy: str, tile: int, vectors: bool) -> dict[str, float]:
    command = [
        str(binary),
        "--mode",
        "preflight",
        "--backend",
        backend,
        "--solver-strategy",
        solver_strategy,
        "--fixture",
        fixture,
        "--input",
        "input.nml",
        "--mesh",
        "%dx%dx%d" % mesh,
        "--tile-size",
        str(tile),
        "--eigenvectors",
        "1" if vectors else "0",
        "--L",
        str(length),
    ]
    completed = subprocess.run(command, cwd=fixture_dir, capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    if completed.returncode:
        raise RuntimeError(f"ACC-P0 preflight failed ({completed.returncode}):\n{text}")
    line = next((line for line in text.splitlines() if line.startswith("ACCP0_PREFLIGHT ")), None)
    if line is None:
        raise RuntimeError(f"ACC-P0 preflight emitted no record:\n{text}")
    values: dict[str, float] = {}
    for token in line.split()[1:]:
        key, value = token.split("=", 1)
        try:
            values[key] = float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
    values["output"] = text  # type: ignore[assignment]
    return values


def run_validation(
    *, binary: Path, fixture_dir: Path, fixture: str, length: int,
    mesh: tuple[int, int, int], solver_strategy: str, tile: int,
) -> dict[str, Any]:
    command = [
        str(binary), "--mode", "validate", "--backend", "cuda",
        "--solver-strategy", solver_strategy, "--fixture", fixture,
        "--input", "input.nml", "--mesh", "%dx%dx%d" % mesh,
        "--tile-size", str(tile), "--L", str(length),
    ]
    completed = subprocess.run(command, cwd=fixture_dir, capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    if completed.returncode:
        raise RuntimeError(f"ACC-P1 validation failed ({completed.returncode}):\n{text}")
    line = next((line for line in text.splitlines() if line.startswith("ACCP1_VALIDATION ")), None)
    if line is None:
        raise RuntimeError(f"ACC-P1 validation emitted no record:\n{text}")
    values: dict[str, Any] = {}
    for token in line.split()[1:]:
        key, value = token.split("=", 1)
        try:
            values[key] = float(value.replace("D", "E").replace("d", "e"))
        except ValueError:
            values[key] = value
    values.update({"fixture": fixture, "L": length, "solver_strategy": solver_strategy, "mesh": "x".join(map(str, mesh))})
    return values


def matched_mesh(base: int, length: int) -> tuple[int, int, int]:
    value = max(1, int(round(base / length)))
    return (value, value, value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/benchmarks/accp0")
    parser.add_argument("--quick", action="store_true", help="small validation campaign; still exercises both workloads")
    parser.add_argument("--skip-cuda", action="store_true")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--tiles", type=parse_int_list, default=[1, 8, 16])
    parser.add_argument("--meshes", type=parse_int_list, default=[1, 2, 4, 8])
    parser.add_argument("--vectors", action="store_true", help="include eigenvector requests in addition to values-only")
    args = parser.parse_args()

    binary = args.binary.resolve()
    build_dir = args.build_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    si_source = ROOT / "tests/scf/cases/bulk/diamondSi"
    fe_source = ROOT / "tests/scf/cases/bulk/bccFe"
    if args.quick:
        meshes = [1, 2]
        tiles = [1, 8]
        vectors = [False, True] if args.vectors else [False]
        fe_lengths = [1, 2]
    else:
        meshes = args.meshes
        tiles = args.tiles
        vectors = [False, True] if args.vectors else [False]
        fe_lengths = [1, 2, 3, 4, 5]

    rows: list[dict[str, Any]] = []
    folding: list[dict[str, Any]] = []
    validation: list[dict[str, Any]] = []
    configurations = [("lapack", "lapack")]
    if not args.skip_cuda:
        configurations.extend([("cuda", "zheevd_serial"), ("cuda", "zheevj_batched")])
    with tempfile.TemporaryDirectory(prefix="rslmto-accp0-", dir="/tmp") as temporary:
        scratch = Path(temporary)
        si_dir = scratch / "diamondSi"
        fixture_input(si_source, si_dir)
        for backend, solver_strategy in configurations:
            if backend == "cuda":
                validation.append(run_validation(
                    binary=binary, fixture_dir=si_dir, fixture="diamondSi", length=1,
                    mesh=(1, 1, 1), solver_strategy=solver_strategy, tile=max(1, tiles[0]),
                ))
        for mesh_value in meshes:
            for tile in tiles:
                for vectors_value in vectors:
                    for backend, solver_strategy in configurations:
                        output = raw_dir / f"si_{mesh_value}_{backend}_{solver_strategy}_{tile}_{int(vectors_value)}.json"
                        rows.append(
                            run_measurement(
                                binary=binary,
                                build_dir=build_dir,
                                fixture_dir=si_dir,
                                output=output,
                                fixture="diamondSi",
                                source="tests/scf/cases/bulk/diamondSi",
                                workload="crossover",
                                length=1,
                                mesh=(mesh_value, mesh_value, mesh_value),
                                backend=backend,
                                solver_strategy=solver_strategy,
                                tile=tile,
                                vectors=vectors_value,
                                warmups=args.warmups,
                                repetitions=args.repetitions,
                            )
                        )

        for length in fe_lengths:
            fe_dir = scratch / f"bccFe_L{length}"
            if length == 1:
                fixture_input(fe_source, fe_dir)
                source_name = "tests/scf/cases/bulk/bccFe"
            else:
                make_bcc_supercell(fe_source, fe_dir, length)
                source_name = "tests/scf/cases/bulk/bccFe:production_supercell"

            oracle_info: dict[str, Any] = {"fixture": "bccFe", "L": length}
            if length > 1:
                primitive_dir = scratch / f"bccFe_primitive_for_L{length}"
                fixture_input(fe_source, primitive_dir)
                primitive_dump = run_oracle(binary, primitive_dir, raw_dir / f"oracle_primitive_L{length}.log", (length, length, length), 1)
                super_dump = run_oracle(binary, fe_dir, raw_dir / f"oracle_supercell_L{length}.log", (1, 1, 1), length)
                oracle_info.update(compare_oracle(primitive_dump, super_dump))
                folding.append(oracle_info)
                if not oracle_info["passed"]:
                    raise RuntimeError(f"ACC-P0 band-folding oracle failed for bccFe L={length}: {oracle_info}")

            if length <= 3:
                for backend, solver_strategy in configurations:
                    if backend == "cuda":
                        validation.append(run_validation(
                            binary=binary, fixture_dir=fe_dir, fixture="bccFe", length=length,
                            mesh=(1, 1, 1), solver_strategy=solver_strategy, tile=max(1, tiles[0]),
                        ))

            workload_meshes = [("crossover", (mesh_value,) * 3) for mesh_value in meshes]
            workload_meshes.extend(("matched-density", matched_mesh(max(meshes), length)) for _ in [0])
            for workload, mesh in workload_meshes:
                for vectors_value in vectors:
                    for tile in tiles:
                        for backend, solver_strategy in configurations:
                            if backend == "cuda":
                                # The driver repeats this check immediately
                                # before solving; doing it here makes the
                                # campaign log explicit and keeps large jobs
                                # from being launched blindly.
                                preflight(binary, fe_dir, "bccFe", length, mesh, backend, solver_strategy, tile, vectors_value)
                            output = raw_dir / f"fe_L{length}_{workload}_{mesh[0]}_{backend}_{solver_strategy}_{tile}_{int(vectors_value)}.json"
                            rows.append(
                                run_measurement(
                                    binary=binary,
                                    build_dir=build_dir,
                                    fixture_dir=fe_dir,
                                    output=output,
                                    fixture="bccFe",
                                    source=source_name,
                                    workload=workload,
                                    length=length,
                                    mesh=mesh,
                                    backend=backend,
                                    solver_strategy=solver_strategy,
                                    tile=tile,
                                    vectors=vectors_value,
                                    warmups=args.warmups,
                                    repetitions=args.repetitions,
                                )
                            )

    rows.sort(key=lambda row: str(row.get("name", "")))
    with (output_dir / "accp0_table.csv").open("w", newline="", encoding="utf-8") as stream:
        fieldnames = REQUIRED_COLUMNS + ["name", "vectors"]
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "rslmto.accp0-real-material.v1",
        "policy": {
            "persistent_process": True,
            "warmups_inside_process": args.warmups,
            "measured_repetitions_inside_process": args.repetitions,
            "cuda_solver_algorithm_changed": True,
            "solver_strategies": sorted({str(row.get("solver_strategy")) for row in rows}),
        },
        "rows": rows,
        "folding_oracle": folding,
        "validation": validation,
        "summary": {
            "rows": len(rows),
            "configuration_totals": {
                f"{backend}:{strategy}": sum(
                    1 for row in rows
                    if row.get("backend") == backend and row.get("solver_strategy") == strategy
                )
                for backend, strategy in configurations
            },
            "steady_total_medians_by_backend": {
                f"{backend}:{strategy}": statistics.median(
                    float(row["total_steady_s"]) for row in rows
                    if row.get("backend") == backend
                    and row.get("solver_strategy") == strategy
                    and row.get("support_status") == "supported"
                )
                for backend, strategy in configurations
                if any(
                    row.get("backend") == backend
                    and row.get("solver_strategy") == strategy
                    and row.get("support_status") == "supported"
                    for row in rows
                )
            },
        },
    }
    (output_dir / "accp0_results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE {output_dir / 'accp0_table.csv'} ({len(rows)} rows)")
    print(f"WROTE {output_dir / 'accp0_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
