#!/usr/bin/env python3
"""Benchmark the supplied explicit bcc-Fe supercells with ACC-P0.

The source fixtures are never modified. Each fixture is copied to a temporary
run directory with its corrected ``crystal_sym='file'`` input unchanged, so
the explicit ``lattice.nml`` (with ``nrec=L^3``) is consumed by the production
lattice path. The measured copies use one canonical ``Fe1.nml`` potential for
every site and every selected supercell; only the site label is changed. This
keeps the scaling campaign a uniform-potential comparison even when the source
directories contain site-dependent self-consistent ``FeX.nml`` files.

The campaign records persistent CPU/GPU timings at the scheduled k meshes and performs two
independent correctness checks:

* CPU LAPACK versus CUDA eigenvalues for the identical assembled Hamiltonian;
* a CPU band-folding comparison using the same canonical ``Fe1.nml`` potential,
  which isolates geometry/assembler correctness.

The default k-point policy is a geometric workload schedule anchored at a
32x32x32 primitive-cell mesh and halved per added supercell level: 16x16x16,
8x8x8, 4x4x4, and 2x2x2 for L=2,3,4,5 respectively.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from accp0_real_material import (  # noqa: E402
    REQUIRED_COLUMNS,
    compare_oracle,
    preflight,
    run_measurement,
    run_oracle,
)


FIXTURE_PATTERN = re.compile(r"^(?P<length>[2-5])x(?P=length)x(?P=length)$")
CRYSTAL_SYM_PATTERN = re.compile(
    r"(?m)^(?P<prefix>\s*crystal_sym\s*=\s*)['\"](?P<value>[^'\"]+)['\"](?P<suffix>\s*(?:!.*)?)$"
)


def parse_int_list(text: str) -> list[int]:
    values = [int(part) for part in text.split(",") if part]
    if not values or any(value < 1 for value in values):
        raise argparse.ArgumentTypeError("expected a comma-separated list of positive integers")
    return values


def density_mesh(base_mesh: int, length: int) -> tuple[int, int, int]:
    value = max(1, base_mesh // (2 ** (length - 1)))
    return (value, value, value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_potential_digest(path: Path) -> str:
    """Hash a potential after removing only the site label from ``symbol``."""

    text = path.read_text(encoding="utf-8")
    text = re.sub(r"(symbol\s*=\s*)['\"][^'\"]+['\"]", r"\1'FeX'", text, count=1)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_manifest(path: Path, source_root: Path) -> list[dict[str, Any]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "rslmto.accp0-supercell-fe-examples.v1":
        raise ValueError(f"unsupported supercell fixture manifest: {path}")
    fixtures = document.get("fixtures")
    if not isinstance(fixtures, list) or len(fixtures) != 4:
        raise ValueError("the supercell manifest must contain the four supplied fixtures")

    validated: list[dict[str, Any]] = []
    for fixture in fixtures:
        directory = str(fixture["directory"])
        match = FIXTURE_PATTERN.fullmatch(directory)
        if match is None:
            raise ValueError(f"invalid supercell directory in manifest: {directory}")
        length = int(fixture["length"])
        atom_count = int(fixture["atom_count"])
        if length != int(match.group("length")) or atom_count != length**3:
            raise ValueError(f"inconsistent manifest entry: {fixture}")

        source = source_root / directory
        input_path = source / "input.nml"
        lattice_path = source / "lattice.nml"
        if not source.is_dir() or not input_path.is_file() or not lattice_path.is_file():
            raise FileNotFoundError(f"incomplete supplied fixture: {source}")
        potentials = sorted(source.glob("Fe*.nml"))
        if len(potentials) != atom_count:
            raise ValueError(f"{source}: expected {atom_count} Fe*.nml files, found {len(potentials)}")
        if any(path.name in {"Fe.nml", "Fe_out.nml", "Fe_scf.nml"} for path in potentials):
            raise ValueError(f"{source}: fixture potential set is not site-labelled")
        potential_digests = {normalized_potential_digest(path) for path in potentials}

        input_text = input_path.read_text(encoding="utf-8")
        ntype_match = re.search(r"(?m)^\s*ntype\s*=\s*(\d+)", input_text)
        if ntype_match is None or int(ntype_match.group(1)) != atom_count:
            raise ValueError(f"{input_path}: ntype does not match the fixture atom count")
        nstep_match = re.search(r"(?m)^\s*nstep\s*=\s*(\d+)", input_text)
        if nstep_match is None or int(nstep_match.group(1)) < 5:
            raise ValueError(f"{input_path}: nstep must be at least 5 for benchmark timing policy")
        crystal_match = CRYSTAL_SYM_PATTERN.search(input_text)
        if crystal_match is None:
            raise ValueError(f"{input_path}: crystal_sym is missing")

        validated.append(
            {
                **fixture,
                "source": source,
                "input_sha256": sha256(input_path),
                "lattice_sha256": sha256(lattice_path),
                "source_crystal_sym": crystal_match.group("value"),
                "potential_count": len(potentials),
                "potential_unique_normalized": len(potential_digests),
                "potential_uniform_after_symbol_normalization": len(potential_digests) == 1,
                "input_nstep": int(nstep_match.group(1)),
            }
        )
    return validated


def stage_fixture(
    source: Path, destination: Path, canonical_potential: Path
) -> dict[str, Any]:
    """Stage one fixture with a common potential and assert its lattice route."""

    import shutil

    shutil.copytree(source, destination)
    input_path = destination / "input.nml"
    original = input_path.read_text(encoding="utf-8")
    match = CRYSTAL_SYM_PATTERN.search(original)
    if match is None:
        raise ValueError(f"{input_path}: crystal_sym is missing")
    if match.group("value").lower() != "file":
        raise ValueError(
            f"{input_path}: expected corrected crystal_sym='file', found {match.group('value')!r}"
        )
    canonical_text = canonical_potential.read_text(encoding="utf-8")
    for potential in sorted(destination.glob("Fe*.nml")):
        text = re.sub(
            r"(symbol\s*=\s*)['\"][^'\"]+['\"]",
            rf"\1'{potential.stem}'",
            canonical_text,
            count=1,
        )
        potential.write_text(text, encoding="utf-8")
    return {
        "source_crystal_sym": match.group("value"),
        "staged_crystal_sym": match.group("value"),
        "input_sha256_staged": sha256(input_path),
        "canonical_potential_source": str(canonical_potential),
        "staged_potential_unique_normalized": len(
            {normalized_potential_digest(path) for path in destination.glob("Fe*.nml")}
        ),
        "adapter_applied": True,
    }


def make_uniform_folding_fixtures(
    staged_supercell: Path, primitive_destination: Path, uniform_destination: Path
) -> None:
    """Make a same-potential primitive/supercell pair for a folding check."""

    import shutil

    primitive_input = (staged_supercell / "input.nml").read_text(encoding="utf-8")
    primitive_input = re.sub(
        r"(?m)^(\s*crystal_sym\s*=\s*)['\"][^'\"]+['\"]",
        r"\1'bcc'",
        primitive_input,
        count=1,
    )
    primitive_input = re.sub(r"(?m)^(\s*ntype\s*=\s*)\d+", r"\g<1>1", primitive_input, count=1)
    primitive_input = re.sub(r"(?m)^\s*label\(\d+\)\s*=.*\n", "", primitive_input)
    primitive_input = primitive_input.replace("  label(1) = 'Fe'", "  label(1) = 'Fe'")
    atoms_marker = re.search(r"(?m)^&atoms\s*$", primitive_input)
    if atoms_marker is None:
        raise ValueError(f"{staged_supercell / 'input.nml'}: missing &atoms section")
    section_end = re.search(r"(?m)^/\s*$", primitive_input[atoms_marker.end():])
    if section_end is None:
        raise ValueError(f"{staged_supercell / 'input.nml'}: unterminated &atoms section")
    atoms_start = atoms_marker.end()
    atoms_end = atoms_start + section_end.start()
    atoms_section = primitive_input[atoms_start:atoms_end]
    database_match = re.search(r"(?m)^\s*database\s*=.*\n", atoms_section)
    database_line = database_match.group(0) if database_match else "  database = './'\n"
    primitive_input = primitive_input[:atoms_start] + database_line + "  label(1) = 'Fe'\n" + primitive_input[atoms_end:]
    (primitive_destination / "input.nml").write_text(primitive_input, encoding="utf-8")
    first_potential = staged_supercell / "Fe1.nml"
    primitive_potential = first_potential.read_text(encoding="utf-8")
    primitive_potential = re.sub(
        r"(symbol\s*=\s*)['\"][^'\"]+['\"]", r"\1'Fe'", primitive_potential, count=1
    )
    (primitive_destination / "Fe.nml").write_text(primitive_potential, encoding="utf-8")

    shutil.copytree(staged_supercell, uniform_destination)
    for potential in sorted(uniform_destination.glob("Fe*.nml")):
        text = first_potential.read_text(encoding="utf-8")
        text = re.sub(
            r"(symbol\s*=\s*)['\"][^'\"]+['\"]",
            rf"\1'{potential.stem}'",
            text,
            count=1,
        )
        potential.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--example-root", type=Path, default=ROOT / "example/bulk/supercellFe")
    parser.add_argument(
        "--manifest", type=Path, default=Path(__file__).resolve().with_name("supercellFe_accp0_examples.json")
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/benchmarks/accp0-supercellFe")
    parser.add_argument("--skip-cuda", action="store_true")
    parser.add_argument("--skip-folding", action="store_true")
    parser.add_argument("--lengths", type=parse_int_list, default=[2, 3, 4, 5])
    parser.add_argument(
        "--base-mesh",
        type=int,
        default=32,
        help="mesh density for L=1; each increment in L halves the per-axis mesh",
    )
    parser.add_argument("--tiles", type=parse_int_list, default=[1])
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--vectors", action="store_true", help="also benchmark eigenvector requests")
    args = parser.parse_args()

    binary = args.binary.resolve()
    build_dir = args.build_dir.resolve()
    source_root = args.example_root.resolve()
    manifest_path = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)
    if args.warmups < 0 or args.repetitions < 5 or args.base_mesh < 1:
        raise ValueError("warmups must be non-negative, repetitions must be at least 5, and base_mesh must be positive")
    if any(length not in {2, 3, 4, 5} for length in args.lengths):
        raise ValueError("lengths must be drawn from the supplied 2x2x2 through 5x5x5 fixtures")

    fixtures = load_manifest(manifest_path, source_root)
    fixtures_by_length = {int(fixture["length"]): fixture for fixture in fixtures}
    selected = [fixtures_by_length[length] for length in args.lengths]
    canonical_potential = fixtures_by_length[2]["source"] / "Fe1.nml"
    if not canonical_potential.is_file():
        raise FileNotFoundError(canonical_potential)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    backends = ["lapack"] if args.skip_cuda else ["lapack", "cuda"]
    vectors = [False, True] if args.vectors else [False]
    rows: list[dict[str, Any]] = []
    backend_agreement: list[dict[str, Any]] = []
    folding_oracle: list[dict[str, Any]] = []
    adaptations: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="rslmto-accp0-supercellFe-", dir="/tmp") as temporary:
        scratch = Path(temporary)
        for fixture in selected:
            length = int(fixture["length"])
            mesh = density_mesh(args.base_mesh, length)
            primitive_mesh = tuple(length * value for value in mesh)
            fixture_dir = scratch / str(fixture["name"])
            adaptation = stage_fixture(
                Path(fixture["source"]), fixture_dir, canonical_potential
            )
            adaptations.append(
                {
                    "fixture": fixture["name"],
                    "directory": fixture["directory"],
                    "source_input_sha256": fixture["input_sha256"],
                    "potential_unique_normalized": fixture["potential_unique_normalized"],
                    "potential_uniform_after_symbol_normalization": fixture[
                        "potential_uniform_after_symbol_normalization"
                    ],
                    "input_nstep": fixture["input_nstep"],
                    **adaptation,
                }
            )

            cpu_values = run_oracle(
                binary, fixture_dir, raw_dir / f"{fixture['name']}_oracle_cpu.log", mesh, length, "lapack"
            )
            if not args.skip_cuda:
                gpu_values = run_oracle(
                    binary, fixture_dir, raw_dir / f"{fixture['name']}_oracle_cuda.log", mesh, length, "cuda"
                )
                # For very dense meshes, tiny backend-dependent perturbations
                # can split a nominally degenerate group despite sub-1e-12
                # eigenvalue agreement. State degeneracy matching as a
                # diagnostic, but make the backend correctness gate the
                # sorted eigenvalue multiset and state count.
                agreement = compare_oracle(cpu_values, gpu_values, require_degeneracy=False)
                agreement.update(
                    {
                        "fixture": fixture["name"],
                        "L": length,
                        "check": "cpu_lapack_vs_cuda",
                        "tolerance_energy": 1.0e-7,
                    }
                )
                backend_agreement.append(agreement)

            if not args.skip_folding:
                primitive_dir = scratch / f"{fixture['name']}_uniform_primitive"
                uniform_dir = scratch / f"{fixture['name']}_uniform_supercell"
                primitive_dir.mkdir()
                make_uniform_folding_fixtures(fixture_dir, primitive_dir, uniform_dir)
                primitive_values = run_oracle(
                    binary,
                    primitive_dir,
                    raw_dir / f"{fixture['name']}_uniform_folding_primitive.log",
                    primitive_mesh,
                    1,
                    "lapack",
                )
                uniform_values = run_oracle(
                    binary,
                    uniform_dir,
                    raw_dir / f"{fixture['name']}_uniform_folding_supercell.log",
                    mesh,
                    length,
                    "lapack",
                )
                folding = compare_oracle(primitive_values, uniform_values)
                folding.update(
                    {
                        "fixture": fixture["name"],
                        "L": length,
                        "check": "uniform_Fe1_cpu_band_folding",
                        "interpretation": (
                            "pass: explicit supercell geometry folds to the same-potential primitive spectrum"
                            if folding["passed"]
                            else "fail: explicit supercell geometry or assembler is not folding consistently"
                        ),
                    }
                )
                folding_oracle.append(folding)

            for backend in backends:
                for tile in args.tiles:
                    for vectors_value in vectors:
                        if backend == "cuda":
                            preflight(
                                binary,
                                fixture_dir,
                                str(fixture["name"]),
                                length,
                                mesh,
                                backend,
                                tile,
                                vectors_value,
                            )
                        output = raw_dir / (
                            f"{fixture['name']}_k_density_{backend}_tile{tile}_v{int(vectors_value)}.json"
                        )
                        rows.append(
                            run_measurement(
                                binary=binary,
                                build_dir=build_dir,
                                fixture_dir=fixture_dir,
                                output=output,
                                fixture=str(fixture["name"]),
                                source=f"example/bulk/supercellFe/{fixture['directory']}",
                                workload="k-density",
                                length=length,
                                mesh=mesh,
                                backend=backend,
                                tile=tile,
                                vectors=vectors_value,
                                warmups=args.warmups,
                                repetitions=args.repetitions,
                            )
                        )

    rows.sort(key=lambda row: (int(row.get("L", 0)), str(row.get("backend", "")), int(row.get("vectors", 0))))
    with (output_dir / "supercellFe_accp0_table.csv").open("w", newline="", encoding="utf-8") as stream:
        fieldnames = REQUIRED_COLUMNS + ["name", "vectors"]
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    speedups: dict[str, float] = {}
    for length in args.lengths:
        cpu = [float(row["total_steady_s"]) for row in rows if row.get("L") == length and row.get("backend") == "lapack"]
        cuda = [float(row["total_steady_s"]) for row in rows if row.get("L") == length and row.get("backend") == "cuda"]
        if cpu and cuda:
            speedups[str(length)] = statistics.median(cpu) / statistics.median(cuda)

    report = {
        "schema": "rslmto.accp0-supercell-fe.v1",
        "manifest": str(manifest_path),
        "source_root": str(source_root),
        "policy": {
            "persistent_process": True,
            "warmups_inside_process": args.warmups,
            "measured_repetitions_inside_process": args.repetitions,
            "minimum_input_nstep": 5,
            "minimum_measured_repetitions": 5,
            "mesh_policy": f"floor({args.base_mesh}/2**(L-1)) per axis",
            "base_mesh": f"{args.base_mesh}x{args.base_mesh}x{args.base_mesh}",
            "cuda_solver_algorithm_changed": False,
            "source_inputs_staged_unchanged": True,
            "source_potentials_replaced_in_staging": True,
            "canonical_potential_source": str(canonical_potential),
        },
        "adaptations": adaptations,
        "rows": rows,
        "correctness": {
            "backend_agreement": backend_agreement,
            "uniform_band_folding_oracle": folding_oracle,
        },
        "summary": {
            "rows": len(rows),
            "backend_totals": {backend: sum(1 for row in rows if row.get("backend") == backend) for backend in backends},
            "cpu_over_cuda_speedup_by_L": speedups,
            "source_potentials_uniform_after_symbol_normalization": all(
                bool(fixture["potential_uniform_after_symbol_normalization"]) for fixture in selected
            ),
            "staged_potentials_uniform": all(
                int(item["staged_potential_unique_normalized"]) == 1 for item in adaptations
            ),
            "backend_agreement_passed": all(item["passed"] for item in backend_agreement) if backend_agreement else None,
            "uniform_folding_passed": all(item["passed"] for item in folding_oracle) if folding_oracle else None,
        },
    }
    (output_dir / "supercellFe_accp0_results.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"WROTE {output_dir / 'supercellFe_accp0_table.csv'} ({len(rows)} rows)")
    print(f"WROTE {output_dir / 'supercellFe_accp0_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
