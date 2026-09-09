#!/usr/bin/env python3
"""Run the reproducible TDVAL-01 Fe/Ni mesh and q-path ladder.

The source material decks remain the provenance anchor.  Each execution gets
an isolated effective deck and q file, so a failed or interrupted run cannot
leave a mixed-prefix result that looks converged.  Raw output is deliberately
written below ``results/validation/TDVAL-01_FE_NI`` and is not a checked-in
golden response.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


MATERIALS = {
    "fe": {
        "label": "bcc Fe",
        "directory": Path("bccFe"),
        "deck": "input_eigenpairs.nml",
        "alat": 2.86120,
        "direction": "[011]",
    },
    "ni": {
        "label": "fcc Ni",
        "directory": Path("fccNi"),
        "deck": "input_eigenpairs.nml",
        "alat": 3.520,
        "direction": "[-111]",
    },
}


def replace_value(text: str, key: str, value: str, quoted: bool = False) -> str:
    replacement = f"'{value}'" if quoted else value
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=\s*)[^!\n]*(.*)$")
    updated, count = pattern.subn(rf"\g<1>{replacement}\g<2>", text, count=1)
    if count != 1:
        raise ValueError(f"{key} is missing from the source deck")
    return updated


def namelist_value(text: str, key: str) -> str | None:
    match = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*([^!\n]+)", text)
    return match.group(1).strip().strip("'\"") if match else None


def q_text(values: list[tuple[float, float, float]]) -> str:
    lines = [str(len(values))]
    lines.extend("%.16g %.16g %.16g" % value for value in values)
    return "\n".join(lines) + "\n"


def prepare_run(source_deck: Path, run_dir: Path, mesh: int, q_values: list[tuple[float, float, float]],
                prefix: str, q_set: str) -> Path:
    text = source_deck.read_text(encoding="utf-8")
    database = namelist_value(text, "database")
    if database is None:
        raise ValueError(f"{source_deck}: database is missing")
    database_path = Path(database)
    if not database_path.is_absolute():
        database_path = (source_deck.parent / database_path).resolve()
    if not database_path.is_dir():
        raise FileNotFoundError(database_path)

    relative_database = Path(os.path.relpath(database_path, run_dir)).as_posix()
    if not relative_database.endswith("/"):
        relative_database += "/"
    text = replace_value(text, "database", relative_database, quoted=True)
    text = replace_value(text, "n1", str(mesh))
    text = replace_value(text, "n2", str(mesh))
    text = replace_value(text, "n3", str(mesh))
    text = replace_value(text, "nk1", str(mesh))
    text = replace_value(text, "nk2", str(mesh))
    text = replace_value(text, "nk3", str(mesh))
    text = replace_value(text, "q_file", "q_points.dat", quoted=True)
    text = replace_value(text, "output_prefix", prefix, quoted=True)
    text = replace_value(text, "output_xi", ".true.")
    text = replace_value(text, "output_chi", ".true.")
    text = replace_value(text, "omega_min", "-0.002" if q_set == "covariance" else "0.0")
    text = replace_value(text, "omega_max", "0.002" if q_set == "covariance" else "0.02")
    text = replace_value(text, "nomega", "9" if q_set == "covariance" else "101")
    text = replace_value(text, "eta", "0.0002")
    text = replace_value(text, "output_modes", ".false." if q_set == "covariance" else ".true.")
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "input.nml").write_text(text, encoding="utf-8")
    (run_dir / "q_points.dat").write_text(q_text(q_values), encoding="utf-8")
    shutil.copy2(source_deck, run_dir / "source_input.nml")
    return run_dir / "input.nml"


def run_one(binary: Path, run_dir: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        [str(binary), "input.nml"],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    (run_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    return {
        "status": "PASS" if completed.returncode == 0 else "FAIL",
        "returncode": completed.returncode,
        "directory": str(run_dir),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument("--binary", type=Path, default=None)
    parser.add_argument("--mesh", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--material", choices=["fe", "ni", "both"], default="both")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--q-set", choices=["both", "commensurate", "arbitrary", "covariance"], default="both")
    args = parser.parse_args()
    repo = args.repo.resolve()
    binary = (args.binary or repo / "build" / "bin" / "rslmto.x").resolve()
    output_root = repo / "results" / "validation" / "TDVAL-01_FE_NI" / "runs"
    selected = ["fe", "ni"] if args.material == "both" else [args.material]
    records: list[dict[str, Any]] = []

    for material in selected:
        spec = MATERIALS[material]
        source_dir = repo / "tests" / "regression" / "tddft_validation" / "materials" / spec["directory"]
        source_deck = source_dir / spec["deck"]
        for mesh in args.mesh:
            # q=(0,1/N,2/N) is commensurate with the mesh under the direct
            # reciprocal-coordinate convention.
            commensurate = [(0.0, 0.0, 0.0), (1.0 / mesh, 0.0, 0.0), (2.0 / mesh, 0.0, 0.0)]
            arbitrary = [(0.0, 0.0, 0.0), (0.01375, 0.0, 0.0)]
            q_sets = [("commensurate", commensurate), ("arbitrary", arbitrary),
                      ("covariance", [(0.0, 0.0, 0.0), (0.01375, 0.0, 0.0), (-0.01375, 0.0, 0.0)])]
            if args.q_set != "both":
                q_sets = [item for item in q_sets if item[0] == args.q_set]
            for label, q_values in q_sets:
                run_dir = output_root / material / f"{label}_nk{mesh:02d}"
                prefix = f"tdval01_{material}_{label}_nk{mesh:02d}"
                prepare_run(source_deck, run_dir, mesh, q_values, prefix, label)
                record = {"material": spec["label"], "mesh": [mesh, mesh, mesh], "q_set": label, "q_direct": q_values}
                record.update(run_one(binary, run_dir))
                records.append(record)
                print(f"{record['status']:4s} {material} {label} N={mesh}")
                if record["status"] == "FAIL" and not args.continue_on_error:
                    summary = {"binary": str(binary), "runs": records}
                    (output_root / "summary.json").parent.mkdir(parents=True, exist_ok=True)
                    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
                    return 1

    summary = {"binary": str(binary), "runs": records}
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0 if all(record["status"] == "PASS" for record in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
