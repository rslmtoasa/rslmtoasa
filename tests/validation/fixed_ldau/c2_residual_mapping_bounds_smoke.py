#!/usr/bin/env python3
"""Run a two-step bounds-checked smoke of the C2 impurity residual mapping.

The default Fe fixture generates ntype=2 with a recursion-site representative
whose physical cluster index exceeds ntype. This reaches the production
``update_fixed_u_residuals`` path and catches using that physical index as a
type id. Supply a Debug/bounds-enabled executable; this is a diagnostic smoke,
not a convergence or physics acceptance test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "c2_fe"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_substitution(pattern: str, replacement: str, value: str, label: str) -> str:
    result, count = re.subn(pattern, replacement, value, count=1, flags=re.IGNORECASE | re.MULTILINE)
    if count != 1:
        raise ValueError(f"expected exactly one {label}, found {count}")
    return result


def prepare_input(template: Path) -> str:
    text = template.read_text()
    text = one_substitution(r"^\s*nstep\s*=\s*[^\n]+$", "nstep = 2", text, "nstep field")
    hamiltonian = "&hamiltonian\nhubbard_u_impurity(1,3) = 3.8\nhubbard_j_impurity(1,3) = 1.1\n/"
    pattern = r"(?ims)^\s*&hamiltonian\b.*?^\s*/\s*$"
    text, count = re.subn(pattern, hamiltonian, text, count=1)
    if count != 1:
        raise ValueError("expected exactly one &hamiltonian block")
    return text


def validate_outnewclu(path: Path) -> dict[str, int | list[int]]:
    lines = path.read_text().splitlines()
    header_index = next((i for i, line in enumerate(lines) if "NTYPE=" in line and "NREC=" in line), None)
    if header_index is None:
        raise ValueError("outnewclu has no NTYPE/NREC control header")
    match = re.search(r"NTYPE=\s*(\d+)\s+NMAX=\s*(\d+)\s+NBAS=\s*(\d+)\s+NREC=\s*(\d+)", lines[header_index])
    if match is None:
        raise ValueError("could not parse outnewclu control header")
    ntype, nmax, nbas, nrec = map(int, match.groups())
    # The source writes nnmax(1..nbulk), then ibulk(1..nbulk), immediately
    # after this header. The latter are physical cluster indices, not types.
    numeric_lines: list[list[int]] = []
    for line in lines[header_index + 1 :]:
        if not line.strip():
            continue
        fields = line.split()
        if not all(re.fullmatch(r"[+-]?\d+", field) for field in fields):
            break
        numeric_lines.append([int(field) for field in fields])
        if len(numeric_lines) == 2:
            break
    if len(numeric_lines) != 2:
        raise ValueError("could not parse bulk representative indices from outnewclu")
    bulk_indices = numeric_lines[1]
    if ntype != 2 or nrec != 1:
        raise ValueError(f"fixture no longer exercises ntype=2,nrec=1 (got {ntype},{nrec})")
    if not bulk_indices or max(bulk_indices) <= ntype:
        raise ValueError(f"physical representative index does not exceed ntype: {bulk_indices}")
    return {"ntype": ntype, "nmax": nmax, "nbas": nbas, "nrec": nrec,
            "bulk_representative_indices": bulk_indices}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="Debug executable built with Fortran bounds checking")
    parser.add_argument("--workdir", type=Path, required=True, help="new directory for this isolated smoke")
    parser.add_argument("--input-template", type=Path, default=FIXTURE_DIR / "input_impurity.nml")
    parser.add_argument("--plus-u-seed", type=Path, default=FIXTURE_DIR / "Fe_plus_u_seed.nml")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()

    binary = args.binary.resolve()
    template = args.input_template.resolve()
    seed = args.plus_u_seed.resolve()
    root = args.workdir.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError(f"binary is not executable: {binary}")
    root.mkdir(parents=True, exist_ok=False)

    input_path = root / "input.nml"
    input_path.write_text(prepare_input(template))
    labels = re.findall(r"(?im)^\s*label\(\d+\)\s*=\s*'([^']+)'", input_path.read_text())
    if labels != ["Fe-bulk", "Fe-imp"]:
        raise ValueError(f"fixture labels changed; expected Fe-bulk/Fe-imp, got {labels}")
    for label in labels:
        shutil.copyfile(seed, root / f"{label}.nml")

    provenance = {
        "purpose": "two-step diagnostic smoke only; no convergence/physics acceptance",
        "binary": str(binary),
        "binary_sha256": sha(binary),
        "template": str(template),
        "template_sha256": sha(template),
        "plus_u_seed": str(seed),
        "plus_u_seed_sha256": sha(seed),
        "input_sha256": sha(input_path),
        "staged_seed_sha256": {label: sha(root / f"{label}.nml") for label in labels},
        "nstep": 2,
        "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                        "fixed_u_trace": "unset", "fixed_u_trace_ql_only": "unset"},
        "command": [str(binary)],
        "status": "PREPARED",
    }
    (root / "PRELAUNCH.json").write_text(json.dumps(provenance, indent=2) + "\n")
    env = os.environ.copy()
    for key in list(env):
        if (key.startswith("RSLMTO_FIXED_U_TRACE") or key.startswith("RSLMTO_C2_HUBBARD_AUDIT") or
                key == "RSLMTO_FIXED_U_RUN_ID"):
            env.pop(key)
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    try:
        with (root / "console.log").open("w") as log:
            completed = subprocess.run([str(binary)], cwd=root, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        provenance.update({"status": "TIMEOUT", "timeout_seconds": args.timeout})
        (root / "RESULT.json").write_text(json.dumps(provenance, indent=2) + "\n")
        raise ValueError(f"smoke timed out; artifacts retained at {root}") from exc

    console = (root / "console.log").read_text(errors="replace")
    residuals = re.findall(r"FIXED_U_RESIDUAL\b[^\n]*", console)
    site_type_rows = re.findall(r"HUBBARD_MOMENT_ENERGY_SITE\s+site=\s*1\s+type=\s*2\b[^\n]*", console)
    geometry = validate_outnewclu(root / "outnewclu")
    errors = re.findall(r"(?im)^.*(?:Fortran runtime error|ERROR STOP|Index .* above upper bound|index .* above upper bound).*$", console)
    finished = "Calculation finished." in console

    failures = []
    if completed.returncode != 0:
        failures.append(f"exit code {completed.returncode}")
    if errors:
        failures.append("bounds/runtime errors present")
    if len(residuals) != 2:
        failures.append(f"expected two production residual updates, found {len(residuals)}")
    if len(site_type_rows) != 2:
        failures.append(f"expected two site=1/type=2 records, found {len(site_type_rows)}")
    if not finished:
        failures.append("missing normal calculation completion marker")

    provenance.update({"status": "PASS" if not failures else "FAIL", "returncode": completed.returncode,
                       "completed_marker": finished, "residual_updates": residuals,
                       "site_type_records": site_type_rows, "bounds_errors": errors,
                       "geometry": geometry, "console_sha256": sha(root / "console.log"),
                       "outnewclu_sha256": sha(root / "outnewclu"),
                       "report_sha256": sha(root / "report.out") if (root / "report.out").exists() else None,
                       "failures": failures})
    (root / "RESULT.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if failures:
        raise ValueError(f"residual mapping bounds smoke failed: {'; '.join(failures)}; artifacts: {root}")

    print(f"PASS: two residual updates; recursion site 1 maps to type 2; physical representative index "
          f"{max(geometry['bulk_representative_indices'])} exceeds ntype={geometry['ntype']}")
    print(f"BINARY_SHA256={sha(binary)}")
    print(f"ARTIFACTS={root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"RESULT: FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
