#!/usr/bin/env python3
"""Exercise the production bulk/local Hubbard insertion paths on matched Fe decks.

The caller supplies one-step bcc-Fe bulk/impurity input templates and the
accepted +U potential used as the common seed.  The test runs D1, J=0, host-only,
impurity-only via the general array, and impurity-only via the override array.
All runs are isolated and retain their input hashes, audit streams, traces, and
logs under --workdir.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from c2_hubbard_operator_audit import choose_record, read_records

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "c2_fe"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def labels(path: Path) -> list[str]:
    values = re.findall(r"(?im)^\s*label\(\d+\)\s*=\s*'([^']+)'", path.read_text())
    if not values:
        raise ValueError(f"no atom labels found in {path}")
    return values


def replace_hamiltonian(text: str, fields: list[str]) -> str:
    block = "&hamiltonian\n" + "\n".join(fields) + "\n/"
    pattern = re.compile(r"(?ims)^\s*&hamiltonian\b.*?^\s*/\s*$")
    if pattern.search(text):
        return pattern.sub(block, text, count=1)
    return text.rstrip() + "\n\n" + block + "\n"


def zero_seed(source: Path, target: Path) -> None:
    text = source.read_text()
    for key in ("hubbard_u", "hubbard_j"):
        text, count = re.subn(rf"(?im)^\s*{key}\s*=.*$", f" {key} = 0.0, 0.0, 0.0", text)
        if count != 1:
            raise ValueError(f"expected one {key} field in {source}, found {count}")
    target.write_text(text)


def matrix_trace(path: Path, potential_type: int) -> dict[tuple[int, int, int], complex]:
    values: dict[tuple[int, int, int], complex] = {}
    for line in path.read_text().splitlines():
        fields = line.split(",")
        if len(fields) != 19 or fields[0] != "ldm" or fields[4] != "P_OUT":
            continue
        if int(fields[5]) != 0 or int(fields[6]) != potential_type or int(fields[7]) != 2:
            continue
        key = (int(fields[8]), int(fields[11]), int(fields[12]))
        values[key] = complex(float(fields[17]), float(fields[18]))
    if len(values) != 50:
        raise ValueError(f"{path} has {len(values)} d-shell P_OUT entries, expected 50")
    return values


def run_case(binary: Path, root: Path, case: str, role: str, template: Path,
             fields: list[str], plus_seed: Path, impurity_seed_u0: bool = False,
             host_seed_u0: bool = False) -> Path:
    directory = root / case / role
    directory.mkdir(parents=True, exist_ok=False)
    text = replace_hamiltonian(template.read_text(), fields)
    input_path = directory / "input.nml"
    input_path.write_text(text)
    atom_labels = labels(input_path)
    for index, label in enumerate(atom_labels):
        seed = directory / f"{label}.nml"
        use_zero = (role == "imp" and index > 0 and impurity_seed_u0) or (role == "imp" and index == 0 and host_seed_u0)
        if use_zero:
            zero_seed(plus_seed, seed)
        else:
            shutil.copyfile(plus_seed, seed)
    manifest = [f"{sha(input_path)}  input.nml\n"]
    manifest += [f"{sha(directory / f'{label}.nml')}  {label}.nml\n" for label in atom_labels]
    (directory / "PRELAUNCH_SHA256SUMS.txt").write_text("".join(manifest))
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "RSLMTO_C2_HUBBARD_AUDIT": "1",
        "RSLMTO_C2_HUBBARD_AUDIT_FILE": str(directory / "c2_hubbard_audit.dat"),
        "RSLMTO_FIXED_U_TRACE": str(directory / "fixed_u_trace.csv"),
        "RSLMTO_FIXED_U_TRACE_ITERATION": "1",
        "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1",
        "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
    })
    with (directory / "testrun.log").open("w") as log:
        subprocess.run([str(binary)], cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=300)
    return directory


def record(directory: Path, assembler: str, site: int, potential_type: int):
    matches = read_records(directory / "c2_hubbard_audit.dat")
    return choose_record(matches, assembler, site, 1)


def maxdiff(a: dict, b: dict, field: str) -> float:
    return max(abs(a["block"][key][field] - b["block"][key][field]) for key in a["block"])


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def assert_pair(bulk_dir: Path, imp_dir: Path, label: str) -> tuple[dict, dict]:
    bulk = record(bulk_dir, "BULK_EE", 1, 1)
    local = record(imp_dir, "LOCAL_HALL", 1, 2)
    require(bulk["symbol"] == local["symbol"] == "Fe", f"{label}: target symbols differ")
    require(bulk["count"] == local["count"] == 1, f"{label}: insertion count is not one")
    require(abs(bulk["u_ry"] - local["u_ry"]) <= 1e-12, f"{label}: U mismatch")
    require(abs(bulk["j_ry"] - local["j_ry"]) <= 1e-12, f"{label}: J mismatch")
    for field in ("before", "correction", "after"):
        require(maxdiff(bulk, local, field) <= 1e-12, f"{label}: {field} block mismatch")
    return bulk, local


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--bulk-template", type=Path, default=FIXTURE_DIR / "input_bulk.nml")
    parser.add_argument("--impurity-template", type=Path, default=FIXTURE_DIR / "input_impurity.nml")
    parser.add_argument("--plus-u-seed", type=Path, default=FIXTURE_DIR / "Fe_plus_u_seed.nml")
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    bulk_template = args.bulk_template.resolve()
    imp_template = args.impurity_template.resolve()
    seed = args.plus_u_seed.resolve()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / "RUN_PROVENANCE.txt").write_text(
        f"binary {sha(binary)}  {binary}\n"
        f"bulk_template {sha(bulk_template)}  {bulk_template}\n"
        f"impurity_template {sha(imp_template)}  {imp_template}\n"
        f"plus_u_seed {sha(seed)}  {seed}\n"
        f"runner {sha(Path(__file__).resolve())}  {Path(__file__).resolve()}\n"
    )
    uplus = ["hubbard_u_general(1,3)=3.8", "hubbard_j_general(1,3)=1.1"]
    both_uplus = uplus + ["hubbard_u_general(2,3)=3.8", "hubbard_j_general(2,3)=1.1"]
    jzero = ["hubbard_u_general(1,3)=3.8", "hubbard_j_general(1,3)=0.0"]
    both_jzero = jzero + ["hubbard_u_general(2,3)=3.8", "hubbard_j_general(2,3)=0.0"]

    bulk_plus = run_case(binary, root, "d1_u38j11", "bulk", bulk_template, uplus, seed)
    imp_plus = run_case(binary, root, "d1_u38j11", "imp", imp_template, both_uplus, seed)
    bulk, local = assert_pair(bulk_plus, imp_plus, "D1 U=3.8 J=1.1")
    pout_b = matrix_trace(bulk_plus / "fixed_u_trace.csv", 1)
    pout_i = matrix_trace(imp_plus / "fixed_u_trace.csv", 2)
    pout_diff = max(abs(pout_b[key] - pout_i[key]) for key in pout_b)
    require(pout_diff <= 1e-5, f"D1 P_OUT mismatch {pout_diff:g} exceeds 1e-5 control bound")

    bulk_j0 = run_case(binary, root, "d2a_u38j0", "bulk", bulk_template, jzero, seed)
    imp_j0 = run_case(binary, root, "d2a_u38j0", "imp", imp_template, both_jzero, seed)
    assert_pair(bulk_j0, imp_j0, "D2(a) U=3.8 J=0")

    host_only = run_case(binary, root, "d2b_host_only", "imp", imp_template,
                         uplus + ["hubbard_u_general(2,3)=0.0", "hubbard_j_general(2,3)=0.0"],
                         seed, impurity_seed_u0=True)
    hrec = record(host_only, "LOCAL_HALL", 1, 2)
    require(abs(hrec["u_ry"]) <= 1e-12 and abs(hrec["j_ry"]) <= 1e-12, "D2(b) target impurity should have U=J=0")
    require(max(abs(value["correction"]) for value in hrec["block"].values()) <= 1e-12,
            "D2(b) host-only case applied an impurity correction")

    for case, general in (("d2c_imp_only_general", True), ("d2c_imp_only_override", False)):
        fields = ["hubbard_u_general(1,3)=0.0", "hubbard_j_general(1,3)=0.0"]
        if general:
            fields += ["hubbard_u_general(2,3)=3.8", "hubbard_j_general(2,3)=1.1"]
        else:
            fields += ["hubbard_u_impurity(1,3)=3.8", "hubbard_j_impurity(1,3)=1.1"]
        idir = run_case(binary, root, case, "imp", imp_template, fields, seed,
                        host_seed_u0=True, impurity_seed_u0=True)
        irec = record(idir, "LOCAL_HALL", 1, 2)
        require(abs(irec["u_ry"] - bulk["u_ry"]) <= 1e-12 and abs(irec["j_ry"] - bulk["j_ry"]) <= 1e-12,
                f"{case}: impurity target U/J mismatch")
        require(max(abs(value["correction"]) for value in irec["block"].values()) > 1e-6,
                f"{case}: target impurity correction is missing")
        require(maxdiff(bulk, irec, "correction") <= 1e-12, f"{case}: impurity correction differs from bulk")

    print(f"D1 block pre/correction/post maxdiff: {maxdiff(bulk, local, 'before'):.3e} / "
          f"{maxdiff(bulk, local, 'correction'):.3e} / {maxdiff(bulk, local, 'after'):.3e}")
    print(f"D1 P_OUT maxdiff: {pout_diff:.6e}")
    print("D2(a) J=0, D2(b) host-only, D2(c) impurity general/override: PASS")
    print(f"BINARY_SHA256={sha(binary)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"RESULT: FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
