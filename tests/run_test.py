#!/usr/bin/env python3
"""
Unified runner for RSLMTO example tests.

Reads a test case from a cases.json manifest, sets up a scratch directory,
patches input.nml with the test parameters, runs the binary, checks for
fatal errors, and optionally compares or saves reference data.

Usage (smoke):
  run_test.py --binary BIN --cases-json JSON --case-name NAME --scratch-root DIR

Usage (reference comparison):
  run_test.py ... --compare-ref REF_DIR [--abs-tol 1e-6] [--rel-tol 1e-6]

Usage (save reference):
  run_test.py ... --gen-ref REF_DIR

Reference data is driven by a "checks" dict in cases.json:

  "checks": {
    "nml": [
      {
        "file": "Fe_out.nml",
        "scalars": ["etot", "ws_r", "vmad"],
        "arrays": { "moment": [1, 2], "dos": [1, 5, 10] }
      }
    ],
    "text": [
      {
        "file": "Fe_dos.out", "rows": [50, 100], "cols": [1, 2],
        "tolerances": { "1": { "abs_tol": 1e-3 }, "2": { "abs_tol": 1e-4 } }
      }
    ],
    "dimensions": [
      { "file": "band_structure.dat", "data_rows": 201, "data_columns": 19 }
    ],
    "finite": [
      { "file": "totaldos.out", "cols": [1, 2] }
    ],
    "assertions": [
      { "file": "testrun.log", "pattern": "N\\(Emax\\)-N\\(Emin\\) =\\s*([-+0-9.EeDd]+)", "value": 16.0, "abs_tol": 1e-6 }
    ],
    "log": [
      {
        "file": "testrun.log",
        "patterns": { "fermi_level": "Canonical k-space occupations: EF=\\s*([-+0-9.EeDd]+)" },
        "tolerances": { "fermi_level": { "abs_tol": 1e-3 } }
      }
    ],
    "outputs": [
      { "file": "totaldos.out", "min_lines": 50 }
    ]
  }

If a case has no "checks" key the test still runs as a smoke test.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys

import f90nml

from manifest_defaults import apply_manifest_defaults


# DOS text files are written by the Fortran code with five digits after the
# decimal point (for example, ``(2f16.5)``). Their references therefore carry
# a one-last-printed-digit quantization uncertainty. Keep this floor local to
# text DOS output; namelist and log comparisons retain their requested
# tolerances.
DOS_TEXT_ABS_TOL = 1.1e-5


# ---------------------------------------------------------------------------
# Case loading
# ---------------------------------------------------------------------------

def load_case(cases_json: str, case_name: str) -> dict:
    with open(cases_json) as fh:
        data = json.load(fh)
    for case in data["cases"]:
        if case["name"] == case_name:
            return apply_manifest_defaults(data, case)
    raise KeyError(f"Case '{case_name}' not found in {cases_json}")


# ---------------------------------------------------------------------------
# Scratch directory setup
# ---------------------------------------------------------------------------

def setup_scratch(case_dir: str, workdir: str, preserve_outputs: bool = False) -> None:
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    shutil.copytree(case_dir, workdir)
    patterns = ["data.nml"]
    if not preserve_outputs:
        patterns.append("*_out.nml")
    for pattern in patterns:
        for path in glob.glob(os.path.join(workdir, pattern)):
            os.remove(path)


# ---------------------------------------------------------------------------
# NML patching
# ---------------------------------------------------------------------------

def patch_input_nml(workdir: str, case: dict) -> None:
    input_path = os.path.join(workdir, "input.nml")
    tmp_path = input_path + ".tmp"
    f90nml.patch(input_path, case["namelists"], tmp_path)
    os.replace(tmp_path, input_path)


# ---------------------------------------------------------------------------
# Binary invocation
# ---------------------------------------------------------------------------

def run_binary(binary: str, workdir: str, mpi_procs: int = 1, serial_omp_threads: int | None = None) -> None:
    run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_binary.sh")
    # MSYS Python uses the native Windows process API, where the POSIX path
    # /bin/bash is not executable.  Resolve bash through its MSYS-provided
    # PATH there; Unix keeps the explicit system interpreter.
    bash = shutil.which("bash") if os.name == "nt" else "/bin/bash"
    if bash is None:
        raise SystemExit("ERROR: bash is required to run test binaries on Windows/MSYS2")
    cmd = [bash, run_script, binary]
    if mpi_procs > 1:
        cmd.append(str(mpi_procs))
    # propagate controlled serial OMP thread count to the wrapper via env
    env = os.environ.copy()
    if mpi_procs <= 1 and serial_omp_threads is not None:
        env["RSLMTO_OMP_THREADS_SERIAL"] = str(serial_omp_threads)
    result = subprocess.run(cmd, cwd=workdir, env=env)
    if result.returncode != 0:
        print(f"ERROR: binary returned non-zero exit code {result.returncode}")
        log_path = os.path.join(workdir, "testrun.log")
        if os.path.exists(log_path):
            with open(log_path) as fh:
                content = fh.read()
            tail = content[-3000:] if len(content) > 3000 else content
            print("--- testrun.log (last 3000 chars) ---")
            print(tail)
            print("--------------------------------------")
        else:
            print("  (testrun.log not found)")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Post-run checks
# ---------------------------------------------------------------------------

def check_log(workdir: str, case_name: str) -> None:
    log_path = os.path.join(workdir, "testrun.log")
    if not os.path.exists(log_path):
        print(f"ERROR [{case_name}]: testrun.log missing in {workdir}")
        sys.exit(1)
    with open(log_path) as fh:
        content = fh.read()
    if "fatal" in content.lower():
        print(f"ERROR [{case_name}]: fatal message in testrun.log")
        lines = content.splitlines()
        print("\n".join(lines[-50:]))
        sys.exit(1)


def check_output_files(workdir: str, case_name: str, output_checks: list[dict]) -> None:
    """Require one complete shared output for each root-owned file contract."""
    for output_check in output_checks:
        filename = output_check["file"]
        matches = glob.glob(os.path.join(workdir, filename))
        if len(matches) != 1:
            print(f"ERROR [{case_name}]: expected exactly one '{filename}', found {len(matches)}")
            sys.exit(1)
        with open(matches[0]) as fh:
            line_count = sum(1 for _ in fh)
        min_lines = output_check.get("min_lines", 1)
        if line_count < min_lines:
            print(f"ERROR [{case_name}]: '{filename}' incomplete ({line_count} < {min_lines} lines)")
            sys.exit(1)


def check_text_dimensions(workdir: str, case_name: str, dimension_checks: list[dict]) -> None:
    """Check exact dimensions of the numeric records in a text output."""
    for dimension_check in dimension_checks:
        filename = dimension_check["file"]
        filepath = os.path.join(workdir, filename)
        data_rows = []
        with open(filepath) as fh:
            for line_number, line in enumerate(fh, start=1):
                fields = line.split("#", 1)[0].split()
                if not fields:
                    continue
                data_rows.append((line_number, fields))

        expected_rows = dimension_check["data_rows"]
        expected_columns = dimension_check["data_columns"]
        if len(data_rows) != expected_rows:
            print(
                f"ERROR [{case_name}]: '{filename}' has {len(data_rows)} data rows "
                f"(expected {expected_rows})"
            )
            sys.exit(1)
        wrong_columns = [(line_number, len(fields)) for line_number, fields in data_rows
                         if len(fields) != expected_columns]
        if wrong_columns:
            line_number, actual_columns = wrong_columns[0]
            print(
                f"ERROR [{case_name}]: '{filename}' line {line_number} has "
                f"{actual_columns} data columns (expected {expected_columns})"
            )
            sys.exit(1)


def check_text_finite(workdir: str, case_name: str, finite_checks: list[dict]) -> None:
    """Require selected numeric columns in text outputs to remain finite."""
    for finite_check in finite_checks:
        filename = finite_check["file"]
        filepath = os.path.join(workdir, filename)
        requested_cols = finite_check.get("cols")
        with open(filepath) as fh:
            for line_number, line in enumerate(fh, start=1):
                fields = line.split("#", 1)[0].split()
                if not fields:
                    continue
                columns = requested_cols or list(range(1, len(fields) + 1))
                for column in columns:
                    if column < 1 or column > len(fields):
                        print(
                            f"ERROR [{case_name}]: '{filename}' line {line_number} "
                            f"has no requested column {column}"
                        )
                        sys.exit(1)
                    try:
                        value = float(fields[column - 1].replace("D", "E").replace("d", "e"))
                    except ValueError:
                        print(
                            f"ERROR [{case_name}]: '{filename}' line {line_number} "
                            f"column {column} is not numeric"
                        )
                        sys.exit(1)
                    if not math.isfinite(value):
                        print(
                            f"ERROR [{case_name}]: '{filename}' line {line_number} "
                            f"column {column} is non-finite ({fields[column - 1]!r})"
                        )
                        sys.exit(1)


def check_log_assertions(workdir: str, case_name: str, assertions: list[dict]) -> None:
    """Check explicit numerical invariants captured from a log file."""
    contents: dict[str, str] = {}
    for assertion in assertions:
        filename = assertion["file"]
        if filename not in contents:
            with open(os.path.join(workdir, filename)) as fh:
                contents[filename] = fh.read()
        match = re.search(assertion["pattern"], contents[filename])
        if match is None:
            print(
                f"ERROR [{case_name}]: log assertion pattern did not match "
                f"{filename}: {assertion['pattern']}"
            )
            sys.exit(1)
        try:
            actual = float(match.group(1).replace("D", "E").replace("d", "e"))
        except ValueError:
            print(f"ERROR [{case_name}]: log assertion value is not numeric in {filename}")
            sys.exit(1)
        expected = float(assertion["value"])
        tolerance = float(assertion.get("abs_tol", 0.0))
        if not math.isfinite(actual) or abs(actual - expected) > tolerance:
            print(
                f"ERROR [{case_name}]: {filename} invariant failed: "
                f"actual={actual!r}, expected={expected!r}, tolerance={tolerance!r}"
            )
            sys.exit(1)


# ---------------------------------------------------------------------------
# Value extraction
# ---------------------------------------------------------------------------

def _nml_flat(filepath: str) -> dict:
    """Read a namelist file and return a flat key->value dict (all sections merged)."""
    nml = f90nml.read(filepath)
    flat: dict = {}
    for section in nml.values():
        if isinstance(section, dict):
            flat.update(section)
    return flat


def extract_nml_values(workdir: str, nml_check: dict) -> dict:
    """Extract scalars and array elements from a namelist output file."""
    filepath = os.path.join(workdir, nml_check["file"])
    flat = _nml_flat(filepath)
    result: dict = {}

    for key in nml_check.get("scalars", []):
        result[key] = float(flat[key])

    for key, indices in nml_check.get("arrays", {}).items():
        arr = flat[key]  # f90nml returns a Python list (0-based)
        result[key] = {str(i): float(arr[i - 1]) for i in indices}

    return result


def extract_text_values(workdir: str, text_check: dict) -> dict:
    """Extract values at specified (row, col) positions from a space-separated file."""
    filepath = os.path.join(workdir, text_check["file"])
    rows = text_check["rows"]
    cols = text_check["cols"]
    row_set = set(rows)

    result: dict = {}
    with open(filepath) as fh:
        for lineno, line in enumerate(fh, start=1):
            if lineno in row_set:
                parts = line.split()
                result[str(lineno)] = {str(c): float(parts[c - 1]) for c in cols}

    return result


def extract_log_values(workdir: str, log_check: dict) -> dict:
    """Extract named scalar values from the first match of each log regex."""
    filepath = os.path.join(workdir, log_check["file"])
    with open(filepath) as fh:
        content = fh.read()

    result: dict = {}
    for key, pattern in log_check.get("patterns", {}).items():
        match = re.search(pattern, content)
        if match is None:
            raise ValueError(f"Log pattern '{key}' did not match {filepath}: {pattern}")
        result[key] = float(match.group(1).replace("D", "E").replace("d", "e"))
    return result


def build_ref_data(workdir: str, checks: dict) -> dict:
    """Build a complete reference data dict from the workdir outputs."""
    ref: dict = {}

    if "nml" in checks:
        ref["nml"] = {}
        for nml_check in checks["nml"]:
            ref["nml"][nml_check["file"]] = extract_nml_values(workdir, nml_check)

    if "text" in checks:
        ref["text"] = {}
        for text_check in checks["text"]:
            ref["text"][text_check["file"]] = extract_text_values(workdir, text_check)

    if "log" in checks:
        ref["log"] = {}
        for log_check in checks["log"]:
            ref["log"][log_check["file"]] = extract_log_values(workdir, log_check)

    return ref


# ---------------------------------------------------------------------------
# Reference comparison
# ---------------------------------------------------------------------------

def _check_value(
    failures: list,
    label: str,
    run_v: float | None,
    ref_v: float,
    abs_tol: float,
    rel_tol: float,
) -> None:
    if run_v is None:
        failures.append(f"  {label}: missing in run output")
        return
    if not math.isfinite(run_v):
        failures.append(f"  {label}: non-finite run value ({run_v!r})")
        return
    if not math.isfinite(ref_v):
        failures.append(f"  {label}: non-finite reference value ({ref_v!r})")
        return
    abs_diff = abs(run_v - ref_v)
    scale = max(abs(ref_v), 1.0)
    rel_diff = abs_diff / scale
    if abs_diff > abs_tol and rel_diff > rel_tol:
        failures.append(
            f"  {label}  run={run_v:.12e}  ref={ref_v:.12e}"
            f"  |diff|={abs_diff:.3e}  rel={rel_diff:.3e}"
        )


def _text_abs_tol(filename: str, abs_tol: float) -> float:
    """Apply the DOS text-output resolution floor without weakening other checks."""
    if "dos" in filename.lower():
        return max(abs_tol, DOS_TEXT_ABS_TOL)
    return abs_tol


def _check_tolerances(
    check: dict | None,
    key: str | None,
    default_abs_tol: float,
    default_rel_tol: float,
) -> tuple[float, float]:
    """Resolve a check-level or key-level tolerance override.

    A check may specify ``abs_tol``/``rel_tol`` for all values, or a
    ``tolerances`` mapping for selected values.  For text checks the key is
    normally the column number; a ``row:column`` entry takes precedence when
    present.  This keeps backend-sensitive columns (for example an EF-shifted
    energy grid) from forcing a loose tolerance on unrelated observables.
    """
    if not check:
        return default_abs_tol, default_rel_tol

    spec = check
    overrides = check.get("tolerances", {})
    if key is not None and isinstance(overrides, dict):
        spec = overrides.get(key, check)

    if not isinstance(spec, dict):
        spec = check
    return (
        float(spec.get("abs_tol", default_abs_tol)),
        float(spec.get("rel_tol", default_rel_tol)),
    )


def _check_for_file(checks: dict, kind: str, filename: str) -> dict | None:
    """Return the manifest check describing *filename*, if one exists."""
    for check in checks.get(kind, []):
        if check.get("file") == filename:
            return check
    return None


def _reference_variant() -> str | None:
    """Return the native reference variant for this runner.

    A manual override is useful when reproducing a CI comparison locally.
    Unknown platforms use the unqualified ``ref.json`` baseline.
    """
    override = os.environ.get("RSLMTO_REFERENCE_VARIANT", "").strip()
    if override:
        return None if override.lower() in {"default", "none"} else override

    machine = platform.machine().lower()
    if sys.platform == "darwin":
        arch = "arm64" if machine in {"arm64", "aarch64"} else machine
        return f"macos-{arch}"
    if sys.platform.startswith("linux"):
        arch = "x86_64" if machine in {"x86_64", "amd64"} else machine
        return f"linux-{arch}"
    return None


def _reference_path(ref_dir: str, case_name: str) -> str:
    """Select a native reference when one exists, otherwise the baseline."""
    case_dir = os.path.join(ref_dir, case_name)
    variant = _reference_variant()
    if variant:
        native_path = os.path.join(case_dir, f"ref.{variant}.json")
        if os.path.exists(native_path):
            return native_path
    return os.path.join(case_dir, "ref.json")


def compare_ref(
    workdir: str,
    case_name: str,
    ref_dir: str,
    checks: dict,
    abs_tol: float,
    rel_tol: float,
) -> None:
    if not checks:
        print(f"PASS [{case_name}]: no checks defined (smoke only)")
        return

    ref_path = _reference_path(ref_dir, case_name)
    if not os.path.exists(ref_path):
        print(f"PASS [{case_name}]: no reference found, smoke only (run --gen-ref to create)")
        return

    with open(ref_path) as fh:
        ref_data = json.load(fh)

    run_data = build_ref_data(workdir, checks)
    failures: list[str] = []
    n_checked = 0

    # NML comparisons
    for filename, ref_vals in ref_data.get("nml", {}).items():
        run_vals = run_data.get("nml", {}).get(filename, {})
        nml_check = _check_for_file(checks, "nml", filename)
        for key, ref_val in ref_vals.items():
            if isinstance(ref_val, dict):
                run_val = run_vals.get(key, {})
                for idx, ref_v in ref_val.items():
                    run_v = run_val.get(idx) if isinstance(run_val, dict) else None
                    check_abs_tol, check_rel_tol = _check_tolerances(
                        nml_check, f"{key}[{idx}]", abs_tol, rel_tol
                    )
                    _check_value(
                        failures, f"{filename}:{key}[{idx}]", run_v, ref_v,
                        check_abs_tol, check_rel_tol,
                    )
                    n_checked += 1
            else:
                check_abs_tol, check_rel_tol = _check_tolerances(
                    nml_check, key, abs_tol, rel_tol
                )
                _check_value(
                    failures, f"{filename}:{key}", run_vals.get(key), ref_val,
                    check_abs_tol, check_rel_tol,
                )
                n_checked += 1

    # Text comparisons
    for filename, ref_rows in ref_data.get("text", {}).items():
        run_rows = run_data.get("text", {}).get(filename, {})
        text_check = _check_for_file(checks, "text", filename)
        for row_str, ref_cols in ref_rows.items():
            run_cols = run_rows.get(row_str, {})
            for col_str, ref_v in ref_cols.items():
                check_abs_tol, check_rel_tol = _check_tolerances(
                    text_check, f"{row_str}:{col_str}", abs_tol, rel_tol
                )
                if text_check:
                    check_abs_tol, check_rel_tol = _check_tolerances(
                        text_check, col_str, check_abs_tol, check_rel_tol
                    )
                text_abs_tol = _text_abs_tol(filename, check_abs_tol)
                _check_value(
                    failures,
                    f"{filename}:row{row_str}:col{col_str}",
                    run_cols.get(col_str),
                    ref_v,
                    text_abs_tol,
                    check_rel_tol,
                )
                n_checked += 1

    # Log comparisons
    for filename, ref_vals in ref_data.get("log", {}).items():
        run_vals = run_data.get("log", {}).get(filename, {})
        log_check = _check_for_file(checks, "log", filename)
        for key, ref_v in ref_vals.items():
            check_abs_tol, check_rel_tol = _check_tolerances(
                log_check, key, abs_tol, rel_tol
            )
            _check_value(
                failures, f"{filename}:{key}", run_vals.get(key), ref_v,
                check_abs_tol, check_rel_tol,
            )
            n_checked += 1

    if failures:
        print(f"FAIL [{case_name}]: {len(failures)} of {n_checked} value(s) out of tolerance")
        for msg in failures:
            print(msg)
        sys.exit(1)

    print(f"PASS [{case_name}]: reference OK ({n_checked} values checked)")


# ---------------------------------------------------------------------------
# Reference generation
# ---------------------------------------------------------------------------

def save_ref(workdir: str, case_name: str, ref_dir: str, case: dict) -> None:
    checks = case.get("checks", {})
    dest_dir = os.path.join(ref_dir, case_name)
    os.makedirs(dest_dir, exist_ok=True)

    ref_data = build_ref_data(workdir, checks)
    with open(os.path.join(dest_dir, "ref.json"), "w") as fh:
        json.dump(ref_data, fh, indent=2)
    meta = dict(case)
    provenance = os.environ.get("RSLMTO_REFERENCE_PROVENANCE")
    if provenance:
        try:
            meta["reference_provenance"] = json.loads(provenance)
        except json.JSONDecodeError as exc:
            raise ValueError("RSLMTO_REFERENCE_PROVENANCE must contain valid JSON") from exc
    with open(os.path.join(dest_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    n = sum(
        (len(v) if isinstance(v, dict) else 1)
        for nml_vals in ref_data.get("nml", {}).values()
        for v in nml_vals.values()
    ) + sum(
        len(cols)
        for rows in ref_data.get("text", {}).values()
        for cols in rows.values()
    ) + sum(
        len(values)
        for values in ref_data.get("log", {}).values()
    )
    print(f"REF  [{case_name}]: {n} values saved to {dest_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Run one RSLMTO example test case.")
    parser.add_argument("--binary", required=True, help="Path to rslmto.x binary")
    parser.add_argument("--cases-json", required=True, help="Path to cases.json manifest")
    parser.add_argument("--case-name", required=True, help="Name of the case to run")
    parser.add_argument("--scratch-root", required=True, help="Root directory for scratch runs")
    parser.add_argument("--scratch-name", default=None,
                        help="Scratch subdirectory name (defaults to --case-name); "
                             "use to avoid collisions between launch-mode variants of one case")
    parser.add_argument("--ref-name", default=None,
                        help="Reference subdirectory name (defaults to --case-name); "
                             "set to the base case name so a serial/mpi pair shares one reference")
    parser.add_argument("--force-serial", action="store_true",
                        help="Force mpi_procs=1 regardless of case mpi_procs or --mpi-enabled")
    parser.add_argument("--mpi-procs", type=int, default=None,
                        help="Override the case MPI rank count (requires --mpi-enabled)")

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--compare-ref", metavar="REF_DIR", help="Compare output against stored reference")
    mode.add_argument("--gen-ref", metavar="REF_DIR", help="Save output as new reference")

    parser.add_argument("--abs-tol", type=float, default=1e-6)
    parser.add_argument("--rel-tol", type=float, default=1e-6)
    parser.add_argument("--mpi-enabled", action="store_true",
                        help="Honor case mpi_procs entries. Without this, cases run serially.")
    args = parser.parse_args()

    case = load_case(args.cases_json, args.case_name)
    binary = os.path.abspath(args.binary)
    cases_dir = os.path.join(os.path.dirname(os.path.abspath(args.cases_json)), "cases")
    case_dir = os.path.join(cases_dir, case["case"])
    scratch_name = args.scratch_name or args.case_name
    ref_name = args.ref_name or args.case_name
    workdir = os.path.join(args.scratch_root, scratch_name)

    setup_scratch(case_dir, workdir, preserve_outputs=case.get("preserve_outputs", False))
    patch_input_nml(workdir, case)
    serial_omp = case.get("omp_threads", None)
    if args.mpi_procs is not None and args.mpi_procs < 1:
        parser.error("--mpi-procs must be a positive integer")
    if args.mpi_procs is not None and not args.mpi_enabled:
        parser.error("--mpi-procs requires --mpi-enabled")
    mpi_procs = 1 if args.force_serial else (
        args.mpi_procs if args.mpi_enabled and args.mpi_procs is not None else
        (case.get("mpi_procs", 1) if args.mpi_enabled else 1)
    )
    run_binary(binary, workdir, mpi_procs, serial_omp)
    check_log(workdir, args.case_name)
    check_output_files(workdir, args.case_name, case.get("outputs", []))
    check_text_dimensions(workdir, args.case_name, case.get("checks", {}).get("dimensions", []))
    check_text_finite(workdir, args.case_name, case.get("checks", {}).get("finite", []))
    check_log_assertions(workdir, args.case_name, case.get("checks", {}).get("assertions", []))

    if args.compare_ref:
        abs_tol = case.get("abs_tol", args.abs_tol)
        rel_tol = case.get("rel_tol", args.rel_tol)
        compare_ref(
            workdir, ref_name, args.compare_ref,
            case.get("checks", {}), abs_tol, rel_tol,
        )
    elif args.gen_ref:
        save_ref(workdir, ref_name, args.gen_ref, case)
    else:
        print(f"PASS [{args.case_name}]")


if __name__ == "__main__":
    main()
