#!/usr/bin/env python3
"""Focused LDA+U energy-routing, restart, scope, and MPI receipt runner.

The positive path starts from the accepted KS U=2 eV,J=0 endpoint and omits
``hubbard_radial_moments`` to exercise the native-KS default. It checks the
logged route identity from the same production helper and records a direct
P_OUT re-evaluation from the lightweight trace. The serial continuation and
two-rank run each start from their documented restart state. One-step scope
controls cover explicit spectral, no-U default, and accepted RS-bulk auto
routing. The one-step controls check routing identities, not convergence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from decimal import Decimal
from pathlib import Path

ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_hashes() -> dict[str, str]:
    root = Path(__file__).resolve().parents[3]
    files = ("source/hamiltonian.f90", "source/hamiltonian_build.f90",
             "source/hamiltonian_hubbard.f90", "source/include_codes/namelists/hamiltonian.f90",
             "source/self.f90")
    return {name: sha(root / name) for name in files}


def assert_route_eligibility_source_contract() -> dict[str, bool]:
    root = Path(__file__).resolve().parents[3]
    source = re.sub(r"\s+", " ", (root / "source/self.f90").read_text())
    predicates = {
        "calctype_B": "this%control%calctype == 'B'" in source,
        "constraints_disabled": ".not. this%control%constraints_enable" in source,
        "clean_bulk_no_impurity": "this%lattice%nbulk == 0 .and. this%lattice%nrec == this%lattice%ntype" in source,
        "single_spin_d_shell": "this%control%nsp == 1 .and. this%control%lmax == 2" in source,
        "fixed_liechtenstein": "this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein'" in source,
    }
    failed = [name for name, present in predicates.items() if not present]
    if failed:
        raise AssertionError("local-routing eligibility source contract missing: " + ", ".join(failed))
    return predicates


def set_group(text: str, group: str, values: dict[str, str | None]) -> str:
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines)
                  if line.strip().lower() == f"&{group}"), None)
    if start is None:
        lines.extend([f"&{group}", "/"])
        start = len(lines) - 2
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == "/"), None)
    if end is None:
        raise ValueError(f"unterminated &{group} namelist")
    pending = dict(values)
    kept = []
    for line in lines[start + 1:end]:
        stripped = line.lstrip()
        key = stripped.split("=", 1)[0].strip().lower() if "=" in stripped and not stripped.startswith("!") else ""
        if key not in pending:
            kept.append(line)
            continue
        value = pending.pop(key)
        if value is not None:
            kept.append(f"{key} = {value}")
    kept.extend(f"{key} = {value}" for key, value in pending.items() if value is not None)
    lines[start + 1:end] = kept
    return "\n".join(lines) + "\n"


def configure(base: str, *, backend: str, u: float, j: float,
              mode: str | None, test_lambda: float = 0.0,
              explicit_uj: bool = True) -> str:
    deck = set_group(base, "self", {"use_kspace": ".true." if backend == "ks" else ".false.",
                                     "nstep": "1"})
    ham: dict[str, str | None] = {
        "hubbard_u_general(1,3)": f"{u:.16g}" if explicit_uj else None,
        "hubbard_j_general(1,3)": f"{j:.16g}" if explicit_uj else None,
        "hubbard_u_potential_form": "'liechtenstein'",
        "hubbard_radial_moments": f"'{mode}'" if mode else None,
        "hubbard_test_operator": "'w1'",
        "hubbard_test_lambda": f"{test_lambda:.16g}",
        "hubbard_test_spin": "0",
        "hubbard_log_level": "'detailed'",
        "fixed_u_test_tighten_r05": ".false.",
    }
    return set_group(deck, "hamiltonian", ham)


def kv_fields(log: str) -> tuple[dict[str, str], str]:
    rows = re.findall(r"HUBBARD_MOMENT_ENERGY_TOTAL ([^\r\n]*)", log)
    if not rows:
        raise ValueError("missing HUBBARD_MOMENT_ENERGY_TOTAL row")
    fields = dict(re.findall(rf"([A-Za-z0-9_]+)=\s*({FLOAT}|UNAVAILABLE(?:_[A-Z0-9_]+)?)", rows[-1]))
    return fields, rows[-1]


def site_base_energy(log: str) -> tuple[Decimal, str]:
    rows = re.findall(r"HUBBARD_MOMENT_ENERGY_SITE ([^\r\n]*)", log)
    if len(rows) != 1:
        raise ValueError(f"clean single-site test expected one energy-site row, observed {len(rows)}")
    values = dict(re.findall(rf"([A-Za-z0-9_]+)=\s*({FLOAT}|UNAVAILABLE(?:_[A-Z0-9_]+)?)", rows[0]))
    if "ETOT_selected" not in values:
        raise ValueError("missing uncorrected per-site ETOT_selected")
    return Decimal(values["ETOT_selected"].replace("D", "E").replace("d", "e")), rows[0]


def dec(field: dict[str, str], key: str) -> Decimal:
    if key not in field or field[key].startswith("UNAVAILABLE"):
        raise ValueError(f"missing numeric log field: {key}")
    return Decimal(field[key].replace("D", "E").replace("d", "e"))


def comparison_tolerance(*values: Decimal) -> float:
    """Allow 16 binary64 ULPs at the largest compared magnitude.

    This covers the 16-digit E-format diagnostics plus ordinary accumulation
    roundoff. At the ~2.5e3 Ry total-energy scale it is below 8e-12 Ry, over
    seven orders smaller than the smallest registered ASA allowance (~5.8e-4
    Ry/eV). It is only an implementation/reproducibility check, not a physics
    or convergence tolerance.
    """
    scale = max(1.0, *(abs(float(x)) for x in values))
    return 16.0 * math.ulp(scale)


def check_sum(actual: Decimal, a: Decimal, b: Decimal) -> dict:
    residual = actual - a - b
    tolerance = comparison_tolerance(actual, a, b)
    return {"actual": str(actual), "sum_operands": str(a + b),
            "residual": str(residual), "tolerance_ry": tolerance,
            "pass": abs(float(residual)) <= tolerance}


def check_equal(actual: Decimal, expected: Decimal, label: str) -> dict:
    residual = actual - expected
    tolerance = comparison_tolerance(actual, expected)
    return {"label": label, "actual": str(actual), "expected": str(expected),
            "residual": str(residual), "tolerance_ry": tolerance,
            "pass": abs(float(residual)) <= tolerance}


def extract_p(records: Path, out: Path) -> int:
    matrix: dict[tuple[int, int, int], tuple[str, str]] = {}
    with records.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if (row.get("kind") != "test_density_matrix" or row.get("stage") != "P_OUT" or
                    row.get("site") != "0" or row.get("l") != "2"):
                continue
            spin0, i, j = int(row["spin"]), int(row["row"]), int(row["col"])
            if spin0 in (0, 1) and 0 <= i < 5 and 0 <= j < 5:
                matrix[(spin0 + 1, i, j)] = (row["re"].strip(), row["im"].strip())
    if len(matrix) != 50:
        raise ValueError(f"expected 50 direct P_OUT d-shell rows, got {len(matrix)}")
    with out.open("w") as stream:
        for (spin, i, j), (re_value, im_value) in sorted(matrix.items()):
            stream.write(f"{spin} {i + 1} {j + 1} {re_value} {im_value}\n")
    return len(matrix)


def compare_logged_decimal(a: str, b: str) -> dict:
    da = Decimal(a.replace("D", "E").replace("d", "e"))
    db = Decimal(b.replace("D", "E").replace("d", "e"))
    delta = da - db
    tolerance = comparison_tolerance(da, db)
    return {"serial": str(da), "mpi": str(db), "difference": str(delta),
            "tolerance": tolerance, "pass": abs(float(delta)) <= tolerance}


def compare_matrix_files(serial_file: Path, mpi_file: Path) -> dict:
    def read(path):
        values = {}
        with path.open() as stream:
            for line in stream:
                spin, i, j, re_value, im_value = line.split()
                values[(int(spin), int(i), int(j))] = (re_value, im_value)
        return values
    left, right = read(serial_file), read(mpi_file)
    if left.keys() != right.keys():
        return {"pass": False, "reason": "P_OUT matrix entry sets differ"}
    comparisons = []
    for index in sorted(left):
        for component in range(2):
            comparisons.append(compare_logged_decimal(left[index][component], right[index][component]))
    failed = [x for x in comparisons if not x["pass"]]
    max_difference = max((abs(Decimal(x["difference"])) for x in comparisons), default=Decimal(0))
    return {"pass": not failed, "entry_component_count": len(comparisons),
            "max_abs_difference": str(max_difference),
            "failed_components": len(failed),
            "comparison_basis": "16 binary64 ULPs at each component's comparison scale"}


def run_one(name: str, args, *, binary: Path, launcher: list[str], base_input: Path,
            seed: Path, backend: str, u: float, j: float, mode: str | None,
            expected_mode: str, restart_from: Path | None = None,
            test_lambda: float = 0.0, explicit_uj: bool = True,
            prepared_input: Path | None = None) -> dict:
    case = args.outdir / name
    case.mkdir(parents=False, exist_ok=False)
    input_path = case / "input.nml"
    if prepared_input:
        shutil.copyfile(prepared_input, input_path)
    else:
        deck = configure(base_input.read_text(), backend=backend, u=u, j=j, mode=mode,
                         test_lambda=test_lambda, explicit_uj=explicit_uj)
        input_path.write_text(deck)
    restart = case / "Fe_out.nml"
    shutil.copyfile(restart_from or seed, restart)
    for asset in ASSETS:
        shutil.copyfile(args.assets / asset, case / asset)
    input_hash_before = sha(input_path)
    restart_hash_before = sha(restart)
    restart_source_hash_before = sha(restart_from or seed)
    asset_hashes_before = {asset: sha(args.assets / asset) for asset in ASSETS}
    trace = case / "records.csv"
    log_path = case / "stdout.log"
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
                "RSLMTO_FIXED_U_TRACE": str(trace), "RSLMTO_FIXED_U_TRACE_RUN_ID": name,
                "RSLMTO_FIXED_U_TRACE_ITERATION": "0", "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
                "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1"})
    command = launcher + [str(binary)]
    with log_path.open("wb") as log:
        proc = subprocess.run(command, cwd=case, env=env, stdin=subprocess.DEVNULL,
                              stdout=log, stderr=subprocess.STDOUT, check=False)
    text = log_path.read_text(errors="replace")
    if proc.returncode != 0:
        raise RuntimeError(f"{name}: executable returned {proc.returncode}; inspect {log_path}")
    ks_marker = "run_recursion: use_kspace=.true., skipping recursion solver stage" in text
    if backend == "ks" and not ks_marker:
        raise AssertionError(f"{name}: input requested KS but runtime did not confirm reciprocal branch")
    if backend == "rs" and ks_marker:
        raise AssertionError(f"{name}: RS scope control unexpectedly entered reciprocal KS branch")
    fields, raw_energy = kv_fields(text)
    base_site, raw_site = site_base_energy(text)
    if backend == "ks":
        entropy = dec(fields, "entropy")
        ts = dec(fields, "TS")
        if not entropy.is_finite() or not ts.is_finite():
            raise AssertionError(f"{name}: KS entropy/TS is non-finite")
    mode_rows = re.findall(r"HUBBARD_RUN ([^\r\n]*)", text)
    mode_fields = dict(re.findall(r"([A-Za-z0-9_]+)=([^\s]+)", mode_rows[-1])) if mode_rows else {}
    actual_mode = mode_fields.get("radial_mode")
    if actual_mode not in {"local", "spectral"}:
        raise ValueError(f"{name}: missing/invalid HUBBARD_RUN radial_mode field")
    adopted = mode_fields.get("energy_routing", "") == "adopted"
    if actual_mode != expected_mode:
        raise AssertionError(f"{name}: expected mode {expected_mode}, observed {actual_mode}")
    expected_adopted = expected_mode == "local" and test_lambda == 0.0
    if adopted != expected_adopted:
        raise AssertionError(f"{name}: unexpected local-energy adoption marker {mode_fields}")
    physical_energy = dec(fields, "physical_total_energy")
    hubbard_logged = dec(fields, "E_Hub_P_OUT")
    route = check_sum(physical_energy, base_site,
                      hubbard_logged if adopted else Decimal(0))
    base_field_check = check_sum(dec(fields, "base_physical_total_energy"), base_site, Decimal(0))
    route["base_field_matches_unmodified_site_ETOT"] = base_field_check
    route["pass"] = route["pass"] and base_field_check["pass"]
    candidate_checks = {}
    ts = dec(fields, "TS") if backend == "ks" else Decimal(0)
    if expected_mode == "local" and adopted:
        candidate_checks["candidate_b"] = check_equal(
            dec(fields, "candidate_b"), physical_energy, "local candidate_b equals adopted physical total")
        if backend == "ks":
            candidate_checks["candidate_b_F"] = check_equal(
                dec(fields, "candidate_b_F"), physical_energy - ts,
                "local candidate_b_F equals adopted physical total minus TS")
    elif expected_mode == "spectral":
        e_hub_diag = dec(fields, "E_Hub")
        t_diag = dec(fields, "t_in_out")
        candidate_checks["candidate_a"] = check_equal(
            dec(fields, "candidate_a"), physical_energy - t_diag + e_hub_diag,
            "spectral candidate_a equals base ETOT - t + EHub")
        candidate_checks["candidate_c"] = check_equal(
            dec(fields, "candidate_c"), physical_energy + e_hub_diag,
            "spectral candidate_c equals base ETOT + EHub")
        if backend == "ks":
            candidate_checks["candidate_a_F"] = check_equal(
                dec(fields, "candidate_a_F"), physical_energy - t_diag + e_hub_diag - ts,
                "spectral candidate_a_F subtracts finite TS")
            candidate_checks["candidate_c_F"] = check_equal(
                dec(fields, "candidate_c_F"), physical_energy + e_hub_diag - ts,
                "spectral candidate_c_F subtracts finite TS")
    else:
        candidate_checks["nonadopted_EHub"] = check_equal(
            hubbard_logged, Decimal(0), "nonadopted local test field leaves production EHub unapplied")
        candidate_checks["nonadopted_physical_energy"] = check_equal(
            physical_energy, base_site,
            "nonzero external test field is excluded from adopted physical energy")
    route["candidate_energy_identities"] = candidate_checks
    route["pass"] = route["pass"] and all(item["pass"] for item in candidate_checks.values())
    if not route["pass"]:
        raise AssertionError(f"{name}: physical energy routing identity failed: {route}")
    matrix_file = case / "P_OUT_d.txt"
    matrix_rows = extract_p(trace, matrix_file) if backend == "ks" else 0
    helper_result = None
    helper_delta = None
    if adopted and backend == "ks":
        if not args.energy_helper:
            raise ValueError("--energy-helper is required to verify adopted local production EHub")
        helper = subprocess.run([str(args.energy_helper), str(matrix_file), str(u), str(j)],
                                check=True, text=True, capture_output=True)
        helper_result = helper.stdout.strip()
        match = re.search(rf"(?:E_?Hub_Ry|EHub_Ucenter_Ry)=\s*({FLOAT})", helper_result)
        if not match:
            raise ValueError(f"{name}: direct-P helper did not return E_Hub_Ry")
        helper_delta = str(Decimal(match.group(1).replace("D", "E")) - hubbard_logged)
        tolerance = comparison_tolerance(Decimal(match.group(1)), hubbard_logged)
        if abs(float(Decimal(helper_delta))) > tolerance:
            raise AssertionError(f"{name}: direct-P production helper disagrees with logged EHub: {helper_delta}")
    runtime_hashes = {p.name: sha(p) for p in case.iterdir() if p.is_file() and
                      p.name in {"Fe_out.nml", "Fe_scf.nml", "records.csv", "stdout.log"}}
    receipt = {
        "case": name, "backend": backend, "U_eV": u, "J_eV": j,
        "fixed_UJ_fields_explicit": explicit_uj,
        "prepared_input_source": str(prepared_input) if prepared_input else None,
        "requested_mode": mode or "auto_default", "expected_resolved_mode": expected_mode,
        "test_lambda_Ry": test_lambda,
        "observed_resolved_mode": actual_mode, "mode_runtime_fields": mode_fields,
        "local_energy_adopted": adopted, "expected_local_energy_adopted": expected_adopted,
        "returncode": proc.returncode,
        "native_kspace_runtime_marker": ks_marker,
        "KS_entropy_finite": backend != "ks" or entropy.is_finite(),
        "KS_TS_finite": backend != "ks" or ts.is_finite(),
        "binary_sha256": sha(binary), "input_sha256_pre_run": input_hash_before,
        "restart_sha256_pre_run": restart_hash_before,
        "assets_sha256_pre_run": asset_hashes_before,
        "restart_source_sha256": restart_source_hash_before,
        "route_fields_raw": raw_energy, "site_energy_fields_raw": raw_site,
        "unmodified_site_ETOT_Ry": str(base_site), "route_identity": route,
        "P_OUT_rows": matrix_rows, "P_OUT_sha256": sha(matrix_file) if matrix_file.exists() else None,
        "direct_P_energy_helper": helper_result, "helper_minus_logged_E_Hub_Ry": helper_delta,
        "runtime_output_sha256": runtime_hashes,
        "log_path": str(log_path), "trace_path": str(trace),
    }
    (case / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--serial-binary", type=Path, required=True)
    p.add_argument("--mpi-binary", type=Path, required=True)
    p.add_argument("--mpi-launcher", default="mpirun")
    p.add_argument("--mpi-ranks", type=int, default=2)
    p.add_argument("--accepted-input", type=Path, required=True)
    p.add_argument("--accepted-u2-restart", type=Path, required=True)
    p.add_argument("--spectral-u2-restart", type=Path, required=True)
    p.add_argument("--no-u-restart", type=Path, required=True)
    p.add_argument("--rs-restart", type=Path, required=True)
    p.add_argument("--assets", type=Path, required=True)
    p.add_argument("--energy-helper", type=Path, required=True,
                   help="production fixed-P EHub evaluator; required for adopted local cases")
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    required = [args.serial_binary, args.mpi_binary, args.accepted_input,
                args.accepted_u2_restart, args.spectral_u2_restart,
                args.no_u_restart, args.rs_restart,
                args.energy_helper, *(args.assets / x for x in ASSETS)]
    missing = [str(f) for f in required if not f.is_file()]
    if missing:
        p.error("missing file(s): " + ", ".join(missing))
    if args.mpi_ranks < 2:
        p.error("--mpi-ranks must be at least 2")
    args.outdir.mkdir(parents=True, exist_ok=False)
    launcher = [args.mpi_launcher, "--oversubscribe", "-n", str(args.mpi_ranks)]
    rows = []
    first = run_one("ks_auto_local_serial", args, binary=args.serial_binary, launcher=[],
                    base_input=args.accepted_input, seed=args.accepted_u2_restart,
                    backend="ks", u=2.0, j=0.0, mode=None, expected_mode="local")
    rows.append(first)
    first_restart = args.outdir / "ks_auto_local_serial" / "Fe_out.nml"
    if not first_restart.is_file():
        raise FileNotFoundError("serial run did not produce restart Fe_out.nml")
    restart_row = run_one("ks_auto_local_restart", args, binary=args.serial_binary, launcher=[],
                          base_input=args.accepted_input, seed=args.accepted_u2_restart,
                          restart_from=first_restart, backend="ks", u=2.0, j=0.0,
                          mode=None, expected_mode="local")
    restart_row["restart_chain_matches_previous_output"] = (
        restart_row["restart_sha256_pre_run"] == first["runtime_output_sha256"].get("Fe_out.nml"))
    rows.append(restart_row)
    mpi_row = run_one("ks_auto_local_mpi", args, binary=args.mpi_binary, launcher=launcher,
                      base_input=args.accepted_input, seed=args.accepted_u2_restart,
                      backend="ks", u=2.0, j=0.0, mode=None, expected_mode="local")
    rows.append(mpi_row)
    rows.append(run_one("ks_explicit_spectral", args, binary=args.serial_binary, launcher=[],
                        base_input=args.accepted_input, seed=args.spectral_u2_restart,
                        backend="ks", u=2.0, j=0.0, mode="spectral", expected_mode="spectral"))
    rows.append(run_one("ks_explicit_local_nonzero_test_field", args,
                        binary=args.serial_binary, launcher=[], base_input=args.accepted_input,
                        seed=args.accepted_u2_restart, backend="ks", u=2.0, j=0.0,
                        mode="local", expected_mode="local", test_lambda=0.005))
    rows.append(run_one("ks_no_u_auto", args, binary=args.serial_binary, launcher=[],
                        base_input=args.accepted_input, seed=args.no_u_restart,
                        backend="ks", u=0.0, j=0.0, mode=None, expected_mode="spectral",
                        explicit_uj=False))
    rows.append(run_one("ks_explicit_zero_auto", args, binary=args.serial_binary, launcher=[],
                        base_input=args.accepted_input, seed=args.no_u_restart,
                        backend="ks", u=0.0, j=0.0, mode=None, expected_mode="spectral"))
    rows.append(run_one("rs_u2_auto", args, binary=args.serial_binary, launcher=[],
                        base_input=args.accepted_input, seed=args.rs_restart,
                        backend="rs", u=2.0, j=0.0, mode=None, expected_mode="local"))
    mpi_case = args.outdir / "ks_auto_local_mpi"
    serial_case = args.outdir / "ks_auto_local_serial"
    serial_fields, _ = kv_fields((serial_case / "stdout.log").read_text(errors="replace"))
    mpi_fields, _ = kv_fields((mpi_case / "stdout.log").read_text(errors="replace"))
    serialized_comparison = {
        key: compare_logged_decimal(serial_fields[key], mpi_fields[key])
        for key in ("physical_total_energy", "base_physical_total_energy", "E_Hub_P_OUT")
    }
    serialized_comparison["unmodified_site_ETOT"] = compare_logged_decimal(
        rows[0]["unmodified_site_ETOT_Ry"], mpi_row["unmodified_site_ETOT_Ry"])
    serialized_comparison["direct_P_OUT"] = compare_matrix_files(
        serial_case / "P_OUT_d.txt", mpi_case / "P_OUT_d.txt")
    serialized_comparison["same_input_hash"] = (
        rows[0]["input_sha256_pre_run"] == mpi_row["input_sha256_pre_run"])
    serialized_comparison["same_initial_restart_hash"] = (
        rows[0]["restart_sha256_pre_run"] == mpi_row["restart_sha256_pre_run"])
    serialized_comparison["exact_restart_output_match"] = (
        rows[0]["runtime_output_sha256"].get("Fe_out.nml") ==
        mpi_row["runtime_output_sha256"].get("Fe_out.nml"))
    serialized_comparison["numerical_fields_pass"] = all(
        serialized_comparison[key]["pass"] for key in
            ("physical_total_energy", "base_physical_total_energy", "E_Hub_P_OUT",
             "unmodified_site_ETOT", "direct_P_OUT"))
    omitted_zero = args.outdir / "ks_no_u_auto"
    explicit_zero = args.outdir / "ks_explicit_zero_auto"
    omitted_fields, _ = kv_fields((omitted_zero / "stdout.log").read_text(errors="replace"))
    explicit_fields, _ = kv_fields((explicit_zero / "stdout.log").read_text(errors="replace"))
    zero_noop = {
        key: compare_logged_decimal(omitted_fields[key], explicit_fields[key])
        for key in ("physical_total_energy", "base_physical_total_energy", "E_Hub_P_OUT")
    }
    zero_noop["direct_P_OUT"] = compare_matrix_files(
        omitted_zero / "P_OUT_d.txt", explicit_zero / "P_OUT_d.txt")
    zero_noop["all_equal_within_tolerance"] = all(
        v["pass"] for v in zero_noop.values())
    summary = {
        "schema": "a4_4_energy_route_repro_v1",
        "status": "PASS" if all(row["route_identity"]["pass"] for row in rows) and
            serialized_comparison["numerical_fields_pass"] and
            serialized_comparison["same_input_hash"] and
            serialized_comparison["same_initial_restart_hash"] and
            restart_row["restart_chain_matches_previous_output"] and
            zero_noop["all_equal_within_tolerance"] else "FAIL",
        "run_count": len(rows), "MPI_tool": subprocess.run([args.mpi_launcher, "--version"],
             text=True, capture_output=True, check=False).stdout.splitlines()[:2],
        "serial_binary_sha256": sha(args.serial_binary), "mpi_binary_sha256": sha(args.mpi_binary),
        "source_sha256": source_hashes(),
        "eligibility_source_contract": assert_route_eligibility_source_contract(),
        "runner_sha256": sha(Path(__file__).resolve()),
        "energy_helper_sha256": sha(args.energy_helper) if args.energy_helper else None,
        "runs": rows,
        "serial_mpi_reproducibility": serialized_comparison,
        "zero_UJ_noop_comparison": zero_noop,
        "serial_restart_chain": {
            "first_run_output_Fe_out_sha256": first["runtime_output_sha256"].get("Fe_out.nml"),
            "restart_run_input_Fe_out_sha256": restart_row["restart_sha256_pre_run"],
            "same_restart_bytes": restart_row["restart_chain_matches_previous_output"],
            "identity_rechecked_on_restart": restart_row["route_identity"]["pass"],
            "pass": restart_row["restart_chain_matches_previous_output"] and
                    restart_row["route_identity"]["pass"],
        },
        "comparison_note": "Route and serial/MPI comparisons use 16 binary64 ULPs at the compared magnitude, covering 16-digit E-format output and ordinary accumulation roundoff. The energy tolerance stays more than seven orders below the smallest registered ASA allowance; actual residuals and restart hashes are reported. No SCF convergence or physics tolerance is inferred.",
    }
    out = args.outdir / "summary.json"
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": summary["status"], "run_count": len(rows), "summary": str(out)}, sort_keys=True))
    if summary["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
