#!/usr/bin/env python3
"""Run one RC2 deck while archiving each complete rolling P_OUT trace.

The production writer resets its CSV each SCF iteration.  This wrapper samples
that file during the run and preserves iterations only after the direct d-P
matrix, all six ql records, and both W1 traces are complete.  It does not alter
solver source or accept partial CSV rows as evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import signal
import shutil
import subprocess
import time
from pathlib import Path

HEADER = "kind,run_id,lambda,iteration,stage,site,type,l,spin,k,band,row,col,order,weight,energy,scalar,re,im"
THREADS = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
DIFF_RE = re.compile(r"(?:Not converged!\s*)Diff=\s*([-+0-9.eEdD]+)")


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


# Freeze wrapper provenance at process start so later poll-rate updates cannot
# be attributed to already-running monitor processes.
RUNNER_SOURCE_PATH = Path(__file__).resolve()
RUNNER_SOURCE_SHA256 = sha(RUNNER_SOURCE_PATH)


def complete_records(raw: bytes) -> dict[int, list[dict[str, str]]]:
    try:
        rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    except (UnicodeDecodeError, csv.Error):
        return {}
    by_iteration: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        try:
            it = int(row["iteration"])
            if it < 1 or row["stage"] != "P_OUT":
                continue
            by_iteration.setdefault(it, []).append(row)
        except (KeyError, TypeError, ValueError):
            continue
    complete: dict[int, list[dict[str, str]]] = {}
    for it, subset in by_iteration.items():
        qls = [r for r in subset if r.get("kind") == "ql" and r.get("site") == "0"]
        matrices = [r for r in subset if r.get("kind") == "test_density_matrix" and r.get("site") == "0" and r.get("l") == "2"]
        wp = [r for r in subset if r.get("kind") == "test_trace_wp" and r.get("site") == "0"]
        try:
            qkeys = {(int(r["l"]), int(r["spin"])) for r in qls}
            mkeys = {(int(r["spin"]), int(r["row"]), int(r["col"])) for r in matrices}
            wkeys = {int(r["spin"]) for r in wp}
        except (KeyError, TypeError, ValueError):
            continue
        if len(qls) == 6 and qkeys == {(l, s) for l in range(3) for s in (1, 2)} and \
           len(matrices) == 50 and mkeys == {(s, i, j) for s in (0, 1) for i in range(5) for j in range(5)} and \
           len(wp) == 2 and wkeys == {0, 1}:
            complete[it] = subset
    return complete


def matrix_checks(rows: list[dict[str, str]]) -> dict:
    matrix = {}
    ql = {}
    wp = {}
    nonfinite = []
    malformed = []
    for row in rows:
        kind = row.get("kind")
        if kind == "test_density_matrix" and row.get("site") == "0" and row.get("l") == "2":
            try:
                key = (int(row["spin"]), int(row["row"]), int(row["col"]))
                real, imag = float(row["re"]), float(row["im"])
            except (KeyError, TypeError, ValueError) as exc:
                malformed.append({"kind": kind, "error": str(exc), "row": row})
                continue
            if not math.isfinite(real) or not math.isfinite(imag):
                nonfinite.append({"kind": kind, "key": key, "re": row["re"], "im": row["im"]})
            matrix[key] = complex(real, imag)
        elif kind == "ql" and row.get("site") == "0":
            try:
                val = float(row["scalar"])
            except (KeyError, TypeError, ValueError) as exc:
                malformed.append({"kind": kind, "error": str(exc), "row": row})
                continue
            if not math.isfinite(val):
                nonfinite.append({"kind": kind, "l": row["l"], "spin": row["spin"], "scalar": row["scalar"]})
            if row.get("l") == "2":
                ql[int(row["spin"]) - 1] = val
        elif kind == "test_trace_wp" and row.get("site") == "0":
            try:
                val = float(row["scalar"])
            except (KeyError, TypeError, ValueError) as exc:
                malformed.append({"kind": kind, "error": str(exc), "row": row})
                continue
            if not math.isfinite(val):
                nonfinite.append({"kind": kind, "spin": row["spin"], "scalar": row["scalar"]})
            wp[int(row["spin"])] = val
    if nonfinite or malformed:
        return {"pass_1e_11": False, "nonfinite_values": nonfinite, "malformed_values": malformed}
    herm = 0.0
    traces = {}
    for spin in (0, 1):
        traces[spin] = sum(matrix[(spin, i, i)].real for i in range(5))
        for i in range(5):
            for j in range(5):
                herm = max(herm, abs(matrix[(spin, i, j)] - matrix[(spin, j, i)].conjugate()))
    ql_errors = [abs(traces[s] - ql[s]) for s in (0, 1)]
    wp_errors = [abs(traces[s] - wp[s]) for s in (0, 1)]
    return {
        "max_hermiticity_residual": herm,
        "spin_trace_differences_vs_ql_l2": ql_errors,
        "spin_trace_differences_vs_direct_W1_trace": wp_errors,
        "pass_1e_11": herm <= 1e-11 and max(ql_errors + wp_errors) <= 1e-11,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=1.0)
    parser.add_argument("--stop-on-c1", action="store_true",
                        help="freeze and capture the first accepted C1 endpoint instead of running to cap")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    trace = run / "fixed_u_trace.csv"
    log = run / "console.log"
    samples = run / "trace_snapshots"
    samples.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(THREADS)
    env.update({
        "RSLMTO_FIXED_U_TRACE": str(trace),
        "RSLMTO_FIXED_U_TRACE_ITERATION": "0",
        "RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY": "1",
        "RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY": "1",
        "RSLMTO_FIXED_U_TRACE_QL_ONLY": "0",
        "RSLMTO_FIXED_U_LAMBDA": "0.0",
        "RSLMTO_FIXED_U_RUN_ID": run.name,
    })
    input_initial_sha256 = sha(run / "input.nml")
    seed_initial_sha256 = sha(run / "Fe_out.nml")
    captured: dict[int, bytes] = {}
    per_iteration_checks: dict[int, dict] = {}
    early_checkpoint = None
    start = time.time()
    with log.open("wb") as stdout:
        process = subprocess.Popen([str(args.binary.resolve())], cwd=run, env=env,
                                   stdin=subprocess.DEVNULL, stdout=stdout, stderr=subprocess.STDOUT)
        while process.poll() is None:
            try:
                raw = trace.read_bytes()
            except FileNotFoundError:
                time.sleep(args.poll_seconds)
                continue
            for iteration, rows in complete_records(raw).items():
                if iteration not in captured:
                    buf = io.StringIO(newline="")
                    writer = csv.DictWriter(buf, fieldnames=HEADER.split(","), lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
                    data = buf.getvalue().encode()
                    captured[iteration] = data
                    per_iteration_checks[iteration] = matrix_checks(rows)
                    (samples / f"P_OUT_iter{iteration:04d}.csv").write_bytes(data)
            if args.stop_on_c1 and captured:
                ids = sorted(captured)
                if len(ids) >= 20 and ids[-20:] == list(range(ids[-1] - 19, ids[-1] + 1)):
                    log_text = log.read_text(errors="replace") if log.exists() else ""
                    diffs = [float(x.replace("D", "E").replace("d", "e")) for x in DIFF_RE.findall(log_text)]
                    if len(diffs) == ids[-1] and all(x["pass_1e_11"] for x in per_iteration_checks.values()):
                        rows_by_it = {}
                        for iteration, data in captured.items():
                            rows_by_it[iteration] = list(csv.DictReader(io.StringIO(data.decode())))
                        d_up, d_down, moments = [], [], []
                        valid = True
                        for iteration in ids[-20:]:
                            qrows = [r for r in rows_by_it[iteration]
                                     if r.get("kind") == "ql" and r.get("site") == "0"]
                            try:
                                q = {(int(r["l"]), int(r["spin"])): float(r["scalar"]) for r in qrows}
                                valid = valid and set(q) == {(l, s) for l in range(3) for s in (1, 2)}
                                valid = valid and all(math.isfinite(v) for v in q.values())
                                d_up.append(q[(2, 1)]); d_down.append(q[(2, 2)])
                                moments.append(abs(sum(q[(l, 1)] - q[(l, 2)] for l in range(3))))
                            except (KeyError, TypeError, ValueError):
                                valid = False
                                break
                        spreads = (max(d_up) - min(d_up), max(d_down) - min(d_down), max(moments) - min(moments)) if valid else (math.inf,)*3
                        if valid and diffs[-1] <= 3e-8 and spreads[0] <= 5e-7 and spreads[1] <= 5e-7 and spreads[2] <= 1e-6:
                            # Stop only after the matching Diff marker is visible. The
                            # P_OUT trace and post-mix *_scf checkpoint are copied
                            # while the solver is frozen, then the owned process ends.
                            process.send_signal(signal.SIGSTOP)
                            time.sleep(0.05)
                            checkpoint = run / "C1_ACCEPTED_LIVE_CHECKPOINT"
                            checkpoint.mkdir()
                            checkpoint_files = ("Fe_scf.nml", "Fe_out.nml", "Fe_dos.out", "Fe_orbital_dos.out", "totaldos.out", "report.out")
                            for name in checkpoint_files:
                                if (run / name).is_file():
                                    shutil.copyfile(run / name, checkpoint / name)
                            checkpoint_trace = checkpoint / "fixed_u_trace_history.csv"
                            with checkpoint_trace.open("wb") as stream:
                                stream.write((HEADER + "\n").encode())
                                for it in ids:
                                    stream.writelines(captured[it].splitlines(keepends=True)[1:])
                            early_checkpoint = {
                                "iteration": ids[-1], "window_iterations": ids[-20:], "final_diff": diffs[-1],
                                "spreads": {"d_up": spreads[0], "d_down": spreads[1], "sphere_moment_muB": spreads[2]},
                                "files_sha256": {p.name: sha(p) for p in checkpoint.iterdir() if p.is_file()},
                                "checkpoint_dir": str(checkpoint), "trace_sha256": sha(checkpoint_trace),
                                "stop_method": "SIGSTOP, copy frozen solver state, SIGTERM owned process",
                            }
                            process.send_signal(signal.SIGTERM)
                            break
            time.sleep(args.poll_seconds)
    returncode = process.wait()
    if trace.exists():
        raw = trace.read_bytes()
        for iteration, rows in complete_records(raw).items():
            if iteration not in captured:
                buf = io.StringIO(newline="")
                writer = csv.DictWriter(buf, fieldnames=HEADER.split(","), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                data = buf.getvalue().encode()
                captured[iteration] = data
                per_iteration_checks[iteration] = matrix_checks(rows)
                (samples / f"P_OUT_iter{iteration:04d}.csv").write_bytes(data)
    aggregate = run / "fixed_u_trace_history.csv"
    with aggregate.open("wb") as stream:
        stream.write((HEADER + "\n").encode())
        for iteration in sorted(captured):
            lines = captured[iteration].splitlines(keepends=True)
            stream.writelines(lines[1:])
    log_text = log.read_text(errors="replace") if log.exists() else ""
    nonfinite_console = []
    for line_no, line in enumerate(log_text.splitlines(), 1):
        relevant = ("HUBBARD_MOMENT_ENERGY" in line) or ("Diff=" in line)
        if not relevant:
            continue
        for key, token in re.findall(r"([A-Za-z0-9_]+)=\s*([^\s]+)", line):
            if token.upper().startswith("UNAVAILABLE"):
                continue
            try:
                number = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if not math.isfinite(number):
                nonfinite_console.append({"line": line_no, "field": key, "value": token})
    captured_ids = sorted(captured)
    missing_ids = [i for i in range(captured_ids[0], captured_ids[-1] + 1) if i not in captured] if captured_ids else []
    suffix = []
    if captured_ids:
        for i in range(captured_ids[-1], captured_ids[0] - 1, -1):
            if i not in captured:
                break
            suffix.append(i)
        suffix.reverse()
    endpoint_files = ("Fe_out.nml", "Fe_scf.nml", "Fe_dos.out", "Fe_orbital_dos.out", "totaldos.out", "input.nml")
    manifest = {
        "schema": "c3-rc2-monitored-full-pout-history-v1",
        "runner_source_path": str(RUNNER_SOURCE_PATH), "runner_source_sha256": RUNNER_SOURCE_SHA256,
        "poll_seconds": args.poll_seconds,
        "binary_path": str(args.binary.resolve()), "binary_sha256": sha(args.binary),
        "run_dir": str(run), "input_initial_sha256": input_initial_sha256,
        "seed_initial_sha256": seed_initial_sha256,
        "input_terminal_sha256": sha(run / "input.nml"),
        "Fe_out_terminal_sha256": sha(run / "Fe_out.nml") if (run / "Fe_out.nml").is_file() else None,
        "Fe_scf_terminal_sha256": sha(run / "Fe_scf.nml") if (run / "Fe_scf.nml").is_file() else None,
        "environment": THREADS,
        "trace_settings": {k: env[k] for k in env if k.startswith("RSLMTO_FIXED_U_TRACE") or k in ("RSLMTO_FIXED_U_LAMBDA", "RSLMTO_FIXED_U_RUN_ID")},
        "source_trace_behavior": "full writer rolls per iteration; wrapper archives only complete P_OUT groups with six ql, 50 direct d-matrix entries and two spin W1 traces",
        "captured_iterations": captured_ids, "captured_count": len(captured),
        "missing_iteration_ids_between_first_and_last": missing_ids,
        "contiguous_terminal_suffix": suffix, "contiguous_terminal_suffix_count": len(suffix),
        "per_iteration_matrix_checks": per_iteration_checks,
        "all_matrix_checks_pass": bool(per_iteration_checks) and all(x["pass_1e_11"] for x in per_iteration_checks.values()),
        "nonfinite_console_fields": nonfinite_console,
        "nonfinite_supported_diagnostic_observed": bool(nonfinite_console) or any(
            x.get("nonfinite_values") for x in per_iteration_checks.values()),
        "captured_files_sha256": {f"P_OUT_iter{i:04d}.csv": sha_bytes(b) for i, b in sorted(captured.items())},
        "aggregate_path": str(aggregate), "aggregate_sha256": sha(aggregate),
        "final_trace_path": str(trace), "final_trace_sha256": sha(trace) if trace.exists() else None,
        "console_path": str(log), "console_sha256": sha(log), "returncode": returncode,
        "endpoint_files_sha256": {name: sha(run / name) for name in endpoint_files if (run / name).is_file()},
        "elapsed_seconds": time.time() - start,
        "early_c1_checkpoint": early_checkpoint,
        "intentional_accepted_stop": early_checkpoint is not None,
    }
    (run / "TRACE_ARCHIVE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2))
    if returncode and early_checkpoint is None:
        raise SystemExit(returncode)


if __name__ == "__main__":
    main()
