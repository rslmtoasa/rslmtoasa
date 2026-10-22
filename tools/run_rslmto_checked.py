#!/usr/bin/env python3
"""Run one solver process only after the repository C5.4 preflight passes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT = ROOT / "docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py"
DEFAULT_GUARD = ROOT / "tests/validation/fixed_ldau/c5_7_launch.sh"
THREAD_CAPS = {
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "RSLMTO_OMP_THREADS_SERIAL": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".preflight.", suffix=".json", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def cli() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="existing writable SCF run directory")
    parser.add_argument("--binary", required=True, type=Path, help="rslmto.x executable")
    parser.add_argument("--deck", type=Path, help="deck to preflight (default: RUN_DIR/input.nml)")
    parser.add_argument("--restart", type=Path, help="optional writable restart potential for C5.4 checks")
    parser.add_argument("--guard", type=Path, default=DEFAULT_GUARD,
                        help=f"thread-capped bash launcher checked and used (default: {DEFAULT_GUARD})")
    parser.add_argument("--host-ef-trace", type=Path,
                        help="optional host fixed_u_trace.csv with one fermi_level row")
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    deck = (args.deck if args.deck is not None else run_dir / "input.nml").resolve()
    binary = args.binary.resolve()
    guard = args.guard.resolve()
    restart = args.restart.resolve() if args.restart else None
    host_trace = args.host_ef_trace.resolve() if args.host_ef_trace else None
    record_path = run_dir / "preflight.json"
    record: dict[str, Any] = {
        "schema": "run_rslmto_checked_v1",
        "status": "PREFLIGHT_RUNNING",
        "run_dir": str(run_dir),
        "deck": str(deck),
        "deck_sha256": None,
        "restart": str(restart) if restart else None,
        "restart_sha256": None,
        "host_ef_trace": str(host_trace) if host_trace else None,
        "host_ef_trace_sha256": None,
        "binary": str(binary),
        "binary_sha256": None,
        "guard": str(guard),
        "preflight_script": str(PREFLIGHT),
        "preflight_python": sys.executable,
        "thread_caps": THREAD_CAPS,
        "launch": {"status": "NOT_STARTED"},
    }

    # Invalidate any earlier PASS record before even attempting preflight. If a
    # new record cannot be written, the solver is never invoked.
    try:
        write_json_atomic(record_path, record)
    except OSError as exc:
        print(f"cannot create {record_path}: {exc}", file=sys.stderr)
        return 2

    try:
        if Path(sys.prefix).resolve() != (ROOT / ".venv").resolve():
            raise RuntimeError("run this command with the repository .venv Python interpreter")
        if not run_dir.is_dir():
            raise RuntimeError(f"run directory does not exist: {run_dir}")
        if not PREFLIGHT.is_file():
            raise RuntimeError(f"C5.4 preflight script missing: {PREFLIGHT}")
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise RuntimeError(f"solver binary missing or not executable: {binary}")
        for name, path in (("deck", deck), ("restart", restart), ("host_ef_trace", host_trace), ("guard", guard), ("binary", binary)):
            if path is not None and path.is_file():
                record[f"{name}_sha256"] = sha256(path)

        command = [sys.executable, str(PREFLIGHT), str(run_dir), "--deck", str(deck), "--guard", str(guard)]
        if restart:
            command += ["--restart", str(restart)]
        if host_trace:
            command += ["--host-ef-trace", str(host_trace)]
        record["preflight_command"] = command
        preflight = subprocess.run(command, text=True, capture_output=True, check=False, cwd=ROOT)
        try:
            parsed = json.loads(preflight.stdout)
        except json.JSONDecodeError:
            parsed = None
        record["preflight_returncode"] = preflight.returncode
        record["preflight_stdout"] = preflight.stdout
        record["preflight_stderr"] = preflight.stderr
        record["preflight_result"] = parsed
        matched = isinstance(parsed, dict) and parsed.get("status") == "PASS"
        if matched:
            matched = Path(parsed.get("deck", "")).resolve() == deck
            matched = matched and Path(parsed.get("fresh_run_dir", "")).resolve() == run_dir
            if host_trace:
                host_check = parsed.get("host_ef_check")
                matched = matched and isinstance(host_check, dict)
                matched = matched and Path(host_check.get("host_ef_trace", "")).resolve() == host_trace
        if preflight.returncode != 0 or not matched:
            record["status"] = "PREFLIGHT_FAILED"
            record["launch"] = {"status": "REFUSED"}
            write_json_atomic(record_path, record)
            print(f"preflight failed; solver was not launched (see {record_path})", file=sys.stderr)
            if preflight.stderr:
                print(preflight.stderr.rstrip(), file=sys.stderr)
            if preflight.stdout:
                print(preflight.stdout.rstrip(), file=sys.stderr)
            return 1

        # Bind the pass to the exact inputs that will be used. Refuse any input
        # changed between the C5.4 check and launch.
        for name, path in (("deck", deck), ("restart", restart), ("host_ef_trace", host_trace), ("guard", guard), ("binary", binary)):
            if path is not None:
                current_hash = sha256(path)
                if current_hash != record.get(f"{name}_sha256"):
                    raise RuntimeError(f"{name} changed after preflight; solver was not launched")

        record["status"] = "PREFLIGHT_PASSED"
        record["launch"] = {"status": "STARTED", "command": [str(guard), str(binary)]}
        write_json_atomic(record_path, record)
        env = os.environ.copy()
        env.update(THREAD_CAPS)
        with (run_dir / "checked_launcher.log").open("w", encoding="utf-8") as log:
            launched = subprocess.run(["/bin/bash", str(guard), str(binary)], cwd=run_dir,
                                      env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        record["launch"] = {"status": "COMPLETED", "returncode": launched.returncode,
                            "log": str(run_dir / "checked_launcher.log")}
        record["status"] = "RUN_COMPLETED" if launched.returncode == 0 else "RUN_FAILED"
        write_json_atomic(record_path, record)
        return launched.returncode
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        record["status"] = "PREFLIGHT_FAILED"
        record["error"] = str(exc)
        record["launch"] = {"status": "REFUSED"}
        try:
            write_json_atomic(record_path, record)
        except OSError:
            pass
        print(f"{exc}; solver was not launched", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(cli())
