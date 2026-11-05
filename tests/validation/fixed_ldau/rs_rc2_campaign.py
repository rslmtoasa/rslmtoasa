#!/usr/bin/env python3
"""Gated, bounded C3 RC2 endpoint campaign with the C1 fallback ladder.

The driver refuses to launch unless the canonical RC JSON contains an
    finite RC1 residual evidence at or below the registered 1% limit. A C1
    provisional endpoint may support the C3 RC1 numerical record, but does not
    turn the RC1 gate into PASS; its status is propagated to RC2 provenance.
    It runs no more than ``workers`` SCFs at
once. Each case gets an initial 200-step cap, then at most the two C1
continuations (beta/2, then Broyden beta=.01), each capped at 200 steps.
Nonzero exit or a non-finite supported diagnostic is returned as a C2 triage
event and halts that case; it is never silently treated as ordinary
nonconvergence.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / "tests/validation/fixed_ldau/rs_rc2_monitored_runner.py"
THREADS = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
DIFF_RE = re.compile(r"(?:Not converged!\s*)?Diff=\s*([-+0-9.eEdD]+)")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def potential_ql(path: Path) -> dict[int, list[float]]:
    text = path.read_text(errors="replace")
    out = {}
    for spin, payload in re.findall(r"(?im)^\s*ql\(\s*1\s*,\s*:\s*,\s*([12])\s*\)\s*=\s*([^\n]+)", text):
        vals = [float(x.strip().replace("D", "E").replace("d", "e")) for x in payload.split(",")]
        out[int(spin)] = vals
    if set(out) != {1, 2} or any(len(x) < 3 or not all(math.isfinite(v) for v in x[:3]) for x in out.values()):
        raise ValueError(f"missing/nonfinite saved ql in {path}")
    return out


def gate_check(path: Path) -> dict:
    data = json.loads(path.read_text())
    gate = data.get("rc1", {}).get("final_gate", {})
    status = gate.get("status")
    residuals = gate.get("relative_residuals", [])
    vals = [abs(float(x)) for x in residuals]
    if not vals or not all(math.isfinite(x) for x in vals) or max(vals) > 0.01:
        raise RuntimeError(f"RC2 blocked: final RC1 numeric residual evidence missing/invalid/outside 1%: {residuals!r}")
    if status not in {"PASS", "UNESTABLISHED_PROVISIONAL_RECORD"}:
        raise RuntimeError(f"RC2 blocked: RC1 status is not a measured in-limit result: {status!r}")
    return {"status": status, "relative_residuals": vals, "max_abs_relative_residual": max(vals),
            "source_path": str(path), "source_sha256": sha(path),
            "gate_verdict": "PASS" if status == "PASS" else "UNESTABLISHED",
            "provisional_record_used": status == "UNESTABLISHED_PROVISIONAL_RECORD"}


def c1_acceptance(run: Path, runner_manifest: dict) -> dict:
    if runner_manifest["returncode"] != 0 and not runner_manifest.get("intentional_accepted_stop"):
        return {"accepted": False, "classification": "C2_TRIAGE_REQUIRED_NONZERO_EXIT"}
    if runner_manifest.get("intentional_accepted_stop"):
        checkpoint = runner_manifest.get("early_c1_checkpoint") or {}
        checkpoint_dir = Path(checkpoint.get("checkpoint_dir", ""))
        expected = checkpoint.get("files_sha256", {}).get("Fe_scf.nml")
        actual_path = checkpoint_dir / "Fe_scf.nml"
        if not checkpoint_dir.is_dir() or not expected or not actual_path.is_file() or sha(actual_path) != expected:
            return {"accepted": False, "classification": "INTENTIONAL_STOP_CHECKPOINT_HASH_FAILURE"}
    if runner_manifest["nonfinite_supported_diagnostic_observed"]:
        return {"accepted": False, "classification": "C2_TRIAGE_REQUIRED_NONFINITE"}
    if not runner_manifest["all_matrix_checks_pass"]:
        return {"accepted": False, "classification": "INCOMPLETE_OR_INVALID_CAPTURE"}

    trace = run / "fixed_u_trace_history.csv"
    rows = list(csv.DictReader(trace.open(newline="")))
    by_it: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        if row.get("stage") == "P_OUT":
            by_it.setdefault(int(row["iteration"]), []).append(row)
    ids = sorted(by_it)
    suffix = runner_manifest["contiguous_terminal_suffix"]
    checkpoint = runner_manifest.get("early_c1_checkpoint") or {}
    checkpoint_iteration = checkpoint.get("iteration") if runner_manifest.get("intentional_accepted_stop") else None
    selected_iteration = int(checkpoint_iteration) if checkpoint_iteration is not None else max(ids)
    eligible_ids = [iteration for iteration in ids if iteration <= selected_iteration]
    if (len(eligible_ids) < 20 or eligible_ids[-1] != selected_iteration or
            eligible_ids[-20:] != list(range(selected_iteration - 19, selected_iteration + 1))):
        return {"accepted": False, "classification": "INSUFFICIENT_CONTIGUOUS_20_ITERATION_HISTORY",
                "terminal_suffix_count": len(suffix), "selected_iteration": selected_iteration}
    window = eligible_ids[-20:]
    d = {1: [], 2: []}
    moment = []
    for it in window:
        ql = [r for r in by_it[it] if r.get("kind") == "ql" and r.get("site") == "0"]
        keyed = {(int(r["l"]), int(r["spin"])): float(r["scalar"]) for r in ql}
        if set(keyed) != {(l, s) for l in range(3) for s in (1, 2)} or not all(math.isfinite(v) for v in keyed.values()):
            return {"accepted": False, "classification": "INCOMPLETE_OR_NONFINITE_QL_WINDOW", "iteration": it}
        d[1].append(keyed[(2, 1)])
        d[2].append(keyed[(2, 2)])
        moment.append(abs(sum(keyed[(l, 1)] - keyed[(l, 2)] for l in range(3))))

    log = (run / "console.log").read_text(errors="replace")
    diffs = [float(x.replace("D", "E").replace("d", "e")) for x in DIFF_RE.findall(log)]
    # A live accepted checkpoint is frozen as soon as its Diff and P_OUT
    # record agree. The solver may complete one more P_OUT before the monitor
    # notices the checkpoint condition, leaving a trailing trace row without
    # a corresponding Diff marker. That trailing row is not part of the
    # checkpoint window. Keep strict association through the selected
    # checkpoint, and reject missing/interior Diff markers; ignore only rows
    # strictly after the checkpoint iteration.
    if selected_iteration not in ids or len(diffs) < selected_iteration:
        return {"accepted": False, "classification": "DIFF_TRACE_ITERATION_ASSOCIATION_UNPROVEN",
                "captured_max_iteration": max(ids), "selected_iteration": selected_iteration,
                "diff_marker_count": len(diffs)}
    final_diff = diffs[selected_iteration - 1]
    spread = lambda xs: max(xs) - min(xs)
    spreads = {"d_up": spread(d[1]), "d_down": spread(d[2]), "sphere_moment_muB": spread(moment)}
    accepted = final_diff <= 3e-8 and spreads["d_up"] <= 5e-7 and spreads["d_down"] <= 5e-7 and spreads["sphere_moment_muB"] <= 1e-6
    selected = (Path(runner_manifest["early_c1_checkpoint"]["checkpoint_dir"]) / "Fe_scf.nml"
                if runner_manifest.get("intentional_accepted_stop") else run / "Fe_out.nml")
    endpoint = {}
    if accepted:
        try:
            q = potential_ql(selected)
            q_up, q_dn = q[1][2], q[2][2]
            sphere = abs(sum(q[1][l] - q[2][l] for l in range(3)))
            bounds = {
                "d_up": [min(d[1]), max(d[1])], "d_down": [min(d[2]), max(d[2])],
                "sphere_moment_muB": [min(moment), max(moment)],
            }
            endpoint_in_ranges = (bounds["d_up"][0] <= q_up <= bounds["d_up"][1] and
                           bounds["d_down"][0] <= q_dn <= bounds["d_down"][1] and
                           bounds["sphere_moment_muB"][0] <= sphere <= bounds["sphere_moment_muB"][1])
            endpoint = {"path": str(selected), "sha256": sha(selected),
                        "d_up": q_up, "d_down": q_dn, "sphere_moment_muB": sphere,
                        "inside_last20_ranges_diagnostic_only": endpoint_in_ranges, "last20_ranges": bounds}
        except (OSError, ValueError) as exc:
            accepted = False
            endpoint = {"path": str(selected), "verification_error": str(exc)}
    return {"accepted": accepted, "classification": "C1_ACCEPTED" if accepted else "NOT_C1_ACCEPTED_AT_CAP",
            "window_iterations": window, "final_diff": final_diff, "spreads": spreads,
            "endpoint": endpoint,
            "gate_limits": {"diff": 3e-8, "d_up": 5e-7, "d_down": 5e-7, "sphere_moment_muB": 1e-6}}


def set_mix(text: str, beta: str, mixtype: str | None = None) -> str:
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().lower() == "&mix")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "/")
    found = set()
    for i in range(start + 1, end):
        m = re.match(r"\s*(beta|mixtype)\s*=", lines[i], re.I)
        if m:
            key = m.group(1).lower()
            if key == "beta":
                lines[i] = f"beta = {beta}"
            elif mixtype:
                lines[i] = f"mixtype = '{mixtype}'"
            found.add(key)
    if "beta" not in found:
        lines.insert(end, f"beta = {beta}")
        end += 1
    if mixtype and "mixtype" not in found:
        lines.insert(end, f"mixtype = '{mixtype}'")
    return "\n".join(lines) + "\n"


def prepare_phase(source: Path, dest: Path, beta: str, mixtype: str | None, seed: Path | None) -> None:
    shutil.copytree(source, dest, ignore=shutil.ignore_patterns("trace_snapshots", "fixed_u_trace*", "console.log", "TRACE_ARCHIVE_MANIFEST.json", "Fe_*out*", "Fe_*scf*", "Fe_dos.out", "Fe_orbital_dos.out", "totaldos.out"))
    deck = (dest / "input.nml").read_text()
    (dest / "input.nml").write_text(set_mix(deck, beta, mixtype))
    if seed is not None:
        shutil.copyfile(seed, dest / "Fe_out.nml")


PHASE_DIRS = ("phase0_initial", "phase1_beta_halved", "phase2_broyden_beta001")
PHASE_SETTINGS = (("0.01", None), ("0.005", None), ("0.01", "broyden"))


def runner_processes_for(run_dir: Path) -> list[int]:
    """Return exact monitored-runner PIDs whose argv names this run directory."""
    try:
        proc = subprocess.run(["ps", "-axo", "pid=,command="], check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    except (OSError, subprocess.CalledProcessError):
        return []
    needle = str(run_dir.resolve())
    found = []
    for line in proc.stdout.splitlines():
        if "rs_rc2_monitored_runner.py" not in line or needle not in line:
            continue
        match = re.match(r"\s*(\d+)\s+", line)
        if match:
            found.append(int(match.group(1)))
    return found


def wait_existing_phase(phase_dir: Path, trusted_active_dirs: set[str]) -> dict:
    """Adopt an existing active phase without ever launching into its directory."""
    manifest = phase_dir / "TRACE_ARCHIVE_MANIFEST.json"
    while True:
        if manifest.exists():
            try:
                return json.loads(manifest.read_text())
            except (OSError, json.JSONDecodeError):
                if str(phase_dir.resolve()) not in trusted_active_dirs and not runner_processes_for(phase_dir):
                    raise RuntimeError(f"existing phase manifest is incomplete after runner exit: {phase_dir}")
        elif str(phase_dir.resolve()) not in trusted_active_dirs and not runner_processes_for(phase_dir):
            raise RuntimeError(f"existing phase has no manifest and no exact active runner: {phase_dir}")
        time.sleep(1.0)


def run_one(case: dict, binary: Path, campaign: Path, runner: Path, gate: dict,
            trusted_active_dirs: set[str] | None = None) -> dict:
    name = case["case"]
    base = Path(case["directory"])
    case_dir = campaign / name
    case_dir.mkdir(parents=True, exist_ok=True)
    phases = []
    accepted_phase = None
    for phase_idx, (beta, mix) in enumerate(PHASE_SETTINGS):
        phase_dir = case_dir / PHASE_DIRS[phase_idx]
        if phase_dir.exists():
            adopted_existing = True
            manifest_path = phase_dir / "TRACE_ARCHIVE_MANIFEST.json"
            if manifest_path.exists():
                m = json.loads(manifest_path.read_text())
            else:
                m = wait_existing_phase(phase_dir, trusted_active_dirs or set())
                manifest_path = phase_dir / "TRACE_ARCHIVE_MANIFEST.json"
            proc_returncode = m.get("returncode")
            proc_stdout = "adopted existing run; no solver relaunched"
        else:
            adopted_existing = False
            if phase_idx == 0:
                seed = base / "Fe_out.nml"
            else:
                previous_record = phases[-1]
                previous_phase = Path(previous_record["run_dir"])
                previous_solver = previous_record.get("solver_manifest", {})
                previous_checkpoint = previous_solver.get("early_c1_checkpoint")
                if previous_solver.get("intentional_accepted_stop") and previous_checkpoint:
                    seed = Path(previous_checkpoint["checkpoint_dir"]) / "Fe_scf.nml"
                    expected_sha = previous_checkpoint.get("files_sha256", {}).get("Fe_scf.nml")
                    if not seed.is_file() or sha(seed) != expected_sha:
                        phases.append({"phase": phase_idx, "status": "CANNOT_CONTINUE_FROZEN_CHECKPOINT_HASH_MISMATCH"})
                        break
                else:
                    seed = previous_phase / "Fe_scf.nml"
                if not seed.is_file():
                    phases.append({"phase": phase_idx, "status": "CANNOT_CONTINUE_TERMINAL_FE_SCF_MISSING"})
                    break
            prepare_phase(base, phase_dir, beta, mix, seed)
            proc = subprocess.run([sys.executable, str(runner), "--binary", str(binary), "--run-dir", str(phase_dir), "--stop-on-c1"],
                                  cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            manifest_path = phase_dir / "TRACE_ARCHIVE_MANIFEST.json"
            if not manifest_path.exists():
                phases.append({"phase": phase_idx, "run_dir": str(phase_dir), "status": "RUNNER_FAILED_NO_MANIFEST",
                               "runner_returncode": proc.returncode, "stdout": proc.stdout[-2000:]})
                break
            m = json.loads(manifest_path.read_text())
            proc_returncode = proc.returncode
            proc_stdout = proc.stdout[-2000:]
        acceptance = c1_acceptance(phase_dir, m)
        record = {"phase": phase_idx, "run_dir": str(phase_dir), "beta": beta,
                  "mixtype": mix or "inherited", "runner_returncode": proc_returncode,
                  "adopted_existing_phase": adopted_existing,
                  "trace_manifest_path": str(manifest_path), "trace_manifest_sha256": sha(manifest_path),
                  "acceptance": acceptance, "solver_manifest": m}
        phases.append(record)
        if acceptance.get("accepted") and accepted_phase is None:
            accepted_phase = phase_idx
        if proc_returncode not in (0, None) and not m.get("intentional_accepted_stop"):
            break
        if acceptance["classification"].startswith("C2_TRIAGE_REQUIRED"):
            break
        # If an accepted state already exists, only adopt any parser-debug
        # continuations that were launched by the old in-memory campaign.
        # Never launch an additional continuation from an accepted phase.
        if accepted_phase is not None:
            next_phase_exists = phase_idx < 2 and (case_dir / PHASE_DIRS[phase_idx + 1]).exists()
            if not next_phase_exists:
                break
    if accepted_phase is not None:
        final_status = "C1_ACCEPTED"
    elif len(phases) == 3:
        final_status = "PROVISIONAL_AFTER_TWO_C1_CONTINUATIONS"
    else:
        final_status = phases[-1].get("acceptance", {}).get("classification", phases[-1].get("status"))
    result = {"case": name, "family": case["family"], "mode": case["mode"], "u_eV": case["u_eV"],
              "j_eV": case["j_eV"], "status": final_status,
              "phases": phases, "gate_provenance": gate}
    (case_dir / "CAMPAIGN_CASE.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rc-json", type=Path, required=True)
    ap.add_argument("--staged-manifest", type=Path, required=True)
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--resume", action="store_true", help="adopt existing phase directories and continue only missing phases")
    ap.add_argument("--adopt-active-run-dir", type=Path, action="append", default=[],
                    help="exact existing runner directory independently verified live before scheduler replacement")
    ap.add_argument("--launch", action="store_true", help="required to actually launch the gated campaign")
    args = ap.parse_args()
    if not args.launch:
        raise SystemExit("refusing to launch without --launch")
    gate = gate_check(args.rc_json)
    staged = json.loads(args.staged_manifest.read_text())
    if staged.get("case_count") != 30 or len(staged.get("cases", [])) != 30:
        raise SystemExit("staged RC2 campaign must contain exactly the registered 30 cases")
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("workers must be in 1..4")
    if args.outdir.exists() and not args.resume:
        raise SystemExit(f"refusing to overwrite campaign output {args.outdir}")
    if not args.binary.is_file():
        raise SystemExit(f"binary missing: {args.binary}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.outdir / "RC2_CAMPAIGN.json"
    if args.resume:
        if not manifest_path.is_file():
            raise SystemExit(f"cannot resume without existing campaign manifest: {manifest_path}")
        campaign = json.loads(manifest_path.read_text())
        if campaign.get("staged_manifest_sha256") != sha(args.staged_manifest) or campaign.get("binary_sha256") != sha(args.binary):
            raise SystemExit("resume provenance mismatch: staged manifest or solver binary differs")
        campaign.setdefault("resume_history", []).append({"unix": time.time(), "prior_manifest_sha256": sha(manifest_path),
                                                            "runner_path": str(RUNNER), "runner_sha256": sha(RUNNER),
                                                            "adopt_existing_phase_dirs": True})
        campaign["runner_path"] = str(RUNNER)
        campaign["runner_sha256"] = sha(RUNNER)
        campaign["status"] = "RUNNING_RESUMED"
    else:
        campaign = {"schema": "c3-rc2-campaign-v1", "status": "RUNNING", "created_unix": time.time(),
                    "rc1_gate": gate, "staged_manifest_path": str(args.staged_manifest.resolve()),
                    "staged_manifest_sha256": sha(args.staged_manifest), "binary_path": str(args.binary.resolve()),
                    "binary_sha256": sha(args.binary), "runner_path": str(RUNNER), "runner_sha256": sha(RUNNER),
                    "workers": args.workers, "thread_caps": THREADS, "cases": []}
    manifest_path.write_text(json.dumps(campaign, indent=2, sort_keys=True) + "\n")
    failures = []
    trusted_active_dirs = {str(p.resolve()) for p in args.adopt_active_run_dir}
    existing_by_case = {entry.get("case"): entry for entry in campaign.get("cases", [])}
    campaign["cases"] = [existing_by_case[k] for k in sorted(existing_by_case)]
    # Fail closed before scheduling any new SCF if an existing directory is
    # incomplete and was not independently verified live for this resume.
    # This prevents an adoption/inspection error from freeing a worker slot
    # and silently launching a second solver into an occupied workflow.
    if args.resume:
        unknown_incomplete = []
        for case in staged["cases"]:
            case_dir = args.outdir / case["case"]
            for phase in PHASE_DIRS:
                phase_dir = case_dir / phase
                if phase_dir.is_dir() and not (phase_dir / "TRACE_ARCHIVE_MANIFEST.json").is_file():
                    if str(phase_dir.resolve()) not in trusted_active_dirs:
                        unknown_incomplete.append(str(phase_dir.resolve()))
        if unknown_incomplete:
            raise SystemExit("resume refused: incomplete existing phase(s) lack exact live-run receipt: " +
                             ", ".join(unknown_incomplete))
    cases_to_run = list(staged["cases"])
    if args.resume:
        def resume_priority(case: dict) -> tuple[int, str]:
            case_dir = args.outdir / case["case"]
            active_existing = any((case_dir / phase / "TRACE_ARCHIVE_MANIFEST.json").exists() is False and
                                  (case_dir / phase).is_dir() and str((case_dir / phase).resolve()) in trusted_active_dirs
                                  for phase in PHASE_DIRS)
            any_existing = any((case_dir / phase).is_dir() for phase in PHASE_DIRS)
            return (0 if active_existing else (1 if any_existing else 2), case["case"])
        cases_to_run.sort(key=resume_priority)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, c, args.binary.resolve(), args.outdir.resolve(), RUNNER, gate, trusted_active_dirs): c for c in cases_to_run}
        for future in concurrent.futures.as_completed(futures):
            case = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"case": case["case"], "status": "CAMPAIGN_DRIVER_ERROR", "error": repr(exc)}
                failures.append(result)
            by_case = {entry.get("case"): entry for entry in campaign["cases"]}
            by_case[result["case"]] = result
            campaign["cases"] = [by_case[k] for k in sorted(by_case)]
            campaign["status"] = "RUNNING"
            manifest_path.write_text(json.dumps(campaign, indent=2, sort_keys=True) + "\n")
    campaign["status"] = "COMPLETE_WITH_ERRORS" if failures else "COMPLETE"
    campaign["completed_unix"] = time.time()
    manifest_path.write_text(json.dumps(campaign, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": campaign["status"], "cases": len(campaign["cases"]), "manifest": str(manifest_path)}, indent=2))
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
