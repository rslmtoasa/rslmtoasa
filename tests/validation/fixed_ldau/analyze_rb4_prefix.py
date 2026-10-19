#!/usr/bin/env python3
"""Read-only RB4 restart-prefix evaluator using C1 full-precision QL traces.

Usage: analyze_rb4_prefix.py ACCEPTED_PREVIOUS_RUN RESTART_RUN [--site N]

Each run directory must contain fixed_u_trace.csv and a console log named
c1_console.log, run.log, or testrun.log. Only QL rows at P_OUT are consumed.
The previous run must itself pass C1 RS acceptance on its final 20 completed
iterations. For each k=1..10 the script joins previous last (20-k) samples
with restart first k samples, and reports all metrics and the earliest pass.
This is analysis tooling: it does not inspect or modify inputs, run binaries,
or certify that the two runs share a valid matched recipe/restart state.
"""
import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

DIFF_LIMIT = 3.0e-8
D_SPREAD_LIMIT = 5.0e-7
MOM_SPREAD_LIMIT = 1.0e-6
LOG_EVENT = re.compile(
    r"Perform recursion at step\s+(\d+)|"
    r"Not converged! Diff=\s*([-+0-9.eEdD]+)|"
    r"Converged!\s*([-+0-9.eEdD]+)"
)


def number(s):
    x = float(s.replace("D", "E").replace("d", "e"))
    if not math.isfinite(x):
        raise ValueError(f"non-finite numeric value: {s!r}")
    return x


def read_diffs(run_dir):
    logs = [run_dir / n for n in ("c1_console.log", "run.log", "testrun.log")]
    log = next((p for p in logs if p.is_file()), None)
    if log is None:
        raise ValueError(f"no console log found in {run_dir}")
    diffs = {}
    pending_step = None
    for m in LOG_EVENT.finditer(log.read_text(errors="replace")):
        if m.group(1) is not None:
            if pending_step is not None:
                raise ValueError(f"step {pending_step} has no Diff before next step in {log}")
            pending_step = int(m.group(1))
        else:
            if pending_step is None:
                # Ignore matching text outside a recursion (e.g. summaries).
                continue
            if pending_step in diffs:
                raise ValueError(f"duplicate Diff for step {pending_step} in {log}")
            diffs[pending_step] = number(m.group(2) or m.group(3))
            pending_step = None
    # A final currently-running/uncompleted step may lack a Diff. It must not
    # be mistaken for a completed iteration; QL/log intersection handles it.
    return diffs, str(log)


def read_run(path, site):
    run_dir = Path(path).resolve()
    trace = run_dir / "fixed_u_trace.csv"
    if not trace.is_file():
        raise ValueError(f"missing {trace}")
    diffs, log_path = read_diffs(run_dir)
    by_iter = {}
    required = {(l, s) for l in range(3) for s in (1, 2)}
    with trace.open(newline="") as f:
        reader = csv.DictReader(f)
        needed = {"kind", "iteration", "stage", "site", "l", "spin", "scalar"}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise ValueError(f"{trace}: header missing {sorted(needed - set(reader.fieldnames or []))}")
        for line, row in enumerate(reader, 2):
            if row["kind"].strip() != "ql" or row["stage"].strip() != "P_OUT":
                continue
            if int(row["site"]) != site:
                continue
            it, l, spin = int(row["iteration"]), int(row["l"]), int(row["spin"])
            if l not in range(3) or spin not in (1, 2):
                raise ValueError(f"{trace}:{line}: unexpected QL indices l={l}, spin={spin}")
            vals = by_iter.setdefault(it, {})
            key = (l, spin)
            if key in vals:
                raise ValueError(f"{trace}:{line}: duplicate P_OUT QL for iteration={it}, site={site}, l={l}, spin={spin}")
            vals[key] = number(row["scalar"])
    if not by_iter:
        raise ValueError(f"{trace}: no P_OUT QL rows for site {site}")
    complete = {}
    unpaired_trace = []
    for it, vals in by_iter.items():
        missing = required - vals.keys()
        if missing:
            raise ValueError(f"{trace}: incomplete iteration {it}, missing {sorted(missing)}")
        if it not in diffs:
            unpaired_trace.append(it)
        else:
            complete[it] = vals
    # Some historical traces were captured just after P_OUT but before the
    # console emitted the corresponding Diff (e.g. accepted RB2 step 32 has
    # an extra P_OUT row at step 33). Preserve those bytes for provenance but
    # exclude them from iteration counts/windows. Only a contiguous suffix is
    # admissible; an interior unpaired P_OUT sample still signals corruption.
    if unpaired_trace:
        paired_max = max(complete) if complete else -1
        if any(it <= paired_max for it in unpaired_trace):
            raise ValueError(f"{trace}: unpaired P_OUT iterations are not trailing: {sorted(unpaired_trace)}")
        all_ids = sorted(by_iter)
        if any(b != a + 1 for a, b in zip(all_ids, all_ids[1:])):
            raise ValueError(f"{trace}: P_OUT iteration sequence is not contiguous: {all_ids[:3]}..{all_ids[-3:]}")
    orphan = set(diffs) - set(complete)
    # Ignore only a single trailing log Diff without QL if an incomplete
    # iteration was not recorded; interior orphan events indicate corruption.
    if orphan:
        max_trace = max(complete)
        if any(i <= max_trace for i in orphan):
            raise ValueError(f"{log_path}: Diff steps without matching P_OUT QL: {sorted(orphan)}")
    its = sorted(complete)
    if any(b != a + 1 for a, b in zip(its, its[1:])):
        raise ValueError(f"{trace}: P_OUT iteration sequence is not contiguous: {its[:3]}..{its[-3:]}")
    if len(its) < 1:
        raise ValueError(f"{trace}: no complete iterations")
    return {"path": str(run_dir), "trace": str(trace), "log": log_path,
            "iterations": its, "ql": complete, "diff": diffs}


def observables(sample):
    up = sample[(2, 1)]
    down = sample[(2, 2)]
    moment = abs(sum(sample[(l, 1)] for l in range(3)) -
                 sum(sample[(l, 2)] for l in range(3)))
    return up, down, moment


def summarize(samples, diff, refs):
    obs = [observables(s) for s in samples]
    ups, downs, moments = ([x[i] for x in obs] for i in range(3))
    spreads = [max(x) - min(x) for x in (ups, downs, moments)]
    ok = (len(samples) == 20 and diff <= DIFF_LIMIT and
          spreads[0] <= D_SPREAD_LIMIT and spreads[1] <= D_SPREAD_LIMIT and
          spreads[2] <= MOM_SPREAD_LIMIT)
    return {
        "n": len(samples), "iteration_refs": refs,
        "diff": diff, "d_up_final": ups[-1], "d_down_final": downs[-1],
        "sphere_moment_muB_final": moments[-1],
        "spread_d_up": spreads[0], "spread_d_down": spreads[1],
        "spread_sphere_moment_muB": spreads[2], "qualifies": ok,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("accepted_previous_run")
    ap.add_argument("restart_run")
    ap.add_argument("--site", type=int, default=0)
    args = ap.parse_args()
    if args.site < 0:
        ap.error("--site must be nonnegative")
    prev = read_run(args.accepted_previous_run, args.site)
    cur = read_run(args.restart_run, args.site)
    pi, ci = prev["iterations"], cur["iterations"]
    if len(pi) < 20:
        raise ValueError(f"preceding history has only {len(pi)} complete iterations; need >=20")
    prev_ids = pi[-20:]
    prev_tail = [prev["ql"][i] for i in prev_ids]
    prev_accept = summarize(prev_tail, prev["diff"][prev_ids[-1]],
                            [f"previous:{i}" for i in prev_ids])
    if not prev_accept["qualifies"]:
        raise ValueError("preceding run is not C1 RS-accepted on its final 20 samples; "
                         "refusing RB4 join: " + json.dumps(prev_accept, sort_keys=True))
    candidates = []
    for k in range(1, min(10, len(ci)) + 1):
        pids = prev_ids[-(20-k):] if 20-k else []
        cids = ci[:k]
        joined = [prev["ql"][i] for i in pids] + [cur["ql"][i] for i in cids]
        refs = [f"previous:{i}" for i in pids] + [f"restart:{i}" for i in cids]
        item = summarize(joined, cur["diff"][cids[-1]], refs)
        item["k"] = k
        item["current_iteration"] = cids[-1]
        candidates.append(item)
    first = next((x["k"] for x in candidates if x["qualifies"]), None)
    result = {
        "analyzer": "RB4 joined-prefix C1 RS analyzer",
        "thresholds": {"diff_max": DIFF_LIMIT, "d_spin_spread_max": D_SPREAD_LIMIT,
                       "sphere_moment_spread_muB_max": MOM_SPREAD_LIMIT,
                       "window_n": 20, "max_restart_steps": 10},
        "site_index": args.site,
        "previous": {k: prev[k] for k in ("path", "trace", "log")},
        "restart": {k: cur[k] for k in ("path", "trace", "log")},
        "previous_final20_acceptance": prev_accept,
        "restart_complete_iterations": len(ci),
        "candidate_prefixes": candidates,
        "first_qualifying_k": first,
        "RB4_restart_pass": first is not None,
        "warning": "Does not validate matched recipe, input/seed identity, executable identity, or serial/MPI pairing.",
    }
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if first is not None else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, csv.Error) as e:
        print(json.dumps({"error": str(e)}, indent=2), file=sys.stderr)
        sys.exit(2)
