#!/usr/bin/env python3
"""Write tests/test_inventory.csv: one row per registered test.

usage: make_inventory.py BUILD_DIR RUNTIMES.txt TRACE_DIR [OUT.csv]

BUILD_DIR   configured build directory (`ctest --show-only=json-v1` gives names, labels, commands)
RUNTIMES    lines `<test name> <seconds>` (from a ctest run)
TRACE_DIR   one `<test name>.paths` per test: the repository files the test opened, from
            `strace -f -e trace=%file -o raw ctest -R '^<name>$'` filtered to this repository
            (a test without a trace file lists the path arguments of its command that exist)
Columns: name, tier, labels (without the tier), runtime_s, family, reads (`;`-separated repo-relative paths).
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TIERS = ("quick", "standard", "validation")


def family(name: str, labels: list[str]) -> str:
    if "lda_u" in labels or name.startswith("LdaU"):
        return "LDA+U"
    if "gbt" in labels:
        return "GBT"
    if "spin-response" in labels:
        return "spin response"
    if "xc" in labels:
        return "XC"
    if name.startswith(("Regression_", "Triad_")):
        return "backend regression"
    if name.startswith("Val"):
        return "validation"
    if name.startswith(("Example_", "App")):
        return "SCF example"
    return "unit core"


def reads_of(name: str, command: list[str], trace_dir: Path) -> list[str]:
    paths: set[str] = set()
    trace = trace_dir / f"{name}.paths"
    if trace.exists():
        for line in trace.read_text().split():
            paths.add(line)
    else:
        for arg in command:
            if arg.startswith(str(ROOT)) and Path(arg).exists():
                paths.add(arg)
    out = set()
    for p in paths:
        try:
            rel = Path(p).resolve().relative_to(ROOT)
        except ValueError:
            continue
        if rel.parts and rel.parts[0] in {".git", "build", "__pycache__"} or "__pycache__" in rel.parts:
            continue
        if (ROOT / rel).exists():
            out.add(str(rel))
    return sorted(out)


def main() -> int:
    build, runtimes, trace_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out = Path(sys.argv[4]) if len(sys.argv) > 4 else ROOT / "tests/test_inventory.csv"
    times = dict(line.split() for line in runtimes.read_text().splitlines() if line.strip())
    listing = json.loads(subprocess.check_output(["ctest", "--test-dir", str(build), "--show-only=json-v1"], text=True))
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "tier", "labels", "runtime_s", "family", "reads"])
        for test in listing["tests"]:
            name = test["name"]
            labels = next((p["value"] for p in test.get("properties", []) if p["name"] == "LABELS"), [])
            tier = [t for t in labels if t in TIERS]
            if len(tier) != 1:
                raise SystemExit(f"{name}: expected exactly one tier label, got {tier}")
            rest = [l for l in labels if l not in TIERS]
            writer.writerow([name, tier[0], ";".join(rest), f"{float(times.get(name, 0)):.2f}", family(name, rest),
                             ";".join(reads_of(name, test.get("command", []), trace_dir))])
    return 0


if __name__ == "__main__":
    sys.exit(main())
