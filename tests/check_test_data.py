#!/usr/bin/env python3
"""Every file or directory that a registered test reads (tests/test_inventory.csv, column `reads`) must exist.

Bug caught: a data or deck deletion that a retained test still needs, found at the next quick run instead of
when that test runs. The inventory is written by tests/make_inventory.py.

usage: check_test_data.py [inventory.csv]
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    inventory = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "tests/test_inventory.csv"
    missing: list[tuple[str, str]] = []
    n_tests = n_paths = 0
    with inventory.open(newline="") as handle:
        for row in csv.DictReader(handle):
            n_tests += 1
            for rel in filter(None, row["reads"].split(";")):
                n_paths += 1
                if not (ROOT / rel).exists():
                    missing.append((row["name"], rel))
    for name, rel in missing:
        print(f"MISSING {rel}  (read by {name})")
    print(f"{n_tests} tests, {n_paths} paths, {len(missing)} missing")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
